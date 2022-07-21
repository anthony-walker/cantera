//! @file IdealGasConstPressureMoleReactor.cpp A constant pressure
//! zero-dimensional reactor with moles as the state

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasConstPressureMoleReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/utilities.h"
#include <limits>

using namespace std;

namespace Cantera
{

void IdealGasConstPressureMoleReactor::setThermoMgr(ThermoPhase& thermo)
{
    if (thermo.type() != "IdealGas") {
        throw CanteraError("IdealGasConstPressureMoleReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    MoleReactor::setThermoMgr(thermo);
}

void IdealGasConstPressureMoleReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("IdealGasConstPressureMoleReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);
    // get mass for calculations
    m_mass = m_thermo->density() * m_vol;
    // set the first component to the temperature
    y[0] = m_thermo->temperature();
    // Use inverse molecular weights
    const double* Y = m_thermo->massFractions();
    const vector_fp& imw = m_thermo->inverseMolecularWeights();
    double *ys = y + m_sidx;
    for (size_t i = 0; i < m_nsp; i++) {
        ys[i] = m_mass * imw[i] * Y[i];
    }
    // set the remaining components to the surface species moles on the walls
    getSurfaceInitialConditions(y + m_nsp + m_sidx);
}

void IdealGasConstPressureMoleReactor::initialize(double t0)
{
    MoleReactor::initialize(t0);
    m_nv -= 1; // const pressure system loses 1 more variable from MoleReactor
    m_hk.resize(m_nsp, 0.0);
}

void IdealGasConstPressureMoleReactor::updateState(double* y)
{
    // the components of y are: [0] the temperature, [1...K+1) are the
    // moles of each species, and [K+1...] are the moles of surface
    // species on each wall.
    m_thermo->setMolesNoTruncate(y + m_sidx);
    m_thermo->setState_TP(y[0], m_pressure);
    // calculate mass from moles
    const vector_fp& mw = m_thermo->molecularWeights();
    m_mass = 0;
    for (size_t i = 0; i < m_nsp; i++) {
        m_mass += y[i + m_sidx] * mw[i];
    }
    m_vol = m_mass / m_thermo->density();
    updateConnected(false);
    updateSurfaceState(y + m_nsp + m_sidx);
}

void IdealGasConstPressureMoleReactor::eval(double time, double* LHS, double* RHS)
{
    double& mcpdTdt = RHS[0]; // m * c_p * dT/dt
    double* dndt = RHS + m_sidx; // kmol per s

    evalWalls(time);

    m_thermo->restoreState(m_state);

    m_thermo->getPartialMolarEnthalpies(&m_hk[0]);
    const vector_fp& imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // evaluate reactor surfaces
    evalSurfaces(LHS + m_nsp + m_sidx, RHS + m_nsp + m_sidx, m_sdot.data());

    // external heat transfer
    mcpdTdt += m_Qdot;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcpdTdt -= m_wdot[n] * m_hk[n] * m_vol;
        mcpdTdt -= m_sdot[n] * m_hk[n];
        // production in gas phase and from surfaces
        dndt[n] = m_wdot[n] * m_vol + m_sdot[n];
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species into system and dilution by other species
            dndt[n] -= outlet->outletSpeciesMassFlowRate(n) * imw[n];
        }
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        mcpdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dndt[n] += inlet->outletSpeciesMassFlowRate(n) * imw[n];
            mcpdTdt -= m_hk[n] * imw[n] * mdot_spec;
        }
    }

    if (m_energy) {
        LHS[0] = m_mass * m_thermo->cp_mass();
    } else {
        RHS[0] = 0.0;
    }
}

Eigen::SparseMatrix<double> IdealGasConstPressureMoleReactor::jacobian(double t,
    double* y)
{
    // clear former jacobian elements
    m_jac_trips.clear();
    // Determine Species Derivatives
    // volume / moles * rates portion of equation
    size_t nspecies = m_thermo->nSpecies();
    // create sparse structures for rates and volumes
    Eigen::SparseMatrix<double> rates(nspecies, 1);
    Eigen::SparseMatrix<double> volumes(1, nspecies);
    // reserve space for data
    rates.reserve(nspecies);
    volumes.reserve(nspecies);
    // fill sparse structures
    m_kin->getNetProductionRates(rates.valuePtr()); // "omega dot"
    std::fill(volumes.valuePtr(), volumes.valuePtr() + nspecies, m_vol);
    // get ROP derivatives
    double scalingFactor = m_vol/accumulate(y + m_sidx, y + m_nv, 0.0);
    Eigen::SparseMatrix<double> speciesDervs = m_kin->netProductionRates_ddX();
    // sum parts
    speciesDervs = scalingFactor * speciesDervs + rates * volumes;
    // add elements to jacobian triplets
    for (int k=0; k<speciesDervs.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(speciesDervs, k); it; ++it) {
            m_jac_trips.emplace_back(it.row() + m_sidx, it.col() + m_sidx, it.value());
        }
    }
    // Temperature Derivatives
    if (m_analytical_temperature) {
        analyticalDerivs(speciesDervs, t, y);
    } else {
        finiteDiffDerivs(speciesDervs, t, y);
    }
    Eigen::SparseMatrix<double> jac (m_nv, m_nv);
    jac.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    return jac;
}

void IdealGasConstPressureMoleReactor::finiteDiffDerivs(Eigen::SparseMatrix<double>& speciesDervs, double t, double* y)
{
    if (m_energy) {
        // getting perturbed state for finite difference
        double deltaTemp = y[0] * std::numeric_limits<double>::epsilon();
        // finite difference temperature derivatives
        vector_fp yNext(m_nv);
        vector_fp ydotNext(m_nv);
        vector_fp yCurrent(m_nv);
        vector_fp ydotCurrent(m_nv);
        // copy LHS to current and next
        copy(y, y + m_nv, yCurrent.begin());
        copy(y, y + m_nv, yNext.begin());
        // perturb temperature
        yNext[0] += deltaTemp;
        // getting perturbed state
        updateState(yNext.data());
        eval(t, yNext.data(), ydotNext.data());
        ydotNext[0] /= yNext[0];
        // reset and get original state
        updateState(yCurrent.data());
        eval(t, yCurrent.data(), ydotCurrent.data());
        ydotCurrent[0] /= yCurrent[0];
        // d T_dot/dT
        m_jac_trips.emplace_back(0, 0, (ydotNext[0] - ydotCurrent[0]) / deltaTemp);
        // d omega_dot_j/dT
        for (size_t j = m_sidx; j < m_nv; j++) {
            m_jac_trips.emplace_back(j, 0, (ydotNext[j] - ydotCurrent[j]) / deltaTemp);
        }
        // d T_dot/dnj
        vector_fp specificHeat(m_nsp);
        vector_fp netProductionRates(m_nsp);
        vector_fp enthalpy(m_nsp);
        // getting physical quantities
        m_thermo->getPartialMolarCp(specificHeat.data());
        m_thermo->getPartialMolarEnthalpies(enthalpy.data());
        m_kin->getNetProductionRates(netProductionRates.data());
        // getting perturbed changes w.r.t temperature
        double hkndotksum = 0;
        double NtotalCp = accumulate(y + m_sidx, y + m_nv, 0.0) * m_thermo->cp_mole();
        // scale net production rates by  volume to get molar rate
        scale(netProductionRates.begin(), netProductionRates.end(),
            netProductionRates.begin(), m_vol);
        // determine a sum in derivative
        for (size_t i = 0; i < m_nsp; i++) {
            hkndotksum += enthalpy[i] * netProductionRates[i];
        }
        // determine derivatives
        // spans columns
        for (size_t j = 0; j < m_nsp; j++) {
            double hkdnkdnjSum = 0;
            // spans rows
            for (size_t k = 0; k < m_nsp; k++) {
                hkdnkdnjSum += enthalpy[k] * speciesDervs.coeff(k, j);
            }
            // add elements to jacobian triplets
            m_jac_trips.emplace_back(0, j + m_sidx, (-hkdnkdnjSum * NtotalCp +
            specificHeat[j] * hkndotksum) / (NtotalCp * NtotalCp));
        }
    }
}

void IdealGasConstPressureMoleReactor::analyticalDerivs(Eigen::SparseMatrix<double>& speciesDervs, double t, double* y)
{
    // analytical temperature derivatives
    // domega_dot_j/dT
    vector_fp dndotdT(m_nsp);
    m_kin->getNetProductionRates_ddT(dndotdT.data());
    // scale to get molar production rate derivatives
    scale(dndotdT.begin(), dndotdT.end(), dndotdT.begin(), m_vol);
    for (size_t j = m_sidx; j < m_nv; j++) {
        m_jac_trips.emplace_back(j, 0, dndotdT[j]);
    }
    // Energy derivatives
    if (m_energy) {
        // d T_dot/dnj
        vector_fp specificHeat (m_nsp);
        vector_fp netProductionRates (m_nsp);
        vector_fp enthalpy (m_nsp);
        vector_fp dcpdT (m_nsp);
        // getting physical quantities
        m_thermo->specific_heats_ddT(dcpdT.data());
        m_thermo->getPartialMolarCp(specificHeat.data());
        m_thermo->getPartialMolarEnthalpies(enthalpy.data());
        m_kin->getNetProductionRates(netProductionRates.data());
        // getting perturbed changes w.r.t temperature
        double nkcpksum = 0;
        double hkndotksum = 0;
        double ndotkcpksum = 0;
        double nkdcpdTksum = 0;
        double hkdndotkdTsum = 0;
        double Ntotal = accumulate(y + m_sidx, y + m_nsp, 0.0);
        double NtotalCp = Ntotal * m_thermo->cp_mole();
        double cpksum = accumulate(specificHeat.begin(), specificHeat.end(), 0.0);
        // scale net production rates by inverse of volume to get molar rate
        scale(netProductionRates.begin(), netProductionRates.end(), netProductionRates.begin(), m_vol);
        // determine a sum in derivative
        for (size_t i = 0; i < m_nsp; i++) {
            nkdcpdTksum += y[m_sidx + i] * dcpdT[i];
            hkdndotkdTsum += enthalpy[i] * dndotdT[i];
            nkcpksum += y[m_sidx + i] * specificHeat[i];
            hkndotksum += enthalpy[i] * netProductionRates[i];
            ndotkcpksum += netProductionRates[i] * specificHeat[i];
        }
        // dTdot/dT
        double dTdotdT = -(ndotkcpksum + hkdndotkdTsum) * nkcpksum;
        dTdotdT += hkndotksum * (Ntotal / temperature() * cpksum + nkdcpdTksum);
        dTdotdT /= nkcpksum * nkcpksum;
        m_jac_trips.emplace_back(0, 0, dTdotdT);
        // determine derivatives
        for (size_t j = 0; j < m_nsp; j++) { // spans columns
            double hkdnkdnjSum = 0;
            for (size_t k = 0; k < m_nsp; k++) { // spans rows
                hkdnkdnjSum += enthalpy[k] * speciesDervs.coeff(k, j);
            }
            // set appropriate column of preconditioner
            m_jac_trips.emplace_back(0, j + m_sidx, (-hkdnkdnjSum * NtotalCp +
            specificHeat[j] * hkndotksum) / (NtotalCp * NtotalCp));
        }
    }
}

size_t IdealGasConstPressureMoleReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + m_sidx;
    } else if (nm == "temperature") {
        return 0;
    } else {
        return npos;
    }
}

std::string IdealGasConstPressureMoleReactor::componentName(size_t k) {
    if (k == 0) {
        return "temperature";
    } else if (k >= m_sidx && k < neq()) {
        k -= m_sidx;
        if (k < m_thermo->nSpecies()) {
            return m_thermo->speciesName(k);
        } else {
            k -= m_thermo->nSpecies();
        }
        for (auto& S : m_surfaces) {
            ThermoPhase* th = S->thermo();
            if (k < th->nSpecies()) {
                return th->speciesName(k);
            } else {
                k -= th->nSpecies();
            }
        }
    }
    throw CanteraError("IdealGasConstPressureMoleReactor::componentName",
                       "Index is out of bounds.");
}

}
