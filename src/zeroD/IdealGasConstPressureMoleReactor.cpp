//! @file IdealGasConstPressureMoleReactor.cpp A constant pressure
//! zero-dimensional reactor with moles as the state

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#include "cantera/zeroD/IdealGasConstPressureMoleReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{

void IdealGasConstPressureMoleReactor::setThermoMgr(ThermoPhase& thermo)
{
    if (thermo.type() != "IdealGas") {
        throw CanteraError("IdealGasConstPressureMoleReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void IdealGasConstPressureMoleReactor::getState(double* N)
{
    if (m_thermo == 0) {
        throw CanteraError("IdealGasConstPressureMoleReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);
    // get mass for calculations
    m_mass = m_thermo->density() * m_vol;
    // set the second component to the temperature
    N[0] = m_thermo->temperature();
    // Use inverse molecular weights
    const double* Y = m_thermo->massFractions();
    const vector_fp& imw = m_thermo->inverseMolecularWeights();
    double *Ns = N + m_sidx;
    for (size_t i = 0; i < m_nsp; i++)
    {
        Ns[i] = m_mass * imw[i] * Y[i];
    }
    // set the remaining components to the surface species coverages on
    // the walls
    getSurfaceInitialConditions(N + m_nsp + m_sidx);
}

void IdealGasConstPressureMoleReactor::initialize(double t0)
{
    MoleReactor::initialize(t0);
    m_nv -= 1; // const pressure system loses 1 more variable from MoleReactor
    m_hk.resize(m_nsp, 0.0);
}

void IdealGasConstPressureMoleReactor::updateState(double* N)
{
    // the components of N are: [0] the temperature, [1...K+1) are the
    // moles of each species, and [K+1...] are the coverages of surface
    // species on each wall.
    m_thermo->setMolesNoNorm(N + m_sidx);
    m_thermo->setState_TP(N[0], m_pressure);
    // get mass
    vector_fp mass(m_nv-m_sidx);
    const vector_fp& mw = m_thermo->molecularWeights();
    copy(N + m_sidx, N + m_nv, mass.begin());
    transform(mass.begin(), mass.end(), mw.begin(),
              mass.begin(), multiplies<double>());
    m_mass = accumulate(mass.begin(), mass.end(), 0.0);
    m_vol = m_mass / m_thermo->density();
    updateSurfaceState(N + m_nsp + m_sidx);
    updateConnected(false);
}

void IdealGasConstPressureMoleReactor::eval(double time, double* LHS,
                                   double* RHS)
{
    double mcpdTdt = 0.0; // m * c_p * dT/dt
    double* dNdt = RHS + m_sidx; // kmol per s

    evalWalls(time);

    m_thermo->restoreState(m_state);

    m_thermo->getPartialMolarEnthalpies(&m_hk[0]);
    const vector_fp& imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // evaluate surfaces
    evalSurfaces(LHS + m_nsp + m_sidx, RHS + m_nsp + m_sidx, m_sdot.data());

    // external heat transfer
    mcpdTdt -= m_Q;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcpdTdt -= m_wdot[n] * m_hk[n] * m_vol;
        mcpdTdt -= m_sdot[n] * m_hk[n];
        // production in gas phase and from surfaces
        dNdt[n] = (m_wdot[n] * m_vol + m_sdot[n]);
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species into system and dilution by other species
            dNdt[n] -= outlet->outletSpeciesMolarFlowRate(n);
        }
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        mcpdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dNdt[n] += inlet->outletSpeciesMolarFlowRate(n);
            mcpdTdt -= m_hk[n] * imw[n] * mdot_spec;
        }
    }

    if (m_energy) {
        RHS[0] = mcpdTdt / (m_mass * m_thermo->cp_mass());
    } else {
        RHS[0] = 0.0;
    }

}

void IdealGasConstPressureMoleReactor::reactorPreconditionerSetup(AdaptivePreconditioner& preconditioner, double t, double* LHS, double* RHS)
{
    // strictly positive composition
    vector_fp LHSCopy(m_nv);
    preconditioner.getStrictlyPositiveComposition(m_nv, LHS, LHSCopy.data());
    updateState(LHSCopy.data());
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
    double scalingFactor = m_vol/accumulate(LHS + m_sidx, LHS + m_nv, 0.0);
    Eigen::SparseMatrix<double> speciesDervs = m_kin->netProductionRates_ddX();
    // sum parts
    speciesDervs = scalingFactor * speciesDervs + rates * volumes;
    // add to preconditioner
    for (int k=0; k<speciesDervs.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(speciesDervs, k); it; ++it)
        {
            preconditioner(it.row() + m_sidx, it.col() + m_sidx, it.value());
        }
    }
    // Temperature Derivatives
    if (m_energy)
    {
        // getting perturbed state for finite difference
        double deltaTemp = LHSCopy[0] * preconditioner.getPerturbationConst();
        // finite difference temperature derivatives
        vector_fp NNext (m_nv);
        vector_fp NdotNext (m_nv);
        vector_fp NCurrent (m_nv);
        vector_fp NdotCurrent (m_nv);
        // copy LHS to current and next
        copy(LHSCopy.begin(), LHSCopy.end(), NCurrent.begin());
        copy(LHSCopy.begin(), LHSCopy.end(), NNext.begin());
        // perturb temperature
        NNext[0] += deltaTemp;
        // getting perturbed state
        updateState(NNext.data());
        eval(t, NNext.data(), NdotNext.data());
        // reset and get original state
        updateState(NCurrent.data());
        eval(t, NCurrent.data(), NdotCurrent.data());
        // d T_dot/dT
        preconditioner(0, 0, (NdotNext[0] - NdotCurrent[0]) / deltaTemp);
        // d omega_dot_j/dT
        for (size_t j = m_sidx; j < m_nv; j++)
        {
            preconditioner(j, 0, (NdotNext[j] - NdotCurrent[j]) / deltaTemp);
        }
        // d T_dot/dnj
        vector_fp specificHeat (m_nsp);
        vector_fp netProductionRates (m_nsp);
        vector_fp enthalpy (m_nsp);
        // getting physical quantities
        m_thermo->getPartialMolarCp(specificHeat.data());
        m_thermo->getPartialMolarEnthalpies(enthalpy.data());
        m_kin->getNetProductionRates(netProductionRates.data());
        // getting perturbed changes w.r.t temperature
        double hkndotksum = 0;
        double inverseVolume = 1/volume();
        double NtotalCp = accumulate(LHSCopy.begin() + m_sidx, LHSCopy.end(), 0.0) * m_thermo->cp_mole();
        // scale net production rates by inverse of volume to get molar rate
        scale(netProductionRates.begin(), netProductionRates.end(), netProductionRates.begin(), inverseVolume);
        // determine a sum in derivative
        for (size_t i = 0; i < m_nsp; i++)
        {
            hkndotksum += enthalpy[i] * netProductionRates[i];
        }
        // determine derivatives
        for (size_t j = 0; j < m_nsp; j++) // spans columns
        {
            double hkdnkdnjSum = 0;
            for (size_t k = 0; k < m_nsp; k++) // spans rows
            {
                hkdnkdnjSum += enthalpy[k] * speciesDervs.coeff(k, j);
            }
            // set appropriate column of preconditioner
            preconditioner(0, j + m_sidx, (-hkdnkdnjSum * NtotalCp + specificHeat[j] * hkndotksum) / (NtotalCp * NtotalCp));
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
    } else if (k >= 1 && k < neq()) {
        k -= 1;
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
