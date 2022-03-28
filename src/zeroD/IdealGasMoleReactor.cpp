//! @file IdealGasMoleReactor.cpp A constant volume zero-dimensional
//! reactor with moles as the state

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#include "cantera/zeroD/IdealGasMoleReactor.h"
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

void IdealGasMoleReactor::setThermoMgr(ThermoPhase& thermo)
{
    if (thermo.type() != "IdealGas") {
        throw CanteraError("IdealGasMoleReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void IdealGasMoleReactor::getState(double* N)
{
    if (m_thermo == 0) {
        throw CanteraError("IdealGasMoleReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // get mass for calculations
    m_mass = m_thermo->density() * m_vol;

    // set the first component to the temperature
    N[0] = m_thermo->temperature();

    // set the second component to the volume
    N[1] = m_vol;

    // use inverse molecular weights
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

void IdealGasMoleReactor::initialize(double t0)
{
   MoleReactor::initialize(t0);
    m_uk.resize(m_nsp, 0.0);
}

void IdealGasMoleReactor::updateState(double* N)
{
    // the components of N are: [0] the temperature, [1] the volume, [2...K+1) are the
    // moles of each species, and [K+1...] are the coverages of surface
    // species on each wall. get mass
    vector_fp mass(m_nv-m_sidx);
    const vector_fp& mw = m_thermo->molecularWeights();
    copy(N + m_sidx, N + m_nv, mass.begin());
    transform(mass.begin(), mass.end(), mw.begin(),
              mass.begin(), multiplies<double>());
    m_mass = accumulate(mass.begin(), mass.end(), 0.0);
    m_vol = N[1];
    // set state
    m_thermo->setMolesNoNorm(N + m_sidx);
    m_thermo->setState_TR(N[0], m_mass / m_vol);
    updateSurfaceState(N + m_nsp + m_sidx);
    updateConnected(true);
}

void IdealGasMoleReactor::eval(double time, double* LHS,
                                   double* RHS)
{
    double mcvdTdt = RHS[0]; // m * c_v * dT/dt
    double* dNdt = RHS + m_sidx; // kmol per s

    evalWalls(time);

    m_thermo->restoreState(m_state);

    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
    const vector_fp& imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // evaluate surfaces
    evalSurfaces(LHS + m_nsp + m_sidx, RHS + m_nsp + m_sidx, m_sdot.data());

    // external heat transfer
    mcvdTdt += - m_pressure * m_vdot - m_Q;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcvdTdt -= m_wdot[n] * m_uk[n] * m_vol;
        mcvdTdt -= m_sdot[n] * m_uk[n];
        // production in gas phase and from surfaces
        dNdt[n] = (m_wdot[n] * m_vol + m_sdot[n]);
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species into system and dilution by other species
            dNdt[n] -= outlet->outletSpeciesMolarFlowRate(n);
        }
        double mdot = outlet->massFlowRate();
        mcvdTdt -= mdot * m_pressure * m_vol / m_mass; // flow work
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        mcvdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dNdt[n] += inlet->outletSpeciesMolarFlowRate(n);
            mcvdTdt -= m_uk[n] * imw[n] * mdot_spec;
        }
    }

    RHS[1] = m_vdot;
    if (m_energy) {
        RHS[0] = mcvdTdt / (m_mass * m_thermo->cv_mass());
    } else {
        RHS[0] = 0.0;
    }

}

void IdealGasMoleReactor::reactorPreconditionerSetup(AdaptivePreconditioner& preconditioner, double t, double* LHS, double* RHS)
{
    // strictly positive composition
    vector_fp LHSCopy(m_nv);
    preconditioner.getStrictlyPositiveComposition(m_nv, LHS, LHSCopy.data());
    updateState(LHSCopy.data());
    // Determine Species Derivatives
    // get ROP derivatives
    double scalingFactor = m_vol/accumulate(LHS + m_sidx, LHS + m_nv, 0.0);
    Eigen::SparseMatrix<double> speciesDervs = scalingFactor * m_kin->netProductionRates_ddX();
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
        // find derivatives d T_dot/dNj
        vector_fp specificHeat (m_nsp);
        vector_fp netProductionRates (m_nsp);
        vector_fp internal_energy (m_nsp);
        vector_fp concentrations (m_nsp);
        // getting species concentrations
        m_thermo->getConcentrations(concentrations.data());
        m_thermo->getPartialMolarCp(specificHeat.data());
        m_thermo->getPartialMolarIntEnergies(internal_energy.data());
        m_kin->getNetProductionRates(netProductionRates.data());
        // convert Cp to Cv for ideal gas as Cp - Cv = R
        for (size_t i = 0; i < specificHeat.size(); i++)
        {
            specificHeat[i] -= GasConstant;
        }
        double uknkSum = 0;
        double NtotalCv = accumulate(LHSCopy.begin() + m_sidx, LHSCopy.end(), 0.0) * m_thermo->cv_mole();
        for (size_t i = 0; i < m_nsp; i++)
        {
            uknkSum += internal_energy[i] * netProductionRates[i];
        }
        for (size_t j = 0; j < m_nsp; j++) // spans columns
        {
            double ukdnkdnjSum = 0;
            for (size_t k = 0; k < m_nsp; k++) // spans rows
            {
                ukdnkdnjSum += internal_energy[k] * speciesDervs.coeff(k, j);
            }
            // set appropriate column of preconditioner
            preconditioner(0, j + m_sidx, (-ukdnkdnjSum * NtotalCv + specificHeat[j] *  uknkSum) / (NtotalCv * NtotalCv));
        }
    }
}

}
