//! @file IdealGasReactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"

using namespace std;

namespace Cantera
{

void IdealGasReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @TODO: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (thermo.type() != "IdealGas") {
        throw CanteraError("IdealGasReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void IdealGasReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("IdealGasReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // set the first component to the total mass
    m_mass = m_thermo->density() * m_vol;
    y[0] = m_mass;

    // set the second component to the total volume
    y[1] = m_vol;

    // Set the third component to the temperature
    y[2] = m_thermo->temperature();

    // set components y+3 ... y+K+2 to the mass fractions of each species
    m_thermo->getMassFractions(y+3);

    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(y + m_nsp + 3);
}

void IdealGasReactor::initialize(doublereal t0)
{
    Reactor::initialize(t0);
    m_uk.resize(m_nsp, 0.0);
}

void IdealGasReactor::updateState(doublereal* y)
{
    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the temperature, [3...K+3] are the mass fractions of each species,
    // and [K+3...] are the coverages of surface species on each wall.
    m_mass = y[0];
    m_vol = y[1];
    m_thermo->setMassFractions_NoNorm(y+3);
    m_thermo->setState_TR(y[2], m_mass / m_vol);
    updateSurfaceState(y + m_nsp + 3);
    updateConnected(true);
}

void IdealGasReactor::evalEqs(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{
    double dmdt = 0.0; // dm/dt (gas phase)
    double* dYdt = ydot + 3;
<<<<<<< HEAD

    evalWalls(time);
    applySensitivity(params);
    m_thermo->restoreState(m_state);
    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
=======
    /*
    Evaluating energy equation where
    flow devices, walls, and etc are also evaluated
    */
    evaluateEnergyEquation(time,y,ydot,params);
    //Getting mass fractions and molecular weights
>>>>>>> c9b17a9b9 (Adding evaluation of energy equation as separate function so as not to repeat work)
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();
    //Evaluation of surfaces
    double mdot_surf = evalSurfaces(time, ydot + m_nsp + 3);
    dmdt += mdot_surf;

    for (size_t n = 0; n < m_nsp; n++) {
        // production in gas phase and from surfaces
        dYdt[n] = (m_wdot[n] * m_vol + m_sdot[n]) * mw[n] / m_mass;
        // dilution by net surface mass flux
        dYdt[n] -= Y[n] * mdot_surf / m_mass;
    }

    // add terms for outlets
<<<<<<< HEAD
    for (auto outlet : m_outlet) {
        double mdot = outlet->massFlowRate();
        dmdt -= mdot; // mass flow out of system
        mcvdTdt -= mdot * m_pressure * m_vol / m_mass; // flow work
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        dmdt += mdot; // mass flow into system
        mcvdTdt += inlet->enthalpy_mass() * mdot;
=======
    for (size_t i = 0; i < m_outlet.size(); i++) {
        // double mdot_out = m_outlet[i]->massFlowRate(time);
        dmdt -= m_mdot_out[i]; // mass flow out of system
    }

    // add terms for inlets
    for (size_t i = 0; i < m_inlet.size(); i++) {
        // double mdot_in = m_inlet[i]->massFlowRate(time);
        dmdt += m_mdot_in[i]; // mass flow into system
>>>>>>> c9b17a9b9 (Adding evaluation of energy equation as separate function so as not to repeat work)
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
<<<<<<< HEAD
            dYdt[n] += (mdot_spec - mdot * Y[n]) / m_mass;

            // In combination with h_in*mdot_in, flow work plus thermal
            // energy carried with the species
            mcvdTdt -= m_uk[n] / mw[n] * mdot_spec;
=======
            dYdt[n] += (mdot_spec - m_mdot_in[i] * Y[n]) / m_mass;
>>>>>>> c9b17a9b9 (Adding evaluation of energy equation as separate function so as not to repeat work)
        }
    }
    ydot[0] = dmdt;
    ydot[1] = m_vdot;
    ydot[2] = m_dEdt / (m_mass * m_thermo->cv_mass()); //m_dEdt set to zero in evaluateEnergyEquation if not m_energy // m * c_v * dT/dt
    resetSensitivity(params);
}

void IdealGasReactor::evaluateEnergyEquation(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{ 
    evalFlowDevices(time);
    evalWalls(time);
    applySensitivity(params);
    m_thermo->restoreState(m_state);
    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
    const vector_fp& mw = m_thermo->molecularWeights();
    if (m_chem) 
    {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }
    if (m_energy)
    {
        // compression work and external heat transfer
        m_dEdt += - m_pressure * m_vdot - m_Q;

        for (size_t n = 0; n < m_nsp; n++) {
            // heat release from gas phase and surface reactions
            m_dEdt -= m_wdot[n] * m_uk[n] * m_vol;
            m_dEdt -= m_sdot[n] * m_uk[n];
        }

        // add terms for outlets
        for (size_t i = 0; i < m_outlet.size(); i++) 
        {
            m_dEdt -= m_mdot_out[i] * m_pressure * m_vol / m_mass; // flow work
        }

        // add terms for inlets
        for (size_t i = 0; i < m_inlet.size(); i++) 
        {
            m_dEdt += m_inlet[i]->enthalpy_mass() * m_mdot_in[i];
            for (size_t n = 0; n < m_nsp; n++) {
                double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
                // In combination with h_in*mdot_in, flow work plus thermal
                // energy carried with the species
                m_dEdt -= m_uk[n] / mw[n] * mdot_spec;
            }
        }
    }
    else
    {
        this->m_dEdt = 0.0; // m * c_v * dT/dt
    }
}


void IdealGasReactor::reactorPrecSetup(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params)
{   
    //Setting up preconditioner
    this->m_preconditioner.setDimensions(this->m_nv,this->m_nv); //setting number of dimensions for preconditioner
    //Defining index variables
    size_t speciesStart = 3; //starting index for species
    //Getting derivative of temp w.r.t time - a lot of extra evaluation here, potentially add boolean to prevent recalculation in RHS function if possible
    this->evaluateEnergyEquation(t,y,ydot,params); //Solving energy equation for preconditioner
    double dTdt =this->m_dEdt/(this->m_mass*(this->m_thermo)->cv_mass()); //adjusting energy to get dTdt - K/s
    //Filling preconditioner based on type of preconditioner
    switch (this->m_preconditioner_type)
    {
    case PRECONDITIONER_NOT_SET:
        throw CanteraError("Reactor::reactorPrecSetup", "preconditioner type not set");
        break;
    case ADAPTIVE_MECHANISM_PRECONDITIONER:
        std::cout<<"IdealGasReactor"<<std::endl;
        // Cantera::AMP::printReactorComponents(this);
        //Species derivatives
        Cantera::AMP::SpeciesSpeciesDerivatives<SundialsSparseMatrix>(&(this->m_preconditioner),this,speciesStart);
        Cantera::AMP::TemperatureDerivatives<SundialsSparseMatrix>(&(this->m_preconditioner),this,ydot,dTdt,2,speciesStart); //Temperature is index location 2
        break;
    default:
        throw CanteraError("Reactor::reactorPrecSetup", "unknown preconditioner type");
        break;
    }
    resetSensitivity(params); //Not quite sure if this is needed
}

void IdealGasReactor::reactorPrecSolve(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params)
{

}

size_t IdealGasReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + 3;
    } else if (nm == "mass") {
        return 0;
    } else if (nm == "volume") {
        return 1;
    } else if (nm == "temperature") {
        return 2;
    } else {
        return npos;
    }
}

std::string IdealGasReactor::componentName(size_t k) {
    if (k == 2) {
        return "temperature";
    } else {
        return Reactor::componentName(k);
    }
}


}
