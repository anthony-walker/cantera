//! @file IdealGasConstPressureMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#ifndef CT_IDEALGASCONSTPRESSMOLE_REACTOR_H
#define CT_IDEALGASCONSTPRESSMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/*!
 * IdealGasConstPressureMoleReactor is a class for ideal gas
 * constant-pressure reactors which use a state of moles. The reactor
 * may have an arbitrary number of inlets and outlets, each of which may
 * be connected to a "flow device" such as a mass flow controller, a
 * pressure regulator, etc. Additional reactors may be connected to the
 * other end of the flow device, allowing construction of arbitrary
 * reactor networks.
 */
class IdealGasConstPressureMoleReactor : public MoleReactor
{
public:
    IdealGasConstPressureMoleReactor(){};

    //! Deprecated function for returning type as a string
    virtual std::string typeStr() const {
        warn_deprecated("IdealGasConstPressureMoleReactor::typeStr",
                        "To be removed after Cantera 2.6. Use type() instead.");
        return "IdealGasConstPressureMoleReactor";
    };

    //! Use this function to return a string of reactor type
    virtual std::string type() const {
        return "IdealGasConstPressureMoleReactor";
    };

    //! Return the index in the solution vector for this MoleReactor of
    //! the component named *nm*. Possible values for *nm* are "mass",
    //! "volume", "int_energy", the name of a homogeneous phase species,
    //! or the name of a surface species.
    virtual size_t componentIndex(const std::string& nm) const;

    //! Return the name of the solution component with index *i*.
    //! @see componentIndex()
    virtual std::string componentName(size_t k);

    //! Use this function to set the thermo manager of this reactor
    virtual void setThermoMgr(ThermoPhase& thermo);

    //! Use this function to get the state in moles
    virtual void getState(double* N);

    //! Use this function to initialize the reactor
    virtual void initialize(double t0 = 0.0);

    //! Right hand side function used to integrate by CVODES
    //! @param t current time of the simulation
    //! @param LHS state vector in moles
    //! @param RHS derivative vector in moles per second
    virtual void eval(double t, double* LHS, double* RHS);

    //! Use to update state vector N
    virtual void updateState(double* N);

    //! This function is the next level of preconditioner setup used in
    //! the visitor design pattern. This is necessary for determining
    //! specific types of both the reactor and preconditioner object
    //! @param preconditioner the preconditioner being used by cvodes
    //! @param t current time of the simulation
    //! @param LHS state vector in moles
    //! @param RHS derivative vector in moles per second
    virtual void reactorPreconditionerSetup(AdaptivePreconditioner& preconditioner, double t, double* LHS, double* RHS);

protected:
    vector_fp m_hk; //!< Species molar enthalpies

    //! const value for the species start index
    const int m_sidx = 1;
};

}

#endif
