//! @file IdealGasConstPressureMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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
    IdealGasConstPressureMoleReactor() {}

    //! Return a string of reactor type
    virtual std::string type() const {
        return "IdealGasConstPressureMoleReactor";
    };

    //! Return the index in the solution vector the component named
    //! *nm*. Possible values for *nm* are "temperature", the name of a
    //! homogeneous phase species, or the name of a surface species.
    virtual size_t componentIndex(const std::string& nm) const;

    //! Return the name of the solution component with index *k*.
    //! @see componentIndex()
    virtual std::string componentName(size_t k);

    //! Set the thermo manager of this reactor
    virtual void setThermoMgr(ThermoPhase& thermo);

    //! Get the state in moles
    virtual void getState(double* y);

    //! Initialize the reactor
    virtual void initialize(double t0 = 0.0);

    //! Right hand side function used to integrate by CVODES
    //! @param t current time of the simulation
    //! @param LHS state vector in moles
    //! @param RHS derivative vector in moles per second
    virtual void eval(double t, double* LHS, double* RHS);

    //! Use to update state vector y
    virtual void updateState(double* y);

    //! Method to calculate the reactor specific jacobian
    //! @param t current time of the simulation
    //! @param LHS state vector in moles
    //! @param RHS derivative vector in moles per second
    //! @warning  This method is an experimental part of the %Cantera
    //! API and may be changed or removed without notice.
    virtual Eigen::SparseMatrix<double> jacobian(double t, double* LHS, double* RHS);

protected:
    vector_fp m_hk; //!< Species molar enthalpies

    //! const value for the species start index
    const int m_sidx = 1;
};

}

#endif
