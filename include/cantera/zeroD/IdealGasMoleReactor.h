//! @file IdealGasMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASMOLE_REACTOR_H
#define CT_IDEALGASMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/*!
 * IdealGasMoleReactor is a class for ideal gas
 * constant-volume reactors which use a state of moles. The reactor
 * may have an arbitrary number of inlets and outlets, each of which may
 * be connected to a "flow device" such as a mass flow controller, a
 * pressure regulator, etc. Additional reactors may be connected to the
 * other end of the flow device, allowing construction of arbitrary
 * reactor networks.
 */
class IdealGasMoleReactor : public MoleReactor
{
public:
    IdealGasMoleReactor() {}

    //! Return a string of reactor type
    virtual std::string type() const {
        return "IdealGasMoleReactor";
    }

    //! Return the index in the solution vector the component named
    //! *nm*. Possible values for *nm* are "temperature", "volume", the
    //! name of a homogeneous phase species, or the name of a surface
    //! species.
    virtual size_t componentIndex(const std::string& nm) const;

    //! Set the thermo manager of this reactor
    virtual void setThermoMgr(ThermoPhase& thermo);

    //! Get the state in moles
    virtual void getState(double* y);

    //! Initialize the reactor
    virtual void initialize(double t0 = 0.0);

    //! Right hand side function used to integrate by CVODES
    //! @param t current time of the simulation
    //! @param[out] LHS pointer to start of vector of left-hand side
    //! coefficients for governing equations, length m_nv, default values 1
    //! @param[out] RHS pointer to start of vector of right-hand side
    //! coefficients for governing equations, length m_nv, default values 0
    virtual void eval(double t, double* LHS, double* RHS);

    //! Use to update state vector y
    virtual void updateState(double* y);

    //! Method to calculate the reactor specific jacobian
    //! @param t current time of the simulation
    //! @param[out] LHS pointer to start of vector of left-hand side
    //! coefficients for governing equations, length m_nv, default values 1
    //! @param[out] RHS pointer to start of vector of right-hand side
    //! coefficients for governing equations, length m_nv, default values 0
    //! @warning  This method is an experimental part of the %Cantera
    //! API and may be changed or removed without notice.
    virtual Eigen::SparseMatrix<double> jacobian(double t, double* LHS, double* RHS);

protected:
    vector_fp m_uk; //!< Species molar internal energies
};

}

#endif
