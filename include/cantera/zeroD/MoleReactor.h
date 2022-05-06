//! @file MoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MOLEREACTOR_H
#define CT_MOLEREACTOR_H

#include "Reactor.h"

namespace Cantera
{

/*!
 * MoleReactor is a class base class similar which use a state of moles.
 * The reactor may have an arbitrary number of inlets and outlets, each
 * of which may be connected to a "flow device" such as a mass flow
 * controller, a pressure regulator, etc. Additional reactors may be
 * connected to the other end of the flow device, allowing construction
 * of arbitrary reactor networks.
 */
class MoleReactor : public Reactor
{
public:
    MoleReactor() {}

    //! Return a string of reactor type
    virtual std::string type() const {
        return "MoleReactor";
    }

    //! Run reactor specific initialization
    virtual void initialize(double t0 = 0.0);

    //! Right hand side function used to integrate by CVODES
    //! @param t current time of the simulation
    //! @param LHS state vector in moles
    //! @param RHS derivative vector in moles per second
    virtual void eval(double t, double* LHS, double* RHS) {
        throw NotImplementedError("MoleReactor::eval()");
    }

    //! Return the index in the solution vector the component named
    //! *nm*. Possible values for *nm* are "int_energy", "volume", the
    //! name of a homogeneous phase species, or the name of a surface
    //! species.
    virtual size_t componentIndex(const std::string& nm) const;

    //! Return the name of the solution component with index *k*.
    //! @see componentIndex()
    virtual std::string componentName(size_t k);

protected:
    //! Evaluate terms related to surface reactions.
    //! @param[out] LHS   Multiplicative factor on the left hand side of ODE for surface
    //!                   species moles
    //! @param[out] RHS   Right hand side of ODE for surface species moles
    //! @param[out] sdot  array of production rates of bulk phase species on surfaces
    //!                   [kmol/s]
    virtual void evalSurfaces(double* LHS, double* RHS, double* sdot);

    //! Update the state of SurfPhase objects attached to this
    //! MoleReactor
    virtual void updateSurfaceState(double* y);

    //! Get initial conditions for SurfPhase objects attached to this
    //! MoleReactor
    virtual void getSurfaceInitialConditions(double* y);

    //! const value for the species start index
    const int m_sidx = 2;
};

}

#endif
