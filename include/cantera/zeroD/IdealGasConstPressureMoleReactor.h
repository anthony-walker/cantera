//! @file ConstPressureReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASCONSTPRESSMOLE_REACTOR_H
#define CT_IDEALGASCONSTPRESSMOLE_REACTOR_H

#include "MoleReactor.h"

namespace Cantera
{

/**
 * Class ConstPressureReactor is a class for constant-pressure reactors. The
 * reactor may have an arbitrary number of inlets and outlets, each of which may
 * be connected to a "flow device" such as a mass flow controller, a pressure
 * regulator, etc. Additional reactors may be connected to the other end of the
 * flow device, allowing construction of arbitrary reactor networks.
 */
class IdealGasConstPressureMoleReactor : public MoleReactor
{
public:
    IdealGasConstPressureMoleReactor() {}

    virtual std::string typeStr() const {
        warn_deprecated("IdealGasConstPressureMoleReactor::typeStr",
                        "To be removed after Cantera 2.6. Use type() instead.");
        return "IdealGasConstPressureMoleReactor";
    }
    virtual std::string type() const {
        return "IdealGasConstPressureMoleReactor";
    }
    virtual void setThermoMgr(ThermoPhase& thermo);
    virtual void getState(double* N);
    virtual void initialize(double t0 = 0.0);
    virtual void evalEqs(double t, double* N,
                         double* Ndot, double* params);
    //! Use to update state vector N
    virtual void updateState(double* N);

    virtual void acceptPreconditioner(PreconditionerBase *preconditioner, double t, double* N, double* Ndot, double* params);

protected:
    std::vector<double> m_hk; //!< Species molar enthalpies
};
}

#endif
