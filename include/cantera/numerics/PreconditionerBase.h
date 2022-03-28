/**
 *  @file PreconditionerBase.h Declarations for the class
 *   PreconditionerBase which is a virtual base class for
 *   preconditioning systems.
 */

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#ifndef PRECONDITIONERBASE_H
#define PRECONDITIONERBASE_H

#include "cantera/base/ctexceptions.h"
#include "cantera/numerics/Integrator.h"

namespace Cantera
{

//! Forward declarations
class ReactorNet;
class Reactor;
class MoleReactor;
class IdealGasMoleReactor;
class IdealGasConstPressureMoleReactor;

//! Flag to indicate preconditioner is not set
const int PRECONDITIONER_NOT_SET = 0;

//! PreconditionerBase serves as an abstract type to extend different
//! preconditioners from.
class PreconditionerBase
{
public:
    PreconditionerBase(/* args */){}
    ~PreconditionerBase(){}

    //! Use this function to return zero for preconditioner not set
    virtual size_t getPreconditionerMethod(){return PRECONDITIONER_NOT_SET;};

    //! Use this function to the "type" for CVODES
    virtual PreconditionerType getPreconditionerType(){return NO_PRECONDITION;};

    //! Use this function to solve a linear system Ax=b where A is the
    //! preconditioner contained in this matrix
    //! @param[in] state_len length of the rhs and output vectors
    //! @param[in] rhs_vector right hand side vector used in linear
    //! system
    //! @param[out] output guess vector used by GMRES
    virtual void solve(const size_t state_len, double *rhs_vector, double* output){
        throw NotImplementedError("PreconditionerBase::solve");
    };

    //! Use this function to perform preconditioner specific post-reactor
    //! setup operations such as factorize.
    virtual void setup()
    {
        throw NotImplementedError("PreconditionerBase::setup");
    };

    //! This function performs preconditioner specific post-reactor
    //! setup operations such as factorize.
    virtual void reset()
    {
        throw NotImplementedError("PreconditionerBase::reset");
    };

    //! This function is called during setup for any processes that
    //! need to be completed prior to setup functions
    //! @param network A pointer to the reactor net object
    //! associated with the integration
    virtual void initialize(ReactorNet& network){
        throw NotImplementedError("PreconditionerBase::initialize");
    };

    //! Use this function to print preconditioner contents
    virtual void printPreconditioner(){
        throw NotImplementedError("PreconditionerBase::printPreconditioner");
    };

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with Reactor
    virtual void acceptReactor(Reactor& reactor, double t, double* LHS, double* RHS){
        throw NotImplementedError("PreconditionerBase::acceptReactor");
    };

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with MoleReactor
    virtual void acceptReactor(MoleReactor& reactor, double t, double* LHS, double* RHS){
        throw NotImplementedError("PreconditionerBase::acceptReactor");
    };

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with IdealGasMoleReactor
    virtual void acceptReactor(IdealGasMoleReactor& reactor, double t, double* LHS, double* RHS){
        throw NotImplementedError("PreconditionerBase::acceptReactor");
    };

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with IdealGasConstPressureMoleReactor
    virtual void acceptReactor(IdealGasConstPressureMoleReactor& reactor, double t, double* LHS, double* RHS){
        throw NotImplementedError("PreconditionerBase::acceptReactor");
    };

    //! Use this function to set gamma
    //! @param gamma used in M = I - gamma*J
    virtual void setGamma(double gamma){
        m_gamma = gamma;
    };

    //! Use this function to get gamma
    virtual double getGamma(){
        return m_gamma;
    };

    //! Use this function to set the absolute tolerance in the
    //! solver outside of the network initialization
    //! @param atol the specified tolerance
    void setAbsoluteTolerance(double atol){m_atol = atol;};

    //! Use this function to set dimensions of the preconditioner
    //! @param dims A pointer to a vector of the dimensions
    void setDimensions(std::vector<size_t> *dims){
        this->m_dimensions.clear();
        for (auto it = dims->begin(); it != dims->end(); ++it)
        {
            this->m_dimensions.push_back(*it);
        }
    };

    //! Use this function to return pointer to dimensions
    std::vector<size_t>* getDimensions(){
        return &(this->m_dimensions);
    };

    //! Use this function to get the absolute tolerance from the
    //! preconditioner
    double getAbsoluteTolerance(){return m_atol;};

    //! A counter variable used by network to identify the current reactor index during setup
    size_t m_rctr = 0;

protected:
    //! a size_t vector of dimensions
    std::vector<size_t> m_dimensions;

    //! gamma value used in M = I - gamma*J
    double m_gamma = 1.0;

    //! bool saying whether or not the preconditioner is initialized
    bool m_init = false;

    //! Absolute tolerance of the ODE solver
    double m_atol = 0;
};

}
#endif
