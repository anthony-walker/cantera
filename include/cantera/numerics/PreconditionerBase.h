/**
 *  @file PreconditionerBase.h Declarations for the class
 *   PreconditionerBase which is a virtual base class for
 *   preconditioning systems.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef PRECONDITIONERBASE_H
#define PRECONDITIONERBASE_H

#include "cantera/base/ctexceptions.h"
#include "cantera/base/AnyMap.h"
#include "float.h"

namespace Cantera
{

/**
 * Specifies the preconditioner type used in the integration if any. Not all methods are
 * supported by all integrators.
 */
enum class PreconditionerType {
    NO_PRECONDITION, //! No preconditioning
    LEFT_PRECONDITION, //! Left side preconditioning
    RIGHT_PRECONDITION, //! Right side preconditioning
    BOTH_PRECONDITION //! Left and right side preconditioning
};

//! PreconditionerBase serves as an abstract type to extend different preconditioners
class PreconditionerBase
{
public:
    PreconditionerBase() {}

    virtual const double& operator() (size_t row, size_t col) {
        throw NotImplementedError("PreconditionerBase::operator()");
    }

    virtual void operator() (size_t row, size_t col, double value) {
        throw NotImplementedError("PreconditionerBase::operator()");
    }

    //! Get preconditioner type for CVODES
    virtual PreconditionerType preconditionerType() { return PreconditionerType::NO_PRECONDITION; };

    //! Solve a linear system Ax=b where A is the preconditioner
    //! @param[in] stateSize length of the rhs and output vectors
    //! @param[in] rhs_vector right hand side vector used in linear system
    //! @param[out] output output vector for solution
    virtual void solve(const size_t stateSize, double *rhs_vector, double* output) {
        throw NotImplementedError("PreconditionerBase::solve");
    };

    //! Perform preconditioner specific post-reactor setup operations such as factorize.
    virtual void setup() {
        throw NotImplementedError("PreconditionerBase::setup");
    };

    //! Function to reset preconditioner parameters as needed
    virtual void reset() {
        throw NotImplementedError("PreconditionerBase::reset");
    };

    //! This function is called during setup for any processes that need to be completed
    //! prior to setup functions
    //! @param network A pointer to the reactor net object associated with the
    //! integration
    virtual void initialize(size_t networkSize) {
        throw NotImplementedError("PreconditionerBase::initialize");
    };

    //! Print preconditioner contents
    virtual void printPreconditioner() {
        throw NotImplementedError("PreconditionerBase::printPreconditioner");
    };


    //! Set gamma used in preconditioning
    //! @param gamma used in M = I - gamma*J
    virtual void setGamma(double gamma) {
        m_gamma = gamma;
    };

    //! Get gamma used in preconditioning
    virtual double gamma() {
        return m_gamma;
    };

    //! Set the absolute tolerance in the solver outside of the network initialization
    //! @param atol the specified tolerance
    virtual void setAbsoluteTolerance(double atol) {
        m_atol = atol;
    }

    //! Set dimensions of the preconditioner
    //! @param dims A pointer to a vector of the dimensions
    virtual void setDimensions(std::vector<size_t> *dims) {
        this->m_dimensions.clear();
        for (auto it = dims->begin(); it != dims->end(); ++it) {
            this->m_dimensions.push_back(*it);
        }
    };

    //! Return pointer to dimensions
    virtual std::vector<size_t>* dimensions() { return &(this->m_dimensions); }

    //! Get the absolute tolerance from the preconditioner
    virtual double absoluteTolerance() { return m_atol; };

    //! Set derivative settings
    virtual void setPreconSettings(AnyMap& settings) {
        m_settings = settings;
    }

    //! Get derivative settings from kinetics objects
    virtual void preconSettings(AnyMap& settings) {
        settings = m_settings;
    }

    //! Set the perturbation constant used in finite difference calculations
    //! @param perturb the new pertubation constant
    virtual void setPerturbation(double perturb) {
        m_perturb = perturb;
    }

    //! Get the pertubation constant
    virtual double perturbation() { return m_perturb; }

    //! Update counter variable for use with identifying current reactor
    virtual size_t& counter() { return m_rctr; }

protected:

    //! A counter variable used by network to identify the current reactor during setup
    size_t m_rctr = 0;

    //! a size_t vector of dimensions
    std::vector<size_t> m_dimensions;

    //! gamma value used in M = I - gamma*J
    double m_gamma = 1.0;

    //! bool saying whether or not the preconditioner is initialized
    bool m_init = false;

    //! Absolute tolerance of the ODE solver
    double m_atol = 0;

    //! Perturbation for numerical calculations
    double m_perturb = std::sqrt(DBL_EPSILON);

    //! Derivative settings to be passed to reactors
    AnyMap m_settings;
};

}
#endif
