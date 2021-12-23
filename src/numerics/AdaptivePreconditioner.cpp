//! @file AdaptivePreconditioner.cpp

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/zeroD/MoleReactor.h"
#include "cantera/zeroD/IdealGasConstPressureMoleReactor.h"
#include "cantera/zeroD/IdealGasMoleReactor.h"
#include "cantera/base/global.h"
#include <iostream>

namespace Cantera
{

bool AdaptivePreconditioner::operator== (const AdaptivePreconditioner &externalPrecon)
{
    return m_precon_matrix.isApprox(externalPrecon.m_precon_matrix);
}

void AdaptivePreconditioner::operator= (const AdaptivePreconditioner &externalPrecon)
{
    // copy all variables
    m_jac_trips = externalPrecon.m_jac_trips;
    m_identity = externalPrecon.m_identity;
    m_precon_matrix = externalPrecon.m_precon_matrix;
    m_threshold = externalPrecon.m_threshold;
    m_perturb = externalPrecon.m_perturb;
    m_dimensions = externalPrecon.m_dimensions;
    m_atol = externalPrecon.m_atol;
    m_gamma = externalPrecon.m_gamma;
    m_init = externalPrecon.m_init;
    m_rctr = 0;
    return;
}

void AdaptivePreconditioner::operator() (size_t row, size_t col, double value)
{
    m_jac_trips.push_back(Eigen::Triplet<double>(row + m_rctr, col + m_rctr, value));
}

void AdaptivePreconditioner::initialize(ReactorNet& network)
{
    //! don't use legacy rate constants
    use_legacy_rate_constants(false);
    // reset arrays in case of re-initialization
    m_dimensions.clear();
    m_jac_trips.clear();
    // set dimensions of preconditioner from network
    size_t totalColLen = network.neq();
    m_dimensions.push_back(totalColLen);
    m_dimensions.push_back(totalColLen);
    // Derivative settings
    AnyMap m_settings;
    m_settings["skip-third-bodies"] = true;
    m_settings["skip-falloff"] = true;
    // Loop through reactors
    for (size_t i = 0; i < network.nreactors(); i++)
    {
        // apply settings to each reactor
        Reactor& currReactor = network.reactor(i);
        currReactor.setKineticsDerivativeSettings(m_settings);
    }
    // reserve maximum space for vectors making up SparseMatrix
    m_jac_trips.reserve(totalColLen * totalColLen);
    // reserve space for preconditioner
    m_precon_matrix.resize(m_dimensions[0], m_dimensions[1]);
    // creating sparse identity matrix
    m_identity.resize(m_dimensions[0], m_dimensions[1]);
    m_identity.setIdentity();
    m_identity.makeCompressed();
    // setting ILUT parameters
    setDropTolILUT();
    setFillFactorILUT();
    // update initialized status
    m_init = true;
}

void AdaptivePreconditioner::acceptReactor(MoleReactor& reactor, double t, double* LHS, double* RHS)
{
    reactor.reactorPreconditionerSetup(*this, t, LHS, RHS);
}

void AdaptivePreconditioner::acceptReactor(IdealGasMoleReactor& reactor, double t, double* LHS, double* RHS)
{
     reactor.reactorPreconditionerSetup(*this, t, LHS, RHS);
}

void AdaptivePreconditioner::acceptReactor(IdealGasConstPressureMoleReactor& reactor, double t, double* LHS, double* RHS)
{
     reactor.reactorPreconditionerSetup(*this, t, LHS, RHS);
}

void AdaptivePreconditioner::setup()
{
    // make into preconditioner as P = (I - gamma * J_bar)
    transformJacobianToPreconditioner();
    // compressing sparse matrix structure
    m_precon_matrix.makeCompressed();
    // analyze and factorize
    m_solver.compute(m_precon_matrix);
    // check for errors
    preconditionerErrorCheck();
}

void AdaptivePreconditioner::transformJacobianToPreconditioner()
{
    // set precon to jacobian
    m_precon_matrix.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
    // convert to preconditioner
    m_precon_matrix = m_identity - m_gamma * m_precon_matrix;
    // prune by threshold if desired
    if (m_prune_precon)
    {
        prunePreconditioner();
    }
}

void AdaptivePreconditioner::prunePreconditioner()
{
    for (int k=0; k<m_precon_matrix.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_precon_matrix, k); it; ++it)
        {
            if (std::abs(it.value()) < m_threshold && it.row() != it.col())
            {
                m_precon_matrix.coeffRef(it.row(), it.col()) = 0;
            }
        }
    }
}

void AdaptivePreconditioner::preconditionerErrorCheck()
{
    if (m_solver.info() != Eigen::Success)
    {
        throw CanteraError("AdaptivePreconditioner::solve",
                           "error code: {}", m_solver.info());
    }
}

void AdaptivePreconditioner::solve(const size_t state_len, double *rhs_vector, double* output)
{
    // creating vectors in the form of Ax=b
    Eigen::Map<Eigen::VectorXd> bVector(rhs_vector, state_len);
    Eigen::Map<Eigen::VectorXd> xVector(output, state_len);
    // solve for xVector
    xVector = m_solver.solve(bVector);
    preconditionerErrorCheck();
}

}
