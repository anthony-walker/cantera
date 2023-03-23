//! @file AdaptivePreconditioner.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/global.h"
#include "cantera/numerics/StateDiagonalPreconditioner.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/ReactorNet.h"

namespace Cantera
{

void StateDiagonalPreconditioner::initialize(size_t networkSize)
{
    AdaptivePreconditioner::initialize(networkSize);
    addReactor(&(m_rnet->reactor(0)));
}

void StateDiagonalPreconditioner::setup()
{
    // make jacobian from simple jacobian
    for (size_t i = 0; i < m_reactors.size(); i++) {
        Eigen::SparseMatrix<double> rJac;
        if (m_reactors[i]->isMoleReactor()) {
            rJac = m_reactors[i]->jacobian();

        } else {
            rJac = m_reactors[i]->finiteDifferenceJacobian(true);
        }
        // add to preconditioner
        for (int k=0; k<rJac.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(rJac, k); it; ++it) {
                if (it.row() == 0 || it.col() == 0 || it.row() == it.col()) {
                    setValue(it.row(), it.col(), it.value());
                }
            }
        }
    }
    // make into preconditioner as P = (I - gamma * J_bar)
    updatePreconditioner();
    // compressing sparse matrix structure
    m_precon_matrix.makeCompressed();
    // analyze and factorize
    m_solver.compute(m_precon_matrix);
    // check for errors
    if (m_solver.info() != Eigen::Success) {
        throw CanteraError("StateDiagonalPreconditioner::setup",
                           "error code: {}", static_cast<int>(m_solver.info()));
    }
}

}
