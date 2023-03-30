//! @file AdaptivePreconditioner.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/global.h"
#include "cantera/numerics/SubmodelPreconditioner.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/ReactorNet.h"

namespace Cantera
{

void SubmodelPreconditioner::initialize(size_t networkSize)
{
    AdaptivePreconditioner::initialize(networkSize);
    // create map
    for (size_t i = 0; i < networkSize; i++) {
        m_submodel_map[i] = true;
    }
    size_t ctr = 0;
    for (auto subr : m_reactors) {
        for (size_t i = 0; i < subr->neq(); i++) {
            string comp = subr->componentName(i);
            size_t netIdx = m_rnet->reactor(ctr).componentIndex(comp);
            m_submodel_map[netIdx] = false;
        }
        ctr += 1;
    }
}

void SubmodelPreconditioner::prunePreconditioner()
{
    for (int k=0; k<m_precon_matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_precon_matrix, k); it;
            ++it) {
            bool subrow = m_submodel_map[it.row()];
            bool subcol = m_submodel_map[it.col()];
            if (std::abs(it.value()) < m_threshold && it.row() != it.col() && subrow && subcol) {
                it.valueRef() = 0;
            }
        }
    }
}

// void SubmodelPreconditioner::stateAdjustment(vector_fp& state) {
//     // copy initial reactors states over
//     for (auto subr : m_reactors) {
//         vector_fp sub_state(subr->neq());
//         for (size_t i = 0; i < subr->neq(); i++) {
//             string comp = subr->componentName(i);
//             sub_state[i] = state[m_submodel_map[i]];
//         }
//         subr->updateState(sub_state.data());
//     }
// }

// void SubmodelPreconditioner::setup()
// {
//     // make jacobian from simple jacobian
//     for (size_t i = 0; i < m_reactors.size(); i++) {
//         Eigen::SparseMatrix<double> rJac;
//         if (m_reactors[i]->isMoleReactor()) {
//             rJac = m_reactors[i]->jacobian();
//         } else {
//             rJac = m_reactors[i]->finiteDifferenceJacobian();
//         }
//         // add to preconditioner
//         for (int k=0; k<rJac.outerSize(); ++k) {
//             for (Eigen::SparseMatrix<double>::InnerIterator it(rJac, k); it; ++it) {
//                 size_t row = m_submodel_map[it.row()];
//                 size_t col = m_submodel_map[it.col()];
//                 setValue(row, col, it.value());
//             }
//         }
//     }
//     // make into preconditioner as P = (I - gamma * J_bar)
//     updatePreconditioner();
//     // compressing sparse matrix structure
//     m_precon_matrix.makeCompressed();
//     // analyze and factorize
//     m_solver.compute(m_precon_matrix);
//     // check for errors
//     if (m_solver.info() != Eigen::Success) {
//         throw CanteraError("SubmodelPreconditioner::setup",
//                            "error code: {}", static_cast<int>(m_solver.info()));
//     }
// }

}
