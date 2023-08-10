/**
 *  @file MultiPhaseEquilSolver.cpp
 *    Driver routine for the VCSnonideal equilibrium solver package
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/equil/MultiPhaseEquilSolver.h"

namespace Cantera
{

MultiPhaseEquilSolver::MultiPhaseEquilSolver(MultiPhase* mix) :
    m_mix(mix)
{
    // get complete formula matrix matrix (elements x species)
    // data for sparse matrix
    m_data.reserve(mix->nSpecies());
    m_columns.reserve(mix->nSpecies());
    m_row_starts.reserve(mix->nElements()+1);
    // loop over all species - columns
    for (size_t j = 0; j < mix->nElements(); j++) {
        m_row_starts.push_back(m_data.size());
        for (size_t k = 0; k < mix->nSpecies(); k++) {
        // loop over all elements - rows
            double atoms = mix->nAtoms(k, j);
            if (atoms > 0) {
                m_data.push_back(atoms);
                m_columns.push_back(k);
            }
        }
    }
    m_row_starts.push_back(m_data.size());
    // create the necessary nonlinear algebraic solver in Sundials
    m_sp_formula_mat = SUNSparseMatrix(mix->nElements(), mix->nSpecies(), m_data.size(), CSR_MAT, m_sundials_ctx.get());
    // assign pointers to  SUNMatrix Content
    auto indexvals = SUNSparseMatrix_IndexValues(m_sp_formula_mat); // column indices
    indexvals = m_columns.data();
    auto colvals = SUNSparseMatrix_Data(m_sp_formula_mat);
    colvals = m_data.data();
    auto rowptrs = SUNSparseMatrix_IndexPointers(m_sp_formula_mat);
    rowptrs = m_row_starts.data();

    // free memory if it is not free
    if (m_kin_mem) {
        KINFree(&m_kin_mem);
    }

    #if CT_SUNDIALS_VERSION < 40
        m_kin_mem = KINCreate(m_method, CV_NEWTON);
    #elif CT_SUNDIALS_VERSION < 60
        m_kin_mem = KINCreate(m_method);
    #else
        m_kin_mem = KINCreate(m_sundials_ctx.get());
    #endif
    // if (!m_kin_mem) {
    //     throw CanteraError("kinsIntegrator::initialize",
    //                        "kinCreate failed.");
    // }
    // m_kin_mem = KINCreate(m_sundials_ctx);
}

int MultiPhaseEquilSolver::equilibrate_TV(int XY, double xtarget, int estimateEquil,
                                        double err, int maxsteps, int loglevel)
{

}

int MultiPhaseEquilSolver::equilibrate_HP(double Htarget, int XY, double Tlow,
    double Thigh, int estimateEquil, double err, int maxsteps,
    int loglevel)
{
     return 0;
}

int MultiPhaseEquilSolver::equilibrate_SP(double Starget, double Tlow, double Thigh,
    int estimateEquil, double err, int maxsteps, int loglevel)
{
    return 0;
}

int MultiPhaseEquilSolver::equilibrate(int XY, int estimateEquil, double err, int maxsteps, int loglevel)
{
    return 0;
}

int MultiPhaseEquilSolver::equilibrate_TP(int estimateEquil, double err,
                                        int maxsteps, int loglevel)
{
    return 0;
}

}
