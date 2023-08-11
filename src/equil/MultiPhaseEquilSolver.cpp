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

extern "C" {

static int kin_rhs(N_Vector x, N_Vector fx, void* user_data)
{
    MultiPhaseEquilSolver* solver = (MultiPhaseEquilSolver*) user_data;
    solver->equilibrate(ConstantState::TP, NV_DATA_S(fx), NV_DATA_S(x), user_data);
    return 0;
}

}

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

    #if CT_SUNDIALS_VERSION < 40
        m_sp_formula_mat = SUNSparseMatrix(mix->nElements(), mix->nSpecies(), m_data.size(), CSR_MAT);
    #elif CT_SUNDIALS_VERSION < 60
        m_sp_formula_mat = SUNSparseMatrix(mix->nElements(), mix->nSpecies(), m_data.size(), CSR_MAT);
    #else
        m_sp_formula_mat = SUNSparseMatrix(mix->nElements(), mix->nSpecies(), m_data.size(), CSR_MAT, m_sundials_ctx.get());
    #endif

    // assign pointers to  SUNMatrix Content
    auto indexvals = SUNSparseMatrix_IndexValues(m_sp_formula_mat); // column indices
    indexvals = m_columns.data();
    auto colvals = SUNSparseMatrix_Data(m_sp_formula_mat);
    colvals = m_data.data();
    auto rowptrs = SUNSparseMatrix_IndexPointers(m_sp_formula_mat);
    rowptrs = m_row_starts.data();

    // resize n vector and get state
    m_state = N_VNew_Serial(2 * mix->nSpecies() + mix->nElements());
    m_constraints = N_VNew_Serial(2 * mix->nSpecies() + mix->nElements());
    mix->getMoles(NV_DATA_S(m_state));

    // free memory if it is not free
    if (m_kin_mem) {
        KINFree(&m_kin_mem);
    }
    #if CT_SUNDIALS_VERSION < 40
        // TODO: Update for lower versions
        m_kin_mem = KINCreate();
    #elif CT_SUNDIALS_VERSION < 60
        m_kin_mem = KINCreate();
        if (!m_kin_mem) {
            throw CanteraError("kinsIntegrator::initialize", "kinCreate failed.");
        }
        // May need to set constraints
        // flag = KINSetConstraints(kmem, c);
        KINInit(m_kin_mem, kin_rhs, m_state);
        /* Create sparse SUNMatrix */
        m_linsol_matrix = (void *) SUNDenseMatrix(mix->nSpecies(), mix->nSpecies());
        // if(check_flag(m_linsol_matrix, "SUNSparseMatrix", 0)) return(1);
        /* Create linear solver object */
        #if CT_SUNDIALS_USE_LAPACK
            m_linsol = (void *) SUNLinSol_LapackDense(m_state,
                                                      (SUNMatrix) m_linsol_matrix);
        #else
            m_linsol = (void *) SUNLinSol_Dense(m_state, (SUNMatrix) m_linsol_matrix);
        #endif
        /* Attach dense linear solver */
        KINSetLinearSolver(m_kin_mem, (SUNLinearSolver) m_linsol, (SUNMatrix) m_linsol_matrix);
    #else
        //TODO: Update for higher versions
        m_kin_mem = KINCreate(m_sundials_ctx.get());
    #endif

}

int MultiPhaseEquilSolver::equilibrate(ConstantState cs, double* LHS, double* RHS, void* f_data)
{
    if (cs == ConstantState::TP) {
        equilibrate_TP(LHS, RHS, f_data);
    }

    return 0;
}

int MultiPhaseEquilSolver::equilibrate_TP(double* LHS, double* RHS, void* f_data)
{
    vector<double> mu;
    double* n = RHS;
    double* y = RHS + m_mix->nSpecies();
    double* z = RHS + m_mix->nSpecies() + m_mix->nElements();
    m_mix->getChemPotentials(mu.data());

    return 0;
}

int MultiPhaseEquilSolver::equilibrate_TV(int XY, double xtarget, int estimateEquil,
                                        double err, int maxsteps, int loglevel)
{
    return 0;
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

}
