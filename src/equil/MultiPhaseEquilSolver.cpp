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

N_Vector newNVector(size_t N, Cantera::SundialsContext& context)
{
#if CT_SUNDIALS_VERSION >= 60
    return N_VNew_Serial(static_cast<sd_size_t>(N), context.get());
#else
    return N_VNew_Serial(static_cast<sd_size_t>(N));
#endif
}


MultiPhaseEquilSolver::MultiPhaseEquilSolver(MultiPhase* mix) :
    m_mix(mix)
{
    // Sizes needed
    size_t nsp = m_mix->nSpecies();
    size_t nele = m_mix->nElements();
    m_laplace_size = 2 * nsp + nele;
    // get complete formula matrix matrix (elements x species)
    vector<Eigen::Triplet<double>> trips;
    // reserve some initial space
    trips.reserve(nsp);
    // loop over all elements - rows
    for (size_t j = 0; j < nele; j++) {
        // loop over all species - columns
        for (size_t k = 0; k < nsp; k++) {
            double atoms = mix->nAtoms(k, j);
            if (atoms > 0) {
                trips.emplace_back(j, k, atoms);
            }
        }
    }
    // resize eigen formula matrix
    m_formula_mat.resize(nele, nsp);
    m_formula_mat.setFromTriplets(trips.begin(), trips.end());
    // resize n vector and get state
    m_state = newNVector(m_laplace_size, m_sundials_ctx);
    m_constraints = newNVector(m_laplace_size, m_sundials_ctx);
    // define system constraints
    for (size_t i = 0; i < nsp; i++) {
        NV_Ith_S(m_constraints, i) = 1.0;
        NV_Ith_S(m_constraints, i + nsp + nele) = 1.0;
    }
    for (size_t i = nsp; i < nsp + nele; i++) {
        NV_Ith_S(m_constraints, i) = 0.0;
    }
    // create the necessary nonlinear algebraic solver in Sundials
    // free memory if it has been created
    KINFree(&m_kin_mem);
    SUNLinSolFree((SUNLinearSolver) m_linsol);
    SUNMatDestroy((SUNMatrix) m_linsol_matrix);
    // create memory
    #if CT_SUNDIALS_VERSION < 40
        // Create KINSOL Memory
        m_kin_mem = KINCreate();
        if (!m_kin_mem) {
            throw CanteraError("MultiPhaseEquilSolver::KINCreate", "KINCreate failed.");
        }
        // Create SUNMatrix
        m_linsol_matrix =  (void *) SUNDenseMatrix(m_laplace_size, m_laplace_size);
        // Create linear solver object
        #if CT_SUNDIALS_USE_LAPACK
            m_linsol = (void *) SUNLapackDense(m_state, (SUNMatrix) m_linsol_matrix);
        #else
            m_linsol = (void *) SUNDenseLinearSolver(m_state, (SUNMatrix) m_linsol_matrix);
        #endif
        // Set the solver as user data
        int flag = KINSetUserData(m_kin_mem, (void *) this);
        if (flag != 0) {
            throw CanteraError("MultiPhaseEquilSolver::KINSetUserData", "KINSetUserData failed.");
        }
        // Initialize KINSOL memory
        flag = KINInit(m_kin_mem, kin_rhs, m_state);
        if (flag != 0) {
            throw CanteraError("MultiPhaseEquilSolver::KINInit", "KINInit failed.");
        }
        // set constraints
        flag = KINSetConstraints(m_kin_mem, m_constraints);
        if (flag != 0) {
            throw CanteraError("MultiPhaseEquilSolver::KINSetConstraints", "KINSetConstraints failed.");
        }
        // Attach linear solver to memory
        flag = KINDlsSetLinearSolver(m_kin_mem, (SUNLinearSolver) m_linsol, (SUNMatrix) m_linsol_matrix);
        if (flag != 0) {
            throw CanteraError("MultiPhaseEquilSolver::KINSetLinearSolver", "KINSetLinearSolver failed.");
        }
    #else
        #if CT_SUNDIALS_VERSION < 60
            // Create KINSOL Memory
            m_kin_mem = KINCreate();
            if (!m_kin_mem) {
                throw CanteraError("MultiPhaseEquilSolver::KINCreate", "KINCreate failed.");
            }
            // Create SUNMatrix
            m_linsol_matrix =  (void *) SUNDenseMatrix(m_laplace_size, m_laplace_size);
            // Create linear solver object
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = (void *) SUNLinSol_LapackDense(m_state, (SUNMatrix) m_linsol_matrix);
            #else
                m_linsol = (void *) SUNLinSol_Dense(m_state, (SUNMatrix) m_linsol_matrix);
            #endif

        #else
            m_kin_mem = KINCreate(m_sundials_ctx.get());
            if (!m_kin_mem) {
                throw CanteraError("MultiPhaseEquilSolver::KINCreate", "KINCreate failed.");
            }
            // Create SUNMatrix
            m_linsol_matrix =  (void *) SUNDenseMatrix(m_laplace_size, m_laplace_size, m_sundials_ctx.get());
            // Create linear solver object
            #if CT_SUNDIALS_USE_LAPACK
                m_linsol = (void *) SUNLinSol_LapackDense(m_state, (SUNMatrix) m_linsol_matrix, m_sundials_ctx.get());
            #else
                m_linsol = (void *) SUNLinSol_Dense(m_state, (SUNMatrix) m_linsol_matrix, m_sundials_ctx.get());
            #endif
        #endif
            // Set the solver as user data
            int flag = KINSetUserData(m_kin_mem, (void *) this);
            if (flag != 0) {
                throw CanteraError("MultiPhaseEquilSolver::KINSetUserData", "KINSetUserData failed.");
            }
            // Initialize KINSOL memory
            flag = KINInit(m_kin_mem, kin_rhs, m_state);
            if (flag != 0) {
                throw CanteraError("MultiPhaseEquilSolver::KINInit", "KINInit failed.");
            }
            // set constraints
            flag = KINSetConstraints(m_kin_mem, m_constraints);
            if (flag != 0) {
                throw CanteraError("MultiPhaseEquilSolver::KINSetConstraints", "KINSetConstraints failed.");
            }
            // Attach linear solver to memory
            flag = KINSetLinearSolver(m_kin_mem, (SUNLinearSolver) m_linsol, (SUNMatrix) m_linsol_matrix);
            if (flag != 0) {
                throw CanteraError("MultiPhaseEquilSolver::KINSetLinearSolver", "KINSetLinearSolver failed.");
            }
    #endif
}

int MultiPhaseEquilSolver::evalEquilibrium()
{
    // set scaling vector - no scaling
    N_Vector scale = newNVector(m_laplace_size, m_sundials_ctx);
    N_VConst(1.0, scale);
    // set initial guess to current state
    N_VConst(0, m_state);
    m_mix->getMoles(NV_DATA_S(m_state));
    // Run solver
    int flag = KINSol(m_kin_mem, m_state, KIN_NONE, scale, scale);
    // destroy scaling vector
    N_VDestroy_Serial(scale);
    return flag;
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
    // right hand sides
    size_t nsp = m_mix->nSpecies();
    size_t nele = m_mix->nElements();
    // RHS vectors
    Eigen::Map<Eigen::VectorXd> n(RHS, nsp);
    Eigen::Map<Eigen::VectorXd> y(RHS + nsp, nele);
    Eigen::Map<Eigen::VectorXd> z(RHS + nsp + nele, nsp);
    // LHS vectors
    Eigen::Map<Eigen::VectorXd> fn(LHS, nsp);
    Eigen::Map<Eigen::VectorXd> fy(LHS + nsp, nele);
    Eigen::Map<Eigen::VectorXd> fz(LHS + nsp + nele, nsp);
    // getting chemical potentials and moles
    Eigen::VectorXd mu(nsp);
    Eigen::VectorXd b(nsp);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(nsp);
    m_mix->getMoles(b.data());
    m_mix->getChemPotentials(mu.data());
    // assignments
    fn = mu + m_formula_mat.transpose() * y - z;
    fy = m_formula_mat * n - b;
    fz = Eigen::MatrixXd(n.asDiagonal()) * Eigen::MatrixXd(z.asDiagonal()) * ones;
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
