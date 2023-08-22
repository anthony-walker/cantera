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

    // function to setup hessian approximation as the Jacobian
    static int kin_jac(N_Vector x, N_Vector fx, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2) {
        MultiPhaseEquilSolver* solver = (MultiPhaseEquilSolver*) user_data;
        solver->setupHessianJac(J, NV_DATA_S(x), NV_DATA_S(fx));
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
    m_mix(mix) {
}

void MultiPhaseEquilSolver::initialize()
{
    // Sizes needed
    size_t N = m_mix->nSpecies();
    size_t C = m_mix->nElements();
    m_laplace_size = 2 * N + C;
    // get complete formula matrix matrix (elements x species)
    vector<Eigen::Triplet<double>> trips;
    // reserve some initial space
    trips.reserve(N);
    // loop over all elements - rows
    for (size_t j = 0; j < C; j++) {
        // loop over all species - columns
        for (size_t k = 0; k < N; k++) {
            double atoms = m_mix->nAtoms(k, j);
            if (atoms > 0) {
                trips.emplace_back(j, k, atoms);
            }
        }
    }
    // resize eigen formula matrix
    m_formula_mat.resize(C, N);
    m_formula_mat.setFromTriplets(trips.begin(), trips.end());
    // free matrix memory if it exists
    N_VDestroy_Serial(m_state);
    N_VDestroy_Serial(m_constraints);
    // resize n vector and get state
    m_state = newNVector(m_laplace_size, m_sundials_ctx);
    m_constraints = newNVector(m_laplace_size, m_sundials_ctx);
    // define system constraints
    for (size_t i = 0; i < N; i++) {
        NV_Ith_S(m_constraints, i) = 1.0;
        NV_Ith_S(m_constraints, i + N + C) = 1.0;
    }
    for (size_t i = N; i < N + C; i++) {
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
            // Create the new SUNMatrix
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
            // Attach jac function
            flag = KINSetJacFn(m_kin_mem, kin_jac);
            if (flag != 0) {
                throw CanteraError("MultiPhaseEquilSolver::KINSetJacFn", "KINSetJacFn failed.");
            }
            // Setting print level to be verbose
            // FIXME: It seems that the initial guess may be the problem?
            KINSetPrintLevel(m_kin_mem, 3);
    #endif
}

int MultiPhaseEquilSolver::evalEquilibrium()
{
    for (size_t i = 0; i < m_mix->nSpecies(); i++) {
        std::cout<<m_mix->speciesName(i)<<" ";
    }
    std::cout<<std::endl;

    std::cout<<m_formula_mat<<std::endl;
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

void MultiPhaseEquilSolver::setupHessianJac(SUNMatrix J, double* u, double* fu)
{
    // constants
    size_t N = m_mix->nSpecies();
    size_t C = m_mix->nElements();
    // set first diagonal to approximate Hessian
    for (size_t i = 0; i < N; i++) {
        SM_ELEMENT_D(J, i, i) = u[i] > 0 ? 1 / u[i] : 0;            /* code */
    }
    // add formula matrix to Jacobian
    for (int k = 0; k < m_formula_mat.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_formula_mat, k); it; ++it) {
            SM_ELEMENT_D(J, it.row() + N, it.col()) = it.value();
            SM_ELEMENT_D(J, it.col() + N, it.row()) = it.value();
        }
    }
    // upper right quadrant is negative identity
    // lower left quadrant is Z matrix
    // lower right quadrant is N matrix
    for (size_t i = 0; i < N; i++) {
        SM_ELEMENT_D(J, i, i + N + C) = -1;
        SM_ELEMENT_D(J, i + N + C, i) = u[i+N+C];
        SM_ELEMENT_D(J, i + N + C, i + N + C) = u[i];
    }
}

int MultiPhaseEquilSolver::equilibrate_TP(double* LHS, double* RHS, void* f_data)
{
    // right hand sides
    size_t N = m_mix->nSpecies();
    size_t C = m_mix->nElements();
    // RHS vectors
    Eigen::Map<Eigen::VectorXd> n(RHS, N);
    Eigen::Map<Eigen::VectorXd> y(RHS + N, C);
    Eigen::Map<Eigen::VectorXd> z(RHS + N + C, N);
    // LHS vectors
    Eigen::Map<Eigen::VectorXd> fn(LHS, N);
    Eigen::Map<Eigen::VectorXd> fy(LHS + N, C);
    Eigen::Map<Eigen::VectorXd> fz(LHS + N + C, N);
    // getting chemical potentials and moles
    Eigen::VectorXd mu(N);
    Eigen::VectorXd b(C);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(N);
    m_mix->getElemAbundances(b.data());
    m_mix->getChemPotentials(mu.data());
    // assignments
    std::cout<<y<<std::endl;
    // std::cout<<n<<std::endl;
    fn = m_formula_mat.transpose() * y + z - mu;
    fy = b - m_formula_mat * n;
    fz = - Eigen::MatrixXd(n.asDiagonal()) * Eigen::MatrixXd(z.asDiagonal()) * ones;
    // std::cout<<fn<<std::endl;
    // std::cout<<"----------"<<std::endl;
    // std::cout<<fy<<std::endl;
    // std::cout<<fz<<std::endl;
    // std::cout<<"----------"<<std::endl;
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
