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
        solver->equilibrate(NV_DATA_S(fx), NV_DATA_S(x), user_data);
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
    m_sys_size = N + C;
    m_temp_jac.resize(m_sys_size, m_sys_size);
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
    // setup perturbation
    vector<double> elemAbundence(C);
    m_mix->getElemAbundances(elemAbundence.data());
    int min = npos;
    for (auto val : elemAbundence) {
        if (val < min && val > 0) {
            min = val;
        }
    }
    m_perturb *= min;

    // lagrange vector should be size of number of species
    m_lagrange_z.resize(N);
    m_mix->getMoles(m_lagrange_z.data());
    // determining an initial z vector
    for (size_t i = 0; i < N; i++) {
        m_lagrange_z[i] = m_lagrange_z[i] > 0 ? m_perturb / m_lagrange_z[i] : 1;
    }
    // free matrix memory if it exists
    N_VDestroy_Serial(m_state);
    N_VDestroy_Serial(m_constraints);
    // resize n vector and get state
    m_state = newNVector(m_sys_size, m_sundials_ctx);
    m_constraints = newNVector(m_sys_size, m_sundials_ctx);
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
        m_linsol_matrix =  (void *) SUNDenseMatrix(m_sys_size, m_sys_size);
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
            m_linsol_matrix =  (void *) SUNDenseMatrix(m_sys_size, m_sys_size);
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
            m_linsol_matrix =  (void *) SUNDenseMatrix(m_sys_size, m_sys_size, m_sundials_ctx.get());
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
            // // Attach jac function
            // flag = KINSetJacFn(m_kin_mem, kin_jac);
            // if (flag != 0) {
            //     throw CanteraError("MultiPhaseEquilSolver::KINSetJacFn", "KINSetJacFn failed.");
            // }
            // // set scaled step tolerance
            // flag = KINSetScaledStepTol(m_kin_mem, 1e-20);
            // if (flag != 0) {
            //     throw CanteraError("MultiPhaseEquilSolver::KINSetJacFn", "KINSetJacFn failed.");
            // }
            // // set scaled step tolerance
            // flag = KINSetFuncNormTol(m_kin_mem, 1e8);
            // if (flag != 0) {
            //     throw CanteraError("MultiPhaseEquilSolver::KINSetJacFn", "KINSetJacFn failed.");
            // }
            // Setting print level to be verbose
            // FIXME: It seems that the initial guess may be the problem?
            KINSetPrintLevel(m_kin_mem, 3);
    #endif
}

int MultiPhaseEquilSolver::evalEquilibrium(const string xy)
{
    if (xy == "TP") {
        m_xy = ConstantState::TP;
        m_const_one = m_mix->temperature();
        m_const_two = m_mix->pressure();
    } else if (xy == "TV") {
        m_xy = ConstantState::TV;
        m_const_one = m_mix->temperature();
        m_const_two = m_mix->volume();
    } else if (xy == "SP") {
        m_xy = ConstantState::SP;
        m_const_one = m_mix->entropy();
        m_const_two = m_mix->pressure();
    } else if (xy == "HP") {
        m_xy = ConstantState::HP;
        m_const_one = m_mix->enthalpy();
        m_const_two = m_mix->pressure();
    } else if (xy == "UV") {
        m_xy = ConstantState::UV;
        m_const_one = m_mix->IntEnergy();
        m_const_two = m_mix->volume();
    } else {
        throw CanteraError("MultiPhaseEquilSolver::evalEquilibrium", "Invalid constant state given.");
    }
    // set scaling vector - no scaling
    N_Vector scale = newNVector(m_sys_size, m_sundials_ctx);
    N_VConst(1.0, scale);
    // set initial guess to current state
    N_VConst(0, m_state);
    m_mix->getMoles(NV_DATA_S(m_state));
    // set state of lagrange multiplier y
    for (size_t i = m_mix->nSpecies(); i < m_sys_size; i++) {
        NV_Ith_S(m_state, i) = 0;
    }
    // Run solver
    int flag = KINSol(m_kin_mem, m_state, KIN_NONE, scale, scale);
    // destroy scaling vector
    N_VDestroy_Serial(scale);
    return flag;
    return 0;
}

int MultiPhaseEquilSolver::equilibrate(double* LHS, double* RHS, void* f_data)
{
    // sizes used in computation
    size_t N = m_mix->nSpecies();
    size_t C = m_mix->nElements();
    // RHS vectors
    Eigen::Map<Eigen::VectorXd> n(RHS, N);
    Eigen::Map<Eigen::VectorXd> y(RHS + N, C);
    // vector for the change in n
    // LHS vectors
    Eigen::Map<Eigen::VectorXd> fn(LHS, N);
    Eigen::Map<Eigen::VectorXd> fy(LHS + N, C);
    // Vectors for getting chemical potentials and moles
    Eigen::VectorXd mu(N);
    Eigen::MatrixXd nInv(N, N);
    Eigen::VectorXd b(C);
    Eigen::VectorXd dn(N);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(N);
    // get former state
    m_mix->getMoles(dn.data());
    // create inverse matrix
    nInv = dn.asDiagonal();
    for (size_t i = 0; i < N; i++) {
         nInv(i, i) += (nInv(i, i) == 0) ? m_perturb : 0;
    }
    nInv = nInv.inverse();
    // set state according to m_xy
    if (m_xy == ConstantState::TP) {
        m_mix->setState_TP(m_const_one, m_const_two);
    } else if (m_xy == ConstantState::TV) {
        // FIXME: These set state functions don't exist
        // m_mix->setState_TV(m_const_one, m_const_two);
    } else if (m_xy == ConstantState::SP) {
        // m_mix->setState_SP(m_const_one, m_const_two);
    } else if (m_xy == ConstantState::HP) {
        // m_mix->setState_HP(m_const_one, m_const_two);
    } else {
        throw CanteraError("MultiPhaseEquilSolver::equilibrate", "Invalid constant state given.");
    }
    // set moles of the system
    m_mix->setMoles(RHS);
    // get system parameters
    m_mix->getElemAbundances(b.data());
    m_mix->getChemPotentials(mu.data());
    // make dn
    dn = n - dn;
    // assignments
    fn = mu - m_formula_mat.transpose() * y - m_perturb * nInv * ones;
    fy = m_formula_mat * n - b;
    std::cout<<fn<<std::endl;
    std::cout<<fy<<std::endl;
    m_lagrange_z = (m_perturb * nInv * ones - nInv * m_lagrange_z.asDiagonal() * dn);
    return 0;
}

void MultiPhaseEquilSolver::setupHessianJac(SUNMatrix J, double* u, double* fu)
{
    // constants
    size_t N = m_mix->nSpecies();
    size_t C = m_mix->nElements();
    // set first diagonal to approximate Hessian
    for (size_t i = 0; i < N; i++) {
        double n_inv = u[i] != 0 ? 1 / u[i] : 1 / m_perturb;
        SM_ELEMENT_D(J, i, i) = n_inv + n_inv * m_lagrange_z[i];
        m_temp_jac(i, i) = n_inv + n_inv * m_lagrange_z[i];
    }
    // add formula matrix to Jacobian
    for (int k = 0; k < m_formula_mat.outerSize(); k++) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_formula_mat, k); it; ++it) {
            SM_ELEMENT_D(J, it.row() + N, it.col()) = it.value();
            SM_ELEMENT_D(J, it.col(), it.row() + N) = -it.value();
            m_temp_jac(it.row() + N, it.col()) = it.value();
            m_temp_jac(it.col(), it.row() + N) = -it.value();
        }
    }
}

}
