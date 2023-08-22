/**
 *  @file  MultiPhaseEquilSolver.h
 *  Interface class for the vcsnonlinear solver
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef MULTIPHASEEQUILSOLVER_H
#define MULTIPHASEEQUILSOLVER_H

#include "MultiPhase.h"
#include "cantera/numerics/sundials_headers.h"
#include "cantera/numerics/SundialsContext.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/numerics/eigen_dense.h"

namespace Cantera
{
/**
 * Specifies the side of the system on which the preconditioner is applied. Not all
 * methods are supported by all integrators.
 */
enum class ConstantState {
    XY, //! Invalid constant state which an invalid string is provided
    TV, //! Temperature and volume are held constant
    TP, //! Temperature and pressure are held constant
    UV, //! Internal energy and volume are held constant
    HP, //! Enthalpy and pressure are held constant
    SP,  //! Entropy and pressure are held constant
};

//! %Cantera's Interface to the Multiphase chemical equilibrium solver.
/*!
 * Class MultiPhaseEquilSolver is designed to be used to set a mixture containing
 * one or more phases to a state of chemical equilibrium.
 *
 * Note, as currently constructed, the underlying ThermoPhase objects are shared
 * between the MultiPhase object and this object. Therefore, mix is not a const
 * argument, and the return parameters are contained in underlying ThermoPhase
 * objects.
 *
 * @ingroup equilGroup
 */
class MultiPhaseEquilSolver
{
public:
    //! Constructor for the multiphase equilibrium solver
    /*!
     * This constructor will initialize the object with a MultiPhase object,
     * setting up the internal equilibration problem. Note, as currently
     * constructed, the underlying ThermoPhase objects are shared between the
     * MultiPhase object and this object. Therefore, mix is not a const
     * argument, and the return parameters are contained in underlying
     * ThermoPhase objects.
     *
     * @param mix Object containing the MultiPhase object
     */
    MultiPhaseEquilSolver(MultiPhase* mix);

    virtual ~MultiPhaseEquilSolver() {}

    //! return the number of nonlinear solver iterations
    int iterations() const {
        long int iters = 0;
        KINGetNumNonlinSolvIters(m_kin_mem, &iters);
        return iters;
    }

    /*! Equilibrate the solution using the current element abundances stored in the
     * MultiPhase object
     * @param LHS left hand side vector F(u)
     * @param RHS right hand side vector u
     * @param f_data a pointer to data passed into the function for solving the problem
     */
    int equilibrate(double* LHS, double* RHS, void* f_data);

    /*! Interface function to evaluate equilibrium based on the given state
     * @param xy a string TP, TV, etc specifying variables to be held constant.
     */
    int evalEquilibrium(const string xy);

    /*! Internal function to setup the Jacobian approximation using approximate Hessian
    * @param J SUNMatrix that is the output Jacobian in the calculation
    * @param u current iterate of u within the solver
    * @param fu current evaluation of F(u) from equilibrate function
    */
    void setupHessianJac(SUNMatrix J, double* u, double* fu);

    //! A function to initialize the solver based on `m_mix`.
    void initialize();

    //! A function to check if the solver has been initialized.
    bool initialized() {
        return m_init;
    }

protected:
    // structures for sparse kinsol solution
    void* m_kin_mem = nullptr; //!< Pointer to the KINSOL memory for the problem
    void* m_linsol = nullptr; //!< Sundials linear solver object
    void* m_linsol_matrix = nullptr; //!< matrix used by Sundials
    SundialsContext m_sundials_ctx; //!< SUNContext object for Sundials>=6.0
    N_Vector m_state; //! State of system in an n-vector for solver
    N_Vector m_constraints; //! Constraints on system
    Eigen::SparseMatrix<double> m_formula_mat; //! Formula matrix used in equilibrate
    Eigen::MatrixXd m_temp_jac; //! Formula matrix used in equilibrate
    Eigen::VectorXd m_lagrange_z; //! Lagrange multiplier z vector
    double m_perturb = 1e-25; //! default perturbation parameter
    ConstantState m_xy = ConstantState::XY; //! default constant variables set to invalid state
    // // An index map to dynamically restructure the order for better Gaussian elimination
    // AnyMap m_index_map;
    //! Pointer to the MultiPhase mixture that will be equilibrated.
    /*!
     *  Equilibrium solutions will be returned via this variable.
     */
    MultiPhase* m_mix;
    //! Initialization bool
    bool m_init = false;
    //! Size of state with laplace transform
    size_t m_sys_size;
    // constant variable one
    double m_const_one;
    //! constant variable two
    double m_const_two;
};

}

#endif
