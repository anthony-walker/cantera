/**
 *  @file AdaptivePreconditioner.h Declarations for the class
 *   AdaptivePreconditioner which is a child class of PreconditionerBase
 *   for preconditioners used by sundials
 */

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#ifndef ADAPTIVEPRECONDITIONER_H
#define ADAPTIVEPRECONDITIONER_H

#include "cantera/numerics/eigen_sparse.h"
#include "cantera/numerics/PreconditionerBase.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "float.h"
#include <unordered_map>
#include <iostream>

namespace Cantera
{

//! Flag to indicate adaptive preconditioner is set
const int ADAPTIVE_MECHANISM_PRECON_MATRIX = 1;

//! AdaptivePreconditioner a preconditioner designed for use with large
//! mechanisms that leverages sparse solvers. It does this by pruning
//! the preconditioner by a threshold value. It also neglects pressure
//! dependence and thirdbody contributions in its formation and has a
//! finite difference approximation for temperature.
class AdaptivePreconditioner : public PreconditionerBase
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // Required for mis-alignment of EIGEN matrix
    AdaptivePreconditioner(){};
    ~AdaptivePreconditioner(){};
    AdaptivePreconditioner(const AdaptivePreconditioner &externalPrecon){};

    //! This function is called during setup for any processes that need
    //! to be completed prior to setup functions used in sundials.
    //! @param network A pointer to the reactor net object associated
    //! with the integration
    void initialize(ReactorNet& network);

    //! Use this function to reset arrays within preconditioner object
    void reset(){
        m_precon_matrix.setZero();
        m_jac_trips.clear();
    };

    //! Use this function to perform preconditioner specific post-reactor
    //! setup operations such as factorize.
    void setup();

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with MoleReactor
    void acceptReactor(MoleReactor& reactor, double t, double* LHS, double* RHS);

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with IdealGasMoleReactor
    void acceptReactor(IdealGasMoleReactor& reactor, double t, double* LHS, double* RHS);

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with IdealGasConstPressureMoleReactor
    void acceptReactor(IdealGasConstPressureMoleReactor& reactor, double t, double* LHS, double* RHS);

    //! Use this function to check if there was an error in eigen methods and throw
    //! it if so.
    void preconditionerErrorCheck();

    //! Use this function to transform Jacobian vector and write into
    //! preconditioner
    void transformJacobianToPreconditioner();

    //! Use this function to prune preconditioner elements
    void prunePreconditioner();

    //! Function to solve a linear system Ax=b where A is the
    //! preconditioner contained in this matrix
    //! @param[in] state_len length of vectors supplied
    //! @param[in] rhs_vector right hand side vector supplied by cvodes
    //! @param[out] output output vector "z" sent back to cvodes
    void solve(const size_t state_len, double *rhs_vector, double* output);

    //! Use this function to return the preconditioning method as an integer
    size_t getPreconditionerMethod(){return ADAPTIVE_MECHANISM_PRECON_MATRIX;};

    //! Use this function to return the preconditioning type as an integer
    PreconditionerType getPreconditionerType(){return LEFT_PRECONDITION;};

    //! Use this function to return pointer to the preconditioner matrix
    Eigen::SparseMatrix<double>* getMatrix(){return &(m_precon_matrix);};

    //! Function used to return semi-analytical jacobian matrix
    Eigen::SparseMatrix<double> getJacobian(){
        Eigen::SparseMatrix<double> jacobian(m_dimensions[0], m_dimensions[1]);
        jacobian.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
        return jacobian;
    };

    //! Use this function to get the threshold value for setting
    //! elements
    double getThreshold(){return m_threshold;};

    //! Use this function to get the pertubation constant
    double getPerturbationConst(){return m_perturb;};

    //! Use this function to get a strictly positive composition
    void getStrictlyPositiveComposition(size_t vlen, double* in, double* out){
        for (size_t i = 0; i < vlen; i++)
        {
            out[i] = std::max(in[i], m_atol);
        }
    };

    //! Use this function to set the threshold value to compare elements
    //! against
    //! @param threshold double value used in setting by threshold
    void setThreshold(double threshold)
    {
        m_threshold = threshold;
        m_prune_precon = (threshold <= 0) ? false : true;
    };

    //! Use this function to set drop tolerance for ILUT
    //! @param droptol double value used in setting solver drop tolerance
    void setDropTolILUT(double droptol = 1e-10){m_solver.setDroptol(droptol);};

    //! Use this function to set the fill factor for ILUT
    void setFillFactorILUT(int fillfactor = -1)
    {
        int newfillfactor = (fillfactor < 0) ? m_dimensions[0]/4 : fillfactor;
        m_solver.setFillfactor(newfillfactor);
    }

    //! Use this function to set the perturbation constant used in
    //! finite difference calculations.
    //! @param perturb the new pertubation constant
    void setPerturbationConst(double perturb){m_perturb = perturb;};

    //! Overloading of the == operator to compare values strictly inside
    //! preconditioner matrix
    //! @param externalPrecon == comparison with this object
    bool operator== (const AdaptivePreconditioner &externalPrecon);
    //! Overloading of the = operator to copy one preconditioner to
    //! another
    //! @param externalPrecon the preconditioner becoming this object
    void operator= (const AdaptivePreconditioner &externalPrecon);

    //! Overloading of the () operator to assign values to the jacobian
    //! this function does not assume that index is index map
    //! @param row row index of jacobian
    //! @param col column index of jacobian
    //! @param value to place in jacobian vector
    void operator() (size_t row, size_t col, double value);

    //! Use this function to print preconditioner contents
    void printPreconditioner(){
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        std::cout<<Eigen::MatrixXd(m_precon_matrix).format(HeavyFmt)<<std::endl;
    };

    //! Use this function to print jacobian contents
    void printJacobian(){
        Eigen::SparseMatrix<double> jacobian(m_dimensions[0], m_dimensions[1]);
        jacobian.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
        std::cout<<Eigen::MatrixXd(jacobian)<<std::endl;
    };

    //! Use this function to set the fill factor for factorizing the ILUT preconditioner
    void setFillFactor(int n) {
        m_solver.setFillfactor(n);
    }

    //! Use this function to set the tolerance for dropping elements when factorizing the ILUT preconditioner
    void setDropTol(double tol) {
        m_solver.setDroptol(tol);
    }

protected:
    //! Vector of triples representing the jacobian used in preconditioning
    std::vector<Eigen::Triplet<double>> m_jac_trips;

    //! Storage of appropriately sized identity matrix for making the preconditioner
    Eigen::SparseMatrix<double> m_identity;

    //! Container that is the sparse preconditioner
    Eigen::SparseMatrix<double> m_precon_matrix;

    //! Solver used in solving the linear system
    Eigen::IncompleteLUT<double> m_solver;

    //! Minimum value a non-diagonal element must be to be included in
    //! the preconditioner
    double m_threshold = DBL_EPSILON; // default

    //! Perturbation constant that is multiplied by temperature for
    //! perturbation and finite difference calculations
    double m_perturb = std::sqrt(DBL_EPSILON);

    //! Bool set whether to prune the matrix or not
    double m_prune_precon = true;
};

}

#endif
