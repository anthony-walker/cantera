/**
 *  @file AdaptivePreconditioner.h Declarations for the class
 *   AdaptivePreconditioner which is a child class of PreconditionerBase
 *   for preconditioners used by sundials
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef ADAPTIVEPRECONDITIONER_H
#define ADAPTIVEPRECONDITIONER_H

#include "cantera/numerics/PreconditionerBase.h"
#include "cantera/numerics/eigen_sparse.h"
#include <iostream>

namespace Cantera
{

//! AdaptivePreconditioner a preconditioner designed for use with large
//! mechanisms that leverages sparse solvers. It does this by pruning
//! the preconditioner by a threshold value. It also neglects pressure
//! dependence and thirdbody contributions in its formation and has a
//! finite difference approximation for temperature.
class AdaptivePreconditioner : public PreconditionerBase
{
public:
    AdaptivePreconditioner() {}

    //! This function is called during setup for any processes that need
    //! to be completed prior to setup functions used in sundials.
    //! @param network A pointer to the reactor net object associated
    //! with the integration
    void initialize(size_t networkSize);

    //! Reset arrays within preconditioner object
    void reset() {
        m_precon_matrix.setZero();
        m_jac_trips.clear();
    };

    //! Perform preconditioner specific post-reactor
    //! setup operations such as factorize.
    void setup();

    //! Transform Jacobian vector and write into
    //! preconditioner
    void transformJacobianToPreconditioner();

    //! Prune preconditioner elements
    void prunePreconditioner();

    //! Solve a linear system Ax=b where A is the preconditioner
    //! @param[in] stateSize length of the rhs and output vectors
    //! @param[in] rhs_vector right hand side vector used in linear system
    //! @param[out] output output vector for solution
    void solve(const size_t stateSize, double *rhs_vector, double* output);

    //! Return the preconditioning type as an integer
    PreconditionerType preconditionerType() { return PreconditionerType::LEFT_PRECONDITION; }

    //! Function used to return semi-analytical jacobian matrix
    Eigen::SparseMatrix<double> getJacobian() {
        Eigen::SparseMatrix<double> jacobian(m_dimensions[0], m_dimensions[1]);
        jacobian.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
        return jacobian;
    }

    //! Get the threshold value for setting
    //! elements
    double threshold() { return m_threshold; }

    //! Get ilut fill factor
    double ilutFillFactor() { return m_fill_factor; }

    //! Get ilut drop tolerance
    double ilutDropTol() { return m_drop_tol; }

    //! Set the threshold value to compare elements against
    //! @param threshold double value used in setting by threshold
    void setThreshold(double threshold) {
        m_threshold = threshold;
        m_prune_precon = (threshold <= 0) ? false : true;
    }

    //! Set drop tolerance for ILUT
    //! @param droptol double value used in setting solver drop tolerance
    void setIlutDropTol(double droptol) {
        m_drop_tol = droptol;
        m_solver.setDroptol(droptol);
        }

    //! Set the fill factor for ILUT solver
    //! @param fillFactor fill in factor for ILUT solver
    void setIlutFillFactor(int fillFactor) {
        m_fill_factor = fillFactor;
        m_solver.setFillfactor(fillFactor);
    }

    //! Overloading of the () operator to assign values to the jacobian
    //! this function does not assume that index is index map
    //! @param row row index of jacobian
    //! @param col column index of jacobian
    //! @param value to place in jacobian vector
    const double& operator() (size_t row, size_t col);

    //! Overloading of the () operator to assign values to the jacobian
    //! this function does not assume that index is index map
    //! @param row row index of jacobian
    //! @param col column index of jacobian
    //! @param value to place in jacobian vector
    void operator() (size_t row, size_t col, double value);

    //! Print preconditioner contents
    void printPreconditioner() {
        std::stringstream ss;
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        ss << Eigen::MatrixXd(m_precon_matrix).format(HeavyFmt);
        writelog(ss.str());
    }

    //! Print jacobian contents
    void printJacobian() {
        std::stringstream ss;
        Eigen::SparseMatrix<double> jacobian(m_dimensions[0], m_dimensions[1]);
        jacobian.setFromTriplets(m_jac_trips.begin(), m_jac_trips.end());
        ss << Eigen::MatrixXd(jacobian);
        writelog(ss.str());
    }

protected:
    //! ilut fill factor
    double m_fill_factor = 0;

    //! ilut drop tolerance
    double m_drop_tol = 0;

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
    double m_threshold = 1e-8;

    //! Bool set whether to prune the matrix or not
    double m_prune_precon = true;
};

}

#endif
