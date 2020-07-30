/*
Programmer: Anthony Walker
This is the this file contains functions to adaptively precondition the sparse matrix class
*/
#ifndef ADAPTIVESPARSEPRECONDITIONER_H
#define ADAPTIVESPARSEPRECONDITIONER_H

//Cantera imports
#include "cantera/numerics/SparseMatrix.h"
//Sundials imports
#include "sundials/sundials_matrix.h"
#include "sunmatrix/sunmatrix_sparse.h"
//Eigen Imports
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif
//Type definitions
typedef Eigen::SparseMatrix<double> EigenSparseMatrix;
typedef SUNMatrix SundialsSparseMatrix;

//Function declarations
namespace Cantera //Making ASP apart of Cantera namespace
{

template<class MATTYPE> void AdaptivelyPrecondition(MATTYPE *preconditioner);

template<class MATTYPE> void TemperatureSpeciesDerivative(MATTYPE *preconditioner);

}

#endif