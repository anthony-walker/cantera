/**
 * @Author: Anthony Walker
 * @Date:   2020-04-16T14:30:22-07:00
 * @Email:  walkanth@oregonstate.edu
 * @Filename: AdaptivePreconditioner.h
 * @Last modified time: 2020-04-16T15:33:18-07:00
 */

#ifndef ADAPTIVEPRECONDITIONER_H
#define ADAPTIVEPRECONDITIONER_H
#endif

//Cantera imports
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"

//Eigen Imports
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

template<typename _Scalar, int _Options, typename _StorageIndex> void AdaptivelyPrecondition(Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex> *sparse_matrix)
{
  /*
    This is the main preconditioner function which takes a SparseMatrix created by Eigen of the appropriate size.
  */
  // sparse_matrix->insert(0,0) = 1;
  for (int i = 0; i < sparse_matrix->cols(); i++) {
    sparse_matrix->insert(i,i) = i;
  }
}
