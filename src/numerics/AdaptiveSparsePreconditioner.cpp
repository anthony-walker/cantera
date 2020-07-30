/**
 * @Author: Anthony Walker
 * @Date:   2020-04-16T14:30:22-07:00
 * @Email:  walkanth@oregonstate.edu
 * @Filename: AdaptivePreconditioner.h
 * @Last modified time: 2020-04-16T15:33:18-07:00
 */

//Cantera imports
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
//Adaptive imports
#include "cantera/numerics/AdaptiveSparsePreconditioner.h"
#include "cantera/numerics/SparseMatrix.h"
//Sundials imports
#include "sunmatrix/sunmatrix_sparse.h"
//Eigen Imports
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

template<class MATTYPE> void AdaptivelyPrecondition(MATTYPE *preconditioner)
{
  /*
    This is the main preconditioner function which takes a SparseMatrix created by Eigen of the appropriate size.
  */
}



template<class MATTYPE> void TemperatureSpeciesDerivative(MATTYPE *preconditioner)
{
  /*
    This is the main preconditioner function which takes a SparseMatrix created by Eigen of the appropriate size.
  */
}

