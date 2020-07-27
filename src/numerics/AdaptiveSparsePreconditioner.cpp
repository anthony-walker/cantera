/**
 * @Author: Anthony Walker
 * @Date:   2020-04-16T14:30:22-07:00
 * @Email:  walkanth@oregonstate.edu
 * @Filename: AdaptivePreconditioner.h
 * @Last modified time: 2020-04-16T15:33:18-07:00
 */


//Cantera imports
#include "cantera/numerics/AdaptiveSparsePreconditioner.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"

template<class MATTYPE> void Cantera::AdaptivelyPrecondition(SparseMatrix<MATTYPE> *sparseMatrix)
{
  /*
    This is the main preconditioner function which takes a SparseMatrix created by Eigen of the appropriate size.
  */
}

template<class MATTYPE> void Cantera::TemperatureSpeciesDerivative(SparseMatrix<MATTYPE> *sparseMatrix)
{
  /*
    This is the main preconditioner function which takes a SparseMatrix created by Eigen of the appropriate size.
  */
}

