/*
Programmer: Anthony Walker
This is the this file contains functions to adaptively precondition the sparse matrix class
*/
#ifndef ADAPTIVESPARSEPRECONDITIONER_H
#define ADAPTIVESPARSEPRECONDITIONER_H

//Cantera imports
#include "cantera/numerics/SparseMatrix.h"

//Function Library for preconditioniong
namespace Cantera //Making ASP apart of Cantera namespace
{


template<class MATTYPE> void AdaptivelyPrecondition(SparseMatrix<MATTYPE> *preconditioner)
{
  /*
    This is the main preconditioner function which takes a SparseMatrix created by Eigen of the appropriate size.
  */
}



template<class MATTYPE> void TemperatureSpeciesDerivative(SparseMatrix<MATTYPE> *preconditioner)
{
  /*
    This is the main preconditioner function which takes a SparseMatrix created by Eigen of the appropriate size.
  */
}

}
#endif