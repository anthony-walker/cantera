/*
Programmer: Anthony Walker
This is the this file contains functions to adaptively precondition the sparse matrix class
*/
#ifndef ADAPTIVE_SPARSE_PRECONDITIONER_H
#define ADAPTIVE_SPARSE_PRECONDITIONER_H

const int ADAPTIVE_MECHANISM_PRECONDITIONER = 1;
const int PRECONDITIONER_NOT_SET = 0;


//Cantera imports
#include "cantera/numerics/SparseMatrix.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/kinetics/Reaction.h"
//Sundials imports
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "nvector/nvector_serial.h"

#if CT_SUNDIALS_VERSION < 25
typedef int sd_size_t;
#else
typedef long int sd_size_t;
#endif
//Other imports

//Function Library for preconditioniong
namespace Cantera::AMP //Making ASP apart of Cantera namespace
{
  inline void printReactorComponents(Reactor* reactor);

  inline void speciesDerivative(std::map<std::string, double> comp,std::map<std::string,size_t> indexMap, double* omega, double* concentrations, double k_direction, double volume);

  void TemperatureDerivatives(SparseMatrix *preconditioner,Reactor* reactor, double* ydot, double dTdt, size_t index, size_t speciesStart);
  
  void SpeciesSpeciesDerivatives(SparseMatrix *preconditioner,Reactor* reactor, size_t speciesStart);

}

#endif