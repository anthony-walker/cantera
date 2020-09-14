/*
Programmer: Anthony Walker
This is the this file contains functions to adaptively precondition the sparse matrix class
*/
#ifndef PRECONDITIONING_LIBRARY_H
#define PRECONDITIONING_LIBRARY_H

//Cantera imports
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/kinetics/Reaction.h"
#include "Preconditioners.h"

//Function Library for preconditioniong
namespace Cantera::AMP //Making ASP apart of Cantera namespace
{
  inline void printReactorComponents(Reactor* reactor);

  inline void speciesDerivative(std::map<std::string, double> comp,std::map<std::string,size_t> indexMap, double* omega, double* concentrations, double k_direction, double volume);

  void TemperatureDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, double* ydot, double dTdt, size_t index, size_t speciesStart);
  
  void SpeciesSpeciesDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, size_t speciesStart);

}

#endif