/**
 *  @file MultiPhaseEquilSolver.cpp
 *    Driver routine for the VCSnonideal equilibrium solver package
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/equil/MultiPhaseEquilSolver.h"

#include <cstdio>

namespace Cantera
{

MultiPhaseEquilSolver::MultiPhaseEquilSolver(MultiPhase* mix) :
    m_mix(mix)
{
    // vector of all element compositions in all phases
    std::vector<Eigen::Triplet<double>> formulaVector;
    // get complete formula matrix matrix (elements x species)
    size_t formerElements = 0;
    size_t formerSpecies = 0;
    for (size_t i=0; i<mix->nPhases(); i++) {
        ThermoPhase& currPhase = mix->phase(i);
        // loop over all elements
        for (size_t j = 0; j < currPhase.nElements(); j++) {
            // loop over all species
            for (size_t k = 0; k < currPhase.nSpecies(); k++) {
                formulaVector.emplace_back(j+formerElements, k+formerSpecies,
                                           currPhase.nAtoms(k, j));
            }
        }
        // add to former terms for indexing
        formerElements += currPhase.nElements();
        formerSpecies += currPhase.nSpecies();
    }
    // set formula matrix from triplets
    m_formula_matrix.resize(mix->nElements(), mix->nSpecies());
    m_formula_matrix.setFromTriplets(formulaVector.begin(), formulaVector.end());
}

int MultiPhaseEquilSolver::equilibrate_TV(int XY, double xtarget, int estimateEquil,
                                        double err, int maxsteps, int loglevel)
{

}

int MultiPhaseEquilSolver::equilibrate_HP(double Htarget, int XY, double Tlow,
    double Thigh, int estimateEquil, double err, int maxsteps,
    int loglevel)
{
     return 0;
}

int MultiPhaseEquilSolver::equilibrate_SP(double Starget, double Tlow, double Thigh,
    int estimateEquil, double err, int maxsteps, int loglevel)
{
    return 0;
}

int MultiPhaseEquilSolver::equilibrate(int XY, int estimateEquil, double err, int maxsteps, int loglevel)
{
    return 0;
}

int MultiPhaseEquilSolver::equilibrate_TP(int estimateEquil, double err,
                                        int maxsteps, int loglevel)
{
    return 0;
}

}
