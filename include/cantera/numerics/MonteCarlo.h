/**
 *  @file MonteCarlo.h Declarations for the class
 *  MonteCarlo which is a novel implementation and solver type so it does
 *  not inherit from other types.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "cantera/numerics/eigen_sparse.h"
#include "cantera/base/global.h"

namespace Cantera
{

//! MonteCarlo is a solver designed for simulating particle phase physics with Monte
//! Carlo Method.
class MonteCarlo
{
private:
    /* data */
public:
    MonteCarlo(/* args */) {};
    ~MonteCarlo() {};

};

}
