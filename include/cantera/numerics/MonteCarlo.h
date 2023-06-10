/**
 *  @file MonteCarlo.h Declarations for the class
 *  MonteCarlo which is a novel implementation and solver type so it does
 *  not inherit from other types.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "cantera/numerics/eigen_dense.h"
#include "cantera/base/global.h"

namespace Cantera
{

//! MonteCarlo is a solver designed for simulating particle phase physics with Monte
//! Carlo Method.
class MonteCarlo {
public:
    MonteCarlo(/* args */) {};
    ~MonteCarlo() {};
    void initialize(size_t nsp, size_t bins);

protected:
    /* data */
    size_t m_bins;
    size_t m_nsp;
    double m_stdev;
    double m_mean;
    Eigen::MatrixXd m_probabilities;
};

}

#endif
