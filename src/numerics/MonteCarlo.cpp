//! @file MonteCarlo.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/global.h"
#include "cantera/numerics/MonteCarlo.h"

namespace Cantera
{

    void MonteCarlo::initialize(size_t nsp, size_t bins) {
        // set initial variables
        m_bins = bins;
        m_nsp = nsp;
        // resize matrix

    }


}
