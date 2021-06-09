//! @file PreconditionerBase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/PreconditionerBase.h"

namespace Cantera
{
    double PreconditionerBase::getThreshold()
    {
        return this->m_threshold;
    }

    void PreconditionerBase::setThreshold(double threshold)
    {
        this->m_threshold = threshold;
    }

    void PreconditionerBase::setElementByThreshold(size_t row,size_t col, double element)
    {
        if (this->m_threshold<element)
        {
            this->setElement(row,col,element);
        }
    }

    void PreconditionerBase::setDimensions(size_t nrows,size_t ncols)
    {
        this->m_dimensions[0] = nrows;
        this->m_dimensions[1] = ncols;
    }

    size_t* PreconditionerBase::getDimensions()
    {
        return &(this->m_dimensions[0]);
    }
}
