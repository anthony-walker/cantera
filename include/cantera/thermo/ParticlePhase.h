/**
 *  @file ParticlePhase.h
 * Header file for class ParticlePhase, the base class for phases with
 * thermodynamic properties, and the text for the Module thermoprops
 * (see \ref thermoprops and class \link Cantera::ParticlePhase ParticlePhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PARTICLEPHASE_H
#define CT_PARTICLEPHASE_H

#include "cantera/base/Units.h"
#include "cantera/base/AnyMap.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

//! Base class for a phase with thermodynamic properties.
/*!
 * Class ParticlePhase is a class for phase representation
 * @ingroup thermoprops
 * @ingroup phases
 */
class ParticlePhase : public ThermoPhase
{
public:
    ParticlePhase() = default;

    //! @name  Information Methods
    //! @{
    virtual std::string type() const {
        return "ParticlePhase";
    }
    //! @}

    //! The number of bins used in a MonteCarlo simulation
    size_t bins() {
        return m_nbins;
    }

    //! Set the number of bins used in a MonteCarlo simulation
    void setBins(size_t nbins) {
        m_nbins = nbins;
    }

protected:
    size_t m_nbins;

    };
}

#endif
