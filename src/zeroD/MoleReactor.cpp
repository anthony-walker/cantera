//! @file MoleReactor.cpp A zero-dimensional reactor with a moles as the
//! state

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#include "cantera/zeroD/MoleReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{

void MoleReactor::getSurfaceInitialConditions(double* N)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        auto currPhase = S->thermo();
        currPhase->getConcentrations(N + loc);
        double area = S->area();
        for (size_t i = loc; i < loc + currPhase->nSpecies(); i++)
        {
            N[i] *= area;
        }
        loc += currPhase->nSpecies();
    }
}

void MoleReactor::initialize(double t0)
{
    Reactor::initialize(t0);
    m_nv -= 1; // moles gives the state one fewer variables
}

void MoleReactor::updateSurfaceState(double* N)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        auto currPhase = S->thermo();
        currPhase->setMolesNoNorm(N+loc);
        loc += currPhase->nSpecies();
    }
}

void MoleReactor::evalSurfaces(double* LHS, double* RHS, double* sdot)
{
    if (!m_surfaces.empty())
    {
        warn_user("MoleReactor::evalSurfaces", "evalSurfaces not currently implemented for MoleReactors");
    }
}

size_t MoleReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + m_sidx;
    } else if (nm == "temperature") {
        return 0;
    } else if (nm == "volume") {
        return 1;
    } else {
        return npos;
    }
}

std::string MoleReactor::componentName(size_t k) {
    if (k == 0) {
        return "temperature";
    }
    else if (k == 1){
        return "volume";
    }
    else if (k >= 2 && k < neq()) {
        k -= 1;
        if (k < m_thermo->nSpecies()) {
            return m_thermo->speciesName(k);
        } else {
            k -= m_thermo->nSpecies();
        }
        for (auto& S : m_surfaces) {
            ThermoPhase* th = S->thermo();
            if (k < th->nSpecies()) {
                return th->speciesName(k);
            } else {
                k -= th->nSpecies();
            }
        }
    }
    throw CanteraError("MoleReactor::componentName",
                       "Index is out of bounds.");
}

}
