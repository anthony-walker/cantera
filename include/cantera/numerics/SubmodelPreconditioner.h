/**
 *  @file SubmodelPreconditioner.h Declarations for the class
 *   SubmodelPreconditioner which is a child class of PreconditionerBase
 *   for preconditioners used by sundials
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef SUBMODELPRECONDITIONER_H
#define SUBMODELPRECONDITIONER_H

#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/base/global.h"
#include <iostream>
#include "cantera/base/AnyMap.h"

namespace Cantera
{


// Forward Declaration of Reactor
class Reactor;

//! SubmodelPreconditioner a preconditioner designed for use with large
//! mechanisms that leverages sparse solvers. It does this by pruning
//! the preconditioner by a threshold value. It also neglects pressure
//! dependence and thirdbody contributions in its formation and has a
//! finite difference approximation for temperature.
class SubmodelPreconditioner : public AdaptivePreconditioner
{
public:
    SubmodelPreconditioner() = default;

    string type() override { return "SubmodelPreconditioner"; }

    void initialize(size_t networkSize) override;

    // void stateAdjustment(vector_fp& state) override;

    // void setup() override;

    void prunePreconditioner() override;

    void addReactor(Reactor* reactor) { m_reactors.push_back(reactor); }

protected:
    //! map to translate from detailed model to submodel Jacobian
    std::map<size_t, bool> m_submodel_map;

    //! Vector of reactors
    std::vector<Reactor*> m_reactors;
};

}

#endif
