/**
 *  @file StateDiagonalPreconditioner.h Declarations for the class
 *   StateDiagonalPreconditioner which is a child class of PreconditionerBase
 *   for preconditioners used by sundials
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef STATEDIAGONALPRECONDITIONER_H
#define STATEDIAGONALPRECONDITIONER_H

#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/base/global.h"
#include <iostream>
#include "cantera/base/AnyMap.h"

namespace Cantera
{

// Forward Declaration of Reactor
class Reactor;

//! StateDiagonalPreconditioner a preconditioner designed for use with large
//! mechanisms that leverages sparse solvers. It does this by pruning
//! the preconditioner by a threshold value. It also neglects pressure
//! dependence and thirdbody contributions in its formation and has a
//! finite difference approximation for temperature.
class StateDiagonalPreconditioner : public AdaptivePreconditioner
{
public:
    StateDiagonalPreconditioner() = default;

    string type() override { return "SubmodelPreconditioner"; }

    void initialize(size_t networkSize) override;

    virtual void stateAdjustment(vector_fp& state) override {};

    void setup() override;

    void addReactor(Reactor* reactor) { m_reactors.push_back(reactor); };

protected:
    //! Vector of reactors
    std::vector<Reactor*> m_reactors;
};

}

#endif
