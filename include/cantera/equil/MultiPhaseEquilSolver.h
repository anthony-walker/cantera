/**
 *  @file  MultiPhaseEquilSolver.h
 *  Interface class for the vcsnonlinear solver
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef MULTIPHASEEQUILSOLVER_H
#define MULTIPHASEEQUILSOLVER_H

#include "MultiPhase.h"
#include "cantera/numerics/sundials_headers.h"
#include "cantera/numerics/SundialsContext.h"

namespace Cantera
{
/**
 * Specifies the side of the system on which the preconditioner is applied. Not all
 * methods are supported by all integrators.
 */
enum class ConstantState {
    TV, //! Temperature and volume are held constant
    TP, //! Temperature and pressure are held constant
    UV, //! Internal energy and volume are held constant
    HP, //! Enthalpy and pressure are held constant
    SP  //! Entropy and pressure are held constant
};

//! %Cantera's Interface to the Multiphase chemical equilibrium solver.
/*!
 * Class MultiPhaseEquilSolver is designed to be used to set a mixture containing
 * one or more phases to a state of chemical equilibrium.
 *
 * Note, as currently constructed, the underlying ThermoPhase objects are shared
 * between the MultiPhase object and this object. Therefore, mix is not a const
 * argument, and the return parameters are contained in underlying ThermoPhase
 * objects.
 *
 * @ingroup equilGroup
 */
class MultiPhaseEquilSolver
{
public:
    //! Constructor for the multiphase equilibrium solver
    /*!
     * This constructor will initialize the object with a MultiPhase object,
     * setting up the internal equilibration problem. Note, as currently
     * constructed, the underlying ThermoPhase objects are shared between the
     * MultiPhase object and this object. Therefore, mix is not a const
     * argument, and the return parameters are contained in underlying
     * ThermoPhase objects.
     *
     * @param mix Object containing the MultiPhase object
     */
    MultiPhaseEquilSolver(MultiPhase* mix);

    virtual ~MultiPhaseEquilSolver() {}

    //! return the number of iterations
    int iterations() const {
        // FIXME: returns 0
        return 0;
    }

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object
    /*!
     */
    int equilibrate(ConstantState cs, double* LHS, double* RHS, void* f_data);

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object using constant T and P
    /*!
     */
    int equilibrate_TP(double* LHS, double* RHS, void* f_data);

    //! Equilibrate the solution using the current element abundances
    //! stored in the MultiPhase object using either constant H and P
    //! or constant U and P.
    /*!
     * Use the vcs algorithm to equilibrate the current multiphase mixture. The
     * pressure of the calculation is taken from the current pressure stored
     * with the MultiPhase object.
     *
     * @param Htarget Value of the total mixture enthalpy or total internal
     *     energy that will be kept constant. Note, this is and must be an
     *     extensive quantity.  units = Joules
     * @param XY      Integer flag indicating what is held constant. Must be
     *     either HP or UP.
     * @param Tlow    Lower limit of the temperature. It's an error condition
     *     if the temperature falls below Tlow.
     * @param Thigh   Upper limit of the temperature. It's an error condition
     *     if the temperature goes higher than Thigh.
     * @param estimateEquil integer indicating whether the solver
     *     should estimate its own initial condition.
     *     - If 0, the initial mole fraction vector in the ThermoPhase object is
     *       used as the initial condition.
     *     - If 1, the initial mole fraction vector is used if the element
     *       abundances are satisfied.
     *     - if -1, the initial mole fraction vector is thrown out, and an
     *       estimate is formulated.
     * @param err     Internal error level
     * @param maxsteps max steps allowed.
     * @param loglevel Determines the amount of printing to the output file.
     */
    int equilibrate_HP(double Htarget, int XY, double Tlow, double Thigh,
                       int estimateEquil = 0,
                       double err = 1.0E-6,
                       int maxsteps = 1000, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances stored in
    //! the MultiPhase object using constant S and P.
    /*!
     * Use the vcs algorithm to equilibrate the current multiphase mixture. The
     * pressure of the calculation is taken from the current pressure stored
     * with the MultiPhase object.
     *
     * @param Starget Value of the total mixture entropy that will be kept
     *     constant. Note, this is and must be an extensive quantity.
     *     units = Joules/K
     * @param Tlow    Lower limit of the temperature. It's an error condition if
     *     the temperature falls below Tlow.
     * @param Thigh   Upper limit of the temperature. It's an error condition if
     *     the temperature goes higher than Thigh.
     * @param estimateEquil integer indicating whether the solver should
     *     estimate its own initial condition.
     *     - If 0, the initial mole fraction vector in the ThermoPhase object is
     *       used as the initial condition.
     *     - If 1, the initial mole fraction vector is used if the element
     *       abundances are satisfied.
     *     - If -1, the initial mole fraction vector is thrown out, and an
     *       estimate is formulated.
     * @param err     Internal error level
     * @param maxsteps max steps allowed.
     * @param loglevel Determines the amount of printing to the output file.
     */
    int equilibrate_SP(double Starget, double Tlow, double Thigh,
                       int estimateEquil = 0,
                        double err = 1.0E-6,
                       int maxsteps = 1000, int loglevel=-99);

    //! Equilibrate the solution using the current element abundances stored
    //! in the MultiPhase object using constant V and constant T, H, U or S.
    /*!
     * Use the vcs algorithm to equilibrate the current multiphase mixture. The
     * pressure of the calculation is taken from the current pressure stored
     * with the MultiPhase object.
     *
     * @param XY      Integer flag indicating what is held constant.
     *     Must be either TV, HV, UV, or SV.
     * @param xtarget Value of the total thermodynamic parameter to be held
     *     constant in addition to V. Note, except for T, this must be an
     *     extensive quantity.  units = Joules/K or Joules
     * @param estimateEquil integer indicating whether the solver should
     *     estimate its own initial condition.
     *     - If 0, the initial mole fraction vector in the ThermoPhase object is
     *       used as the initial condition.
     *     - If 1, the initial mole fraction vector is used if the element
     *       abundances are satisfied.
     *     - if -1, the initial mole fraction vector is thrown out, and an
     *       estimate is formulated.
     * @param err      Internal error level
     * @param maxsteps max steps allowed.
     * @param logLevel Determines the amount of printing to the output file.
     */
    int equilibrate_TV(int XY, double xtarget,
                       int estimateEquil = 0,
                        double err = 1.0E-6,
                       int maxsteps = 1000, int logLevel = -99);

protected:
    //! Vector that takes into account of the current sorting of the species
    /*!
     * The index of m_order is the original k value of the species in the
     * multiphase.  The value of m_order, k_sorted, is the current value of the
     * species index.
     *
     * `m_order[korig] = k_sorted`
     */
    vector<int> m_order;


    // structures for sparse kinsol solution
    void* m_kin_mem = nullptr; //!< Pointer to the KINSOL memory for the problem
    void* m_linsol = nullptr; //!< Sundials linear solver object
    void* m_linsol_matrix = nullptr; //!< matrix used by Sundials
    SundialsContext m_sundials_ctx; //!< SUNContext object for Sundials>=6.0
    SUNMatrix m_sp_formula_mat; //! Sparse Matrix Object for sundials
    vector<double> m_data; //! Vector for data in sparse matrix object
    vector<sunindextype> m_columns; //! Column indices for sparse matrix object
    vector<sunindextype> m_row_starts; //! Starts of rows within m_data
    N_Vector m_state; //! State of system in an n-vector for solver
    N_Vector m_constraints; //! Constraints on system
    //! Pointer to the MultiPhase mixture that will be equilibrated.
    /*!
     *  Equilibrium solutions will be returned via this variable.
     */
    MultiPhase* m_mix;

    //! Vector of indices for species that are included in the calculation. This
    //! is used to exclude pure-phase species with invalid thermo data
    vector<int> m_species;

};

}

#endif
