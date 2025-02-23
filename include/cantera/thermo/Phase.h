/**
 * @file Phase.h
 * Header file for class Phase.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PHASE_H
#define CT_PHASE_H

#include "cantera/base/ctexceptions.h"
#include "cantera/thermo/Elements.h"
#include "cantera/base/ValueCache.h"

namespace Cantera
{

class Solution;
class Species;
class Kinetics;

//! Class Phase is the base class for phases of matter, managing the species and
//! elements in a phase, as well as the independent variables of temperature,
//! mass density (compressible substances) or pressure (incompressible
//! substances), species mass/mole fraction, and other generalized forces and
//! intrinsic properties (such as electric potential) that define the
//! thermodynamic state.
/*!
 *
 * Class Phase provides information about the elements and species in a
 * phase - names, index numbers (location in arrays), atomic or molecular
 * weights, etc. The set of elements must include all those that compose the
 * species, but may include additional elements.
 *
 * It also stores an array of species molecular weights, which are used to
 * convert between mole and mass representations of the composition. For
 * efficiency in mass/mole conversion, the vector of mass fractions divided
 * by molecular weight @f$ Y_k/M_k @f$ is also stored.
 *
 * Class Phase is not usually used directly. Its primary use is as a base class
 * for class ThermoPhase. It is not generally necessary to overloaded any of
 * class Phase's methods, which handles both compressible and incompressible
 * phases. For incompressible phases, the density is replaced by the pressure
 * as the independent variable, and can no longer be set directly. In this case,
 * the density needs to be calculated from a suitable equation of state, and
 * assigned to the object using the assignDensity() method. This also applies
 * for nearly-incompressible phases or phases which utilize standard states
 * based on a T and P, in which case they need to overload these functions too.
 *
 * Class Phase contains a number of utility functions that will set the state
 * of the phase in its entirety, by first setting the composition, and then
 * temperature and pressure. An example of this is the function
 * Phase::setState_TPY(double t, double p, const double* y).
 *
 * For bulk (3-dimensional) phases, the mass density has units of kg/m^3, and the molar
 * density and concentrations have units of kmol/m^3, and the units listed in the
 * methods of the Phase class assume a bulk phase. However, for surface (2-dimensional)
 * phases have units of kg/m^2 and kmol/m^2, respectively. And for edge (1-dimensional)
 * phases, these units kg/m and kmol/m.
 *
 * Class Phase contains methods for saving and restoring the full internal state
 * of a given phase. These are saveState() and restoreState(). These functions
 * operate on a state vector, which by default uses the first two entries for
 * temperature and density (compressible substances) or temperature and
 * pressure (incompressible substances). If the substance is not pure in a
 * thermodynamic sense (that is, it may contain multiple species), the state also
 * contains nSpecies() entries that specify the composition by corresponding
 * mass fractions. Default definitions can be overloaded by derived classes.
 * For any phase, the native definition of its thermodynamic state is defined
 * the method nativeState(), with the length of the state vector returned by
 * by stateSize(). In addition, methods isPure() and isCompressible() provide
 * information on the implementation of a Phase object.
 *
 * A species name is referred to via speciesName(), which is unique within a
 * given phase. Note that within multiphase mixtures (MultiPhase()), both a
 * phase name/index as well as species name are required to access information
 * about a species in a particular phase. For surfaces, the species names are
 * unique among the phases.
 *
 * @todo
 *   - Specify that the input mole, mass, and volume fraction vectors must sum
 *     to one on entry to the set state routines. Non-conforming mole/mass
 *     fraction vectors are not thermodynamically consistent. Moreover, unless
 *     we do this, the calculation of Jacobians will be altered whenever the
 *     treatment of non- conforming mole fractions is changed. Add setState
 *     functions corresponding to specifying mole numbers, which is actually
 *     what is being done (well one of the options, there are many) when non-
 *     conforming mole fractions are input. Note, we realize that most numerical
 *     Jacobian and some analytical Jacobians use non-conforming calculations.
 *     These can easily be changed to the set mole number setState functions.
 *
 * @ingroup thermoprops
 */
class Phase
{
public:
    Phase() = default; //!< Default constructor.
    virtual ~Phase() = default;

    // Phase objects are not copyable or assignable
    Phase(const Phase&) = delete;
    Phase& operator=(const Phase&) = delete;

    /**
     * @name Name
     * Class Phase uses the string name to identify a phase. For phases instantiated
     * from YAML input files, the name is the value of the corresponding key in the
     * phase map.
     *
     * However, the name field may be changed to another value during the
     * course of a calculation. For example, if duplicates of a phase object
     * are instantiated and used in multiple places (such as a ReactorNet), they
     * will have the same constitutive input, that is, the names of the phases will
     * be the same. Note that this is not a problem for %Cantera internally;
     * however, a user may want to rename phase objects in order to clarify.
     */
    //!@{

    //! Return the name of the phase.
    /*!
     *   Names are unique within a %Cantera problem.
     */
    string name() const;

    //! Sets the string name for the phase.
    //!     @param nm String name of the phase
    void setName(const string& nm);

    //! String indicating the thermodynamic model implemented. Usually
    //! corresponds to the name of the derived class, less any suffixes such as
    //! "Phase", TP", "VPSS", etc.
    //! @since Starting in %Cantera 3.0, the name returned by this method corresponds
    //!     to the canonical name used in the YAML input format.
    virtual string type() const {
        return "Phase";
    }

    //! @} end group Name

    //! @name Element and Species Information
    //! @{

    //! Name of the element with index m.
    //!     @param m  Element index.
    string elementName(size_t m) const;

    //! Return the index of element named 'name'. The index is an integer
    //! assigned to each element in the order it was added. Returns @ref npos
    //! if the specified element is not found.
    //!     @param name Name of the element
    size_t elementIndex(const string& name) const;

    //! Return a read-only reference to the vector of element names.
    const vector<string>& elementNames() const;

    //! Atomic weight of element m.
    //!     @param m  Element index
    double atomicWeight(size_t m) const;

    //! Entropy of the element in its standard state at 298 K and 1 bar.
    //! If no entropy value was provided when the phase was constructed,
    //! returns the value `ENTROPY298_UNKNOWN`.
    //!     @param m  Element index
    double entropyElement298(size_t m) const;

    //! Atomic number of element m.
    //!     @param m Element index
    int atomicNumber(size_t m) const;

    //! Return the element constraint type
    //! Possible types include:
    //!
    //!   - `CT_ELEM_TYPE_TURNEDOFF        -1`
    //!   - `CT_ELEM_TYPE_ABSPOS            0`
    //!   - `CT_ELEM_TYPE_ELECTRONCHARGE    1`
    //!   - `CT_ELEM_TYPE_CHARGENEUTRALITY  2`
    //!   - `CT_ELEM_TYPE_LATTICERATIO      3`
    //!   - `CT_ELEM_TYPE_KINETICFROZEN     4`
    //!   - `CT_ELEM_TYPE_SURFACECONSTRAINT 5`
    //!   - `CT_ELEM_TYPE_OTHERCONSTRAINT   6`
    //!
    //! The default is `CT_ELEM_TYPE_ABSPOS`.
    //!     @param m  Element index
    //!     @returns the element type
    int elementType(size_t m) const;

    //! Change the element type of the mth constraint
    //! Reassigns an element type.
    //!     @param m  Element index
    //!     @param elem_type New elem type to be assigned
    //!     @returns the old element type
    int changeElementType(int m, int elem_type);

    //! Return a read-only reference to the vector of atomic weights.
    const vector<double>& atomicWeights() const;

    //! Number of elements.
    size_t nElements() const;

    //! Check that the specified element index is in range.
    //! Throws an exception if m is greater than nElements()-1
    void checkElementIndex(size_t m) const;

    //! Check that an array size is at least nElements().
    //! Throws an exception if mm is less than nElements(). Used before calls
    //! which take an array pointer.
    void checkElementArraySize(size_t mm) const;

    //! Number of atoms of element @c m in species @c k.
    //!     @param k    species index
    //!     @param m    element index
    double nAtoms(size_t k, size_t m) const;

    //! Returns the index of a species named 'name' within the Phase object.
    //! The first species in the phase will have an index 0, and the last one
    //! will have an index of nSpecies() - 1.
    //!     @param name String name of the species. It may also be in the form
    //!            phaseName:speciesName
    //!     @return The index of the species. If the name is not found,
    //!             the value @ref npos is returned.
    size_t speciesIndex(const string& name) const;

    //! Name of the species with index k
    //!     @param k index of the species
    string speciesName(size_t k) const;

    //! Return a const reference to the vector of species names
    const vector<string>& speciesNames() const;

    //! Returns the number of species in the phase
    size_t nSpecies() const {
        return m_kk;
    }

    //! Check that the specified species index is in range.
    //! Throws an exception if k is greater than nSpecies()-1
    void checkSpeciesIndex(size_t k) const;

    //! Check that an array size is at least nSpecies().
    //! Throws an exception if kk is less than nSpecies(). Used before calls
    //! which take an array pointer.
    void checkSpeciesArraySize(size_t kk) const;

    //! @} end group Element and Species Information

    //! Return whether phase represents a pure (single species) substance
    virtual bool isPure() const {
        return false;
    }

    //! Return whether phase represents a substance with phase transitions
    virtual bool hasPhaseTransition() const {
        return false;
    }

    //! Return whether phase represents a compressible substance
    virtual bool isCompressible() const {
        return true;
    }

    //! Return a map of properties defining the native state of a substance.
    //! By default, entries include "T", "D", "Y" for a compressible substance
    //! and "T", "P", "Y" for an incompressible substance, with offsets 0, 1 and
    //! 2, respectively. Mass fractions "Y" are omitted for pure species.
    //! In all cases, offsets into the state vector are used by saveState()
    //! and restoreState().
    virtual map<string, size_t> nativeState() const;

    //! Return string acronym representing the native state of a Phase.
    //! Examples: "TP", "TDY", "TPY".
    //! @see nativeState
    //! @since New in %Cantera 3.0
    string nativeMode() const;

    //! Return a vector containing full states defining a phase.
    //! Full states list combinations of properties that allow for the
    //! specification of a thermodynamic state based on user input.
    //! Properties and states are represented by single letter acronyms, and
    //! combinations of letters, respectively (for example, "TDY", "TPX", "SVX").
    //! Supported property acronyms are:
    //!    "T": temperature
    //!    "P": pressure
    //!    "D": density
    //!    "X": mole fractions
    //!    "Y": mass fractions
    //!    "T": temperature
    //!    "U": specific internal energy
    //!    "V": specific volume
    //!    "H": specific enthalpy
    //!    "S": specific entropy
    //!    "Q": vapor fraction
    virtual vector<string> fullStates() const;

    //! Return a vector of settable partial property sets within a phase.
    //! Partial states encompass all valid combinations of properties that allow
    //! for the specification of a state while ignoring species concentrations
    //! (such as "TD", "TP", "SV").
    virtual vector<string> partialStates() const;

    //! Return size of vector defining internal state of the phase.
    //! Used by saveState() and restoreState().
    virtual size_t stateSize() const;

    //! Save the current internal state of the phase.
    //! Write to vector 'state' the current internal state.
    //!     @param state output vector. Will be resized to stateSize().
    void saveState(vector<double>& state) const;

    //! Write to array 'state' the current internal state.
    //!     @param lenstate length of the state array. Must be >= stateSize()
    //!     @param state    output vector. Must be of length stateSizes() or
    //!                     greater.
    virtual void saveState(size_t lenstate, double* state) const;

    //! Restore a state saved on a previous call to saveState.
    //!     @param state State vector containing the previously saved state.
    void restoreState(const vector<double>& state);

    //! Restore the state of the phase from a previously saved state vector.
    //!     @param lenstate   Length of the state vector
    //!     @param state      Vector of state conditions.
    virtual void restoreState(size_t lenstate, const double* state);

    //! @name Set Thermodynamic State
    //!
    //! Set the internal thermodynamic state by setting the internally stored
    //! temperature, density and species composition. Note that the composition
    //! is always set first.
    //!
    //! Temperature and density are held constant if not explicitly set.
    //! @{

    //! Set the species mole fractions by name.
    //! Species not listed by name in @c xMap are set to zero.
    //!     @param xMap map from species names to mole fraction values.
    void setMoleFractionsByName(const Composition& xMap);

    //! Set the mole fractions of a group of species by name. Species which
    //! are not listed by name in the composition map are set to zero.
    //!     @param x string x in the form of a composition map
    void setMoleFractionsByName(const string& x);

    //! Set the species mass fractions by name.
    //! Species not listed by name in @c yMap are set to zero.
    //!     @param yMap map from species names to mass fraction values.
    void setMassFractionsByName(const Composition& yMap);

    //! Set the species mass fractions by name.
    //! Species not listed by name in @c x are set to zero.
    //!     @param x String containing a composition map
    void setMassFractionsByName(const string& x);

    //! Set the internally stored temperature (K) and density (kg/m^3)
    //!     @param t     Temperature in kelvin
    //!     @param rho   Density (kg/m^3)
    //!     @since New in %Cantera 3.0.
    void setState_TD(double t, double rho);

    //! @} end group set thermo state

    //! Molecular weight of species @c k.
    //!     @param k   index of species @c k
    //!     @returns the molecular weight of species @c k.
    double molecularWeight(size_t k) const;

    //! Copy the vector of molecular weights into array weights.
    //!     @param weights  Output array of molecular weights (kg/kmol)
    void getMolecularWeights(double* weights) const;

    //! Return a const reference to the internal vector of molecular weights.
    //! units = kg / kmol
    const vector<double>& molecularWeights() const;

    //! Return a const reference to the internal vector of molecular weights.
    //! units = kmol / kg
    const vector<double>& inverseMolecularWeights() const;

    //! Copy the vector of species charges into array charges.
    //!     @param charges Output array of species charges (elem. charge)
    void getCharges(double* charges) const;

    //! @name Composition
    //! @{

    //! Get the mole fractions by name.
    //!     @param threshold   Exclude species with mole fractions less than or
    //!                        equal to this threshold.
    //!     @return Map of species names to mole fractions
    Composition getMoleFractionsByName(double threshold=0.0) const;

    //! Return the mole fraction of a single species
    //!     @param  k  species index
    //!     @return Mole fraction of the species
    double moleFraction(size_t k) const;

    //! Return the mole fraction of a single species
    //!     @param  name  String name of the species
    //!     @return Mole fraction of the species
    double moleFraction(const string& name) const;

    //! Get the mass fractions by name.
    //!     @param threshold   Exclude species with mass fractions less than or
    //!                        equal to this threshold.
    //!     @return Map of species names to mass fractions
    Composition getMassFractionsByName(double threshold=0.0) const;

    //! Return the mass fraction of a single species
    //!     @param  k species index
    //!     @return Mass fraction of the species
    double massFraction(size_t k) const;

    //! Return the mass fraction of a single species
    //!     @param  name  String name of the species
    //!     @return Mass Fraction of the species
    double massFraction(const string& name) const;

    //! Get the species mole fraction vector.
    //!     @param x On return, x contains the mole fractions. Must have a
    //!          length greater than or equal to the number of species.
    void getMoleFractions(double* const x) const;

    //! Set the mole fractions to the specified values.
    //! There is no restriction on the sum of the mole fraction vector.
    //! Internally, the Phase object will normalize this vector before storing
    //! its contents.
    //!     @param x Array of unnormalized mole fraction values (input). Must
    //! have a length greater than or equal to the number of species, m_kk.
    virtual void setMoleFractions(const double* const x);

    //! Set the mole fractions to the specified values without normalizing.
    //! This is useful when the normalization condition is being handled by
    //! some other means, for example by a constraint equation as part of a
    //! larger set of equations.
    //!     @param x  Input vector of mole fractions. Length is m_kk.
    virtual void setMoleFractions_NoNorm(const double* const x);

    //! Get the species mass fractions.
    //!     @param[out] y Array of mass fractions, length nSpecies()
    void getMassFractions(double* const y) const;

    //! Return a const pointer to the mass fraction array
    const double* massFractions() const {
        return &m_y[0];
    }

    //! Set the mass fractions to the specified values and normalize them.
    //!     @param[in] y Array of unnormalized mass fraction values. Length
    //!                  must be greater than or equal to the number of
    //!                  species. The Phase object will normalize this vector
    //!                  before storing its contents.
    virtual void setMassFractions(const double* const y);

    //! Set the mass fractions to the specified values without normalizing.
    //! This is useful when the normalization condition is being handled by
    //! some other means, for example by a constraint equation as part of a
    //! larger set of equations.
    //!     @param y  Input vector of mass fractions. Length is m_kk.
    virtual void setMassFractions_NoNorm(const double* const y);

    //! Get the species concentrations (kmol/m^3).
    /*!
     *    @param[out] c The vector of species concentrations. Units are
     *                  kmol/m^3. The length of the vector must be greater than
     *                  or equal to the number of species within the phase.
     */
    virtual void getConcentrations(double* const c) const;

    //! Concentration of species k.
    //! If k is outside the valid range, an exception will be thrown.
    /*!
     *    @param[in] k Index of the species within the phase.
     *
     *    @returns the concentration of species k (kmol m-3).
     */
    virtual double concentration(const size_t k) const;

    //! Set the concentrations to the specified values within the phase.
    //! We set the concentrations here and therefore we set the overall density
    //! of the phase. We hold the temperature constant during this operation.
    //! Therefore, we have possibly changed the pressure of the phase by
    //! calling this routine.
    //!     @param[in] conc Array of concentrations in dimensional units. For
    //!                     bulk phases c[k] is the concentration of the kth
    //!                     species in kmol/m3. For surface phases, c[k] is the
    //!                     concentration in kmol/m2. The length of the vector
    //!                     is the number of species in the phase.
    virtual void setConcentrations(const double* const conc);

    //! Set the concentrations without ignoring negative concentrations
    virtual void setConcentrationsNoNorm(const double* const conc);
    //! @}

    //! Set the state of the object with moles in [kmol]
    virtual void setMolesNoTruncate(const double* const N);

    //! Elemental mass fraction of element m
    /*!
     *  The elemental mass fraction @f$ Z_{\mathrm{mass},m} @f$ of element @f$ m @f$
     *  is defined as
     *  @f[
     *      Z_{\mathrm{mass},m} = \sum_k \frac{a_{m,k} M_m}{M_k} Y_k
     *  @f]
     *  with @f$ a_{m,k} @f$ being the number of atoms of element @f$ m @f$ in
     *  species @f$ k @f$, @f$ M_m @f$ the atomic weight of element @f$ m @f$,
     *  @f$ M_k @f$ the molecular weight of species @f$ k @f$, and @f$ Y_k @f$
     *  the mass fraction of species @f$ k @f$.
     *
     *  @param[in] m Index of the element within the phase. If m is outside
     *               the valid range, an exception will be thrown.
     *
     *  @return the elemental mass fraction of element m.
     */
    double elementalMassFraction(const size_t m) const;

    //! Elemental mole fraction of element m
    /*!
     *  The elemental mole fraction @f$ Z_{\mathrm{mole},m} @f$ of element @f$ m @f$
     *  is the number of atoms of element *m* divided by the total number of
     *  atoms. It is defined as:
     *
     *  @f[
     *      Z_{\mathrm{mole},m} = \frac{\sum_k a_{m,k} X_k}
     *                                 {\sum_k \sum_j a_{j,k} X_k}
     *  @f]
     *  with @f$ a_{m,k} @f$ being the number of atoms of element @f$ m @f$ in
     *  species @f$ k @f$, @f$ \sum_j @f$ being a sum over all elements, and
     *  @f$ X_k @f$ being the mole fraction of species @f$ k @f$.
     *
     *  @param[in] m Index of the element within the phase. If m is outside the
     *               valid range, an exception will be thrown.
     *  @return the elemental mole fraction of element m.
     */
    double elementalMoleFraction(const size_t m) const;

    //! Dimensionless electrical charge of a single molecule of species k
    //! The charge is normalized by the the magnitude of the electron charge
    //!     @param k species index
    double charge(size_t k) const {
        return m_speciesCharge[k];
    }

    //! Charge density [C/m^3].
    double chargeDensity() const;

    //! Returns the number of spatial dimensions (1, 2, or 3)
    size_t nDim() const {
        return m_ndim;
    }

    //! Set the number of spatial dimensions (1, 2, or 3). The number of
    //! spatial dimensions is used for vector involving directions.
    //!     @param ndim Input number of dimensions.
    void setNDim(size_t ndim) {
        m_ndim = ndim;
    }

    //! @name Thermodynamic Properties
    //! @{

    //! Temperature (K).
    //!     @return The temperature of the phase
    double temperature() const {
        return m_temp;
    }

    //! Electron Temperature (K)
    //!     @return The electron temperature of the phase
    virtual double electronTemperature() const {
        return m_temp;
    }

    //! Return the thermodynamic pressure (Pa).
    /*!
     *  This method must be overloaded in derived classes. Within %Cantera, the
     *  independent variable is either density or pressure. If the state is
     *  defined by temperature, density, and mass fractions, this method should
     *  use these values to implement the mechanical equation of state @f$ P(T,
     *  \rho, Y_1, \dots, Y_K) @f$. Alternatively, it returns a stored value.
     */
    virtual double pressure() const {
        throw NotImplementedError("Phase::pressure",
            "Not implemented for thermo model '{}'", type());
    }

    //! Density (kg/m^3).
    //!     @return The density of the phase
    virtual double density() const {
        return m_dens;
    }

    //! Molar density (kmol/m^3).
    //!     @return The molar density of the phase
    virtual double molarDensity() const;

    //! Molar volume (m^3/kmol).
    //!     @return The molar volume of the phase
    virtual double molarVolume() const;

    //! Set the internally stored density (kg/m^3) of the phase.
    //! Note the density of a phase is an independent variable.
    //!     @param[in] density_ density (kg/m^3).
    virtual void setDensity(const double density_);

    //! Set the internally stored pressure (Pa) at constant temperature and
    //! composition
    /*!
     *  This method must be reimplemented in derived classes, where it may
     *  involve the solution of a nonlinear equation. Within %Cantera, the
     *  independent variable is either density or pressure. Therefore, this
     *  function may either solve for the density that will yield the desired
     *  input pressure or set an independent variable. The temperature
     *  and composition are held constant during this process.
     *
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(double p) {
        throw NotImplementedError("Phase::setPressure",
            "Not implemented for thermo model '{}'", type());
    }

    //! Set the internally stored temperature of the phase (K).
    //!     @param temp Temperature in Kelvin
    virtual void setTemperature(double temp) {
        if (temp > 0) {
            m_temp = temp;
        } else {
            throw CanteraError("Phase::setTemperature",
                               "temperature must be positive. T = {}", temp);
        }
    }

    //! Set the internally stored electron temperature of the phase (K).
    //!     @param etemp Electron temperature in Kelvin
    virtual void setElectronTemperature(double etemp) {
        throw NotImplementedError("Phase::setElectronTemperature",
            "Not implemented for thermo model '{}'", type());
    }

    //! @}

    //! @name Mean Properties
    //! @{

    //! Evaluate the mole-fraction-weighted mean of an array Q.
    //! @f[ \sum_k X_k Q_k. @f]
    //! Q should contain pure-species molar property values.
    //!     @param[in] Q Array of length m_kk that is to be averaged.
    //!     @return mole-fraction-weighted mean of Q
    double mean_X(const double* const Q) const;

    //! @copydoc Phase::mean_X(const double* const Q) const
    double mean_X(const vector<double>& Q) const;

    //!  The mean molecular weight. Units: (kg/kmol)
    double meanMolecularWeight() const {
        return m_mmw;
    }

    //! Evaluate @f$ \sum_k X_k \ln X_k @f$.
    //! @return The indicated sum. Dimensionless.
    double sum_xlogx() const;

    //! @}
    //! @name Adding Elements and Species
    //!
    //! These methods are used to add new elements or species. These are not
    //! usually called by user programs.
    //!
    //! Since species are checked to insure that they are only composed of
    //! declared elements, it is necessary to first add all elements before
    //! adding any species.
    //! @{

    //! Add an element.
    //!     @param symbol Atomic symbol string.
    //!     @param weight Atomic mass in amu.
    //!     @param atomicNumber Atomic number of the element (unitless)
    //!     @param entropy298 Entropy of the element at 298 K and 1 bar in its
    //!         most stable form. The default is the value ENTROPY298_UNKNOWN,
    //!         which is interpreted as an unknown, and if used will cause
    //!         %Cantera to throw an error.
    //!     @param elem_type Specifies the type of the element constraint
    //!         equation. This defaults to CT_ELEM_TYPE_ABSPOS, that is, an element.
    //!     @return index of the element added
    size_t addElement(const string& symbol, double weight=-12345.0,
                      int atomicNumber=0, double entropy298=ENTROPY298_UNKNOWN,
                      int elem_type=CT_ELEM_TYPE_ABSPOS);

    //! Add a Species to this Phase. Returns `true` if the species was
    //! successfully added, or `false` if the species was ignored.
    //!
    //! Derived classes which need to size arrays according to the number of
    //! species should overload this method. The derived class implementation
    //! should call the base class method, and, if this returns `true`
    //! (indicating that the species has been added), adjust their array sizes
    //! accordingly.
    //!
    //! @see ignoreUndefinedElements addUndefinedElements throwUndefinedElements
    virtual bool addSpecies(shared_ptr<Species> spec);

    //! Modify the thermodynamic data associated with a species.
    /*!
     * The species name, elemental composition, and type of thermo
     * parameterization must be unchanged. If there are Kinetics objects that
     * depend on this phase, Kinetics::invalidateCache() should be called on
     * those objects after calling this function.
     */
    virtual void modifySpecies(size_t k, shared_ptr<Species> spec);

    //! Add a species alias (that is, a user-defined alternative species name).
    //! Aliases are case-sensitive.
    //!     @param name original species name
    //!     @param alias alternate name
    void addSpeciesAlias(const string& name, const string& alias);

    //! Lock species list to prevent addition of new species.
    //! Increments a reference counter used to track whether the Phase is being used by
    //! a Reactor, Domain1D, or MultiPhase object, which require the number of species
    //! to remain constant. Should be called in C++ by the object owning the reference.
    void addSpeciesLock() {
        m_nSpeciesLocks++;
    }

    //! Decrement species lock counter.
    //! Should only be called in C++ by the object owning the reference.
    void removeSpeciesLock();

    //! Return a vector with isomers names matching a given composition map
    //!     @param compMap Composition of the species.
    //!     @return A vector of species names for matching species.
    virtual vector<string> findIsomers(const Composition& compMap) const;

    //! Return a vector with isomers names matching a given composition string
    //!     @param comp String containing a composition map
    //!     @return A vector of species names for matching species.
    virtual vector<string> findIsomers(const string& comp) const;

    //! Return the Species object for the named species. Changes to this object
    //! do not affect the ThermoPhase object until the #modifySpecies function
    //! is called.
    shared_ptr<Species> species(const string& name) const;

    //! Return the Species object for species whose index is *k*. Changes to
    //! this object do not affect the ThermoPhase object until the
    //! #modifySpecies function is called.
    shared_ptr<Species> species(size_t k) const;

    //! Set behavior when adding a species containing undefined elements to just
    //! skip the species.
    void ignoreUndefinedElements();

    //! Set behavior when adding a species containing undefined elements to add
    //! those elements to the phase. This is the default behavior.
    void addUndefinedElements();

    //! Set the behavior when adding a species containing undefined elements to
    //! throw an exception.
    void throwUndefinedElements();

    struct UndefElement { enum behavior {
        error, ignore, add
    }; };

    //! @} end group adding species and elements

    //!  Returns a bool indicating whether the object is ready for use
    /*!
     *  @returns true if the object is ready for calculation, false otherwise.
     */
    virtual bool ready() const;

    //! Return the State Mole Fraction Number
    int stateMFNumber() const {
        return m_stateNum;
    }

    //! Invalidate any cached values which are normally updated only when a
    //! change in state is detected
    virtual void invalidateCache();

    //! Returns `true` if case sensitive species names are enforced
    bool caseSensitiveSpecies() const {
        return m_caseSensitiveSpecies;
    }

    //! Set flag that determines whether case sensitive species are enforced
    //! in look-up operations, for example speciesIndex
    void setCaseSensitiveSpecies(bool cflag = true) {
        m_caseSensitiveSpecies = cflag;
    }

    //! Converts a Composition to a vector with entries for each species
    //! Species that are not specified are set to zero in the vector
    /*!
     * @param[in] comp Composition containing the mixture composition
     * @return vector with length m_kk
     */
    vector<double> getCompositionFromMap(const Composition& comp) const;

    //! Converts a mixture composition from mole fractions to mass fractions
    //!     @param[in] Y mixture composition in mass fractions (length m_kk)
    //!     @param[out] X mixture composition in mole fractions (length m_kk)
    void massFractionsToMoleFractions(const double* Y, double* X) const;

    //! Converts a mixture composition from mass fractions to mole fractions
    //!     @param[in] X mixture composition in mole fractions (length m_kk)
    //!     @param[out] Y mixture composition in mass fractions (length m_kk)
    void moleFractionsToMassFractions(const double* X, double* Y) const;

protected:
    //! Ensure that phase is compressible.
    //! An error is raised if the state is incompressible
    //!     @param setter  name of setter (used for exception handling)
    void assertCompressible(const string& setter) const {
        if (!isCompressible()) {
            throw CanteraError("Phase::assertCompressible",
                               "Setter '{}' is not available. Density is not an "
                               "independent \nvariable for "
                               "'{}' ('{}')", setter, name(), type());
        }
    }

    //! Set the internally stored constant density (kg/m^3) of the phase.
    //! Used for incompressible phases where the density is not an independent
    //! variable, that is, density does not affect pressure in state calculations.
    //!     @param[in] density_ density (kg/m^3).
    void assignDensity(const double density_);

    //! Cached for saved calculations within each ThermoPhase.
    /*!
     *   For more information on how to use this, see examples within the source
     *   code and documentation for this within ValueCache class itself.
     */
    mutable ValueCache m_cache;

    //! Set the molecular weight of a single species to a given value.
    //!
    //! Used by phases where the equation of state is defined for a specific
    //! value of the molecular weight which may not exactly correspond to the
    //! value computed from the chemical formula.
    //!     @param k       id of the species
    //!     @param mw      Molecular Weight (kg kmol-1)
    void setMolecularWeight(const int k, const double mw);

    //! Apply changes to the state which are needed after the composition
    //! changes. This function is called after any call to setMassFractions(),
    //! setMoleFractions(), or similar. For phases which need to execute a
    //! callback after any change to the composition, it should be done by
    //! overriding this function rather than overriding all of the composition-
    //! setting functions. Derived class implementations of compositionChanged()
    //! should call the parent class method as well.
    virtual void compositionChanged();

    size_t m_kk = 0; //!< Number of species in the phase.

    //! Dimensionality of the phase. Volumetric phases have dimensionality 3
    //! and surface phases have dimensionality 2.
    size_t m_ndim = 3;

    //! Atomic composition of the species. The number of atoms of element i
    //! in species k is equal to m_speciesComp[k * m_mm + i]
    //! The length of this vector is equal to m_kk * m_mm
    vector<double> m_speciesComp;

    vector<double> m_speciesCharge; //!< Vector of species charges. length m_kk.

    map<string, shared_ptr<Species>> m_species; //!< Map of Species objects

    size_t m_nSpeciesLocks = 0; //!< Reference counter preventing species addition

    //! Flag determining behavior when adding species with an undefined element
    UndefElement::behavior m_undefinedElementBehavior = UndefElement::add;

    //! Flag determining whether case sensitive species names are enforced
    bool m_caseSensitiveSpecies = false;

    //! Vector of size m_kk, used as a temporary holding area.
    mutable vector<double> m_workS;

private:
    //! Find lowercase species name in m_speciesIndices when case sensitive
    //! species names are not enforced and a user specifies a non-canonical
    //! species name. Raise exception if lowercase name is not unique.
    size_t findSpeciesLower(const string& nameStr) const;

    //! Name of the phase.
    //! Initially, this is the name specified in the YAML input file. It may be changed
    //! to another value during the course of a calculation.
    string m_name;

    double m_temp = 0.001; //!< Temperature (K). This is an independent variable

    //! Density (kg m-3). This is an independent variable except in the case
    //! of incompressible phases, where it has to be changed using the
    //! assignDensity() method. For compressible substances, the pressure is
    //! determined from this variable rather than other way round.
    double m_dens = 0.001;

    double m_mmw = 0.0; //!< mean molecular weight of the mixture (kg kmol-1)

    //! m_ym[k] = mole fraction of species k divided by the mean molecular
    //! weight of mixture.
    mutable vector<double> m_ym;

    //! Mass fractions of the species
    /*!
     *   Note, this vector
     *   Length is m_kk
     */
    mutable vector<double> m_y;

    vector<double> m_molwts; //!< species molecular weights (kg kmol-1)

    vector<double> m_rmolwts; //!< inverse of species molecular weights (kmol kg-1)

    //! State Change variable. Whenever the mole fraction vector changes,
    //! this int is incremented.
    int m_stateNum = -1;

    //! Vector of the species names
    vector<string> m_speciesNames;

    //! Map of species names to indices
    map<string, size_t> m_speciesIndices;

    //! Map of lower-case species names to indices
    map<string, size_t> m_speciesLower;

    size_t m_mm = 0; //!< Number of elements.
    vector<double> m_atomicWeights; //!< element atomic weights (kg kmol-1)
    vector<int> m_atomicNumbers; //!< element atomic numbers
    vector<string> m_elementNames; //!< element names
    vector<int> m_elem_type; //!< Vector of element types

    //! Entropy at 298.15 K and 1 bar of stable state pure elements (J kmol-1)
    vector<double> m_entropy298;
};

}

#endif
