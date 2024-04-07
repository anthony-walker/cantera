//! @file Wall.h Header file for base class WallBase.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_WALL_H
#define CT_WALL_H

#include "cantera/base/ctexceptions.h"
#include "cantera/zeroD/ReactorBase.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera
{

class Func1;

/**
 * Base class for 'walls' (walls, pistons, etc.) connecting reactors.
 * @ingroup wallGroup
 */
class WallBase
{
public:
    WallBase() = default;

    virtual ~WallBase() {}
    WallBase(const WallBase&) = delete;
    WallBase& operator=(const WallBase&) = delete;

    //! String indicating the wall model implemented. Usually
    //! corresponds to the name of the derived class.
    virtual string type() const {
        return "WallBase";
    }

    //! Rate of volume change (m^3/s) for the adjacent reactors.
    /*!
     * This method is called by Reactor::evalWalls(). Base class method
     * does nothing (that is, constant volume), but may be overloaded.
     * @deprecated Still used by traditional MATLAB toolbox; replaceable by expansionRate.
     */
    virtual double vdot(double t) {
        warn_deprecated("WallBase::vdot",
            "To be removed; replaceable by 'expansionRate'.");
        return 0.0;
    }

    //! Rate of volume change (m^3/s) for the adjacent reactors at current reactor
    //! network time.
    /*!
     * This method is called by Reactor::evalWalls(). Base class method
     * does nothing (that is, constant volume), but may be overloaded.
     * @since New in %Cantera 3.0.
     */
    virtual double expansionRate() {
        return 0.0;
    }

    //! Heat flow rate through the wall (W).
    /*!
     * This method is called by Reactor::evalWalls(). Base class method
     * does nothing (that is, an adiabatic wall), but may be overloaded.
     * @deprecated Still used by traditional MATLAB toolbox; replaceable by heatRate.
     */
    virtual double Q(double t) {
        warn_deprecated("WallBase::Q",
            "To be removed; replaceable by 'heatRate'.");
        return 0.0;
    }

    //! Heat flow rate through the wall (W) at current reactor network time.
    /*!
     * This method is called by Reactor::evalWalls(). Base class method
     * does nothing (that is, an adiabatic wall), but may be overloaded.
     * @since New in %Cantera 3.0.
     */
    virtual double heatRate() {
        return 0.0;
    }

    //! Area in (m^2).
    double area() {
        return m_area;
    }

    //! Set the area [m^2].
    virtual void setArea(double a);

    //! Install the wall between two reactors or reservoirs
    bool install(ReactorBase& leftReactor, ReactorBase& rightReactor);

    //! Called just before the start of integration
    virtual void initialize() {}

    //! True if the wall is correctly configured and ready to use.
    virtual bool ready() {
        return (m_left != 0 && m_right != 0);
    }

    //! Return a reference to the Reactor or Reservoir to the left of the wall.
    ReactorBase& left() const {
        return *m_left;
    }

    //! Return a reference to the Reactor or Reservoir to the right of the wall.
    const ReactorBase& right() {
        return *m_right;
    }

    //! Set current reactor network time
    /*!
     * @since New in %Cantera 3.0.
     */
    void setSimTime(double time) {
        m_time = time;
    }

    //! Build the Jacobian terms specific to the flow device for the given connected
    //! reactor.
    //! @param r a pointer to the calling reactor
    //! @param jacVector a vector of triplets to be added to the Jacobian for the
    //! reactor
    //! @warning This function is an experimental part of the %Cantera API and may be
    //! changed or removed without notice.
    //! @since New in %Cantera 3.1.
    //!
    virtual void buildReactorJacobian(ReactorBase* r, vector<Eigen::Triplet<double>>& jacVector) {
        throw NotImplementedError("WallBase::buildReactorJacobian");
    }

    //! Build the Jacobian terms specific to the flow device for the network. These
    //! terms
    //! will be adjusted to the networks indexing system outside of the reactor.
    //! @param jacVector a vector of triplets to be added to the Jacobian for the
    //! reactor
    //! @warning This function is an experimental part of the %Cantera API and may be
    //! changed or removed without notice.
    //! @since New in %Cantera 3.1.
    //!
    virtual void buildNetworkJacobian(vector<Eigen::Triplet<double>>& jacVector) {
        throw NotImplementedError("WallBase::buildNetworkJacobian");
    }

    //! Specify the Jacobian terms have been calculated and should not be recalculated.
    //! @warning This function is an experimental part of the %Cantera API and may be
    //! changed or removed without notice.
    //! @since New in %Cantera 3.1.
    //!
    void calculatedJacobian() { m_jac_calculated = true; };

    //! Specify that Jacobian terms have not been calculated and should be recalculated.
    //! @warning This function is an experimental part of the %Cantera API and may be
    //! changed or removed without notice.
    //! @since New in %Cantera 3.1.
    //!
    void notCalculatedJacobian() { m_jac_calculated = false; };

protected:
    ReactorBase* m_left = nullptr;
    ReactorBase* m_right = nullptr;

    //! current reactor network time
    double m_time = 0.0;

    double m_area = 1.0;

    //! a variable to switch on and off so calculations are not doubled by the calling
    //! reactor or network
    bool m_jac_calculated = false;
};

//! Represents a wall between between two ReactorBase objects.
/*!
 * Walls can move (changing the volume of the adjacent reactors) and allow heat
 * transfer between reactors.
 * @ingroup wallGroup
 */
class Wall : public WallBase
{
public:
    Wall() = default;

    //! String indicating the wall model implemented. Usually
    //! corresponds to the name of the derived class.
    string type() const override {
        return "Wall";
    }

    //! Wall velocity @f$ v(t) @f$ at current reactor network time.
    //! @since New in %Cantera 3.0.
    double velocity() const;

    //! Set the wall velocity to a specified function of time, @f$ v(t) @f$.
    void setVelocity(Func1* f=0) {
        if (f) {
            m_vf = f;
        }
    }

    //! Rate of volume change (m^3/s) for the adjacent reactors.
    /*!
     * The volume rate of change is given by
     * @f[
     *     \dot V = K A (P_{left} - P_{right}) + F(t)
     * @f]
     * where *K* is the specified expansion rate coefficient, *A* is the wall
     * area, and *F(t)* is a specified function of time. Positive values for
     * `vdot` correspond to increases in the volume of reactor on left, and
     * decreases in the volume of the reactor on the right.
     * @deprecated Still used by traditional MATLAB toolbox; replaceable by
     *      expansionRate.
     */
    double vdot(double t) override;

    //! Rate of volume change (m^3/s) for the adjacent reactors.
    /*!
     * The volume rate of change is given by
     * @f[
     *     \dot V = K A (P_{left} - P_{right}) + F(t)
     * @f]
     * where *K* is the specified expansion rate coefficient, *A* is the wall area,
     * and and *F(t)* is a specified function evaluated at the current network time.
     * Positive values for `expansionRate` correspond to increases in the volume of
     * reactor on left, and decreases in the volume of the reactor on the right.
     * @since New in %Cantera 3.0.
     */
    double expansionRate() override;

    //! Heat flux function @f$ q_0(t) @f$ evaluated at current reactor network time.
    //! @since New in %Cantera 3.0.
    double heatFlux() const;

    //! Specify the heat flux function @f$ q_0(t) @f$.
    void setHeatFlux(Func1* q) {
        m_qf = q;
    }

    //! Heat flow rate through the wall (W).
    /*!
     * The heat flux is given by
     * @f[
     *     Q = h A (T_{left} - T_{right}) + A G(t)
     * @f]
     * where *h* is the heat transfer coefficient, *A* is the wall area, and
     * *G(t)* is a specified function of time. Positive values denote a flux
     * from left to right.
     * @deprecated Still used by traditional MATLAB toolbox; replaceable by heatRate.
     */
    double Q(double t) override;

    //! Heat flow rate through the wall (W).
    /*!
     * The heat flux is given by
     * @f[
     *     Q = h A (T_{left} - T_{right}) + A G(t)
     * @f]
     * where *h* is the heat transfer coefficient, *A* is the wall area, and
     * *G(t)* is a specified function of time evaluated at the current network
     * time. Positive values denote a flux from left to right.
     * @since New in %Cantera 3.0.
     */
    double heatRate() override;

    void setThermalResistance(double Rth) {
        m_rrth = 1.0/Rth;
    }

    //! Set the overall heat transfer coefficient [W/m^2/K].
    void setHeatTransferCoeff(double U) {
        m_rrth = U;
    }

    //! Get the overall heat transfer coefficient [W/m^2/K].
    double getHeatTransferCoeff() const {
        return m_rrth;
    }

    //! Set the emissivity.
    void setEmissivity(double epsilon) {
        if (epsilon > 1.0 || epsilon < 0.0) {
            throw CanteraError("WallBase::setEmissivity",
                               "emissivity must be between 0.0 and 1.0");
        }
        m_emiss = epsilon;
    }

    //! Get the emissivity.
    double getEmissivity() const {
        return m_emiss;
    }

    //! Set the expansion rate coefficient.
    void setExpansionRateCoeff(double k) {
        m_k = k;
    }

    //! Get the expansion rate coefficient
    double getExpansionRateCoeff() const {
        return m_k;
    }

    virtual void buildReactorJacobian(ReactorBase* r, vector<Eigen::Triplet<double>>& jacVector) override;

    virtual void buildNetworkJacobian(vector<Eigen::Triplet<double>>& jacVector) override;

protected:

    //! expansion rate coefficient
    double m_k = 0.0;

    //! heat transfer coefficient
    double m_rrth = 0.0;

    //! emissivity
    double m_emiss = 0.0;

    //! Velocity function
    Func1* m_vf = nullptr;

    //! Heat flux function
    Func1* m_qf = nullptr;
};

}

#endif
