#include <cstdio>
#include <string>
#include "gtest/gtest.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/zerodim.h"
#include "cantera/numerics/AdaptivePreconditioner.h"

using namespace Cantera;

// This test ensures that prior reactor initialization of a reactor does
// not affect later integration within a network. This example was
// adapted from test_reactor.py::test_equilibrium_HP.
TEST(ZeroDim, test_individual_reactor_initialization)
{
    // Initial conditions
    double T0 = 1100.0;
    double P0 = 10 * OneAtm;
    double tol = 1e-7;
    std::string X0 = "H2:1.0, O2:0.5, AR:8.0";
    // Reactor solution, phase, and kinetics objects
    std::shared_ptr<Solution> sol1 = newSolution("h2o2.yaml");
    sol1->thermo()->setState_TPX(T0, P0, X0);
    // Set up reactor object
    Reactor reactor1;
    reactor1.insert(sol1);
    // Initialize reactor prior to integration to ensure no impact
    reactor1.initialize();
    // Setup reactor network and integrate
    ReactorNet network;
    network.addReactor(reactor1);
    network.advance(1.0);
    // Secondary gas for comparison
    std::shared_ptr<Solution> sol2 = newSolution("h2o2.yaml");
    sol2->thermo()->setState_TPX(T0, P0, X0);
    sol2->thermo()->equilibrate("UV");
    // Secondary reactor for comparison
    Reactor reactor2;
    reactor2.insert(sol2);
    reactor2.initialize(0.0);
    // Get state of reactors
    std::vector<double> state1 (reactor1.neq());
    std::vector<double> state2 (reactor2.neq());
    reactor1.getState(state1.data());
    reactor2.getState(state2.data());
    // Compare the reactors.
    EXPECT_EQ(reactor1.neq(), reactor2.neq());
    for (size_t i = 0; i < reactor1.neq(); i++)
    {
        EXPECT_NEAR(state1[i], state2[i], tol);
    }
}

TEST(ZeroDim, test_get_reaction_from_reactor)
{
    // Setting up solution to insert in reactor
    auto sol = newSolution("methane_twostep.yaml");
    sol->thermo()->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
    // Set up reactor object
    Reactor reactor;
    reactor.insert(sol);
    // Get reactions
    auto kinetics = reactor.getKineticsMgr();
    auto reaction0 = kinetics->getReaction(0);
    auto reaction1 = kinetics->getReaction(1);
    // Compare reaction strings to manual input
    EXPECT_EQ(reaction0->equation(), "CH4 + 1.5 O2 => CO + 2 H2O");
    EXPECT_EQ(reaction1->equation(), "CO + 0.5 O2 <=> CO2");
}

TEST(ZeroDim, test_get_thermo_from_reactor)
{
    // Setting up solution to insert in reactor
    auto sol = newSolution("h2o2.yaml");
    auto thermo1 = sol->thermo();
    // Set up reactor object
    Reactor reactor;
    reactor.setThermoMgr(*thermo1);
    // Get thermo from reactor
    auto thermo2 = reactor.getThermoMgr();
    // Compare two thermo objects
    EXPECT_EQ(thermo1->nSpecies(), thermo2->nSpecies());
    EXPECT_EQ(thermo1->enthalpy_mass(), thermo2->enthalpy_mass());
    EXPECT_EQ(thermo1->intEnergy_mass(), thermo2->intEnergy_mass());
    EXPECT_EQ(thermo1->pressure(), thermo2->pressure());
}

TEST(MoleReactorTestSet, test_preconditioned_ideal_const_vol)
{
    // Reactor 1
    auto sol1 = newSolution("h2o2.yaml");
    sol1->thermo()->setState_TPX(1100.0, 10 * OneAtm, "H2:1.0, O2:0.5, AR:8.0");
    IdealGasMoleReactor r1;
    r1.insert(sol1);
    // Reactor 2
    auto sol2 = newSolution("h2o2.yaml");
    sol2->thermo()->setState_TPX(1100.0, 10 * OneAtm, "H2:1.0, O2:0.5, AR:8.0");
    IdealGasMoleReactor r2;
    r2.insert(sol2);
    // Network 1
    ReactorNet net1;
    net1.addReactor(r1);
    AdaptivePreconditioner precon;
    net1.setIntegratorType(&precon, GMRES);
    net1.initialize();
    // Network 2
    ReactorNet net2;
    net2.addReactor(r2);
    net2.setIntegratorType(GMRES);
    net2.initialize();
    // Advancing
    net1.advance(0.5);
    net2.advance(0.5);
    // Comparison
    std::vector<double> state1(r1.neq(), 0.0);
    std::vector<double> state2(r2.neq(), 0.0);
    for (size_t i = 0; i < r1.neq(); i++)
    {
        EXPECT_NEAR(state1[i], state2[i], net1.atol());
    }
}

TEST(MoleReactorTestSet, test_preconditioned_ideal_const_pressure)
{
    // Reactor 1
    auto sol1 = newSolution("h2o2.yaml");
    sol1->thermo()->setState_TPX(1100.0, 10 * OneAtm, "H2:1.0, O2:0.5, AR:8.0");
    IdealGasConstPressureMoleReactor r1;
    r1.insert(sol1);
    // Reactor 2
    auto sol2 = newSolution("h2o2.yaml");
    sol2->thermo()->setState_TPX(1100.0, 10 * OneAtm, "H2:1.0, O2:0.5, AR:8.0");
    IdealGasConstPressureMoleReactor r2;
    r2.insert(sol2);
    // Network 1
    ReactorNet net1;
    net1.addReactor(r1);
    AdaptivePreconditioner precon;
    net1.setIntegratorType(&precon, GMRES);
    net1.initialize();
    // Network 2
    ReactorNet net2;
    net2.addReactor(r2);
    net2.setIntegratorType(GMRES);
    net2.initialize();
    // Advancing
    net1.advance(0.5);
    net2.advance(0.5);
    // Comparison
    std::vector<double> state1(r1.neq(), 0.0);
    std::vector<double> state2(r2.neq(), 0.0);
    for (size_t i = 0; i < r1.neq(); i++)
    {
        EXPECT_NEAR(state1[i], state2[i], net1.atol());
    }
}

TEST(MoleReactorTestSet, test_preconditioned_ideal_network_integrations)
{
    // Network 1
    auto sol1 = newSolution("h2o2.yaml");
    sol1->thermo()->setState_TPX(1100.0, 10 * OneAtm, "H2:1.0, O2:0.5, AR:8.0");
    IdealGasMoleReactor r1;
    IdealGasConstPressureMoleReactor r2;
    r1.insert(sol1);
    r2.insert(sol1);
    ReactorNet net1;
    net1.addReactor(r1);
    net1.addReactor(r2);
    AdaptivePreconditioner precon;
    net1.setIntegratorType(&precon, GMRES);
    net1.initialize();
    // Network 2
    ReactorNet net2;
    auto sol2 = newSolution("h2o2.yaml");
    sol2->thermo()->setState_TPX(1100.0, 10 * OneAtm, "H2:1.0, O2:0.5, AR:8.0");
    IdealGasMoleReactor r3;
    IdealGasConstPressureMoleReactor r4;
    r3.insert(sol1);
    r4.insert(sol1);
    net2.addReactor(r3);
    net2.addReactor(r4);
    net2.setIntegratorType(GMRES);
    net2.initialize();
    // Advancing
    net1.advance(0.1);
    net2.advance(0.1);
    // Comparison
    std::vector<double> state1(net1.neq(), 0.0);
    std::vector<double> state2(net2.neq(), 0.0);
    for (size_t i = 0; i < net1.neq(); i++)
    {
        EXPECT_NEAR(state1[i], state2[i], net1.atol());
    }
}

TEST(MoleReactorTestSet, test_mole_reactor_get_state)
{
    // Setting up solution object and thermo/kinetics pointers
    auto sol = newSolution("h2o2.yaml");
    sol->thermo()->setState_TPY(1000.0, 101325, "H2:0.5, O2:0.5");
    IdealGasConstPressureMoleReactor reactor;
    reactor.insert(sol);
    reactor.setInitialVolume(0.5);
    reactor.setEnergy(false);
    reactor.initialize();
    std::vector<double> state(reactor.neq());
    std::vector<double> updatedState(reactor.neq());
    // test get state
    auto thermo = reactor.getThermoMgr();
    const std::vector<double>& imw = thermo->inverseMolecularWeights();
    // prescribed state
    double mass = reactor.volume() * thermo->density();
    size_t H2I = reactor.componentIndex("H2")-1;
    size_t O2I = reactor.componentIndex("O2")-1;
    double O2_Moles = imw[O2I] * 0.5 * mass;
    double  H2_Moles = imw[H2I] * 0.5 * mass;
    // test getState
    reactor.getState(state.data());
    EXPECT_EQ(state[reactor.componentIndex("H2")], H2_Moles);
    EXPECT_EQ(state[reactor.componentIndex("O2")], O2_Moles);
    // Test updateState
    reactor.updateState(state.data());
    reactor.getState(updatedState.data());
    EXPECT_EQ(reactor.volume(), 0.5);
    EXPECT_EQ(reactor.pressure(), 101325);
    for (size_t i = 0; i < reactor.neq(); i++)
    {
       EXPECT_EQ(state[i], updatedState[i]);
    }
}

int main(int argc, char** argv)
{
    printf("Running main() from test_zeroD.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    Cantera::make_deprecation_warnings_fatal();
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
