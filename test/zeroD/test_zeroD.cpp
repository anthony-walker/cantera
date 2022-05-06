#include "gtest/gtest.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/zerodim.h"
#include "cantera/base/Interface.h"
#include "cantera/numerics/AdaptivePreconditioner.h"

using namespace Cantera;

// This test ensures that prior reactor initialization of a reactor does
// not affect later integration within a network. This example was
// adapted from test_reactor.py::test_equilibrium_HP.
TEST(ZeroDim, test_individual_reactor_initialization)
{
    // initial conditions
    double T0 = 1100.0;
    double P0 = 10 * OneAtm;
    double tol = 1e-7;
    std::string X0 = "H2:1.0, O2:0.5, AR:8.0";
    // reactor solution, phase, and kinetics objects
    std::shared_ptr<Solution> sol1 = newSolution("h2o2.yaml");
    sol1->thermo()->setState_TPX(T0, P0, X0);
    // set up reactor object
    Reactor reactor1;
    reactor1.insert(sol1);
    // initialize reactor prior to integration to ensure no impact
    reactor1.initialize();
    // setup reactor network and integrate
    ReactorNet network;
    network.addReactor(reactor1);
    network.initialize();
    network.advance(1.0);
    // secondary gas for comparison
    std::shared_ptr<Solution> sol2 = newSolution("h2o2.yaml");
    sol2->thermo()->setState_TPX(T0, P0, X0);
    sol2->thermo()->equilibrate("UV");
    // secondary reactor for comparison
    Reactor reactor2;
    reactor2.insert(sol2);
    reactor2.initialize(0.0);
    // get state of reactors
    vector_fp state1 (reactor1.neq());
    vector_fp state2 (reactor2.neq());
    reactor1.getState(state1.data());
    reactor2.getState(state2.data());
    // compare the reactors.
    EXPECT_EQ(reactor1.neq(), reactor2.neq());
    for (size_t i = 0; i < reactor1.neq(); i++) {
        EXPECT_NEAR(state1[i], state2[i], tol);
    }
}

TEST(MoleReactorTestSet, test_mole_reactor_get_state)
{
    // setting up solution object and thermo/kinetics pointers
    double tol = 1e-8;
    auto sol = newSolution("h2o2.yaml");
    sol->thermo()->setState_TPY(1000.0, OneAtm, "H2:0.5, O2:0.5");
    IdealGasConstPressureMoleReactor reactor;
    reactor.insert(sol);
    reactor.setInitialVolume(0.5);
    reactor.setEnergy(false);
    reactor.initialize();
    vector_fp state(reactor.neq());
    vector_fp updatedState(reactor.neq());
    // test get state
    const ThermoPhase& thermo = reactor.contents();
    const vector_fp& imw = thermo.inverseMolecularWeights();
    // prescribed state
    double mass = reactor.volume() * thermo.density();
    size_t H2I = reactor.componentIndex("H2")-1;
    size_t O2I = reactor.componentIndex("O2")-1;
    double O2_Moles = imw[O2I] * 0.5 * mass;
    double  H2_Moles = imw[H2I] * 0.5 * mass;
    // test getState
    reactor.getState(state.data());
    EXPECT_NEAR(state[reactor.componentIndex("H2")], H2_Moles, tol);
    EXPECT_NEAR(state[reactor.componentIndex("O2")], O2_Moles, tol);
    // test updateState
    reactor.updateState(state.data());
    reactor.getState(updatedState.data());
    EXPECT_NEAR(reactor.volume(), 0.5, tol);
    EXPECT_NEAR(reactor.pressure(), OneAtm, tol);
    for (size_t i = 0; i < reactor.neq(); i++) {
       EXPECT_NEAR(state[i], updatedState[i], tol);
    }
}

bool isclose(double a, double b, double tol)
{
    return std::abs(a-b)<tol;
}

TEST(AdaptivePreconditionerTestSet, test_is_close_to_jacobian)
{
    // constants
    double vol = 0.2;
    double sharedThreshold = 0;
    double gamma = 1;
    // constant pressure reactor
    IdealGasConstPressureMoleReactor r1;
    auto sol1 = newSolution("gri30.yaml");
    auto thermo1 = sol1->thermo();
    thermo1->setEquivalenceRatio(1, "CH4", "O2:1, N2:3.76");
    thermo1->setState_TP(1000, 101325);
    r1.insert(sol1);
    r1.setInitialVolume(vol);
    // constant volume reactor
    IdealGasMoleReactor r2;
    auto sol2 = newSolution("gri30.yaml");
    auto thermo2 = sol2->thermo();
    thermo2->setEquivalenceRatio(1, "CH4", "O2:1, N2:3.76");
    thermo2->setState_TP(1000, 101325);
    r2.insert(sol2);
    r2.setInitialVolume(vol);
    // create and add preconditioner
    AdaptivePreconditioner precon;
    precon.setThreshold(sharedThreshold);
    // create derivative settings
    AnyMap settings;
    settings["skip-third-bodies"] = true;
    settings["skip-falloff"] = true;
    // network
    ReactorNet network;
    network.addReactor(r1);
    network.addReactor(r2);
    network.setDerivativeSettings(settings);
    network.setLinSolverType(GMRES);
    network.setPreconditioner(precon);
    network.initialize();
    // take a step
    double ctime = network.step();
    // state produced within CVODES for this example
    vector_fp Ndot(network.neq(), 0.0);
    vector_fp N(network.neq(), 0.0);
    vector_fp rhs(network.neq(), 0.0);
    vector_fp output(network.neq(), 0.0);
    network.getState(N.data());
    network.getState(rhs.data());
    // two array for jacobian
    Array2D jacobian(network.neq(), network.neq(), 0.0);
    network.evalJacobian(ctime, N.data(), Ndot.data(), nullptr, &jacobian);
    network.preconditioner_setup_nothrow(ctime, N.data(), Ndot.data(), gamma);
    Eigen::SparseMatrix<double> preconJac = precon.getJacobian();
    size_t same = 0;
    size_t notsame = 0;
    for (size_t i = 0; i < network.neq(); i++) {
        for (size_t j = 0; j < network.neq(); j++) {
            if (isclose(jacobian.value(i, j), preconJac.coeffRef(i, j), 1e-8)) {
                same += 1;
            } else {
                notsame += 1;
            }
        }
    }
    // expect matches to be 10 times greater than non-matches
    EXPECT_GT(same, notsame * 10);
}

TEST(AdaptivePreconditionerTestSet, test_preconditioned_hydrogen_auto_ignition)
{
    // create an ideal gas mixture that corresponds to GRI-Mech 3.0
    auto sol = newSolution("gri30.yaml", "gri30", "None");
    auto gas = sol->thermo();
    // set the state
    gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
    // create a reactor
    IdealGasConstPressureMoleReactor r;
    // 'insert' the gas into the reactor and environment.
    r.insert(sol);
    // create preconditioner
    AdaptivePreconditioner precon;
    // create derivative settings
    AnyMap settings;
    settings["skip-third-bodies"] = true;
    settings["skip-falloff"] = true;
    // create reactor network and set to use preconditioner
    ReactorNet sim;
    sim.addReactor(r);
    sim.setDerivativeSettings(settings);
    sim.setLinSolverType(GMRES);
    sim.setPreconditioner(precon);
    // main loop
    double dt = 1.e-5; // interval at which output is written
    int nsteps = 100; // number of intervals
    for (int i = 1; i <= nsteps; i++) {
        double tm = i*dt;
        sim.advance(tm);
    }
}

TEST(AdaptivePreconditionerTestSet, test_utilities_get_set)
{
    // seed random number generator
    std::srand(std::time(NULL));
    // generate random test dimensions between 10 and 50
    size_t base = 10;
    size_t limit = 50;
    size_t randomNum = std::rand() % (limit-base+1) + base;
    std::vector<size_t> dims{randomNum, randomNum};
    // create preconditioner object
    AdaptivePreconditioner precon;
    // set threshold
    double thresh = (double) base+2;
    precon.setThreshold(thresh);
    // set dimensions randomly
    precon.setDimensions(&dims);
    // get dimensions newly
    std::vector<size_t> *preconDims = precon.dimensions();
    // test that the dimensions are set properly via setDimensions
    EXPECT_EQ (dims.at(0), preconDims->at(0));
    EXPECT_EQ (dims.at(1), preconDims->at(1));
    // testing some get/set utilities
    double randSetGet = std::rand();
    precon.setThreshold(randSetGet);
    EXPECT_EQ(precon.threshold(), randSetGet);
    precon.setAbsoluteTolerance(randSetGet);
    EXPECT_EQ(precon.absoluteTolerance(), randSetGet);
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
