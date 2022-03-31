#include "gtest/gtest.h"
// #include "cantera/base/Solution.h"
// #include "cantera/base/utilities.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/zerodim.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
#include <stdlib.h>
#include <time.h>

using namespace Cantera;

bool isclose(double a, double b, double tol)
{
    return std::abs(a-b)<tol;
}


TEST(AdaptivePreconditionerTestSet, test_is_close_to_jacobian)
{
    // Constants
    double vol = 0.2;
    double sharedThreshold = 0;
    double gamma = 1;
    // Constant pressure reactor
    IdealGasConstPressureMoleReactor r1;
    auto sol1 = newSolution("gri30.yaml");
    auto thermo1 = sol1->thermo();
    thermo1->setEquivalenceRatio(1, "CH4", "O2:1, N2:3.76");
    thermo1->setState_TP(1000, 101325);
    r1.insert(sol1);
    r1.setInitialVolume(vol);
    // Constant volume reactor
    IdealGasMoleReactor r2;
    auto sol2 = newSolution("gri30.yaml");
    auto thermo2 = sol2->thermo();
    thermo2->setEquivalenceRatio(1, "CH4", "O2:1, N2:3.76");
    thermo2->setState_TP(1000, 101325);
    r2.insert(sol2);
    r2.setInitialVolume(vol);
    // Network
    ReactorNet network;
    network.addReactor(r1);
    network.addReactor(r2);
    // Create and add preconditioner
    AdaptivePreconditioner precon;
    precon.setThreshold(sharedThreshold);
    network.setProblemType(GMRES);
    network.setPreconditioner(precon);
    network.initialize();
    // take a step
    double ctime = network.step();
    // State produced within CVODES for this example
    vector_fp Ndot(network.neq(), 0.0);
    vector_fp N(network.neq(), 0.0);
    vector_fp rhs(network.neq(), 0.0);
    vector_fp output(network.neq(), 0.0);
    network.getState(N.data());
    network.getState(rhs.data());
    // Two array for jacobian
    Array2D jacobian(network.neq(), network.neq(), 0.0);
    network.evalJacobian(ctime, N.data(), Ndot.data(), nullptr, &jacobian);
    network.preconditioner_setup_nothrow(ctime, N.data(), Ndot.data(), gamma);
    Eigen::SparseMatrix<double> preconJac = precon.getJacobian();
    size_t same = 0;
    size_t notsame = 0;
    for (size_t i = 0; i < network.neq(); i++)
    {
        for (size_t j = 0; j < network.neq(); j++)
        {
            if (isclose(jacobian.value(i, j), preconJac.coeffRef(i, j), 1e-8))
            {
                same += 1;
            }
            else
            {
                notsame += 1;
            }
        }
    }
    // expect matches to be 10 times greater than non-matches
    EXPECT_GT(same, notsame * 10);
}

TEST(AdaptivePreconditionerTestSet, test_run_sim)
{
    // Setting up solution object and thermo/kinetics pointers
    auto sol = newSolution("methane_twostep.yaml");
    sol->thermo()->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
    // Set up reactor object
    IdealGasConstPressureMoleReactor reactor;
    reactor.insert(sol);
    reactor.setInitialVolume(1.0);
    // Creating inlet reservoir object and adding gas
    Reservoir inlet;
    inlet.insert(sol);
    //Creating exhaust reservoir object and adding gas
    Reservoir exhaust;
    exhaust.insert(sol);
    // Creating mass flow controllers
    MassFlowController inletMassFlowController;
    PressureController outletMassFlowController;
    //Connecting reactors
    inletMassFlowController.install(inlet,reactor);
    outletMassFlowController.install(reactor,exhaust);
    outletMassFlowController.setMaster(&inletMassFlowController);
    outletMassFlowController.setPressureCoeff(0.01);
    // Set constant massflow rate
    inletMassFlowController.setMassFlowRate(1.0);
    // Creating reactor network
    ReactorNet network;
    network.addReactor(reactor); // Adding reactor to network
    // Create and add preconditioner
    AdaptivePreconditioner precon;
    network.setProblemType(GMRES);
    network.setPreconditioner(precon);
    // Running simulation
    network.initialize();
    network.advance(0.1);
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
    // create reactor network and set to use preconditioner
    ReactorNet sim;
    sim.addReactor(r);
    sim.setProblemType(GMRES);
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
    // Seed random number generator
    std::srand(std::time(NULL));
    // Generate random test dimensions between 10 and 50
    size_t base = 10;
    size_t limit = 50;
    size_t randomNum = std::rand() % (limit-base+1) + base;
    std::vector<size_t> dims{randomNum, randomNum};
    // Create preconditioner object
    AdaptivePreconditioner precon;
    auto precon_mat = precon.getMatrix();
    precon_mat->resize(dims[0], dims[1]);
    // Set threshold
    double thresh = (double) base+2;
    precon.setThreshold(thresh);
    // Set dimensions randomly
    precon.setDimensions(&dims);
    // Get dimensions newly
    std::vector<size_t> *preconDims = precon.getDimensions();
    // Test that the dimensions are set properly via setDimensions
    EXPECT_EQ (dims.at(0), preconDims->at(0));
    EXPECT_EQ (dims.at(1), preconDims->at(1));
    // Testing some get/set utilities
    double randSetGet = std::rand();
    precon.setThreshold(randSetGet);
    EXPECT_EQ(precon.getThreshold(), randSetGet);
    precon.setAbsoluteTolerance(randSetGet);
    EXPECT_EQ(precon.getAbsoluteTolerance(), randSetGet);
}

TEST(PreconditionerBaseTestSet, test_preconditioner_base)
{
    // Required variables
    PreconditionerBase precon;
    IdealGasConstPressureMoleReactor pressureReactor;
    Reactor generalReactor;
    ReactorNet network;
    // Test error throwing of preconditioner base
    EXPECT_THROW(precon.solve(0, nullptr, nullptr), NotImplementedError);
    EXPECT_THROW(precon.setup(), NotImplementedError);
    EXPECT_EQ(PRECONDITIONER_NOT_SET, precon.getPreconditionerMethod());
}

TEST(PreconditionerBaseTestSet, test_funceval_precon_funcs)
{
    // Setting up solution to insert in reactor
    auto sol = newSolution("h2o2.yaml");
    // Set up reactor object
    Reactor reactor;
    reactor.insert(sol);
    ReactorNet network;
    network.addReactor(reactor);
    network.initialize();
    EXPECT_THROW(network.FuncEval::preconditionerSetup(0, nullptr, nullptr, 0.0), NotImplementedError);
    EXPECT_THROW(network.FuncEval::preconditionerSolve(0, nullptr, nullptr, nullptr, nullptr), NotImplementedError);
}

int main(int argc, char** argv)
{
    printf("Running main() from test_preconditioners.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    Cantera::make_deprecation_warnings_fatal();
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
