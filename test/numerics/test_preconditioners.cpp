#include "gtest/gtest.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/transport.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/base/utilities.h"
#include "cantera/zerodim.h"
#include <stdlib.h>
#include <time.h>
#include <queue>

using namespace Cantera;

bool isclose(double a, double b, double tol)
{
    return std::abs(a-b)<tol;
}


// void twoStepManualPrecondition(AdaptivePreconditioner& externalPrecon, MoleReactor& reactor, double* N, double* Ndot, double t=0.0, size_t loc=0)
//     {
//         const ThermoPhase& thermo = reactor.contents();
//         auto kinetics = reactor.getKineticsMgr();
//         double volInv = 1/reactor.volume();
//         double vol = reactor.volume();
//         // Getting data for manual calculations
//         size_t numberOfReactions = kinetics->nReactions();
//         size_t numberOfSpecies = kinetics->nTotalSpecies();
//         double moles = std::accumulate(N + 1, N + numberOfSpecies, 0.0);
//         // volume portion of constant pressure system
//         vector_fp omega(numberOfSpecies, 0.0);
//         kinetics->getNetProductionRates(omega.data());
//         size_t sidx = reactor.species_start();
//         for (size_t i = 0; i < numberOfSpecies; i++)
//         {
//             for (size_t j = 0; j < numberOfSpecies; j++)
//             {
//                 externalPrecon(i + sidx, j + sidx, vol / moles * omega[i]);
//             }
//         }
//         // vectors for manual calcs
//         vector_fp kf(numberOfReactions, 0.0);
//         vector_fp kr(numberOfReactions, 0.0);
//         // getting actual data
//         vector_fp C(thermo.nSpecies(), 0.0);
//         thermo.getConcentrations(C.data());
//         kinetics->getFwdRateConstants(kf.data());
//         kinetics->getRevRateConstants(kr.data());
//         // Assign concentrations to species
//         double O2 = C[0];
//         double CH4 = C[2];
//         double CO = C[4];
//         // Setting elements with manual rop derivatives
//         // O2
//         externalPrecon(1 + loc, 1 + loc, (- 2.25 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.25 * kf[1] * CO * std::pow(O2, -0.5)) * volInv); // dO2/dO2
//         externalPrecon(1 + loc, 3 + loc, - 1.5 * kf[0] * std::pow(O2, 1.5) * volInv); // dO2/dCH4
//         externalPrecon(1 + loc, 4 + loc, 0.5 * kr[1] * volInv); // dO2/dCO2
//         externalPrecon(1 + loc, 5 + loc, - 0.5 * kf[1] * std::pow(O2, 0.5) * volInv); // dO2/dCO
//         // H2O
//         externalPrecon(2 + loc, 1 + loc, 3 * kf[0] * CH4 * std::pow(O2, 0.5) * volInv); // dH2O/dO2
//         externalPrecon(2 + loc, 3 + loc, 2 * kf[0] * std::pow(O2, 1.5) * volInv); // dH2O/dCH4
//         // CH4
//         externalPrecon(3 + loc, 1 + loc, - 1.5 * kf[0] * CH4 * std::pow(O2, 0.5) * volInv); // dCH4/dO2
//         externalPrecon(3 + loc, 3 + loc, - kf[0] * std::pow(O2, 1.5) * volInv); // dCH4/dCH4
//         // CO2
//         externalPrecon(4 + loc, 1 + loc, (0.5 * kf[1] * CO * std::pow(O2, -0.5)) * volInv); // dCO2/dO2
//         externalPrecon(4 + loc, 4 + loc, -kr[1] * volInv); // dCO2/dCO2
//         externalPrecon(4 + loc, 5 + loc, kf[1] * std::pow(O2, 0.5) * volInv); // dCO2/CO
//         //CO
//         externalPrecon(5 + loc, 1 + loc, (1.5 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.5 * kf[1] * CO * std::pow(O2, -0.5)) * volInv); // dCO/dO2
//         externalPrecon(5 + loc, 3 + loc, kf[0] * std::pow(O2, 1.5) * volInv); // dCO/dCH4
//         externalPrecon(5 + loc, 4 + loc, kr[1] * volInv); // dCO/dCO2
//         externalPrecon(5 + loc, 5 + loc, -kf[1] * std::pow(O2,0.5) * volInv); // dCO/CO
//         // temperature derivatives
//         // rop derivatives w.r.t temperature
//         // vector_fp dwdot(numberOfSpecies, 0.0);
//         // kinetics->getNetProductionRates_ddT(dwdot.data());
//         // for (size_t i = 0; i < numberOfSpecies; i++)
//         // {
//         //     externalPrecon(i + sidx, 0, dwdot[i] * vol);
//         // }
//         // temperature derivatives w.r.t species

//     };


// TEST(AdaptivePreconditionerTestSet, test_manual_two_step_mechanism)
// {
//     // Constants
//     double volume = 0.2;
//     double gamma = 0.1;
//     double threshold = 0;
//     // Setting up solution object and thermo/kinetics pointers
//     auto sol = newSolution("methane_twostep.yaml");
//     auto thermo = sol->thermo();
//     auto kinetics = sol->kinetics();
//     thermo->setEquivalenceRatio(1, "CH4", "O2:1");
//     thermo->setState_TP(1000, 101325);
//     // Set up reactor object
//     IdealGasConstPressureMoleReactor reactor;
//     reactor.insert(sol);
//     reactor.setInitialVolume(volume);
//     reactor.initialize();
//     // Setup network
//     ReactorNet network;
//     // Internal preconditioner
//     AdaptivePreconditioner internalPrecon;
//     internalPrecon.setGamma(gamma);
//     internalPrecon.setThreshold(threshold);
//     // Create and add preconditioner
//     network.addReactor(reactor); //Adding reactor to network
//     network.setProblemType(GMRES);
//     network.setPreconditioner(internalPrecon);
//     network.initialize();
//     // State produced within CVODES for this example
//     vector_fp N(reactor.neq(), 0.0);
//     vector_fp Ndot(reactor.neq(), 0.0);
//     vector_fp NCopy(reactor.neq(), 0.0);
//     // Creating external preconditioner for comparison
//     AdaptivePreconditioner externalPrecon;
//     externalPrecon.setGamma(gamma);
//     externalPrecon.setThreshold(threshold);
//     externalPrecon.initialize(network);
//     // Compare the preconditioners while stepping
//     double ct = 0.0;
//     for (size_t i = 0; i < 1; i++)
//     {
//         // Get reactor state
//         reactor.getState(N.data());
//         // reset preconditioner
//         internalPrecon.reset();
//         // precondition from reactor object
//         network.preconditionerSetup(ct, N.data(), Ndot.data(), gamma);
//         // manual precondition
//         // set state to strictly positive composition for manual comparison
//         externalPrecon.getStrictlyPositiveComposition(reactor.neq(), N.data(), NCopy.data());
//         reactor.updateState(NCopy.data());
//         // reset external preconditioner
//         externalPrecon.reset();
//         twoStepManualPrecondition(externalPrecon, reactor, NCopy.data(), Ndot.data(), ct);
//         internalPrecon.printJacobian();
//         std::cout<< "----------------" <<std::endl;
//         externalPrecon.printJacobian();
//         // // post setup processes
//         // externalPrecon.setup();
//         // // check that the two are equal
//         // EXPECT_EQ(externalPrecon == internalPrecon, true);
//         // // reset state
//         // reactor.updateState(N.data());
//         // // step the network
//         // ct = network.step();
//     }
// }

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
    IdealGasMoleReactor reactor;
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
    // Create and add preconditioner
    AdaptivePreconditioner precon;
    network.addReactor(reactor); //Adding reactor to network
    network.setProblemType(GMRES);
    network.setPreconditioner(precon);
    // Setting up simulation
    network.initialize();
    network.step();
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
