#include "gtest/gtest.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/transport.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/zerodim.h"
#include <stdlib.h>
#include <time.h>
#include <queue>

using namespace Cantera;

#if CT_SUNDIALS_VERSION >= 40
    TEST(AdaptivePreconditioning,test_run_sim)
    {
        // Setting up solution object and thermo/kinetics pointers
        std::shared_ptr<Solution> sol = newSolution("methane_twostep.yaml");
        sol->thermo()->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
        // Set up reactor object
        IdealGasConstPressureReactor reactor;
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
        network.setIntegratorType(&precon,GMRES);
        // network->setVerbose(); //Setting verbose to be true
        network.addReactor(reactor); //Adding reactor to network
        // Setting up simulation
        network.setInitialTime(0.0);
        network.setMaxTimeStep(0.1);
        network.setMaxSteps(10000);
        network.setTolerances(1e-6,1e-6);
        network.setSensitivityTolerances(1e-6,1e-6);
        network.step();
    }

    TEST(AdaptivePreconditioning,test_two_step_mechanism)
    {
        // Constants
        double volume = 1.0;
        double startTime = 0.0;
        size_t reactorStart = 0;
        // Setting up solution object and thermo/kinetics pointers
        std::shared_ptr<Solution> sol = newSolution("methane_twostep.yaml");
        sol->thermo()->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0, CO:0");
        // Set up reactor object
        IdealGasConstPressureReactor reactor;
        reactor.insert(sol);
        reactor.setInitialVolume(volume);
        reactor.initialize();
        // State produced within CVODES for this example
        std::vector<double> y(reactor.neq(), 0.0);
        std::vector<double> ydot(reactor.neq(), 0.0);
        reactor.getState(y.data());
        // Internal preconditioner
        AdaptivePreconditioner internalPrecon;
        std::vector<size_t> preconDims{reactor.neq(), reactor.neq()};
        internalPrecon.initialize(&preconDims);
        internalPrecon.setReactorStart(reactorStart);
        internalPrecon.reactorLevelSetup(&reactor, reactorStart, startTime, y.data(), ydot.data(), nullptr);
        // Creating external preconditioner for comparison
        AdaptivePreconditioner externalPrecon;
        externalPrecon.initialize(internalPrecon.getDimensions());
        externalPrecon.setReactorStart(reactorStart);
        // Getting data for manual calculations
        std::vector<double> concs(sol->thermo()->nSpecies(), 0.0);
        std::vector<double> kf(sol->kinetics()->nReactions(), 0.0);
        std::vector<double> kr(sol->kinetics()->nReactions(), 0.0);
        sol->thermo()->getConcentrations(concs.data());
        sol->kinetics()->getFwdRateConstants(kf.data());
        sol->kinetics()->getRevRateConstants(kr.data());
        // Assign concentrations to species
        double O2 = concs[0];
        double CH4 = concs[2];
        double CO = concs[4];
        // Setting elements
        // Mass
        externalPrecon.setElement(0, 0, 1);
        // O2
        externalPrecon.setElement(2, 2, (- 2.25 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.25 * kf[1] * CO * std::pow(O2, -0.5)) / volume); // dO2/dO2
        externalPrecon.setElement(2, 3, 0); // dO2/dH20
        externalPrecon.setElement(2, 4, - 1.5 * kf[0] * std::pow(O2, 1.5) / volume); // dO2/dCH4
        externalPrecon.setElement(2, 5, 0.5 * kr[1]); // dO2/dCO2
        externalPrecon.setElement(2, 6, - 0.5 * kf[1] * std::pow(O2, 0.5)); // dO2/dCO
        // H2O
        externalPrecon.setElement(3, 2, 3 * kf[0] * CH4 * std::pow(O2, 0.5) / volume); // dH2O/dO2
        externalPrecon.setElement(3, 3, 0); // dH2O/dH20
        externalPrecon.setElement(3, 4, 2 * kf[0] * std::pow(O2, 1.5) / volume); // dH2O/dCH4
        externalPrecon.setElement(3, 5, 0); // dH2O/dCO2
        externalPrecon.setElement(3, 6, 0); // dH2O/dCO
        // CH4
        externalPrecon.setElement(4, 2, - 1.5 * kf[0] * CH4 * std::pow(O2, 0.5) / volume); // dCH4/dO2
        externalPrecon.setElement(4, 3, 0); // dCH4/dH20
        externalPrecon.setElement(4, 4, - kf[0] * std::pow(O2, 1.5) / volume); // dCH4/dCH4
        externalPrecon.setElement(4, 5, 0); // dCH4/dCO2
        externalPrecon.setElement(4, 6, 0); // dCH4/dCO
        // CO2
        externalPrecon.setElement(5, 2, (0.5 * kf[1] * CO * std::pow(O2, -0.5)) / volume); // dCO2/dO2
        externalPrecon.setElement(5, 3, 0); // dCO2/dH20
        externalPrecon.setElement(5, 4, 0); // dCO2/dCH4
        externalPrecon.setElement(5, 5, - kr[1]); // dCO2/dCO2
        externalPrecon.setElement(5, 6, kf[1] * std::pow(O2, 0.5)); // dCO2/CO
        //CO
        externalPrecon.setElement(6, 2, 1.5 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.5 * kf[1] * CO * std::pow(O2, -0.5) / volume); // dCO/dO2
        externalPrecon.setElement(6, 3, 0); // dCO/dH20
        externalPrecon.setElement(6, 4, kf[0] * std::pow(O2, 1.5) / volume); // dCO/dCH4
        externalPrecon.setElement(6, 5, kr[1]); // dCO/dCO2
        externalPrecon.setElement(6, 6, - kf[1] * std::pow(O2,0.5)); // dCO/CO
        // Temperature Derivatives
        externalPrecon.TemperatureDerivatives(&reactor, startTime, y.data(), ydot.data(), nullptr);
        // Check that the two are equal
        ASSERT_EQ(externalPrecon==internalPrecon, true);
    }

    TEST(AdaptivePreconditioning,test_one_step_mechanism_network)
    {
        // Constants
        double volume = 1.0;
        double startTime = 0.0;
        double sharedThreshold = 1e-16;
        // Setting up solution object and thermo/kinetics pointers one
        std::shared_ptr<Solution> sol1 = newSolution("methane_onestep.yaml");
        sol1->thermo()->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
        // Set up reactor object
        IdealGasConstPressureReactor reactor1;
        reactor1.insert(sol1);
        reactor1.setInitialVolume(volume);
        // Setting up solution object and thermo/kinetics pointers two
        std::shared_ptr<Solution> sol2 = newSolution("methane_onestep.yaml");
        sol2->thermo()->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
        // Set up reactor object
        IdealGasConstPressureReactor reactor2;
        reactor2.insert(sol2);
        reactor2.setInitialVolume(volume);
        // Network
        ReactorNet network;
        network.addReactor(reactor1);
        network.addReactor(reactor2);
        //Create and add preconditioner
        AdaptivePreconditioner internalPrecon;
        internalPrecon.setThreshold(sharedThreshold);
        network.setIntegratorType(&internalPrecon,GMRES);
        network.initialize();
        // State produced within CVODES for this example
        std::vector<double> ydot(network.neq(), 0.0);
        std::vector<double> y(network.neq(), 0.0);
        network.getState(y.data());
        network.preconditionerSetup(startTime, y.data(), ydot.data(), nullptr);
        // Creating external preconditioner for comparison
        AdaptivePreconditioner externalPrecon;
        externalPrecon.initialize(internalPrecon.getDimensions());
        externalPrecon.setThreshold(sharedThreshold);
        externalPrecon.setReactorStart(0);
        // Getting data for manual calculations
        std::vector<double> concs(sol1->thermo()->nSpecies(), 0.0);
        std::vector<double> kf(sol1->kinetics()->nReactions(), 0.0);
        sol1->thermo()->getConcentrations(concs.data());
        sol1->kinetics()->getFwdRateConstants(kf.data());
        // Assign concentrations to species
        double O2 = concs[0];
        double CH4 = concs[2];
        // Setting elements
        // Mass
        externalPrecon.setElement(0, 0, 1);
        // O2
        externalPrecon.setElement(2, 2, -4 * kf[0] * O2 * CH4 / volume); // dO2/dO2
        externalPrecon.setElement(2, 3, 0); // dO2/dH20
        externalPrecon.setElement(2, 4, -2 * kf[0] * O2 * O2 / volume); // dO2/dCH4
        externalPrecon.setElement(2, 5, 0); // dO2/dCO2
        // H2O
        externalPrecon.setElement(3, 2, 4 * kf[0]*O2*CH4/volume); // dH2O/dO2
        externalPrecon.setElement(3, 3, 0); // dH2O/dH20
        externalPrecon.setElement(3, 4, 2 * kf[0] * O2 * O2 / volume); // dH2O/dCH4
        externalPrecon.setElement(3, 5, 0); // dH2O/dCO2
        // CH4
        externalPrecon.setElement(4, 2, -2 * kf[0] * O2 * CH4 / volume); // dCH4/dO2
        externalPrecon.setElement(4, 3, 0); // dCH4/dH20
        externalPrecon.setElement(4, 4, -kf[0] * O2 * O2 / volume); // dCH4/dCH4
        externalPrecon.setElement(4, 5, 0); // dCH4/dCO2
        // CO2
        externalPrecon.setElement(5, 2, 2 * kf[0] * O2 * CH4 / volume); // dCO2/dO2
        externalPrecon.setElement(5, 3, 0); // dCO2/dH20
        externalPrecon.setElement(5, 4, kf[0] * O2 * O2 / volume); // dCO2/dCO2
        externalPrecon.setElement(5, 5, 0); // dCO2/dCO2
        externalPrecon.TemperatureDerivatives(&reactor1, startTime, y.data(), ydot.data(), nullptr);
        size_t neq = reactor1.neq();
        for (size_t i = 0; i < neq; i++)
        {
            for (size_t j = 0; j < neq; j++)
            {
                externalPrecon.setElement(i+neq,j+neq,externalPrecon.getElement(i,j));
            }
        }
        // Check that the two are equal
        ASSERT_EQ(externalPrecon==internalPrecon, true);
    }
#endif

// This tests initialization and dimension setting
TEST(AdaptivePreconditioning,test_initialize_set_dimensions)
{
    //Seed random number generator
    std::srand(std::time(NULL));
    // Generate random test dimensions between 5 and 15
    size_t base = 5;
    size_t limit = 15;
    std::vector<size_t> dims{std::rand() % (limit-base+1) + base, std::rand() % (limit-base+1) + base};
    // Create preconditioner object
    AdaptivePreconditioner precon;
    // Initialize matrix
    precon.initialize(&dims);
    // Test get dimensions
    std::vector<size_t> *preconDims = precon.getDimensions();
    // Test that the dimensions are set properly via initialize
    ASSERT_EQ (dims.at(0),preconDims->at(0));
    ASSERT_EQ (dims.at(1),preconDims->at(1));
    // Checking inside of sparse matrix
    Eigen::SparseMatrix<double>* SpMat = precon.getMatrix();
    ASSERT_EQ (dims.at(0),(size_t) SpMat->innerSize());
    ASSERT_EQ (dims.at(1),(size_t) SpMat->outerSize());
    // Test Set Dimensions
    std::vector<size_t> newDims{std::rand() % (limit-base+1) + base, std::rand() % (limit-base+1) + base};
    // Initialize with different dimensions
    precon.setDimensions(&newDims);
    // Test that the dimensions are set properly via setDimensions
    ASSERT_EQ (newDims.at(0),preconDims->at(0));
    ASSERT_EQ (newDims.at(1),preconDims->at(1));
}

TEST(AdaptivePreconditioning,test_get_set_copy_assignment_compare)
{
    // Seed random number generator
    std::srand(std::time(NULL));
    // Generate random test dimensions between 10 and 50
    size_t base = 10;
    size_t limit = 50;
    std::vector<size_t> dims{std::rand() % (limit-base+1) + base, std::rand() % (limit-base+1) + base};
    // Create preconditioner object
    AdaptivePreconditioner precon;
    // Initialize matrix
    precon.initialize(&dims);
    // Set threshold
    double thresh = (double) base+2;
    precon.setThreshold(thresh);
    // Random values to put in matrix
    std::queue<size_t> values;
    for(size_t i=0; i<base; i++)
    {
        values.push(std::rand() % 100);
        values.push(std::rand() % dims[1]);
        values.push(std::rand() % dims[0]);
    }
    // Check set and get elements
    for(size_t i=0; i<base; i++)
    {
        double currElement = (double) values.front();
        values.pop();
        size_t col = values.front();
        values.pop();
        size_t row = values.front();
        values.pop();
        precon.setElement(row,col,currElement);
        double returnedElement = precon.getElement(row,col);
        if(std::abs(currElement) >= thresh || row==col)
        {
            ASSERT_EQ (returnedElement,currElement);
        }
        else
        {
            ASSERT_EQ (returnedElement,0.0);
        }
    }
    // Create preconditioner object for copy and compare
    AdaptivePreconditioner preconCopy = precon; //call copy constructor
    // Call comparison
    ASSERT_EQ(preconCopy==precon,true);
    // Reset preconditioner then compare again
    precon.reset();
    ASSERT_EQ(preconCopy==precon,false);
    // Call assignment then compare again
    precon=preconCopy;
    ASSERT_EQ(preconCopy==precon,true);
}

TEST(DISABLED_AdaptivePreconditioning, hydrogen_auto_ignition)
{
    // See C++ examples for this - plan to complete this soon.
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
