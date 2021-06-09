#include "gtest/gtest.h"
//Creating kinetics and thermo objects
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/transport.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
//Creating reactor and reactor network
#include "cantera/zerodim.h"
//Other imports
#include <stdlib.h>
#include <time.h>
#include <queue>

using namespace Cantera;

TEST(AdaptivePreconditioning,single_step_mechanism_test)
{
    //Declaring ReactorNet related components
    IdealGasConstPressureReactor* combustor;
    Reservoir* inlet;
    Reservoir* exhaust;
    MassFlowController* inletMassFlowController;
    PressureController* outletMassFlowController;
    ReactorNet* network;
    //Setting up solution object and thermo/kinetics pointers
    std::shared_ptr<Solution> sol = newSolution("methane_onestep_multiphase.yaml");
    std::shared_ptr<ThermoPhase> gas= sol->thermo();
    std::shared_ptr<Kinetics> kin= sol->kinetics();
    gas->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, N2:3.76");
    //Set up reactor object
    combustor = new IdealGasConstPressureReactor();
    combustor->setKineticsMgr(*kin);
    combustor->setThermoMgr(*gas);
    combustor->setInitialVolume(1.0);
    //Creating inlet reservoir object and adding gas
    inlet = new Reservoir();
    inlet->insert(*gas);
    //equilibrate gas with enthalpy and pressure
    gas->equilibrate("HP");
    //Creating exhaust reservoir object and adding gas
    exhaust = new Reservoir();
    exhaust->insert(*gas);
    //Creating mass flow controllers
    inletMassFlowController = new MassFlowController();
    outletMassFlowController = new PressureController();
    //Connecting reactors
    inletMassFlowController->install(*inlet,*combustor);
    outletMassFlowController->install(*combustor,*exhaust);
    outletMassFlowController->setMaster(inletMassFlowController);
    outletMassFlowController->setPressureCoeff(0.01);
    //Set constant massflow rate
    double mf = 1.0;
    inletMassFlowController->setMassFlowRate(mf);
    //Creating preconditioner
    AdaptivePreconditioner* precon = new AdaptivePreconditioner();
    //Creating reactor network
    network = new ReactorNet();
    // network->setVerbose(); //Setting verbose to be true
    network->addReactor(*combustor); //Adding combustor to network
    //Setting up simulation
    network->setIntegratorType(precon,GMRES);
    network->setInitialTime(0.0);
    network->setMaxTimeStep(0.1);
    network->setMaxSteps(10000);
    network->setTolerances(1e-6,1e-6);
    network->setSensitivityTolerances(1e-6,1e-6);
    network->step();
    delete precon;
}

//This tests initialization and dimension setting of sparse matrix inside preconditioner.
TEST(AdaptivePreconditioning,initialize_dimensions_sparse_matrix)
{
    //Seed random number generator
    std::srand(std::time(NULL));
    //Generate random test dimensions between 5 and 15
    size_t base = 5;
    size_t limit = 15;
    size_t matrixDimOne = std::rand() % (limit-base+1) + base;
    size_t matrixDimTwo = std::rand() % (limit-base+1) + base;
    //Create preconditioner object
    AdaptivePreconditioner* precon = new AdaptivePreconditioner();
    //Initialize matrix
    precon->initialize(matrixDimOne,matrixDimTwo);
    //Test get dimensions
    size_t* preconDims = precon->getDimensions();
    //Test that the dimensions are set properly via initialize
    ASSERT_EQ (matrixDimOne,preconDims[0]);
    ASSERT_EQ (matrixDimTwo,preconDims[1]);
    //Checking inside of sparse matrix
    Eigen::SparseMatrix<double>* SpMat = precon->getMatrix();
    ASSERT_EQ (matrixDimOne,SpMat->innerSize());
    ASSERT_EQ (matrixDimTwo,SpMat->outerSize());
    //Test Set Dimensions
    size_t newDimOne = std::rand() % (limit-base+1) + base;
    size_t newDimTwo = std::rand() % (limit-base+1) + base;
    //Initialize with different dimensions
    precon->setDimensions(newDimOne,newDimTwo);
    //Test that the dimensions are set properly via setDimensions
    ASSERT_EQ (newDimOne,preconDims[0]);
    ASSERT_EQ (newDimTwo,preconDims[1]);
}

//Test get and set of sparse matrix and adaptive preconditioner
TEST(AdaptivePreconditioning,get_set_sparse_matrix_elements)
{
    //Seed random number generator
    std::srand(std::time(NULL));
    //Generate random test dimensions between 10 and 50
    size_t base = 10;
    size_t limit = 50;
    size_t matrixDimOne = std::rand() % (limit-base+1) + base;
    size_t matrixDimTwo = std::rand() % (limit-base+1) + base;
    //Create preconditioner object
    AdaptivePreconditioner* precon = new AdaptivePreconditioner();
    //Initialize matrix
    precon->initialize(matrixDimOne,matrixDimTwo);
    std::cout<<matrixDimOne<<" and "<<matrixDimTwo<<std::endl;
    // Set threshold
    size_t thresh = base+2;
    precon->setThreshold((double) thresh);
    //Random values to put in matrix
    std::queue<size_t> values;
    for(size_t i=0; i<base; i++)
    {
        values.push(std::rand());
        values.push(std::rand() % matrixDimTwo);
        values.push(std::rand() % matrixDimOne);
    }
    // Check set and get elements
    for(size_t i=0; i<base; i++)
    {
        size_t currElement = values.front();
        values.pop();
        size_t col = values.front();
        values.pop();
        size_t row = values.front();
        values.pop();
        precon->setElementByThreshold(row,col,currElement);
        size_t returnedElement = (size_t) precon->getElement(row,col);
        EXPECT_EQ (returnedElement,currElement);
    }
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
