/*
Programmer: Anthony Walker
This is the this file contains functions to adaptively precondition the sparse matrix class
*/
#ifndef PRECONDITIONERS_H
#define PRECONDITIONERS_H

const int ADAPTIVE_MECHANISM_PRECONDITIONER = 1;
const int PRECONDITIONER_NOT_SET = 0;

//PreconditionerBase imports
//Cantera Imports
#include "cantera/base/ctexceptions.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/FlowReactor.h"

//Eigen Imports
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

//Adaptive preconditioning imports
//Cantera imports
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/Reaction.h"

//Other imports

//Function Library for preconditioning
namespace Cantera//Making ASP apart of Cantera namespace
{
  class PreconditionerBase
  { 
    private:
        
    protected:
        //@param threshold a double value to selectively fill the matrix structure based on this threshold
        double threshold=10e-16; //default 
        //@param dimensions an unsigned int pointer to store dimensions
        unsigned long dimensions[2];
    public:
        PreconditionerBase(/* args */){}
        PreconditionerBase(const PreconditionerBase &precBase){*this=precBase;} //Copy constructor
        virtual ~PreconditionerBase(){} //destructor
        virtual unsigned long getPreconditionerType(){return PRECONDITIONER_NOT_SET;};
        /*

          Reactor Setup & Solve Functions

        */

       //Pure virtual functions

       //!Function to solve a linear system Ax=b where A is the preconditioner contained in this matrix
        //@param x a double pointer to the vector (array) to store inv(A)*b
        //@param b a double pointer to the vector (array) multiplied by inv(A)
        virtual void solve(double* x, double *b,unsigned long size)=0;
        //! This function performs the setup of the preconditioner for the specified reactor type and should be overloaded for each different reactor time
        //!@param reactor A Reactor object pointer
        //!@param reactorStart an unsigned long providing the index location in which the state of the given reactor starts
        virtual void setup(Reactor *reactor, double t, double* y, double* ydot, double* params, unsigned long reactorStart)=0;
        //!This function is called during setup for any processes that need to be completed prior to setup functions
        //! e.g. dynamic memory allocation
        virtual void initialize(unsigned long nrows,unsigned long ncols)=0;
        //!This function is called during setup for any processes that need to be completed post to setup functions
        //! e.g. dynamic memory allocation
        virtual void reset()=0;
        //!Function used to set a specific element of the matrix structure
        //!@param row unsigned long specifying the row location
        //!@param col unsigned long specifying the column location
        //!@param element double value to be inserted into matrix structure
        virtual void setElement(unsigned long row, unsigned long col, double element)=0; //set element
        //!Function used to get a specific element of the matrix structure
        //!@param row unsigned long specifying the row location
        //!@param col unsigned long specifying the column location
        virtual double getElement(unsigned long row, unsigned long col)=0; //get element

        //Other preconditioner functions

        //!Use this function to get the threshold value for setting elements
        virtual double getThreshold();
        //!Use this function to set the threshold value to compare elements against
        //!@param threshold double value used in setting by threshold
        virtual void setThreshold(double threshold);
        //!Use this function to set an element by the threshold
        //!@param row unsigned long specifying the row location
        //!@param col unsigned long specifying the column location
        //!@param element double value to be inserted into matrix structure
        virtual void setElementByThreshold(unsigned long row,unsigned long col, double element);
        //!Function used to set the dimensions of and construct the matrix structure - required for initialization and use of the class
        //!@param nrows unsigned long number of rows in the structure
        //!@param ncols unsigned long nubmer of columns in the structure
        //!@param otherData void* for passing other data necessary for subclasses to initialize the matrix structure
        virtual void setDimensions(unsigned long nrows,unsigned long ncols);
        //!Function to return the dimensions of the matrix structure
        virtual unsigned long* getDimensions();
        
  };
}

/**
 * 
 * 
 * Adaptive Mechanism Preconditioning Namespace
 * 
 * 
 * */
namespace Cantera::AMP //Making ASP apart of Cantera namespace
{

 //!Typedef used for getting indices based on strings

  typedef std::map<std::string,unsigned long> StateMap;
  typedef void (*AdaptiveFunction)(PreconditionerBase *preconditioner,Reactor* reactor, double** inputs,StateMap indexMap, std::string key);
  typedef std::map<std::string,AdaptiveFunction> FunctionMap;

 class AdaptivePreconditioner : public PreconditionerBase
  {
    protected:
        Eigen::SparseMatrix<double> matrix;
        FunctionMap functionMap;
        unsigned long nonzeros;
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        AdaptivePreconditioner(/* args */){};
        AdaptivePreconditioner(int rtype);
        ~AdaptivePreconditioner(){};
        AdaptivePreconditioner(const AdaptivePreconditioner &preconditioner){*this=preconditioner;} //Copy constructor
        virtual unsigned long getPreconditionerType(){return ADAPTIVE_MECHANISM_PRECONDITIONER;};
        //! This function performs the setup of the preconditioner for the specified reactor type and should be overloaded for each different reactor time
        //!@param reactor A Reactor object pointer
        //!@param reactorStart an unsigned long providing the index location in which the state of the given reactor starts
        virtual void setup(Reactor *reactor, double t, double* y, double* ydot, double* params, unsigned long reactorStart);
        //!This function is called during setup for any processes that need to be completed prior to setup functions
        //! e.g. dynamic memory allocation
        virtual void initialize(unsigned long nrows,unsigned long ncols);
        //!This function is called during setup for any processes that need to be completed post to setup functions
        //! e.g. dynamic memory allocation
        virtual void reset();
        //!Function used to get a specific element of the matrix structure
        //!@param row unsigned long specifying the row location
        //!@param col unsigned long specifying the column location
        virtual double getElement(unsigned long row, unsigned long col); //get element
        //!Function used to return compressed version of the matrix structure
        virtual Eigen::SparseMatrix<double>* getMatrix();
        //!Function used to set the dimensions of and construct the matrix structure - required for initialization and use of the class
        //!@param nrows unsigned long number of rows in the structure
        //!@param ncols unsigned long nubmer of columns in the structure
        //!@param otherData void* for passing other data necessary for subclasses to initialize the matrix structure
        virtual void setDimensions(unsigned long nrows,unsigned long ncols);
        //!Function used to set a specific element of the matrix structure
        //!@param row unsigned long specifying the row location
        //!@param col unsigned long specifying the column location
        //!@param element double value to be inserted into matrix structure
        virtual void setElement(unsigned long row, unsigned long col, double element);//set element
        //!Function used to set compressed version of the matrix structure
        //!@param sparseMatrix a SUNMatrix pointer to a type of SUNMatrix
        //!@param compress a bool dictating whether or not the set matrix needs compressed or not
        virtual void setMatrix(Eigen::SparseMatrix<double>* sparseMatrix);  
        //!Function to solve a linear system Ax=b where A is the preconditioner contained in this matrix
        //@param x a double pointer to the vector (array) to store inv(A)*b
        //@param b a double pointer to the vector (array) multiplied by inv(A)
        virtual void solve(double* x, double *b,unsigned long size);
        //!Function to add extra functions to function map
        //@param component a string type used as the map key
        //@param newFunction AdaptiveFunction type for the function to be added
        virtual void addToFunctionMap(std::string component, AdaptiveFunction newFunction);
        //!Function to remove undesirable functions from function map
        //@param component a string type used as the map key
        virtual void removeFromFunctionMap(std::string component);

    };

  //!This function returns an index map of nonspecies
  //!@param reactor the current Reactor object
  StateMap getStateMap(Reactor *reactor,unsigned long start);

  //! This function determines the rate of progress derivatives given a composition of reactants or products
  void checkEigenError(std::string method, unsigned long info);

  //! Use this function to print and check reactor components
  inline void printReactorComponents(Reactor* reactor);

  //!This function is a subfunction of SpeciesDerivatives that gets the species w.r.t species derivatives for each reaction
  inline void speciesSpeciesDerivative(std::map<std::string, double> comp,std::map<std::string,unsigned long> indexMap, double* omega, double* concentrations, double k_direction, double volume, unsigned long numberOfSpecies);

  //!This function is a subfunction of SpeciesDerivatives that gets the temperature w.r.t species derivatives
  inline void temperatureSpeciesDerivatives(PreconditionerBase *preconditioner, double *concentrations, double *dwdnj, double volume, ThermoPhase *thermo,Kinetics *kinetics, StateMap stateMap);
  
  //!This function is a subfunction of TemperatureDerivatives that gets the species w.r.t temperature derivatives
  inline void speciesTemperatureDerivatives(PreconditionerBase *preconditioner, double *netProductionRatesCurrent, double *netProductionRatesNext, double deltaTemp,unsigned long numberOfSpecies,unsigned long speciesStart, unsigned long tempIndex);

  //!This function is a subfunction of TemperatureDerivatives that gets the temperature dot w.r.t temperature derivatives
  inline void temperatureTemperatureDerivatives(PreconditionerBase *preconditioner, double TDotCurrent, double TDotNext, double deltaTemp, unsigned long tempIndex);

  
  //!This function does not precondition the associated equation by assigning it's preconditioner value to a value of 1
  //!@param row the row index of the variable
  //!@param col the column index of the variable
  void NoPrecondition(PreconditionerBase *preconditioner,Reactor* reactor, double** inputs,StateMap indexMap, std::string key);

  //! This function determines derivatives of Species and Temperature with respect to Temperature for jacobian preconditioning with a finite difference.
    //! @param *preconditioner A pointer to a PreconditionerBase Object for preconditioning the system and storing preconditioner values
    //! @param *reactor A pointer to the current reactor being precondition
    //! @param *ydot A pointer to the current data of ydot passed from CVODES
    //! @param meanSpecificHeat The mean specific heat used based on reactor type
    //! @param index The index location of temperature in the state vector
  void TemperatureDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, double** inputs,StateMap indexMap, std::string key);

  //! This function determines derivatives of Species with respect to species for jacobian preconditioning;
  //! specifically it determines the derivatives of the rate laws of all species with respect to other species in terms of moles.
  //! @param *preconditioner A pointer to a PreconditionerBase Object for preconditioning the system and storing preconditioner values
  //! @param *reactor A pointer to the current reactor being precondition
  void SpeciesDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, double** inputs,StateMap indexMap);

  //! This function is used to convert the system from mass fraction to mole fraction for solving the linear system with a mole based jacobian.
  //! @param *reactor A pointer to the current reactor being converted
  //! @param *tempState A double pointer to the temporary state used to solve the linear system
  //! @param *rhs A double pointer provided to preconditioner solve of the initial rhs state
  //! @param m_start An unsigned long for the global index of each species
  void ForwardConversion(Reactor *currReactor, double *tempState, double *rhs, unsigned long m_start);

  //! This function is used to convert the system from moles back to mass fraction after being solved with AMP preconditioner
  //! @param *reactor A pointer to the current reactor being converted
  //! @param *output A double pointer to the output vector of mole values
  //! @param m_start An unsigned long for the global index of each species
  void BackwardConversion(Reactor *currReactor, double *output, unsigned long m_start);
} 

#endif