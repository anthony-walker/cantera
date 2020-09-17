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
    protected:
        //@param threshold a double value to selectively fill the matrix structure based on this threshold
        double threshold=10e-8; //default 
        //@param dimensions an unsigned int pointer to store dimensions
        unsigned long dimensions[2];
    public:
        PreconditionerBase(/* args */){}
        PreconditionerBase(const PreconditionerBase &precBase){*this=precBase;} //Copy constructor
        virtual ~PreconditionerBase(){} //destructor

        /*

          Reactor Setup & Solve Functions

        */
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
        //!Function used to set a specific element of the matrix structure
        //!@param row unsigned long specifying the row location
        //!@param col unsigned long specifying the column location
        //!@param element double value to be inserted into matrix structure
        virtual void setElement(unsigned long row, unsigned long col, double element); //set element
        //!Function used to get a specific element of the matrix structure
        //!@param row unsigned long specifying the row location
        //!@param col unsigned long specifying the column location
        virtual double getElement(unsigned long row, unsigned long col); //get element
        //!Function used to set the dimensions of and construct the matrix structure - required for initialization and use of the class
        //!@param nrows unsigned long number of rows in the structure
        //!@param ncols unsigned long nubmer of columns in the structure
        //!@param otherData void* for passing other data necessary for subclasses to initialize the matrix structure
        virtual void setDimensions(unsigned long nrows,unsigned long ncols);
        //!Function to return the dimensions of the matrix structure
        virtual unsigned long* getDimensions();
        //!Function to solve a linear system Ax=b where A is the preconditioner contained in this matrix
        //@param x a double pointer to the vector (array) to store inv(A)*b
        //@param b a double pointer to the vector (array) multiplied by inv(A)
        virtual void solve(double* x, double *b,unsigned long size);
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
 class AdaptivePreconditioner : public PreconditionerBase
  {
    protected:
        Eigen::SparseMatrix<double> matrix;
        unsigned long nonzeros;
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        AdaptivePreconditioner(/* args */){};
        ~AdaptivePreconditioner(){};
        AdaptivePreconditioner(const AdaptivePreconditioner &preconditioner){*this=preconditioner;} //Copy constructor
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
    };

  void checkEigenError(std::string method, unsigned long info);

  inline void printReactorComponents(Reactor* reactor);

  inline void speciesDerivative(std::map<std::string, double> comp,std::map<std::string,size_t> indexMap, double* omega, double* concentrations, double k_direction, double volume);

  void NoPrecondition(PreconditionerBase *preconditioner,unsigned long row, unsigned long col);

  void TemperatureDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, double* ydot, double dTdt, size_t index, size_t speciesStart);
  
  void SpeciesSpeciesDerivatives(PreconditionerBase *preconditioner,Reactor* reactor, size_t speciesStart);

}

#endif