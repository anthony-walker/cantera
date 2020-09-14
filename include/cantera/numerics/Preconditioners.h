/*
Programmer: Anthony Walker
This is the this file contains functions to adaptively precondition the sparse matrix class
*/
#ifndef PRECONDITIONERS_H
#define PRECONDITIONERS_H

const int ADAPTIVE_MECHANISM_PRECONDITIONER = 1;
const int PRECONDITIONER_NOT_SET = 0;

//Cantera Imports
#include "cantera/base/ctexceptions.h"
//Eigen Imports
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

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
        virtual void setElement(unsigned long row, unsigned long col, double element)=0; //set element
        //!Function used to get a specific element of the matrix structure
        //!@param row unsigned long specifying the row location
        //!@param col unsigned long specifying the column location
        virtual double getElement(unsigned long row, unsigned long col)=0; //get element
        //!Function used to set the dimensions of and construct the matrix structure - required for initialization and use of the class
        //!@param nrows unsigned long number of rows in the structure
        //!@param ncols unsigned long nubmer of columns in the structure
        //!@param otherData void* for passing other data necessary for subclasses to initialize the matrix structure
        virtual void setDimensions(unsigned long nrows,unsigned long ncols, void* otherData=NULL)=0;
        //!Function to return the dimensions of the matrix structure
        virtual unsigned long* getDimensions();
        //!Function to solve a linear system Ax=b where A is the preconditioner contained in this matrix
        //@param x a double pointer to the vector (array) to store inv(A)*b
        //@param b a double pointer to the vector (array) multiplied by inv(A)
        virtual void solveLinearSystem(double* x, double *b)=0;
  };

  class Preconditioner : public PreconditionerBase
  {
    protected:
        //@param matrix a Eigen::SparseMatrix<double> type of structure used for storing data
        Eigen::SparseMatrix<double> *matrix;
    public:
        Preconditioner(/* args */);
        ~Preconditioner();
        Preconditioner(const Preconditioner &preconditioner){*this=preconditioner;} //Copy constructor
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
        virtual void setDimensions(unsigned long nrows,unsigned long ncols, void* otherData=NULL);
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
        virtual void solveLinearSystem(double* x, double *b);
    };
}

#endif