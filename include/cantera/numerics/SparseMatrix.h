/*
Programmer: Anthony Walker

This is the SparseMatrix class which serves as a wrapper for Sundials and eventually Eigen Sparse Matrices.
*/
#ifndef CANTERA_SPARSE_MATRIX_H
#define CANTERA_SPARSE_MATRIX_H
//Cantera Imports
#include "cantera/base/ctexceptions.h"
//Sundials imports
#include "sunmatrix/sunmatrix_sparse.h"
//Eigen Imports
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

typedef Eigen::SparseMatrix<double> SparseMatrix_Eigen;

namespace Cantera
{

class SparseMatrix
{
protected:
    double threshold=10e-8; //default 
    size_t dimensions[2];
public:
    SparseMatrix(/* args */){}
    SparseMatrix(const SparseMatrix &sparseMat){*this=sparseMat;} //Copy constructor
    virtual ~SparseMatrix(){} //destructor
    //!Use this function to get the threshold value for setting elements
    virtual double getThreshold();
    //!Use this function to set the threshold value to compare elements against
    virtual void setThreshold(double threshold);
    //!Use this function to set an element by the threshold
    virtual void setElementByThreshold(size_t row,size_t col, double element);
    virtual void setElement(size_t row, size_t col, double element)=0; //set element
    virtual double getElement(size_t row, size_t col)=0; //get element
    virtual void setDimensions(size_t nrows,size_t ncols, void* otherData=NULL)=0;
    virtual size_t* getDimensions();
};

class SundialsSparseMatrix : public SparseMatrix
{
protected:
    SUNMatrix matrix;
    size_t dimensions[3]; //Array pointer for storage of number of  rows and columns and max non zero elements
    size_t datactr=0; //counter for data array
    size_t colctr=0; //counter for column ptrs
public:
    SundialsSparseMatrix(/* args */){}
    ~SundialsSparseMatrix(){}
    SundialsSparseMatrix(const SundialsSparseMatrix &sparseMat){*this=sparseMat;} //Copy constructor
    //Getter declarations
    virtual double getElement(size_t row, size_t col); //get element
    virtual SUNMatrix* getMatrix();
    //Setter declarations
    virtual void setDimensions(size_t nrows,size_t ncols, void* otherData=NULL);
    virtual void setElement(size_t row, size_t col, double element);//set element
    virtual void setMatrix(SUNMatrix *sparseMatrix);
};

} //Namespace end bracket

#endif