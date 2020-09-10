/*
Programmer: Anthony Walker

This is the SparseMatrix class which serves as a wrapper for Sundials and eventually Eigen Sparse Matrices.
*/
#include "cantera/numerics/SparseMatrix.h"

using namespace Cantera; //Need this because SparseMatrix.h is defined in the Cantera namespace
/**
 * 
 * SparseMatrix Implementation
 * 
 * **/
double SparseMatrix::getThreshold()
{
    return this->threshold;
}

void SparseMatrix::setThreshold(double threshold)
{
    this->threshold = threshold;
}

void SparseMatrix::setElementByThreshold(size_t row,size_t col, double element)
{
    if (this->threshold<element)
    {
        this->setElement(row,col,element);
    }
}


/**
 * 
 * SundialsSparseMatrix Implementation
 * 
 * **/

void SundialsSparseMatrix::setDimensions(size_t nrows,size_t ncols, void* otherData)
{
    this->dimensions[0] = nrows;
    this->dimensions[1] = ncols;
    this->dimensions[2] = (nrows*ncols)/2 ? !(otherData):(size_t)otherData;
    // this->matrix = SUNSparseMatrix(nrows,ncols,this->dimensions[2],CSC_MAT);
}

void SundialsSparseMatrix::setElement(size_t row,size_t col,double element)
{
    // SUNSparseMatrix_Data(this->matrix)[row]=element;
}

double SundialsSparseMatrix::getElement(size_t row,size_t col)
{
    return 1.0;//SUNSparseMatrix_Data(this->matrix)[row]; //FIXME
}

