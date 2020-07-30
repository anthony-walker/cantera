/*
Programmer: Anthony Walker

This is the SparseMatrix class which serves as a wrapper for Sundials and eventually Eigen Sparse Matrices.
*/
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

//Sundials imports
#include "sunmatrix/sunmatrix_sparse.h"
//Eigen Imports
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

typedef Eigen::SparseMatrix<double> EigenSparseMatrix;
typedef SUNMatrix SundialsSparseMatrix;

namespace Cantera
{
template<class MATTYPE> class SparseMatrix
{
    private:
        MATTYPE *matrix; //Matrix stored by sparse matrix
        double threshold=10e-8; //default threshold value
    public:
        ~SparseMatrix(); //destructor
        SparseMatrix(); //default constructor
        SparseMatrix(int nrows, int ncols); //overloaded
        SparseMatrix(int nrows, int ncols, int maxNonZero); //overloaded
        SparseMatrix(MATTYPE *sparseMatrix); //overloaded
        SparseMatrix(SparseMatrix *OtherSparseMat); //copy constructor
        double getElement(int row, int col); //get element
        MATTYPE getMatrix();
        double getThreshold();
        void setElement(int row, int col, double element);//set element
        void setMatrix(MATTYPE *sparseMatrix);
        void setThreshold(double threshold);
    };

/*
    TEMPLATE FUNCTIONS 
*/

//Destructor
template<class MATTYPE> SparseMatrix<MATTYPE>::~SparseMatrix()
{
    delete this->matrix;
}

//Default Constructor
template<class MATTYPE> SparseMatrix<MATTYPE>::SparseMatrix()
{
    //Do nothing
}

//Overloaded Constructor
template<class MATTYPE> SparseMatrix<MATTYPE>::SparseMatrix(MATTYPE *sparseMatrix)
{
    this->matrix = sparseMatrix;
}

//Threshold setter
template<class MATTYPE> void SparseMatrix<MATTYPE>::setThreshold(double threshold)
{
    this->threshold = threshold;
}

//Matrix setter
template<class MATTYPE> void SparseMatrix<MATTYPE>::setMatrix(MATTYPE *sparseMatrix)
{
    this->matrix = sparseMatrix;
}

//Threshold getter
template<class MATTYPE> double SparseMatrix<MATTYPE>::getThreshold()
{
    return this->threshold;
}

//Matrix getter
template<class MATTYPE> MATTYPE SparseMatrix<MATTYPE>::getMatrix()
{
    return this->matrix;
}

/*
    SUNDIALS SPECIALIZED FUNCTIONS 
*/
//Overloaded Constructor
template<> SparseMatrix<SundialsSparseMatrix>::SparseMatrix(int nrows, int ncols, int maxNonZero)
{   
    SundialsSparseMatrix temporaryMatrix = SUNSparseMatrix(nrows,ncols,maxNonZero,CSC_MAT);
    this->matrix = &temporaryMatrix;
}


//Specialized SundialsSparseMatrix SETTER
template<> void SparseMatrix<SundialsSparseMatrix>::setElement(int row, int col, double element)
{
    
}

//Specialized SundialsSparseMatrix GETTER
template<> double SparseMatrix<SundialsSparseMatrix>::getElement(int row, int col)
{
    //IMPLEMENT ME
    return 1.0;
}

/*
    EIGEN SPECIALIZED FUNCTIONS 
*/
//Overloaded Constructor
template<> SparseMatrix<EigenSparseMatrix>::SparseMatrix(int nrows, int ncols)
{   
    EigenSparseMatrix temporaryMatrix = EigenSparseMatrix(nrows,ncols);
    this->matrix = &temporaryMatrix;
}


}

#endif