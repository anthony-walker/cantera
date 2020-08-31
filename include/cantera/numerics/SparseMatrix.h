/*
Programmer: Anthony Walker

This is the SparseMatrix class which serves as a wrapper for Sundials and eventually Eigen Sparse Matrices.
*/
#ifndef CANTERA_SPARSE_MATRIX_H
#define CANTERA_SPARSE_MATRIX_H
//Cantera imports
#include "cantera/zeroD/ReactorNet.h"
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
        int nrows; //number of rows
        int ncols; //number of columns
        void* cvode_mem; //Cvode memory pointer for access to time step data
    public:
        SparseMatrix(); //default constructor
        ~SparseMatrix();      // destructor 
        SparseMatrix(int nrows, int ncols); //overloaded
        SparseMatrix(int nrows, int ncols, int maxNonZero); //overloaded
        SparseMatrix(MATTYPE *sparseMatrix); //overloaded
        SparseMatrix(SparseMatrix *OtherSparseMat); //copy constructor
        //Getter declaration
        double getElement(int row, int col); //get element
        MATTYPE getMatrix();
        double getThreshold();
        void* getCvodeMemoryPtr();
        //Setter declaration
        void setElement(int row, int col, double element);//set element
        void setMatrix(MATTYPE *sparseMatrix);
        void setThreshold(double threshold);
        void setCvodeMemoryPtr(void* cv_mem_ptr);
        //Threshold set
        void setElement(int row,int col, double element, bool checkThreshold);
    };

/*
    TEMPLATE FUNCTIONS 
*/

//Default Constructor
template<class MATTYPE> SparseMatrix<MATTYPE>::SparseMatrix()
{
    //Do nothing
}

//Default Destructor
template<class MATTYPE> SparseMatrix<MATTYPE>::~SparseMatrix()
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

//Cvode memory pointer setter
template<class MATTYPE> void SparseMatrix<MATTYPE>::setCvodeMemoryPtr(void* cv_mem_ptr)
{
    this->cvode_mem = cv_mem_ptr;
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

//Cvode memory pointer getter
template<class MATTYPE> void* SparseMatrix<MATTYPE>::getCvodeMemoryPtr()
{
    return this->cvode_mem;
}

} //Namespace end bracket

#endif