/*
Programmer: Anthony Walker

This is the SparseMatrix class which serves as a wrapper for Sundials and eventually Eigen Sparse Matrices.
*/
#ifndef CANTERA_SPARSE_MATRIX_H
#define CANTERA_SPARSE_MATRIX_H
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
        size_t dimensions[2]; //Array pointer for storage of number of  rows and columns
        void* cvode_mem; //Cvode memory pointer for access to time step data
    public:
        SparseMatrix(); //default constructor
        ~SparseMatrix();      // destructor 
        SparseMatrix(SparseMatrix *OtherSparseMat); //copy constructor
        //Getter declarations
        size_t* getDimensions();
        double getElement(size_t row, size_t col); //get element
        MATTYPE getMatrix();
        double getThreshold();
        void* getCvodeMemoryPtr();
        //Setter declarations
        void setDimensions(size_t nrows,size_t ncols);
        void setElement(size_t row, size_t col, double element);//set element
        void setMatrix(MATTYPE *sparseMatrix);
        void setThreshold(double threshold);
        void setCvodeMemoryPtr(void* cv_mem_ptr);
        void setElementByThreshold(size_t row,size_t col, double element);
        //! Use this function to construct the object for use by
        // void buildSparseRepresentation(); 
    };

/*
    TEMPLATE FUNCTIONS 
*/

//Default Constructor - do nothing
template<class MATTYPE> SparseMatrix<MATTYPE>::SparseMatrix(){}

//Default Destructor - do nothing
template<class MATTYPE> SparseMatrix<MATTYPE>::~SparseMatrix(){}

/*

Setter functions

*/

//Set element by threshold
template<class MATTYPE> void SparseMatrix<MATTYPE>::setElementByThreshold(size_t row, size_t col, double element)
{   
    if (element > this->threshold)
    {
        this->setElement(row,col,element);
    }
}

//Dimensions setter
template<class MATTYPE> void SparseMatrix<MATTYPE>::setDimensions(size_t nrows, size_t ncols)
{
    this->dimensions[0] = nrows;
    this->dimensions[1] = ncols;
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
/*

Getter functions

*/

//dimensions getter
template<class MATTYPE> size_t* SparseMatrix<MATTYPE>::getDimensions()
{
    return this->dimensions;
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