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
        int typeIndex; //This is used to run a certain switch case.
    public:
        SparseMatrix(); //default constructor
        SparseMatrix(int nrows, int ncols); //overloaded
        SparseMatrix(MATTYPE *sparseMatrix); //overloaded
        SparseMatrix(SparseMatrix *OtherSparseMat); //copy constructor
        double getElement(int row, int col); //get element
        int getTypeIndex(); // returns type index
        void setElement(int row, int col, double element);//set element
        void setMatrix(MATTYPE *SparseMatrix);
    };

/*
    TEMPLATE FUNCTIONS 
*/

//Default Constructor
template<class MATTYPE> SparseMatrix<MATTYPE>::SparseMatrix()
{
    //Matrix is not allocated
    if (std::is_same<MATTYPE,SUNMatrix>::value)
    {   
        this->typeIndex = 0; // Type id is zero for sunmatrix
    }
    else if (std::is_same<MATTYPE,EigenSparseMatrix>::value)
    {
        this->typeIndex = 1; // Type id is zero for Eigen sparse matrix
    }
    else //Throw error if type isn't defined properly
    {
        throw std::invalid_argument( "The type supplied to the template class is not supported." );
    }
}

//Overloaded Constructor
template<class MATTYPE> SparseMatrix<MATTYPE>::SparseMatrix(int nrows, int ncols)
{
    if (std::is_same<MATTYPE,SUNMatrix>::value)
    {   
        this->typeIndex = 0; // Type id is zero for sunmatrix
        MATTYPE temporaryMatrix = SUNSparseMatrix(nrows,ncols,nrows*ncols,CSC_MAT);
        this->matrix = &temporaryMatrix;
    }
    else if (std::is_same<MATTYPE,EigenSparseMatrix>::value)
    {
        this->typeIndex = 1; // Type id is zero for Eigen sparse matrix
        // MATTYPE temporaryMatrix = EigenSparseMatrix(nrows,ncols);
        // this->matrix = &temporaryMatrix;
    }
    else //Throw error if type isn't defined properly
    {
        throw std::invalid_argument( "The type supplied to the template class is not supported." );
    }   
}

//Overloaded Constructor
template<class MATTYPE> SparseMatrix<MATTYPE>::SparseMatrix(MATTYPE *sparseMatrix)
{
    //Matrix is not allocated
    if (std::is_same<MATTYPE,SUNMatrix>::value)
    {   
        this->typeIndex = 0; // Type id is zero for sunmatrix
    }
    else if (std::is_same<MATTYPE,EigenSparseMatrix>::value)
    {
        this->typeIndex = 1; // Type id is zero for Eigen sparse matrix
    }
    else //Throw error if type isn't defined properly
    {
        throw std::invalid_argument( "The type supplied to the template class is not supported." );
    }
    this->matrix = sparseMatrix;
}


//This function returns the typeIndex
template<class MATTYPE> int SparseMatrix<MATTYPE>::getTypeIndex()
{
    return this->typeIndex;
}

/*
    SUNDIALS SPECIALIZED FUNCTIONS 
*/
// //Specialized Overloaded Constructor - SUNMATRIX
// SparseMatrix<SundialsSparseMatrix>::SparseMatrix(int32_t nrows, int32_t ncols, int32_t maxNonzero, int storageType)
// {   
//     this->typeIndex = 0; // Type id is zero for sunmatrix
//     MATTYPE temporaryMatrix = SUNSparseMatrix(nrows,ncols,maxNonzero,storageType);
//     this->matrix = &temporaryMatrix;
// }

//Specialized SundialsSparseMatrix GETTER
template<> double SparseMatrix<SundialsSparseMatrix>::getElement(int row, int col)
{
    //IMPLEMENT ME
    return 1.0;
}

//Specialized SundialsSparseMatrix GETTER
template<> void SparseMatrix<SundialsSparseMatrix>::setElement(int row, int col, double element)
{
    //IMPLEMENT ME
}

/*
    EIGEN SPECIALIZED FUNCTIONS 
*/
//Specialized Overloaded Constructor - EIGEN
template<> SparseMatrix<EigenSparseMatrix>::SparseMatrix(EigenSparseMatrix *sparseMatrix)
{
    //IMPLEMENT ME
}

/*
Explicit instantiation allows the template header to be separate from the implementation.
*/


}

#endif