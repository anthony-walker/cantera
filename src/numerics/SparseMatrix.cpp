/*
Programmer: Anthony Walker

This is the SparseMatrix class which serves as a wrapper for Sundials and eventually Eigen Sparse Matrices.
*/
#include "cantera/numerics/SparseMatrix.h"

using namespace Cantera; //Need this because SparseMatrix.h is defined in the Cantera namespace
/*
    SUNDIALS SPECIALIZED FUNCTIONS 

    --fully specialized functions have to be placed inside of cpp or causes multiple definition error.
*/

//Specialized SundialsSparseMatrix SETTER
template<> void SparseMatrix<SundialsSparseMatrix>::setElement(size_t row, size_t col, double element)
{   
    
}

//Specialized SundialsSparseMatrix GETTER
template<> double SparseMatrix<SundialsSparseMatrix>::getElement(size_t row, size_t col)
{
    //IMPLEMENT ME
    return 1.0;
}


/*
    EIGEN SPECIALIZED FUNCTIONS 
*/

//Specialized EigenSparseMatrix SETTER
template<> void SparseMatrix<EigenSparseMatrix>::setElement(size_t row, size_t col, double element)
{   
    this->matrix->insert(row,col) = element;
}


//Specialized EigenSparseMatrix GETTER
template<> double SparseMatrix<EigenSparseMatrix>::getElement(size_t row, size_t col)
{
    return 1.0;
}

