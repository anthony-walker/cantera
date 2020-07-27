//Cantera imports
#include "cantera/numerics/SparseMatrix.h"
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

/*
    TEMPLATE FUNCTIONS 
*/
//Overloaded Constructor
template<class MATTYPE> Cantera::SparseMatrix::SparseMatrix(MATTYPE *sparseMatrix)
{
    
}

/*
    SUNDIALS SPECIALIZED FUNCTIONS 
*/
//Specialized Overloaded Constructor - SUNMATRIX
Cantera::SparseMatrix<SundialsSparseMatrix>::SparseMatrix(SundialsSparseMatrix *sparseMatrix)
{
    //IMPLEMENT ME
}

//Specialized SundialsSparseMatrix GETTER
double Cantera::SparseMatrix<SundialsSparseMatrix>::get(int row, int col)
{
    //IMPLEMENT ME
}

//Specialized SundialsSparseMatrix GETTER
double Cantera::SparseMatrix<SundialsSparseMatrix>::set(int row, int col, double element)
{
    //IMPLEMENT ME
}

/*
    EIGEN SPECIALIZED FUNCTIONS 
*/
//Specialized Overloaded Constructor - EIGEN
Cantera::SparseMatrix<EigenSparseMatrix>::SparseMatrix(EigenSparseMatrix *sparseMatrix)
{
    //IMPLEMENT ME
}

