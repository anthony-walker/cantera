/*
Programmer: Anthony Walker

This is the SparseMatrix class which serves as a wrapper for Sundials and eventually Eigen Sparse Matrices.
*/
#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H


namespace Cantera
{
template<class MATTYPE> class SparseMatrix
{
    private:
        MATTYPE matrix; //Matrix stored by sparse matrix
    public:
        SparseMatrix(); //default constructor
        SparseMatrix(MATTYPE *sparseMatrix); //overloaded
        double get(int row, int col); //get element
        void set(int row, int col, double element);//set element
    };
}

#endif