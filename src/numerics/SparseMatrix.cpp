#include "cantera/numerics/SparseMatrix.h"
#include "cantera/numerics/AdaptiveSparsePreconditioner.h"
using namespace Cantera; //Need this because SparseMatrix.h is defined in the Cantera namespace
/*
    SUNDIALS SPECIALIZED FUNCTIONS 

    --fully specialized functions have to be placed inside of cpp or causes multiple definition error.
*/
//Overloaded Constructor
template<> SparseMatrix<SundialsSparseMatrix>::SparseMatrix(int nrows, int ncols, int maxNonZero)
{   
    SundialsSparseMatrix temporaryMatrix = SUNSparseMatrix(nrows,ncols,maxNonZero,CSC_MAT);
    this->matrix = &temporaryMatrix;
    this->nrows = nrows;
    this->ncols = ncols;
}

//Specialized SundialsSparseMatrix SETTER
template<> void SparseMatrix<SundialsSparseMatrix>::setElement(int row, int col, double element)
{   
    printf("%s","Setting by threshold\n");
}

template<> void SparseMatrix<SundialsSparseMatrix>::setElement(int row, int col, double element,bool checkThreshold)
{   
    if (element > this->threshold)
    {
        printf("%s","Setting by threshold\n");
    }
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
    this->nrows = nrows;
    this->ncols = ncols;
}

//Specialized EigenSparseMatrix SETTER
template<> void SparseMatrix<EigenSparseMatrix>::setElement(int row, int col, double element)
{   
    this->matrix->insert(row,col) = element;
}

//Specialized EigenSparseMatrix SETTER
template<> void SparseMatrix<EigenSparseMatrix>::setElement(int row, int col, double element,bool checkThreshold)
{   
    if (element > this->threshold)
    {
        this->matrix->insert(row,col) = element;
    }
}

//Specialized EigenSparseMatrix GETTER
template<> double SparseMatrix<EigenSparseMatrix>::getElement(int row, int col)
{
    // int *innerIndexList = this->matrix->innerIndexPtr();
    // int *outerIndexList = this->matrix->outerIndexPtr();
    // int innerLen = sizeof(&innerIndexList)/sizeof(innerIndexList[0]);
    // int outerLen = sizeof(&outerIndexList)/sizeof(outerIndexList[0]);
    // for (int i = 0; i < innerLen; i++)
    // {
        
    //     // std::cout << innerIndexList[i] << std::endl;
    //     // if (innerIndexList[i] == row)
    //     // {
    //     //     for (int j = 0; j < outerLen; j++)
    //     //     {
    //     //         if (outerIndexList[j] == col)
    //     //         {
                    
    //     //         }
    //     //     }
            
    //     // }
    // }
    
    return 1.0;
}

