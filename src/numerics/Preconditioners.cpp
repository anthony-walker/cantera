#include "cantera/numerics/Preconditioners.h"

namespace Cantera
{
    /**
     * 
     * Preconditioner Implementation
     * 
     * **/
    double PreconditionerBase::getThreshold()
    {
        return this->threshold;
    }

    void PreconditionerBase::setThreshold(double threshold)
    {
        this->threshold = threshold;
    }

    void PreconditionerBase::setElementByThreshold(size_t row,size_t col, double element)
    {
        if (this->threshold<element)
        {
            // this->setElement(row,col,element);
        }
    }

    size_t* PreconditionerBase::getDimensions()
    {
        return &(this->dimensions[0]);
    }
    

    /**
     * 
     * Preconditioner implementations
     * 
     * **/
    Preconditioner::Preconditioner()
    {
        this->matrix =  new Eigen::SparseMatrix<double>;
    }

    Preconditioner::~Preconditioner()
    {
        delete this->matrix;
    }

    void Preconditioner::setDimensions(size_t nrows,size_t ncols, void* otherData)
    {
        this->dimensions[0] = nrows;
        this->dimensions[1] = ncols;
        this->matrix->resize(nrows,ncols);
        this->matrix->reserve(nrows*ncols);
    }

    void Preconditioner::setElement(size_t row,size_t col,double element)
    {
        this->matrix->insert(row,col)=element;
    }

    double Preconditioner::getElement(size_t row,size_t col)
    {
        return 1.0;//FIXME
    }

    Eigen::SparseMatrix<double>* Preconditioner::getMatrix()
    {
        return this->matrix;
    }

    void Preconditioner::setMatrix(Eigen::SparseMatrix<double> *sparseMat)
    {
        this->matrix=sparseMat;
    }

    void Preconditioner::solveLinearSystem(double* x, double* b)
    {

    }
}
