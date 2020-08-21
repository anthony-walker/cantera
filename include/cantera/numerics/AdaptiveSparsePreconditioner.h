/*
Programmer: Anthony Walker
This is the this file contains functions to adaptively precondition the sparse matrix class
*/
#ifndef ADAPTIVESPARSEPRECONDITIONER_H
#define ADAPTIVESPARSEPRECONDITIONER_H

//Cantera imports
#include "cantera/numerics/SparseMatrix.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
//Function Library for preconditioniong
namespace Cantera //Making ASP apart of Cantera namespace
{

extern "C" 
{
    
    /**
    If the user’s Jacobian-times-vector routine requires that any Jacobian-related data be preprocessed
    or evaluated, then this needs to be done in a user-supplied function of type CVLsJacTimesSetupFn,
    defined as follows:
    t is the current value of the independent variable.
    y is the current value of the dependent variable vector.
    fy is the current value of the vector f(t, y).
    user data is a pointer to user data, the same as the user data parameter passed to
    CVodeSetUserData.

     **/
    static int adaptiveMatLinSolSetup(realtype t, N_Vector y, N_Vector fy, booleantype jok, booleantype *jcurPtr, realtype gamma, void *user_data)
    {
        /*
            This is a function implemented to setup the preconditioner during integration.
        */
       
       
       return 0; //Success, return negative value for unrecoverable error or positive for recoverable error
    }

    static int adaptiveMatLinSolSolve(realtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, realtype gamma, realtype delta, int lr, void *user_data) 
    {
        /*
            This is a function implemented to solve with preconditioner during integration.
        */
       


       return 0; //Success, return negative value for unrecoverable error or positive for recoverable error
    }
}


template<class MATTYPE> void AdaptivelyPrecondition(SparseMatrix<MATTYPE> *preconditioner)
{
  /*
    This is the main preconditioner function which takes a SparseMatrix created by Eigen of the appropriate size.
  */
 
 


}



template<class MATTYPE> void TemperatureSpeciesDerivative(SparseMatrix<MATTYPE> *preconditioner)
{
  /*
    This is the main preconditioner function which takes a SparseMatrix created by Eigen of the appropriate size.
  */
}

}
#endif