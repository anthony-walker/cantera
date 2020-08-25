#include "cantera/numerics/AdaptiveSparsePreconditioner.h"
using namespace Cantera;
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
       ReactorNet *network = (ReactorNet*)user_data;
       SparseMatrix<SundialsSparseMatrix> *preconditioner = (SparseMatrix<SundialsSparseMatrix> *)network->getNetworkPreconditioner();

      if (jok)
      {
        /** jok = SUNTRUE means that the Jacobian data, if saved from the previous call to this function, can be reused
         * **/
        printf("%s","Jacobian does not need recomputed.\n");
        return 0;
      }
      else
      {
        printf("%s","Jacobian needs recomputed.\n");
        AdaptivelyPrecondition<SundialsSparseMatrix>(preconditioner,network);
        return 0; //Success, return negative value for unrecoverable error or positive for recoverable error
      }
    }

    static int adaptiveMatLinSolSolve(realtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, realtype gamma, realtype delta, int lr, void *user_data) 
    {
        /*
            This is a function implemented to solve with preconditioner during integration.
        */
       


       return 0; //Success, return negative value for unrecoverable error or positive for recoverable error
    }
}