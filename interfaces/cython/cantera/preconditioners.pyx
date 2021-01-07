# distutils: language = c++
from libcpp.map cimport map
from libcpp.string  cimport string

# Implementation of python wrapper's for preconditioner functions
cdef extern from "cantera/numerics/Preconditioners.h" namespace "Cantera":
    cdef cppclass PreconditionerBase:
        PreconditionerBase() except +


# ctypedef map<string,unsigned long> StateMap
# ctypedef void (*AdaptiveFunction)(PreconditionerBase *preconditioner,Reactor* reactor, double** inputs,StateMap indexMap, string key)
# ctypedef map<string,vector<AdaptiveFunction>> FunctionMap

# cdef extern from "cantera/numerics/Preconditioners.h" namespace "Cantera::AMP":
#     cdef cppclass AdaptivePreconditioner:
#         AdaptivePreconditioner() except +
#         # AdaptivePreconditioner(FunctionMap) except +
#         double threshold
#         double getThreshold()

# cdef class PyPreconditionerBase:
#     cdef PreconditionerBase preconditioner
#     def __cinit__(self):
#         self.preconditioner = PreconditionerBase()

    # def __del__(self):
    #     del self.preconditioner

# cdef extern from "Preconditioners.h" namespace Cantera::AMP:
#     cppclass AdaptivePreconditioner:
#         AdaptivePreconditioner()
#         AdaptivePrecondition(FunctionMap fmap)

