cdef class PreconditionerBase:
    """
    Common base class for reactors and reservoirs.
    """
    precon_type = "PreconditionerBase"
    def __cinit__(self, *args, **kwargs):
        self.pbase = newPreconditioner(stringify(self.precon_type))

    def __dealloc__(self):
        del self.pbase

cdef class AdaptivePreconditioner(PreconditionerBase):
    precon_type = "AdaptivePreconditioner"

    def __cinit__(self, *args, **kwargs):
        self.preconditioner = <CxxAdaptivePreconditioner*>(self.pbase)

    def get_threshold(self):
        return self.preconditioner.getThreshold()

    def set_threshold(self, val):
        self.preconditioner.setThreshold(val)

    def set_fill_factor(self, val=-1):
        self.preconditioner.setFillFactorILUT(val)

    def set_drop_tol(self, val=1e-10):
        self.preconditioner.setDropTolILUT(val)

    def print_contents(self):
        self.preconditioner.printPreconditioner()
