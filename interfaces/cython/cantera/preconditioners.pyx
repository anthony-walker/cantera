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

    property threshold:
        """
        Property for setting behavior of preconditioner pruning

        Update the threshold to a desired value as:

            >>> precon.threshold = 1e-8
        """
        def __get__(self):
            return self.preconditioner.getThreshold()

        def __set__(self, val):
            self.preconditioner.setThreshold(val)

    property perturb_const:
        """
        Property for setting the perturbation constant used in finite difference calculations for temperature

        Update the threshold to a desired value as:

            >>> precon.perturb_const = 1e-8

        Default is __DBL_EPSILON__
        """
        def __get__(self):
            return self.preconditioner.getPerturbationConst()

        def __set__(self, val):
            self.preconditioner.setPerturbationConst(val)

    property ilut_fill_factor:
        """
        Property setting the linear solvers fill factor.

        This must be called after initialization of the network because a default fill factor is set during network initialization.

        Update the ILUT fill factor to a desired value as:

            >>> precon.ilut_fill_factor = 2
        """
        def __set__(self, val):
            self.preconditioner.setFillFactorILUT(val)

        def __get__(self):
            pass

    property ilut_drop_tol:
        """
        Property setting the linear solvers drop tolerance.

        This must be called after initialization of the network because a default is set during network initialization.

        Update the ILUT drop tolerance to a desired value as:

            >>> precon.ilut_drop_tol = 1e-10
        """
        def __set__(self, val):
            self.preconditioner.setDropTolILUT(val)

        def __get__(self):
            pass

    property precon_derv_settings:
        """
        Property setting behavior of preconditioner derivative evaluation.

        the following keyword/value pairs are supported:

        -  ``skip-third-bodies`` (boolean) ... if `False` (default), third body
           concentrations are considered for the evaluation of derivatives

        -  ``skip-falloff`` (boolean) ... if `True` (default), third-body effects
           on reaction rates are not considered.

        -  ``rtol-delta`` (double) ... relative tolerance used to perturb properties
           when calculating numerical derivatives. The default value is 1e-8.

        Derivative settings are updated using a dictionary::

            >>> precon.derivative_settings = {"skip-falloff": True}

        Passing an empty dictionary will reset all values to their defaults.
        """
        def __get__(self):
            cdef CxxAnyMap settings
            self.preconditioner.getPreconditionerDerivativeSettings(settings)
            return anymap_to_dict(settings)

        def __set__(self, settings):
            self.preconditioner.setPreconditionerDerivativeSettings(dict_to_anymap(settings))

    def print_contents(self):
        self.preconditioner.printPreconditioner()
