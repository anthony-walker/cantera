# This file is part of Cantera. See License.txt in the top-level directory or
# at https://cantera.org/license.txt for license and copyright information.

cdef class PreconditionerBase:
    """
    Common base class for preconditioners.
    """
    precon_type = "PreconditionerBase"
    precon_prob_type = "DENSE + NOJAC"

    def __cinit__(self, *args, **kwargs):
        self.pbase = newPreconditioner(stringify(self.precon_type))

cdef class AdaptivePreconditioner(PreconditionerBase):
    precon_type = "Adaptive"
    precon_prob_type = "GMRES"

    def __cinit__(self, *args, **kwargs):
        self.preconditioner = <CxxAdaptivePreconditioner*>(self.pbase)

    property threshold:
        """
        Property for setting behavior of preconditioner pruning

        Update the threshold to a desired value as:
            >>> precon.threshold = 1e-8

        Default is 1e-8.
        """
        def __get__(self):
            return self.preconditioner.threshold()

        def __set__(self, val):
            self.preconditioner.setThreshold(val)

    property perturb_const:
        """
        Property for setting the perturbation constant used in finite difference calculations for temperature.

        Update the threshold to a desired value as:
            >>> precon.perturb_const = 1e-8

        Default is DBL_EPSILON.
        """
        def __get__(self):
            return self.preconditioner.perturbation()

        def __set__(self, val):
            self.preconditioner.setPerturbation(val)

    property ilut_fill_factor:
        """
        Property setting the linear solvers fill factor.

        During factorization, after row elimination, only some of the largest elements in the L and U in addition to the diagonal element are kept. The number of elements kept is computed from the fill factor (a ratio) relative to the initial number of nonzero elements.

        Update the ILUT fill factor to a desired value as:
            >>> precon.ilut_fill_factor = 2

        Default is the state size divided by 4.
        """
        def __set__(self, val):
            self.preconditioner.setIlutFillFactor(val)

        def __get__(self):
            return self.preconditioner.ilutFillFactor()

    property ilut_drop_tol:
        """
        Property setting the linear solvers drop tolerance.

        During factorization any element below the product of the drop tolerance and average magnitude is dropped.

        Update the ILUT drop tolerance to a desired value as:
            >>> precon.ilut_drop_tol = 1e-10

        Default is 1e-10.
        """
        def __set__(self, val):
            self.preconditioner.setIlutDropTol(val)

        def __get__(self):
            return self.preconditioner.ilutDropTol()

    property derv_settings:
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
            self.preconditioner.preconSettings(settings)
            return anymap_to_dict(settings)

        def __set__(self, settings):
            self.preconditioner.setPreconSettings(dict_to_anymap(settings))

    def print_contents(self):
        self.preconditioner.printPreconditioner()