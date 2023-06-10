import numpy as np

import cantera as ct
from . import utilities
from .utilities import allow_deprecated
import pytest


class TestMonteCarlo(utilities.CanteraTest):
    def setUp(self):
        self.solver = ct.MonteCarlo()

    def test_monte_carlo(self):
        self.solver.initialize(10, 10)
