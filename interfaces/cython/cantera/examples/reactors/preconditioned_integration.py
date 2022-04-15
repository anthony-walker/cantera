# -*- coding: utf-8 -*-
"""
Ideal gas, constant-pressure, adiabatic kinetics simulation that compares preconditioned and non-preconditioned integration of nDodecane.

Requires: cantera >= 2.6.0
Keywords: combustion, reactor network, preconditioner
"""
import cantera as ct
import numpy as np
from timeit import default_timer

def compare_preconditioned_integration(preconditioner=True):
    # Use a reduced n-dodecane mechanism with PAH formation pathways
    gas = ct.Solution('nDodecane_Reitz.yaml', 'nDodecane_IG')
    # Create Reactor and set initial contents to be products of lean combustion
    gas.TP = 1000, ct.one_atm
    gas.set_equivalence_ratio(0.30, 'c12h26', 'n2:3.76, o2:1.0')
    gas.equilibrate('TP')
    r = ct.IdealGasConstPressureMoleReactor(gas)
    r.volume = 0.001
    # Create reactor network
    sim = ct.ReactorNet([r])
    # Add preconditioner
    if preconditioner:
        sim.preconditioner = ct.AdaptivePreconditioner()
    # Advance to final time
    integ_time = default_timer()
    sim.advance(10)
    integ_time = default_timer() - integ_time
    # Return time to integrate
    if preconditioner:
        print("Preconditioned Integration Time: {:f}".format(integ_time))
    else:
        print("Non-preconditioned Integration Time: {:f}".format(integ_time))
    # Get and output solver stats
    lin_stats = sim.linear_stats
    nonlin_stats = sim.nonlinear_stats
    for key in lin_stats:
        print(key, lin_stats[key])
    for key in nonlin_stats:
        print(key, nonlin_stats[key])
    print("\n")

if __name__ == "__main__":
    compare_preconditioned_integration()
    compare_preconditioned_integration(False)
