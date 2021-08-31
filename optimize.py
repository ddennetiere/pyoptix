from scipy.stats import pearsonr
from scipy.optimize import minimize
import pandas as pd
import numpy as np


def focus(beamline, variable_oe, variable, screen, dimension="y", method=None):
    def correlation(value):
        variable_oe.__setattr__(variable, value)
        print(variable_oe.__getattribute__(variable))
        beamline.clear_impacts()
        print(beamline.radiate())
        spots = screen.get_diagram(int(beamline.active_chain[0].nrays), show_first_rays=True)
        if dimension.lower() == "xy":
            return abs(pearsonr(spots["X"], spots["Y"])[0])
        elif dimension.lower() == "x":
            return abs(pearsonr(spots["X"], spots["dX"])[0])
        elif dimension.lower() == "y":
            return abs(pearsonr(spots["Y"], spots["dY"])[0])
        else:
            raise AttributeError("Unknown dimension, should be 'x', 'y' or 'xy'")
    solution = minimize(correlation, variable_oe.__getattribute__(variable), method=method)
    print(f"Minimization success: {solution.success}, converged to {variable_oe.name}.{variable} = {solution.x}")
    if solution.success:
        return solution.x
    else:
        raise RuntimeError("Unable to reach an optimum")
