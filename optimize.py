from scipy.stats import pearsonr
from scipy.optimize import minimize
import pandas as pd
import numpy as np


def focus(beamline, variable_oe, variable, wavelength, screen, dimension="y", nrays=None, method=None,
          show_progress=False, tol=1e-3):
    """
    Function to be called for minimizing a focused spot diagram placed at "screen" by varying the parameter "variable"
    of the optical element "variable_oe" of the "beamline" beamline. Each iteration is realigned at wavelength
    "wavelength". One can either try to focus horizontally, vertically or in both dimension by specifying the
    "dimension" parameter with respectively "x", "y" or "xy". Number of rays for the computation can be specified or
    the nrays parameter of the beamline source will be used. Minimization method can be specified, see
    scipy.optimize.minimize documentation for available algorithms.
    :param beamline: Beamline along which to propagate rays
    :type beamline: pyoptix.Beamline
    :param variable_oe: Optical element whose parameter must be varied
    :type variable_oe: any class inheriting pyoptix.OpticalElement
    :param variable: name of the parameter to be varied
    :type variable: str
    :param wavelength: wavelength at which beamline must be aligned in m.
    :type wavelength: float
    :param screen: recording surface where beam must be focused
    :type screen: any class inheriting pyoptix.OpticalElement
    :param dimension: dimension along which focusing is desired. Must be "x", "y" or "xy"
    :type dimension: str
    :param nrays: Number of rays to propagate
    :type nrays: int
    :param method: Method to be used for minimisation
    :type method: str
    :param show_progress: If True, each iteration will print the current variable value and value of the function
                          to be minimized
    :type show_progress: bool
    :param tol: Tolerance for the optimizer, See scipy.optimize.minimize documentation
    :type tol: float
    :return: optimal value of the variable to achieve focusing
    :rtype: float
    :raises RuntimeError: if variable has no effect or minimum cannot be reached with asked tolerance
    """
    old_nrays = int(beamline.active_chain[0].nrays)
    old_link = screen.next
    screen.next = None
    if nrays is None:
        nrays = int(beamline.active_chain[0].nrays)
    else:
        beamline.active_chain[0].nrays = int(nrays)
    bounds = None
    if variable_oe.get_whole_parameter(variable)["bounds"] != (0, 0):
        bounds = variable_oe.get_whole_parameter(variable)["bounds"]

    def correlation(value):
        variable_oe.__setattr__(variable, value)
        beamline.clear_impacts()
        beamline.align(wavelength)
        beamline.radiate()
        spots = screen.get_diagram(nrays, show_first_rays=False)
        if dimension.lower() == "xy":
            ret = np.std(spots["X"]**2 + spots["Y"]**2)
        elif dimension.lower() == "x":
            ret = abs(pearsonr(spots["X"], spots["dX"])[0])
        elif dimension.lower() == "y":
            ret = abs(pearsonr(spots["Y"], spots["dY"])[0])
        else:
            raise AttributeError("Unknown dimension, should be 'x', 'y' or 'xy'")
        if show_progress:
            print(value, ret)
        return ret
    solution = minimize(correlation, variable_oe.__getattribute__(variable), method=method, tol=tol, bounds=bounds)
    print(f"Minimization success: {solution.success}, converged to {variable_oe.name}.{variable} = {solution.x}")
    beamline.active_chain[0].nrays = int(old_nrays)
    screen.next = old_link
    if solution.success:
        return solution.x[0]
    else:
        raise RuntimeError("Unable to reach an optimum")
