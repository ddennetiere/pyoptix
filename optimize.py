from scipy.stats import pearsonr
from scipy.optimize import minimize
import pandas as pd
import numpy as np
from ipywidgets import interact, Layout
import ipywidgets as widgets
from bokeh.models.widgets import Slider, TextInput
from bokeh.application import Application
from bokeh.application.handlers import FunctionHandler
from bokeh.layouts import row, column
from bokeh.io import  show
from bokeh.io import push_notebook


def slider_optimizer(variable_oe=None, variable="", variable_bounds=(), variable_step=0.1, screen=None,
                     wavelength=6e-9, nrays=None, display="yyp", light_spd=False):
    """
    Prints out the spot diagram kind asked in display on the surface of 'screen' and a slider object which when moved
    will update the spot diagram.
    The slider sets the value of the variable 'variable' of the optical element 'variable_oe' between
    'variable_bounds' by increment of 'variable_step'.
    All propagated rays are at wavelength 'wavelength'.

    :param variable_oe: Optical element whose parameter must be varied
    :type variable_oe: any class inheriting pyoptix.OpticalElement
    :param variable: name of the parameter to be varied
    :type variable: str
    :param wavelength: wavelength at which beamline must be aligned in m.
    :type wavelength: float
    :param nrays: Number of rays to propagate
    :type nrays: int
    :param screen: recording surface where beam must be focused
    :type screen: any class inheriting pyoptix.OpticalElement
    :param variable_bounds: Bounds of the variable and thus the slider
    :type variable_bounds: tuple fo float
    :param variable_step: Minimum increment for the slider
    :type variable_step: float
    :param display: Indicates which representation is to be shown, can be "xy", "xxp" or "yyp"
    :type display: str
    :param light_spd: set to True for quick monochromatic rendering of the scatter plot
    :type light_spd: bool
    :return: None
    :rtype: NoneType
    """
    assert display != "all"
    assert " " not in display
    if nrays is None:
        nrays = int(screen.beamline.active_chain[0].nrays)
    else:
        screen.beamline.active_chain[0].nrays = int(nrays)
    screen.beamline.clear_impacts(clear_source=True)
    screen.beamline.align(wavelength)
    screen.beamline.generate(wavelength)
    screen.beamline.radiate()
    datasource, handles = screen.show_diagram(display=display, light_yyp=light_spd, light_xy=light_spd,
                                              light_xxp=light_spd, show_spd=False)
    v0 = variable_oe.__getattribute__(variable)

    # def f(x):
    def update_data(attrname, old, new):
        screen.beamline.clear_impacts(clear_source=False)
        variable_oe.__setattr__(variable, value_slider.value)
        screen.beamline.align(wavelength)
        screen.beamline.radiate()
        spd = screen.get_diagram(nrays)
        datasource[display].data.update(spd)
        # push_notebook(handle=handles[0])

    def update_title(attrname, old, new):
        handles[0][0].title.text = text.value

    # interact(f, x=widgets.FloatSlider(min=variable_bounds[0], max=variable_bounds[1], step=variable_step,
    #                                   value=v0, continuous_update=False, layout=Layout(width='90%'),
    #                                   description=f"{variable}", readout_format='.3e'))
    value_slider = Slider(start=variable_bounds[0], end=variable_bounds[1], step=variable_step,
                          value=v0, sizing_mode="stretch_width")
    text = TextInput(title="title", value='my sine wave')

    value_slider.on_change('value_throttled', update_data)
    layout = column(handles[0], value_slider, text)

    def modify_doc(doc):
        doc.add_root(row(layout, width=800))
        doc.title = "Sliders"
        text.on_change('value', update_title)

    handler = FunctionHandler(modify_doc)
    app = Application(handler)
    return app
    # show(app)


def focus(beamline, variable_oe, variable, wavelength, screen, dimension="y", nrays=None, method="Nelder-Mead",
          show_progress=False, tol=1e-3, options=None):
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
    :param method: Method to be used for minimisation, default "Nedler-Mead"
    :type method: str
    :param show_progress: If True, each iteration will print the current variable value and value of the function
        to be minimized
    :type show_progress: bool
    :param tol: Tolerance for the optimizer, See scipy.optimize.minimize documentation
    :type tol: float
    :param options: Method-specific options, see scipy.optimize.minimize for details
    :type options: dict
    :return: optimal value of the variable to achieve focusing
    :rtype: float
    :raises RuntimeError: if variable has no effect or minimum cannot be reached with asked tolerance
    """
    if options is None:
        options = {}
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
    beamline.clear_impacts(clear_source=True)
    beamline.generate(wavelength)

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

    if "fatol" not in options.keys() and method == "Nedler-Mead":
        options["fatol"] = tol
    if "xatol" not in options.keys() and method == "Nedler-Mead":
        options["xatol"] = tol
    solution = minimize(correlation, variable_oe.__getattribute__(variable), method=method, tol=tol, bounds=bounds,
                        options=options)
    print(f"Minimization success: {solution.success}, converged to {variable_oe.name}.{variable} = {solution.x}")
    beamline.active_chain[0].nrays = int(old_nrays)
    screen.next = old_link
    if solution.success:
        return solution.x[0]
    else:
        raise RuntimeError("Unable to reach an optimum")


def find_focus(beamline, screen, wavelength, dimension="y", nrays=None, method="Nelder-Mead",
               show_progress=False, tol=1e-3, options=None, adjust_distance=True):
    """
    Function to be called for minimizing a focused spot diagram placed at "screen" by varying it distance to the
    previous oe in the "beamline" beamline. One can either try to focus horizontally, vertically or in both dimension by
    specifying the "dimension" parameter with respectively "x", "y" or "xy". Number of rays for the computation can be
    specified or the nrays parameter of the beamline source will be used. Minimization method can be specified, see
    scipy.optimize.minimize documentation for available algorithms. If adjust_distance is set to True, when convergence
    is reached, the distance between screen and previous OE is set to the calculated optimal value.

    :param beamline: Beamline along which to propagate rays
    :type beamline: pyoptix.Beamline
    :param screen:  recording surface where beam must be focused
    :type screen: any class inheriting pyoptix.OpticalElement
    :param wavelength: wavelength at which beamline must be aligned in m.
    :type wavelength: float
    :param dimension: dimension along which focusing is desired. Must be "x", "y" or "xy"
    :type dimension: str
    :param nrays: Number of rays to propagate
    :type nrays: int
    :param method: Method to be used for minimisation, default "Nedler-Mead"
    :type method: str
    :param show_progress: If True, each iteration will print the current variable value and value of the function
        to be minimized
    :type show_progress: bool
    :param tol: Tolerance for the optimizer, See scipy.optimize.minimize documentation
    :type tol: float
    :param options: Method-specific options, see scipy.optimize.minimize for details
    :type options: dict
    :param adjust_distance: if True, sets the screen at the optimal distance
    :type adjust_distance: bool
    :return: Optimal distance from the screen to get focus in given dimension
    :rtype: float
    :raises RuntimeError: if distance has no effect or minimum cannot be reached with asked tolerance
    """

    if options is None:
        options = {}
    old_nrays = int(beamline.active_chain[0].nrays)
    if nrays is None:
        nrays = int(beamline.active_chain[0].nrays)
    else:
        beamline.active_chain[0].nrays = int(nrays)
    bounds = None
    if screen.get_whole_parameter("distance")["bounds"] != (0, 0):
        bounds = screen.get_whole_parameter("distance")["bounds"]
    beamline.clear_impacts(clear_source=True)
    beamline.generate(wavelength)
    beamline.align(wavelength)
    beamline.radiate()

    def correlation(value):
        spots = screen.get_diagram(nrays, show_first_rays=False, distance_from_oe=value)
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

    if "fatol" not in options.keys() and method == "Nedler-Mead":
        options["fatol"] = tol
    if "xatol" not in options.keys() and method == "Nedler-Mead":
        options["xatol"] = tol
    solution = minimize(correlation, 0, method=method, tol=tol, bounds=bounds,
                        options=options)
    print(f"Minimization success: {solution.success}, converged to a distance of {solution.x}")
    beamline.active_chain[0].nrays = int(old_nrays)
    if solution.success:
        if adjust_distance:
            screen.distance_from_previous = solution.x[0]
        return solution.x[0]
    else:
        raise RuntimeError("Unable to reach an optimum")
