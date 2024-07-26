from scipy.stats import pearsonr
from scipy.optimize import minimize
import pandas as pd
import numpy as np
from ipywidgets import interact, Layout
import ipywidgets as widgets
from IPython.display import display as ipy_display
from bokeh.models.widgets import Slider, TextInput
from bokeh.application import Application
from bokeh.application.handlers import FunctionHandler
from bokeh.layouts import row, column
import plotly.graph_objs as go
from bokeh.io import show
from bokeh.io import push_notebook


class AttrDict(dict):
    """
    Magic class that let's the user access a dictionnary keys as attribute
    """
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


fig_keys = AttrDict(dict(xy=AttrDict({"x": "X", "y": "Y", "xunit": "m", "yunit": "m"}),
                         xxp=AttrDict({"x": "X", "y": "dX", "xunit": "m", "yunit": "rad"}),
                         yyp=AttrDict({"x": "Y", "y": "dY", "xunit": "m", "yunit": "rad"})))


def slider_optimizer_bokeh(variable_oe=None, variable="", variable_bounds=(), variable_step=0.1, screen=None,
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
    print(screen.beamline.active_chain[0].nrays)
    screen.beamline.clear_impacts(clear_source=True)
    screen.beamline.align(wavelength, wavelength)
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
        spd = screen.get_diagram()
        datasource[display].data.update(spd)
        # push_notebook(handle=handles[0])

    # interact(f, x=widgets.FloatSlider(min=variable_bounds[0], max=variable_bounds[1], step=variable_step,
    #                                   value=v0, continuous_update=False, layout=Layout(width='90%'),
    #                                   description=f"{variable}", readout_format='.3e'))
    value_slider = Slider(start=variable_bounds[0], end=variable_bounds[1], step=variable_step,
                          value=v0, sizing_mode="stretch_width")

    value_slider.on_change('value_throttled', update_data)
    layout = column(handles[0], value_slider)

    def modify_doc(doc):
        doc.add_root(row(layout, width=800))
        doc.title = "Sliders"

    handler = FunctionHandler(modify_doc)
    app = Application(handler)
    return app
    # show(app)


def slider_optimizer_plotly(variable_oe=None, variable="", variable_bounds=(), variable_step=0.1, screen=None,
                            wavelength=6e-9, nrays=1000, display="yyp", inverse_value=False, run_func=None):
    v0 = variable_oe.__getattribute__(variable)
    if inverse_value:
        v0 = 1/v0
    layout = widgets.Layout(width='auto', height='40px')
    slider = widgets.FloatSlider(
        value=v0,
        min=variable_bounds[0],
        max=variable_bounds[1],
        step=variable_step,
        description=f"{variable_oe.name} {variable}:",
        disabled=False,
        continuous_update=True,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',
        layout=layout
    )
    output = widgets.Output()

    screen.beamline.active_chain[0].nrays = nrays
    if run_func:
        run_func()
    else:
        screen.beamline.clear_impacts(clear_source=False)
        screen.beamline.align(wavelength)
        screen.beamline.radiate()
    df = screen.get_diagram(0)
    fig = go.FigureWidget(go.Scatter(x=df[fig_keys[display].x], y=df[fig_keys[display].y], mode="markers"))
    fig.layout.xaxis.title = f"{fig_keys[display].x} ({fig_keys[display].xunit})"
    fig.layout.yaxis.title = f"{fig_keys[display].y} ({fig_keys[display].yunit})"
    points = fig.data[0]
    points = fig.data[0]

    def on_slider_moved(_):
        if inverse_value:
            value = 1/slider.value
        else:
            value = slider.value
        variable_oe.__setattr__(variable, value)
        if run_func:
            run_func()
        else:
            screen.beamline.clear_impacts(clear_source=False)
            screen.beamline.align(wavelength)
            screen.beamline.radiate()
        with output:
            points.x, points.y = (screen.get_diagram(0)[fig_keys[display].x],
                                  screen.get_diagram(0)[fig_keys[display].y])

    slider.observe(on_slider_moved, names="value")
    ipy_display(fig, slider, output)


def multi_slider_optimizer_plotly(variable_oes=None, variables=None, variable_bounds=None, variable_steps=None,
                                  screen=None,
                                  wavelength=6e-9, nrays=1000, display="yyp", inverse_values=None, run_func=None):
    """
    Create interactive Plotly sliders for optimizing multiple beamline variables.

    This function sets up Jupyter widget sliders to adjust specified variables in multiple optical
    elements (OEs) and dynamically updates Plotly scatter plots based on the changes. It is used for
    optimizing multiple beamline parameters interactively.

    Parameters
    ----------
    variable_oes : list of objects, optional
        The optical elements (OEs) containing the variables to be optimized. Default is None.
    variables : list of str, optional
        The names of the variables in `variable_oes` to be optimized. Default is None.
    variable_bounds : list of tuples, optional
        The bounds (min, max) for the sliders corresponding to `variables`. Default is None.
    variable_steps : list of floats, optional
        The step sizes for the sliders. Default is None.
    screen : object, optional
        The screen object used for beamline simulations and plotting. Default is None.
    wavelength : float, optional
        The wavelength used for beamline alignment and radiation, in meters. Default is 6e-9.
    nrays : int, optional
        The number of rays used in the beamline simulation. Default is 1000.
    display : str, optional
        The keys for the x and y data to be displayed in the Plotly scatter plot. Default is "yyp".
    inverse_values : list of bool, optional
        If True for a variable, the slider value is used as the inverse of the variable. Default is None.
    run_func : function, optional
        A custom function to run for each slider update instead of the default beamline simulation
        process. Default is None.

    Returns
    -------
    None
        Displays interactive sliders and dynamically updating Plotly scatter plots in a Jupyter
        notebook.

    Notes
    -----
    - The sliders update the specified `variables` in `variable_oes` within the provided bounds.
    - If `inverse_values` is True for a variable, the slider value is taken as the inverse of the variable.
    - The function either uses a custom `run_func` or performs default beamline simulation steps
      (clear impacts, align, radiate) to update the plots.
    - The scatter plots update dynamically based on the changes in `variables`, displaying the
      results of the beamline simulation or custom function output.
    """
    if variable_oes is None or variables is None or variable_bounds is None or variable_steps is None:
        raise ValueError("variable_oes, variables, variable_bounds, and variable_steps must be provided")

    if inverse_values is None:
        inverse_values = [False] * len(variables)

    sliders = []

    for i, (variable_oe, variable, bounds, step, inverse_value) in enumerate(
            zip(variable_oes, variables, variable_bounds, variable_steps, inverse_values)):
        v0 = variable_oe.__getattribute__(variable)
        if inverse_value:
            v0 = 1 / v0
        layout = widgets.Layout(width='auto', height='40px')
        slider = widgets.FloatSlider(
            value=v0,
            min=bounds[0],
            max=bounds[1],
            step=step,
            description=f"{variable_oe.name} {variable}:",
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='.1f',
            layout=layout
        )
        sliders.append(slider)

    screen.beamline.active_chain[0].nrays = nrays
    if run_func:
        run_func()
    else:
        screen.beamline.clear_impacts(clear_source=False)
        screen.beamline.align(wavelength)
        screen.beamline.radiate()

    df = screen.get_diagram(0)
    fig = go.FigureWidget(go.Scatter(x=df[fig_keys[display].x], y=df[fig_keys[display].y], mode="markers"))
    fig.layout.xaxis.title = f"{fig_keys[display].x} ({fig_keys[display].xunit})"
    fig.layout.yaxis.title = f"{fig_keys[display].y} ({fig_keys[display].yunit})"
    points = fig.data[0]

    def on_slider_moved(change):
        for slider, variable_oe, variable, inverse_value in zip(sliders, variable_oes, variables, inverse_values):
            if inverse_value:
                value = 1 / slider.value
            else:
                value = slider.value
            variable_oe.__setattr__(variable, value)

        if run_func:
            run_func()
        else:
            screen.beamline.clear_impacts(clear_source=False)
            screen.beamline.align(wavelength)
            screen.beamline.radiate()

        points.x, points.y = (screen.get_diagram(0)[fig_keys[display].x],
                              screen.get_diagram(0)[fig_keys[display].y])

    for slider in sliders:
        slider.observe(on_slider_moved, names="value")

    widgets_list = [fig] + sliders
    ipy_display(*widgets_list)


def focus(beamline, variable_oe, variable, wavelength, screen, dimension="y", nrays=None, method="Nelder-Mead",
          show_progress=False, tol=1e-3, options=None, verbose=1):
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
    :param verbose: Verbose level control. 0 is silent
    :type verbose: int
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
        try:
            assert spots["Y"].std() != 0
            assert spots["X"].std() != 0
            assert spots["dX"].std() != 0
            assert spots["dY"].std() != 0
        except AssertionError:
            print(beamline.active_chain)
            for oe in beamline.active_chain[1:]:
                diag = oe.get_diagram()
                print(oe.name," : ", diag["Y"].std(), " -> ", (oe.next.name if oe.next is not None else None))
            raise AssertionError("Unexploitable rays")
        if dimension.lower() == "xy":
            ret = np.std(spots["X"]**2 + spots["Y"]**2)
        elif dimension.lower() == "x":
            # ret = abs(pearsonr(spots["X"], spots["dX"])[0])
            ret = spots["X"].std()
        elif dimension.lower() == "y":
            # ret = abs(pearsonr(spots["Y"], spots["dY"])[0])
            ret = spots["X"].std()
        else:
            raise AttributeError("Unknown dimension, should be 'x', 'y' or 'xy'")
        if show_progress:
            print(variable_oe.__getattribute__(variable), ret)
        return ret

    if "fatol" not in options and method == "Nedler-Mead":
        options["fatol"] = tol
    if "xatol" not in options and method == "Nedler-Mead":
        options["xatol"] = tol
    solution = minimize(correlation, variable_oe.__getattribute__(variable), method=method, tol=tol, bounds=bounds,
                        options=options)
    if verbose:
        print(f"Minimization success: {solution.success}, converged to {variable_oe.name}.{variable} = {solution.x}")
    beamline.active_chain[0].nrays = int(old_nrays)
    screen.next = old_link
    if solution.success:
        return solution.x[0]
    else:
        raise RuntimeError("Unable to reach an optimum")


def find_focus(beamline, screen, wavelength, dimension="y", nrays=None, method="Nelder-Mead",
               show_progress=False, tol=1e-3, options=None, adjust_distance=True, verbose=1):
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
    :param verbose: Verbose level control. 0 is silent
    :type verbose: int
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
        ret = correlation_on_screen(screen, nrays, dimension=dimension, distance_to_screen=value)
        if show_progress:
            print(value, ret)
        return ret

    if "fatol" not in options and method == "Nedler-Mead":
        options["fatol"] = tol
    if "xatol" not in options and method == "Nedler-Mead":
        options["xatol"] = tol
    solution = minimize(correlation, 0, method=method, tol=tol, bounds=bounds,
                        options=options)
    if verbose:
        print(f"Minimization success: {solution.success}, converged to a distance of {solution.x}")
    beamline.active_chain[0].nrays = int(old_nrays)
    if solution.success:
        if adjust_distance:
            screen.distance_from_previous += solution.x[0]
        return solution.x[0]
    else:
        raise RuntimeError("Unable to reach an optimum")


def correlation_on_screen(screen, nrays, dimension="y", distance_to_screen=0):
    spots = screen.get_diagram(nrays, show_first_rays=False, distance_from_oe=distance_to_screen)
    if dimension.lower() == "xy":
        ret = np.std(spots["X"]**2 + spots["Y"]**2)
    elif dimension.lower() == "x":
        ret = abs(pearsonr(spots["X"], spots["dX"])[0])
    elif dimension.lower() == "y":
        ret = abs(pearsonr(spots["Y"], spots["dY"])[0])
    else:
        raise AttributeError("Unknown dimension, should be 'x', 'y' or 'xy'")
    return ret


def custom_optimizer(beamline, screen, wavelengths, oes, attributes, dimension="y", nrays=None, method="Nelder-Mead",
                     show_progress=False, tol=1e-3, options=None, move_screen=False, norm=2, special_align=lambda:None,
                     verbose=1):
    """
    Optimize the optical system based on the figure of merit (FOM) using the specified optimization method.

    Parameters
    ----------
    beamline : pyoptix.Beamline
        Beamline object containing the optical elements and their properties
    screen : pyoptix.OpticalElement
        Screen object containing the position and orientation of the screen in the beamline
    wavelengths : float or list of floats
        Wavelengths to be used in the optimization
    oes : object or list of objects
        Optical element(s) to be optimized, length has to match the number of attibutes.
    attributes : str or list of str
        Attribute(s) of the optical element(s) to be optimized, length has to match the number or optical elements
    dimension : str, optional
        Dimension of the screen diagram to be optimized (default is "y")
    nrays : int, optional
        Number of rays used to generate the diagram (default is None)
    method : str, optional
        Optimization method to be used (default is "Nelder-Mead")
    show_progress : bool, optional
        If True, print the current value of the FOM during optimization (default is False)
    tol : float, optional
        Tolerance for termination (default is 1e-3)
    options : dict, optional
        Options for the optimization method (default is None)
    move_screen : bool, optional
        If True, move the screen to the focal position during optimization (default is False)
    norm : int, optional
        Order of the norm used to calculate the FOM (default is 2)
    special_align : callable, optional
        function to be called before running simulation, can be used for alignment
    verbose : int, optional
        If 1, print the optimized values of the attributes (default is 1)

    Returns
    -------
    float
        Optimized value of the first attribute of the first optical element

    Raises
    ------
    RuntimeError
        If the optimization fails to converge to a minimum

    Notes
    -----
    The function uses the figure of merit (FOM) defined as the norm of the standard deviations of the screen
    diagrams for the specified wavelengths. The optimization method minimizes the FOM with respect to the
    specified attribute(s) of the optical element(s).
    """

    old_nrays = beamline.active_chain[0].nrays
    old_link = screen.next
    screen.next = None
    if nrays is not None:
        beamline.active_chain[0].nrays = int(nrays)
    if type(oes) == list:
        assert len(oes) == len(attributes), "Number of optical element must match number of attributes"
    else:
        oes = [oes]
        attributes = [attributes]
    if type(wavelengths) == float:
        wavelengths = [wavelengths]

    def FOM(attributes_values):
        for oe, attrib, attrib_value in zip(oes, attributes, attributes_values):
            attributes_values_0.append(oe.__setattr__(attrib, attrib_value))
        contributions = []
        for wavelength in wavelengths:
            beamline.clear_impacts(clear_source=True)
            beamline.align(wavelength)
            special_align()
            beamline.generate(wavelength)
            beamline.radiate()
            if move_screen:
                find_focus(beamline, screen, wavelength, dimension=dimension, nrays=nrays, method="Nelder-Mead",
                           show_progress=False, tol=1e-3)
            contributions.append(screen.get_diagram()[dimension].std())
            # contributions.append(correlation_on_screen(screen, nrays, dimension=dimension))
        ret = np.linalg.norm(contributions, norm)
        if show_progress:
            for oe, attrib in zip(oes, attributes):
                print(f"{oe.name}.{attrib} = {oe.__getattribute__(attrib)}")
            print(f"-> FOM = {ret}")
        return ret

    attributes_values_0 = []
    for oe, attrib in zip(oes, attributes):
        attributes_values_0.append(oe.__getattribute__(attrib))
    solution = minimize(FOM, np.array(attributes_values_0), method=method, tol=tol,
                        options=options)
    if verbose:
        print(f"Minimization success: {solution.success}, converged to :")
        for oe, attrib in zip(oes, attributes):
            print(f"{oe.name}.{attrib} = {oe.__getattribute__(attrib)}")
        print(f"-> FOM = {solution.x}")
    beamline.active_chain[0].nrays = int(old_nrays)
    screen.next = old_link
    if solution.success:
        return solution.x[0]
    else:
        raise RuntimeError("Unable to reach an optimum")
