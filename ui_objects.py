import re
import pandas as pd
import pyoptix.classes
from bokeh.plotting import figure, show
from bokeh.io import export_png
import numpy as np
from numpy.polynomial.legendre import legval2d
from numpy.polynomial.polynomial import polyval2d
from scipy.stats import gaussian_kde
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer, Label
from bokeh.layouts import row, column
from bokeh.util.hex import hexbin
from bokeh.transform import linear_cmap
from bokeh.models import PolyAnnotation, ColumnDataSource, LabelSet
import ipysheet as ipys  # sheet, cell, s_row, s_column, cell_range
import ipywidgets
from IPython.display import display, clear_output
import plotly.express as px
import pandas as pd
import plotly.graph_objs as go
from scipy.constants import h, c, eV, degree
import inspect
from numpy.polynomial import Polynomial


# Definition des fonctions d'affichage
def blue(grayscale):
    if grayscale < 0.33:
        return 1
    elif grayscale < 0.5:
        return (0.5 - grayscale) / (0.5 - 0.33)
    else:
        return 0


def green(grayscale):
    if grayscale < 0.33:
        return grayscale / 0.33
    elif grayscale < 0.66:
        return 1
    else:
        return (1 - grayscale) / (1 - 0.66)


def red(grayscale):
    if grayscale < 0.5:
        return 0
    elif grayscale < 0.66:
        return (grayscale - 0.5) / (0.66 - 0.5)
    else:
        return 1


def to_jet(grayscale):
    return (int(255 * red(grayscale)),
            int(255 * green(grayscale)),
            int(255 * blue(grayscale)))


TOOLS = "pan,wheel_zoom,box_select,lasso_select,reset,save, box_zoom"


def format_prefix(value, precision=2, unit="m"):
    """
   Formats a numerical value with an appropriate SI prefix.

   Args:
   - value: A float representing the numerical value to be formatted.
   - precision: An integer representing the number of decimal places to include in the formatted value (default: 2).
   - unit: A string representing the unit of measurement to be appended to the formatted value (default: "m").

   Returns:
   - A string representing the formatted value with the appropriate SI prefix and unit.

   Raises:
   - ValueError: If the exponent is out of range of the available prefixes.
    """
    SI_PREFIX_UNITS = u"yzafpnµm kMGTPEZY"
    negative = False
    digits = precision + 2

    if value == 0.:
        expof10 = 0
    else:
        if value < 0.:
            value = -value
            negative = True
        expof10 = int(np.log10(value))
        if expof10 > 0:
            expof10 = (expof10 // 3) * 3
        else:
            expof10 = (-expof10 + 3) // 3 * (-3)

        value *= 10 ** (-expof10)

        if value >= 1000.:
            value /= 1000.0
            expof10 += 3
        elif value >= 100.0:
            digits -= 2
        elif value >= 10.0:
            digits -= 1

        if negative:
            value *= -1
    expof10 = int(expof10)
    prefix_levels = (len(SI_PREFIX_UNITS) - 1) // 2
    si_level = expof10 // 3

    if abs(si_level) > prefix_levels:
        raise ValueError("Exponent out range of available prefixes.")
    return f"{round(value * 10 ** digits) / 10 ** digits} " + SI_PREFIX_UNITS[si_level + prefix_levels].strip() + unit


def display_progress_bar(max_count):
    """
    Displays and a progressbar jupyterlab widget for long calculations. To update, add 1 to the returned object.value.

    :param max_count: maximum expected value for progressbar
    :type max_count: int
    :return: progress bar widget
    :rtype: ipywidgets.IntProgress
    """
    f = ipywidgets.IntProgress(min=0, max=max_count)  # instantiate the bar
    display(f)  # display the bar
    return f


def display_parameter_sheet(oe):
    """
    Displays the parameters of an optical element in a table jupyterlab widget

    :param oe: optical element the parameters of which to display
    :type oe: pyoptix.OpticalElement or inherited
    :return: table widget
    :rtype: ipysheet.sheet
    """
    properties = oe.dump_properties(verbose=0)
    params = properties["oe_params"]
    params_sheet = ipys.sheet(rows=len(params.keys()), columns=8, column_headers=["parameter", "value", "min", "max",
                                                                                  "multiplier", "type", "group",
                                                                                  "flags"])
    for i, param in enumerate(params):
        ipys.row(i, [param, params[param]["value"], params[param]["bounds"][0], params[param]["bounds"][1],
                     params[param]["multiplier"], params[param]["type"], params[param]["group"],
                     params[param]["flags"]])
    params_sheet.column_width = [15, 20, 10, 10, 10, 10, 10, 10]
    params_sheet.layout = ipywidgets.Layout(width='750px', height='100%')
    for k, c in enumerate(params_sheet.cells):
        c.style['textAlign'] = 'center'
        c.send_state()

    display(params_sheet)
    return params_sheet


def general_FWHM(x):
    # x_dist = x - x.mean()
    # x_dist = np.sort(x_dist)
    # x_integral_dist_pos = np.cumsum(np.where(x_dist > 0, 1, 0))
    # x_fwhm_pos = x_dist[x_integral_dist_pos > 0.87 * x_integral_dist_pos.max()][0]
    # x_integral_dist_neg = np.cumsum(np.where(x_dist[::-1] < 0, 1, 0))
    # x_fwhm_neg = x_dist[::-1][x_integral_dist_neg > 0.87 * x_integral_dist_neg.max()][0]
    # return x_fwhm_pos - x_fwhm_neg
    hist, bin_edges = np.histogram(x, bins="auto", density=True)
    top_half = bin_edges[:-1][hist > hist.max() / 2]
    return np.abs(top_half[-1] - top_half[0])


def plot_spd(columndatasource, x_key="x", y_key="y", oe_name="", **kwargs):
    """
    Wrapper function for the PyOptix object of the scatter_plot_2d function. Data is expected as first parameter
    in a pandas.Dataframe, the column names must contain X_key and y_key.

    :param columndatasource: ColumnDataSource to be plotted with keys given by xkey and ykey
    :type columndatasource: bokeh.models.ColumnDataSource
    :param x_key: name of the Dataframe column containing the X value
    :type x_key: str
    :param y_key:  name of the Dataframe column containing the Y value
    :type y_key: str
    :param kwargs: See scatter_plot_2d doc for additional parameters
    :param oe_name: Name of the optical element where the diagram is captured
    :type oe_name: str
    :return: layout of the plot to be used as parameter of bokeh.plotting.show
    :rtype: bokeh.models.layouts.LayoutDOM
    """
    if x_key[0] == "d":
        x_unit = "rad"
    else:
        x_unit = "m"
    if y_key[0] == "d":
        y_unit = "rad"
    else:
        y_unit = "m"
    title = f"{y_key} vs {x_key}"
    if oe_name != "":
        title += f" {oe_name}"
    beamline_name = kwargs.get("beamline_name")
    chain_name = kwargs.get("chain_name")
    if beamline_name is not None:
        title += f" of {beamline_name}"
    if chain_name is not None:
        title += f" in config. {chain_name}"
    title += f" at E = {1239.842e-9 / columndatasource.data['Lambda'].mean():.1f} eV"
    layout = scatter_plot_2d(columndatasource, x_key, y_key, title=title,
                             x_label=x_key, x_unit=x_unit, y_label=y_key, y_unit=y_unit, **kwargs)
    return layout


def scatter_plot_2d(cds, xkey, ykey, title="", x_unit="", y_unit="", show_map=False, light_plot=False, orthonorm=False,
                    radius=5, x_label="", y_label="", save_in_file="", return_fwhm=False, color_scale="density"):
    """
    Function that draws a scatter plot using bokeh. Scatter points are colored as a function of local density so as
    to better grip spot shape when density gets high.

    :param cds: ColumnDataSource to be plotted with keys given by xkey and ykey
    :type cds: bokeh.models.ColumnDataSource
    :param xkey: X data coordinate
    :type xkey: numpy.ndarray
    :param ykey: Y data coordinate
    :type ykey: numpy.ndarray
    :param title: Title of the plot
    :type title: str
    :param x_unit: Unit of the X coordinate
    :type x_unit: str
    :param y_unit: Unit of the Y coordinate
    :type y_unit: str
    :param show_map: Plot an hexagonal map instead of scatter plot
    :type show_map: bool
    :param light_plot: Does not color scatter points for quick plot
    :type light_plot: bool
    :param orthonorm: Sets X range = Y range so as to have a real life spot shape
    :type orthonorm: bool
    :param radius: Radius of scatter point, default = 5
    :type radius: float
    :param x_label: Label of the X coordinate
    :type x_label: str
    :param y_label: Label of the y coordinate
    :type y_label: str
    :param save_in_file: Name of the file to save the layout to as a png. Default = "", which means not saving to file
    :type save_in_file: str
    :param return_fwhm: If True, returns the Full Width at Half Maximum of the spot is both dimensions as well as the
        plot layout
    :type return_fwhm: bool
    :param color_scale: Value to color the points with, default = "density", acceptable values : any spot diagram
        column name or "density"
    :type color_scale: str
    :return: layout of the plot to be used as parameter of bokeh.plotting.show or tuple (layout, FWHMs)
    :rtype: bokeh.models.layouts.LayoutDOM or (LayoutDOM, (float, float))
    """
    x = np.copy(cds.data[xkey])
    y = np.copy(cds.data[ykey])
    if not light_plot and not show_map:
        if color_scale == "density":
            # Calculate the point density
            xy = np.vstack([x, y])
            z = gaussian_kde(xy)(xy)
            z = z / z.max()
        elif color_scale in cds.data:
            z = (cds.data[color_scale] - cds.data[color_scale].min()) / cds.data[color_scale].ptp()
        else:
            raise AttributeError("Unknown color scale")
        colors_jet = ["#%02x%02x%02x" % (to_jet(grayscale)) for grayscale in z]
    else:
        colors_jet = ["blue"] * x.shape[0]
    cds.data[f"color{xkey}{ykey}"] = colors_jet
    # create the scatter plot

    if orthonorm:
        x_ptp = np.ptp(x)
        y_ptp = np.ptp(y)
        x_range = (x.mean() - max(x_ptp, y_ptp) / 2, x.mean() + max(x_ptp, y_ptp) / 2)
        y_range = (y.mean() - max(x_ptp, y_ptp) / 2, y.mean() + max(x_ptp, y_ptp) / 2)
    else:
        x_range = (x.min(), x.max())
        y_range = (y.min(), y.max())

    p = figure(tools=TOOLS, width=600, height=600, min_border=10, min_border_left=50,
               toolbar_location="left", x_axis_location=None, y_axis_location=None,
               title=title, x_range=x_range, y_range=y_range)
    p.background_fill_color = "#fafafa"
    p.select(BoxSelectTool).select_every_mousemove = False
    p.select(LassoSelectTool).select_every_mousemove = False

    if show_map:
        bins = hexbin(x, y, min(np.ptp(x), np.ptp(y)) / 100)

        p = figure(tools="wheel_zoom,reset", match_aspect=True, background_fill_color='#440154')
        p.grid.visible = False

        p.hex_tile(q="q", r="r", size=min(np.ptp(x), np.ptp(y)) / 100, line_color=None, source=bins,
                   fill_color=linear_cmap('counts', 'Viridis256', 0, max(bins.counts)), alpha=1)
    else:
        p.circle(x=xkey, y=ykey, source=cds, size=radius, color=f"color{xkey}{ykey}", alpha=0.1)
        # p.circle(x, y, size=radius, color=colors_jet, alpha=0.1)

    # create the horizontal histogram
    # hhist, hedges = np.histogram(x, bins=max(20, int(x.shape[0] / 100)))
    hhist, hedges = np.histogram(x, bins="doane")
    hzeros = np.zeros(len(hedges) - 1)
    hmax = max(hhist) * 1.1

    line_args = dict(color="#3A5785", line_color=None)

    ph = figure(toolbar_location=None, width=p.width, height=200, x_range=p.x_range,
                y_range=(-hmax, hmax), min_border=10, min_border_left=50, y_axis_location="right")
    ph.xgrid.grid_line_color = None
    ph.yaxis.major_label_orientation = np.pi / 4
    ph.background_fill_color = "#fafafa"

    ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hhist, color="white", line_color="#3A5785")
    ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.5, **line_args)
    ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.1, **line_args)

    # x_dist = x - x.mean()
    # x_dist = np.sort(x_dist)
    # x_integral_dist_pos = np.cumsum(np.where(x_dist > 0, 1, 0))
    # x_fwhm_pos = x_dist[x_integral_dist_pos > 0.87 * x_integral_dist_pos.max()][0]
    # x_integral_dist_neg = np.cumsum(np.where(x_dist[::-1] < 0, 1, 0))
    # x_fwhm_neg = x_dist[::-1][x_integral_dist_neg > 0.87 * x_integral_dist_neg.max()][0]

    # mytext = Label(x=x.mean(), y=-hhist.max() / 2, text='%.2e %s FWHM' % (x_fwhm_pos - x_fwhm_neg, x_unit))
    mytext = Label(x=x.mean() + 0.75 * (x.min() - x.mean()), y=-hhist.max() / 2,
                   text=f'{format_prefix(general_FWHM(x), precision=2, unit=x_unit)} FWHM, '
                        f'{format_prefix(x.std(), precision=2, unit=x_unit)} RMS')
    # text=f'{general_FWHM(x):.2e} {x_unit} FWHM, {x.std():.2e} {x_unit} RMS')
    ph.add_layout(mytext)

    # create the vertical histogram
    # vhist, vedges = np.histogram(y, bins=max(20, int(y.shape[0] / 100)))
    vhist, vedges = np.histogram(y, bins="doane", density=True)
    vzeros = np.zeros(len(vedges) - 1)
    vmax = max(vhist) * 1.1

    pv = figure(toolbar_location=None, width=200, height=p.height, x_range=(-vmax, vmax),
                y_range=p.y_range, min_border=10, y_axis_location="right")
    pv.ygrid.grid_line_color = None
    pv.xaxis.major_label_orientation = np.pi / 4
    pv.background_fill_color = "#fafafa"

    pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vhist, color="white", line_color="#3A5785")
    pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vzeros, alpha=0.5, **line_args)
    pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vzeros, alpha=0.1, **line_args)

    # y_dist = y - y.mean()
    # y_dist = np.sort(y_dist)
    # y_integral_dist_pos = np.cumsum(np.where(y_dist > 0, 1, 0))
    # y_fwhm_pos = y_dist[y_integral_dist_pos > 0.87 * y_integral_dist_pos.max()][0]
    # y_integral_dist_neg = np.cumsum(np.where(y_dist[::-1] < 0, 1, 0))
    # y_fwhm_neg = y_dist[::-1][y_integral_dist_neg > 0.87 * y_integral_dist_neg.max()][0]

    mytext = Label(x=-vhist.max() / 2, y=y.mean() + 0.75 * (y.max() - y.mean()),
                   text=f'{format_prefix(general_FWHM(y), precision=2, unit=y_unit)} FWHM, '
                        f'{format_prefix(y.std(), precision=2, unit=y_unit)} RMS',
                   # text=f'{general_FWHM(y):.2e} {y_unit} FWHM, {y.std():.2e} {y_unit} RMS',
                   angle=-np.pi / 2)
    pv.add_layout(mytext)

    layout = column(row(p, pv), row(ph, Spacer(width=200, height=200)))

    ph.xaxis.axis_label = x_label
    pv.yaxis.axis_label = y_label

    if save_in_file != "":
        export_png(layout, filename=save_in_file)

    if return_fwhm:
        return layout, (2.35 * x.std(), 2.35 * y.std())
    else:
        return layout


def plot_spd_plotly(df, x_key="x", y_key="y", oe_name="", show_map=False, light_plot=False, orthonorm=False,
                    save_in_file="", return_fwhm=False, **kwargs):
    """
    Wrapper function for the PyOptix object of the scatter_plot_2d function. Data is expected as first parameter
    in a pandas.Dataframe, the column names must contain X_key and y_key.

    :param df: DataFrame to be plotted with keys given by xkey and ykey or OpticalElement the diagram of which to plot
    :type df: pandas.DataFrame or pyoptix.OpticalElement
    :param light_plot: Does not color scatter points for quick plot
    :type light_plot: bool
    :param orthonorm: Sets X range = Y range so as to have a real life spot shape
    :type orthonorm: bool
    :param save_in_file: Name of the file to save the layout to as a png. Default = "", which means not saving to file
    :type save_in_file: str
    :param return_fwhm: If True, returns the Full Width at Half Maximum of the spot is both dimensions as well as the
        plot layout
    :type return_fwhm: bool
    :param show_map: Plot an hexagonal map instead of scatter plot
    :type show_map: bool
    :param x_key: name of the Dataframe column containing the X value
    :type x_key: str
    :param y_key:  name of the Dataframe column containing the Y value
    :type y_key: str
    :param kwargs: See scatter_plot_2d doc for additional parameters
    :param oe_name: Name of the optical element where the diagram is captured
    :type oe_name: str
    :return: layout of the plot to be used as parameter of bokeh.plotting.show
    :rtype: bokeh.models.layouts.LayoutDOM
    """
    if isinstance(df, pyoptix.classes.OpticalElement):
        kwargs["chain_name"] = df.beamline.active_chain_name
        kwargs["beamline_name"] = df.beamline.name
        oe_name = df.name
        df = df.get_diagram()
    elif isinstance(df, ColumnDataSource):
        df = df.data
    else:
        assert isinstance(df, pd.DataFrame)
    if x_key[0] == "d":
        x_unit = "rad"
    else:
        x_unit = "m"
    if y_key[0] == "d":
        y_unit = "rad"
    else:
        y_unit = "m"
    title = f"{y_key} vs {x_key}"
    if oe_name != "":
        title += f" {oe_name}"
    beamline_name = kwargs.get("beamline_name")
    chain_name = kwargs.get("chain_name")
    if beamline_name is not None:
        title += f" of {kwargs.pop('beamline_name')}"
    if chain_name is not None:
        title += f" in config. {kwargs.pop('chain_name')}"
    if "width" not in kwargs:
        kwargs["width"] = 800
    if "height" not in kwargs:
        kwargs["height"] = 800
    title += f" at E = {1239.842e-9 / df['Lambda'].mean():.1f} eV"
    x = np.copy(df[x_key])
    y = np.copy(df[y_key])
    if show_map:
        plot = px.density_heatmap
    else:
        plot = px.scatter
    if not light_plot and not show_map:
        # Calculate the point density
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        z = z / z.max()
        kwargs["color"] = z
    if orthonorm:
        x_ptp = np.ptp(x)
        y_ptp = np.ptp(y)
        x_range = (x.mean() - max(x_ptp, y_ptp) / 2, x.mean() + max(x_ptp, y_ptp) / 2)
        y_range = (y.mean() - max(x_ptp, y_ptp) / 2, y.mean() + max(x_ptp, y_ptp) / 2)
    else:
        x_range = (x.min(), x.max())
        y_range = (y.min(), y.max())
    fig = plot(df, x_key, y_key, marginal_x="histogram", marginal_y="histogram", range_x=x_range, range_y=y_range,
               title=title, **kwargs)
    fig.layout["xaxis"]["title"] = f"{x_key}: {format_prefix(x.std(), precision=2, unit=x_unit)} RMS, " \
                                   f"{format_prefix(general_FWHM(x), precision=2, unit=x_unit)} FWHM"
    fig.layout["yaxis"]["title"] = f"{y_key}: {format_prefix(y.std(), precision=2, unit=y_unit)} RMS, " \
                                   f"{format_prefix(general_FWHM(y), precision=2, unit=y_unit)} FWHM"
    fig.update_coloraxes(showscale=False)
    if save_in_file:
        fig.write_image(save_in_file)
    if return_fwhm:
        return fig, general_FWHM(x), general_FWHM(y)
    else:
        return fig


def plot_polynomial_surface(coeffs, xy_limits, legendre=False, mesh=100, probe_stat=0.002):
    """
    Plots a 3D surface of a polynomial described by its coefficients on a specified XY plane.

    Parameters
    ----------
    coeffs : list or array-like
        Coefficients of the polynomial in descending order of degree.
    xy_limits : list or array-like
        A list or array of length 4 that specifies the limits of the surface on the XY plane.
        The first and second elements correspond to the minimum and maximum values of `x`, respectively,
        while the third and fourth elements correspond to the minimum and maximum values of `y`, respectively.
    probe_stat : float
        A size over which to compute the statistical form and slope errors as with an optical profiler
    mesh : int
        number of sampling points over the surface size for each axis
    legendre : bool
        if True the coefficients are applied to the legendre serie of polynomials, if False over the natural polynomials

    Returns
    -------
    plotly.graph_objects.Figure

    Example
    -------
    >>> coeffs = [1, 2, 3, 4, 5, 6]
    >>> xy_limits = [-5, 5, -10, 10]
    >>> fig = plot_polynomial_surface(coeffs, xy_limits)
    >>> fig.show()

    """
    # Generate a grid of values for x and y
    x = np.linspace(-1, 1, mesh)  # axes are scaled later because the polynomials of legval2d are defined over [-1,1]
    y = np.linspace(-1, 1, mesh)
    X, Y = np.meshgrid(x, y)

    # Calculate z-values from the polynomial coefficients
    if not legendre:
        Z = polyval2d(X, Y, coeffs)
    else:
        Z = legval2d(X, Y, coeffs)

    x_scaled = (x + xy_limits[0]) * (xy_limits[1] - xy_limits[0])
    y_scaled = (y + xy_limits[2]) * (xy_limits[3] - xy_limits[2])
    X, Y = np.meshgrid(x_scaled, y_scaled)
    # Create the 3D surface using Plotly
    fig = go.Figure(data=[go.Surface(x=X, y=Y, z=Z)])

    # Set the figure parameters
    fig.update_layout(title='Polynomial 3D Surface',
                      scene=dict(xaxis_title='X',
                                 yaxis_title='Y',
                                 zaxis_title='Z'))
    if probe_stat:
        Z_y = Z[np.abs(Y) < probe_stat / 2]
        Z_x = Z[np.abs(X) < probe_stat / 2]

        Zp_y = np.diff(Z_y) / np.diff(Y[np.abs(X) < probe_stat / 2])
        Zp_x = np.diff(Z_x) / np.diff(X[np.abs(Y) < probe_stat / 2])

        print(f"Along X : \n\tshape error = {Z_x.std()} m RMS\n\tslope error = {Zp_x.std()} rad RMS")
        print(f"Along Y : \n\tshape error = {Z_y.std()} m RMS\n\tslope error = {Zp_y.std()} rad RMS")
    return fig


def plot_beamline(spots, plot_3D=False, beamline_walls=None, orthonorm=False):
    spots = spots.loc[spots['Intensity'] != 0]
    if "size" not in spots.columns:
        spots.assign(size=0.1)
    if len(list(spots["configuration"].unique())) == 1:
        color_by = "name"
    else:
        color_by = "configuration"
    if plot_3D:
        fig = px.scatter_3d(spots, x="Z", y="X", z="Y", color=color_by,
                            labels={"Z": "S", "X": "X", "Y": "Z"},
                            hover_data=['name', "configuration", "center_s", "center_x"], title="3D View",
                            height=800, size="size"
                            )
        fig.update_layout(scene={'aspectmode': 'data'})
        fig.show()
    else:
        fig = px.scatter(spots, x="Z", y="X", color=color_by,
                         labels={
                             "Z": "S",
                             "X": "X",
                             "Y": "Z"},
                         hover_data=['name', "configuration", "center_s", "center_x"], title="Top view")
        if beamline_walls is not None:
            beamline_walls_fig = px.line({"Z": beamline_walls[:, 1],
                                          "X": beamline_walls[:, 0]}, x="Z", y="X")
            fig.add_trace(beamline_walls_fig.data[0])
        fig.update_layout(scene=dict(xaxis_title='X',
                                     yaxis_title='Z', ), title="Top view")
        if orthonorm:
            fig.update_layout(xaxis=dict(scaleanchor='y', scaleratio=1),
                              yaxis=dict(scaleanchor='x', scaleratio=1))
        fig.show()
        fig = px.scatter(spots, x="Z", y="Y", color=color_by,
                         labels={
                             "Z": "S",
                             "X": "X",
                             "Y": "Z"},
                         hover_data=['name', "configuration", "center_s", "center_z"], title="Side view")
        fig.update_layout(scene=dict(xaxis_title='S',
                                     yaxis_title='Z', ))
        if orthonorm:
            fig.update_layout(xaxis=dict(scaleanchor='y', scaleratio=1),
                              yaxis=dict(scaleanchor='x', scaleratio=1))
        fig.show()


def ellipse(x_center=0, y_center=0, angle=0, x_axis=1, y_axis=1, N=100):
    # x_center, y_center the coordinates of ellipse center
    # ax1 ax2 two orthonormal vectors representing the ellipse axis directions
    # a, b the ellipse parameters
    ax1 = [np.cos(angle), np.sin(angle)]
    ax2 = [-np.sin(angle), np.cos(angle)]
    try:
        assert np.linalg.norm(ax1) == 1
        assert np.linalg.norm(ax2) == 1
    except AssertionError:
        raise ValueError('ax1, ax2 must be unit vectors')
    assert abs(np.dot(ax1, ax2)) < 1e-06, ValueError('ax1, ax2 must be orthogonal vectors')
    t = np.linspace(0, 2 * np.pi, N)
    # ellipse parameterization with respect to a system of axes of directions a1, a2
    xs = x_axis * np.cos(t)
    ys = y_axis * np.sin(t)
    # rotation matrix
    rotation_matrix = np.array([ax1, ax2]).T
    # coordinate of the  ellipse points with respect to the system of axes [1, 0], [0,1] with origin (0,0)
    xp, yp = np.dot(rotation_matrix, [xs, ys])
    x = xp + x_center
    y = yp + y_center
    return x, y


def plot_aperture(stops, title=""):
    fig = go.Figure()
    for stop_number, stop_details in stops.items():
        fill_color = "white"
        if stop_details["opacity"]:
            fill_color = "black"
        if stop_details["kind"] == "Polygon":
            vertex = np.array(stop_details["vertex"])
            fig.add_trace(go.Scatter(x=vertex[:, 0], y=vertex[:, 1], fill="toself", fillcolor=fill_color,
                                     name=f"Polygon #{stop_number}"))
        else:
            ellipse_coords = ellipse(x_center=stop_details["x_center"], y_center=stop_details["y_center"],
                                     angle=stop_details["angle"],
                                     x_axis=stop_details["x_axis"], y_axis=stop_details["y_axis"])
            fig.add_trace(go.Scatter(x=ellipse_coords[0], y=ellipse_coords[1], fill="toself", fillcolor=fill_color,
                                     name=f"Ellipse #{stop_number}"))
    fig.update_layout(title=title)
    if stops["0"]["opacity"] == 0:
        fig.update_layout(plot_bgcolor='black')
    return fig


def plot_beamline_normals(beamline, enlarge_normals=True, orthonormal=True):
    """
    Plots the normal to each optical element surface along the beamline
    :param orthonormal: makes the axes orthonormal (angles appear as in reality, but small offsets are hard to see)
    :type orthonormal: bool
    :param enlarge_normals: multiply the normals by a big factor to make them more visible
    :type enlarge_normals: bool
    :param beamline: Aligned beamline the optical element of which to plot
    :type beamline: pyoptix.Beamline
    :return: None
    :rtype: NoneType
    """
    fig = go.Figure()
    length_bl = beamline.active_chain[-1].get_local_frame()["Center_soleil"][0]
    for oe in beamline.active_chain:
        loc_fr = oe.get_local_frame()
        if not oe.get_transmissive() and enlarge_normals:
            loc_fr["Z_soleil"] *= length_bl / 2
        fig.add_trace(go.Scatter3d(x=[loc_fr["Center_soleil"][0], loc_fr["Center_soleil"][0] + loc_fr["Z_soleil"][0]],
                                   y=[loc_fr["Center_soleil"][1], loc_fr["Center_soleil"][1] + loc_fr["Z_soleil"][1]],
                                   z=[loc_fr["Center_soleil"][2], loc_fr["Center_soleil"][2] + loc_fr["Z_soleil"][2]],
                                   mode="lines+markers",
                                   marker=dict(symbol="diamond", size=2),
                                   name=oe.name
                                   ))
    fig.update_layout(scene_camera_eye=dict(x=-0.76, y=1.8, z=0.92), height=800,
                      title="Optical element normals in " + beamline.active_chain_name,
                      scene=dict(xaxis_title='S', yaxis_title='X', zaxis_title="Z"))
    if orthonormal:
        # fig.update_layout(scene=dict(xaxis_range=[0, length_bl],
        #                              yaxis_range=[-length_bl / 2, length_bl / 2],
        #                              zaxis_range=[-length_bl / 2, length_bl / 2]))
        fig.update_layout(scene={'aspectmode': 'data'}, )
    fig.show()


def generate_widget_for_callable_argument(callable_func):
    """
    Generate JupyterLab widgets to set the argument of a callable function.

    Parameters:
        callable_func (callable): The function for which you want to create widgets.

    Returns:
        widgets.VBox: A VBox widget containing input widgets for the callable's arguments.
    """

    # Get the names of the callable's arguments and their default values
    argspec = inspect.getfullargspec(callable_func)
    arg_names = argspec.args
    defaults = argspec.defaults
    print((arg_names, defaults))

    # Create input widgets for each argument
    input_widgets = []

    for arg_name, default_value in zip(arg_names[2:], defaults):
        # Use text widgets for simplicity, but you can use other widgets as needed
        print(arg_name, default_value)
        if callable(default_value):
            print(default_value)
            lambda_pattern = r'lambda\s+([\w\s]+):[\w\s\.\+\-\/\*]+'
            # Search for lambda functions in the code
            match = re.search(lambda_pattern, inspect.getsource(default_value))
            if match:
                print(match[0])
                default_value = match[0]
            else:
                default_value = inspect.getsource(default_value)
            widget = ipywidgets.Text(
                value=str(default_value),
                description=arg_name,
                disabled=False,
                style={'description_width': 'initial'}
            )

        else:
            widget = ipywidgets.Text(
                value=str(default_value),
                description=arg_name,
                disabled=False,
                style={'description_width': 'initial'}
            )
        input_widgets.append(widget)

    return input_widgets, arg_names[2:]


def config_forge(beamline, lambda_align, lambda_radiate, names=None):
    """
    Creates interactive Jupyterlab widgets that sets up the simulation using alignment parameters and entropy from
    beamline.chains names so that a user friendly selection can be made. Since the meaning of the names of the
    configuration are user dependent, the names can be given in argument. For example, if configuration are names
    in a pattern "<source name>_<grating name>_<endstation name>", names should be
    ["source name","grating","endstation"], otherwise widget will be names "Configuration element 1" and so on.

    Use example :

    .. code-block:: python

        >> display(config_forge(beamline=Hermes, lambda_align=lambda_align, lambda_radiate=lambda_radiate,
            names=["Ondulator","End station","Grating"]))

    :param beamline: beamline to be configured
    :type beamline: pyoptix.classes.Beamline
    :param lambda_align: wavelength at which to align the optical elements in m
    :type lambda_align: float
    :param lambda_radiate: wavelength to radiate from source in m
    :type lambda_radiate: float
    :param names: list of names indexing beamline.chains names, see above.
    :type names: list
    :return: layout containing widgets to display
    :rtype: ipywidgets.VBox
    """
    config_variables = [list(set([config.split("_")[i] for config in beamline.chains.keys()])) for i in
                        range(len(beamline.active_chain_name.split("_")))]
    hc = h * c / eV
    output = ipywidgets.Output()
    widgets_config = []
    invariable = []
    if names is not None:
        names = iter(names)
    else:
        names = iter([f"Configuration element {i}:" for i in range(len(config_variables))])
    for i, variable in enumerate(config_variables):
        if len(variable) > 1:
            widgets_config.append(ipywidgets.Dropdown(options=variable,
                                                      value=beamline.active_chain_name.split("_")[i],
                                                      description=next(names),
                                                      disabled=False,
                                                      style={'description_width': 'initial'})
                                  )
        else:
            invariable.append((i, variable[0]))

    def run_simulation():
        clear_output()
        conf_list = [widget.value for widget in widgets_config]
        for invar in invariable:
            conf_list.insert(*invar)

        lambda_align = hc / E_align.value
        lambda_radiate = hc / E_radiate.value
        beamline.active_chain[0].nrays = N_rays.value
        print(f"Running {N_rays.value} rays of {E_radiate.value:.2f} eV through ", "_".join(conf_list),
              f" aligned for {E_align.value:.2f} eV ")
        try:
            assert "_".join(conf_list) in beamline.chains.keys()
        except AssertionError:
            raise AttributeError("This combination of parameters leads to a non existant configuration")
        beamline.active_chain = "_".join(conf_list)

        align_args = dict()
        for name, widget in zip(align_arg_names, align_widgets):
            try:
                align_args[name] = eval(widget.value)
            except NameError:  # case of string value parameter
                align_args[name] = widget.value
        beamline.align(lambda_align, lambda_radiate,
                       **align_args)
        beamline.clear_impacts(clear_source=True)
        beamline.generate(lambda_radiate)
        beamline.radiate()

    button = ipywidgets.Button(description="Rerun")

    E_align = ipywidgets.FloatText(
        value=hc / lambda_align,
        description='E_align:',
        disabled=False
    )
    E_radiate = ipywidgets.FloatText(
        value=hc / lambda_radiate,
        description='E_radiate:',
        disabled=False
    )
    N_rays = ipywidgets.IntText(
        value=beamline.active_chain[0].nrays,
        description='N rays:',
        disabled=False
    )
    align_widgets, align_arg_names = generate_widget_for_callable_argument(beamline.align_steps)

    def on_button_clicked(b):
        with output:
            run_simulation()

    button.on_click(on_button_clicked)
    show_widget = widgets_config + align_widgets + [E_align, E_radiate, N_rays, button, output]
    return ipywidgets.VBox(show_widget)
