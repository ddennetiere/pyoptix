from bokeh.plotting import figure, show
from bokeh.io import export_png
import numpy as np
from scipy.stats import gaussian_kde
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer, Label
from bokeh.layouts import row, column
from bokeh.util.hex import hexbin
from bokeh.transform import linear_cmap
from bokeh.models import PolyAnnotation, ColumnDataSource, LabelSet
import ipysheet as ipys  # sheet, cell, s_row, s_column, cell_range
import ipywidgets
from IPython.display import display


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
        title += f" on {oe_name}"
    if "beamline_name" in kwargs and "chain_name" in kwargs:
        beamline_name = kwargs.pop("beamline_name")
        chain_name = kwargs.pop("chain_name")
        if beamline_name is not None:
            title += f" of {beamline_name}"
        if chain_name is not None:
            title += f" in config. {chain_name}"
    title += f" at E = {1239.842e-9/columndatasource.data['Lambda'].mean():.1f} eV"
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
    hhist, hedges = np.histogram(x, bins=max(20, int(x.shape[0] / 100)))
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

    x_dist = x - x.mean()
    x_dist = np.sort(x_dist)
    x_integral_dist_pos = np.cumsum(np.where(x_dist > 0, 1, 0))
    x_fwhm_pos = x_dist[x_integral_dist_pos > 0.87 * x_integral_dist_pos.max()][0]
    x_integral_dist_neg = np.cumsum(np.where(x_dist[::-1] < 0, 1, 0))
    x_fwhm_neg = x_dist[::-1][x_integral_dist_neg > 0.87 * x_integral_dist_neg.max()][0]

    # mytext = Label(x=x.mean(), y=-hhist.max() / 2, text='%.2e %s FWHM' % (x_fwhm_pos - x_fwhm_neg, x_unit))
    mytext = Label(x=3 * x.min() / 4, y=-hhist.max() / 2,
                   text=f'{x_fwhm_pos - x_fwhm_neg:.2e} {x_unit} FWHM, {x.std():.2e} {x_unit} RMS')
    ph.add_layout(mytext)

    # create the vertical histogram
    vhist, vedges = np.histogram(y, bins=max(20, int(y.shape[0] / 100)))
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

    y_dist = y - y.mean()
    y_dist = np.sort(y_dist)
    y_integral_dist_pos = np.cumsum(np.where(y_dist > 0, 1, 0))
    y_fwhm_pos = y_dist[y_integral_dist_pos > 0.87 * y_integral_dist_pos.max()][0]
    y_integral_dist_neg = np.cumsum(np.where(y_dist[::-1] < 0, 1, 0))
    y_fwhm_neg = y_dist[::-1][y_integral_dist_neg > 0.87 * y_integral_dist_neg.max()][0]

    mytext = Label(x=-vhist.max() / 2, y=3 * y.max() / 4,
                   text=f'{y_fwhm_pos - y_fwhm_neg:.2e} {y_unit} FWHM, {y.std():.2e} {y_unit} RMS',
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
