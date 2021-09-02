from bokeh.plotting import figure, show
from bokeh.io import export_png
import numpy as np
from scipy.stats import gaussian_kde
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer, Label
from bokeh.layouts import row, column
from bokeh.util.hex import hexbin
from bokeh.transform import linear_cmap
from bokeh.models import PolyAnnotation, ColumnDataSource, LabelSet


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


def plot_spd(dataframe, x_key="x", y_key="y", oe_name="", **kwargs):
    """
    Wrapper function for the PyOptix object of the scatter_plot_2d function. Data is expected as first parameter
    in a pandas.Dataframe, the column names must contain X_key and y_key.

    :param dataframe: Data to be plotted.
    :type dataframe: pandas.Dataframe
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
    if "beamline_name" in kwargs.keys() and "chain_name" in kwargs.keys():
        beamline_name = kwargs.pop("beamline_name")
        chain_name = kwargs.pop("chain_name")
        if beamline_name is not None:
            title += f" of {beamline_name}"
        if chain_name is not None:
            title += f" in config. {chain_name}"
    layout = scatter_plot_2d(dataframe[x_key], dataframe[y_key], title=title,
                             x_label=x_key, x_unit=x_unit, y_label=y_key, y_unit=y_unit, **kwargs)
    return layout


def scatter_plot_2d(x, y, title="", x_unit="", y_unit="", show_map=False, light_plot=False, orthonorm=False, radius=5,
                    x_label="", y_label="", save_in_file="", return_fwhm=False):
    """
    Function that draws a scatter plot using bokeh. Scatter points are colored as a function of local density so as
    to better grip spot shape when density gets high.

    :param x: X data coordinate
    :type x: numpy.ndarray
    :param y: Y data coordinate
    :type y: numpy.ndarray
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
    :return: layout of the plot to be used as parameter of bokeh.plotting.show or tuple (layout, FWHMs)
    :rtype: bokeh.models.layouts.LayoutDOM or (LayoutDOM, (float, float))
    """
    if not light_plot and not show_map:
        # Calculate the point density
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        z = z / z.max()
        colors_jet = ["#%02x%02x%02x" % (to_jet(grayscale)) for grayscale in z]
    else:
        colors_jet = ["blue"]*x.shape[0]
    # create the scatter plot

    if orthonorm:
        x_ptp = np.ptp(x)
        y_ptp = np.ptp(y)
        x_range = (x.mean() - max(x_ptp, y_ptp) / 2, x.mean() + max(x_ptp, y_ptp) / 2)
        y_range = (y.mean() - max(x_ptp, y_ptp) / 2, y.mean() + max(x_ptp, y_ptp) / 2)
    else:
        x_range = (x.min(), x.max())
        y_range = (y.min(), y.max())

    p = figure(tools=TOOLS, plot_width=600, plot_height=600, min_border=10, min_border_left=50,
               toolbar_location="above", x_axis_location=None, y_axis_location=None,
               title=title, x_range=x_range, y_range=y_range)
    p.background_fill_color = "#fafafa"
    p.select(BoxSelectTool).select_every_mousemove = False
    p.select(LassoSelectTool).select_every_mousemove = False

    if show_map:
        bins = hexbin(x, y, min(np.ptp(x), np.ptp(y))/100)

        p = figure(tools="wheel_zoom,reset", match_aspect=True, background_fill_color='#440154')
        p.grid.visible = False

        p.hex_tile(q="q", r="r", size=min(np.ptp(x), np.ptp(y))/100, line_color=None, source=bins,
                   fill_color=linear_cmap('counts', 'Viridis256', 0, max(bins.counts)), alpha=1)
    else:
        p.circle(x, y, size=radius, color=colors_jet, alpha=0.1)

    # create the horizontal histogram
    hhist, hedges = np.histogram(x, bins=max(20, int(x.shape[0] / 1000)))
    hzeros = np.zeros(len(hedges) - 1)
    hmax = max(hhist) * 1.1

    line_args = dict(color="#3A5785", line_color=None)

    ph = figure(toolbar_location=None, plot_width=p.plot_width, plot_height=200, x_range=p.x_range,
                y_range=(-hmax, hmax), min_border=10, min_border_left=50, y_axis_location="right")
    ph.xgrid.grid_line_color = None
    ph.yaxis.major_label_orientation = np.pi / 4
    ph.background_fill_color = "#fafafa"

    ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hhist, color="white", line_color="#3A5785")
    ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.5, **line_args)
    ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.1, **line_args)

    mytext = Label(x=x.mean(), y=-hhist.max() / 2, text='%.2e %s FWHM' % (2.35 * x.std(), x_unit))
    ph.add_layout(mytext)

    # create the vertical histogram
    vhist, vedges = np.histogram(y, bins=max(20, int(y.shape[0] / 1000)))
    vzeros = np.zeros(len(vedges) - 1)
    vmax = max(vhist) * 1.1

    pv = figure(toolbar_location=None, plot_width=200, plot_height=p.plot_height, x_range=(-vmax, vmax),
                y_range=p.y_range, min_border=10, y_axis_location="right")
    pv.ygrid.grid_line_color = None
    pv.xaxis.major_label_orientation = np.pi / 4
    pv.background_fill_color = "#fafafa"

    pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vhist, color="white", line_color="#3A5785")
    pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vzeros, alpha=0.5, **line_args)
    pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vzeros, alpha=0.1, **line_args)

    mytext = Label(x=-vhist.max() / 2, y=y.mean(), text='%.2e %s FWHM' % (2.35 * y.std(), y_unit), angle=-np.pi / 2)
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
