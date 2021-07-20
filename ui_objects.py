from bokeh.plotting import figure, show
from bokeh.io import export_png
import numpy as np
from scipy.constants import *
hc = h*c/eV/nano

from scipy.stats import gaussian_kde
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer, Label
from bokeh.layouts import row, column


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


def plot_spd(dataframe, x_key="x", y_key="y", **kwargs):
    layout = scatter_plot_2d(dataframe[x_key], dataframe[y_key], title=f"{y_key} vs {x_key}",
                             x_label=x_key, x_unit="m", y_label=y_key, y_unit="m", **kwargs)
    return layout


def scatter_plot_2d(x, y, **kwargs):
    if 'title' in kwargs.keys():
        title = kwargs['title']
    else:
        title = ""
    if "x_unit" in kwargs.keys():
        x_unit = kwargs["x_unit"]
    else:
        x_unit = "m"
    if "y_unit" in kwargs.keys():
        y_unit = kwargs["y_unit"]
    else:
        y_unit = "m"
    # Calculate the point density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    z = z / z.max()
    colors_jet = ["#%02x%02x%02x" % (to_jet(grayscale)) for grayscale in z]
    # create the scatter plot

    if 'othonorm' in kwargs.keys():
        if kwargs['orthonorm']:
            x_range = (x.mean() - max(x.ptp(), y.ptp()) / 2, x.mean() + max(x.ptp(), y.ptp()) / 2)
            y_range = (y.mean() - max(x.ptp(), y.ptp()) / 2, y.mean() + max(x.ptp(), y.ptp()) / 2)
        else:
            x_range = (x.min(), x.max())
            y_range = (y.min(), y.max())
    else:
        x_range = (x.min(), x.max())
        y_range = (y.min(), y.max())

    if 'radius' in kwargs.keys():
        radius = kwargs['radius']
    else:
        # radius = min(x.ptp(), y.ptp())/1000
        # radius =  4e-6/z.mean()
        radius = y.std() / 100
        print(radius)

    p = figure(tools=TOOLS, plot_width=600, plot_height=600, min_border=10, min_border_left=50,
               toolbar_location="above", x_axis_location=None, y_axis_location=None,
               title=title, x_range=x_range, y_range=y_range)
    p.background_fill_color = "#fafafa"
    p.select(BoxSelectTool).select_every_mousemove = False
    p.select(LassoSelectTool).select_every_mousemove = False

    r = p.circle(x, y, size=5, color=colors_jet, alpha=0.6)

    # create the horizontal histogram
    hhist, hedges = np.histogram(x, bins=max(20, int(x.shape[0] / 1000)))
    hzeros = np.zeros(len(hedges) - 1)
    hmax = max(hhist) * 1.1

    LINE_ARGS = dict(color="#3A5785", line_color=None)

    ph = figure(toolbar_location=None, plot_width=p.plot_width, plot_height=200, x_range=p.x_range,
                y_range=(-hmax, hmax), min_border=10, min_border_left=50, y_axis_location="right")
    ph.xgrid.grid_line_color = None
    ph.yaxis.major_label_orientation = np.pi / 4
    ph.background_fill_color = "#fafafa"

    ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hhist, color="white", line_color="#3A5785")
    hh1 = ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.5, **LINE_ARGS)
    hh2 = ph.quad(bottom=0, left=hedges[:-1], right=hedges[1:], top=hzeros, alpha=0.1, **LINE_ARGS)

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
    vh1 = pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vzeros, alpha=0.5, **LINE_ARGS)
    vh2 = pv.quad(left=0, bottom=vedges[:-1], top=vedges[1:], right=vzeros, alpha=0.1, **LINE_ARGS)

    mytext = Label(x=-vhist.max() / 2, y=y.mean(), text='%.2e %s FWHM' % (2.35 * y.std(), y_unit), angle=-np.pi / 2)
    pv.add_layout(mytext)

    layout = column(row(p, pv), row(ph, Spacer(width=200, height=200)))

    if "x_label" in kwargs.keys():
        ph.xaxis.axis_label = kwargs["x_label"]
    if "y_label" in kwargs.keys():
        pv.yaxis.axis_label = kwargs["y_label"]

    if "save_in_file" in kwargs.keys():
        export_png(layout, filename=kwargs["save_in_file"])

    if "return_FWHM" in kwargs.keys():
        return layout, (2.35 * x.std(), 2.35 * y.std())
    else:
        return layout