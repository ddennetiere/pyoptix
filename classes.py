# coding: utf-8
import ctypes
import numpy as np
from ctypes import (Union, Structure, c_int, c_double, c_uint32, c_int32, POINTER, c_void_p, cast, c_ulong,
                    create_string_buffer, c_int64, pointer)
from ctypes.wintypes import BYTE, INT, HMODULE, LPCSTR, HANDLE, DOUBLE, UINT
from .mcpl import PyOptixMCPLWriter
from .exposed_functions import (get_parameter, set_parameter, align, generate, radiate, enumerate_parameters,
                                get_element_name, set_recording, get_next_element, get_previous_element,
                                get_spot_diagram, chain_element_by_id, create_element, clear_impacts, get_impacts_data,
                                set_transmissive, get_transmissive, get_hologram_pattern, fit_surface_to_slopes,
                                fit_surface_to_heights, get_array_parameter, get_array_parameter_size,
                                get_array_parameter_dims, dump_parameter, dump_arg_parameter,
                                get_parameter_flags, Bounds, Parameter, ParamArray, FrameID, RecordingMode, Diagram,
                                GratingPatternInformation, save_as_xml, generate_polarized, load_from_xml,
                                enumerate_stops, set_aperture_activity, get_aperture_activity,
                                add_rectangular_stop, get_polygon_parameters, add_polygonal_stop, insert_polygonal_stop,
                                replace_stop_by_polygon, insert_rectangular_stop, replace_stop_by_rectangle,
                                get_ellipse_parameters, add_elliptical_stop, insert_elliptical_stop,
                                replace_stop_by_ellipse, add_circular_stop, insert_circular_stop,
                                replace_stop_by_circle, get_surface_frame, find_element_id, get_exit_frame)
from scipy.constants import degree, milli
from lxml import etree
import pandas as pd
from .ui_objects import show, plot_spd, figure, PolyAnnotation, ColumnDataSource, LabelSet, display_parameter_sheet, \
    plot_beamline, plot_aperture, general_FWHM
from numpy import pi, cos, sin, tan, arccos, arcsin, arctan
from scipy.signal import find_peaks, peak_widths
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
import pickle

# dictionnary for optix to pyoptix attribute import
optix_dictionnary = {
    "DX": "d_x",
    "DY": "d_y",
    "DZ": "d_z",
    "Dphi": "d_phi",
    "Dpsi": "d_psi",
    "Dtheta": "d_theta",
    "distance": "distance_from_previous",
    "sigmaX": "sigma_x",
    "sigmaY": "sigma_y",
    "sigmaXdiv": "sigma_x_div",
    "sigmaYdiv": "sigma_y_div",
    "nRays": "nrays",
    "azimuthAngle1": "azimuth_angle1",
    "azimuthAngle2": "azimuth_angle2",
    "elevationAngle1": "elevation_angle1",
    "inverseDist1": "inverse_distance1",
    "inverseDist2": "inverse_distance2",
    "recordingWavelength": "recording_wavelength",
    "lineDensity": "line_density",
    "invp": "inverse_p",
    "invq": "inverse_q",
    "waistX": "waist_x",
    "waistY": "waist_y",
    "surfaceLimits": "surface_limits",
    "apertureX": "aperture_x",
    "trajectoryRadius": "trajectory_radius",

}


class PostInitMeta(type):
    """
    Metaclass for freezing definition of attributes outside init
    """

    def __call__(cls, *args, **kw):
        instance = super().__call__(*args, **kw)  # < runs __new__ and __init__
        instance.__post_init__()
        return instance


array_parameter_list = ['coefficients', "surfaceLimits"]


class ChainList(list):
    """
    Class inheriting list allowing clear and concise text prints of beamlines
    """

    def __repr__(self):
        ret_str = ''
        for oe in self:
            ret_str += f"{oe.name} -> "
        ret_str += "\n"
        return ret_str[:-4]


class Point(object):
    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"\t{self.x}\n\t{self.y}\n\t{self.z}"

    def get_spherical_coord(self):
        r = np.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
        theta = arccos(self.z / r)
        if self.y >= 0:
            phi = arccos(self.x / np.sqrt(self.x ** 2 + self.y ** 2))
        else:
            phi = 2 * pi - arccos(self.x / np.sqrt(self.x ** 2 + self.y ** 2))
        return r, theta, phi

    def __add__(self, point2):
        return Point(self.x + point2.x, self.y + point2.y, self.z + point2.z)

    def __sub__(self, point2):
        return Point(self.x - point2.x, self.y - point2.y, self.z - point2.z)

    def get_coord_horiba(self):
        r, theta, phi = self.get_spherical_coord()
        return r, theta / degree, phi / degree

    def get_coord_solemio(self, point1=None):
        r, theta, phi = self.get_spherical_coord()
        if point1:
            r1, theta1, phi1 = point1.get_spherical_coord()
            return 1 / r, cos(pi / 2 - theta) - cos(pi / 2 - theta1), phi
        else:
            return 1 / r, cos(pi / 2 - theta), phi

    def get_coord_pyoptix(self):
        r, theta, phi = self.get_spherical_coord()
        return 1 / r, pi / 2 - theta, phi


class SphericalPoint(Point):
    def __init__(self, r, elevation, azimuth):
        super().__init__()
        # from convention elevation is 90-\theta and azimuth is \phi
        self.x = np.longdouble(r * cos(azimuth) * sin(90 * degree - elevation))
        self.y = np.longdouble(r * sin(azimuth) * sin(90 * degree - elevation))
        self.z = np.longdouble(r * cos(90 * degree - elevation))


class BeamlineChainDict(dict):
    """
    Class inheriting from dict for holding and representing all possible optical elements chains in a beamline.
    Allows for convenient commutations between chains.
    """

    def __setitem__(self, key, value):
        assert isinstance(value, list)
        super(BeamlineChainDict, self).__setitem__(key, ChainList(value))

    def __repr__(self):
        ret_str = ""
        for key in self.keys():
            ret_str += f"Chaîne {key}:\n\t"
            chain = self.__getitem__(key)
            ret_str += chain.__repr__() + "\n"
        return ret_str


class Beamline(object):
    """
    Class for describing a beamline in optix
    """

    def __init__(self, name="Beamline"):
        """
        Constructor of the class Beamline which holds the name of the beamline and all configurations of it called
        chains.

        :param name: Name of the beamline
        :type name: str
        """
        super().__init__()
        self._elements = []
        self._chains = BeamlineChainDict({})
        self._active_chain = None
        self.active_chain_name = None
        self.name = name
        self.optical_distances = None
        self.align_steps = lambda lambda_align: None
        self.n_generated_rays = 0

    def save_configuration(self, filename=None):
        """
        Saves all the parameters of a beamline active chain and returns them as a list of dictionaries.
        If filename is supplied, saves the configuration object returned in a file using pickle

        :param filename: name of the file where to save configuration, default: None
        :type filename: str
        :return: list of parameters as dict
        :rtype: list of dict
        """
        config = []
        for oe in self.active_chain:
            config.append(oe.dump_properties(verbose=0))
        if filename is not None:
            with open(filename, "wb") as filout:
                pickle.dump(config, filout)
        return config

    def load_configuration(self, configuration=None, filename=None):
        """
        Loads all the parameters of a beamline active chain as a list of dictionaries.
        If filename is supplied, loads the configuration object from pickle file and ignores configuration param.

        :param configuration: list of dictionnaries containing the parameters of the optical element. See output
            of save_configuration for more details.
        :type configuration: list of dict
        :param filename: name of the file from which to load configuration, default: None
        :type filename: str
        :return: None
        :rtype: Nonetype
        """
        if filename is not None:
            with open(filename, "rb") as filin:
                configuration = pickle.load(filin)
        for oe in self.active_chain:
            for description in configuration:
                if oe.name == description["oe_name"]:
                    params = description["oe_params"]
                    for param in params:
                        oe._set_parameter(param, params[param])

    def save_beamline(self, filename=None):
        if filename:
            save_as_xml(filename)

    def load_beamline(self, filename=None):
        if filename:
            load_from_xml(filename)

    @property
    def chains(self):
        """
        :return: all the beamline chains as BeamlineChainDict
        """
        return self._chains

    @property
    def active_chain(self):
        """
        :return: the chain along which rays will be propagated
        """
        return self._active_chain

    @active_chain.setter
    def active_chain(self, chain_name):
        """
        Sets the chain whose key is `chain_name` as the active chain and links all optical elements accordingly.

        :param chain_name: name of the chain to become active
        :type chain_name: str
        :return: None
        """
        assert chain_name in self._chains
        self.active_chain_name = chain_name
        self.optical_distances = None
        self._active_chain = self.chains[chain_name]
        for i, oe in enumerate(self._active_chain):
            oe.beamline = self
            try:
                oe.next = self._active_chain[i + 1]
            except IndexError:
                pass
            if i:
                oe.previous = self._active_chain[i - 1]
        self._active_chain[-1].next = None
        ret_str = f"Chaîne {chain_name}:\n\t"
        ret_str += self._active_chain.__repr__()
        print(ret_str)

    def align(self, lambda_align, from_element=None, **kwargs):
        """
        Computes the absolute positions of the optics using the optics parameters. To be called before radiate.

        :param lambda_align: Wavelength to be used for coordinate calculations in m. Can be different from actual
            radiated wavelength
        :type lambda_align: float
        :param from_element: Source element from which to compute the coordinates. Default is first element of
            active_chain
        :type from_element: OpticalElement or inherited
        :param kwargs: Keyword arguments passed to beamline.align_steps
        :type kwargs: Keyword arguments
        :return: code result from relative optix function
        """

        self.get_distance_between_oe(None, None)
        try:
            self.align_steps(lambda_align, **kwargs)
        except AttributeError:
            raise AttributeError("Beamline's method align_steps must have one positional argument which is the "
                                 "wavelength of alignment and can have as many keyword argument as needed."
                                 "Default is lambda wavelength: None. ")
        if from_element is not None:
            ret = align(from_element.element_id, lambda_align)
        else:
            ret = align(self.active_chain[0].element_id, lambda_align)
        return ret

    def get_distance_between_oe(self, oe1, oe2):
        """
        Returns the optical distance traveled by the chief ray between oe1 and oe2 if oe1 and oe2 are not None.

        :param oe1: First optical element
        :type oe1: pyoptix.OpticalElement instance or inherited class
        :param oe2: Second optical element
        :type oe2: pyoptix.OpticalElement instance or inherited class
        :return: None if oe1 or oe2 is None, distance between oe1 and oe2 otherwise.
        :rtype: float
        """
        totdist = 0
        self.optical_distances = {}
        for oe in self.active_chain:
            totdist += oe.distance_from_previous
            self.optical_distances[oe.name] = totdist

        if oe1 is None:
            return
        if oe2 is None:
            return
        try:
            assert oe1 in self.active_chain
            assert oe2 in self.active_chain
        except AssertionError:
            raise AssertionError(f"Both {oe1.name} ({oe1 in self.active_chain}) and "
                                 f"{oe2.name} ({oe2 in self.active_chain}) must be in active chain or be None")
        return self.optical_distances[oe1.name] - self.optical_distances[oe2.name]

    def clear_impacts(self, clear_source=False):
        """
        Removes any impact on the spot diagrams computed on all element following the source object at the exception
        of the source itself, as it needs its impact to repropagate rays from the previous distribution.

        :param clear_source: also clears the source
        :type clear_source: bool
        :return: code result from the optix function
        """
        from_element = 1
        if clear_source:
            from_element = 0
            self.n_generated_rays = 0
        return clear_impacts(self.active_chain[from_element].element_id)

    def radiate(self, from_element=None):
        """
        Propagates rays from a source element.

        :param from_element: Source element from which to propagate rays
        :type from_element: Source element or inherited
        :return:
        """
        if from_element is not None:
            return radiate(from_element.element_id)
        else:
            return radiate(self.active_chain[0].element_id)

    def generate(self, lambda_radiate, polarization=None):
        """
        Generates rays at `lambda_radiate` wavelength.

        :param lambda_radiate: Wavelength of the rays in m
        :param polarization: Polarization of the generated rays which can be 'S', 'P', 'R' (circular right),
         'L'(circular left)
        :return:
        """
        if polarization:
            n_generated_rays = generate_polarized(self.active_chain[0].element_id, lambda_radiate, polarization)
        else:
            n_generated_rays = generate(self.active_chain[0].element_id, lambda_radiate)
        self.n_generated_rays += n_generated_rays

    def _add_element(self, new_element):
        if new_element not in self._elements:
            self._elements.append(new_element)

    def show_active_chain_orientation(self):
        """
        Prints the optical element in the active chain with the orientation of their local vertical axis.
        For example :

        - a vertical deflecting mirror with a reflective surface pointed up (deflection towards +Z) will
          appear will an "up" arrow,
        - a screen should appear with an "up" arrow, meaning that its "y" axis is vertical and pointed upwards,
        - an horizontal deflecting mirror with a reflective surface on the left side (deflection towards +X) will
          appear with a "left" arrow.

        :return: None
        :rtype: Nonetype
        """
        phi = 0
        arrows = {"up": "\u2191", "left": "\u2190", "right": "\u2192", "down": "\u2193"}
        for oe in self.active_chain:
            phi += oe.phi
            local_phi = (phi + 2 * pi) % (2 * pi)
            if local_phi < pi / 4:
                arrow = arrows["up"]
            elif local_phi < 3 * pi / 4:
                arrow = arrows["left"]
            elif local_phi < 5 * pi / 4:
                arrow = arrows["down"]
            elif local_phi < 7 * pi / 4:
                arrow = arrows["right"]
            else:
                arrow = arrows["up"]
            print(oe.name, arrow)

    def draw_active_chain(self, top_only=False, side_only=False):
        """
        Draws a *not to scale* diagram of the beamline top and side views for debug purposes.

        :param top_only: If True only draws the beamline viewed from the top
        :type top_only: bool
        :param side_only: If True only draws the beamline viewed from the side
        :type side_only: bool
        :return: None
        """

        def rotate(x, y, theta=0):
            xp = x * np.cos(theta) + y * np.sin(theta)
            yp = -x * np.sin(theta) + y * np.cos(theta)
            return xp, yp

        def make_box(oe_type="film", center=(0, 0), angle=0, fig=None, direction="straight", height=10,
                     color="blue"):
            if oe_type == "film":
                width = 1
            else:
                width = 10
            direction_offset = 0
            if direction == "up":
                direction_offset = -height / 2
            if direction == "down":
                direction_offset = height / 2

            x1, x2 = -width / 2, -width / 2
            x3, x4 = width / 2, width / 2
            y2, y3 = -height / 2 + direction_offset, -height / 2 + direction_offset
            y1, y4 = height / 2 + direction_offset, height / 2 + direction_offset
            x1, y1 = rotate(x1, y1, angle)
            x2, y2 = rotate(x2, y2, angle)
            x3, y3 = rotate(x3, y3, angle)
            x4, y4 = rotate(x4, y4, angle)
            polygon = PolyAnnotation(
                fill_color=color,
                fill_alpha=0.3,
                xs=[x1 + center[0], x2 + center[0], x3 + center[0], x4 + center[0]],
                ys=[y1 + center[1], y2 + center[1], y3 + center[1], y4 + center[1]],
            )
            fig.add_layout(polygon)

        print(self.active_chain)
        total_phi = 0
        total_theta_side = 0
        total_theta_top = 0
        top_points = [(0, 0)]
        side_points = [(0, +50)]
        length = 25
        p = figure()
        for oe in self.active_chain:
            total_phi += oe.phi
            if oe.theta == 0:
                make_box(oe_type="film", center=top_points[-1], angle=-total_theta_top * pi / 180, fig=p,
                         direction="straight", height=10)
                make_box(oe_type="film", center=side_points[-1], angle=-total_theta_side * pi / 180, fig=p,
                         direction="straight", height=10)
            elif abs(abs(total_phi) % pi) < 1e-5:
                if abs(total_phi % (2 * pi) - pi) < 1e-5:
                    make_box(oe_type="mirror", center=side_points[-1], angle=-(total_theta_side - 15) * pi / 180,
                             fig=p, direction="down", height=10)
                    total_theta_side -= 30
                else:
                    make_box(oe_type="mirror", center=side_points[-1], angle=-(total_theta_side + 15) * pi / 180,
                             fig=p, direction="up", height=10)
                    total_theta_side += 30
                make_box(oe_type="mirror", center=top_points[-1], angle=-total_theta_top * pi / 180, fig=p,
                         direction="straight", height=10)
            elif abs(abs(total_phi) % (pi / 2)) < 1e-5:
                if abs(total_phi % (pi * 2) - pi / 2) < 1e-5:
                    make_box(oe_type="mirror", center=top_points[-1], angle=-(total_theta_top + 15) * pi / 180,
                             fig=p, direction="up", height=10)
                    total_theta_top += 30
                else:
                    make_box(oe_type="mirror", center=top_points[-1], angle=-(total_theta_top - 15) * pi / 180,
                             fig=p, direction="down", height=10)
                    total_theta_top -= 30
                make_box(oe_type="mirror", center=side_points[-1], angle=-total_theta_side * pi / 180, fig=p,
                         direction="straight", height=10)
            else:
                raise Exception("unable to parse oe", oe.name)
            top_points.append((top_points[-1][0] + length,
                               top_points[-1][1] + length * np.tan(total_theta_top * pi / 180)))
            side_points.append((side_points[-1][0] + length,
                                side_points[-1][1] + length * np.tan(total_theta_side * pi / 180)))
        side_points = np.array(side_points)
        top_points = np.array(top_points)
        source_top = ColumnDataSource(data=dict(X=top_points[:-1, 0],
                                                Y=top_points[:-1, 1],
                                                names=[element.name for element in self.active_chain]))
        source_side = ColumnDataSource(data=dict(X=side_points[:-1, 0],
                                                 Y=side_points[:-1, 1],
                                                 names=[element.name for element in self.active_chain]))
        if not side_only:
            p.line(top_points[:-1, 0], top_points[:-1, 1], color="blue", legend_label="top view")
            labels_top = LabelSet(x='X', y='Y', text='names', angle=30,
                                  x_offset=0, y_offset=-15, source=source_top)  # , render_mode='canvas')
            p.add_layout(labels_top)
        if not top_only:
            p.line(side_points[:-1, 0], side_points[:-1, 1], color="red", legend_label="side view")
            labels_side = LabelSet(x='X', y='Y', text='names', angle=-30,
                                   x_offset=0, y_offset=15, source=source_side)  # , render_mode='canvas')
            p.add_layout(labels_side)
        show(p)

    def get_resolution(self, mono_slit=None, wavelength=None, orientation="vertical", dlambda_over_lambda=1 / 5000,
                       show_spd=False, verbose=0, nrays=5000, criterion="fwhm", return_all=False, distance_from_slit=0,
                       **kwargs):
        """
        Computes the resolution of a beamline in its `mono_slit` plane at a given `wavelength`. An a priori resolution
        must be given as `dlambda_over_lambda` for calculation purposes and the orientation of deviation relative
        to the slit plane must also be given. THe computation is as follows : two lambdas are  genreated and propagated,
        their distance in the slit plane is computed and their width (either 2.35*rms or fwhm from histogram
        depending on parameter criterion).
        Resolution is then the lambda/dlambda such as the two spots don't overlap.

        :param distance_from_slit: offset along incoming beam at which to measure resolution from mono_slit position
        :type distance_from_slit: float
        :param criterion: Criterion for resolution computation either rms or fwhm.
        :type criterion: str
        :param mono_slit: Slit plane
        :type mono_slit: pyoptix.OpticalElement
        :param wavelength: wavelength at which to compute resolution
        :type wavelength: float
        :param orientation: Orientation of the grating deviation "vertical" or "horizontal"
        :type orientation: str
        :param dlambda_over_lambda: estimation of inverse of resolution
        :type dlambda_over_lambda: float
        :param show_spd: If True displays the calculated spots and the Y projection of spot at wavelength.
            Default: False
        :type show_spd: bool
        :param verbose: If > 0, prints details of the resolution calculations. Default: 0
        :type verbose: int
        :param nrays: Number of rays per spot to be propagated in order to compute the resolution. Default: 5000
        :type nrays: int
        :param return_all: If True returns a dictionary containing spot size according to criterion, dispersion (defined
            as the distance between spots at lambda and at lambda+dlambda in µm), lambda (in m),
            dlambda, if False returns the resolution in lambda/dlambda
        :return: Resolution in lambda/dlambda of dict (see return_all)
        :rtype: float
        """
        self.clear_impacts(clear_source=True)
        stored_nrays = self.active_chain[0].nrays
        slit_next_OE = mono_slit.next
        mono_slit.next = None
        self.active_chain[0].nrays = nrays
        self.align(wavelength, **kwargs)
        self.generate(wavelength)
        self.generate(wavelength + wavelength * dlambda_over_lambda)
        self.radiate()
        if orientation == "vertical":
            dim = "Y"
        elif orientation == "horizontal":
            dim = "X"
        else:
            raise AttributeError("Unknown orientation")
        spd = mono_slit.get_diagram(distance_from_slit)
        spd = spd[spd["Intensity"] != 0]
        if show_spd:
            print(spd)
            mono_slit.show_diagram(distance_from_slit)
        projection = np.array(spd.where(spd["Lambda"] == wavelength).dropna()[dim])
        projection_dl = np.array(
            spd.where(spd["Lambda"] == (wavelength + wavelength * dlambda_over_lambda)).dropna()[dim])
        if criterion == "fwhm":
            # vhist, vedges = np.histogram(spd.where(spd["Lambda"] == wavelength).dropna()[dim], bins=100)
            # peaks, _ = find_peaks(vhist, height=vhist.max())
            # res_half = peak_widths(vhist, peaks, rel_height=0.5)
            # mono_chr_fwhm = (res_half[0] * (vedges[1] - vedges[0]))[0]
            mono_chr_fwhm = general_FWHM(spd.where(spd["Lambda"] == wavelength).dropna()[dim])
        else:
            mono_chr_fwhm = np.array(spd.where(spd["Lambda"] == wavelength).dropna()[dim]).std()
        if show_spd and criterion == "fwhm":
            import matplotlib.pyplot as plt
            vhist, vedges = np.histogram(spd.where(spd["Lambda"] == wavelength).dropna()[dim], bins="auto")
            plt.plot(vedges[:-1], vhist)
            # plt.plot(vedges[1:], vhist)
            # plt.plot(vedges[peaks], vhist[peaks], "x")
            # h, left, right = res_half[1:]
            # widths = (h, vedges[int(left[0])], vedges[int(right[0])])
            # plt.hlines(*widths, color="C2")
            plt.show()
        distance = abs(np.mean(projection) - np.mean(projection_dl))
        resolution = (1 / dlambda_over_lambda) * distance / mono_chr_fwhm
        if verbose:
            if verbose >1:
                print(f"FWHM monochromatique : {mono_chr_fwhm * 1e6:.2f} µm")
                print(f"dispersion dans le plan (m) : {distance * 1e6:.2f} µm")
                print(f"Lambda = {wavelength * 1e9:.5f} nm")
                print("lambda_over_dlambda for calculation :", 1 / dlambda_over_lambda)
                print("calculated resolution :", resolution)
            else:
                print(f"FWHM @ {1239.842e-9/wavelength:.2f} eV : {mono_chr_fwhm * 1e6:.2f} µm, dispersion :"
                      f"{distance * 1e6:.2f} µm, resolution : {resolution}")
        self.active_chain[0].nrays = stored_nrays
        mono_slit.next = slit_next_OE

        if return_all:
            return {criterion: mono_chr_fwhm * 1e6, "dispersion_on_screen": distance * 1e6, "wavelength": wavelength,
                    "dlambda": wavelength * dlambda_over_lambda}
        else:
            return resolution

    def draw_to_scale(self, wavelength=300e-9, radiate=False, configurations=[], plot3D=False, beamline_walls=None,
                      **kwargs):
        """
        Generate and plot the beamline spot diagram to scale for each recording optical element
        for the specified wavelength(s) and configurations.


        Parameters
        ----------
        wavelength : float or list of floats, optional
            Wavelength(s) for which to generate the beamline diagram (default is 300e-9). If list, length must match
            the number of configuration.
        radiate : bool, optional
            If True, radiate the beamline for each configuration (default is False).
            Set to True if multiple configurations are provided
        configurations : list of str, optional
            Configurations to be plotted (default is []). If [], the current configuration of the beamline is used
        plot3D : bool
            If True, the spots are plotted in 3 dimensions. Beware if the beamline is long,
            the plot is hard to navigate
        beamline_walls : numpy.array
            Array of (X,S) pair of coordinates representing the physical beamline corners that gets plotted in the top
            view diagram
        kwargs : keyword arguments, optional
            Keyword arguments passed to beamline align_steps method

        Returns
        -------
        pandas DataFrame
            DataFrame containing the recorded impact data for the plotted configurations

        Notes
        -----
        - If configurations is provided, the beamline will be modified and radiated, regardless of the radiate flag.
        - If radiate is True or multiple configurations are specified, a warning will be printed.
        - If wavelength is a list, it must have the same length as configurations.
        - The beamline diagram is generated to scale for the specified wavelength and configurations.
        - The recorded impact data for each configuration is collected and returned as a DataFrame.
        - The X-coordinate values in the impact data are flipped (multiplied by -1) to match the plotting convention.
        - The beamline diagram is plotted using the plot_beamline function.
        - The PyOptix reference frame (X,Y,Z) translates in the SOLEIL frame as (-X, Z, S)
        """

        if configurations:
            radiate = True
        else:
            configurations = [self.active_chain_name]
        if radiate:
            print("WARNING : for radiate = True or multiple configurations, beamline will be modified after plot")
        if isinstance(wavelength, list):
            assert len(wavelength) == len(configurations)
        else:
            wavelength = [wavelength] * len(configurations)
        diags = []
        for config, lamda in zip(configurations, wavelength):
            if radiate:
                self.active_chain = config
                for oe in self.active_chain:
                    oe.recording_mode = RecordingMode.recording_output
                self.clear_impacts(clear_source=True)
                self.align(lamda, **kwargs)
                self.generate(lamda)
                self.radiate()
            for oe in self.active_chain:
                if oe.recording_mode != RecordingMode.recording_none:
                    impacts = oe.get_impacts(reference_frame="general_frame")
                    center = oe.get_local_frame()["Center_pyoptix"]
                    impacts["name"] = oe.name
                    impacts["center_x"] = -center[0]
                    impacts["center_z"] = center[1]
                    impacts["center_s"] = center[2]
                    impacts["configuration"] = config
                    impacts["X"] *= -1
                    diags.append(impacts)
        spots = pd.concat(diags)
        plot_beamline(spots, plot_3D=plot3D, beamline_walls=beamline_walls)
        return spots


class OpticalElement(metaclass=PostInitMeta):
    """
    Base class for all pyoptix optical element
    """
    _frozen = False

    def __init__(self, name="", phi=0, psi=0, theta=0, d_phi=0, d_psi=0, d_theta=0,
                 d_x=0, d_y=0, d_z=0, next_element=None, previous=None, distance_from_previous=0, element_id=None,
                 element_type="", beamline=None):
        """
        Constructor for all optical elements. Parameters can be set in constructor or at a later time.

        Underlying OpticalElement parameters are complex c structure for optimization ond different magic handling.
        See class Parameter docstring.
        Two methods of setting parameters are provided: implicit and explicit. Implicit method
        `OpticalElement.parameter = value` only sets the value of the parameter without changing any other field.
        Explicit method `OpticalElement.parameter = {"key":value,...}` sets any parameter field without changing those
        not in the provided dictionary keys.

        Positioning referential of OpticalElement is linked to its surface : Z is normal to the surface,
        X is along it width, Z along its height. For optical elements in normal incidence such as films usually, this
        convention makes Y vertical in absence of phi modifier. Warning : referentials are inherited from an element
        to the next. A horizontally deflecting mirror should have phi=pi/2, but the next mirror if deflecting
        horizontally towards the other side should have phi=pi. See provided examples and use method
        Beamline.draw_active_chain profusely.

        If parameter `element_id` is set to an preexisting pyoptix object, all parameters will be imported from
        pyoptix

        :param name: Name of the OpticalElement as will be stored in optix and displayed by Beamline
        :type name: str (32 char long max)
        :param phi: Angle of rotation around Y (rad) (roll)
        :type phi: float
        :param psi: Angle of rotation around Z (rad) (yaw)
        :type psi: float
        :param theta: Angle of rotation around X (rad) (incidence or pitch)
        :type theta: float
        :param d_phi: Misalignment of angle of rotation around Y (rad) (roll), not transferred to the next element
        :type d_phi: float
        :param d_psi: Misalignment of angle of rotation around Z (rad) (yaw), not transferred to the next element
        :type d_psi: float
        :param d_theta: Misalignment of angle of rotation around X (rad) (pitch), not transferred to the next element
        :type d_theta: float
        :param d_x: Misalignment of the position along X (not transferred to the next element) in m
        :type d_x: float
        :param d_y: Misalignment of the position along Y (not transferred to the next element) in m
        :type d_y: float
        :param d_z: Misalignment of the position along Z (not transferred to the next element) in m
        :type d_z: float
        :param next_element: Next element in the active chain
        :type next_element: OpticalElement or inherited
        :param previous: Previous element in the chain
        :type previous: OpticalElement or inherited
        :param distance_from_previous: distance from the previous optical element in m
        :type distance_from_previous: float
        :param element_id: Handle of underlying optix object if object already exists
        :type element_id: wintypes.INT
        :param element_type: class of the OpticalElement (if unsure of which to chose, use an inherited class)
        :type element_type: str
        :param beamline: Beamline of the element
        :type beamline: pyoptix.Beamline
        """
        super().__init__()
        self._recording_mode = RecordingMode.recording_none
        self._element_id = find_element_id(name)
        self._element_type = element_type
        if element_id is not None:
            self.from_element_id(element_id)
            self._element_id = element_id
        elif self._element_id is not None:  # the element already exists
            pass
        elif element_type != "" and name != "":
            self._element_id = create_element(element_type, name)
            self._name = name
        else:
            raise AttributeError("please provide either element_type and name or element_id")
        self._name = name
        self.phi = phi
        self.psi = psi
        self.theta = theta
        self.d_phi = d_phi
        self.d_psi = d_psi
        self.d_theta = d_theta
        self.d_x = d_x
        self.d_y = d_y
        self.d_z = d_z
        self.next = next_element
        self.previous = previous
        self.distance_from_previous = distance_from_previous
        self.beamline = beamline

    def __setattr__(self, name, value):  # enable the frozen decorator for child classes
        if not self._frozen or hasattr(self, name):
            super().__setattr__(name, value)
        else:
            raise AttributeError(f"Forbidden attempt to define new attribute {name} of a frozen object")

    def __post_init__(self):
        self._frozen = True

    def get_whole_parameter(self, param_name):
        param = Parameter()
        inv_dict_optix = {v: k for k, v in optix_dictionnary.items()}
        if param_name in inv_dict_optix:
            param_name = inv_dict_optix[param_name]
        get_parameter(self._element_id, param_name, param)
        # return {"value": param.value, "bounds": (param.bounds.min, param.bounds.max),
        #         "multiplier": param.multiplier, "type": param.type, "group": param.group, "flags": param.flags}
        return {"value": param.value, "bounds": (param.bounds.min, param.bounds.max),
                "multiplier": param.multiplier, "type": param.type, "group": param.group, "flags": param.flags}

    def _get_parameter(self, param_name):
        param = self._get_c_parameter(param_name)
        if param.flags & 0x8:
            return param.array
        else:
            return param.value

    def _get_c_parameter(self, param_name):
        param = Parameter()
        flags = c_uint32()
        get_parameter_flags(self.element_id, param_name, flags)
        if flags.value & 0x8:
            # the paracmeter contains an array
            dims = (c_int64 * 2)()
            if get_array_parameter_dims(self.element_id, param_name, dims) != 0:
                param.array = np.empty([dims[1], dims[0]])
                get_array_parameter(self.element_id, param_name, param, dims[0] * dims[1])
        else:
            # release the array data and unset flags array bit, if needed
            param.flags &= ~(1 << 3)
            get_parameter(self.element_id, param_name, param)
        return param

    def _set_parameter(self, param_name, value):
        if param_name in optix_dictionnary:
            pyoptix_param_name = optix_dictionnary[param_name]
        else:
            pyoptix_param_name = param_name
        param = self._get_c_parameter(param_name)
        if isinstance(value, dict):  # value is a dictionnary
            for key in value:
                assert key in ("value", "bounds", "multiplier", "type", "group", "flags")
                if key == "value":
                    param.value = DOUBLE(value[key])
                elif key == "multiplier":
                    param.multiplier = DOUBLE(value[key])
                elif key == "type":
                    param.type = INT(value[key])
                elif key == "group":
                    param.group = INT(value[key])
                elif key == "flags":
                    param.flags = c_ulong(value[key])
                else:
                    bounds = Bounds()
                    bounds.min = DOUBLE(value[key][0])
                    bounds.max = DOUBLE(value[key][1])
                    param.bounds = bounds
            set_parameter(self._element_id, param_name, param)
        elif isinstance(value, np.ndarray):  # value is an array
            assert param_name in array_parameter_list
            # print(f"setting {param_name} with array {value}")
            param.array = value
            # dump_arg_parameter(param)
        elif isinstance(value, (float, int)):  # value is a float
            try:
                value = float(value)
                param.value = DOUBLE(value)
            except TypeError:
                raise AttributeError(
                    "value of parameter must be castable in a float, an numpy.ndarray or a dictionnary")
        else:
            raise AttributeError("Unknown value type")
        if param_name in array_parameter_list:
            # print(f"Setting in _set_parameter {param_name} with {param} ")
            pass
        set_parameter(self._element_id, param_name, param)
        if "lineDensityCoeff_" in pyoptix_param_name:  # cas particulier des classes filles de Poly1D
            deg = int(pyoptix_param_name.split("_")[-1])
            max_deg = max(len(self._line_density_coeffs), deg)
            if max_deg > self.degree:
                self.degree = max_deg
            get_parameter(self._element_id, param_name, param)
        else:
            self.__getattribute__(pyoptix_param_name)  # update of internal variable
        updated_parameter = self._get_c_parameter(param_name)
        if updated_parameter.flags & (1 << 3):  # parameter is an array parameter
            # print(f"updated {param_name}:\n", updated_parameter)
            return updated_parameter.array
        else:
            return updated_parameter.value

    @property
    def element_id(self):
        return self._element_id

    @property
    def name(self):
        element_name = create_string_buffer(32)
        get_element_name(self._element_id, element_name, confirm=False)
        self._name = element_name.value.decode()
        return self._name

    @property
    def phi(self):
        self._phi = self._get_parameter("phi")
        return self._phi

    @phi.setter
    def phi(self, value):
        self._phi = self._set_parameter("phi", value)

    @property
    def psi(self):
        self._psi = self._get_parameter("psi")
        return self._psi

    @psi.setter
    def psi(self, value):
        self._psi = self._set_parameter("psi", value)

    @property
    def theta(self):
        self._theta = self._get_parameter("theta")
        return self._theta

    @theta.setter
    def theta(self, value):
        self._theta = self._set_parameter("theta", value)

    @property
    def d_phi(self):
        self._d_phi = self._get_parameter("Dphi")
        return self._d_phi

    @d_phi.setter
    def d_phi(self, value):
        self._d_phi = self._set_parameter("Dphi", value)

    @property
    def d_psi(self):
        self._d_psi = self._get_parameter("Dpsi")
        return self._d_psi

    @d_psi.setter
    def d_psi(self, value):
        self._d_psi = self._set_parameter("Dpsi", value)

    @property
    def d_theta(self):
        self._d_theta = self._get_parameter("Dtheta")
        return self._d_theta

    @d_theta.setter
    def d_theta(self, value):
        self._d_theta = self._set_parameter("Dtheta", value)

    @property
    def d_x(self):
        self._d_x = self._get_parameter("DX")
        return self._d_x

    @d_x.setter
    def d_x(self, value):
        self._d_x = self._set_parameter("DX", value)

    @property
    def d_y(self):
        self._d_y = self._get_parameter("DY")
        return self._d_y

    @d_y.setter
    def d_y(self, value):
        self._d_y = self._set_parameter("DY", value)

    @property
    def d_z(self):
        self._d_z = self._get_parameter("DZ")
        return self._d_z

    @d_z.setter
    def d_z(self, value):
        self._d_z = self._set_parameter("DZ", value)

    @property
    def distance_from_previous(self):
        self._distance_from_previous = self._get_parameter("distance")
        return self._distance_from_previous

    @distance_from_previous.setter
    def distance_from_previous(self, value):
        self._distance_from_previous = self._set_parameter("distance", value)

    @property
    def previous(self):
        return self._previous

    @previous.setter
    def previous(self, previous_oe):
        if previous_oe is not None:
            chain_element_by_id(previous_oe.element_id, self._element_id)
            self._previous = previous_oe
        else:
            self._previous = None

    @property
    def next(self):
        return self._next

    @next.setter
    def next(self, next_oe):
        self._next = next_oe
        if next_oe is not None:
            chain_element_by_id(self._element_id, next_oe.element_id)
            next_oe.previous = self
        else:
            chain_element_by_id(self._element_id, 0)

    @property
    def recording_mode(self):
        return self._recording_mode

    @recording_mode.setter
    def recording_mode(self, value):
        assert value in [RecordingMode.recording_output, RecordingMode.recording_input,
                         RecordingMode.recording_none]
        set_recording(self.element_id, value)
        self._recording_mode = value

    def __repr__(self):
        description = f"Element {self._name} of class {self.__class__}"
        if self._previous is not None:
            description += f"\n\t at {self._distance_from_previous} m from {self._previous.name}"
        else:
            description += f"\n\t at {self._distance_from_previous} m from None"
        if self._next is not None:
            description += f"\n\t pointing to {self._next.name}"
        else:
            description += f"\n\t pointing to None"
        description += f"\n\t oriented in pitch at {self._theta / degree} deg " \
                       f"(deviation {180 - 2 * self._theta / degree} deg)"
        description += f"\n\t oriented in roll at {self._phi / degree} deg"
        description += f"\n\t oriented in yaw at {self._psi / degree} deg\n"
        return description

    def get_local_frame(self, verbose=0, degrees=True):
        """
        Returns the surface frame orientation and position of the given element in the absolute frame

        This method return an array of four vectors of size 3, in the location pointed by frame_vectors.
        They are respectively and in this order, the unit vectors in X, Y, and Z, followed by the position of the
        frame origin in the absolute frame, where Z in normal to the surface at it center, Y tangent to the length
        of the surface, X tangent to the width of the surface.
        The PyOptix reference frame (X,Y,Z) translates in the SOLEIL
        frame as (-X, Z, S), ie the SOLEIL frame (S, X, Z) has a PyOptix matrix of [[0,0,1], [-1,0,0], [0,1,0]]

        :return: dictionnary containing the coordinate of each local axis of the surface in the reference frame,
            a 3 float array describing its center and the pitch, roll and yaw angle of the surface frame relative to
            rotations around X, Y and Z respectively
        :rtype: dict
        """
        vectors = get_surface_frame(self.element_id)
        x_vector, y_vector, z_vector, pos_vector = vectors
        mirror_frame = Rotation.from_matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]]).inv()
        surface_frame = Rotation.from_matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        mat = np.array([x_vector, y_vector, z_vector])
        if "film" in self._element_type.lower() or "source" in self._element_type.lower():
            pitch, roll, yaw = Rotation.from_matrix(mat).as_euler("xyz", degrees=degrees)
            mat = surface_frame.apply(mat)
        else:
            mat = mirror_frame.apply(mat)
            roll, pitch, yaw = Rotation.from_matrix(mat).as_euler("xyz", degrees=degrees)
        mat = Rotation.from_matrix(mat)
        if verbose:
            print("Orientation of", self.name)
            print(list(zip(["roll", "pitch", "yaw"], [roll, pitch, yaw])))
            print("X:", list(zip("SXZ", mat.as_matrix()[0])))
            print("Y:", list(zip("SXZ", mat.as_matrix()[1])))
            print("Z:", list(zip("SXZ", mat.as_matrix()[2])))
        return {"X_pyoptix": x_vector, "Y_pyoptix": y_vector, "Z_pyoptix": z_vector, "Center_pyoptix": pos_vector,
                "X_soleil": mat.as_matrix()[0], "Y_soleil": mat.as_matrix()[1], "Z_soleil": mat.as_matrix()[2],
                "Center_soleil": surface_frame.apply(pos_vector),
                "pitch": pitch, "yaw": yaw, "roll": roll}

    def get_exit_frame(self, verbose=0, degrees=True):
        """
        Returns the beam frame orientation and position after the given element in the absolute frame

        This method return an array of four vectors of size 3, in the location pointed by frame_vectors.
        They are respectively and in this order, the unit vectors in X, Y, and Z, followed by the position of the
        frame origin in the absolute frame, where Z in normal to the surface at it center, Y tangent to the length
        of the surface, X tangent to the width of the surface.
        The PyOptix reference frame (X,Y,Z) translates in the SOLEIL
        frame as (-X, Z, S), ie the SOLEIL frame (S, X, Z) has a PyOptix matrix of [[0,0,1], [-1,0,0], [0,1,0]]

        :return: dictionnary containing the coordinate of each local axis of the exit frame in the reference frame,
            a 3 float array describing its center and the pitch, roll and yaw angle of the surface frame relative to
            rotations around X, Y and Z respectively
        :rtype: dict
        """
        vectors = get_exit_frame(self.element_id)
        x_vector, y_vector, z_vector, pos_vector = vectors
        # mirror_frame = Rotation.from_matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]]).inv()
        surface_frame = Rotation.from_matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
        mat = np.array([x_vector, y_vector, z_vector])
        # if "film" in self._element_type.lower() or "source" in self._element_type.lower():
        pitch, roll, yaw = Rotation.from_matrix(mat).as_euler("xyz", degrees=degrees)
        mat = surface_frame.apply(mat)
        # else:
        #     mat = mirror_frame.apply(mat)
        #     roll, pitch, yaw = Rotation.from_matrix(mat).as_euler("xyz", degrees=degrees)
        mat = Rotation.from_matrix(mat)
        if verbose:
            print("Orientation exit frame after", self.name)
            print(list(zip(["roll", "pitch", "yaw"], [roll, pitch, yaw])))
            print("X:", list(zip("SXZ", mat.as_matrix()[0])))
            print("Y:", list(zip("SXZ", mat.as_matrix()[1])))
            print("Z:", list(zip("SXZ", mat.as_matrix()[2])))
        return {"X_pyoptix": x_vector, "Y_pyoptix": y_vector, "Z_pyoptix": z_vector, "Center_pyoptix": pos_vector,
                "X_soleil": mat.as_matrix()[0], "Y_soleil": mat.as_matrix()[1], "Z_soleil": mat.as_matrix()[2],
                "Center_soleil": surface_frame.apply(pos_vector),
                "pitch": pitch, "yaw": yaw, "roll": roll}

    def show_diagram(self, distance_from_oe=0, map_xy=False, light_xy=False,
                     light_xxp=False, light_yyp=False, show_first_rays=False, display="all", show_spd=True, **kwargs):
        """
        If recording_mode is set to any other value than RecordingMode.recording_none, displays X vs Y,
        X' vs X and Y' vs Y scatter plots at a distance `distance_from_oe` from the element.
        For quick display, light_xy, light_xxp and light_yyp can be set to True.
        For more realistic rendering, map_xy can be set to True.

        :param distance_from_oe: distance from the element at which to draw the spot diagram in m
        :type distance_from_oe: str
        :param map_xy: set to True for hexagonal pixel rendering of X vs Y diagram
        :type map_xy: bool
        :param light_xy: set to True for quick monochromatic rendering of X vs Y scatter plot
        :type light_xy: bool
        :param light_xxp: set to True for quick monochromatic rendering of X' vs X scatter plot
        :type light_xxp: bool
        :param light_yyp: set to True for quick monochromatic rendering of Y' vs Y scatter plot
        :type light_yyp: bool
        :param show_first_rays: set to True for a table display of the parameter of the first rays in diagram
        :param display: Indicates which representation is to be shown, can be "xy", "xy, xxp"... or "all"
        :type display: str
        :param kwargs: See pyoptix.ui_objects.plot_spd doc for additional parameters
        :return: tuple of (spd_data, handle(s) of the figures)
        :rtype: tuple
        """
        beamline_name = None
        chain_name = None
        if self.beamline is not None:
            chain_name = self.beamline.active_chain_name
            beamline_name = self.beamline.name
        print(beamline_name, chain_name)
        spots = self.get_diagram(distance_from_oe=distance_from_oe, show_first_rays=show_first_rays)
        spots = spots.loc[spots['Intensity'] != 0]
        if not spots.loc[spots['Intensity'] == 0].empty:
            print("Warning : la colorisation présente des effets de bords avec les ouvertures")
        if spots["X"].std() == 0:
            print("Warning empty diagram:\n", spots)
            return None, None
        if self.recording_mode == RecordingMode.recording_output:
            oe_name = "after " + self._name
        else:
            oe_name = "before " + self._name
        datasources = {"xy": ColumnDataSource(spots), "xxp": ColumnDataSource(spots), "yyp": ColumnDataSource(spots)}
        figs = []
        if display == "all":
            display = "xy, xxp, yyp"
        if "xy" in display:
            figs.append(plot_spd(datasources["xy"], x_key="X", y_key="Y", light_plot=light_xy, show_map=map_xy,
                                 beamline_name=beamline_name, chain_name=chain_name, oe_name=oe_name, **kwargs))
        if "xxp" in display:
            figs.append(
                plot_spd(datasources["xxp"], x_key="X", y_key="dX", light_plot=light_xxp, beamline_name=beamline_name,
                         chain_name=chain_name, oe_name=oe_name, **kwargs))

            def fun(x):
                return (spots["X"] + x * spots["dX"]).std()

            try:
                res = minimize(fun, 0)
                print(f"Optimal focalisation along X at {res.x[0]} m from this plane")
            except RuntimeError:
                print("Unable to find the optimal focalisation along X automatically")
        if "yyp" in display:
            figs.append(
                plot_spd(datasources["yyp"], x_key="Y", y_key="dY", light_plot=light_yyp, beamline_name=beamline_name,
                         chain_name=chain_name, oe_name=oe_name, **kwargs))

            def fun(x):
                return (spots["Y"] + x * spots["dY"]).std()

            try:
                res = minimize(fun, 0)
                print(f"Optimal focalisation along Y at {res.x[0]} m from this plane")
            except RuntimeError:
                print("Unable to find the optimal focalisation along Y automatically")
        handles = []
        if show_spd:
            for fig in figs:
                handles.append(show(fig, notebook_handle=True))
        return datasources, figs

    def get_diagram(self, distance_from_oe=0, show_first_rays=False):
        """
        If recording_mode is set to any other value than RecordingMode.recording_none, returns a
        (X, Y, dX, dY, Lambda) pandas dataframe where each row is a computed ray

        :param distance_from_oe: distance from the element at which to draw the spot diagram in m
        :type distance_from_oe: str
        :param show_first_rays: set to True for a table display of the parameter of the first rays in diagram
        :type show_first_rays: bool
        :return: pandas.Dataframe containing all rays intercept on the optics surface
        """
        assert self.recording_mode != RecordingMode.recording_none
        nrays = self.beamline.n_generated_rays
        diagram = Diagram(ndim=6, nreserved=int(nrays))
        get_spot_diagram(self.element_id, diagram, distance_from_oe)
        spots = pd.DataFrame(np.copy(np.ctypeslib.as_array(diagram.spots, shape=(diagram.reserved, diagram.dim))),
                             columns=("X", "Y", "dX", "dY", "Lambda", "Intensity"))
        if show_first_rays:
            print(spots.head())
        return spots

    def get_impacts(self, reference_frame="surface_frame", show_first_rays=False):
        """
        If recording_mode is set to any other value than RecordingMode.recording_none, returns a
        (X, Y, Z, dX, dY, dZ, Lambda, Intensity) pandas dataframe where each row is a computed ray in the frame
        of reference given in parameter reference_frame either :

        - "general_frame" : Absolute laboratory frame
        - "local_absolute_frame" : Absolute frame with origin on the surface
        - "aligned_local_frame" : Local frame, with origin on the surface, axe OZ is along the chief ray and OY is in
          the deviation plane of the last preceding reflective element. Transmissive elements do not change the
          AlignedLocalFrame
        - "surface_frame" : Local frame used to describe a surface. Origin is at surface intercept with the chief ray.
          Oz is along the surface normal (at origin). OX is the tangential axis for reflective elements.

        Method to be used for computing footprints on a mirror.


        :param reference_frame: reference frame for coordinates see above
        :type reference_frame: str
        :param show_first_rays: set to True for a table display of the parameter of the first rays in diagram
        :type show_first_rays: bool
        :return: pandas.Dataframe containing all rays intercept on the optics surface
        """
        assert self.recording_mode != RecordingMode.recording_none
        nrays = self.beamline.n_generated_rays
        diagram = Diagram(ndim=8, nreserved=int(nrays))
        frame = {"general_frame": FrameID.general_frame,
                 "local_absolute_frame": FrameID.local_absolute_frame,
                 "aligned_local_frame": FrameID.aligned_local_frame,
                 "surface_frame": FrameID.surface_frame}
        get_impacts_data(self.element_id, diagram, frame[reference_frame])
        spots = pd.DataFrame(np.copy(np.ctypeslib.as_array(diagram.spots, shape=(diagram.reserved, diagram.dim))),
                             columns=("X", "Y", "Z", "dX", "dY", "dZ", "Lambda", "Intensity"))
        if show_first_rays:
            print(spots.head())
        return spots

    def show_impacts(self, reference_frame="local_absolute_frame", plots=None, **kwargs):
        """
        If recording_mode is set to any other value than RecordingMode.recording_none, shows by default the (X,Y) and
        (X,Z) cross sections of the (X, Y, Z, dX, dY, dZ, Lambda, Intensity) pandas dataframe where each row is a
        computed ray in the frame of reference given in parameter reference_frame either :

        - "general_frame" : Absolute laboratory frame
        - "local_absolute_frame" : Absolute frame with origin on the surface
        - "aligned_local_frame" : Local frame, with origin on the surface, axe OZ is along the chief ray and OY is in
          the deviation plane of the last preceding reflective element. Transmissive elements do not change the
          AlignedLocalFrame
        - "surface_frame" : Local frame used to describe a surface. Origin is at surface intercept with the chief ray.
          Oz is along the surface normal (at origin). OX is the tangential axis for reflective elements.

        Method to be used for computing footprints on a mirror.

        :param plots: Plots to be displayed. Default is [("X", "Y"), ("X", "Z")]
        :type plots: list of tuple
        :param nrays: number of expected rays (default: source_oe.nrays). Only use if generate is called multiple times.
        :type nrays: int
        :param reference_frame: reference frame for coordinates see above
        :type reference_frame: str
        :param kwargs: See pyoptix.ui_objects.plot_spd doc for additional parameters
        :return: tuple of (spd_data, handle(s) of the figures)
        :rtype: tuple
        """
        if plots is None:
            plots = [('X', 'Y'), ('X', 'Z')]
        impacts = self.get_impacts(reference_frame=reference_frame, show_first_rays=False)
        impacts = impacts.loc[impacts['Intensity'] != 0]
        if not impacts.loc[impacts['Intensity'] == 0].empty:
            print("Warning : plot colorization will show artefacts around borders if used with apertures")

        figs = []
        coldatasource = ColumnDataSource(impacts)
        for plot in plots:
            x_axis, y_axis = plot
            figs.append(plot_spd(coldatasource, x_key=x_axis, y_key=y_axis,
                                 beamline_name=self.beamline.name, chain_name=self.beamline.active_chain_name,
                                 oe_name="on " + self._name, **kwargs))
        handles = []
        for fig in figs:
            handles.append(show(fig, notebook_handle=True))
        return coldatasource, figs

    def export_beam_mcpl(self, reference_frame):
        """
        If recording_mode is set to any other value than RecordingMode.recording_none, exports a
        MCPL (see https://mctools.github.io/mcpl/) in the frame
        of reference given in parameter reference_frame (see doc of OpticalElement.get_impacts)

        :param reference_frame:reference frame for coordinates see  OpticalElement.get_impacts
        :type reference_frame: str
        :return: None
        :rtype: Nonetype
        """
        mcpl_file_out = PyOptixMCPLWriter("my_mcpl_file.mcpl")
        mcpl_file_out.add_comment(f"File generated with beamline {self.beamline.name}")
        mcpl_file_out.add_comment(f"Test MCPL file from MCPL_interfacing.ipynb")
        mcpl_file_out.dump_diagram(self.get_impacts(reference_frame=reference_frame))
        mcpl_file_out.write_to_file()

    def from_element_id(self, element_id, print_all=False):
        """
        Sets all the optical element parameters using a preexisting optix object which handle is passed as `element_id`.

        :param element_id: Handle of the preexisting optix objects
        :type element_id: wintypes.INT
        :param print_all: set to true to print all values of parameters
        :type print_all: bool
        :return: None
        """
        hparam = HANDLE(0)
        self._element_id = element_id
        print("initializing ", self.name)
        param = Parameter()
        param_name = create_string_buffer(48)
        enumerate_parameters(element_id, hparam, param_name, param, confirm=False)
        while hparam:
            if print_all:
                print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
                            f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
            if param_name.value.decode() in optix_dictionnary:
                self.__dict__["_" + optix_dictionnary[param_name.value.decode()]] = param.value
            else:
                self.__dict__["_" + param_name.value.decode()] = param.value
            enumerate_parameters(element_id, hparam, param_name, param, confirm=False)
        if print_all:
            print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
                        f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
        if param_name.value.decode() in optix_dictionnary:
            self.__dict__["_" + optix_dictionnary[param_name.value.decode()]] = param.value
        else:
            self.__dict__["_" + param_name.value.decode()] = param.value
        enumerate_parameters(element_id, hparam, param_name, param, confirm=False)
        next_id = get_next_element(element_id, confirm=False)
        if next_id is not None:
            next_name = create_string_buffer(32)
            get_element_name(next_id, next_name, confirm=False)
            self._next = next_name.value.decode()
        else:
            self._next = None
        previous_id = get_previous_element(element_id, confirm=False)
        if previous_id is not None:
            previous_name = create_string_buffer(32)
            get_element_name(previous_id, previous_name, confirm=False)
            self._previous = previous_name.value.decode()
        else:
            self._previous = None

    def set_recording(self, recording_mode="output"):
        """
        Sets the recording mode for the spot diagram associated with the element.
        Recording mode can be "output" (records rays affected by the element), "input" (records ray just before
        reaching the element), "not_recording" (does not record the rays).

        :param recording_mode: When to record the rays in reference to the element (see above)
        :type recording_mode: str
        :return: None
        """
        self._recording_mode = {"output": RecordingMode.recording_output,
                                "input": RecordingMode.recording_input,
                                "not_recording": RecordingMode.recording_none}[recording_mode]

    def set_transmissive(self, is_transmissive):
        """
        Sets an optical element as transmissive or reflective depending on the value of is_transmissive.
        Only authorized for gratings

        :param is_transmissive: if True, makes the optical element transmissive, if False, reflective
        :type is_transmissive: bool
        :return: None
        :rtype: Nonetype
        """
        assert "grating" in self._element_type.lower()
        set_transmissive(self._element_id, is_transmissive)

    def get_transmissive(self):
        """
        Returns True if optical element is transmissive, False if reflective.

        :return: True if optical element is transmissive, False if reflective
        :rtype: bool
        """
        return bool(get_transmissive(self._element_id))

    def dump_properties(self, verbose=1):
        """
        Prints all stored optix parameters
        """
        hparam = HANDLE(0)
        param_name = create_string_buffer(48)
        param = Parameter()
        enumerate_parameters(self.element_id, hparam, param_name, param, confirm=False)
        if verbose:
            print(f"Dump of all optix parameter for {self.name}")
        param_dict = {"oe_name": self.name, "oe_params": {}, "oe_class": self.__class__}
        while hparam:
            if param_name.value.decode() not in array_parameter_list:
                value = param.value
            else:
                pa = ParamArray()

                ppa_src = param.p_array
                cast(ppa_src, POINTER(ParamArray))
                pa_src = ppa_src.contents
                pa.dims[0] = pa_src.dims[0]
                pa.dims[1] = pa_src.dims[1]
                size = pa.dims[0] * pa.dims[1]
                databuf = np.fromiter(pa_src.data, np.float64, count=size)
                value = np.ndarray((pa.dims[1], pa.dims[0]), dtype=np.float64, buffer=databuf)
            if verbose:
                print("\t", f"{param_name.value.decode()}: {value} [{param.bounds.min}, {param.bounds.max}],"
                            f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
            param_dict["oe_params"][param_name.value.decode()] = {"value": value,
                                                                  "bounds": [param.bounds.min, param.bounds.max],
                                                                  "multiplier": param.multiplier,
                                                                  "type": param.type, "group": param.group,
                                                                  "flags": param.flags}
            enumerate_parameters(self.element_id, hparam, param_name, param, confirm=False)

        if param_name.value.decode() not in array_parameter_list:
            value = param.value
        else:
            value = "array"
        if verbose:
            print("\t", f"{param_name.value.decode()}: {value} [{param.bounds.min}, {param.bounds.max}],"
                        f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
        param_dict["oe_params"][param_name.value.decode()] = {"value": value,
                                                              "bounds": [param.bounds.min, param.bounds.max],
                                                              "multiplier": param.multiplier,
                                                              "type": param.type, "group": param.group,
                                                              "flags": param.flags}
        return param_dict

    def display_properties(self):
        display_parameter_sheet(self)

    def enumerate_stops(self, max_vertex=8):
        """
        Any Element based on a Surface class has an ApertureStop member, which is composed of an array of superimposed
        stopping Region Regions or Stops;

        - Each Region is a simply connected domain, which can be either a Polygon or an Ellipse.
            Mind that no check is made on the simply connected character of a polygon when created
        - Each region defines the opacity of its inside area and imposes it over all underlying regions. *
        - It doesn't define the opacity of its outside area unless it is the bottom element of the stack. *
            - if the region Opacity=true, the optical transmittance is 0 inside the region
            - if the region Opacity=false, the optical transmittance is 1. inside the region
        - The opacity of the intersection of all outside areas of all region is the opposite of the opacity of the
            bottom region (first define in the stack)

        In other words, the first stop defines an opacity in all space. Adding an opaque Region obstructs the
        transmission under the region area, but does not change the outside
        while adding a transparent Region opens an aperture though all the stacked stops &nd still does not change
        the outside


        - The logical stack should define reducing region sizes from bottom to top. This is intended to allow rounding
        of squared apertures and stops.

        All functions of the Aperture API are returning a size_t value. A negative return value, actually -1, means an
        error occurred and the GetOptiXLastError can be checked for a reason.

        :param max_vertex: number of expected vertexes to be looked up in memory, 8 for a rectangle. Must be even.
            Default: 8
        :type max_vertex: int
        :return: Dictionnary containing the regions describing the aperture and their parameters
        :rtype: dict
        """
        assert max_vertex % 2 == 0, AttributeError("max_vertex must be even")
        stops = {}
        stops_kind = enumerate_stops(self._element_id)
        for i, kind in enumerate(stops_kind):
            stops_details = {}
            if kind[1] == "Polygon":
                stops_details["kind"] = kind[1]
                vertex, opacity = self.get_polygon_parameters(i, max_vertex)
                vertex_x = vertex[::2]
                vertex_y = vertex[1::2]
                stops_details["vertex"] = list(zip(vertex_x, vertex_y))
                stops_details["opacity"] = opacity
            elif kind[1] == "Ellipse":
                stops_details["kind"] = kind[1]
                x_axis, y_axis, opacity, x_center, y_center, angle = self.get_elliptical_stop_parameters(i)
                stops_details["x_axis"] = x_axis
                stops_details["y_axis"] = y_axis
                stops_details["x_center"] = x_center
                stops_details["y_center"] = y_center
                stops_details["angle"] = angle
                stops_details["opacity"] = opacity
            else:
                stops_details["kind"] = kind[1]
            stops[str(i)] = stops_details
        return stops

    def set_aperture_activity(self, active=True):
        """
        Set the activity state of the aperture, i.e. whether the aperture is used for propagation or not.

        Parameters:
        active (bool): The activity state of the aperture. Default is True.
        """
        set_aperture_activity(self.element_id, active)

    def get_aperture_activity(self):
        """
        Get the activity state of the aperture, i.e. whether the aperture is used for propagation or not.

        Returns:
        bool: The activity state of the aperture.
        """
        return get_aperture_activity(self.element_id)

    def add_rectangular_stop(self, x_width, y_width, opacity=0, x_center=0, y_center=0, angle=0):
        """
        Add a rectangular region to the aperture. X axis is along the length of the optical element.

        Parameters:
        x_width (float): Width of the rectangular region.
        y_width (float): Height of the rectangular region.
        opacity (float): Opacity of the rectangular stop. Default is 0 meaning transparent inside the rectangle,
            absorbant outside of it.
        x_center (float): X-coordinate of the center of the rectangular stop. Default is 0.
        y_center (float): Y-coordinate of the center of the rectangular stop. Default is 0.
        angle (float): Rotation angle of the rectangular stop in radians. Default is 0.
        """
        add_rectangular_stop(self.element_id, x_width, y_width, opacity, x_center, y_center, angle)

    def get_polygon_parameters(self, index, array_width):
        """
        Get the parameters of a polygonal stop in the aperture.

        Parameters:
            index (int): The index of the polygonal stop in the list of apertures, see enumerate_stops.
            array_width (int): The width of the array storing the stops.

        Returns:
            tuple: A tuple containing the parameters of the polygonal stop and the opacity of the polygon.
        """
        return get_polygon_parameters(self.element_id, index, array_width)

    def add_polygonal_stop(self, num_sides, vertex_array, opacity=0):
        """
        Add a polygonal stop to the aperture.

        Parameters:
            num_sides (int): The number of sides of the polygon.
            vertex_array (array-like): The vertices of the polygon.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
        """
        add_polygonal_stop(self.element_id, num_sides, vertex_array, opacity)

    def insert_polygonal_stop(self, index, num_sides, vertex_array, opacity=0):
        """
        Insert a polygonal stop at the specified index in the aperture.

        Parameters:
            index (int): The index at which the polygonal stop should be inserted. See enumerate.stops.
            num_sides (int): The number of sides of the polygon.
            vertex_array (array-like): The vertices of the polygon.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
        """
        insert_polygonal_stop(self.element_id, index, num_sides, vertex_array, opacity)

    def replace_stop_by_polygon(self, index, num_sides, vertex_array, opacity=0):
        """
        Replace a stop in the aperture with a polygonal stop.

        Parameters:
            index (int): The index of the stop to be replaced, see enumerate_stops.
            num_sides (int): The number of sides of the polygon.
            vertex_array (array-like): The vertices of the polygon.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
        """
        replace_stop_by_polygon(self.element_id, index, num_sides, vertex_array, opacity)

    def insert_rectangular_stop(self, index, x_width, y_width, opacity=0, x_center=0, y_center=0, angle=0):
        """
        Insert a rectangular stop at the specified index in the aperture.

        Parameters:
            index (int): The index at which the rectangular stop should be inserted, see enumerate_stops.
            x_width (float): The width of the rectangular stop along the x-axis.
            y_width (float): The width of the rectangular stop along the y-axis.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
            x_center (float): The x-coordinate of the center of the stop. Default is 0.
            y_center (float): The y-coordinate of the center of the stop. Default is 0.
            angle (float): The rotation angle of the stop in radians. Default is 0.
        """
        insert_rectangular_stop(self.element_id, index, x_width, y_width, opacity, x_center, y_center, angle)

    def replace_stop_by_rectangle(self, index, x_width, y_width, opacity=0, x_center=0, y_center=0, angle=0):
        """
        Replace a stop in the aperture with a rectangular stop.

        Parameters:
            index (int): The index of the stop to be replaced, see enumerate_stops.
            x_width (float): The width of the rectangular stop along the x-axis.
            y_width (float): The width of the rectangular stop along the y-axis.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
            x_center (float): The x-coordinate of the center of the stop. Default is 0.
            y_center (float): The y-coordinate of the center of the stop. Default is 0.
            angle (float): The rotation angle of the stop in radians. Default is 0.
        """
        replace_stop_by_rectangle(self.element_id, index, x_width, y_width, opacity, x_center, y_center, angle)

    def get_elliptical_stop_parameters(self, index):
        """
        Get the parameters of an elliptical stop in the aperture.

        Parameters:
            index (int): The index of the elliptical stop, see enumerate_stops.

        Returns:
            tuple: A tuple containing the parameters of the elliptical stop : boolean=1 if parameters have been
                successfully retrieved + parameters in the order defined in add_elliptical_stop
        """
        return get_ellipse_parameters(self.element_id, index)

    def add_elliptical_stop(self, x_axis, y_axis, opacity=0, x_center=0, y_center=0, angle=0):
        """
        Add an elliptical stop to the aperture.

        Parameters:
            x_axis (float): The length of the major axis of the ellipse.
            y_axis (float): The length of the minor axis of the ellipse.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
            x_center (float): The x-coordinate of the center of the stop. Default is 0.
            y_center (float): The y-coordinate of the center of the stop. Default is 0.
            angle (float): The rotation angle of the stop in radians. Default is 0.
        """
        add_elliptical_stop(self.element_id, x_axis, y_axis, opacity, x_center, y_center, angle)

    def insert_elliptical_stop(self, index, x_axis, y_axis, opacity=0, x_center=0, y_center=0, angle=0):
        """
        Insert an elliptical stop at the specified index in the aperture.

        Parameters:
            index (int): The index at which the elliptical stop should be inserted, see enumerate_stops.
            x_axis (float): The length of the major axis of the ellipse.
            y_axis (float): The length of the minor axis of the ellipse.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
            x_center (float): The x-coordinate of the center of the stop. Default is 0.
            y_center (float): The y-coordinate of the center of the stop. Default is 0.
            angle (float): The rotation angle of the stop in radians. Default is 0.
        """
        insert_elliptical_stop(self.element_id, index, x_axis, y_axis, opacity, x_center, y_center, angle)

    def replace_stop_by_ellipse(self, index, x_axis, y_axis, opacity=0, x_center=0, y_center=0, angle=0):
        """
        Replace a stop in the aperture with an elliptical stop.

        Parameters:
            index (int): The index of the stop to be replaced, see enumerate_stops.
            x_axis (float): The length of the major axis of the ellipse.
            y_axis (float): The length of the minor axis of the ellipse.
            opacity (float): The opacity of the stop.  Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
            x_center (float): The x-coordinate of the center of the stop. Default is 0.
            y_center (float): The y-coordinate of the center of the stop. Default is 0.
            angle (float): The rotation angle of the stop in radians. Default is 0.
        """
        replace_stop_by_ellipse(self.element_id, index, x_axis, y_axis, opacity, x_center, y_center, angle)

    def add_circular_stop(self, radius, opacity=0, x_center=0, y_center=0):
        """
        Add a circular stop to the aperture. Note a circular stop is of kind "Ellipse"

        Parameters:
            radius (float): The radius of the circular stop.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
            x_center (float): The x-coordinate of the center of the stop. Default is 0.
            y_center (float): The y-coordinate of the center of the stop. Default is 0.
        """

        add_circular_stop(self.element_id, radius, opacity, x_center, y_center)

    def insert_circular_stop(self, index, radius, opacity=0, x_center=0, y_center=0):
        """
        Insert a circular stop at the specified index in the aperture.

        Parameters:
            index (int): The index at which the circular stop should be inserted.
            radius (float): The radius of the circular stop.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
            x_center (float): The x-coordinate of the center of the stop. Default is 0.
            y_center (float): The y-coordinate of the center of the stop. Default is 0.
        """
        insert_circular_stop(self.element_id, index, radius, opacity, x_center, y_center)

    def replace_stop_by_circle(self, index, radius, opacity=0, x_center=0, y_center=0):
        """
        Replace a stop in the aperture with a circular stop.

        Parameters:
            index (int): The index of the stop to be replaced, see enumerate_stops.
            radius (float): The radius of the circular stop.
            opacity (float): The opacity of the stop. Default is 0 meaning transparent inside the rectangle,
                absorbant outside of it.
            x_center (float): The x-coordinate of the center of the stop. Default is 0.
            y_center (float): The y-coordinate of the center of the stop. Default is 0.
        """
        replace_stop_by_circle(self.element_id, index, radius, opacity, x_center, y_center)

    def draw_aperture(self, max_vertex=8):
        """
        Draw all the sub region of the aperture, the ones that are absorbing inside are filled with black.

        Parameters:
            max_vertex (int): The maximum number of vertices to look up in memory for polygons. Default is 8.
        """
        stops = self.enumerate_stops(max_vertex)
        fig = plot_aperture(stops, title=f"Composite aperture of {self.name}")
        fig.show()


class Source(OpticalElement):
    """
    Base class for all Sources elements. Can be used as is with the correct element type or using the inherited
    classes
    """

    def __init__(self, **kwargs):
        """
        Constructor of a Source type element. Inherits OpticalElement class.
        Acceptable element_types are "GaussianSource", "RadialSource" or
        "XYGridSource"
        :param kwargs: See OpticalElement doc for other parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ["GaussianSource", "RadialSource", "XYGridSource",
                                              "AstigmaticGaussianSource", "BMtypeGaussianSource"]
        super().__init__(**kwargs)


class BMSource(OpticalElement):
    """
    Class for Bending magnet sources
    """

    def __init__(self, sigma_x=0, sigma_y=0, aperture_x=0, sigma_y_div=0, trajectory_radius=5.36, nrays=0, **kwargs):
        """
        Constructor of a gaussian source type element of size sigma_x*sigma_y and divergence aperture_x*sigma_y_div.
        Will generate nrays rays when method generate is called.

        :param sigma_x: RMS source size in X direction in m
        :type sigma_x: float
        :param sigma_y: RMS source size in Y direction in m
        :type sigma_y: float
        :param aperture_x: collection angle in the deflection plane
        :type aperture_x: float
        :param sigma_y_div: RMS source divergence in Y direction
        :type sigma_y_div: float
        :param nrays: number of rays to be generated
        :type nrays: int
        :param kwargs: See OpticalElement doc for other parameters
        """
        kwargs["element_type"] = "BMtypeGaussianSource"
        super().__init__(**kwargs)
        self.sigma_x = sigma_x
        self.sigma_y = sigma_y
        self.aperture_x = aperture_x
        self.sigma_y_div = sigma_y_div
        self.nrays = nrays
        self.trajectory_radius = trajectory_radius

    @property
    def nrays(self):
        self._nrays = self._get_parameter("nRays")
        return self._nrays

    @nrays.setter
    def nrays(self, value):
        self._nrays = self._set_parameter("nRays", value)

    @property
    def sigma_x(self):
        self._sigma_x = self._get_parameter("sigmaX")
        return self._sigma_x

    @sigma_x.setter
    def sigma_x(self, value):
        self._sigma_x = self._set_parameter("sigmaX", value)

    @property
    def sigma_y(self):
        self._sigma_y = self._get_parameter("sigmaY")
        return self._sigma_y

    @sigma_y.setter
    def sigma_y(self, value):
        self._sigma_y = self._set_parameter("sigmaY", value)

    @property
    def aperture_x(self):
        self._aperture_x = self._get_parameter("apertureX")
        return self._aperture_x

    @aperture_x.setter
    def aperture_x(self, value):
        self._aperture_x = self._set_parameter("apertureX", value)

    @property
    def sigma_y_div(self):
        self._sigma_y_div = self._get_parameter("sigmaYdiv")
        return self._sigma_y_div

    @sigma_y_div.setter
    def sigma_y_div(self, value):
        self._sigma_y_div = self._set_parameter("sigmaYdiv", value)

    @property
    def trajectory_radius(self):
        self._trajectory_radius = self._get_parameter("trajectoryRadius")
        return self._trajectory_radius

    @trajectory_radius.setter
    def trajectory_radius(self, value):
        self._trajectory_radius = self._set_parameter("trajectoryRadius", value)


class GaussianSource(Source):
    """
    Base class for all Sources elements. Can be used as is with the correct element type or using the inherited
    classes
    """

    def __init__(self, sigma_x=0, sigma_y=0, sigma_x_div=0, sigma_y_div=0, nrays=0, **kwargs):
        """
        Constructor of a gaussian source type element of size sigma_x*sigma_y and divergence sigma_x_div*sigma_y_div.
        Will generate nrays rays when method generate is called.

        :param sigma_x: RMS source size in X direction in m
        :type sigma_x: float
        :param sigma_y: RMS source size in Y direction in m
        :type sigma_y: float
        :param sigma_x_div: RMS source divergence in X direction
        :type sigma_x_div: float
        :param sigma_y_div: RMS source divergence in Y direction
        :type sigma_y_div: float
        :param nrays: number of rays to be generated
        :type nrays: int
        :param kwargs: See OpticalElement doc for other parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ["GaussianSource", "AstigmaticGaussianSource"]
        else:
            kwargs["element_type"] = "GaussianSource"
        super().__init__(**kwargs)
        self.sigma_x = sigma_x
        self.sigma_y = sigma_y
        self.sigma_x_div = sigma_x_div
        self.sigma_y_div = sigma_y_div
        self.nrays = nrays

    @property
    def nrays(self):
        self._nrays = self._get_parameter("nRays")
        return self._nrays

    @nrays.setter
    def nrays(self, value):
        self._nrays = self._set_parameter("nRays", value)

    @property
    def sigma_x(self):
        self._sigma_x = self._get_parameter("sigmaX")
        return self._sigma_x

    @sigma_x.setter
    def sigma_x(self, value):
        self._sigma_x = self._set_parameter("sigmaX", value)

    @property
    def sigma_y(self):
        self._sigma_y = self._get_parameter("sigmaY")
        return self._sigma_y

    @sigma_y.setter
    def sigma_y(self, value):
        self._sigma_y = self._set_parameter("sigmaY", value)

    @property
    def sigma_x_div(self):
        self._sigma_x_div = self._get_parameter("sigmaXdiv")
        return self._sigma_x_div

    @sigma_x_div.setter
    def sigma_x_div(self, value):
        self._sigma_x_div = self._set_parameter("sigmaXdiv", value)

    @property
    def sigma_y_div(self):
        self._sigma_y_div = self._get_parameter("sigmaYdiv")
        return self._sigma_y_div

    @sigma_y_div.setter
    def sigma_y_div(self, value):
        self._sigma_y_div = self._set_parameter("sigmaYdiv", value)


class AstigmaticGaussianSource(GaussianSource):
    """
    Base class for all Sources elements. Can be used as is with the correct element type or using the inherited
    classes
    """

    def __init__(self, waist_x=0, waist_y=0, **kwargs):
        """
        Constructor of a gaussian source type element of size sigma_x*sigma_y and divergence sigma_x_div*sigma_y_div.
        Will generate nrays rays when method generate is called.

        :param waist_x: distance of the horizontal source point from source plane in m
        :type waist_x: float
        :param waist_y: distance of the vertical source point from source plane in m
        :type waist_y: float
        :param kwargs: See GaussianSource doc for other parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "AstigmaticGaussianSource"
        else:
            kwargs["element_type"] = "AstigmaticGaussianSource"
        super().__init__(**kwargs)
        self.waist_x = waist_x
        self.waist_y = waist_y

    @property
    def waist_x(self):
        self._waist_x = self._get_parameter("waistX")
        return self._waist_x

    @waist_x.setter
    def waist_x(self, value):
        self._waist_x = self._set_parameter("waistX", value)

    @property
    def waist_y(self):
        self._waist_y = self._get_parameter("waistY")
        return self._waist_y

    @waist_y.setter
    def waist_y(self, value):
        self._waist_y = self._set_parameter("waistY", value)


class PlaneMirror(OpticalElement):
    """
    Class for plane mirrors. Inherits OpticalElement
    """

    def __init__(self, **kwargs):
        """
        Constructor for the PlaneMirror class
        :param kwargs: See OpticalElement doc for parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "PlaneMirror"
        else:
            kwargs["element_type"] = "PlaneMirror"
        super().__init__(**kwargs)


class PolynomialOpticalElement(OpticalElement):
    def __init__(self, **kwargs):
        """
        Constructor for the PlaneMirror class
        :param kwargs: See OpticalElement doc for parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ["NaturalPolynomialMirror", "LegendrePolynomialMirror"]
        else:
            raise AttributeError("Abstract class not to be instantiated")
        super().__init__(**kwargs)
        self.surface_limits = np.array([[-1, 1], [-1, 1]])
        self.coefficients = np.zeros((1, 1))

    @property
    def coefficients(self):
        self._coefficients = self._get_parameter("coefficients")
        return self._coefficients

    @coefficients.setter
    def coefficients(self, value):
        self._coefficients = self._set_parameter("coefficients", value)

    @property
    def surface_limits(self):
        self._surface_limits = self._get_parameter("surfaceLimits")
        return self._surface_limits

    @surface_limits.setter
    def surface_limits(self, value):
        self._surface_limits = self._set_parameter("surfaceLimits", value)

    def fit_heights(self, order_x, order_y, filename):
        """
        Fitting method to a XYZ file. The file must contain exactly 3 columns of values describing the XYZ surface.
        Values should be in meters.
        :param order_x: order of the polynomial fit in the X direction
        :type order_x: int
        :param order_y: order of the polynomial fit in the Y direction
        :type order_y: int
        :param filename: absolute filename to the file containing the height data
        :type filename: str
        :return: the rms residual of the fit
        :rtype: np.array
        """
        heights = np.loadtxt(filename)
        sigma_h = np.empty(heights.shape)
        fit_surface_to_heights(self.element_id, order_x, order_y, self.surface_limits, heights, sigma_h)
        return sigma_h

    def fit_slopes(self, order_x, order_y, filename):
        """
        Fitting method to a XYX'Y' file. The file must contain exactly 4 columns of values describing the
        X,Y,dZ/dX, dZ/dY slopes.
        Values should be in meters and radians.
        :param order_x: order of the polynomial fit in the X direction
        :type order_x: int
        :param order_y: order of the polynomial fit in the Y direction
        :type order_y: int
        :param filename: absolute filename to the file containing the height data
        :type filename: str
        :return: a tuple containing the maps of residuals of the fits along X and Y respectively
        :rtype: NoneType
        """
        slopes = np.loadtxt(filename)
        sigma_x = np.empty(slopes.shape)
        sigma_y = np.empty(slopes.shape)
        fit_surface_to_heights(self.element_id, order_x, order_y, self.limits, slopes, sigma_x, sigma_y)
        return sigma_x, sigma_y

    def plot_surface(self, nb_points_x=100, nb_points_y=100):
        x = np.linspace(self.limits[0], self.limits[1], nb_points_x)
        y = np.linspace(self.limits[2], self.limits[3], nb_points_y)
        zz = np.polynomial.polynomial.polygrid2d(x, y, self.coeffs)


class NaturalPolynomialMirror(PolynomialOpticalElement):
    """
    Class for mirrors whose surface is defined by a polynomial.
    """

    def __init__(self, **kwargs):
        """
        Constructor for the NaturalPolynomialMirror class
        :param kwargs: See OpticalElement doc for parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "NaturalPolynomialMirror"
        else:
            kwargs["element_type"] = "NaturalPolynomialMirror"
        super().__init__(**kwargs)


class LegendrePolynomialMirror(PolynomialOpticalElement):
    """
    Class for mirrors whose surface is defined by a polynomial.
    """

    def __init__(self, **kwargs):
        """
        Constructor for the LegendrePolynomialMirror class
        :param kwargs: See OpticalElement doc for parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "LegendrePolynomialMirror"
        else:
            kwargs["element_type"] = "LegendrePolynomialMirror"
        super().__init__(**kwargs)


class ConicalMirror(OpticalElement):
    """
        Class for conical mirrors. Inherits OpticalElement.

        Generates all cone-shaped mirrors using parameters curvature, directrix_angle and apex_distance such as:

        curvature is the curvature (1/Rc) of the directric circle passing through the origin point  in m-1
        directrix_angle  is the inclination angle of the directrix circle plane with respect to the YOZ plane in rad
        apex_distance is the distance from apex to the origin point

        The surface is positioned such as the X axis is one generatrix and the apex is the point of coordinates
        (-x0, 0,0), the directrix circle passing by the origin has a curvature 1/r and its plane is inclined by
        theta around the direction Y with 1/r= curvature, -x0=apex_distance.
        The cone is a revolution cone only if the apex is on the circle normal passing by its center,
        ie 1/curvature + apex_distance*sin(directix_angle) = 0

    """

    def __init__(self, curvature=1, directix_angle=0, apex_distance=1, **kwargs):
        """
        Constructor for the ConicalMirror class. For details of convention of conical parameter, see
        ConicalMirror class doc.

        :param curvature: curvature (1/Rc) of the directric circle passing through the origin point  in m-1
        :type curvature: float
        :param directix_angle: inclination angle of the directrix circle plane with respect to the YOZ plane in rad
        :type directix_angle: float
        :param apex_distance: distance from apex to the origin point
        :type apex_distance: float
        :param kwargs: See OpticalElement doc for parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ["ConicalMirror"]
        else:
            kwargs["element_type"] = "ConicalMirror"
        super().__init__(**kwargs)
        self.curvature = curvature
        self.directrix_angle = directix_angle
        self.apex_distance = apex_distance

    @property
    def curvature(self):
        self._curvature = self._get_parameter("curvature")
        return self._curvature

    @curvature.setter
    def curvature(self, value):
        self._curvature = self._set_parameter("curvature", value)

    @property
    def directrix_angle(self):
        self._directrix_angle = self._get_parameter("directrix_angle")
        return self._directrix_angle

    @directrix_angle.setter
    def directrix_angle(self, value):
        self._directrix_angle = self._set_parameter("directrix_angle", value)

    @property
    def apex_distance(self):
        self._apex_distance = self._get_parameter("apex_distance")
        return self._apex_distance

    @apex_distance.setter
    def apex_distance(self, value):
        self._apex_distance = self._set_parameter("apex_distance", value)

    def is_of_revolution(self):
        return 1 / self.curvature + self.apex_distance * np.sin(self.directix_angle) == 0


class RevolutionQuadricMirror(OpticalElement):
    """
    Class for quadric mirrors. Inherits OpticalElement. Inherited by classes "ConicBaseCylindricalMirror" and
    "RevolutionQuadricMirror"

    Generates all quadric surfaced mirrors using quadric parameters p, q and theta such as:

    p is defined at coordinate (p*cos(theta), -p*sin(theta))
    q is defined at coordinate (q*cos(theta), q*sin(theta))
    Due to this convention, the cylinder base (directrix) is

    - an ellipse if p^{−1}*q^{−1}<0
    - an hyperbola if p^{−1}*q^{−1}>0
    - a parabola if either p^{−1}=0 or q^{−1}=0
    - Warning : p^{−1}=q^{−1} is forbidden and will result as an error at any time.

    """

    def __init__(self, inverse_p=0, inverse_q=0.1, theta0=0, **kwargs):
        """
        Constructor for the RevolutionQuadricMirror class. For details of convention of quadric parameter, see
        RevolutionQuadricMirror class doc.

        Note: while p and q are the usual quadric parameter, here 1/p (inverse_p) and 1/q (inverse_q) parameters are
        instanciated so as to have continuity between all quadrics in an optimization problem.

        :param inverse_p: Inverse distance from center of the mirror to object focal point in m-1
        :type inverse_p: float
        :param inverse_q: Inverse distance from center of the mirror to image focal point in m-1
        :type inverse_q: float
        :param theta0: Angle of incidence at the center of the mirror in rad
        :type theta0: float
        :param kwargs: See OpticalElement doc for parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ["ConicBaseCylindricalMirror", "RevolutionQuadricMirror"]
        else:
            kwargs["element_type"] = "RevolutionQuadricMirror"
        super().__init__(**kwargs)
        self.inverse_p = inverse_p
        self.inverse_q = inverse_q
        self.theta0 = theta0

    @property
    def inverse_p(self):
        self._inverse_p = self._get_parameter("invp")
        return self._inverse_p

    @inverse_p.setter
    def inverse_p(self, value):
        self._inverse_p = self._set_parameter("invp", value)

    @property
    def inverse_q(self):
        self._inverse_q = self._get_parameter("invq")
        return self._inverse_q

    @inverse_q.setter
    def inverse_q(self, value):
        self._inverse_q = self._set_parameter("invq", value)

    @property
    def theta0(self):
        self._theta0 = self._get_parameter("theta0")
        return self._theta0

    @theta0.setter
    def theta0(self, value):
        self._theta0 = self._set_parameter("theta0", value)


class ConicCylindricalMirror(RevolutionQuadricMirror):
    """
    Class for conic cylindrical mirrors. Inherits RevolutionQuadricMirror.
    The conic surface is defined as for the RevolutionQuadricMirror class
    except it is invariant along the X axe (width).
    """

    def __init__(self, **kwargs):
        """
        Constructor for the ConicCylindricalMirror class.
        :param kwargs: See RevolutionQuadricMirror doc for parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "ConicBaseCylindricalMirror"
        else:
            kwargs["element_type"] = "ConicBaseCylindricalMirror"
        super().__init__(**kwargs)


class SphericalMirror(OpticalElement):
    """
    Class for the spherical mirrors. Inherits OpticalElement. Inherited by "SphericalFilm",
    "SphericalHoloGrating" and "SphericalPoly1DGrating"

    Defines a 2d spherical mirror of radius 1/curvature.
    """

    def __init__(self, curvature=0, **kwargs):
        """
        Constructor for the SphericalMirror class. Radius of curvature of the surface is defined as 1/curvature.

        :param curvature: inverse of the radius of the sphere in m-1
        :type curvature: float
        :param kwargs: See OpticalElement doc for other parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ("SphericalMirror", "SphericalFilm", "SphericalHoloGrating",
                                              "SphericalPoly1DGrating")
        else:
            kwargs["element_type"] = "SphericalMirror"
        super().__init__(**kwargs)
        self.curvature = curvature

    @property
    def curvature(self):
        self._curvature = self._get_parameter("curvature")
        return self._curvature

    @curvature.setter
    def curvature(self, value):
        self._curvature = self._set_parameter("curvature", value)


class CylindricalMirror(OpticalElement):
    """
    Class for the cylindrical spherical mirrors. Inherits OpticalElement. Inherited by "CylindricalFilm",
    "CylindricalHoloGrating" and "CylindricalPoly1DGrating".

    Defines a 1d spherical mirror of radius 1/curvature invariant along an axes at axis_angle from the X axes of
    the mirror.
    """

    def __init__(self, curvature=0, axis_angle=0, **kwargs):
        """
        Constructor for the CylindricalMirror class. curvature is the inverse of the base circle radius,
        cylinder axis is at an angle axis_angle from the mirror X axis.

        :param curvature: inverse of the base circle radius in m-1
        :type curvature: float
        :param axis_angle: angle between the mirror X axis and the cylinder axis of invariance : 0 = tangential cylinder
        :type axis_angle: float
        :param kwargs: See OpticalElement doc for other parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ("CylindricalMirror", "CylindricalFilm", "CylindricalHoloGrating",
                                              "CylindricalPoly1DGrating")
        else:
            kwargs["element_type"] = "CylindricalMirror"
        super().__init__(**kwargs)
        self.curvature = curvature
        self.axis_angle = axis_angle

    @property
    def curvature(self):
        self._curvature = self._get_parameter("curvature")
        return self._curvature

    @curvature.setter
    def curvature(self, value):
        self._curvature = self._set_parameter("curvature", value)

    @property
    def axis_angle(self):
        self._axis_angle = self._get_parameter("axis_angle")
        return self._axis_angle

    @axis_angle.setter
    def axis_angle(self, value):
        self._axis_angle = self._set_parameter("axis_angle", value)


class ToroidalMirror(OpticalElement):
    """
    Class for the toroidal mirrors. Inherits OpticalElement. Inherited by "ToroidalFilm",
    "ToroidalHoloGrating" and "ToroidalPoly1DGrating".

    Defines a toroidal mirror of major radius 1/major_curvature along its length and minor radius 1/minor_curvature
    along its width.
    """

    def __init__(self, minor_curvature=0, major_curvature=0, **kwargs):
        """
        Constructor for the ToroidalMirror class.

        :param minor_curvature: inverse of radius of curvature in the width of the mirror (around tangential axis)
        :type minor_curvature: float
        :param major_curvature: inverse of radius of curvature in the length of the mirror (around sagittal axis)
        :type major_curvature: float
        :param kwargs: See OpticalElement doc for other parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ("ToroidalMirror", "ToroidalFilm", "ToroidalHoloGrating",
                                              "ToroidalPoly1DGrating")
        else:
            kwargs["element_type"] = "ToroidalMirror"
        super().__init__(**kwargs)
        self.major_curvature = major_curvature
        self.minor_curvature = minor_curvature

    @property
    def minor_curvature(self):
        self._minor_curvature = self._get_parameter("minor_curvature")
        return self._minor_curvature

    @minor_curvature.setter
    def minor_curvature(self, value):
        self._minor_curvature = self._set_parameter("minor_curvature", value)

    @property
    def major_curvature(self):
        self._major_curvature = self._get_parameter("major_curvature")
        return self._major_curvature

    @major_curvature.setter
    def major_curvature(self, value):
        self._major_curvature = self._set_parameter("major_curvature", value)


class PlaneFilm(OpticalElement):
    """
    Class for plane films. Inherits OpticalElement
    """

    def __init__(self, **kwargs):
        """
        Constructor of the PlaneFilm class
        :param kwargs: See OpticalElement doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "PlaneFilm"
        else:
            kwargs["element_type"] = "PlaneFilm"
        super().__init__(**kwargs)


class SphericalFilm(SphericalMirror):
    """
    Class for spherical films. Inherits SphericalMirror
    """

    def __init__(self, **kwargs):
        """
        Constructor of the SphericalFilm class.

        :param kwargs: See SphericalMirror doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "SphericalFilm"
        else:
            kwargs["element_type"] = "SphericalFilm"
        super().__init__(**kwargs)


class CylindricalFilm(CylindricalMirror):
    """
    Class for cylindrical films. Inherits CylindricalMirror
    """

    def __init__(self, **kwargs):
        """
        Constructor of the CylindricalFilm class.

        :param kwargs: See CylindricalMirror doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "CylindricalFilm"
        else:
            kwargs["element_type"] = "CylindricalFilm"
        super().__init__(**kwargs)


class ToroidalFilm(ToroidalMirror):
    """
    Class for toroidal films. Inherits ToroidalMirror
    """

    def __init__(self, **kwargs):
        """
        Constructor of the ToroidalFilm class
        :param kwargs: See ToroidalMirror doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "ToroidalFilm"
        else:
            kwargs["element_type"] = "ToroidalFilm"
        super().__init__(**kwargs)


class Grating(OpticalElement):
    """
    Abstract class to be inherited by all grating type Optical Element
    """

    def __init__(self, line_density=1e6, order_use=1, order_align=1, **kwargs):
        super().__init__(**kwargs)
        self.order_align = order_align
        self.order_use = order_use
        self.line_density = line_density

    @property
    def order_align(self):
        self._order_align = self._get_parameter("order_align")
        return self._order_align

    @order_align.setter
    def order_align(self, value):
        self._order_align = self._set_parameter("order_align", value)

    @property
    def order_use(self):
        self._order_use = self._get_parameter("order_use")
        return self._order_use

    @order_use.setter
    def order_use(self, value):
        self._order_use = self._set_parameter("order_use", value)

    @property
    def line_density(self):
        self._line_density = self._get_parameter("lineDensity")
        return self._line_density

    @line_density.setter
    def line_density(self, value):
        self._line_density = self._set_parameter("lineDensity", value)


class PlaneHoloGrating(Grating):
    """
    Class for plane holographic grating.

    This simulates a grating which pattern was recorded using lights sources at point M and N emitting at
    recording_wavelength wavelength.
    M coordinates are defined in spherical coordinates by (azimuth_angle1, elevation_angle1 and distance1).
    N coordinates are partially defined in spherical coordinates by (azimuth_angle2 and distance2). Elevation angle of
    point N is computed from line density of the grating and the other M and N coordinates.

    Conventions for angles are :
    * angles are defined between -pi/2 and pi/2
    * waves are diverging (real source point) if elevation and distance are of the same sign, converging otherwise
    (virtual source points)
    """

    def __init__(self, azimuth_angle1=0, azimuth_angle2=0, elevation_angle1=0, inverse_distance1=np.inf,
                 inverse_distance2=np.inf, order_align=1, order_use=1, recording_wavelength=351.1e-9,
                 line_density=1e6, **kwargs):
        """
        Constructor for the PlaneHoloGrating class

        :param azimuth_angle1: coordinate in azimuth (around Z axis) of source point M in rad
        :type azimuth_angle1: float
        :param azimuth_angle2: coordinate in azimuth (around Z axis) of source point N in rad
        :type azimuth_angle2: float
        :param elevation_angle1: coordinate in elevation (in (k,Z) plane) of source point M in rad
        :type elevation_angle1: float
        :param inverse_distance1: inverse distance of the M point from the grating center in m-1
        :type inverse_distance1: float
        :param inverse_distance2: inverse distance of the N point from the grating center in m-1
        :type inverse_distance2: float
        :param order_align: order of diffraction at which the grating is aligned
        :type order_align: int
        :param order_use: order of diffraction at which the grating diffracts
        :type order_use: int
        :param recording_wavelength: wavelength of the light emitted from M and N points in m
        :type recording_wavelength: float
        :param line_density: line density at the center of the grating in m-1
        :type line_density: float
        :param kwargs: See OpticalElement doc for other parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ("PlaneHoloGrating", "SphericalHoloGrating", "CylindricalHoloGrating",
                                              "ToroidalHoloGrating")
        else:
            kwargs["element_type"] = "PlaneHoloGrating"
        super().__init__(order_align=order_align, order_use=order_use, line_density=line_density, **kwargs)
        self.azimuth_angle1 = azimuth_angle1
        self.azimuth_angle2 = azimuth_angle2
        self.elevation_angle1 = elevation_angle1
        self.inverse_distance1 = inverse_distance1
        self.inverse_distance2 = inverse_distance2
        self.recording_wavelength = recording_wavelength
        self.construction_point_1 = SphericalPoint(-1 / inverse_distance1, pi / 2 - elevation_angle1, 0)
        self.construction_point_2 = SphericalPoint(-1 / inverse_distance2,
                                                   pi / 2 - arccos(
                                                       cos(elevation_angle1) -
                                                       line_density * recording_wavelength),
                                                   0)

    @property
    def azimuth_angle1(self):
        self._azimuth_angle1 = self._get_parameter("azimuthAngle1")
        return self._azimuth_angle1

    @azimuth_angle1.setter
    def azimuth_angle1(self, value):
        self._azimuth_angle1 = self._set_parameter("azimuthAngle1", value)

    @property
    def azimuth_angle2(self):
        self._azimuth_angle2 = self._get_parameter("azimuthAngle2")
        return self._azimuth_angle2

    @azimuth_angle2.setter
    def azimuth_angle2(self, value):
        self._azimuth_angle2 = self._set_parameter("azimuthAngle2", value)

    @property
    def elevation_angle1(self):
        self._elevation_angle1 = self._get_parameter("elevationAngle1")
        return self._elevation_angle1

    @elevation_angle1.setter
    def elevation_angle1(self, value):
        self._elevation_angle1 = self._set_parameter("elevationAngle1", value)

    @property
    def inverse_distance1(self):
        self._inverse_distance1 = self._get_parameter("inverseDist1")
        return self._inverse_distance1

    @inverse_distance1.setter
    def inverse_distance1(self, value):
        self._inverse_distance1 = self._set_parameter("inverseDist1", value)

    @property
    def inverse_distance2(self):
        self._inverse_distance2 = self._get_parameter("inverseDist2")
        return self._inverse_distance2

    @inverse_distance2.setter
    def inverse_distance2(self, value):
        self._inverse_distance2 = self._set_parameter("inverseDist2", value)

    @property
    def recording_wavelength(self):
        self._recording_wavelength = self._get_parameter("recordingWavelength")
        return self._recording_wavelength

    @recording_wavelength.setter
    def recording_wavelength(self, value):
        self._recording_wavelength = self._set_parameter("recordingWavelength", value)

    def update_construction_points(self):
        """
        Method to be called for recomputation of the construction points of the hologram. In particular, to be called
        after any change applied to the hologram related attribute of the object.

        :return: None
        :rtype: Nonetype
        """
        self.construction_point_1 = SphericalPoint(1 / self.inverse_distance1, self.elevation_angle1, 0)
        self.construction_point_2 = SphericalPoint(1 / self.inverse_distance2,
                                                   arccos(
                                                       cos(self.elevation_angle1) - self.line_density * self.recording_wavelength),
                                                   0)

    def from_solemio(self, inverse_dist1, inverse_dist2, cos1, dcos, lamda=None):
        """
        Definition of the hologram using the SOLEMIO convention where:

        :param inverse_dist1: inverse of the distance between grating and point of construction 1
        :type inverse_dist1: float
        :param inverse_dist2: inverse of the distance between grating and point of construction 2
        :type inverse_dist2: float
        :param cos1: cosine of the grazing angle (elevation) from grating to point of construction 1
        :type cos1: float
        :param dcos: difference between cosines of the grazing angle to points of construction
        :type dcos: float
        :param lamda: wavelength used for recording the hologram
        :type lamda: float
        :return: None
        :rtype: Nonetype
        """
        # if dcos < 0:
        #     self.inverse_distance1 = -inverse_dist1
        #     self.inverse_distance2 = -inverse_dist2
        #     self.elevation_angle1 = arccos(cos1)
        # else:
        #     self.inverse_distance1 = -inverse_dist2
        #     self.inverse_distance2 = -inverse_dist1
        #     self.elevation_angle1 = arccos(cos1 + dcos)
        self.inverse_distance1 = -inverse_dist1
        self.inverse_distance2 = -inverse_dist2
        self.elevation_angle1 = -min(arccos(cos1), arccos(cos1 + dcos))
        if lamda:
            self.recording_wavelength = lamda
        self.line_density = abs(dcos / self.recording_wavelength)
        self.update_construction_points()
        VLS_coeffs = self.get_vls_law(100e-3, 3)
        VLS_poly = np.polynomial.polynomial.Polynomial(VLS_coeffs[::-1])
        print(f"line density computed from solemio definition : {VLS_poly} /mm")
        print("angle ", self.elevation_angle1 / degree, "deg")
        print("inv1 ", self.inverse_distance1)
        print("inv2 ", self.inverse_distance2)

    def from_horiba(self, dist1, dist2, angle1, angle2, lamda=None):
        """
        Definition of the hologram using the Horiba (Jobin-Yvon) convention where:

        :param dist1: distance between grating and point of construction 1
        :type dist1: float
        :param dist2:  distance between grating and point of construction 2
        :type dist2: float
        :param angle1: angle (spherical) between normal and axis from grating to point of construction 1
        :type angle1: float
        :param angle2: angle (spherical) between normal and axis from grating to point of construction 2
        :type angle2: float
        :param lamda: wavelength used for recording the hologram
        :type lamda: float
        :return: None
        :rtype: Nonetype
        """
        if abs(angle1) < abs(angle2):
            self.inverse_distance1 = 1 / dist1
            self.inverse_distance2 = 1 / dist2
            self.elevation_angle1 = -(angle1 + pi / 2)
        else:
            self.inverse_distance1 = 1 / dist2
            self.inverse_distance2 = 1 / dist1
            self.elevation_angle1 = -(angle2 + pi / 2)
        if lamda:
            self.recording_wavelength = lamda
        self.line_density = abs(sin(angle1) - sin(angle2)) / self.recording_wavelength
        print(f"line density at center computed from solemio definition : {self.line_density / 1000} /mm")
        print(self.elevation_angle1)
        print(self.inverse_distance1)
        print(self.inverse_distance2)
        self.update_construction_points()

    def get_horiba_coords(self):
        """
        Method for calculating the hologram parameters in the convention used by Horiba (Jobin Yvon).

        :return: (Coordinates of point 1, coordinates of point 2 , recording wavelength)
        :rtype: tuple
        """
        self.update_construction_points()
        d1, theta1, phi1 = self.construction_point_1.get_coord_horiba()
        d2, theta2, phi2 = self.construction_point_2.get_coord_horiba()

        print("distance1 = ", d1)
        print(f"angle1 =  {theta1} deg")
        print(f"phi1 = {phi1} deg")
        print("distance2 = ", d2)
        print(f"angle2 = {theta2} deg")
        print(f"phi2 = {phi2} deg")
        return self.construction_point_1.get_coord_pyoptix(), self.construction_point_2.get_coord_pyoptix(), \
            self.recording_wavelength

    def get_pyoptix_coords(self):
        """
        Method for calculating the hologram parameters in the convention used by PyOptiX.

        :return: (Coordinates of point 1, coordinates of point 2 , recording wavelength, line density at center)
        :rtype: tuple
        """
        self.update_construction_points()
        inv_d1, elevation1, phi1 = self.construction_point_1.get_coord_pyoptix()
        inv_d2, elevation2, phi2 = self.construction_point_2.get_coord_pyoptix()
        if elevation1 > elevation2:
            tmp = (inv_d1, elevation1, phi1)
            inv_d1, elevation1, phi1 = inv_d2, elevation2, phi2
            inv_d2, elevation2, phi2 = tmp
        print("inverse_distance1 = ", inv_d1)
        print("phi1 = ", phi1)
        print("elevation_angle1 = ", elevation1)
        print("inverse_distance2 = ", inv_d2)
        print("phi2 = ", phi2)
        print("tpm = ", self.get_tpmm(0))
        return self.construction_point_1.get_coord_pyoptix(), self.construction_point_2.get_coord_pyoptix(), \
            self.recording_wavelength, self.get_tpmm(0)

    def get_solemio_coords(self):
        """
        Method for calculating the hologram parameters in the convention used by SOLEMIO.

        :return: (Coordinates of point 1, coordinates of point 2 , recording wavelength)
        :rtype: tuple
        """
        self.update_construction_points()
        inv_d, cos1, phi = self.construction_point_1.get_coord_solemio()
        print("inverse_distance1 = ", inv_d)
        print("cos1 = ", cos1)
        print("phi1 = ", phi)
        inv_d, dcos, phi = self.construction_point_2.get_coord_solemio(point1=self.construction_point_1)
        print("inverse_distance2 = ", inv_d)
        print("dcos = ", dcos)
        print("phi2 = ", phi)
        return self.construction_point_1.get_coord_solemio(), \
            self.construction_point_2.get_coord_solemio(point1=self.construction_point_1), \
            self.recording_wavelength

    def get_tpmm(self, point_m):
        """
        Method for calculating the local line density of an hologram at a given point which coordinates are given
        in point_m.

        :param point_m: point where the local  line density is computed.
        :type point_m: pyoptix.Point
        :return: line density in m^-1
        :rtype: float
        """
        self.update_construction_points()
        if not isinstance(point_m, Point):
            point_m = Point(point_m, 0, 0)
        p1_apparent = self.construction_point_1 - point_m
        p2_apparent = self.construction_point_2 - point_m
        _, theta_apparent_1, _ = p1_apparent.get_spherical_coord()
        _, theta_apparent_2, _ = p2_apparent.get_spherical_coord()
        return np.abs(sin(theta_apparent_2) - sin(theta_apparent_1)) / self.recording_wavelength

    def get_vls_law(self, x_span, order, return_tpm=False, sample_number=20):
        """
        Method for computing the Varied Line Space law of the hologram along its  central axis (azimuth=0).
        Either returns polynomial coefficient (up to degree 'order') or both coefficients and list of local
        line densities at "sample_number" points evenly space over a range of "x_span" meters.

        :param x_span: total range over which to compute the VLS law
        :type x_span: float
        :param order: degree of the polynomial to which the VLS law is fitted
        :type order: int
        :param return_tpm: if True, local list of line densities are returned in a tuple after fit coefficients
        :type return_tpm: bool
        :param sample_number: Number of points along the central axis where local line densities are computed
        :type sample_number: int
        :return: tuple of (array of fit coeffs, list of line densities) or array of fit coeffs)
        :rtype: array or tuple
        """
        self.update_construction_points()
        tpm = []  # in m-1
        for x in np.linspace(-x_span / 2, x_span / 2, sample_number):
            try:
                tpm.append(self.get_tpmm(Point(x, 0, 0)))
            except RuntimeError:
                tpm.append(np.nan)

        gratinfo = GratingPatternInformation()
        get_hologram_pattern(self._element_id, gratinfo, x_span, 10e-3)
        print(f"C hologram info for {self.name}: ")
        print(f"\t groove/m = {gratinfo.axial_line_density.a0} + {gratinfo.axial_line_density.a1} * x + "
              f"{gratinfo.axial_line_density.a2} * x^2 + {gratinfo.axial_line_density.a3} * x^3")
        print(f"\t line radius curvature = {gratinfo.line_curvature}")
        print(f"\t line tilt = {gratinfo.line_tilt}")
        tpm = np.array(tpm, dtype=np.float64)
        if not return_tpm:
            return np.polyfit(np.linspace(-x_span / 2, x_span / 2, sample_number, dtype=np.float64), tpm, order) * \
                np.array([milli ** (n + 1) for n in range(order + 1)][::-1], dtype=np.float64)
        else:
            return (np.polyfit(np.linspace(-x_span / 2, x_span / 2, sample_number, dtype=np.float64), tpm, order) *
                    np.array([milli ** (n + 1) for n in range(order + 1)][::-1], dtype=np.float64), tpm)

    def show_vls_law(self, x_span, order, sample_number=20):
        """
        Displays the local line densities computed by pyoptix.PlaneHoloGrating.get_vls_law as well as the fit of
        those diplaying in legend the coefficients of the polynomial fit.

        :param x_span: total range over which to compute the VLS law
        :type x_span: float
        :param order: degree of the polynomial to which the VLS law is fitted
        :type order: int
        :param sample_number: Number of points along the central axis where local line densities are computed
        :type sample_number: int
        :return: None
        :rtype: Nonetype
        """
        coeffs, tpm = self.get_vls_law(x_span, order, return_tpm=True, sample_number=sample_number)
        coeffs /= np.array([milli ** (n + 1) for n in range(order + 1)][::-1])
        fit = np.poly1d(coeffs)
        plt = figure(width=1000, height=300)
        plt.xaxis.axis_label = "Distance from center (m)"
        plt.yaxis.axis_label = "Number of groove per meter"
        plt.line(np.linspace(-x_span / 2, x_span / 2, sample_number), np.array(tpm), legend_label="local #groove/m")
        plt.line(np.linspace(-x_span / 2, x_span / 2, sample_number),
                 fit(np.linspace(-x_span / 2, x_span / 2, sample_number)),
                 legend_label=f"fit (mm-1) : {coeffs * np.array([milli ** (n + 1) for n in range(order + 1)][::-1])}",
                 color="red")
        show(plt)


class SphericalHoloGrating(SphericalMirror, PlaneHoloGrating):
    """
    Class for spherical holographic gratings. Inherits SphericalMirror and PlaneHoloGrating
    """

    def __init__(self, **kwargs):
        """
        Constructor of the SphericalHoloGrating class.

        :param kwargs: See SphericalMirror and PlaneHoloGrating doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "SphericalHoloGrating"
        else:
            kwargs["element_type"] = "SphericalHoloGrating"
        super().__init__(**kwargs)


class CylindricalHoloGrating(CylindricalMirror, PlaneHoloGrating):
    """
    Class for cylindrical holographic gratings. Inherits CylindricalMirror and PlaneHoloGrating
    """

    def __init__(self, **kwargs):
        """
        Constructor of the CylindricalHoloGrating class.

        :param kwargs: See CylindricalMirror and PlaneHoloGrating doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "CylindricalHoloGrating"
        else:
            kwargs["element_type"] = "CylindricalHoloGrating"
        super().__init__(**kwargs)


class ToroidalHoloGrating(ToroidalMirror, PlaneHoloGrating):
    """
    Class for toroidal holographic gratings. Inherits ToroidalMirror and PlaneHoloGrating
    """

    def __init__(self, **kwargs):
        """
        Constructor of the ToroidalHoloGrating class.

        :param kwargs: See ToroidalMirror and PlaneHoloGrating doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "ToroidalHoloGrating"
        else:
            kwargs["element_type"] = "ToroidalHoloGrating"
        super().__init__(**kwargs)


class PlaneGrating(PlaneHoloGrating):
    """
    Class for plane gratings. Inherits PlaneHoloGrating, inverse_distances are 0 so there is no variation of
    line density
    """

    def __init__(self, **kwargs):
        """
        Constructor of the PlaneGrating class.

        :param kwargs: See PlaneHoloGrating doc for the parameters
        """
        kwargs["inverse_distance1"] = np.inf
        kwargs["inverse_distance2"] = np.inf
        kwargs["elevation_angle1"] = np.arccos(0.93963)
        super().__init__(**kwargs)
        self.inverse_distance1 = np.inf
        self.inverse_distance2 = np.inf
        self.elevation_angle1 = np.arccos(0.93963)

    @property
    def inverse_distance1(self):
        self._inverse_distance1 = super().inverse_distance1
        return self._inverse_distance1

    @inverse_distance1.setter
    def inverse_distance1(self, value):
        if value != np.inf:
            raise Exception("Property 'inverse_distance1' is not settable")

    @property
    def inverse_distance2(self):
        self._inverse_distance2 = super().inverse_distance2
        return self._inverse_distance2

    @inverse_distance2.setter
    def inverse_distance2(self, value):
        if value != np.inf:
            raise Exception("Property 'inverse_distance2' is not settable")

    @property
    def elevation_angle1(self):
        self._elevation_angle1 = super().elevation_angle1
        return self._elevation_angle1

    @elevation_angle1.setter
    def elevation_angle1(self, value):
        if value != np.arccos(0.93963):
            raise Exception("Property 'elevation_angle1' is not settable")


class PlanePoly1DGrating(Grating):
    """
    Class for plane polynomial gratings, which line density varies as a function if only one coordinate in a polynomial
    way. The center of the grating is assumed to be at the 0 value of the coordinate,leading to the parameter
    `line_density` defining the line density at the grating center.

    Inherits OpticalElement. Inherited by SphericalPoly1DGrating, CylindricalPoly1DGrating and
    ToroidalPoly1DGrating.
    """

    def __init__(self, polynomial_degree=0, line_density=1e6, line_density_coeffs=None, order_align=1, order_use=1,
                 **kwargs):
        """
        Constructor method for the class PlanePoly1DGrating.

        :param polynomial_degree: degree of the polynomial law for varying line density. 0 makes an evenly spaced
            classical grating
        :type polynomial_degree: int
        :param line_density: line density at the center of the grating (in m-1)
        :type line_density: float
        :param line_density_coeffs: coefficients of the polynomial law in m-2, m-3, etc. List length must match
            polynomial_degree
        :type line_density_coeffs: list of float
        :param kwargs: See OpticalElement doc for additional parameters
        """
        if line_density_coeffs is None:
            line_density_coeffs = []
        if "element_type" in kwargs:
            assert kwargs["element_type"] in ("PlanePoly1DGrating", "SphericalPoly1DGrating",
                                              "CylindricalPoly1DGrating",
                                              "ToroidalPoly1DGrating")
        else:
            kwargs["element_type"] = "PlanePoly1DGrating"
        super().__init__(line_density=line_density, order_use=order_use, order_align=order_align, **kwargs)
        self.degree = polynomial_degree
        self.line_density_coeffs = line_density_coeffs

    @property
    def degree(self):
        self._degree = int(self._get_parameter("degree"))
        return self._degree

    @degree.setter
    def degree(self, value):
        self._degree = int(self._set_parameter("degree", value))

    @property
    def line_density(self):
        self._line_density = self._get_parameter("lineDensity")
        return self._line_density

    @line_density.setter
    def line_density(self, value):
        self._line_density = self._set_parameter("lineDensity", value)

    @property
    def line_density_coeffs(self):
        self._line_density_coeffs = []
        for i in range(1, self.degree + 1):
            self._line_density_coeffs.append(self._get_parameter(f"lineDensityCoeff_{i}"))
        return self._line_density_coeffs

    @line_density_coeffs.setter
    def line_density_coeffs(self, value):
        self._line_density_coeffs = []
        assert isinstance(value, list)
        for i in range(1, self.degree + 1):
            self._line_density_coeffs.append(self._set_parameter(f"lineDensityCoeff_{i}", value[i - 1]))


class SphericalPoly1DGrating(SphericalMirror, PlanePoly1DGrating):
    """
    Class for spherical polynomial gratings, which line density varies as a function if only one coordinate in a
    polynomial way. The center of the grating is assumed to be at the 0 value of the coordinate,leading to the parameter
    `line_density` defining the line density at the grating center.

    Inherits SphericalMirror and PlanePoly1DGrating.
    """

    def __init__(self, **kwargs):
        """
        Constructor of the SphericalPoly1DGrating class.

        :param kwargs: See SphericalMirror and PlanePoly1DGrating doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "SphericalPoly1DGrating"
        else:
            kwargs["element_type"] = "SphericalPoly1DGrating"
        super().__init__(**kwargs)


class CylindricalPoly1DGrating(CylindricalMirror, PlanePoly1DGrating):
    """
    Class for cylindrical polynomial gratings, which line density varies as a function if only one coordinate in a
    polynomial way. The center of the grating is assumed to be at the 0 value of the coordinate,leading to the parameter
    `line_density` defining the line density at the grating center.

    Inherits CylindricalMirror and PlanePoly1DGrating.
    """

    def __init__(self, **kwargs):
        """
        Constructor of the CylindricalPoly1DGrating class
        :param kwargs: See CylindricalMirror and PlanePoly1DGrating doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "CylindricalPoly1DGrating"
        else:
            kwargs["element_type"] = "CylindricalPoly1DGrating"
        super().__init__(**kwargs)


class ToroidalPoly1DGrating(ToroidalMirror, PlanePoly1DGrating):
    """
    Class for toroidal polynomial gratings, which line density varies as a function if only one coordinate in a
    polynomial way. The center of the grating is assumed to be at the 0 value of the coordinate,leading to the parameter
    `line_density` defining the line density at the grating center.

    Inherits ToroidalMirror and PlanePoly1DGrating.
    """

    def __init__(self, **kwargs):
        """
        Constructor of the ToroidalPoly1DGrating class.

        :param kwargs: See ToroidalMirror and PlanePoly1DGrating doc for the parameters
        """
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "ToroidalPoly1DGrating"
        else:
            kwargs["element_type"] = "ToroidalPoly1DGrating"
        super().__init__(**kwargs)


def parse_xml(filename):
    tree = etree.parse(filename)
    beamline = Beamline()
    for user in tree.xpath("/system/element"):
        new_element = OpticalElement(name=user.get("name"), next_element=user.get("next"),
                                     previous=user.get("previous"))
        beamline.add_element(new_element)
    beamline.chain()

    for chain in beamline.chains:
        desc = ""
        for element in chain:
            desc += element.name + " -> "
        print(desc)


def load_beamline(filename, glob_dict, verbose=1):
    with open(filename, "rb") as filein:
        data = pickle.load(filein)
    glob_dict[data["beamline"]] = Beamline(name=data["beamline"])
    if verbose:
        print("Retrieving beamline", data["beamline"])
    active_chain = []
    for item in data["config"]:
        if item['oe_name'] not in glob_dict:
            glob_dict[item['oe_name']] = item['oe_class'](name=item['oe_name'])
        oe = glob_dict[item['oe_name']]
        active_chain.append(oe)
        for param_name, param in item['oe_params'].items():
            oe._set_parameter(param_name, param)
        if verbose:
            print("Retrieving and setting element", item['oe_name'])
    glob_dict[data["beamline"]].chains[data["active_chain"]] = active_chain
    glob_dict[data["beamline"]].active_chain = data["active_chain"]


def save_beamline(beamline, active_chain_name, filename):
    data = {"beamline": beamline.name, "config": [oe.dump_properties(verbose=0) for oe in beamline.active_chain],
            "active_chain": active_chain_name}
    with open(filename, "wb") as fileout:
        pickle.dump(data, fileout)
