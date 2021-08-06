# coding: utf-8
import numpy as np
from ctypes import Structure, c_int, c_double, c_uint32, c_int32, POINTER, c_void_p, cast, c_int64, create_string_buffer
from ctypes.wintypes import BYTE, INT, HMODULE, LPCSTR, HANDLE, DOUBLE
from .exposed_functions import *
from scipy.constants import degree
from lxml import etree
import pandas as pd
from .ui_objects import show, plot_spd


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
    "NRays": "nrays",
    "azimuthAngle1": "azimuth_angle1",
    "azimuthAngle2": "azimuth_angle2",
    "elevationAngle1": "elevation_angle1",
    "inverseDist1": "inverse_distance1",
    "inverseDist2": "inverse_distance2",
    "recordingWavelength": "recording_wavelength",
    "lineDensity": "line_density",
    "invp": "inverse_p",
    "invq": "inverse_q",
}


class Bounds(Structure):
    _fields_ = [("min", c_double),
                ("max", c_double)]


class Parameter(Structure):
    _fields_ = [("value", c_double),
                ("bounds", Bounds),
                ("multiplier", c_double),
                ("type", c_int32),
                ("group", c_int32),
                ("flags", c_uint32),
                ]


class RecordingMode(object):
    recording_output = c_int32(2)
    recording_input = c_int32(1)
    recording_none = c_int32(0)


class Diagram(Structure):
    _fields_ = [("dim", c_int),
                ("reserved", c_int),
                ("count", c_int),
                ("lost", c_int),
                ("min", POINTER(c_double)),
                ("max", POINTER(c_double)),
                ("mean", POINTER(c_double)),
                ("sigma", POINTER(c_double)),
                ("spots", POINTER(c_double)),
                ]

    def __init__(self, nreserved=0, ndim=5):
        super().__init__()
        self.dim = ndim
        self.reserved = nreserved
        p_min = (c_double * ndim)()
        self.min = cast(p_min, POINTER(c_double))
        p_max = (c_double * ndim)()
        self.max = cast(p_max, POINTER(c_double))
        p_mean = (c_double * ndim)()
        self.mean = cast(p_mean, POINTER(c_double))
        p_sigma = (c_double * ndim)()
        self.sigma = cast(p_sigma, POINTER(c_double))
        p_spots = (c_double * ndim * nreserved)()
        self.spots = cast(p_spots, POINTER(c_double))


class ChainList(list):
    def __repr__(self):
        ret_str = ''
        for oe in self:
            ret_str += f"{oe.name} -> "
        ret_str += "\n"
        return ret_str[:-4]


class BeamlineChainDict(dict):
    def __setitem__(self, key, value):
        assert isinstance(value, list)
        super(BeamlineChainDict, self).__setitem__(key, ChainList(value))

    def __repr__(self):
        ret_str = ""
        for key in self.keys():
            ret_str += f"Chaîne {key}:\n\t"
            chain = self.__getitem__(key)
            ret_str += chain.__repr__()
        return ret_str


class Beamline(object):
    def __init__(self):
        super().__init__()
        self._elements = []
        self._chains = BeamlineChainDict({})
        self._active_chain = None

    @property
    def chains(self):
        return self._chains

    @property
    def active_chain(self):
        return self._active_chain

    @active_chain.setter
    def active_chain(self, chain_name):
        assert chain_name in self._chains.keys()
        self._active_chain = self.chains[chain_name]
        for i, oe in enumerate(self._active_chain):
            try:
                oe.next = self._active_chain[i + 1]
            except IndexError:
                pass
        ret_str = f"Chaîne {chain_name}:\n\t"
        ret_str += self._active_chain.__repr__()
        print(ret_str)

    def align(self, lambda_align, from_element=None):
        if from_element is not None:
            return align(from_element.element_id, lambda_align)
        else:
            return align(self.active_chain[0].element_id, lambda_align)

    def clear_impacts(self):
        return clear_impacts(self.active_chain[0].element_id)

    def radiate(self, from_element=None):
        if from_element is not None:
            return radiate(from_element.element_id)
        else:
            return radiate(self.active_chain[0].element_id)

    def generate(self, lambda_radiate):
        return generate(self.active_chain[0].element_id, lambda_radiate)

    def _add_element(self, new_element):
        if new_element not in self._elements:
            self._elements.append(new_element)


class OpticalElement(object):
    def __init__(self, name="", phi=0, psi=0, theta=0, d_phi=0, d_psi=0, d_theta=0, x=0, y=0, z=0,
                 d_x=0, d_y=0, d_z=0, next=None, previous=None, distance_from_previous=0, element_id=None,
                 element_type=""):
        super().__init__()
        self._element_id = None
        self._element_type = element_type
        if element_id is not None:
            self.from_element_id(element_id)
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
        self.next = next
        self.previous = previous
        self.distance_from_previous = distance_from_previous

    def _get_parameter(self, param_name):
        param = Parameter()
        get_parameter(self._element_id, param_name, param)
        return param.value

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
        param = Parameter()
        get_parameter(self._element_id, "phi", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "phi", param)
        self._phi = self._get_parameter("phi")

    @property
    def psi(self):
        self._psi = self._get_parameter("psi")
        return self._psi

    @psi.setter
    def psi(self, value):
        param = Parameter()
        get_parameter(self._element_id, "psi", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "psi", param)
        self._psi = self._get_parameter("psi")

    @property
    def theta(self):
        self._theta = self._get_parameter("theta")
        return self._theta

    @theta.setter
    def theta(self, value):
        param = Parameter()
        get_parameter(self._element_id, "theta", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "theta", param)
        self._theta = self._get_parameter("theta")

    @property
    def d_phi(self):
        self._d_phi = self._get_parameter("Dphi")
        return self._d_phi

    @d_phi.setter
    def d_phi(self, value):
        param = Parameter()
        get_parameter(self._element_id, "Dphi", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "Dphi", param)
        self._d_phi = self._get_parameter("Dphi")

    @property
    def d_psi(self):
        self._d_psi = self._get_parameter("Dpsi")
        return self._d_psi

    @d_psi.setter
    def d_psi(self, value):
        param = Parameter()
        get_parameter(self._element_id, "Dpsi", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "Dpsi", param)
        self._d_psi = self._get_parameter("Dpsi")

    @property
    def d_theta(self):
        self._d_theta = self._get_parameter("Dtheta")
        return self._d_theta

    @d_theta.setter
    def d_theta(self, value):
        param = Parameter()
        get_parameter(self._element_id, "Dtheta", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "Dtheta", param)
        self._d_theta = self._get_parameter("Dtheta")

    @property
    def distance_from_previous(self):
        self._distance_from_previous = self._get_parameter("distance")
        return self._distance_from_previous

    @distance_from_previous.setter
    def distance_from_previous(self, value):
        param = Parameter()
        get_parameter(self._element_id, "distance", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "distance", param)
        self._distance_from_previous = self._get_parameter("distance")

    @property
    def previous(self):
        self._previous = self._get_parameter("previous")
        return self._previous

    @previous.setter
    def previous(self, previous_oe):
        if previous_oe is not None:
            param = Parameter()
            get_parameter(self._element_id, "previous", param)
            param.value = DOUBLE(previous_oe.element_id)
            set_parameter(self._element_id, "previous", param)
            self._previous = self._get_parameter("previous")
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
        description += f"\n\t at {self._distance_from_previous} m from {self._previous}"
        description += f"\n\t pointing to {self._next.name}"
        description += f"\n\t oriented in pitch at {self._theta / degree} deg " \
                       f"(deviation {180 - 2 * self._theta / degree} deg)"
        description += f"\n\t oriented in roll at {self._phi / degree} deg"
        description += f"\n\t oriented in yaw at {self._psi / degree} deg\n"
        return description

    def show_diagram(self, nrays, distance_from_oe=0, light_xy=False, map_xy=False,
                     light_xxp=False, light_yyp=False, show_first_rays=False):
        assert self.recording_mode != RecordingMode.recording_none
        diagram = Diagram(ndim=5, nreserved=int(nrays))
        get_spot_diagram(self.element_id, diagram, distance=distance_from_oe)
        spots = pd.DataFrame(np.ctypeslib.as_array(diagram.spots, shape=(diagram.reserved, diagram.dim)),
                             columns=("X", "Y", "dX", "dY", "Lambda"))
        if show_first_rays:
            print(spots.head())
        figs = [plot_spd(spots, x_key="X", y_key="Y", light_plot=light_xy, show_map=map_xy),
                plot_spd(spots, x_key="X", y_key="dX", light_plot=light_xxp),
                plot_spd(spots, x_key="Y", y_key="dY", light_plot=light_yyp)]

        for fig in figs:
            show(fig)

    def from_element_id(self, element_id, print_all=False):
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
            if param_name.value.decode() in optix_dictionnary.keys():
                self.__dict__["_" + optix_dictionnary[param_name.value.decode()]] = param.value
            else:
                self.__dict__["_" + param_name.value.decode()] = param.value
            enumerate_parameters(element_id, hparam, param_name, param, confirm=False)
        if print_all:
            print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
                        f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
        if param_name.value.decode() in optix_dictionnary.keys():
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
        set_recording(self._element_id,
                      {"output": RecordingMode.recording_output,
                       "input": RecordingMode.recording_input,
                       "not_recording": RecordingMode.recording_none}[recording_mode])


class Source(OpticalElement):
    def __init__(self, sigma_x=0, sigma_y=0, sigma_x_div=0, sigma_y_div=0, nrays=0, **kwargs):
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
        param = Parameter()
        get_parameter(self._element_id, "nRays", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "nRays", param)
        self._nrays = self._get_parameter("nRays")

    @property
    def sigma_x(self):
        self._sigma_x = self._get_parameter("sigmaX")
        return self._sigma_x

    @sigma_x.setter
    def sigma_x(self, value):
        param = Parameter()
        get_parameter(self._element_id, "sigmaX", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "sigmaX", param)
        self._sigma_x = self._get_parameter("sigmaX")

    @property
    def sigma_y(self):
        self._sigma_y = self._get_parameter("sigmaY")
        return self._sigma_y

    @sigma_y.setter
    def sigma_y(self, value):
        param = Parameter()
        get_parameter(self._element_id, "sigmaY", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "sigmaY", param)
        self._sigma_y = self._get_parameter("sigmaY")

    @property
    def sigma_x_div(self):
        self._sigma_x_div = self._get_parameter("sigmaXdiv")
        return self._sigma_x_div

    @sigma_x_div.setter
    def sigma_x_div(self, value):
        param = Parameter()
        get_parameter(self._element_id, "sigmaXdiv", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "sigmaXdiv", param)
        self._sigma_x_div = self._get_parameter("sigmaXdiv")

    @property
    def sigma_y_div(self):
        self._sigma_y_div = self._get_parameter("sigmaYdiv")
        return self._sigma_y_div

    @sigma_y_div.setter
    def sigma_y_div(self, value):
        param = Parameter()
        get_parameter(self._element_id, "sigmaYdiv", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "sigmaYdiv", param)
        self._sigma_y_div = self._get_parameter("sigmaYdiv")


class ConicCylindricalMirror(OpticalElement):
    """
    p is defined at coordinate (pcos(theta), -psin(theta))
    q is defined at coordinate (qcos(theta), qsin(theta))
    Due to this convention, the cylinder base (directrix) is
    •an ellipse if p−1*q−1<0
    •an hyperbola if p−1*q−1<0
    •a parabola if either p−1=0 or q−1=0
    •Warning : p−1=q−1 is forbidden and will result as an error at any time.

    """

    def __init__(self, inverse_p=0, inverse_q=0.1, theta0=0, **kwargs):
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "ConicBaseCylindricalMirror"
        else:
            kwargs["element_type"] = "ConicBaseCylindricalMirror"
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
        param = Parameter()
        get_parameter(self._element_id, "invp", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "invp", param)
        self._inverse_p = self._get_parameter("invp")

    @property
    def inverse_q(self):
        self._inverse_q = self._get_parameter("invq")
        return self._inverse_q

    @inverse_q.setter
    def inverse_q(self, value):
        param = Parameter()
        get_parameter(self._element_id, "invq", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "invq", param)
        self._inverse_q = self._get_parameter("invq")

    @property
    def theta0(self):
        self._theta0 = self._get_parameter("theta0")
        return self._theta0

    @theta0.setter
    def theta0(self, value):
        param = Parameter()
        get_parameter(self._element_id, "theta0", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "theta0", param)
        self._theta0 = self._get_parameter("theta0")


class CylindricalMirror(OpticalElement):
    def __init__(self, curvature=0, axis_angle=0, **kwargs):
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "CylindricalMirror"
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
        param = Parameter()
        get_parameter(self._element_id, "curvature", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "curvature", param)
        self._curvature = self._get_parameter("curvature")

    @property
    def axis_angle(self):
        self._axis_angle = self._get_parameter("axis_angle")
        return self._axis_angle

    @axis_angle.setter
    def axis_angle(self, value):
        param = Parameter()
        get_parameter(self._element_id, "axis_angle", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "axis_angle", param)
        self._axis_angle = self._get_parameter("axis_angle")


class ToroidalMirror(OpticalElement):
    def __init__(self, minor_curvature=0, major_curvature=0, **kwargs):
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "ToroidalMirror"
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
        param = Parameter()
        get_parameter(self._element_id, "minor_curvature", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "minor_curvature", param)
        self._minor_curvature = self._get_parameter("minor_curvature")

    @property
    def major_curvature(self):
        self._major_curvature = self._get_parameter("major_curvature")
        return self._major_curvature

    @major_curvature.setter
    def major_curvature(self, value):
        param = Parameter()
        get_parameter(self._element_id, "major_curvature", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "major_curvature", param)
        self._major_curvature = self._get_parameter("major_curvature")


class PlaneHoloGrating(OpticalElement):
    def __init__(self, azimuth_angle1=0, azimuth_angle2=0, elevation_angle1=0, inverse_distance1=np.inf,
                 inverse_distance2=np.inf, order_align=1, order_use=1, recording_wavelength=351.1e-9,
                 line_density=1e6, **kwargs):
        if "element_type" in kwargs:
            assert kwargs["element_type"] == "PlaneHoloGrating"
        else:
            kwargs["element_type"] = "PlaneHoloGrating"
        super().__init__(**kwargs)
        self.azimuth_angle1 = azimuth_angle1
        self.azimuth_angle2 = azimuth_angle2
        self.elevation_angle1 = elevation_angle1
        self.inverse_distance1 = inverse_distance1
        self.inverse_distance2 = inverse_distance2
        self.order_align = order_align
        self.order_use = order_use
        self.recording_wavelength = recording_wavelength
        self.line_density = line_density

    @property
    def azimuth_angle1(self):
        self._azimuth_angle1 = self._get_parameter("azimuthAngle1")
        return self._azimuth_angle1

    @azimuth_angle1.setter
    def azimuth_angle1(self, value):
        param = Parameter()
        get_parameter(self._element_id, "azimuthAngle1", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "azimuthAngle1", param)
        self._azimuth_angle1 = self._get_parameter("azimuthAngle1")

    @property
    def azimuth_angle2(self):
        self._azimuth_angle2 = self._get_parameter("azimuthAngle2")
        return self._azimuth_angle2

    @azimuth_angle2.setter
    def azimuth_angle2(self, value):
        param = Parameter()
        get_parameter(self._element_id, "azimuthAngle2", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "azimuthAngle2", param)
        self._azimuth_angle2 = self._get_parameter("azimuthAngle2")

    @property
    def elevation_angle1(self):
        self._elevation_angle1 = self._get_parameter("elevationAngle1")
        return self._elevation_angle1

    @elevation_angle1.setter
    def elevation_angle1(self, value):
        param = Parameter()
        get_parameter(self._element_id, "elevationAngle1", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "elevationAngle1", param)
        self._elevation_angle1 = self._get_parameter("elevationAngle1")

    @property
    def inverse_distance1(self):
        self._inverse_distance1 = self._get_parameter("inverseDist1")
        return self._inverse_distance1

    @inverse_distance1.setter
    def inverse_distance1(self, value):
        param = Parameter()
        get_parameter(self._element_id, "inverseDist1", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "inverseDist1", param)
        self._inverse_distance1 = self._get_parameter("inverseDist1")

    @property
    def inverse_distance2(self):
        self._inverse_distance2 = self._get_parameter("inverseDist2")
        return self._inverse_distance2

    @inverse_distance2.setter
    def inverse_distance2(self, value):
        param = Parameter()
        get_parameter(self._element_id, "inverseDist2", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "inverseDist2", param)
        self._inverse_distance2 = self._get_parameter("inverseDist2")

    @property
    def order_align(self):
        self._order_align = self._get_parameter("order_align")
        return self._order_align

    @order_align.setter
    def order_align(self, value):
        param = Parameter()
        get_parameter(self._element_id, "order_align", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "order_align", param)
        self._order_align = self._get_parameter("order_align")

    @property
    def order_use(self):
        self._order_use = self._get_parameter("order_use")
        return self._order_use

    @order_use.setter
    def order_use(self, value):
        param = Parameter()
        get_parameter(self._element_id, "order_use", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "order_use", param)
        self._order_use = self._get_parameter("order_use")

    @property
    def recording_wavelength(self):
        self._recording_wavelength = self._get_parameter("recordingWavelength")
        return self._recording_wavelength

    @recording_wavelength.setter
    def recording_wavelength(self, value):
        param = Parameter()
        get_parameter(self._element_id, "recordingWavelength", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "recordingWavelength", param)
        self._recording_wavelength = self._get_parameter("recordingWavelength")

    @property
    def line_density(self):
        self._line_density = self._get_parameter("lineDensity")
        return self._line_density

    @line_density.setter
    def line_density(self, value):
        param = Parameter()
        get_parameter(self._element_id, "lineDensity", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "lineDensity", param)
        self._line_density = self._get_parameter("lineDensity")


def parse_xml(filename):
    tree = etree.parse(filename)
    beamline = Beamline()
    for user in tree.xpath("/system/element"):
        new_element = OpticalElement(name=user.get("name"), next=user.get("next"), previous=user.get("previous"))
        beamline.add_element(new_element)
    beamline.chain()

    for chain in beamline.chains:
        desc = ""
        for element in chain:
            desc += element.name + " -> "
        print(desc)
