# coding: utf-8
import numpy as np
from ctypes import Structure, c_int, c_double, c_uint32, c_int32, POINTER, c_void_p, cast, c_int64, create_string_buffer
from ctypes.wintypes import BYTE, INT, HMODULE, LPCSTR, HANDLE, DOUBLE
from exposed_functions import *
from scipy.constants import degree
from lxml import etree


optix_dictionnary = {
                     "DX": "d_x",
                     "DY": "d_y",
                     "DZ": "d_z",
                     "Dphi": "d_phi",
                     "Dpsi": "d_psi",
                     "Dtheta": "d_theta",
                     "distance": "distance_from_previous",
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


class Beamline(object):
    def __init__(self):
        super().__init__()
        self.elements = []
        self.chains = []

    def align(self, from_element):
        pass

    def radiate(self, from_element):
        pass

    def generate(self, n_rays=None):
        pass

    def chain(self):
        for element in self.elements:
            print(element)
            if element.next is None:
                new_chain = True
                for chain in self.chains:
                    if chain[0] == element:
                        new_chain = False
                if new_chain:
                    self.chains.append([element])
                    print("new_chain")

        def get_previous(an_element):
            for element_i in self.elements:
                if element_i.next == an_element.name:
                    return element_i

        for chain in self.chains:
            previous = 1
            while previous is not None:
                previous = get_previous(chain[0])
                if previous:
                    chain.insert(0, previous)

    def add_element(self, new_element):
        if new_element not in self.elements:
            self.elements.append(new_element)
        else:
            print("Warning : element already in beamline")
        for element in self.elements:
            if element.name == new_element.next:
                element.previous = new_element.name


class OpticalElement(object):
    def __init__(self, name="", phi=0, psi=0, theta=0, d_phi=0, d_psi=0, d_theta=0, x=0, y=0, z=0,
                 d_x=0, d_y=0, d_z=0, next=None, previous=None, distance_from_previous=0, element_id=None,
                 element_type=""):
        super().__init__()
        self._name = name
        self._phi = phi
        self._psi = psi
        self._theta = theta
        self._d_phi = d_phi
        self._d_psi = d_psi
        self._d_theta = d_theta
        self._x = x
        self._y = y
        self._z = z
        self._d_x = d_x
        self._d_y = d_y
        self._d_z = d_z
        self._next = next
        self._previous = previous
        self._distance_from_previous = distance_from_previous
        self._element_id = None
        self._element_type = element_type
        if element_id is not None:
            self.from_element_id(element_id)
        elif element_type != "" and name != "":
            self._element_id = create_element(element_type, name)
            self._name = name

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
    def x(self):
        self._x = self._get_parameter("x")
        return self._x

    @x.setter
    def x(self, value):
        param = Parameter()
        get_parameter(self._element_id, "x", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "x", param)
        self._x = self._get_parameter("x")

    @property
    def y(self):
        self._y = self._get_parameter("y")
        return self._y

    @y.setter
    def y(self, value):
        param = Parameter()
        get_parameter(self._element_id, "y", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "y", param)
        self._y = self._get_parameter("y")

    @property
    def z(self):
        self._z = self._get_parameter("z")
        return self._z

    @z.setter
    def z(self, value):
        param = Parameter()
        get_parameter(self._element_id, "z", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "z", param)
        self._z = self._get_parameter("z")

    @property
    def d_x(self):
        self._d_x = self._get_parameter("Dx")
        return self._d_x

    @d_x.setter
    def d_x(self, value):
        param = Parameter()
        get_parameter(self._element_id, "Dx", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "Dx", param)
        self._d_x = self._get_parameter("Dx")

    @property
    def d_y(self):
        self._d_y = self._get_parameter("Dy")
        return self._d_y

    @d_y.setter
    def d_y(self, value):
        param = Parameter()
        get_parameter(self._element_id, "Dy", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "Dy", param)
        self._d_y = self._get_parameter("Dy")

    @property
    def d_z(self):
        self._d_z = self._get_parameter("Dz")
        return self._d_z

    @d_z.setter
    def d_z(self, value):
        param = Parameter()
        get_parameter(self._element_id, "Dz", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "Dz", param)
        self._d_z = self._get_parameter("Dz")

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
    def previous(self, value):
        param = Parameter()
        get_parameter(self._element_id, "previous", param)
        param.value = DOUBLE(value)
        set_parameter(self._element_id, "previous", param)
        self._previous = self._get_parameter("previous")

    @property
    def next(self):
        return self._next

    @next.setter
    def next(self, next_oe):
        self._next = next_oe
        chain_element_by_id(self._element_id, next_oe.element_id)

    def __repr__(self):
        description = f"Element {self._name} of class {self.__class__}"
        description += f"\n\t at {self._distance_from_previous} m from {self._previous}"
        description += f"\n\t pointing to {self._next}"
        description += f"\n\t oriented in pitch at {self._theta/degree} deg (deviation {180-2*self._theta/degree} deg)"
        description += f"\n\t oriented in roll at {self._phi/degree} deg"
        description += f"\n\t oriented in yaw at {self._psi/degree} deg"
        return description

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
                self.__dict__["_"+optix_dictionnary[param_name.value.decode()]] = param.value
            else:
                self.__dict__["_"+param_name.value.decode()] = param.value
            enumerate_parameters(element_id, hparam, param_name, param, confirm=False)
        if print_all:
            print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
                        f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
        if param_name.value.decode() in optix_dictionnary.keys():
            self.__dict__["_"+optix_dictionnary[param_name.value.decode()]] = param.value
        else:
            self.__dict__["_"+param_name.value.decode()] = param.value
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
        self._sigma_x = sigma_x
        self._sigma_y = sigma_y
        self._sigma_x_div = sigma_x_div
        self._sigma_y_div = sigma_y_div
        self._nrays = nrays

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
