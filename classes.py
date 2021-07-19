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

    def align(self):
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
            for element in self.elements:
                if element.next == an_element.name:
                    return element

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
                 d_x=0, d_y=0, d_z=0, next=None, previous=None, distance_from_previous=0):
        super().__init__()
        self.name = name
        self.phi = phi
        self.psi = psi
        self.theta = theta
        self.d_phi = d_phi
        self.d_psi = d_psi
        self.d_theta = d_theta
        self.x = x
        self.y = y
        self.z = z
        self.d_x = d_x
        self.d_y = d_y
        self.d_z = d_z
        self.next = next
        self.previous = previous
        self.distance_from_previous = distance_from_previous

    def __repr__(self):
        description = f"Element {self.name} of class {self.__class__}"
        description += f"\n\t at {self.distance_from_previous} m from {self.previous}"
        description += f"\n\t pointing to {self.next}"
        description += f"\n\t oriented in pitch at {self.theta/degree} deg (deviation {180-2*self.theta/degree} deg)"
        description += f"\n\t oriented in roll at {self.phi/degree} deg"
        description += f"\n\t oriented in yaw at {self.psi/degree} deg"
        return description

    def from_element_id(self, element_id, print_all=False):
        hparam = HANDLE(0)
        element_name = create_string_buffer(32)
        get_element_name(element_id, element_name, confirm=False)
        self.name = element_name.value.decode()
        param = Parameter()
        param_name = create_string_buffer(48)
        enumerate_parameters(element_id, hparam, param_name, param, confirm=False)
        while hparam:
            if print_all:
                print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
                            f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
            if param_name.value.decode() in optix_dictionnary.keys():
                self.__dict__[optix_dictionnary[param_name.value.decode()]] = param.value
            else:
                self.__dict__[param_name.value.decode()] = param.value
            enumerate_parameters(element_id, hparam, param_name, param, confirm=False)
        if print_all:
            print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
                        f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
        if param_name.value.decode() in optix_dictionnary.keys():
            self.__dict__[optix_dictionnary[param_name.value.decode()]] = param.value
        else:
            self.__dict__[param_name.value.decode()] = param.value
        enumerate_parameters(element_id, hparam, param_name, param, confirm=False)


class Source(OpticalElement):
    def __init__(self, name="", phi=0, psi=0, theta=0, d_phi=0, d_psi=0, d_theta=0, x=0, y=0, z=0,
                 d_x=0, d_y=0, d_z=0, next=None, previous=None, distance_from_previous=0, sigma_x=0,
                 sigma_y=0, sigma_x_div=0, sigma_y_div=0, nrays=0):
        super().__init__(name=name, phi=phi, psi=psi, theta=theta, d_phi=d_phi, d_psi=d_psi, d_theta=d_theta,
                         x=x, y=y, z=z, d_x=d_x, d_y=d_y, d_z=d_z, next=next, previous=previous,
                         distance_from_previous=distance_from_previous)
        self.sigma_x = sigma_x
        self.sigma_y = sigma_y
        self.sigma_x_div = sigma_x_div
        self.sigma_x_div = sigma_y_div
        self.nrays = nrays


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
