import functools
import ctypes
from classes import Beamline, OpticalElement
from lxml import etree
from ctypes.wintypes import BYTE, INT, HMODULE

global optix


def init_optix():
    # initialisation auto
    try:
        test = optix
        print(test, "  already initialized")
    except NameError:
        load_optix()
        print("OptiX library initialized")


def catch_c_error(function):
    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        confirm = True
        if "confirm" in kwargs.keys():
            confirm = kwargs["confirm"]
        if function(*args):
            buf = ctypes.create_string_buffer(256)  # create a 128 byte buffer
            optix.GetOptiXLastError(buf, 256)
            print(function.__name__, "error", buf.value)
        else:
            if confirm:
                print(function.__name__, "ok")

    return wrapper


def load_optix():
    global optix
    print("intialzing SR library")
    optix = ctypes.cdll.LoadLibrary(r'D:\Dennetiere\optix\release\OptiX.dll')
    # optix.LoadSolemioFile.restype = BYTE  # gcc bool match to ctypes BYTE
    # optix.LoadSolemioFile.argtypes = [LPCSTR]
    optix.GetOptiXLastError.restype = BYTE
    optix.GetOptiXLastError.argtypes = [ctypes.c_char_p, INT]
    return optix


def release_optix():
    global optix
    ctypes.windll.kernel32.FreeLibrary.argtypes = [HMODULE]
    libhandle = optix._handle
    ctypes.windll.kernel32.FreeLibrary(libhandle)
    # del SRlib


@catch_c_error
def load_solemio_file(name):
    optix.LoadSolemioFile(name)


@catch_c_error
def align(sourceID, lamda_align):
    optix.Align(sourceID, lamda_align)


@catch_c_error
def generate(sourceID, lamda):
    optix.Align(sourceID, lamda)


@catch_c_error
def radiate(sourceID):
    optix.Align(sourceID)


@catch_c_error
def enumerate_elements(hsys, element_id, element_name):
    optix.EnumerateElements(ctypes.byref(hsys), ctypes.byref(element_id), element_name, 32)


@catch_c_error
def enumerate_parameters(element_id, handle_param, parameter_name, parameter, confirm=False):
    optix.EnumerateParameters(element_id, ctypes.byref(handle_param), parameter_name, 48, ctypes.byref(parameter))


@catch_c_error
def get_parameter(element_id, parameter_name, parameter):
    optix.GetParameter(element_id, parameter_name.encode(), ctypes.byref(parameter))


@catch_c_error
def set_parameter(element_id, parameter_name, parameter):
    optix.SetParameter(element_id, parameter_name.encode(), parameter)


@catch_c_error
def version():
    optix.Version()


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
