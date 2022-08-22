import functools
import ctypes
from ctypes.wintypes import BYTE, INT, HMODULE, LPCSTR, HANDLE, DOUBLE, BOOLEAN

global optix


def init_optix():
    # initialisation auto
    try:
        test = optix
        print(test, "already initialized")
    except NameError:
        load_optix()
        print("OptiX library initialized")


def catch_c_error(function):
    @functools.wraps(function)
    def wrapper(*args, **kwargs):
        confirm = True
        confirm_ok = False
        show_return = False
        if "confirm" in kwargs:
            confirm = kwargs["confirm"]
        if "confirm_ok" in kwargs:
            confirm_ok = kwargs["confirm_ok"]
        if "show_return" in kwargs:
            show_return = kwargs["show_return"]
        ret = function(*args)
        if show_return:
            print(function.__name__, "returned", ret)
        if not ret:
            buf = ctypes.create_string_buffer(256)  # create a 128 byte buffer
            optix.GetOptiXLastError(buf, 256)
            print(function.__name__, args, kwargs, "error", buf.value)
        else:
            if confirm and confirm_ok:
                print(function.__name__, "ok")
        return ret

    return wrapper


def load_optix():
    global optix
    print("intialzing SR library")
    optix = ctypes.cdll.LoadLibrary(r'D:\Dennetiere\optix\release\OptiX.dll')
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
    optix.LoadSolemioFile.argtypes = [LPCSTR]
    optix.LoadSolemioFile.restype = INT
    ret = optix.LoadSolemioFile(name)
    return ret


@catch_c_error
def align(source_id, lamda_align):
    optix.Align.argtypes = [HANDLE, DOUBLE]
    optix.Align.restype = INT
    ret = optix.Align(source_id, lamda_align)
    return ret


@catch_c_error
def generate(source_id, lamda):
    optix.Generate.argtypes = [HANDLE, DOUBLE]
    optix.Generate.restype = INT
    ret = optix.Generate(source_id, lamda)
    return ret


@catch_c_error
def radiate(source_id):
    optix.Radiate.argtypes = [HANDLE]
    optix.Radiate.restype = BOOLEAN
    ret = optix.Radiate(source_id)
    return ret


@catch_c_error
def enumerate_elements(hsys, element_id, element_name):
    optix.EnumerateElements.argtypes = [HANDLE, HANDLE, LPCSTR, INT]
    optix.EnumerateElements.restype = INT
    ret = optix.EnumerateElements(ctypes.byref(hsys), ctypes.byref(element_id), element_name, 32)
    return ret


@catch_c_error
def enumerate_parameters(element_id, handle_param, parameter_name, parameter):
    optix.EnumerateParameters.argtypes = [HANDLE, HANDLE, LPCSTR, INT, HANDLE]
    optix.EnumerateParameters.restype = INT
    ret = optix.EnumerateParameters(element_id, ctypes.byref(handle_param), parameter_name, 48, ctypes.byref(parameter))
    return ret


@catch_c_error
def get_parameter(element_id, parameter_name, parameter):
    optix.GetParameter.argtypes = [HANDLE, LPCSTR, HANDLE]
    optix.GetParameter.restype = INT
    ret = optix.GetParameter(element_id, parameter_name.encode(), ctypes.byref(parameter))
    return ret


@catch_c_error
def set_parameter(element_id, parameter_name, parameter):
    optix.SetParameter.argtypes = [HANDLE, LPCSTR, HANDLE]
    optix.SetParameter.restype = INT
    ret = optix.SetParameter(element_id, parameter_name.encode(), ctypes.byref(parameter))
    return ret


@catch_c_error
def create_element(class_name, name):
    optix.CreateElement.argtypes = [LPCSTR, LPCSTR]
    optix.CreateElement.restype = HANDLE
    elem_id = optix.CreateElement(class_name.encode(), name.encode())
    return elem_id


@catch_c_error
def get_element_name(element_id, element_name):
    optix.GetElementName.argtypes = [HANDLE, HANDLE, INT]
    optix.GetElementName.restype = HANDLE
    ret = optix.GetElementName(element_id, ctypes.byref(element_name), 32)
    return ret


@catch_c_error
def set_recording(element_id, recording_mode):
    optix.SetRecording.argtypes = [HANDLE, INT]
    optix.SetRecording.restype = BOOLEAN
    ret = optix.SetRecording(element_id, recording_mode)
    return ret


@catch_c_error
def clear_impacts(element_id):
    optix.ClearImpacts.argtypes = [HANDLE]
    optix.ClearImpacts.restype = INT
    ret = optix.ClearImpacts(element_id)
    return ret


@catch_c_error
def get_next_element(element_id):
    optix.GetNextElement.argtypes = [HANDLE]
    optix.GetNextElement.restype = HANDLE
    next_id = optix.GetNextElement(element_id)
    return next_id


@catch_c_error
def get_previous_element(element_id):
    optix.GetPreviousElement.argtypes = [HANDLE]
    optix.GetPreviousElement.restype = HANDLE
    previous_id = optix.GetPreviousElement(element_id)
    return previous_id


@catch_c_error
def chain_element_by_name(previous_name, next_name):
    optix.ChainElement_byName.argtypes = [LPCSTR, LPCSTR]
    optix.ChainElement_byName.restype = INT
    ret = optix.ChainElement_byName(previous_name.encode(), next_name.encode())
    return ret


@catch_c_error
def chain_element_by_id(previous_id, next_id):
    optix.ChainElement_byID.argtypes = [HANDLE, HANDLE]
    optix.ChainElement_byID.restype = INT
    ret = optix.ChainElement_byID(previous_id, next_id)
    return ret


@catch_c_error
def get_element_id(element_name):
    optix.GetElementID.argtypes = [LPCSTR]
    optix.GetElementID.restype = HANDLE
    ret = optix.GetElementID(element_name.encode())
    return ret


@catch_c_error
def find_element_id(element_name, element_id):
    optix.FindElementID.argtypes = [LPCSTR, HANDLE]
    optix.FindElementID.restype = INT
    ret = optix.FindElementID(element_name.encode(), ctypes.byref(element_id))
    return ret


@catch_c_error
def find_next_element(element_id, next_id):
    optix.FindElementID.argtypes = [LPCSTR, HANDLE]
    optix.FindElementID.restype = INT
    ret = optix.FindNextElement(element_id, ctypes.byref(next_id))
    return ret


@catch_c_error
def get_spot_diagram(element_id, diagram, distance):
    optix.GetSpotDiagram.argtypes = (HANDLE, HANDLE, DOUBLE)
    optix.GetSpotDiagram.restype = INT
    ret = optix.GetSpotDiagram(element_id, ctypes.byref(diagram), DOUBLE(distance))
    return ret


@catch_c_error
def get_impacts_data(element_id, diagram, frame_id):
    optix.GetImpactsData.argtypes = (HANDLE, HANDLE, INT)
    optix.GetImpactsData.restype = INT
    ret = optix.GetImpactsData(element_id, ctypes.byref(diagram), frame_id)
    return ret


@catch_c_error
def set_transmissive(element_id, is_transmissive):
    optix.GetImpactsData.argtypes = (HANDLE,  BOOLEAN)
    optix.GetImpactsData.restype = INT
    ret = optix.SetTransmissive(element_id, is_transmissive)
    return ret


@catch_c_error
def get_transmissive(element_id):
    optix.GetImpactsData.argtypes = HANDLE
    optix.GetImpactsData.restype = BOOLEAN
    ret = optix.GetTransmissive(element_id)
    return ret

@catch_c_error
def version():
    optix.Version()
