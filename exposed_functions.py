import functools
import ctypes
from ctypes import ARRAY, c_voidp, c_double, Structure, c_int64, c_int32, POINTER, Union, pointer, c_int, cast
from ctypes.wintypes import BYTE, INT, HMODULE, LPCSTR, HANDLE, DOUBLE, BOOLEAN, LARGE_INTEGER, LPVOID, UINT
import numpy as np

global optix


class Bounds(Structure):
    """
    C structure to be used in optix parameters
    """
    _pack_ = 1
    _fields_ = [("min", c_double),
                ("max", c_double)]

    def __repr__(self):
        return f"Parameter bounds [{self.min}, {self.max}]"


class ParamArray(Structure):
    _fields_ = [('dims', c_int64 * 2),
                ('data', POINTER(DOUBLE))]


class UData(Union):
    # _pack_ = 8
    _fields_ = [("value", DOUBLE),
                ("p_array", POINTER(ParamArray))]


class Dims(Structure):
    """
    C structure to be used in optix parameters
    """
    _pack_ = 1
    _fields_ = [("x", c_int64),
                ("y", c_int64)]


class Parameter(Structure):
    """
    C structure defining modifiable fields of optix optical element parameters. Note bounds type is Bounds. See Bounds
    docstring.
    """
    # _pack_ = 8
    _anonymous_ = ("u",)
    _fields_ = [("u", UData),
                ("bounds", Bounds),
                ("multiplier", DOUBLE),
                ("type", INT),
                ("group", INT),
                ("flags", UINT),
                ]

    def __init__(self):
        super().__init__()
        self._array = None

    @property
    def array(self):
        return self._array

    @array.setter
    def array(self, np_array):
        """
        np_array doit être un numpy ndarray
        this function sets self.array, the array flag bit, and the required pointer and dims fields of
        the Parameter.paramArray sub structure

        :param np_array: array value
        :type np_array: numpy.ndarray
        :return: None
        :rtype: Nonetype
        """
        if not isinstance(np_array, np.ndarray):
            print("a 'ndarray' argument was expected ,but the received argument type is ", type(np_array).__name__)
            return
        if np_array.dtype.name == 'float64':
            self._array = np_array
        else:
            self._array = np_array.astype('float64')
        pa = ParamArray()
        pa.dims[0] = self._array.shape[1]
        pa.dims[1] = self._array.shape[0]
        pa.data = self._array.ctypes.data_as(POINTER(c_double))
        self.p_array = pointer(pa)
        # pa peut être retrouvé par pa=p_array.contents
        self.flags |= 1 << 3

    def __repr__(self):
        if self.flags & 1 << 3:
            data_length = self.p_array.contents.dims[0]*self.p_array.contents.dims[1]
            data_address = self.p_array.contents.data.contents
            data = np.ctypeslib.as_array((ctypes.c_double * data_length).from_address(ctypes.addressof(data_address)))
            return f"Param object containing:\n" \
                   f"\t ParamArray dimension {self.p_array.contents.dims[0]}*{self.p_array.contents.dims[1]}, " \
                   f"@ {hex(ctypes.addressof(self.p_array))}\n" \
                   f"\t\twith data = {data} @ {hex(ctypes.addressof(data_address))}\n" \
                   f"\t Bounds {self.bounds}\n" \
                   f"\t multiplier {self.multiplier}\n" \
                   f"\t type {self.type}\n" \
                   f"\t group {self.group}\n" \
                   f"\t flags {self.flags}\n"
        else:
            return f"Param object containing:\n" \
                   f"\t Value {self.value}\n" \
                   f"\t Bounds {self.bounds}\n" \
                   f"\t multiplier {self.multiplier}\n" \
                   f"\t type {self.type}\n" \
                   f"\t group {self.group}\n" \
                   f"\t flags {self.flags}\n"


class PolynomialExpansion(Structure):
    """
    C structure to be used in optix parameters
    """
    _pack_ = 1
    _fields_ = [("a0", c_double),
                ("a1", c_double),
                ("a2", c_double),
                ("a3", c_double),
                ]


class GratingPatternInformation(Structure):
    """
    C structure defining the pattern parameters of an hologram used for holographic grating.
    Fields are the polynomial expansion of the line density, and the central line tilt and curvature of an holographic
    grating, where :

    - axial_line_density is an array which will receive the coefficients of the line density approximation by a
      third degree polynomial. The term of degree 0 is the nominal line density at grating center (in lines/m)
    - line_curvature is the angle of the angle of the central line of the grating line with the Y axis (in rad)
    - line_tilt is the angle of the angle of the central line of the grating line with the Y axis (in rad)

    """
    _fields_ = [("axial_line_density", PolynomialExpansion),
                ("line_tilt", c_double),
                ("line_curvature", c_double),
                ]


class FrameID(object):
    general_frame = c_int32(0)
    local_absolute_frame = c_int32(1)
    aligned_local_frame = c_int32(2)
    surface_frame = c_int32(3)


class RecordingMode(object):
    recording_output = c_int32(2)
    recording_input = c_int32(1)
    recording_none = c_int32(0)


class Diagram(Structure):
    """
    C structure for accessing optix spot diagrams
    """
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
    print("intialzing optix library")
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
def generate_polarized(source_id, lamda, pol):
    optix.GeneratePol.argtypes = [HANDLE, DOUBLE, ctypes.c_char]
    optix.GeneratePol.restype = INT
    # ret = optix.GeneratePol(source_id, lamda, 0x53)
    ret = optix.GeneratePol(source_id, lamda, ctypes.c_char(pol.encode()))
    return ret


@catch_c_error
def radiate(source_id):
    optix.Radiate.argtypes = [HANDLE]
    optix.Radiate.restype = BOOLEAN
    ret = optix.Radiate(source_id)
    return ret


@catch_c_error
def save_as_xml(filename):
    optix.SaveSystemAsXml.argtypes = [LPCSTR]
    optix.SaveSystemAsXml.restype = BOOLEAN
    ret = optix.SaveSystemAsXml(filename.encode())
    return ret


@catch_c_error
def load_from_xml(filename):
    optix.LoadSystemFromXml.argtypes = [LPCSTR]
    optix.LoadSystemFromXml.restype = BOOLEAN
    ret = optix.LoadSystemFromXml(filename.encode())
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


def safe_enumerate_parameters(element_id, handle_param, parameter_name, parameter):
    c_param = Parameter()
    ret = enumerate_parameters(element_id, handle_param, parameter_name, parameter)
    if ret:
        parameter.bounds.min = c_param.bounds.min
        parameter.bounds.max = c_param.bounds.max
        parameter.multiplier = c_param.multiplier
        parameter.type = c_param.type
        parameter.group = c_param.group
        parameter.flags = c_param.flags
        if c_param.flags & 1 << 3:
            pa = ParamArray()

            ppa_src = c_param.p_array
            cast(ppa_src, POINTER(ParamArray))
            pa_src = ppa_src.contents
            pa.dims[0] = pa_src.dims[0]
            pa.dims[1] = pa_src.dims[1]
            size = pa.dims[0] * pa.dims[1]
            databuf = np.fromiter(pa_src.data, np.float64, count=size)
            parameter._array = np.ndarray((pa.dims[1], pa.dims[0]), dtype=np.float64, buffer=databuf)
            pa.data = parameter._array.ctypes.data_as(POINTER(c_double))
            parameter.p_array = pointer(pa)
        else:
            parameter.value = c_param.value
    return ret


@catch_c_error
def get_parameter(element_id, parameter_name, parameter):
    optix.GetParameter.argtypes = [HANDLE, LPCSTR, HANDLE]
    optix.GetParameter.restype = INT
    ret = optix.GetParameter(element_id, parameter_name.encode(), ctypes.byref(parameter))
    return ret


@catch_c_error
def get_array_parameter(element_id, parameter_name, parameter_data, size):
    optix.GetArrayParameter.argtypes = [HANDLE, LPCSTR, HANDLE, ctypes.c_ulonglong]
    optix.GetArrayParameter.restype = INT
    ret = optix.GetArrayParameter(element_id, parameter_name.encode(), ctypes.byref(parameter_data), size)
    return ret


@catch_c_error
def get_parameter_flags(element_id, parameter_name, flags):
    optix.GetParameterFlags.argtypes = [HANDLE, LPCSTR, ctypes.POINTER(ctypes.c_uint32)]
    optix.GetParameterFlags.restype = BYTE
    ret = optix.GetParameterFlags(element_id, parameter_name.encode(), ctypes.byref(flags))
    return ret


@catch_c_error
def get_array_parameter_size(element_id, parameter_name, size):
    optix.GetParameterArraySize.argtypes = [HANDLE, LPCSTR, HANDLE]
    optix.GetParameterArraySize.restype = BOOLEAN
    ret = optix.GetParameterArraySize(element_id, parameter_name.encode(), size)
    return ret


@catch_c_error
def get_array_parameter_dims(element_id, parameter_name, dims):
    optix.GetParameterArrayDims.argtypes = [HANDLE, LPCSTR, ctypes.POINTER(ARRAY(ctypes.c_int64, 2))]
    optix.GetParameterArrayDims.restype = BOOLEAN
    ret = optix.GetParameterArrayDims(element_id, parameter_name.encode(), ctypes.byref(dims))
    return ret


@catch_c_error
def set_parameter(element_id, parameter_name, parameter):
    optix.SetParameter.argtypes = [HANDLE, LPCSTR, HANDLE]
    optix.SetParameter.restype = BYTE
    param_ref = ctypes.byref(parameter)
    try:
        ret = optix.SetParameter(element_id, parameter_name.encode(), param_ref)
        return ret
    except OSError as e:
        raise OSError(e.__str__()+f"\nErreur sur parametre, ref {param_ref}\n", parameter)


@catch_c_error
def dump_parameter(element_id, parameter_name):
    optix.DumpParameter.argtypes = [HANDLE, LPCSTR]
    optix.DumpParameter.restype = BOOLEAN
    ret = optix.DumpParameter(element_id, parameter_name.encode())
    return ret


@catch_c_error
def dump_arg_parameter(parameter):
    optix.DumpArgParameter.argtypes = [HANDLE]
    optix.DumpArgParameter.restype = LPVOID
    ret = optix.DumpArgParameter(ctypes.byref(parameter))
    return ret


@catch_c_error
def memory_dump(address, size):
    optix.MemoryDump.argtypes = [c_voidp, ctypes.c_uint64]
    optix.MemoryDump(address, size)


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
def get_hologram_pattern(element_id, gratinfo, half_length, half_width):
    optix.GetHologramPatternInfo.argtypes = (HANDLE, HANDLE, DOUBLE, DOUBLE)
    optix.GetHologramPatternInfo.restype = BOOLEAN
    ret = optix.GetHologramPatternInfo(element_id, ctypes.byref(gratinfo), DOUBLE(half_length), DOUBLE(half_width))
    return ret


@catch_c_error
def fit_surface_to_heights(element_id, order_x, order_y, limits, heights, sigma_h):
    optix.FitSurfaceToHeights.argtypes = (HANDLE, INT, INT, HANDLE, HANDLE, HANDLE)
    optix.FitSurfaceToHeights.restype = BOOLEAN
    ret = optix.FitSurfaceToHeights(element_id, order_x, order_y, ctypes.byref(limits), ctypes.byref(heights),
                                    ctypes.byref(sigma_h))
    return ret


@catch_c_error
def fit_surface_to_slopes(element_id, order_x, order_y, limits, slopes, sigma_x, sigma_y):
    optix.FitSurfaceToSlopes.argtypes = (HANDLE, INT, INT, HANDLE, HANDLE, HANDLE, HANDLE)
    optix.FitSurfaceToSlopes.restype = BOOLEAN
    ret = optix.FitSurfaceToSlopes(element_id, order_x, order_y, ctypes.byref(limits), ctypes.byref(slopes),
                                   ctypes.byref(sigma_x), ctypes.byref(sigma_y))
    return ret


@catch_c_error
def version():
    optix.Version()
