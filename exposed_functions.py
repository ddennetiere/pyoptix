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
    optix.FitSurfaceToSlopes.argtypes = (HANDLE, INT, INT, HANDLE, HANDLE, HANDLE, HANDLE)
    optix.FitSurfaceToSlopes.restype = BOOLEAN
    optix.FitSurfaceToHeights.argtypes = (HANDLE, INT, INT, HANDLE, HANDLE, HANDLE)
    optix.FitSurfaceToHeights.restype = BOOLEAN
    optix.GetHologramPatternInfo.argtypes = (HANDLE, HANDLE, DOUBLE, DOUBLE)
    optix.GetHologramPatternInfo.restype = BOOLEAN
    optix.GetImpactsData.argtypes = (HANDLE,  HANDLE, INT)
    optix.GetImpactsData.restype = INT
    optix.GetSpotDiagram.argtypes = (HANDLE, HANDLE, DOUBLE)
    optix.GetSpotDiagram.restype = INT
    optix.GetSpotDiagram.argtypes = (HANDLE, HANDLE, DOUBLE)
    optix.GetSpotDiagram.restype = INT
    optix.FindElementID.argtypes = [LPCSTR, HANDLE]
    optix.FindElementID.restype = INT
    optix.GetElementID.argtypes = [LPCSTR]
    optix.GetElementID.restype = HANDLE
    optix.ChainElement_byID.argtypes = [HANDLE, HANDLE]
    optix.ChainElement_byID.restype = INT
    optix.ChainElement_byName.argtypes = [LPCSTR, LPCSTR]
    optix.ChainElement_byName.restype = INT
    optix.GetPreviousElement.argtypes = [HANDLE]
    optix.GetPreviousElement.restype = HANDLE
    optix.GetNextElement.argtypes = [HANDLE]
    optix.GetNextElement.restype = HANDLE
    optix.ClearImpacts.argtypes = [HANDLE]
    optix.ClearImpacts.restype = INT
    optix.SetRecording.argtypes = [HANDLE, INT]
    optix.SetRecording.restype = BOOLEAN
    optix.GetElementName.argtypes = [HANDLE, HANDLE, INT]
    optix.GetElementName.restype = HANDLE
    optix.CreateElement.argtypes = [LPCSTR, LPCSTR]
    optix.CreateElement.restype = HANDLE
    optix.MemoryDump.argtypes = [c_voidp, ctypes.c_uint64]
    optix.DumpArgParameter.argtypes = [HANDLE]
    optix.DumpArgParameter.restype = LPVOID
    optix.DumpParameter.argtypes = [HANDLE, LPCSTR]
    optix.DumpParameter.restype = BOOLEAN
    optix.SetParameter.argtypes = [HANDLE, LPCSTR, HANDLE]
    optix.SetParameter.restype = BYTE
    optix.GetParameterArrayDims.argtypes = [HANDLE, LPCSTR, ctypes.POINTER(ARRAY(ctypes.c_int64, 2))]
    optix.GetParameterArrayDims.restype = BOOLEAN
    optix.GetParameterArraySize.argtypes = [HANDLE, LPCSTR, HANDLE]
    optix.GetParameterArraySize.restype = BOOLEAN
    optix.GetParameterFlags.argtypes = [HANDLE, LPCSTR, ctypes.POINTER(ctypes.c_uint32)]
    optix.GetParameterFlags.restype = BYTE
    optix.GetArrayParameter.argtypes = [HANDLE, LPCSTR, HANDLE, ctypes.c_ulonglong]
    optix.GetArrayParameter.restype = INT
    optix.GetParameter.argtypes = [HANDLE, LPCSTR, HANDLE]
    optix.GetParameter.restype = INT
    optix.EnumerateParameters.argtypes = [HANDLE, HANDLE, LPCSTR, INT, HANDLE]
    optix.EnumerateParameters.restype = INT
    optix.EnumerateElements.argtypes = [HANDLE, HANDLE, LPCSTR, INT]
    optix.EnumerateElements.restype = INT
    optix.LoadSystemFromXml.argtypes = [LPCSTR]
    optix.LoadSystemFromXml.restype = BOOLEAN
    optix.SaveSystemAsXml.argtypes = [LPCSTR]
    optix.SaveSystemAsXml.restype = BOOLEAN
    optix.Radiate.argtypes = [HANDLE]
    optix.Radiate.restype = BOOLEAN
    optix.GeneratePol.argtypes = [HANDLE, DOUBLE, ctypes.c_char]
    optix.GeneratePol.restype = INT
    optix.Generate.argtypes = [HANDLE, DOUBLE]
    optix.Generate.restype = INT
    optix.Align.argtypes = [HANDLE, DOUBLE]
    optix.Align.restype = INT
    optix.LoadSolemioFile.argtypes = [LPCSTR]
    optix.LoadSolemioFile.restype = INT
    ctypes.windll.kernel32.FreeLibrary.argtypes = [HMODULE]
    optix.GetStopNumber.argtypes = [HANDLE]
    optix.GetStopNumber.restype = HANDLE
    optix.GetStopType.argtypes = [HANDLE, HANDLE, ctypes.c_char_p, HANDLE]
    optix.GetStopType.restype = HANDLE
    optix.SetApertureActivity.argtypes = [HANDLE, BOOLEAN]
    optix.SetApertureActivity.restype = HANDLE
    optix.GetApertureActivity.argtypes = [HANDLE, POINTER(BOOLEAN)]
    optix.GetApertureActivity.restype = HANDLE
    optix.GetPolygonParameters.argtypes = [HANDLE, HANDLE, HANDLE, POINTER(c_double), POINTER(BOOLEAN)]
    optix.GetPolygonParameters.restype = HANDLE
    optix.AddPolygonalStop.argtypes = [HANDLE, HANDLE, POINTER(c_double), HANDLE]
    optix.AddPolygonalStop.restype = HANDLE
    optix.InsertPolygonalStop.argtypes = [HANDLE, HANDLE, HANDLE, POINTER(c_double), BOOLEAN]
    optix.InsertPolygonalStop.restype = HANDLE
    optix.ReplaceStopByPolygon.argtypes = [HANDLE, HANDLE, HANDLE, POINTER(c_double), BOOLEAN]
    optix.ReplaceStopByPolygon.restype = HANDLE
    optix.AddRectangularStop.argtypes = [HANDLE, c_double, c_double, BOOLEAN, c_double, c_double, c_double]
    optix.AddRectangularStop.restype = HANDLE
    optix.InsertRectangularStop.argtypes = [HANDLE, HANDLE, c_double, c_double, BOOLEAN, c_double, c_double, c_double]
    optix.InsertRectangularStop.restype = HANDLE
    optix.ReplaceStopByRectangle.argtypes = [HANDLE, HANDLE, c_double, c_double, BOOLEAN, c_double, c_double, c_double]
    optix.ReplaceStopByRectangle.restype = HANDLE
    optix.GetEllipseParameters.argtypes = [HANDLE, HANDLE, POINTER(c_double), POINTER(c_double), POINTER(BOOLEAN),
                                           POINTER(c_double), POINTER(c_double), POINTER(c_double)]
    optix.GetEllipseParameters.restype = HANDLE
    optix.AddEllipticalStop.argtypes = [HANDLE, c_double, c_double, BOOLEAN, c_double, c_double, c_double]
    optix.AddEllipticalStop.restype = HANDLE
    optix.InsertEllipticalStop.argtypes = [HANDLE, HANDLE, c_double, c_double, BOOLEAN, c_double, c_double, c_double]
    optix.InsertEllipticalStop.restype = HANDLE
    optix.ReplaceStopByEllipse.argtypes = [HANDLE, HANDLE, c_double, c_double, BOOLEAN, c_double, c_double, c_double]
    optix.ReplaceStopByEllipse.restype = HANDLE
    optix.AddCircularStop.argtypes = [HANDLE, c_double, BOOLEAN, c_double, c_double]
    optix.AddCircularStop.restype = HANDLE
    optix.InsertCircularStop.argtypes = [HANDLE, HANDLE, c_double, BOOLEAN, c_double, c_double]
    optix.InsertCircularStop.restype = HANDLE
    optix.ReplaceStopByCircle.argtypes = [HANDLE, HANDLE, c_double, BOOLEAN, c_double, c_double]
    optix.ReplaceStopByCircle.restype = HANDLE
    optix.SetAperturesActive.argtypes = [BOOLEAN]
    return optix


def release_optix():
    global optix
    libhandle = optix._handle
    ctypes.windll.kernel32.FreeLibrary(libhandle)
    del optix


@catch_c_error
def load_solemio_file(name):
    ret = optix.LoadSolemioFile(name)
    return ret


@catch_c_error
def align(source_id, lamda_align):
    ret = optix.Align(source_id, lamda_align)
    return ret


@catch_c_error
def generate(source_id, lamda):
    ret = optix.Generate(source_id, lamda)
    return ret


@catch_c_error
def generate_polarized(source_id, lamda, pol):
    # ret = optix.GeneratePol(source_id, lamda, 0x53)
    ret = optix.GeneratePol(source_id, lamda, ctypes.c_char(pol.encode()))
    return ret


@catch_c_error
def radiate(source_id):
    ret = optix.Radiate(source_id)
    return ret


@catch_c_error
def save_as_xml(filename):
    ret = optix.SaveSystemAsXml(filename.encode())
    return ret


@catch_c_error
def load_from_xml(filename):
    ret = optix.LoadSystemFromXml(filename.encode())
    return ret


@catch_c_error
def enumerate_elements(hsys, element_id, element_name):
    ret = optix.EnumerateElements(ctypes.byref(hsys), ctypes.byref(element_id), element_name, 32)
    return ret


@catch_c_error
def enumerate_parameters(element_id, handle_param, parameter_name, parameter):
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
    ret = optix.GetParameter(element_id, parameter_name.encode(), ctypes.byref(parameter))
    return ret


@catch_c_error
def get_array_parameter(element_id, parameter_name, parameter_data, size):
    ret = optix.GetArrayParameter(element_id, parameter_name.encode(), ctypes.byref(parameter_data), size)
    return ret


@catch_c_error
def get_parameter_flags(element_id, parameter_name, flags):
    ret = optix.GetParameterFlags(element_id, parameter_name.encode(), ctypes.byref(flags))
    return ret


@catch_c_error
def get_array_parameter_size(element_id, parameter_name, size):
    ret = optix.GetParameterArraySize(element_id, parameter_name.encode(), size)
    return ret


@catch_c_error
def get_array_parameter_dims(element_id, parameter_name, dims):
    ret = optix.GetParameterArrayDims(element_id, parameter_name.encode(), ctypes.byref(dims))
    return ret


@catch_c_error
def set_parameter(element_id, parameter_name, parameter):
    param_ref = ctypes.byref(parameter)
    try:
        ret = optix.SetParameter(element_id, parameter_name.encode(), param_ref)
        return ret
    except OSError as e:
        raise OSError(e.__str__()+f"\nErreur sur parametre, ref {param_ref}\n", parameter)


@catch_c_error
def dump_parameter(element_id, parameter_name):
    ret = optix.DumpParameter(element_id, parameter_name.encode())
    return ret


@catch_c_error
def dump_arg_parameter(parameter):
    ret = optix.DumpArgParameter(ctypes.byref(parameter))
    return ret


@catch_c_error
def memory_dump(address, size):
    optix.MemoryDump(address, size)


@catch_c_error
def create_element(class_name, name):
    elem_id = optix.CreateElement(class_name.encode(), name.encode())
    return elem_id


@catch_c_error
def get_element_name(element_id, element_name):
    ret = optix.GetElementName(element_id, ctypes.byref(element_name), 32)
    return ret


@catch_c_error
def set_recording(element_id, recording_mode):
    ret = optix.SetRecording(element_id, recording_mode)
    return ret


@catch_c_error
def clear_impacts(element_id):
    ret = optix.ClearImpacts(element_id)
    return ret


@catch_c_error
def get_next_element(element_id):
    next_id = optix.GetNextElement(element_id)
    return next_id


@catch_c_error
def get_previous_element(element_id):
    previous_id = optix.GetPreviousElement(element_id)
    return previous_id


@catch_c_error
def chain_element_by_name(previous_name, next_name):
    ret = optix.ChainElement_byName(previous_name.encode(), next_name.encode())
    return ret


@catch_c_error
def chain_element_by_id(previous_id, next_id):
    ret = optix.ChainElement_byID(previous_id, next_id)
    return ret


@catch_c_error
def get_element_id(element_name):
    ret = optix.GetElementID(element_name.encode())
    return ret


@catch_c_error
def find_element_id(element_name, element_id):
    ret = optix.FindElementID(element_name.encode(), ctypes.byref(element_id))
    return ret


@catch_c_error
def find_next_element(element_id, next_id):
    ret = optix.FindNextElement(element_id, ctypes.byref(next_id))
    return ret


@catch_c_error
def get_spot_diagram(element_id, diagram, distance):
    ret = optix.GetSpotDiagram(element_id, ctypes.byref(diagram), DOUBLE(distance))
    return ret


@catch_c_error
def get_impacts_data(element_id, diagram, frame_id):
    ret = optix.GetImpactsData(element_id, ctypes.byref(diagram), frame_id)
    return ret


@catch_c_error
def set_transmissive(element_id, is_transmissive):
    ret = optix.SetTransmissive(element_id, is_transmissive)
    return ret


@catch_c_error
def get_transmissive(element_id):
    ret = optix.GetTransmissive(element_id)
    return ret


@catch_c_error
def get_hologram_pattern(element_id, gratinfo, half_length, half_width):
    ret = optix.GetHologramPatternInfo(element_id, ctypes.byref(gratinfo), DOUBLE(half_length), DOUBLE(half_width))
    return ret


@catch_c_error
def fit_surface_to_heights(element_id, order_x, order_y, limits, heights, sigma_h):
    ret = optix.FitSurfaceToHeights(element_id, order_x, order_y, ctypes.byref(limits), ctypes.byref(heights),
                                    ctypes.byref(sigma_h))
    return ret


@catch_c_error
def fit_surface_to_slopes(element_id, order_x, order_y, limits, slopes, sigma_x, sigma_y):
    ret = optix.FitSurfaceToSlopes(element_id, order_x, order_y, ctypes.byref(limits), ctypes.byref(slopes),
                                   ctypes.byref(sigma_x), ctypes.byref(sigma_y))
    return ret


@catch_c_error
def get_stop_number(element_id):
    """Retrieves the number of stop regions composing the ApertureStop of a given OpticalElement."""
    return optix.GetStopNumber(element_id)


@catch_c_error
def get_stop_type(element_id, index):
    """Gets the runtime class name of the Region object."""
    buffer = ctypes.create_string_buffer(256)
    optix.GetStopType(element_id, index, buffer, 256)
    return buffer.value.decode()


def enumerate_stops(element_id):
    """Enumerates the stops of a given optical element"""
    stops = []
    for index in range(get_stop_number(element_id)):
        stops.append(get_stop_type(element_id, index))
    return stops


@catch_c_error
def set_aperture_activity(element_id, status):
    """Sets the activity status of the aperture stop of the optical element."""
    return optix.SetApertureActivity(element_id, status)


@catch_c_error
def get_aperture_activity(element_id):
    """Gets the activity status of the aperture stop of the optical element."""
    status = BOOLEAN()
    result = optix.GetApertureActivity(element_id, status)
    return result, status.value


@catch_c_error
def get_polygon_parameters(element_id, index, array_width):
    """Returns the parameters defining a polygonal Region."""
    vertex_array = (c_double * array_width)()
    opacity = BOOLEAN()
    result = optix.GetPolygonParameters(element_id, index, array_width, vertex_array, opacity)
    return result, list(vertex_array), opacity.value


@catch_c_error
def add_polygonal_stop(element_id, num_sides, vertex_array, opacity):
    """Append a new polygonal Region in the aperture stop list."""
    vertex_array = (c_double * len(vertex_array))(*vertex_array)
    return optix.AddPolygonalStop(element_id, num_sides, vertex_array, opacity)


@catch_c_error
def insert_polygonal_stop(element_id, index, num_sides, vertex_array, opacity):
    """Insert a new polygonal region in the aperture stop list."""
    vertex_array = (c_double * len(vertex_array))(*vertex_array)
    return optix.InsertPolygonalStop(element_id, index, num_sides, vertex_array, opacity)


@catch_c_error
def replace_stop_by_polygon(element_id, index, num_sides, vertex_array, opacity):
    """Replace a Region in the aperture stop list by a new Polygon."""
    vertex_array = (c_double * len(vertex_array))(*vertex_array)
    return optix.ReplaceStopByPolygon(element_id, index, num_sides, vertex_array, opacity)


@catch_c_error
def add_rectangular_stop(element_id, x_width, y_width, opacity, x_center, y_center, angle):
    """Add a new Rectangle at the end of the aperture stop list."""
    return optix.AddRectangularStop(element_id, x_width, y_width, opacity, x_center, y_center, angle)


@catch_c_error
def insert_rectangular_stop(element_id, index, x_width, y_width, opacity, x_center, y_center, angle):
    """Insert a new Rectangle at the given position in the aperture stop list."""
    return optix.InsertRectangularStop(element_id, index, x_width, y_width, opacity, x_center, y_center, angle)


@catch_c_error
def replace_stop_by_rectangle(element_id, index, x_width, y_width, opacity, x_center, y_center, angle):
    """Replace the Region at the given position in the aperture stop list by a new Rectangle."""
    return optix.ReplaceStopByRectangle(element_id, index, x_width, y_width, opacity, x_center, y_center, angle)


@catch_c_error
def get_ellipse_parameters(element_id, index):
    """Returns the parameters defining an elliptical Region."""
    x_axis = c_double()
    y_axis = c_double()
    opacity = BOOLEAN()
    x_center = c_double()
    y_center = c_double()
    angle = c_double()
    result = optix.GetEllipseParameters(element_id, index, ctypes.byref(x_axis), ctypes.byref(y_axis),
                                        ctypes.byref(opacity), ctypes.byref(x_center), ctypes.byref(y_center),
                                        ctypes.byref(angle))
    return result, x_axis.value, y_axis.value, opacity.value, x_center.value, y_center.value, angle.value


@catch_c_error
def add_elliptical_stop(element_id, x_axis, y_axis, opacity, x_center, y_center, angle):
    """Adds a new elliptical Region at the end of the aperture stop list."""
    return optix.AddEllipticalStop(element_id, x_axis, y_axis, opacity, x_center, y_center, angle)


@catch_c_error
def insert_elliptical_stop(element_id, index, x_axis, y_axis, opacity, x_center, y_center, angle):
    """Inserts a new elliptical Region at the given position in the aperture stop list."""
    return optix.InsertEllipticalStop(element_id, index, x_axis, y_axis, opacity, x_center, y_center, angle)


@catch_c_error
def replace_stop_by_ellipse(element_id, index, x_axis, y_axis, opacity, x_center, y_center, angle):
    """Replace the Region at the given position in the aperture stop list by a new Ellipse."""
    return optix.ReplaceStopByEllipse(element_id, index, x_axis, y_axis, opacity, x_center, y_center, angle)


@catch_c_error
def add_circular_stop(element_id, radius, opacity, x_center, y_center):
    """Adds a new circular Region at the end of the aperture stop list."""
    return optix.AddCircularStop(element_id, radius, opacity, x_center, y_center)


@catch_c_error
def insert_circular_stop(element_id, index, radius, opacity, x_center, y_center):
    """Insert a new circular Region at the given position in the aperture stop list."""
    return optix.InsertCircularStop(element_id, index, radius, opacity, x_center, y_center)


@catch_c_error
def replace_stop_by_circle(element_id, index, radius, opacity, x_center, y_center):
    """Replace the Region at the given position in the aperture stop list by a new circle."""
    return optix.ReplaceStopByCircle(element_id, index, radius, opacity, x_center, y_center)


def set_aperture_active(activity):
    """switch the global aperture activity flag """
    optix.SetAperturesActive(activity)


@catch_c_error
def version():
    optix.Version()
