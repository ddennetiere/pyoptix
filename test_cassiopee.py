# -*- coding: utf-8 -*-
import ctypes, os
# from ctypes import cdll, windll
from ctypes.wintypes import *
from ctypes import *
from ctypes import WINFUNCTYPE, Structure, pointer, byref, POINTER,  c_char_p, c_void_p, c_double, c_float, \
                   create_string_buffer
import numpy as np
from classes import Beamline, OpticalElement, Parameter, Diagram, RecordingMode
from exposed_functions import *


global optix

if __name__ == "__main__":
    # initialisation auto
    try:
        test = optix
        print(test, "  already initialized")
    except NameError:
        optix = load_optix()
        print("OptiX library initialized")
    # parse_xml(r"D:\Dennetiere\optix\bin\test\system.xml")
    # optix.LoadSolemioFile.argtypes = [c_char_p]
    load_solemio_file(create_string_buffer(b"D:\\Dennetiere\\Programmes Python\\optix\\solemio\\CASSIOPEE"))
    hsys, hparam, elemID = c_int64(0), c_int64(0), c_uint64(0)
    elname = create_string_buffer(32, c_char)
    param_name = create_string_buffer(48)
    param = Parameter()
    print("#"*80)
    # enumerate_elements(hsys, elemID, elname)
    sourceID = None
    elements_ID = {}
    elements_ID["screen"] = create_element("PlaneFilm", "screen")
    start = 1
    while hsys or start:
        start = 0
        print("-"*20)
        print("hsys", hsys)
        enumerate_elements(hsys, elemID, elname)
        print(f"element {elname.value.decode()}, ID {elemID.value}")
        elements_ID[elname.value.decode()] = c_int64(elemID.value)
        enumerate_parameters(elemID, hparam, param_name, param, confirm=False)
        while hparam:
            # print("\t", param_name.value, param.value, param.bounds.min, param.bounds.max, param.multiplier, param.type,
            #       param.group, param.flags)
            enumerate_parameters(elemID, hparam, param_name, param, confirm=False)

    print(elements_ID)
    print(elements_ID.keys())
    print("SourceID", elements_ID["S_ONDUL1"])
    sourceID = elements_ID["S_ONDUL1"]
    lamda_align = c_double(2.5e-8)
    nrays = Parameter()
    enumerate_parameters(sourceID, hparam, param_name, param, confirm=False)
    while hparam:
        # print("\t", param_name.value, param.value, param.bounds.min, param.bounds.max, param.multiplier, param.type,
        #       param.group, param.flags)
        enumerate_parameters(sourceID, hparam, param_name, param, confirm=False)
    get_parameter(sourceID, "nRays", nrays)
    print(nrays.value)
    nrays.value = 10000
    set_parameter(sourceID, "nRays", nrays)
    nrays2 = Parameter()
    get_parameter(sourceID, "nRays", nrays2)
    print(nrays2.value)

    # elements_ID["screen"] = create_element("PlaneFilm", "screen")
    print("screen_id", elements_ID["screen"])
    print("pupille_id", elements_ID["pupille"])
    get_element_name(elements_ID["screen"], elname)
    # optix.GetElementName(elements_ID["screen"], ctypes.byref(elname), c_int32(32))
    print("screen name", elname.value.decode())
    set_recording(elements_ID["screen"], RecordingMode.recording_output)
    # optix.SetRecording(elements_ID["screen"], RecordingMode.recording_output)
    diagram = Diagram(ndim=5, nreserved=int(nrays2.value))

    align(sourceID, lamda_align)
    clear_impacts(sourceID)
    generate(sourceID, lamda_align)
    radiate(sourceID)
    optix.GetSpotDiagram(elements_ID["screen"],  ctypes.byref(diagram), c_double(0))
    print(diagram.min.value)
