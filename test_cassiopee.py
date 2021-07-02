# -*- coding: utf-8 -*-
import ctypes, os
# from ctypes import cdll, windll
from ctypes.wintypes import *
from ctypes import *
from ctypes import WINFUNCTYPE, Structure, pointer, byref, POINTER,  c_char_p, c_void_p, c_double, c_float, \
                   create_string_buffer
import numpy as np
from classes import Beamline, OpticalElement, Parameter
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
    hsys, hparam, elemID = c_int64(0), c_int64(0), c_int64(0)
    elname = create_string_buffer(32, c_char)
    param_name = create_string_buffer(48)
    param = Parameter()
    print("#"*80)
    enumerate_elements(hsys, elemID, elname)
    sourceID = None
    while hsys:
        print("hsys", hsys)
        print("-"*20)
        print(f"element {elname.value.decode()}, ID {elemID.value}")
        if elname.value.decode() == "S_ONDUL1":
            sourceID = c_int64(elemID.value)
        enumerate_parameters(elemID, hparam, param_name, param, confirm=False)
        while hparam:
            print("\t", param_name.value, param.value, param.bounds.min, param.bounds.max, param.multiplier, param.type,
                  param.group, param.flags)
            enumerate_parameters(elemID, hparam, param_name, param, confirm=False)

        enumerate_elements(hsys, elemID, elname)
    print("Source ID", sourceID)
    lamda_align = c_double(2.5e-8)
    nrays = Parameter()
    get_parameter(sourceID, "nRays", nrays)
    print(nrays.value)
    nrays.value = 5000
    set_parameter(sourceID, "nRays", nrays)
    get_parameter(sourceID, "nRays", nrays)
    print(nrays.value)
    align(sourceID, lamda_align)
    generate(sourceID, lamda_align)
    radiate(sourceID)
