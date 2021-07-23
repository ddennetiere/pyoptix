# -*- coding: utf-8 -*-
import ctypes, os
# from ctypes import cdll, windll
from ctypes.wintypes import *
from ctypes import *
from ctypes import WINFUNCTYPE, Structure, pointer, byref, POINTER, c_char_p, c_void_p, c_double, c_float, \
    create_string_buffer
import numpy as np
from classes import Beamline, OpticalElement, Parameter, Diagram, RecordingMode
from exposed_functions import *
from ui_objects import scatter_plot_2d, show, plot_spd

global optix

if __name__ == "__main__":
    test_ID = False
    test_enumerate = False
    test_parameter = True
    test_edit_parameter = False
    test_linkage = True
    test_add_element = True
    test_radiate = True
    test_spot_diagram = True
    # initialisation auto
    try:
        test = optix
        print(test, "already initialized")
    except NameError:
        optix = load_optix()
        print("OptiX library initialized")
    # parse_xml(r"D:\Dennetiere\optix\bin\test\system.xml")
    # optix.LoadSolemioFile.argtypes = [c_char_p]
    load_solemio_file(create_string_buffer(b"D:\\Dennetiere\\Programmes Python\\optix\\solemio\\CASSIOPEE"))
    hsys, hparam, elemID = HANDLE(0), HANDLE(0), HANDLE(0)
    elname = create_string_buffer(32, c_char)
    if test_parameter or test_edit_parameter:
        param_name = create_string_buffer(48)
        param = Parameter()
    print("#" * 80)
    enumerate_elements(hsys, elemID, elname)
    sourceID = None
    elements_ID = {}
    start = True
    last = False
    while hsys or start or last:
        start = False
        print("-" * 20)
        print("hsys", hsys)
        if test_enumerate or test_parameter:
            print(f"element {elname.value.decode()}, ID {elemID.value}")
        elements_ID[elname.value.decode()] = HANDLE(elemID.value)
        if test_parameter:
            enumerate_parameters(elemID, hparam, param_name, param, confirm=False)
            while hparam:
                print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
                            f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
                enumerate_parameters(elemID, hparam, param_name, param, confirm=False)
            print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
                        f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
        if not last:
            enumerate_elements(hsys, elemID, elname, confirm=test_enumerate)
            if not hsys:
                last = True
        else:
            last = False

    if test_ID:
        print(elements_ID)
        print(elements_ID.keys())
        print("SourceID", elements_ID["S_ONDUL1"])
        print("S_ONDUL1 has ID", get_element_id("S_ONDUL1"))
        get_element_name(elements_ID['S_ONDUL1'], elname)
        print(f"Stored ID for S_ONDUL1 {elements_ID['S_ONDUL1']} is owned by {elname.value.decode()}")
        an_ID = HANDLE(0)
        find_element_id("S_ONDUL1", an_ID)
        print("find S_ONDUL1 ID:", an_ID)
    sourceID = elements_ID["S_ONDUL1"]
    lamda_align = DOUBLE(2.5e-8)

    if test_edit_parameter:
        nrays = Parameter()
        enumerate_parameters(sourceID, hparam, param_name, param, confirm=True)
        while hparam:
            print("\t", param_name.value, param.value, param.bounds.min, param.bounds.max, param.multiplier, param.type,
                  param.group, param.flags)
            enumerate_parameters(sourceID, hparam, param_name, param, confirm=False)

        get_parameter(sourceID, "nRays", nrays)
        print("Imported nrays", nrays.value)
        nrays.value = 50000
        set_parameter(sourceID, "nRays", nrays)
        nrays2 = Parameter()
        get_parameter(sourceID, "nRays", nrays2)
        print("New nrays", nrays2.value)

    if test_linkage:
        linked_beamline = []
        get_element_name(sourceID, elname)
        print("source name:", elname.value.decode())
        next_ID = sourceID
        this_ID = HANDLE(sourceID.value)
        while next_ID:
            get_element_name(this_ID, elname, show_return=False, confirm=False)
            linked_beamline.append(elname.value.decode())
            next_ID = get_next_element(this_ID, show_return=False, confirm=False)
            this_ID = next_ID
        print("Chained beamline :", linked_beamline)
        for oe_name in linked_beamline:
            oe = OpticalElement(element_id=elements_ID[oe_name])
            print(oe)

    if test_add_element:
        elements_ID["screen"] = create_element("PlaneFilm", "screen")
        print("screen_id", elements_ID["screen"])
        get_element_name(elements_ID["screen"], elname)
        print("screen name", elname.value.decode())
        chain_element_by_id(elements_ID["f"], elements_ID["screen"])
    if test_spot_diagram:
        assert test_add_element
        set_recording(elements_ID["screen"], RecordingMode.recording_output)

    if test_radiate:
        print(lamda_align)
        align(sourceID, lamda_align, show_return=True)
        clear_impacts(sourceID)
        generated_rays = generate(sourceID, lamda_align, show_return=True)
        radiate(sourceID)

    if test_spot_diagram:
        assert test_add_element
        if not test_edit_parameter:
            nrays2 = Parameter()
            get_parameter(sourceID, "nRays", nrays2)
            print("New nrays", nrays2.value)

        # new_distance = Parameter()
        # get_parameter(elements_ID["pupille"], "distance", new_distance)
        # new_distance.value = DOUBLE(20-0.9)
        # set_parameter(elements_ID["pupille"], "distance", new_distance)
        # get_parameter(elements_ID["pupille"], "distance", new_distance)
        # print("New distance", new_distance.value)
        diagram = Diagram(ndim=5, nreserved=int(nrays2.value))
        get_spot_diagram(elements_ID["screen"], diagram, distance=20)
        print(np.ctypeslib.as_array(diagram.min, shape=(diagram.dim,)))
        print(diagram.count)
        import pandas as pd
        spots = pd.DataFrame(np.ctypeslib.as_array(diagram.spots, shape=(diagram.reserved, diagram.dim)),
                             columns=("X", "Y", "dX", "dY", "Lambda"))
        # import matplotlib.pyplot as plt
        # plt.scatter(spots["X"], spots["Y"])
        # plt.show()
        print(spots.head())
        figs = []
        figs.append(plot_spd(spots, x_key="X", y_key="Y", light_plot=True, show_map=True))
        figs.append(plot_spd(spots, x_key="X", y_key="dX", light_plot=False))
        figs.append(plot_spd(spots, x_key="Y", y_key="dY", light_plot=True))

        for fig in figs:
            show(fig)

    if test_ID:
        print("ran test_ID")
    if test_enumerate:
        print("ran test_enumerate")
    if test_parameter:
        print("ran test_parameter")
    if test_edit_parameter:
        print("ran test_edit_parameter")
    if test_linkage:
        print("ran test_linkage")
    if test_add_element:
        print("ran test_add_element")
    if test_radiate:
        print("ran test_radiate")
    if test_spot_diagram:
        print("ran test_spot_diagram")
