# -*- coding: utf-8 -*-
"""
Test file to show the import of a preexisting SOLEMIO file, CASSIOPEE in this instance.

This file should be called with optionnal argument depending on what the user wants to see working:

python test_cassiopee.py [--optionnal_test]

options are :
--test_ID
    Shows internal handles of optix imported objects
--test_enumerate
    Enumerates all optical element in the beamline
--test_parameter
    Enumerates all parameter of all the optical element in the beamline
--test_edit_parameter
    Changes the value of an optix parameter and prints the new stored value
--test_linkage
    Links the optical element as they are linked in the beamline SOLEMIO file
--test_add_element
    Adds an element that does not exist in the SOLEMIO file (a film here)
--test_radiate
    Generates and propagates rays from the source to the EXP1 film
--test_spot_diagram
    Shows the different useful spot diagrams in the EXP1 plane (XY, XX', YY')

"""
from ctypes import create_string_buffer, c_char, sizeof
import numpy as np
import sys
sys.path.append("../..")
from pyoptix import classes
from pyoptix.exposed_functions import (load_optix, load_solemio_file, HANDLE, enumerate_elements, enumerate_parameters,
                                      get_element_id, get_element_name, find_element_id, DOUBLE, get_parameter,
                                      set_parameter, get_next_element, create_element, chain_element_by_id,
                                      set_recording, align , radiate, generate, clear_impacts, get_spot_diagram)
from pyoptix import ui_objects

global optix

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--test_ID", help="shows the elements optix ID", action="store_true")
    parser.add_argument("--test_enumerate", help="shows the elements in the beamline", action="store_true")
    parser.add_argument("--test_parameter", help="shows the elements parameters", action="store_true")
    parser.add_argument("--test_edit_parameter", help="changes the source Nrays", action="store_true")
    parser.add_argument("--test_linkage",
                        help="shows the linkage of the elements through which light will be propagated",
                        action="store_true")
    parser.add_argument("--test_beamline",
                        help="Turns the linked beamline into a pyoptix Beamline object",
                        action="store_true")
    parser.add_argument("--test_add_element", help="adds a screen in the beamline", action="store_true")
    parser.add_argument("--test_radiate", help="runs the raytracing", action="store_true")
    parser.add_argument("--test_spot_diagram",
                        help="shows the spot diagram on the added screen, requires --test_add_element",
                        action="store_true")
    args = parser.parse_args()

    test_ID = args.test_ID
    test_enumerate = args.test_enumerate
    test_parameter = args.test_parameter
    test_edit_parameter = args.test_edit_parameter
    test_linkage = args.test_linkage
    test_beamline = args.test_beamline
    test_add_element = args.test_add_element
    test_radiate = args.test_radiate
    test_spot_diagram = args.test_spot_diagram
    generated_rays = None
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
    elname = create_string_buffer(32, sizeof(c_char))
    param_name = create_string_buffer(48)
    param = classes.Parameter()
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
        element_name = "Reseau_1600H"
        print("Paramètres de l'élément", element_name)
        theta = classes.Parameter()
        enumerate_parameters(elements_ID[element_name], hparam, param_name, param, confirm=True)
        while hparam:
            print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
                        f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")
            enumerate_parameters(elements_ID[element_name], hparam, param_name, param, confirm=False)
        print("\t", f"{param_name.value.decode()}: {param.value} [{param.bounds.min}, {param.bounds.max}],"
              f"x{param.multiplier}, type {param.type}, groupe {param.group}, flags {param.flags}")

        get_parameter(elements_ID["Reseau_1600H"], "theta", theta)
        print("Imported theta", theta.value)
        theta.value = 0.1
        set_parameter(elements_ID["Reseau_1600H"], "theta", theta)
        theta2 = classes.Parameter()
        get_parameter(elements_ID["Reseau_1600H"], "theta", theta2)
        print("New theta", theta2.value)

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
        if test_beamline:
            Cassiopee = classes.Beamline()
            Cassiopee.chains["linked_beamline"] = []
            for oe_name in linked_beamline:
                oe = classes.OpticalElement(element_id=elements_ID[oe_name])
                print(oe)
                Cassiopee.chains["linked_beamline"].append(oe)

    if test_add_element:
        elements_ID["screen"] = create_element("PlaneFilm", "screen")
        print("screen_id", elements_ID["screen"])
        get_element_name(elements_ID["screen"], elname)
        print("screen name", elname.value.decode())
        chain_element_by_id(elements_ID["f"], elements_ID["screen"])
    if test_spot_diagram:
        assert test_add_element
        set_recording(elements_ID["screen"], classes.RecordingMode.recording_output)

    if test_radiate:
        print(lamda_align)
        align(sourceID, lamda_align, show_return=True)
        clear_impacts(sourceID)
        generated_rays = generate(sourceID, lamda_align, show_return=True)
        radiate(sourceID, show_return=True)

    if test_spot_diagram:
        assert test_radiate

        # new_distance = Parameter()
        # get_parameter(elements_ID["pupille"], "distance", new_distance)
        # new_distance.value = DOUBLE(20-0.9)
        # set_parameter(elements_ID["pupille"], "distance", new_distance)
        # get_parameter(elements_ID["pupille"], "distance", new_distance)
        # print("New distance", new_distance.value)
        diagram = classes.Diagram(ndim=5, nreserved=int(generated_rays))
        get_spot_diagram(elements_ID["screen"], diagram, distance=0)
        print(np.ctypeslib.as_array(diagram.min, shape=(diagram.dim,)))
        print(diagram.count)
        import pandas as pd
        spots = pd.DataFrame(np.ctypeslib.as_array(diagram.spots, shape=(diagram.reserved, diagram.dim)),
                             columns=("X", "Y", "dX", "dY", "Lambda"))
        print(spots.head())
        figs = [ui_objects.plot_spd(spots, x_key="X", y_key="Y", light_plot=True, show_map=True),
                ui_objects.plot_spd(spots, x_key="X", y_key="dX", light_plot=False),
                ui_objects.plot_spd(spots, x_key="Y", y_key="dY", light_plot=True)]

        for fig in figs:
            ui_objects.show(fig)

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
