# -*- coding: utf-8 -*-
import ctypes, os
# from ctypes import cdll, windll
from ctypes.wintypes import *
from ctypes import *
from ctypes import WINFUNCTYPE, Structure, pointer, byref, POINTER,  c_char_p, c_void_p, c_double, c_float, \
                   create_string_buffer
import numpy as np
from classes import Beamline, OpticalElement, Parameter
from lxml import etree

global optix


def load_optix():
    global optix
    print("intialzing SR library")
    optix = cdll.LoadLibrary(r'D:\Dennetiere\optix\release\OptiX.dll')
    # optix.LoadSolemioFile.restype = BYTE  # gcc bool match to ctypes BYTE
    # optix.LoadSolemioFile.argtypes = [LPCSTR]
    optix.GetOptiXLastError.restype = BYTE
    optix.GetOptiXLastError.argtypes = [ctypes.c_char_p, INT]


def release_optix():
    global optix
    windll.kernel32.FreeLibrary.argtypes = [HMODULE]
    libhandle = optix._handle
    windll.kernel32.FreeLibrary(libhandle)
    # del SRlib


def version():
    optix.Version()


def load_solemio_file(name):
    if optix.LoadSolemioFile(name) == 0:
        buf = ctypes.create_string_buffer(256)  # create a 128 byte buffer
        optix.GetOptiXLastError(buf, 256)
        print("error loading Solemio file :")
        print(buf.value)


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


if __name__ == "__main__":
    # initialisation auto
    try:
        test = optix
        print(test, "  already initialized")
    except NameError:
        load_optix()
        print("OptiX library initialized")
    # parse_xml(r"D:\Dennetiere\optix\bin\test\system.xml")
    # optix.LoadSolemioFile.argtypes = [c_char_p]
    load_solemio_file(create_string_buffer(b"D:\\Dennetiere\\Programmes Python\\optix\\solemio\\CASSIOPEE"))
    hsys, hparam, elemID = c_int64(0), c_int64(0), c_int64(0)
    elname = create_string_buffer(32, c_char)
    # elemname2 = create_string_buffer(32)
    param_name = create_string_buffer(48)
    param = Parameter()
    print("#"*80)
    optix.EnumerateElements(byref(hsys), byref(elemID), elname, 32)
    while hsys:
        print(hsys)
        optix.EnumerateElements(byref(hsys), byref(elemID), elname, 32)
        # optix.GetElementName(elemID, elemname2, 32)
        print("-"*20)
        print(elname.value)
        if not optix.EnumerateParameters(elemID, byref(hparam), param_name, 48, param):
            buf = ctypes.create_string_buffer(256)  # create a 128 byte buffer
            optix.GetOptiXLastError(buf, 256)
            print("error", buf.value)
        while hparam:
            if not optix.EnumerateParameters(elemID, byref(hparam), param_name, 48, param):
                buf = ctypes.create_string_buffer(256)  # create a 128 byte buffer
                optix.GetOptiXLastError(buf, 256)
                print("error", buf.value)
            else:
                print(param_name.value, param.value, param.bounds.min, param.bounds.max, param.multiplier, param.type,
                      param.group, param.flags)

# system = optix.LoadSolemioFile("D:\\Dennetiere\\Programmes Python\\optix\\solemio\\CASSIOPEE")
# # system_xml = optix.LoadSystemFromXml(r"D:\Dennetiere\optix\bin\test\system.xml")
# # print("system", system_xml)
# hsys, hparam, elemID = c_int64(0), c_int64(0), c_int64(0)
# print(byref(hsys))
# # elname, name2, parmname, errbuf = c_buffer(32*4), c_buffer(32), c_buffer(48), c_buffer(256)
# elname = create_string_buffer(32, c_char)
# # elname_size = c_char * 32
# # elname_buf = elname_size()
# optix.EnumerateElements.argtypes = [c_void_p, c_void_p, c_char_p, c_int]
# optix.EnumerateElements(byref(hsys), byref(elemID), elname, 32)
# print(hsys, elemID, elname)
# optix.Align()
# # optix.EnumerateElements()
