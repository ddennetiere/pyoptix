# coding: utf-8
from .exposed_functions import load_optix
from bokeh.io import output_notebook as bok_on

global optix
try:
    test = optix
    print(test, "already initialized")
except NameError:
    optix = load_optix()
    print("OptiX library initialized")


def output_notebook():
    bok_on()
