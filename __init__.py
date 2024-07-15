# coding: utf-8
from .exposed_functions import load_optix, set_aperture_active
from bokeh.io import output_notebook as bok_on
from .classes import save_beamline, load_beamline
import logging

global optix
logger = logging.getLogger(__name__)
try:
    test = optix
    print(test, "already initialized")
    logging.basicConfig(filename="pyoptix.log", level=logging.DEBUG)
    logger.info("Logging started")
except NameError:
    optix = load_optix()
    print("OptiX library initialized")


def output_notebook():
    bok_on()
