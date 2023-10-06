import os

from win32com.client import Dispatch
from setuptools import setup
from distutils import sysconfig
site_packages_path = sysconfig.get_python_lib()

path = os.path.abspath(os.path.curdir)
parentpath = os.path.abspath(os.path.pardir)
shortcut_filename = "PyOptiX launcher.lnk"

with open(os.path.join(path, "PyOptiX.bat"), "w") as batfile:
    batfile.write(f'start chrome --new-window "{path}\doc\html\index.html"\n')
    batfile.write(r'call C:\Miniconda3\Scripts\activate.bat pyoptix39'+'\n')
    batfile.write('start jupyter lab\n')

with open(os.path.join(site_packages_path, "pyoptix.pth"), "w") as pthfile:
    pthfile.write(f'{path}\n')
    pthfile.write(f'{parentpath}\n')
    pthfile.write(f'# Paths used for running pyoptix')

shell = Dispatch('WScript.Shell')
shortcut = shell.CreateShortCut(os.path.join(path, shortcut_filename))
shortcut.Targetpath = os.path.join(path, "PyOptiX.bat")
shortcut.WorkingDirectory = path
shortcut.IconLocation = os.path.join(path, "PyOptiX_icon.ico")
shortcut.save()

setup(
    name="pyoptix",
    version="1.0",
    url="https://gitlab.synchrotron-soleil.fr/OPTIQUE/optical-simulation/pyoptix",
    license="GPL v3",
    author="David Dennetiere",
    author_email="david.dennetiere@synchrotron-soleil.fr",
    description="Simulation of X-ray synchrotron beamlines",
    packages=['.'],
)
