import os
from win32com.client import Dispatch

path = os.path.abspath(os.path.curdir)
shortcut_filename = "PyOptiX launcher.lnk"

with open(os.path.join(path, "PyOptiX.bat"), "w") as batfile:
    batfile.write(f'start chrome --new-window "{path}\doc\html\index.html"\n')
    batfile.write(r'call C:\Miniconda3\Scripts\activate.bat pyoptix'+'\n')
    batfile.write('start jupyter notebook\n')

shell = Dispatch('WScript.Shell')
shortcut = shell.CreateShortCut(os.path.join(path, shortcut_filename))
shortcut.Targetpath = os.path.join(path, "PyOptiX.bat")
shortcut.WorkingDirectory = path
shortcut.IconLocation = os.path.join(path, "PyOptiX_icon.ico")
shortcut.save()
