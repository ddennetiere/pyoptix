Bibliothèque de simulation optique X pour le synchrotron SOLEIL
===============================================================
This library codes for a simulation tool for designing X-ray beamline. It has specifically
been designed for the SOLEIL Synchrotron, but can be used for any X-ray optical design as shown in examples.
 

Installing
----------

1. Install miniconda

   Installation page pf miniconda: https://docs.conda.io/en/latest/miniconda.html

2. Create a python 3.9 dedicated environment and install required packages

   ```bash
   conda env create -f pyoptix39.yml
   conda activate pyoptix39
   ```

   Create jupyter configuration file

   ```bash
   jupyter notebook --generate-config
   ```

   There is now a file at `C:\Users\<username>\.jupyter\jupyter_notebook_config.py`

   Find the line : `#c.NotebookApp.notebook_dir = ''`

   Replace by `c.NotebookApp.notebook_dir = '/the/path/to/home/folder/'`

   Do not forget to delete the `#` at line start

3. Download optix using git :

   ```bash
   git clone https://github.com/ddennetiere/optix
   ```

   Compile its release version using msys2 cf. https://www.msys2.org/ and codeblocks.

   Follow optix README

Examples
--------
Run README.ipynb in jupyterlab to access the following example notebooks

### Design and simulation

#### Beamline from SOLEIL
**[Herm1](Tests/Hermes%20notebook.ipynb)** - [Hermes beamline](Tests/Hermes%20notebook.ipynb)  
Definition and simulation of the HERMES beamline at SOLEIL
This file shows a complete simulation with chaptering for design, optimisation, visualization, etc.