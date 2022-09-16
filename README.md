Bibliothèque de simulation optique X pour le synchrotron SOLEIL
===============================================================
This library codes for a simulation tool for designing X-ray beamline. It has specifically
been designed for the SOLEIL Synchrotron, but can be used for any X-ray optical design as shown in examples.
 

Installing
----------

1. Install miniconda

   Installation page pf miniconda: https://docs.conda.io/en/latest/miniconda.html

2. Create a python 3.8.10 dedicated environment and install required packages

   ```bash
   conda create PyOptiX python=3.8.10
   conda activate PyOptiX
   conda install lxml
   ```

3. Install Jupyterlab (recommended scripting environment )

   ```bash
   conda install jupyterlab bokeh pandas
   ```

   Create jupyter configuration file

   ```bash
   jupyter notebook --generate-config
   ```

   There is now a file at `C:\Users\<username>\.jupyter\jupyter_notebook_config.py`

   Find the line : `#c.NotebookApp.notebook_dir = ''`

   Replace by `c.NotebookApp.notebook_dir = '/the/path/to/home/folder/'`

   Do not forget to delete the `#` at line start

4. Download optix using git :

   ```bash
   git clone ssh://https://gitlab.synchrotron-soleil.fr/OPTIQUE/optical-simulation/optix
   ```

   Compile its release version using msys2 cf. https://www.msys2.org/ and codeblocks.

   Download libxml2-2 library and link libxml2.a in the codeblocks linker (without its path) and libxml/bin and libxml/lib in the search directories of the project.

5. Download PyOptiX

   ```bash
   git clone ssh://gitlab.synchrotron-soleil.fr/OPTIQUE/optical-simulation/PyOptiX
   ```

6. Install ipysheet

   ```bash
   conda install -c conda-forge ipysheet
   pip install -U "nbclassic>=0.2.8"
   jupyter nbextension enable --py --sys-prefix widgetsnbextension
   ```

Examples
--------
Run README.ipynb in jupyterlab to access the following example notebooks

### Design and simulation

#### Lignes SOLEIL
- **[Herm1](Hermes%20notebook.ipynb)** - [Hermes beamline](Hermes%20notebook.ipynb)  
Definition and simulation of the HERMES beamline at SOLEIL
- **[Herm2](HermesUP%20notebook.ipynb)** - [Hermes beamline](HermesUP%20notebook.ipynb)  
Definition and simulation of the upgraded HERMES beamline at SOLEIL
- **[Deim1](Deimos%20notebook.ipynb)** - [Deimos beamline](Deimos%20notebook.ipynb)  
Definition and simulation of the DEIMOS beamline at SOLEIL
- **[Deim2](DeimosUP%20notebook.ipynb)** - [Deimos beamline](DeimosUP%20notebook.ipynb)  
Definition and simulation of the upgraded DEIMOS beamline at SOLEIL
- **[Anta1](Antares%20notebook.ipynb)** - [Antares beamline](Antares%20notebook.ipynb)  
Definition and simulation of the Antares beamline at SOLEIL
- **[Anta2](AntaresUP%20notebook.ipynb)** - [Antares beamline](AntaresUP%20notebook.ipynb)  
Definition and simulation of the upgraded Antares beamline at SOLEIL
- **[Cass1](Cassiopee%20notebook.ipynb)** - [Cassiopee beamline](Cassiopee%20notebook.ipynb)  
Definition and simulation of the Cassiopee beamline at SOLEIL
- **[Cass2](CassiopeeUP%20notebook.ipynb)** - [Cassiopee beamline](CassiopeeUP%20notebook.ipynb)  
Definition and simulation of the upgraded Cassiopee beamline at SOLEIL


#### Lignes hors SOLEIL
- **[Atto1](Ellipse%20Helimag.ipynb)** - [Helimag beamline V1](Ellipse%20Helimag.ipynb)  
Definition and simulation of a beamline for the Attolab facility using revolution conics mirrors
- **[Atto2](Wolter%20Helimag.ipynb)** - [Helimag beamline V2](Wolter%20Helimag.ipynb)  
Definition and simulation of a beamline for the Attolab facility using toroidal mirrors 
- **[Atto3](KB_HELIMAG.ipynb)** - [Helimag beamline V3](KB_HELIMAG.ipynb)  
Definition and simulation of a beamline for the Attolab facility using conical mirrors in a Kirkpatrick-Baez configuration 


### Tolerancing systems
- **[Atto3b](KB_HELIMAG_tolerancing.ipynb)** - [Helimag beamline V3](KB_HELIMAG_tolerancing.ipynb)  
	Import of a beamline for the Attolab facility using conical mirrors in a Kirkpatrick-Baez configuration defined in
	**[Atto3](KB_HELIMAG.ipynb)** et tolerancing of the beamline

### Optimization
- **[HERM1](Hermes%20notebook.ipynb)** - [Hermes beamline](Hermes%20notebook.ipynb)  
Definition and simulation of the HERMES beamline at SOLEIL
