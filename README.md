Bibliothèque de simulation optique X pour le synchrotron SOLEIL
===============================================================

Install
-------

1. Installer miniconda
2. Créer un environnement dédié python 3.8.10 et y installer les paquets nécessaires

```bash
conda create PyOptiX python=3.8.10
conda activate PyOptiX
conda install lxml
```
3.1 Installer orange 

```bash
conda install orange3
```

Vérifier que l'installation s'est bien déroulée

```bash
python -m Orange.canvas
```

3.2 Installer Jupyter (environnement de script recommandé)

```bash
conda install jupyter
```

Créer le fichier de configuration de jupyter

```bash
jupyter notebook --generate-config
```
Cela crée un fichier dans `C:\Users\<username>\.jupyter\jupyter_notebook_config.py`

Trouver la ligne : `#c.NotebookApp.notebook_dir = ''`

Remplacer par `c.NotebookApp.notebook_dir = '/the/path/to/home/folder/'`

Ne pas oublier d'effacer le `#` en début de ligne

4. Télécharger la bibliothèque optix avec git :

```bash
git clone ssh://https://gitlab.synchrotron-soleil.fr/OPTIQUE/optical-simulation/optix
```

La compiler en version release avec msys2 cf. <https://www.msys2.org/> et codeblocks.

Télécharger la librairie libxml2-2 et linker libxml2.a dans le linker (sans chemin) et libxml/bin et libxml/lib dans les search directories du projet.
5. Télécharger PyOptiX

```bash
git clone ssh://gitlab.synchrotron-soleil.fr/OPTIQUE/optical-simulation/PyOptiX
```

6. Installer bokeh

```bash
conda install bokeh
```

7. Installer pandas

```bash
conda install pandas
```


