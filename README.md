# coding: utf-8
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
3. Installer orange 

```bash
conda install orange3
```

Vérifier que l'installation s'est bien déroulée

```bash
python -m Orange.canvas
```
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



