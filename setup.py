#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""

setup.py
~~~~~~~~

Install the pod_pf package
  $ cd ~/myDev/cloudcast

  Mac
      $ conda env export > environment.yml
      $ /Users/mbauer/miniconda3/envs/stare/bin/python -m pip install -e .
      $ pip install --editable .

  Linux
      $ pip freeze > requirements.txt
      # Edit to remove references to things like pystare
      $ pip install --editable .

$  pip list
      Package                           Version           Editable project location
      --------------------------------- ----------------- ----------------------------
      cloudcast                     0.1               /Users/mbauer/myDev/cloudcast/src

"""
# Standard Imports
from glob import glob
from os.path import basename, splitext
from setuptools import find_packages, setup

##
# Markup Language Specification (see Google Python Style Guide https://google.github.io/styleguide/pyguide.html)
__docformat__ = "Google en"
# ------------------------------------------------------------------------------

setup(name='cloudcast',
      description='Use STARE to POD IMERG precipitation features (PF)',
      author='Mike Bauer',
      author_email='mbauer@bayesics.com',
      version='0.1',
      package_dir = {'': 'src'},
      packages = ['cloudcast', 'cloudcast.cfg', 'cloudcast.util', 'cloudcast.plot'])
      # packages = find_packages('src', exclude=include_mcms),
      # py_modules = [splitext(basename(path))[0] for path in glob('src/*.py')])


# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
