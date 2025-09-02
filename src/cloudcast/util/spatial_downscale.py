#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
spatial_downscale.py
~~~~~~~~~~~~~~~~~~
"""
# Standard Imports
import os

# Third-Party Imports
import numpy as np
import numpy.typing as npt
from scipy.interpolate import griddata

# STARE Imports

# Local Imports

##
# List of Public objects from this module.
__all__ = ['spatial_downscale']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

###############################################################################
# PUBLIC spatial_downscale()
# ---------------------------
def spatial_downscale(indat, file_idx, spoints, flons, flats, spat_file, cheight, cwidth):
    ##
    # Spatial Interpolation
    use_method = ('nearest', 'linear', 'cubic')[2]
    fine_mesh_vals = griddata(spoints, indat.ravel(), (flons, flats), method=use_method)

    ##
    # Reshape
    fine_mesh_vals = np.reshape(fine_mesh_vals, (cheight, cwidth))

    ##
    # Save file
    s_file = spat_file.replace(".npz", f"_{file_idx:04d}.npz")
    np.savez_compressed(s_file, fine_mesh_vals)

    return s_file
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
