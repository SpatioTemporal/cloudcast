#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Read ERA-5 File

read_era5.py
~~~~~~~~~~~~~
"""
# Standard Imports
import os
import typing
from typing import Optional, Union

# Third-Party Imports
import numpy as np
import numpy.typing as npt

# STARE Imports

# Local Imports

##
# List of Public objects from this module.
__all__ = ['read_era5']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

###############################################################################
# PUBLIC read_era5()
# ------------------
def read_era5():


    """
    ccast_crs
        ccast_area_def.corners = [
            (-12.920934886492649, 62.403066090517555), UL
            (33.73865749382469,   60.15059617915855),  UR
            (21.32880156090482,   40.92817004367345),  LR
            (-4.802566482888071,  42.068097533886025)] LL

        ccast_area_def.area_extent  (-855100.436345, -4942000.0, 1448899.563655, -2638000.0)
        ccast_crs.bounds      (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)

        lower_left_xy  = [-855100.436345, -4942000.0] => (-4.816534709314307, 42.053336570266744)
        upper_right_xy = [1448899.563655, -2638000.0] => (33.77742545811928,  60.15622631725923)
    """

    return

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
