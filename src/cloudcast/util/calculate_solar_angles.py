#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Find solar and viewing angles for Meteosat Second Generation (MSG).

Converted from matlab code found in https://github.com/222huan/Calculate_MSG_angles

calculate_solar_angles
~~~~~~~~~~~~~~~~~~~~~~
"""
# Standard Imports
import os
import typing
from typing import Optional, Union, TypeAlias

# Third-Party Imports
import numpy as np
import numpy.typing as npt

# STARE Imports

# Local Imports
#   conda install -c conda-forge pyorbital
##
# List of Public objects from this module.
__all__ = ['calculate_solar_angles']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

###############################################################################
# PUBLIC calculate_solar_angles()
# -------------------------------
def calculate_solar_angles():
    # function [SZA, SAA] = Calculate_solar_angles(lat,lon,doy,hour,minute,timezone)
    # % Reference: "General Solar Position Calculations" by NOAA Global Monitoring Division
    # % Reference: https://en.wikipedia.org/wiki/Solar_azimuth_angle

    # % calculate the solar zenith angle SZA and solar azimuth angle SAA (in degrees)
    # % the hour and minute is in local time system (if in UTC system, timezone=0)
    # % the hour and minute can be a list


    return

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
