#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""Rearrange lon/lat to work with matplotlib pcolormesh.

fix_lon_lat_for_pcolormesh.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

# Standard Imports
import os
import sys

# Third-Party Imports
import numpy as np
import numpy.typing as npt


# Local Imports

##
# List of Public objects from this module.
__all__ = ["fix_lon_lat_for_pcolormesh"]

##
# Markup Language Specification (see Google Python Style Guide https://google.github.io/styleguide/pyguide.html)
__docformat__ = "Google en"
# ------------------------------------------------------------------------------

##
# Converts longitude to +/1 180
LON_180 = lambda x: ((x + 180.0) % 360.0) - 180.0

##
# Converts longitude to 0-360
LON_360 = lambda x: ((x + 360.0) % 360.0)

###############################################################################
# PUBLIC fix_lon_lat_for_pcolormesh()
# -----------------------------------
def fix_lon_lat_for_pcolormesh(in_lons: npt.ArrayLike, in_lats:npt.ArrayLike) -> tuple[npt.ArrayLike, npt.ArrayLike]:
    """https://bairdlangenbrunner.github.io/python-for-climate-scientists/matplotlib/pcolormesh-grid-fix.html"""
    # extend longitude by 2

    lon_extend = np.zeros(in_lons.size + 2)
    # fill in internal values
    lon_extend[1:-1] = in_lons  # fill up with original values
    # fill in extra endpoints
    lon_extend[0] = in_lons[0] - np.diff(in_lons)[0]
    lon_extend[-1] = in_lons[-1] + np.diff(in_lons)[-1]
    # calculate the midpoints
    lon_pcolormesh_midpoints = lon_extend[:-1] + 0.5 * (np.diff(lon_extend))
    # Same as left_lon_edge with cyclic longitude
    # in_lons                  (576): [-180.0000, -179.3750, ... +178.7500, +179.3750]
    # lon_pcolormesh_midpoints (577): [-180.3125, -179.6875, ... +179.0625, +179.6875]
    # print("\nin_lons                  ({}): [{:+9.4f}, {:+9.4f}, ... {:+9.4f}, {:+9.4f}]".format(len(in_lons), *in_lons[:2], *in_lons[-2:]))
    # print("lon_pcolormesh_midpoints ({}): [{:+9.4f}, {:+9.4f}, ... {:+9.4f}, {:+9.4f}]".format(len(lon_pcolormesh_midpoints), *lon_pcolormesh_midpoints[:2], *lon_pcolormesh_midpoints[-2:]))

    # extend latitude by 2
    lat_extend = np.zeros(in_lats.size + 2)
    # fill in internal values
    lat_extend[1:-1] = in_lats
    # fill in extra endpoints
    lat_extend[0] = in_lats[0] - np.diff(in_lats)[0]
    lat_extend[-1] = in_lats[-1] + np.diff(in_lats)[-1]
    #print("\nlat_extend", lat_extend.tolist())
    # calculate the midpoints
    lat_pcolormesh_midpoints = lat_extend[:-1] + 0.5 * (np.diff(lat_extend))
    # Don't extend beyond poles
    lat_pcolormesh_midpoints = np.where(lat_pcolormesh_midpoints > 90.0, 90.0, lat_pcolormesh_midpoints)
    lat_pcolormesh_midpoints = np.where(lat_pcolormesh_midpoints < -90.0, -90.0, lat_pcolormesh_midpoints)
    # in_lats                  (361): [ -90.0000.  -89.5000, ...  +89.5000,  +90.0000]
    # lat_pcolormesh_midpoints (362): [ -90.0000.  -89.7500, ...  +89.7500,  +90.0000]
    # print("\nin_lats                  ({}): [{:+9.4f}. {:+9.4f}, ... {:+9.4f}, {:+9.4f}]".format(len(in_lats), *in_lats[:2], *in_lats[-2:]))
    # print("lat_pcolormesh_midpoints ({}): [{:+9.4f}. {:+9.4f}, ... {:+9.4f}, {:+9.4f}]".format(len(lat_pcolormesh_midpoints), *lat_pcolormesh_midpoints[:2], *lat_pcolormesh_midpoints[-2:]))

    return lon_pcolormesh_midpoints, lat_pcolormesh_midpoints


# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
