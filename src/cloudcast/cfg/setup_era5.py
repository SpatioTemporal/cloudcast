#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Setup geometry for Read ERA-5.

setup_era5.py
~~~~~~~~~~~~
"""
# Standard Imports
import os
import typing

# Third-Party Imports
import numpy as np
import numpy.ma as ma
import numpy.typing as npt
import netCDF4

# STARE Imports

# Local Imports

##
# List of Public objects from this module.
__all__ = ['setup_era5']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

def LON_TO_180(x): return ((x + 180.0) % 360.0) - 180.0
def LON_TO_360(x): return (x + 360.0) % 360.0

###############################################################################
# PUBLIC setup_era5()
# ------------------
def setup_era5(a_file: str) -> tuple[list[float], list[float], list[float], list[float]]:
    """
    Setup geometry for Meteosat Second Generation (MSG).

    Parameters
    ----------
    a_file : str
        Path to a ERA-5 to gather data

    # ('numpy.ma.MaskedArray', 'numpy.ma.MaskedArray', list[float], list)

    Returns
    -------
    era5_lons (161): List of ERA-5 grid center longitudes.
    era5_lats (241): List of ERA-5 grid center latitudes.
    era5_lon_edges (162): List of ERA-5 grid edge longitudes.
    era5_lat_edges (242): List of ERA-5 grid edge latitudes.
    """

    verbose = [False, True][0]
    r"""
    ERA-5 Latitudes  161: [ +70.000, ...  +30.000]
    ERA-5 Longitudes 241: [ -20.000, ...  +40.000]

    Grid Spacing 0.25 x 0.25 deg
                 ~30 x 30 km

    Grid Centers
        (70, -20)----------------(70, 40)
        |                               |
        |                               |
        |                               |
        |                               |
        (30, -20)----------------(30, 40)

    Grid Edges
        (70.125, -20.125)------+------(70.125, -19.875)-----------------(70.125, 39.875)------+------(70.125, 40.125)
                            (70, -20)                                                      (70, 40)
        (69.875, -20.125)------+------(69.875, -19.875)                 (69.875, 39.875)------+------(69.875, 40.125)
        |                                                                                                           |
        |                                                                                                           |
        |                                                                                                           |
        |                                                                                                           |
        |                                                                                                           |
        (30.125, -20.125)------+------(30.125, -19.875)                 (30.125, 39.875)------+------(30.125, 40.125)
                            (30, -20)                                                      (30, 40)
        (29.875, -20.125)------+------(29.875, -19.875)-----------------(29.875, 39.875)------+------(29.875, 40.125)

    For interpolation we will project the data to put it in a linear Cartesian framework.

    Domain based on Grid Edges
        UL_edge(+70.125, -20.125)------------------UR_edge(+70.125, +40.125)
        |                                                                  |
        LL_edge(+29.875, -20.125)------------------LR_edge(+29.875, +40.125)

    Domain based on Grid Centers
        UL(+70.000, -20.000)------------------UR(+70.000, +40.000)
        |                                                        |
        LL(+30.000, -20.000)------------------LR(+30.000, +40.000)

    The MSG CloudCast Domain is (ccast_crs)
        UL(+62.403, -12.921)------------------UR(+60.151, +33.739)
        |
        LL(+42.068,  -4.803)------------------LR(+40.928, +21.329)

    So the MSG domain fits easily within the bounds of the ERA-5 source data.
    """

    ##
    # Read an example NETCDF file for ERA-5
    ncfile = netCDF4.Dataset(a_file, mode='r', format='NETCDF4_CLASSIC')

    ##
    # Get grid center lon/lat
    era5_lats = ncfile.variables["latitude"][:]
    era5_lons = ncfile.variables["longitude"][:]
    ncfile.close()

    nlats = len(era5_lats)
    nlons = len(era5_lons)
    dlat = dlon = 0.25
    ddlat = ddlon = 0.25 * 0.5

    ##
    # Remove Mask
    era5_lons = era5_lons.tolist(fill_value=-999)
    era5_lats = era5_lats.tolist(fill_value=-999)

    ##
    # Find Grid edges
    era5_lat_edges = dict([(round(float(i) + ddlat, 3), 1) for i in era5_lats])
    era5_lat_edges.update([(round(float(i) - ddlat, 3), 1) for i in era5_lats])
    era5_lat_edges = sorted(list(era5_lat_edges.keys()), reverse=True)

    era5_lon_edges = dict([(float(i) + ddlon, 1) for i in era5_lons])
    era5_lon_edges.update([(float(i) - ddlon, 1) for i in era5_lons])
    era5_lon_edges = sorted(list(era5_lon_edges.keys()))

    # Corners of bounding box
    ll_x_edge = era5_lon_edges[0]
    ll_y_edge = era5_lat_edges[-1]

    ul_x_edge = era5_lon_edges[0]
    ul_y_edge = era5_lat_edges[0]

    ur_x_edge = era5_lon_edges[-1]
    ur_y_edge = era5_lat_edges[0]

    lr_x_edge = era5_lon_edges[-1]
    lr_y_edge = era5_lat_edges[-1]

    ll_x = era5_lons[0]
    ll_y = era5_lats[-1]
    ul_x = era5_lons[0]
    ul_y = era5_lats[0]
    ur_x = era5_lons[-1]
    ur_y = era5_lats[0]
    lr_x = era5_lons[-1]
    lr_y = era5_lats[-1]

    if verbose:
        print(f"\tERA-5 Latitudes                           ({nlats}): [{era5_lats[0]:+8.3f}, ... {era5_lats[-1]:+8.3f}]")
        print(f"\tERA-5 Longitudes                          ({len(era5_lons)}): [{era5_lons[0]:+8.3f}, ... {era5_lons[-1]:+8.3f}]")
        print("\n")
        print(f"\tUL_edge({ul_y_edge:+7.3f}, {ul_x_edge:+7.3f})------------------UR_edge({ur_y_edge:+7.3f}, {ur_x_edge:+7.3f})")
        print("\t|                                                                  |")
        print(f"\tLL_edge({ll_y_edge:+7.3f}, {ll_x_edge:+7.3f})------------------LR_edge({lr_y_edge:+7.3f}, {lr_x_edge:+7.3f})")
        print("\n")
        print(f"\tUL({ul_y:+7.3f}, {ul_x:+7.3f})------------------UR({ur_y:+7.3f}, {ur_x:+7.3f})")
        print("\t|                                                        |")
        print(f"\tLL({ll_y:+7.3f}, {ll_x:+7.3f})------------------LR({lr_y:+7.3f}, {lr_x:+7.3f})")

    return (era5_lons, era5_lats, era5_lon_edges, era5_lat_edges)
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
