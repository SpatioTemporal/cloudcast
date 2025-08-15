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
from typing import Optional, Union, TypeAlias

# Third-Party Imports
import numpy as np
import pyresample as pr
import cartopy
import cartopy.crs as ccrs
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

# ##
# # Type Aliases
# SetOut: TypeAlias = tuple[int, int, int, typing.TypeVar('cartopy.crs'),
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')]]

###############################################################################
# PUBLIC setup_era5()
# ------------------
def setup_era5(a_file: str): #  -> SetOut:
    """
    Setup geometry for Meteosat Second Generation (MSG).

    Parameters
    ----------
    a_file : str
        Path to a ERA-5 to gather data

    Returns
    -------
    SetOut
        Tuple of geometry information for the specified projection and variable.
    """
    ##
    # MSG base projection WGS84 - World Geodetic System 1984
    #   https://epsg.io/4326
    #   +proj=longlat +datum=WGS84 +no_defs +type=crs
    MSG_EPSG = 32662
    msg_crs = ccrs.PlateCarree()

    # Lambert Conic Conformal (1SP) EPSG
    LCC_EPSG = 9801

    ##
    # Lambert Conic Conformal (2SP) EPSG
    LCC_EPSG = 9802

    ##
    # WGS-84 Earth equatorial radius at sea level (meters)
    #   STARE uses ccrs.Globe(datum='WGS84', ellipse='WGS84') https://epsg.io/4326
    #   Clarke 1866 ellipsoid https://epsg.io/7008-ellipsoid
    globe = ccrs.Globe(datum='WGS84', ellipse='WGS84')

    # Geodetic:
    #   A 3D/spherical CRS based on latitude and longitude where geographical distance and coordinates are measured in degrees.
    geod_crs = ccrs.Geodetic(globe=globe)

    verbose = [False, True][1]
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
    ncfile = netCDF4.Dataset(a_file, mode='r', format='NETCDF4_CLASSIC')

    era5_lats = ncfile.variables["latitude"][:]
    era5_lons = ncfile.variables["longitude"][:]
    nlats = len(era5_lats)
    nlons = len(era5_lons)
    dlat = dlon = 0.25
    ddlat = ddlon = 0.25 * 0.5

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

    ##
    # Read MSG lon and lats from file (saved in natread.py)
    lonlatfile = "/Users/mbauer/tmp/CloudCast/ccast_lonlat.npy"
    with open(lonlatfile, 'rb') as f:
        ccast_lons = np.load(f)
        ccast_lats = np.load(f)
    if verbose:
        tmp = ccast_lons.flatten()
        tmp_len = len(tmp)
        tmp = tmp[np.abs(tmp) <= 180.0]
        tmp1_len = len(tmp)
        len_frac = 100.0 * (tmp1_len / tmp_len)
        print(f"\n\tccast_lons {ccast_lons.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")

        tmp = ccast_lats.flatten()
        tmp_len = len(tmp)
        tmp = tmp[np.abs(tmp) <= 90.0]
        tmp1_len = len(tmp)
        len_frac = 100.0 * (tmp1_len / tmp_len)
        print(f"\tccast_lats {ccast_lats.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        del tmp

    from scipy.interpolate import griddata

    ##
    # The ERA-5 mesh to interpolate from (38801, 2)
    grid_lons, grid_lats = np.meshgrid(era5_lons, era5_lats)
    sparse_points = np.stack([grid_lons.ravel(), grid_lats.ravel()], -1)  # shape (N, 2) in 2d

    ##
    # The MSG fine mesh to interpolate into (589824, 2)
    # fine_points = np.stack([ccast_lons.ravel(), ccast_lats.ravel()], -1)  # shape (N, 2) in 2d
    ccast_nlons = len(ccast_lons)
    ccast_nlats = len(ccast_lats)

    ##
    # A simple ERA-5 Temperature Field
    era5_temp = np.squeeze(ncfile.variables["t"][:], axis=1)
    grid_centers_z_vals = era5_temp[0, :, :]
    del era5_temp
    print(f"\tERA-5 Temperature {grid_centers_z_vals.shape}: min {np.amin(grid_centers_z_vals):8.3f}, mean {np.mean(grid_centers_z_vals):8.3f}, max {np.amax(grid_centers_z_vals):8.3f} K")

    # z_dense_smooth_griddata = interp.griddata(sparse_points, z_sparse_smooth.ravel(), (x_dense, y_dense), method='cubic')

    # methods = ('nearest', 'linear', 'cubic')
    fine_mesh_vals = griddata(sparse_points, grid_centers_z_vals.ravel(), (ccast_lons.ravel(), ccast_lats.ravel()), method='cubic')
    # Reshape
    fine_mesh_vals = np.reshape(fine_mesh_vals, (ccast_nlats, ccast_nlons))
    print(f"\tInterpolated ERA-5 Temperature {fine_mesh_vals.shape}: min {np.amin(fine_mesh_vals):8.3f}, mean {np.mean(fine_mesh_vals):8.3f}, max {np.amax(fine_mesh_vals):8.3f} K")

    savefile = "/Users/mbauer/tmp/CloudCast/testera52msg.nc"
    ncfile = netCDF4.Dataset(savefile, mode='w', format='NETCDF4_CLASSIC')
    lat_dim = ncfile.createDimension('ny', ccast_nlats) # latitude axis
    lon_dim = ncfile.createDimension('nx', ccast_nlons) # longitude axis
    lat = ncfile.createVariable('lat', np.float32, ('ny', 'nx'))
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lon = ncfile.createVariable('lon', np.float32, ('ny', 'nx'))
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    outdata = ncfile.createVariable('t', np.int32, ('ny', 'nx'), fill_value=0)
    lat[:, :] = ccast_lats
    lon[:, :] = ccast_lons
    outdata[:, :] = fine_mesh_vals
    ncfile.close()

    # ##
    # # Create some information on the reference system
    # CCAST_HEIGHT = 768
    # CCAST_WIDTH = 768
    # lower_left_xy = [-855100.436345, -4942000.0]
    # upper_right_xy = [1448899.563655, -2638000.0]

    # # Define the area
    # #   <class 'pyresample.geometry.AreaDefinition'>
    # ccast_area_def = pr.geometry.AreaDefinition('areaD', 'Europe', 'areaD',
    #                                             {'lat_0': '90.00', 'lat_ts': '50.00', 'lon_0': '5', 'proj': 'stere', 'ellps': 'WGS84'},
    #                                             CCAST_HEIGHT, CCAST_WIDTH,
    #                                             (lower_left_xy[0], lower_left_xy[1], upper_right_xy[0], upper_right_xy[1]))
    # #   +proj=stere +lat_0=90 +lat_ts=50 +lon_0=5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs +type=crs
    # # print(ccast_area_def.proj4_string)

    # ##
    # # Form a cartopy CRS
    # #   <class 'pyresample.utils.cartopy.Projection'>
    # ccast_crs = ccast_area_def.to_cartopy_crs()
    # # print(ccast_merc_crs.bounds)
    # #   ccast_crs.bounds: (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)


    # dtype = "float32"
    # radius = 16000
    # epsilon = 0.5
    # nodata = -3.4E+38
    # ##
    # # Apply a swath definition for our output raster
    # #   <class 'pyresample.geometry.SwathDefinition'>
    # exstract_def = pr.geometry.SwathDefinition(lons=era5_lons, lats=era5_lats)

    # ##
    # # Resample our data to the area of interest
    # ccast_data_vals = pr.kd_tree.resample_nearest(exstract_def, data_vals,
    #                                               ccast_area_def,
    #                                               radius_of_influence=radius, # in meters
    #                                               epsilon=epsilon,
    #                                               fill_value=False)
    # ccast_lons = pr.kd_tree.resample_nearest(exstract_def, lons,
    #                                          ccast_area_def,
    #                                          radius_of_influence=radius, # in meters
    #                                          epsilon=epsilon,
    #                                          fill_value=False)
    # ccast_lats = pr.kd_tree.resample_nearest(exstract_def, lats,
    #                                          ccast_area_def,
    #                                          radius_of_influence=radius, # in meters
    #                                          epsilon=epsilon,
    #                                          fill_value=False)


# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
