#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""

cloudcast_test
~~~~~~~~~~~~~~

$ python cloudcast_test.py

"""
# Standard Imports
import os

# Third-Party Imports
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import BORDERS
from pyresample.geometry import AreaDefinition

# STARE Imports

# Local Imports

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

# Define Global Constants and State Variables
# -------------------------------------------

##
# Change matplotlib colors slightly (optional)
use_original = [False, True][0]
if use_original:
    blues = plt.get_cmap('Blues', 12)
    newcolors = blues(np.linspace(0, 1, 256))
    # Make colors less white
    newcolors[:, :] *= 0.9
    newcmp = matplotlib.colors.ListedColormap(newcolors)
else:
    newcmp = plt.get_cmap('tab20')
# os._exit(1)

##
# Define path to folder
path = ["/Volumes/val/data/CloudCast/full_cropped_cloud/", "/Volumes/val/data/CloudCast/small_cloud/"][0]

##
# Define file name
filename = [f"{path}2017M06.nc", f"{path}TestCloud.nc"][0]

# Load dataset using xarray
# -------------------------------------------
read_data = xr.open_dataarray(filename)

r"""
TestCloud.nc
    <xarray.DataArray (lat: 128, lon: 128, time: 17664)> Size: 289MB
        [289406976 values with dtype=uint8]
    Coordinates:
      * lat      (lat)  float64          1kB -8.461e+05 -8.281e+05 ... 1.422e+06 1.44e+06
      * lon      (lon)  float64          1kB -2.647e+06 -2.665e+06 ... -4.915e+06 -4.933e+06
      * time     (time) datetime64[ns] 141kB 2018-07-01            ... 2018-12-31T23:45:00

2017M06.nc
    <xarray.DataArray (lat: 768, lon: 768, time: 2879)> Size: 2GB
    [1698103296 values with dtype=uint8]
    Coordinates:
      * lat      (lat) float64          6kB -8.536e+05 -8.506e+05 ... 1.444e+06 1.447e+06
      * lon      (lon) float64          6kB -2.64e+06 -2.642e+06  ... -4.938e+06 -4.94e+06
      * time     (time) datetime64[ns] 23kB 2017-06-01T00:09:17   ... 2017-06-30T23...

"""
# print(read_data)
# os._exit(1)

# Pre-processing to match cloud types in paper
# -------------------------------------------

##
# Remove classes 1, 2, 3 and 4, which are cloud-free land, cloud-free sea, snow over land and sea ice.
read_data = read_data.where(read_data > 4)

##
# Subtract 4 to correspond to paper cloud types
read_data = read_data - 4

##
# Set nans to zero
read_data = read_data.fillna(0)

# Descriptive statistics on CloudCast Data
# -------------------------------------------

##
# Frequency of cloud types
r"""
Frequency distribution of various cloud types
    cloud_type  0:    0.31133%
    cloud_type  1:    0.12134%
    cloud_type  2:    0.13055%
    cloud_type  3:    0.10960%
    cloud_type  4:    0.10338%
    cloud_type  5:    0.00968%
    cloud_type  6:    0.07769%
    cloud_type  7:    0.03772%
    cloud_type  8:    0.04774%
    cloud_type  9:    0.02062%
    cloud_type 10:    0.03034%

Minimum            = 0.000
Maximum            = 10.000
Mean               = 2.774
Median             = 2.000
Standard Deviation = 2.870
"""
verbose = 0
if verbose:
    cloud_counts = {i:0 for i in range(0, 11)}
    total_count = read_data.shape[0] * read_data.shape[1] * read_data.shape[2]
    print("\nFrequency distribution of various cloud types")
    for cloud_type in cloud_counts.keys():
        cloud_counts[cloud_type] = ((read_data == cloud_type).sum() / total_count).item()
        print(f"\tcloud_type {cloud_type:2d}: {cloud_counts[cloud_type]:10.5f}%")
    print(f'\nMinimum            = {read_data.min().item():.3f}')
    print(f'Maximum            = {read_data.max().item():.3f}')
    print(f'Mean               = {read_data.mean().item():.3f}')
    print(f'Median             = {read_data.median().item():.3f}')
    print(f'Standard Deviation = {read_data.std().item():.3f}')
    os.exit(1)

# Visualizing CloudCast Data (Map)
# -------------------------------------------

##
# Plot a single time-step
read_data_single = read_data.isel(time = 40)

##
# For plotting purposes, we remove 0 (so we can actually see land)
read_data_single = read_data_single.where(read_data_single > 0)

##
# Set to False if you do not want background image
use_nasa_background = [False, True][1]

if use_nasa_background:
    # Use projection
    height = 768   # 3000 resolution
    width = 768
    lower_left_xy = [-855100.436345, -4942000.0]
    upper_right_xy = [1448899.563655, -2638000.0]
    area_def = AreaDefinition('areaD', 'Europe', 'areaD',
                              {'lat_0': '90.00', 'lat_ts': '50.00',
                               'lon_0': '5', 'proj': 'stere', 'ellps': 'WGS84'},
                              height, width,
                              (lower_left_xy[0], lower_left_xy[1],
                               upper_right_xy[0], upper_right_xy[1]))

    crs = area_def.to_cartopy_crs()
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=crs)
    ax.background_img(name='BM', resolution='low')
    a_image = plt.imshow(read_data_single.values, cmap=newcmp, transform=crs, extent=crs.bounds, origin='upper')
    fig.colorbar(a_image, ax=ax)


else:
    fig = plt.figure(figsize=(10, 8))
    map_proj = ccrs.PlateCarree()
    ax = plt.subplot(projection=map_proj)
    read_data_single.plot.imshow(ax=ax, cmap=newcmp, transform=map_proj, origin='upper')

##
# Show or save plot
plt.show()

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
