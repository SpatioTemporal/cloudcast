#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Read Meteosat Second Generation (MSG) Native Archive Format (.nat) file.

natread.py
~~~~~~~~~~
"""
# Standard Imports
import os
import copy
import typing
from typing import Optional, Union, TypeAlias

# Third-Party Imports
import numpy as np
from osgeo import gdal
from osgeo import osr
import pyresample as pr
#   conda install -c conda-forge satpy
from satpy import Scene
import cartopy
import cartopy.crs as ccrs

# STARE Imports

# Local Imports
from cloudcast.util.readnat import readnat

##
# List of Public objects from this module.
__all__ = ['natread']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------


###############################################################################
# PUBLIC natread()
# ----------------
def natread(fname: str, fvar: str, reader: str, to_euro: bool) -> None:

    verbose = [False, True][1]

    ##
    # Read the file
    #   scn = <class 'satpy.scene.Scene'>
    r"""
    scn.all_dataset_names() =
        ['HRV', 'IR_016', 'IR_039', 'IR_087', 'IR_097', 'IR_108', 'IR_120', 'IR_134', 'VIS006', 'VIS008', 'WV_062', 'WV_073']
    scn.available_dataset_ids() =
        [DataID(name='HRV', wavelength=WavelengthRange(min=0.5, central=0.7, max=0.9, unit='µm'), resolution=1000.134348869, calibration=<calibration.reflectance>, modifiers=()),
         DataID(name='HRV', wavelength=WavelengthRange(min=0.5, central=0.7, max=0.9, unit='µm'), resolution=1000.134348869, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='HRV', wavelength=WavelengthRange(min=0.5, central=0.7, max=0.9, unit='µm'), resolution=1000.134348869, calibration=<calibration.counts>, modifiers=()),
         DataID(name='IR_016', wavelength=WavelengthRange(min=1.5, central=1.64, max=1.78, unit='µm'), resolution=3000.403165817, calibration=<calibration.reflectance>, modifiers=()),
         DataID(name='IR_016', wavelength=WavelengthRange(min=1.5, central=1.64, max=1.78, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='IR_016', wavelength=WavelengthRange(min=1.5, central=1.64, max=1.78, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='IR_039', wavelength=WavelengthRange(min=3.48, central=3.92, max=4.36, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_039', wavelength=WavelengthRange(min=3.48, central=3.92, max=4.36, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='IR_039', wavelength=WavelengthRange(min=3.48, central=3.92, max=4.36, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='IR_087', wavelength=WavelengthRange(min=8.3, central=8.7, max=9.1, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_087', wavelength=WavelengthRange(min=8.3, central=8.7, max=9.1, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='IR_087', wavelength=WavelengthRange(min=8.3, central=8.7, max=9.1, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='IR_097', wavelength=WavelengthRange(min=9.38, central=9.66, max=9.94, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_097', wavelength=WavelengthRange(min=9.38, central=9.66, max=9.94, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='IR_097', wavelength=WavelengthRange(min=9.38, central=9.66, max=9.94, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='IR_108', wavelength=WavelengthRange(min=9.8, central=10.8, max=11.8, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_108', wavelength=WavelengthRange(min=9.8, central=10.8, max=11.8, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='IR_108', wavelength=WavelengthRange(min=9.8, central=10.8, max=11.8, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='IR_120', wavelength=WavelengthRange(min=11.0, central=12.0, max=13.0, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_120', wavelength=WavelengthRange(min=11.0, central=12.0, max=13.0, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='IR_120', wavelength=WavelengthRange(min=11.0, central=12.0, max=13.0, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='IR_134', wavelength=WavelengthRange(min=12.4, central=13.4, max=14.4, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_134', wavelength=WavelengthRange(min=12.4, central=13.4, max=14.4, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='IR_134', wavelength=WavelengthRange(min=12.4, central=13.4, max=14.4, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='VIS006', wavelength=WavelengthRange(min=0.56, central=0.635, max=0.71, unit='µm'), resolution=3000.403165817, calibration=<calibration.reflectance>, modifiers=()),
         DataID(name='VIS006', wavelength=WavelengthRange(min=0.56, central=0.635, max=0.71, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='VIS006', wavelength=WavelengthRange(min=0.56, central=0.635, max=0.71, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='VIS008', wavelength=WavelengthRange(min=0.74, central=0.81, max=0.88, unit='µm'), resolution=3000.403165817, calibration=<calibration.reflectance>, modifiers=()),
         DataID(name='VIS008', wavelength=WavelengthRange(min=0.74, central=0.81, max=0.88, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='VIS008', wavelength=WavelengthRange(min=0.74, central=0.81, max=0.88, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='WV_062', wavelength=WavelengthRange(min=5.35, central=6.25, max=7.15, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='WV_062', wavelength=WavelengthRange(min=5.35, central=6.25, max=7.15, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='WV_062', wavelength=WavelengthRange(min=5.35, central=6.25, max=7.15, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),
         DataID(name='WV_073', wavelength=WavelengthRange(min=6.85, central=7.35, max=7.85, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='WV_073', wavelength=WavelengthRange(min=6.85, central=7.35, max=7.85, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='WV_073', wavelength=WavelengthRange(min=6.85, central=7.35, max=7.85, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=())
        ]
    scn.coarsest_area =
        <bound method Scene.coarsest_area of <satpy.scene.Scene object at 0x1172981f0>>
    scn.images() =
        <generator object Scene.images at 0x14849f920>
    """
    scn = Scene(filenames = {reader:[fname]})
    # print(scn.values())

    # ##
    # # Extract data set names
    # dataset_names = scn.all_dataset_names()
    # """
    #     HRV
    #     IR_016
    #     IR_039
    #     IR_087
    #     IR_097
    #     IR_108
    #     IR_120
    #     IR_134
    #     VIS006
    #     VIS008
    #     WV_062
    #     WV_073
    # """

    ##
    # Extract wanted data
    """
    if fvar == HRV
        lons (11136, 5568)     : [-65.47019958496094, ... 81.21913146972656]
        lats (11136, 5568)     : [-81.13614654541016, ... 81.13614654541016]
        data_vals (11136, 5568): [0.0, ... 26.103036880493164]

        good_lons (50104766): [-57.18206024169922 ... 79.4859390258789]
        good_lats (50104766): [-78.96027374267578 ... 78.3008804321289]

    else:
        lons (3712, 3712)      : [-81.12566375732422, ... 81.12566375732422]
        lats (3712, 3712)      : [-81.0744857788086, ... 81.0744857788086]
        data_vals (3712, 3712) : [0.0, ... 3.5562241077423096]

        good_lons (10213685): [-75.26545715332031 ... 75.56462097167969]
        good_lats (10213685): [-78.95874786376953 ... 78.29975891113281]
    """
    lons, lats, data_vals = readnat(file=fname, calibration="radiance", dataset=fvar, reader=reader, dtype="float32")
    if verbose:
        tmp = lons.flatten()
        tmp = tmp[np.abs(tmp) <= 180.0]
        print(f"\tlons {lons.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        tmp = lats.flatten()
        tmp = tmp[np.abs(tmp) <= 90.0]
        print(f"\tlats {lats.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        tmp = data_vals.flatten()
        tmp = tmp[~np.isnan(tmp)]
        print(f"\tdata_vals {data_vals.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        del tmp
    if to_euro:
        # # lats (928, 3712) (1932430,): [26.656396865844727, ... 81.0744857788086]
        # # lats = lats[3712 - euro_nrows:, :]

        # # row (jj) starts at S and moves N
        # # col (ii) starts in E and moves W
        # # lons[jj, ii] where jj is row and ii is column
        # # 45W - 30E

        # # # jj = 2652 ii = 3146: (lat, lon) (23.81135368347168, -45.049129486083984)
        # # jj = 2652
        # # ii = 3146
        # # print(f"{jj = :4d} {ii = :4d}: (lat, lon) ({lats[jj, ii]}, {lons[jj, ii]})")

        # # # jj = 2652 ii =  910: (lat, lon) (23.114727020263672, 30.13174819946289)
        # # jj = 2652
        # # ii = 910
        # # print(f"{jj = :4d} {ii = :4d}: (lat, lon) ({lats[jj, ii]}, {lons[jj, ii]})")

        # # jj = 2652 ii =  910: (lat, lon) (23.114727020263672, 30.13174819946289)
        # jj = 2652
        # ii = 3146 - euro_ncols
        # print(f"{jj = :4d} {ii = :4d}: (lat, lon) ({lats[jj, ii]}, {lons[jj, ii]})")

        # 3146 - 1530 = 1616
        # # euro_nrows = 928
        # # euro_ncols = 1530
        # # 3146-910 = 2236
        # # 2236 - 1530 = 706
        # os._exit(1)

        # dones = 0
        # for jj in range(3712):
        #     for ii in range(3712):
        #         if np.abs(lons[jj, ii]) <= 180.0 and np.abs(lats[jj, ii] <= 90.0):
        #             if lats[jj, ii] >= 26.0 and lons[jj, ii] <= -45:
        #                 print(f"{jj = :4d} {ii = :4d}: (lat, lon) ({lats[jj, ii]}, {lons[jj, ii]})")
        #                 dones += 1
        #     # if dones > 5:
        #     #     os._exit(1)

        lats = lats[3712 - euro_nrows:, :]
        lons = lons[3712 - euro_nrows:, :]

        # 45W - 30E
        # 3712-1530
        if verbose:
            tmp = lons.flatten()
            tmp = tmp[np.abs(tmp) <= 180.0]
            print(f"\n\tlons {lons.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            tmp = lats.flatten()
            tmp = tmp[np.abs(tmp) <= 90.0]
            print(f"\tlats {lats.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            del tmp

    # ny, nx = lons.shape
    # # for yidx in range(53, 60):
    # for yidx in range(ny):
    #     left_good_idx = -1
    #     right_good_idx = -1
    #     ## print(f"Checking {yidx = :5d}")
    #     for xidx in range(nx):
    #         ## print(f"\tChecking {xidx = :5d} {data_vals[yidx, xidx]}")
    #         if np.isnan(data_vals[yidx, xidx]):
    #             # Invalid data value index
    #             if left_good_idx >= 0:
    #                 # Hit righthand end in the lat/row
    #                 right_good_idx = xidx - 1
    #                 print(f"* row {yidx:5d}: {left_good_idx:5d} {right_good_idx:5d}")
    #                 row_lons = lons[yidx, left_good_idx:right_good_idx + 1]
    #                 row_lats = lats[yidx, left_good_idx:right_good_idx + 1]
    #                 print(f"\trow_lons ({len(row_lons)}): [{np.amin(row_lons):+9.3f} ... {np.amax(row_lons):+9.3f}]")
    #                 print(f"\trow_lats ({len(row_lats)}): [{np.amin(row_lats):+9.3f} ... {np.amax(row_lats):+9.3f}]")
    #                 break
    #             continue
    #         if left_good_idx < 0:
    #             # First valid value in the lat/row on the lefthand side of the array
    #             left_good_idx = xidx
    #             # print(f"\t\tSet {left_good_idx = :5d}")
    #         else:
    #             # Post First value value in lat/row
    #             if xidx == nx - 1:
    #                 # Entire row valid
    #                 right_good_idx = xidx
    #                 # print(f"\t\tSet {right_good_idx = :5d}")
    #                 print(f"+ row {yidx:5d}: {left_good_idx:5d} {right_good_idx:5d}")
    #                 row_lons = lons[yidx, left_good_idx:right_good_idx + 1]
    #                 row_lats = lats[yidx, left_good_idx:right_good_idx + 1]
    #                 print(f"\trow_lons ({len(row_lons)}): [{np.amin(row_lons):+9.3f} ... {np.amax(row_lons):+9.3f}]")
    #                 print(f"\trow_lats ({len(row_lats)}): [{np.amin(row_lats):+9.3f} ... {np.amax(row_lats):+9.3f}]")
    #                 continue
    #             right_good_idx = xidx
    #             # print(f"\t\tSet ** {right_good_idx = :5d}")
    #     if left_good_idx < 0:
    #         # Entire row invalid
    #         print(f"- row {yidx:5d}:")
    # os._exit(1)

    # Find spatial extent of non-missing data
    tmp_vals = data_vals.flatten()
    tmp_lons = lons.flatten()
    tmp_lats = lats.flatten()
    good_lons = []
    good_lats = []
    ll_corner = [None, None]
    ur_corner = [None, None]
    for ii, ival in enumerate(tmp_vals):
        if np.isnan(ival):
            continue
        if ival == 0.0:
            continue
        if np.abs(tmp_lons[ii]) <= 180.0 and np.abs(tmp_lats[ii]) <= 90.0:
            good_lons.append(tmp_lons[ii])
            good_lats.append(tmp_lats[ii])
            # if ll_corner[0] is not None:
            #     if tmp_lats[ii] <= ll_corner[0]:
            #         # New Lower Latitude
            #         if tmp_lons[ii] < ll_corner[1]:
            #             # New ll corner
            #             ll_corner[0] = tmp_lats[ii]
            #             ll_corner[1] = tmp_lons[ii]
            # else:
            #     # New ll corner
            #     ll_corner[0] = tmp_lats[ii]
            #     ll_corner[1] = tmp_lons[ii]
            # if ur_corner[0] is not None:
            #     if tmp_lats[ii] >= ur_corner[0]:
            #         # New higher Latitude
            #         if tmp_lons[ii] > ur_corner[1]:
            #             # New ur corner
            #             ur_corner[0] = tmp_lats[ii]
            #             ur_corner[1] = tmp_lons[ii]
            # else:
            #     # New ur corner
            #     ur_corner[0] = tmp_lats[ii]
            #     ur_corner[1] = tmp_lons[ii]
    if verbose:
        print(f"good_lons ({len(good_lons)}): [{np.amin(good_lons)} ... {np.amax(good_lons)}]")
        print(f"good_lats ({len(good_lats)}): [{np.amin(good_lats)} ... {np.amax(good_lats)}]")
        # print(f"ll_corner: {ll_corner}")
        # print(f"ur_corner: {ur_corner}")

    return

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
