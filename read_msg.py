#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Read/Convert Meteosat Second Generation (MSG) Native Archive Format (.nat) file to GTiff and make image.

$ gdalinfo MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA_HRV.tif

    Driver: GTiff/GeoTIFF
    Files: MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA_HRV.tif
    Size is 768, 768
    Coordinate System is:
    GEOGCRS["WGS 84",
        ENSEMBLE["World Geodetic System 1984 ensemble",
            MEMBER["World Geodetic System 1984 (Transit)"],
            MEMBER["World Geodetic System 1984 (G730)"],
            MEMBER["World Geodetic System 1984 (G873)"],
            MEMBER["World Geodetic System 1984 (G1150)"],
            MEMBER["World Geodetic System 1984 (G1674)"],
            MEMBER["World Geodetic System 1984 (G1762)"],
            MEMBER["World Geodetic System 1984 (G2139)"],
            MEMBER["World Geodetic System 1984 (G2296)"],
            ELLIPSOID["WGS 84",6378137,298.257223563,
            LENGTHUNIT["metre",1]],
            ENSEMBLEACCURACY[2.0]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433]],
        CS[ellipsoidal,2],
            AXIS["geodetic latitude (Lat)",north,
                ORDER[1],
                ANGLEUNIT["degree",0.0174532925199433]],
            AXIS["geodetic longitude (Lon)",east,
                ORDER[2],
                ANGLEUNIT["degree",0.0174532925199433]],
        USAGE[
            SCOPE["Horizontal component of 3D system."],
            AREA["World."],
            BBOX[-90,-180,90,180]],
        ID["EPSG",4326]]
    Data axis to CRS axis mapping: 2,1
    Origin = (-855100.436345000052825,-2638000.000000000000000)
    Pixel Size = (3000.000000000000000,-3000.000000000000000)
    Metadata:
      AREA_OR_POINT=Area
    Image Structure Metadata:
      INTERLEAVE=BAND
    Corner Coordinates:
    Upper Left  ( -855100.436,-2638000.000) (Invalid angle,Invalid angle)
    Lower Left  ( -855100.436,-4942000.000) (Invalid angle,Invalid angle)
    Upper Right ( 1448899.564,-2638000.000) (Invalid angle,Invalid angle)
    Lower Right ( 1448899.564,-4942000.000) (Invalid angle,Invalid angle)
    Center      (  296899.564,-3790000.000) (Invalid angle,Invalid angle)
    Band 1 Block=768x10 Type=Byte, ColorInterp=Gray
      NoData Value=-3.4e+38

read_msg.py
~~~~~~~~~~~~

$ python read_msg.py
"""
# Standard Imports
import os

# Third-Party Imports
import numpy as np
from osgeo import gdal
# from osgeo import osr
import pyresample as pr
#   conda install -c conda-forge satpy
from satpy import Scene
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
# from rasterio.plot import show
# from rio_cogeo import cog_validate, cog_info

# STARE Imports

# Local Imports
from cloudcast.util.setup_msg import setup_msg
from cloudcast.util.readnat import readnat
from cloudcast.util.nat2tif import nat2tif

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

# Define Global Constants and State Variables
# -------------------------------------------
# Convert a nat file to geotif
make_tif = [False, True][0]
# Read the nat file
read_nat = [False, True][0]
# Make a figure
make_fig = [False, True][1]

##
# Select Domain (only one can be True)
# Return MSG full frame data (no reprojection, subsetting or interpolation)
as_full = [False, True][0]
# Return the CloudCast european resolution and domain (reprojection, subsetting and interpolation)
as_euro = [False, True][0]
# Return the CloudCast resolution and domain (reprojection, subsetting and interpolation)
as_ccast = [False, True][1]

if sum([int(as_full), int(as_euro), int(as_ccast)]) > 1:
    raise Exception("Can only use one of as_full as_euro as_ccast.")
if sum([int(as_full), int(as_euro), int(as_ccast)]) == 0:
    raise Exception("Must select one of as_full as_euro as_ccast.")

# Return using a Mercator projection (vs a Sterographic projection) (reprojection and interpolation, no subsetting)
as_merc = [False, True][0]

use_tag = 'full' if as_full else 'ccast' if as_ccast else 'euro' if as_euro else ''
if as_merc:
    use_tag = f"{use_tag}_merc"

# Display region (not data values)
as_region = [False, True][0]

verbose = [False, True][1]

def LON_TO_180(x): return ((x + 180.0) % 360.0) - 180.0
def LON_TO_360(x): return (x + 360.0) % 360.0

euro_nrows = 928
euro_ncols = 1530

##
# From scn.all_dataset_names() below
ds_names = ("HRV", "IR_016", "IR_039", "IR_087", "IR_097", "IR_108", "IR_120",
            "IR_134", "VIS006", "VIS008", "WV_062", "WV_073")

##
# Name(s) of MSG variable to work with.

# These are the base cloudcast fields
# use_dataset = ["VIS006", "IR_039", "IR_108", "IR_120"]

# HRV
use_dataset = ds_names[0]

# IR_039
# use_dataset = ds_names[2]

##
# Define path to folder
FILE_PATH = "/Users/mbauer/tmp/CloudCast/msg/"
BASENAME = ["MSG3-SEVI-MSG15-0100-NA-20170102002740.606000000Z-NA", "MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA"][1]
SUB_PATH = f"{FILE_PATH}{BASENAME}/"
FNAME = f"{SUB_PATH}{BASENAME}.nat"

if isinstance(use_dataset, str):
    TNAME = f"{SUB_PATH}{BASENAME}_{use_dataset}_{use_tag}.tif"
    if verbose:
        print(f"Making: {TNAME}")
else:
    TNAME = []
    for ds in use_dataset:
        tname = f"{SUB_PATH}{BASENAME}_{ds}_{use_tag}.tif"
        TNAME.append(tname)
        if verbose:
            print(f"Making: {TNAME[-1]}")

##
# Color map for imaging
freq_map_cmap = ["plasma", "gist_ncar_r", "bone_r", "nipy_spectral"][2]
if as_region:
    freq_map_cmap = "bone_r"

##
# Define reader (GDAL)
#   https://satpy.readthedocs.io/en/stable/api/satpy.readers.seviri_l1b_native.html
reader = "seviri_l1b_native"

if read_nat:
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
    scn = Scene(filenames = {reader:[FNAME]})
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
    if use_dataset == HRV
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
    lons, lats, data_vals = readnat(file=FNAME, calibration="radiance", dataset=use_dataset, reader=reader, dtype="float32")
    # if verbose:
    #     tmp = lons.flatten()
    #     tmp = tmp[np.abs(tmp) <= 180.0]
    #     print(f"\tlons {lons.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
    #     tmp = lats.flatten()
    #     tmp = tmp[np.abs(tmp) <= 90.0]
    #     print(f"\tlats {lats.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
    #     tmp = data_vals.flatten()
    #     tmp = tmp[~np.isnan(tmp)]
    #     print(f"\tdata_vals {data_vals.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
    #     del tmp
    if as_euro:
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
        os._exit(1)

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
    for ii, ival in enumerate(tmp_vals):
        if np.isnan(ival):
            continue
        if np.abs(tmp_lons[ii]) <= 180.0 and np.abs(tmp_lats[ii]) <= 90.0:
            good_lons.append(tmp_lons[ii])
            good_lats.append(tmp_lats[ii])
    print(f"good_lons ({len(good_lons)}): [{np.amin(good_lons)} ... {np.amax(good_lons)}]")
    print(f"good_lats ({len(good_lats)}): [{np.amin(good_lats)} ... {np.amax(good_lats)}]")

    ##
    # Generate meta data for netcdf files etc.

    ##
    # Save as a netcdf

    os._exit(1)

##
# Run MSG setup
geo_stuff = setup_msg(use_dataset, as_ccast, as_merc, as_euro)
(MSG_EPSG, MERC_EPSG, geod_crs, ccast_area_def, ccast_crs, ccast_merc_area_def,
 ccast_merc_crs, msg_area_def, msg_crs, msg_merc_area_def, msg_merc_crs) = geo_stuff
# labels = ("MSG_EPSG", "MERC_EPSG", "geod_crs", "ccast_area_def", "ccast_crs", "ccast_merc_area_def",
#  "ccast_merc_crs", "msg_area_def", "msg_crs", "msg_merc_area_def", "msg_merc_crs")
# for ii, ival in enumerate(tmp):
#     print(f"\n{labels[ii]}: {ival}")
# os._exit(1)

# merc_stuff = (MERC_EPSG, merc_area_def) if as_merc else (None, None)
# # use_area_def = [area_def, merc_area_def][int(as_merc)]
# # use_epsg = [MSG_EPSG, MERC_EPSG][int(as_merc)]
# # use_crs = [CCAST_CRS, MERC_CRS][int(as_merc)]

if make_tif:
    ##
    # Read the file
    #   scn = <class 'satpy.scene.Scene'>
    r"""
    scn.all_dataset_names() =
        ['HRV', 'IR_016', 'IR_039', 'IR_087', 'IR_097', 'IR_108', 'IR_120', 'IR_134', 'VIS006', 'VIS008', 'WV_062', 'WV_073']
    scn.available_dataset_ids() =
        [DataID(name='HRV', wavelength=WavelengthRange(min=0.5, central=0.7, max=0.9, unit='µm'), resolution=1000.134348869, calibration=<calibration.reflectance>, modifiers=()),
         DataID(name='IR_016', wavelength=WavelengthRange(min=1.5, central=1.64, max=1.78, unit='µm'), resolution=3000.403165817, calibration=<calibration.reflectance>, modifiers=()),
         DataID(name='IR_039', wavelength=WavelengthRange(min=3.48, central=3.92, max=4.36, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_087', wavelength=WavelengthRange(min=8.3, central=8.7, max=9.1, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_097', wavelength=WavelengthRange(min=9.38, central=9.66, max=9.94, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_108', wavelength=WavelengthRange(min=9.8, central=10.8, max=11.8, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_120', wavelength=WavelengthRange(min=11.0, central=12.0, max=13.0, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_134', wavelength=WavelengthRange(min=12.4, central=13.4, max=14.4, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='VIS006', wavelength=WavelengthRange(min=0.56, central=0.635, max=0.71, unit='µm'), resolution=3000.403165817, calibration=<calibration.reflectance>, modifiers=()),
         DataID(name='VIS008', wavelength=WavelengthRange(min=0.74, central=0.81, max=0.88, unit='µm'), resolution=3000.403165817, calibration=<calibration.reflectance>, modifiers=()),
         DataID(name='WV_062', wavelength=WavelengthRange(min=5.35, central=6.25, max=7.15, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='WV_073', wavelength=WavelengthRange(min=6.85, central=7.35, max=7.85, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
        ]
    scn.coarsest_area =
        <bound method Scene.coarsest_area of <satpy.scene.Scene object at 0x1172981f0>>
    scn.images() =
        <generator object Scene.images at 0x14849f920>
    """

    ##
    # Convert to GTIFF
    nat2tif(fname=FNAME, fvar=use_dataset, reader=reader, outdir=SUB_PATH, label=use_dataset, atag=use_tag,
            geo_dat=geo_stuff, to_full=as_full, to_ccast=as_ccast, to_merc=as_merc, to_euro=as_euro)

if make_fig:
    ##
    # Plot the tif made by make_tif or in case of as_merc modified with gdal
    #   $
    #   <class 'osgeo.gdal.Dataset'>
    #   ds.GetGeoTransform()  = (-855100.436345, 3000.0, 0.0, -2638000.0, 0.0, -3000.0)
    #   ds.GetProjection()    = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'

    if as_full:
        if as_merc:
            tif_name = TNAME.replace(".tif", "_merc.tif")
        else:
            tif_name = TNAME
    elif as_ccast:
        if as_merc:
            tif_name = TNAME.replace(".tif", "_merc.tif")
        else:
            tif_name = TNAME
    ds = gdal.Open(tif_name)
    band = ds.GetRasterBand(1)
    data_vals = band.ReadAsArray()
    if as_region:
        data_vals = np.where(data_vals >= 0, 1, 0)

    if as_full:
        if as_merc:
            use_crs = msg_merc_crs
            use_area_def = msg_merc_area_def
        else:
            use_crs = msg_crs
            use_area_def = msg_area_def
    elif as_euro:
        raise Exception("Add as_euro")
    elif as_ccast:
        if as_merc:
            use_crs = ccast_merc_crs
            use_area_def = ccast_merc_area_def
        else:
            use_crs = ccast_crs
            use_area_def = ccast_area_def
    use_extent = use_crs.bounds

    fig = plt.figure(figsize=(10, 8))
    as_geos = False
    if as_full:
        # ax = plt.axes(projection=ccrs.Miller())
        # ax = plt.axes(projection=ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0))
        # ax = plt.axes(projection=ccrs.Mercator(central_longitude=0.0, min_latitude=-80.0, max_latitude=80.0))
        # ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0))
        ax = plt.axes(projection=ccrs.Geostationary(central_longitude=0.0)); as_geos = True

        # proj = ccrs.PlateCarree(central_longitude=0.0, globe=ccrs.Globe(datum='WGS84', ellipse='WGS84'))
        # ax = plt.axes(projection=proj)
        # ax = plt.axes(projection=use_crs)
    elif as_euro:
        raise Exception("Add as_euro")
    elif as_ccast:
        # ax = plt.axes(projection=use_crs)
        ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=0.0, central_latitude=0.0, cutoff=-30))
    else:
        ax = plt.axes(projection=use_crs)

    if as_ccast:
        #   CCAST Sterographic
        #       values: [0.5235565900802612, ... 13.425487518310547]
        #   CCAST Mercator
        #       values: [0.5609534978866577, ... 14.248218536376953]
        the_alpha = 1.0
        if use_dataset == "HRV":
            # bounds = np.linspace(0.0, 15.0, num=16, endpoint=True)
            bounds = np.linspace(0.0, 25.0, num=16, endpoint=True)
            the_min = 0.0
            the_max = 25.0
        else:
            bounds = np.linspace(0.0, 4.0, num=16, endpoint=True)
            the_min = 0.0
            the_max = 4.0
        if as_region:
            bounds = np.linspace(0.0, 1.0, num=2, endpoint=True)
            the_min = 0.0
            the_max = 1.0
            the_alpha = 0.25
        a_image = plt.imshow(data_vals, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=the_alpha)
        ax.set_extent(use_extent, crs=use_crs)
    elif as_merc:
        # ax.set_extent(use_extent, crs=geod_crs)
        a_image = plt.imshow(data_vals, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper')
        ax.set_extent(use_extent, crs=use_crs)
    else:
        # print(use_extent)
        # (-81.12566375732422, -81.0744857788086, 81.12566375732422, 81.0744857788086)
        # lower_left_xy = [-81.12566375732422, 81.12566375732422]
        # upper_right_xy = [-81.0744857788086, 81.0744857788086]

        # x_proj = use_crs.transform_point(use_extent[0], use_extent[2], proj)
        # y_proj = use_crs.transform_point(use_extent[1], use_extent[3], proj)
        # print(f"x_proj {x_proj}")
        # print(f"y_proj: {y_proj}")
        # use_extent = [x_proj[0], x_proj[1], y_proj[0], y_proj[1]]

        if use_dataset == "HRV":
            use_extent = [-2750000, 2750000, -5500000, 5500000]
        else:
            use_extent = [-5500000, 5500000, -5500000, 5500000]
        # use_extent = [-6500000, 6500000, -6500000, 6500000]
        # use_extent = [-6500000, 6500000, -6500000, 6500000]
        # print(use_extent)

        #   CCAST Sterographic
        #       values: [0.5235565900802612, ... 13.425487518310547]
        #   CCAST Mercator
        #       values: [0.5609534978866577, ... 14.248218536376953]
        the_alpha = 1.0
        if use_dataset == "HRV":
            # bounds = np.linspace(0.0, 15.0, num=16, endpoint=True)
            bounds = np.linspace(0.0, 25.0, num=16, endpoint=True)
            the_min = 0.0
            the_max = 25.0
        else:
            bounds = np.linspace(0.0, 4.0, num=16, endpoint=True)
            the_min = 0.0
            the_max = 4.0
        if as_region:
            bounds = np.linspace(0.0, 1.0, num=2, endpoint=True)
            the_min = 0.0
            the_max = 1.0
            the_alpha = 0.25
        a_image = plt.imshow(data_vals, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=the_alpha)
    ax.add_feature(cfeature.COASTLINE, alpha=1)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=2, color='black', alpha=0.5, linestyle='--')

    if not as_region:
        ##
        # Color bar
        use_shrink = 0.9
        if not as_geos:
            use_shrink = 0.5

        cbar = fig.colorbar(a_image, location="right", orientation='vertical', extend="neither", ticks=bounds, shrink=use_shrink,
                            spacing='uniform', fraction=0.15, pad=0.1, aspect=30, drawedges=False,
                            ax=ax)
        cbar.ax.tick_params(labelsize=8)
        cbar.solids.set(alpha=1)

    plt.show()
    # To trim https://www.zackwebster.com/tools/image-trim#_

print("Done")

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
