#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Read/Convert Meteosat Second Generation (MSG) Native Archive Format (.nat) file to GTiff and make image.

read_msg.py
~~~~~~~~~~~~

$ python read_msg.py
"""
# Standard Imports
import os
import pickle

# Third-Party Imports
import numpy as np
from osgeo import gdal
import pyresample as pr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

# STARE Imports

# Local Imports
from cloudcast.cfg.setup_msg import setup_msg
from cloudcast.util.natread import natread
from cloudcast.util.nat2tif import nat2tif

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

# Define Global Constants and State Variables
# -------------------------------------------

##
# Convert a nat file to geotif
make_tif = [False, True][0]

##
# Read the nat file only
read_nat = [False, True][1]

##
# Make a figure
make_fig = [False, True][1]

##
# Select Domain (only one can be True)

##
# Return MSG full frame data (no reprojection, subsetting or interpolation)
as_full = [False, True][0]

##
# Return the CloudCast european resolution and domain (reprojection, subsetting and interpolation)
as_euro = [False, True][1]

##
# Return the CloudCast resolution and domain (reprojection, subsetting and interpolation)
as_ccast = [False, True][0]

if sum([int(as_full), int(as_euro), int(as_ccast)]) > 1:
    raise Exception("Can only use one of as_full as_euro as_ccast.")
if sum([int(as_full), int(as_euro), int(as_ccast)]) == 0:
    raise Exception("Must select one of as_full as_euro as_ccast.")

##
# Return using a Mercator projection (vs a Sterographic projection) (reprojection and interpolation, no subsetting)
as_merc = [False, True][0]

#
# Return using a Lambert Conformal projection (vs a Sterographic projection) (reprojection and interpolation, no subsetting)
#   Doesn't work for full because of latitude span
as_lcc = [False, True][0]
if as_lcc and as_full:
    raise Exception("Due to latitude span as_full can't be used with as_lcc.")

##
# Tag to id new files
use_tag = 'full' if as_full else 'ccast' if as_ccast else 'euro' if as_euro else ''
if as_merc:
    use_tag = f"{use_tag}_merc"
if as_lcc:
    use_tag = f"{use_tag}_lcc"

##
# Display region (not data values)
as_region = [False, True][0]

##
# Obvious
verbose = [False, True][1]

##
# From scn.all_dataset_names() below
ds_names = ("HRV", "IR_016", "IR_039", "IR_087", "IR_097", "IR_108", "IR_120",
            "IR_134", "VIS006", "VIS008", "WV_062", "WV_073")

##
# Name(s) of MSG variable to work with.
#   These are the base cloudcast fields  ["VIS006", "IR_039", "IR_108", "IR_120"]

##
# Just HRV
# use_dataset = ds_names[0]

##
# Just IR_108
use_dataset = ds_names[3]

##
# Define path to folder
FILE_PATH = "/Users/mbauer/tmp/CloudCast/msg/"
BASENAME = ["MSG3-SEVI-MSG15-0100-NA-20170102002740.606000000Z-NA",
            "MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA",
            "MSG3-SEVI-MSG15-0100-NA-20170104005740.099000000Z-NA"][1]
SUB_PATH = f"{FILE_PATH}{BASENAME}/"
FNAME = f"{SUB_PATH}{BASENAME}.nat"
TNAME = f"{SUB_PATH}{BASENAME}_{use_dataset}_{use_tag}.tif"
if verbose:
    print(f"{'Making' if make_tif else 'Reading'}: {TNAME}")

##
# Color map for imaging
freq_map_cmap = ["plasma", "gist_ncar_r", "bone_r", "nipy_spectral"][2]
if as_region:
    freq_map_cmap = "bone_r"

##
# Define reader (GDAL)
#   https://satpy.readthedocs.io/en/stable/api/satpy.readers.seviri_l1b_native.html
reader = "seviri_l1b_native"

if as_euro:
    ##
    # Read raw_lons, raw_lats for CloudCast raw for lon/lat domain matching with MSG
    #   raw_lons (928, 1530): [-69.2706298828125  ... 69.2706298828125]
    #   raw_lats (928, 1530): [ 26.67105484008789 ... 81.09877014160156]
    with open(f"/Volumes/saved/hidden/cloudcast/raw_coords.pkl", 'rb') as f:
        tmp = pickle.load(f)
        raw_lons, raw_lats = tmp
        del tmp
    if verbose:
        tmp = raw_lons.flatten()
        tmp = tmp[np.abs(tmp) <= 180.0]
        print(f"\traw_lons {raw_lons.shape}: [{np.amin(tmp)} ... {np.amax(tmp)}]")
        tmp = raw_lats.flatten()
        tmp = tmp[np.abs(tmp) <= 90.0]
        print(f"\traw_lats {raw_lats.shape}: [{np.amin(tmp)} ... {np.amax(tmp)}]")
else:
    raw_lons = np.zeros((1,1))
    raw_lats = np.zeros((1,1))
if read_nat:
    ##
    # Read the file
    natread(fname=FNAME, fvar=use_dataset, reader=reader, to_euro=as_euro, euro_lons=raw_lons, euro_lats=raw_lats)
    os._exit(1)

##
# Run MSG setup
geo_stuff = setup_msg(use_dataset, as_ccast, as_euro, as_merc, as_lcc)
(MSG_EPSG, MERC_EPSG, LCC_EPSG, geod_crs,
 ccast_area_def, ccast_crs, ccast_merc_area_def, ccast_merc_crs, ccast_lcc_area_def, ccast_lcc_crs,
 msg_area_def, msg_crs, msg_merc_area_def, msg_merc_crs) = geo_stuff
# labels = ("MSG_EPSG", "MERC_EPSG", "geod_crs", "ccast_area_def", "ccast_crs", "ccast_merc_area_def", "ccast_lcc_area_def", "ccast_lcc_crs",
#  "ccast_merc_crs", "msg_area_def", "msg_crs", "msg_merc_area_def", "msg_merc_crs")
# for ii, ival in enumerate(tmp):
#     print(f"\n{labels[ii]}: {ival}")
# os._exit(1)

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
    """

    ##
    # Convert nat to GeoTIF
    nat2tif(fname=FNAME, fvar=use_dataset, reader=reader, outdir=SUB_PATH, label=use_dataset, atag=use_tag,
            geo_dat=geo_stuff, to_full=as_full, to_ccast=as_ccast, to_euro=as_euro, to_merc=as_merc, to_lcc=as_lcc)

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
        elif as_lcc:
            tif_name = TNAME.replace(".tif", "_lcc.tif")
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
        elif as_lcc:
            use_crs = ccast_lcc_crs
            use_area_def = ccast_lcc_area_def
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
        #ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0))
        ax = plt.axes(projection=ccrs.Geostationary(central_longitude=0.0)); as_geos = True

        # proj = ccrs.PlateCarree(central_longitude=0.0, globe=ccrs.Globe(datum='WGS84', ellipse='WGS84'))
        # ax = plt.axes(projection=proj)
        # ax = plt.axes(projection=use_crs)
    elif as_euro:
        raise Exception("Add as_euro")
    elif as_ccast:
        ax = plt.axes(projection=use_crs)
        # ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=0.0, central_latitude=0.0, cutoff=-30))
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

            bounds = np.linspace(0.0, 60.0, num=16, endpoint=True)
            the_min = 0.0
            the_max = 60.0

        if as_region:
            bounds = np.linspace(0.0, 1.0, num=2, endpoint=True)
            the_min = 0.0
            the_max = 1.0
            the_alpha = 0.25
        if as_lcc:
            a_image = plt.imshow(data_vals, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=the_alpha)
            ax.set_extent(use_extent, crs=use_crs)
        else:
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
            bounds = np.linspace(0.0, 25.0, num=11, endpoint=True)
            the_min = 0.0
            the_max = 25.0
        else:
            bounds = np.linspace(0.0, 4.0, num=16, endpoint=True)
            the_min = 0.0
            the_max = 4.0

            bounds = np.linspace(0.0, 100.0, num=26, endpoint=True)
            the_min = 0.0
            the_max = 100.0
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
    # To trim region

print("Done")

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
