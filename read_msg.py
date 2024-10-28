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
from cloudcast.util.readnat import readnat
from cloudcast.util.nat2tif import nat2tif

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

# Define Global Constants and State Variables
# -------------------------------------------
make_tif = [False, True][0]
read_nat = [False, True][0]
make_fig = [False, True][1]
as_spain = [False, True][0]
as_merc = [False, True][1]

##
# From scn.all_dataset_names() below
ds_names = ("HRV", "IR_016", "IR_039", "IR_087", "IR_097", "IR_108", "IR_120",
            "IR_134", "VIS006", "VIS008", "WV_062", "WV_073")

##
# Name(s) of MSG data to work with.
use_dataset = ds_names[0]
# These are the base cloudcast fields
# use_dataset = ["VIS006", "IR_039", "IR_108", "IR_120"]

##
# Define path to folder
FILE_PATH = "/Users/mbauer/tmp/CloudCast/msg/"
BASENAME = ["MSG3-SEVI-MSG15-0100-NA-20170102002740.606000000Z-NA", "MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA"][1]
SUB_PATH = f"{FILE_PATH}{BASENAME}/"
FNAME = f"{SUB_PATH}{BASENAME}.nat"
if isinstance(use_dataset, str):
    TNAME = f"{SUB_PATH}{BASENAME}_{use_dataset}{'_spain' if as_spain else ''}.tif"
else:
    TNAME = []
    for ds in use_dataset:
        TNAME.append(f"{SUB_PATH}{BASENAME}_{ds}{'_spain' if as_spain else ''}.tif")


"""
https://kylebarron.dev/deck.gl-raster/

"""

# ##
# # Check GeoTIF
# gtif_name = "/Volumes/saved/data/CloudCast/msg/MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA/MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA_HRV.tif"
# gtif_cog_name = gtif_name.replace(".tif", "_cog.tif")
# # $ rio cogeo create /Volumes/saved/data/CloudCast/msg/MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA/MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA_HRV.tif /Volumes/saved/data/CloudCast/msg/MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA/MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA_HRV_cog.tif
# # rio_cogeo.create(gtif_name, gtif_cog_name)
# is_valid, errors, warnings = cog_validate(gtif_cog_name)
# print(is_valid, errors, warnings)
# # cog_info(gtif_cog_name)

# os._exit(1)

##
# Color map for imaging
freq_map_cmap = ["plasma", "gist_ncar_r", "bone_r"][1]


##
# Define reader (GDAL)
#   https://satpy.readthedocs.io/en/stable/api/satpy.readers.seviri_l1b_native.html
reader = "seviri_l1b_native"

##
# MSG base projection WGS84 - World Geodetic System 1984
#   https://epsg.io/4326
MSG_EPSG = 4326

# Geodetic:
#   A 3D/spherical CRS based on latitude and longitude where geographical distance and coordinates are measured in degrees.
geod_crs = ccrs.Geodetic()

# The MSG data is provided as Full Disk, meaning that roughly the complete North-South extent of
#   the globe from the Atlantic to the Indian Ocean is present in each file.
if as_spain:
    area_id = "Spain"
    description = "Geographical Coordinate System clipped on Spain"
    proj_id = "Spain"
    # Specify projection parameters
    proj_dict = {"proj": "longlat", "ellps": "WGS84", "datum": "WGS84"}
    ##
    # Calculate the width and height of the aoi in pixels
    llx = -9.5 # lower left x coordinate in degrees
    lly = 35.9 # lower left y coordinate in degrees
    urx = 3.3 # upper right x coordinate in degrees
    ury = 43.8 # upper right y coordinate in degrees
    resolution = 0.005 # target resolution in degrees
    ##
    # Calculate the number of pixels
    width = int((urx - llx) / resolution)
    height = int((ury - lly) / resolution)
    area_extent = (llx, lly, urx, ury)
    ##
    # Define the area
    #   <class 'pyresample.geometry.AreaDefinition'>
    area_def = pr.geometry.AreaDefinition(area_id, proj_id, description, proj_dict, width, height, area_extent)
else:
    ##
    # Create some information on the reference system
    CCAST_HEIGHT = 768
    CCAST_WIDTH = 768
    lower_left_xy = [-855100.436345, -4942000.0]
    upper_right_xy = [1448899.563655, -2638000.0]
    # Define the area
    #   <class 'pyresample.geometry.AreaDefinition'>
    area_def = pr.geometry.AreaDefinition('areaD', 'Europe', 'areaD',
                                          {'lat_0': '90.00', 'lat_ts': '50.00', 'lon_0': '5', 'proj': 'stere', 'ellps': 'WGS84'},
                                          CCAST_HEIGHT, CCAST_WIDTH,
                                          (lower_left_xy[0], lower_left_xy[1], upper_right_xy[0], upper_right_xy[1]))
    #   +proj=stere +lat_0=90 +lat_ts=50 +lon_0=5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs +type=crs
    # print(area_def.proj4_string)

    ##
    # Form a cartopy CRS
    #   <class 'pyresample.utils.cartopy.Projection'>
    CCAST_CRS = area_def.to_cartopy_crs()
    # print(MERC_CRS.bounds)
    #   CCAST_CRS.bounds: (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)

    """
    CCAST_CRS
        area_def.corners = [
            (-12.920934886492649, 62.403066090517555), UL
            (33.73865749382469,   60.15059617915855),  UR
            (21.32880156090482,   40.92817004367345),  LR
            (-4.802566482888071,  42.068097533886025)] LL

        area_def.area_extent  (-855100.436345, -4942000.0, 1448899.563655, -2638000.0)
        CCAST_CRS.bounds      (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)

        lower_left_xy  = [-855100.436345, -4942000.0] => (-4.816534709314307, 42.053336570266744)
        upper_right_xy = [1448899.563655, -2638000.0] => (33.77742545811928,  60.15622631725923)
    """
    # print(area_def.corners)
    # print(area_def.area_extent)
    # print(CCAST_CRS.bounds)
    # ur_xy = geod_crs.transform_point(upper_right_xy[0], upper_right_xy[1], CCAST_CRS)
    # ll_xy = geod_crs.transform_point(lower_left_xy[0], lower_left_xy[1], CCAST_CRS)
    # print(f"ll_xy: {ll_xy}")
    # print(f"ur_xy: {ur_xy}")

    ##
    # Form a Mercator CSR version of CCAST_CRS
    use_pseudo = [False, True][0]
    MERC_EPSG = 3857 if use_pseudo else 3395
    MERC_CRS = ccrs.epsg(MERC_EPSG)
    if use_pseudo:
        ##
        # Mercator version EPSG:3857 for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI
        #   https://epsg.io/3857
        area_id = "EPSG:3857"
        description = "Pseudo-Mercator"
        proj_id = "EPSG:3857"
    else:
        ##
        # World Mercator version EPSG:3395 for WGS 84
        #   https://epsg.io/3395
        area_id = "EPSG:3395"
        description = "Mercator"
        proj_id = "EPSG:3395"

    #   ll_xy: (-536174.1912289965, 5140343.824785849)
    #   ur_xy: (3760085.802305594, 8397504.685448818)
    ur_xy = MERC_CRS.transform_point(upper_right_xy[0], upper_right_xy[1], CCAST_CRS)
    ll_xy = MERC_CRS.transform_point(lower_left_xy[0], lower_left_xy[1], CCAST_CRS)
    # print(f"ll_xy: {ll_xy}")
    # print(f"ur_xy: {ur_xy}")

    ##
    # Specify projection parameters
    if use_pseudo:
        projection = {'proj': 'merc', 'lat_ts': 0, 'lon_0': 0, 'a': 6378137, 'b': 6378137, 'x_0': 0, 'y_0': 0, 'k': 1, 'units': 'm'}
    else:
        projection = {'proj': 'merc','x_0': 0, 'y_0': 0, 'k': 1, 'units': 'm', 'no_defs': None, 'datum': 'WGS84', 'type': 'crs'}

    ##
    # Define the area
    merc_area_def = pr.geometry.AreaDefinition(area_id, 'Europe', proj_id,
                                               projection, CCAST_HEIGHT, CCAST_WIDTH,
                                               (ll_xy[0], ll_xy[1], ur_xy[0], ur_xy[1]))
    #   proj4_string: +proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs
    # print(f"proj4_string: {merc_area_def.proj4_string}")

    MERC_CRS = merc_area_def.to_cartopy_crs()

    """
    MERC_CRS
        merc_area_corners = [(-4.7914084331636335, 60.14672953639852),  UL
                             (33.75229918196862, 60.14672953639852),    UR
                             (33.75229918196862, 42.06753197287029),    LR
                             (-4.7914084331636335, 42.06753197287029)]  LL

        merc_area_def.area_extent (-536174.1912289965, 5140343.824785849, 3760085.802305594, 8397504.685448818)
        MERC_CRS.bounds           (-536174.1912289965, 5140343.824785849, 3760085.802305594, 8397504.685448818)

        ll_xy  = (-536174.1912289965, 5140343.824785849) => (-4.816534709314306, 42.05333657026676)
        ur_xy  = (3760085.802305594, 8397504.685448818)  => (33.77742545811928, 60.15622631725923)
    """
    # print(merc_area_def.corners)
    # print(merc_area_def.area_extent)
    # print(MERC_CRS.bounds)
    # ur_xya = geod_crs.transform_point(ur_xy[0], ur_xy[1], MERC_CRS)
    # ll_xya = geod_crs.transform_point(ll_xy[0], ll_xy[1], MERC_CRS)
    # print(f"ll_xy: {ll_xya}")
    # print(f"ur_xy: {ur_xya}")
    # os._exit(1)

merc_stuff = (MERC_EPSG, merc_area_def) if as_merc else (None, None)

# use_area_def = [area_def, merc_area_def][int(as_merc)]
# use_epsg = [MSG_EPSG, MERC_EPSG][int(as_merc)]
# use_crs = [CCAST_CRS, MERC_CRS][int(as_merc)]

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
    lons (3712, 3712)
    lats (3712, 3712)
    data_vals (4, 3712, 3712)
    """
    lons, lats, data_vals = readnat(file=FNAME, calibration="radiance", dataset=use_dataset, reader=reader, dtype="float32")
    print(f"lons {lons.shape}")
    print(f"lats {lats.shape}")
    print(f"data_vals {data_vals.shape}")

    ##
    # Generate meta data for netcdf files etc.


    ##
    # Save as a netcdf

    #os._exit(1)

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
    scn = Scene(filenames = {reader:[FNAME]})
    # print(scn.values())

    ##
    # Extract data set names
    dataset_names = scn.all_dataset_names()
    """
        HRV
        IR_016
        IR_039
        IR_087
        IR_097
        IR_108
        IR_120
        IR_134
        VIS006
        VIS008
        WV_062
        WV_073
    """

    ##
    # Convert to GTIFF
    nat2tif(file=FNAME, calibration="radiance", area_def=area_def, dataset=use_dataset,
            reader=reader, outdir=SUB_PATH, label=use_dataset, dtype="float32", radius=16000,
            epsilon=0.5, nodata=-3.4E+38, use_epsg=MSG_EPSG, merc_vals=merc_stuff, to_spain=as_spain)

if make_fig:
    ##
    # Plot the tif made by make_tif or in case of as_merc modified with gdal
    #   $
    #   <class 'osgeo.gdal.Dataset'>
    #   ds.GetGeoTransform()  = (-855100.436345, 3000.0, 0.0, -2638000.0, 0.0, -3000.0)
    #   ds.GetProjection()    = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'

    ds = gdal.Open(TNAME.replace(".tif", "_merc.tif") if as_merc else TNAME)
    band = ds.GetRasterBand(1)
    data = band.ReadAsArray()
    use_extent = MERC_CRS.bounds if as_merc else CCAST_CRS.bounds

    if as_spain:
        plt.imshow(data)
    else:
        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes(projection=MERC_CRS if as_merc else CCAST_CRS)
        # a_image = plt.imshow(data, cmap=freq_map_cmap, transform=CCAST_CRS, extent=use_extent, origin='upper')
        if as_merc:
            # ax.set_extent(use_extent, crs=geod_crs)
            a_image = plt.imshow(data, cmap=freq_map_cmap, transform=MERC_CRS, extent=use_extent, origin='upper')
            ax.set_extent(use_extent, crs=MERC_CRS)
        else:
            a_image = plt.imshow(data, cmap=freq_map_cmap, transform=CCAST_CRS, extent=use_extent, origin='upper')
        ax.add_feature(cfeature.COASTLINE, alpha=1)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=2, color='black', alpha=0.5, linestyle='--')

    ##
    # Color bar
    #   CCAST Sterographic
    #       values: [0.5235565900802612, ... 13.425487518310547]
    #   CCAST Mercator
    #       values: [0.5609534978866577, ... 14.248218536376953]
    bounds = np.linspace(0.0, 15.0, num=16, endpoint=True)
    # print(bounds)
    cbar = fig.colorbar(a_image, location="right", orientation='vertical', extend="neither", ticks=bounds, shrink=0.7 if as_merc else 0.9,
                        spacing='uniform', fraction=0.15, pad=0.1, aspect=30, drawedges=False,
                        ax=ax)
    cbar.ax.tick_params(labelsize=8)
    cbar.solids.set(alpha=1)

    plt.show()

print("Done")

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
