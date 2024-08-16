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

# Third-Party Imports
import numpy as np
from osgeo import gdal
from osgeo import osr
import pyresample as pr
#   conda install -c conda-forge satpy
from satpy import Scene
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

# STARE Imports

# Local Imports
from cloudcast.util.readnat import readnat
from cloudcast.util.nat2tif import nat2tif

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

# ###############################################################################
# # PUBLIC nat2tif()
# # ----------------
# def nat2tif(file, calibration, area_def, dataset, reader, outdir, label, dtype, radius, epsilon, nodata):
#     ##
#     # Open the file w/ satpy, which uses Xarray
#     scn = Scene(filenames = {reader: [file]})

#     ##
#     # Check the specified data set is actually available
#     scn_names = scn.all_dataset_names()
#     if dataset not in scn_names:
#         raise Exception("Specified dataset is not available.")

#     ##
#     # Load the data, different calibration can be chosen
#     scn.load([dataset], calibration=calibration)

#     ##
#     # Extract the longitude and latitude data
#     lons, lats = scn[dataset].area.get_lonlats()

#     ##
#     # Apply a swath definition for our output raster
#     swath_def = pr.geometry.SwathDefinition(lons=lons, lats=lats)

#     ##
#     # Extract the data values
#     values = scn[dataset].values

#     ##
#     # Change the datatype of the arrays depending on the present data this can be changed
#     lons = lons.astype(dtype)
#     lats = lats.astype(dtype)
#     values = values.astype(dtype)

#     ##
#     # Resample our data to the area of interest
#     values = pr.kd_tree.resample_nearest(swath_def, values,
#                                              area_def,
#                                              radius_of_influence=radius, # in meters
#                                              epsilon=epsilon,
#                                              fill_value=False)

#     ##
#     # Check if outdir exists, otherwise create it
#     if not os.path.exists(outdir):
#         os.makedirs(outdir)

#     ##
#     # Join our filename based on the input file's basename
#     outname = os.path.join(outdir, os.path.basename(file)[:-4] + "_" + str(label) + ".tif")
#     outname = outname.replace(".tif", '_spain.tif' if as_spain else '.tif')

#     ##
#     # Define some metadata
#     cols = values.shape[1]
#     rows = values.shape[0]
#     pixelWidth = (area_def.area_extent[2] - area_def.area_extent[0]) / cols
#     pixelHeight = (area_def.area_extent[1] - area_def.area_extent[3]) / rows
#     originX = area_def.area_extent[0]
#     originY = area_def.area_extent[3]

#     ##
#     # Create the file
#     driver = gdal.GetDriverByName("GTiff")
#     outRaster = driver.Create(outname, cols, rows, 1)

#     ##
#     # Save the metadata
#     outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

#     ##
#     # Create a new band and write the data
#     outband = outRaster.GetRasterBand(1)
#     outband.SetNoDataValue(nodata) #specified no data value by user
#     outband.WriteArray(np.array(values)) # writing the values
#     outRasterSRS = osr.SpatialReference() # create CRS instance
#     outRasterSRS.ImportFromEPSG(4326) # get info for EPSG 4326
#     outRaster.SetProjection(outRasterSRS.ExportToWkt()) # set CRS as WKT

#     ##
#     # Clean up
#     outband.FlushCache()
#     outband = None
#     outRaster = None

#     return

# Define Global Constants and State Variables
# -------------------------------------------
make_tif = [False, True][0]
read_nat = [False, True][1]
as_spain = [False, True][0]

##
# From scn.all_dataset_names() below
ds_names = ("HRV", "IR_016", "IR_039", "IR_087", "IR_097", "IR_108", "IR_120",
            "IR_134", "VIS006", "VIS008", "WV_062", "WV_073")

##
# Name(s) of MSG data to work with.
use_dataset = ds_names[0]
# These are the base cloudcast fields
use_dataset = ["VIS006", "IR_039", "IR_108", "IR_120"]

##
# Define path to folder
FILE_PATH = "/Volumes/val/data/CloudCast/msg/"
BASENAME = ["MSG3-SEVI-MSG15-0100-NA-20170102002740.606000000Z-NA", "MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA"][1]
SUB_PATH = f"{FILE_PATH}/{BASENAME}/"
FNAME = f"{SUB_PATH}{BASENAME}.nat"
if isinstance(use_dataset, str):
    TNAME = f"{SUB_PATH}{BASENAME}_{use_dataset}{'_spain' if as_spain else ''}.tif"
else:
    TNAME = []
    for ds in use_dataset:
        TNAME.append(f"{SUB_PATH}{BASENAME}_{ds}{'_spain' if as_spain else ''}.tif")

##
# Color map for imaging
freq_map_cmap = "plasma"

##
# Define reader (GDAL)
#   https://satpy.readthedocs.io/en/stable/api/satpy.readers.seviri_l1b_native.html
reader = "seviri_l1b_native"

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
    area_extent = (llx,lly,urx,ury)
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
    ##
    # Form a cartopy CRS
    CCAST_CRS = area_def.to_cartopy_crs()

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

    os._exit(1)

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
    nat2tif(file = FNAME,
            calibration = "radiance",
            area_def = area_def,
            dataset = use_dataset,
            reader = reader,
            outdir = SUB_PATH,
            label = use_dataset,
            dtype = "float32",
            radius = 16000,
            epsilon = .5,
            nodata = -3.4E+38)

##
# Plot
ds = gdal.Open(TNAME)

os._exit(1)

band = ds.GetRasterBand(1)
data = band.ReadAsArray()
if as_spain:
    plt.imshow(data)
else:
    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=CCAST_CRS)
    a_image = plt.imshow(data, cmap=freq_map_cmap, transform=CCAST_CRS, extent=CCAST_CRS.bounds, origin='upper')
    ax.add_feature(cfeature.COASTLINE, alpha=0.5)
plt.show()

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<