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
import matplotlib.pyplot as plt

# STARE Imports

# Local Imports

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

###############################################################################
# PUBLIC nat2tif()
# ----------------
def nat2tif(file, calibration, area_def, dataset, reader, outdir, label, dtype, radius, epsilon, nodata):
    # open the file
    scn = Scene(filenames = {reader: [file]})
    # let us check that the specified data set is actually available
    scn_names = scn.all_dataset_names()
    # raise exception if dataset is not present in available names
    if dataset not in scn_names:
        raise Exception("Specified dataset is not available.")
    # we need to load the data, different calibration can be chosen
    scn.load([dataset], calibration=calibration)
    # let us extract the longitude and latitude data
    lons, lats = scn[dataset].area.get_lonlats()
    # now we can apply a swath definition for our output raster
    swath_def = pr.geometry.SwathDefinition(lons=lons, lats=lats)
    # and finally we also extract the data
    values = scn[dataset].values
    # we will now change the datatype of the arrays
    # depending on the present data this can be changed
    lons = lons.astype(dtype)
    lats = lats.astype(dtype)
    values = values.astype(dtype)
    # now we can already resample our data to the area of interest
    values = pr.kd_tree.resample_nearest(swath_def, values,
                                             area_def,
                                             radius_of_influence=radius, # in meters
                                             epsilon=epsilon,
                                             fill_value=False)
    # we are going to check if the outdir exists and create it if it doesnt
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # let us join our filename based on the input file's basename
    outname = os.path.join(outdir, os.path.basename(file)[:-4] + "_" + str(label) + ".tif")
    # now we define some metadata for our raster file
    cols = values.shape[1]
    rows = values.shape[0]
    pixelWidth = (area_def.area_extent[2] - area_def.area_extent[0]) / cols
    pixelHeight = (area_def.area_extent[1] - area_def.area_extent[3]) / rows
    originX = area_def.area_extent[0]
    originY = area_def.area_extent[3]
    # here we actually create the file
    driver = gdal.GetDriverByName("GTiff")
    outRaster = driver.Create(outname, cols, rows, 1)
    # writing the metadata
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    # creating a new band and writting the data
    outband = outRaster.GetRasterBand(1)
    outband.SetNoDataValue(nodata) #specified no data value by user
    outband.WriteArray(np.array(values)) # writting the values
    outRasterSRS = osr.SpatialReference() # create CRS instance
    outRasterSRS.ImportFromEPSG(4326) # get info for EPSG 4326
    outRaster.SetProjection(outRasterSRS.ExportToWkt()) # set CRS as WKT
    # clean up
    outband.FlushCache()
    outband = None
    outRaster = None


# Define Global Constants and State Variables
# -------------------------------------------
make_tif = [False, True][0]

##
# Define path to folder
FILE_PATH = "/Volumes/val/data/CloudCast/msg/"
BASENAME = ["MSG3-SEVI-MSG15-0100-NA-20170102002740.606000000Z-NA", "MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA"][1]
SUB_PATH = f"{FILE_PATH}/{BASENAME}/"
FNAME = f"{SUB_PATH}{BASENAME}.nat"
TNAME = f"{SUB_PATH}{BASENAME}_HRV.tif"

if make_tif:
    ##
    # Define reader
    reader = "seviri_l1b_native"

    ##
    # Read the file
    scn = Scene(filenames = {reader:[FNAME]})

    ##
    # Extract data set names
    dataset_names = scn.all_dataset_names()

    ##
    # print available datasets
    ## print('\n'.join(map(str, dataset_names)))
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

    # The MSG data is provided as Full Disk, meaning that roughly the complete North-South extent of
    #   the globe from the Atlantic to the Indian Ocean is present in each file.

    ##
    # Create some information on the reference system
    # CCAST_HEIGHT = 768
    # CCAST_WIDTH = 768
    # lower_left_xy = [-855100.436345, -4942000.0]
    # upper_right_xy = [1448899.563655, -2638000.0]
    # area_def = pr.geometry.AreaDefinition('areaD', 'Europe', 'areaD',
    #                                       {'lat_0': '90.00', 'lat_ts': '50.00',
    #                                        'lon_0': '5', 'proj': 'stere', 'ellps': 'WGS84'},
    #                                       CCAST_HEIGHT, CCAST_WIDTH,
    #                                       (lower_left_xy[0], lower_left_xy[1],
    #                                        upper_right_xy[0], upper_right_xy[1]))
    # CCAST_CRS = area_def.to_cartopy_crs()

    area_id = "Spain"
    description = "Geographical Coordinate System clipped on Spain"
    proj_id = "Spain"
    # specifing some parameters of the projection
    proj_dict = {"proj": "longlat", "ellps": "WGS84", "datum": "WGS84"}
    # calculate the width and height of the aoi in pixels
    llx = -9.5 # lower left x coordinate in degrees
    lly = 35.9 # lower left y coordinate in degrees
    urx = 3.3 # upper right x coordinate in degrees
    ury = 43.8 # upper right y coordinate in degrees
    resolution = 0.005 # target resolution in degrees
    # calculating the number of pixels
    width = int((urx - llx) / resolution)
    height = int((ury - lly) / resolution)
    area_extent = (llx,lly,urx,ury)
    # defining the area
    area_def = pr.geometry.AreaDefinition(area_id, proj_id, description, proj_dict, width, height, area_extent)
    # print(area_def)

    nat2tif(file = FNAME,
            calibration = "radiance",
            area_def = area_def,
            dataset = "HRV",
            reader = reader,
            outdir = SUB_PATH,
            label = "HRV",
            dtype = "float32",
            radius = 16000,
            epsilon = .5,
            nodata = -3.4E+38)

##
# Plot
ds = gdal.Open(TNAME)
band = ds.GetRasterBand(1)
data = band.ReadAsArray()
plt.imshow(data)
plt.show()

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
