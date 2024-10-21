#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Read/Convert Meteosat Second Generation (MSG) Native Archive Format (.nat) file to GTiff.

nat2tif.py
~~~~~~~~~~~~
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

# STARE Imports

# Local Imports

##
# List of Public objects from this module.
__all__ = ['nat2tif']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

###############################################################################
# PUBLIC nat2tif()
# ----------------
def nat2tif(file: str, calibration: str, area_def: pr.geometry.AreaDefinition, dataset: str, reader: str, outdir: str, label: str, dtype: str, radius: int, epsilon: float, nodata: float, to_cog_tif: bool, to_spain: bool):
    ##
    # Open the file w/ satpy, which uses Xarray
    scn = Scene(filenames = {reader: [file]})

    ##
    # Check the specified data set is actually available
    scn_names = scn.all_dataset_names()
    if dataset not in scn_names:
        raise Exception("Specified dataset is not available.")

    ##
    # Load the data, different calibration can be chosen
    scn.load([dataset], calibration=calibration)

    ##
    # Extract the longitude and latitude data
    lons, lats = scn[dataset].area.get_lonlats()

    ##
    # Apply a swath definition for our output raster
    swath_def = pr.geometry.SwathDefinition(lons=lons, lats=lats)

    ##
    # Extract the data values
    values = scn[dataset].values

    ##
    # Change the datatype of the arrays depending on the present data this can be changed
    lons = lons.astype(dtype)
    lats = lats.astype(dtype)
    values = values.astype(dtype)

    ##
    # Resample our data to the area of interest
    values = pr.kd_tree.resample_nearest(swath_def, values,
                                             area_def,
                                             radius_of_influence=radius, # in meters
                                             epsilon=epsilon,
                                             fill_value=False)

    ##
    # Check if outdir exists, otherwise create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ##
    # Join our filename based on the input file's basename
    outname = os.path.join(outdir, os.path.basename(file)[:-4] + "_" + str(label) + ".tif")
    outname = outname.replace(".tif", '_spain.tif' if to_spain else '.tif')
    if to_cog_tif:
       outname.replace(".tif", "_cog.tif")

    ##
    # Define some metadata
    cols = values.shape[1]
    rows = values.shape[0]
    pixelWidth = (area_def.area_extent[2] - area_def.area_extent[0]) / cols
    pixelHeight = (area_def.area_extent[1] - area_def.area_extent[3]) / rows
    originX = area_def.area_extent[0]
    originY = area_def.area_extent[3]

    ##
    # Create the file
    if to_cog_tif:
        driver = gdal.GetDriverByName("COG")
    else:
        driver = gdal.GetDriverByName("GTiff")
    print(driver)
    outRaster = driver.Create(outname, cols, rows, 1)
    print(dir(outRaster))

    ##
    # Save the metadata
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

    ##
    # Create a new band and write the data
    outband = outRaster.GetRasterBand(1)
    outband.SetNoDataValue(nodata) #specified no data value by user
    outband.WriteArray(np.array(values)) # writing the values
    outRasterSRS = osr.SpatialReference() # create CRS instance
    outRasterSRS.ImportFromEPSG(4326) # get info for EPSG 4326
    outRaster.SetProjection(outRasterSRS.ExportToWkt()) # set CRS as WKT

    ##
    # Clean up
    outband.FlushCache()
    outband = None
    outRaster = None

    return
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
