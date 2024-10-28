#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Read/Convert Meteosat Second Generation (MSG) Native Archive Format (.nat) file to GTiff.

nat2tif.py
~~~~~~~~~~~~
"""
# Standard Imports
import os
import copy
from typing import Union

# Third-Party Imports
import numpy as np
from osgeo import gdal
from osgeo import osr
import pyresample as pr
#   conda install -c conda-forge satpy
from satpy import Scene
import cartopy.crs as ccrs

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
def nat2tif(file: str, calibration: str, area_def: pr.geometry.AreaDefinition, dataset: str, reader: str, outdir: str, label: str, dtype: str, radius: int, epsilon: float, nodata: float, use_epsg: int, merc_vals: tuple[Union[None, int], Union[None, pr.geometry.AreaDefinition]], to_spain: bool):
    verbose = [False, True][0]
    to_merc = False if merc_vals[0] is None else True

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
    #   lons (11136, 5568): [-65.47019634595965, ... 81.21912937778]
    #   lats (11136, 5568): [-81.13614609568566, ... 81.13614609568566]
    lons, lats = scn[dataset].area.get_lonlats()
    if verbose:
        print("\nCCAST Sterographic")
        tmp = lons.flatten()
        tmp = tmp[np.abs(tmp) <= 180.0]
        print(f"\tlons {lons.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        tmp = lats.flatten()
        tmp = tmp[np.abs(tmp) <= 90.0]
        print(f"\tlats {lats.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        del tmp

    ##
    # Apply a swath definition for our output raster
    #   <class 'pyresample.geometry.SwathDefinition'>
    swath_def = pr.geometry.SwathDefinition(lons=lons, lats=lats)

    ##
    # Extract the data values
    #   values (11136, 5568): [0.0, ... 26.103036880493164]
    values = scn[dataset].values
    if verbose:
        tmp = values.flatten()
        tmp = tmp[~np.isnan(tmp)]
        print(f"\tvalues {values.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        del tmp

    ##
    # Change the datatype of the arrays depending on the present data this can be changed
    lons = lons.astype(dtype)
    lats = lats.astype(dtype)
    values = values.astype(dtype)
    if to_merc:
        values_orig = copy.deepcopy(values)

    ##
    # Resample our data to the area of interest
    #   values* (768, 768): [0.5235565900802612, ... 13.425487518310547]
    values = pr.kd_tree.resample_nearest(swath_def, values,
                                             area_def,
                                             radius_of_influence=radius, # in meters
                                             epsilon=epsilon,
                                             fill_value=False)
    if verbose:
        tmp = values.flatten()
        tmp = tmp[~np.isnan(tmp)]
        print(f"\tvalues* {values.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        del tmp

    ##
    # Check if outdir exists, otherwise create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ##
    # Join our filename based on the input file's basename
    outname = os.path.join(outdir, os.path.basename(file)[:-4] + "_" + str(label) + ".tif")
    outname = outname.replace(".tif", '_spain.tif' if to_spain else '.tif')

    ##
    # Define some metadata
    #   pixelWidth : 3000.0
    #   pixelHeight: -3000.0
    #   originX    : -855100.436345
    #   originY    : -2638000.0
    cols = values.shape[1]
    rows = values.shape[0]
    pixelWidth = (area_def.area_extent[2] - area_def.area_extent[0]) / cols
    pixelHeight = (area_def.area_extent[1] - area_def.area_extent[3]) / rows
    originX = area_def.area_extent[0]
    originY = area_def.area_extent[3]
    if verbose:
        print(f"\tpixelWidth : {pixelWidth}")
        print(f"\tpixelHeight: {pixelHeight}")
        print(f"\toriginX    : {originX}")
        print(f"\toriginY    : {originY}")

    ##
    # Create the file
    driver = gdal.GetDriverByName("GTiff")
    #   outRaster : <osgeo.gdal.Dataset>
    outRaster = driver.Create(outname, cols, rows, 1)

    ##
    # Save the metadata
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

    ##
    # Create a new band and write the data
    #   <class 'osgeo.gdal.Dataset'>
    outband = outRaster.GetRasterBand(1)
    outband.SetNoDataValue(nodata) #specified no data value by user
    outband.WriteArray(np.array(values)) # writing the values
    outRasterSRS = osr.SpatialReference() # create CRS instance
    outRasterSRS.ImportFromEPSG(use_epsg) # get info for EPSG 4326
    outRaster.SetProjection(outRasterSRS.ExportToWkt()) # set CRS as WKT

    ##
    # Clean up
    outband.FlushCache()
    outband = None
    outRaster = None

    """
    CCAST Stereographic
        lons (11136, 5568)  : [-65.47019634595965, ... 81.21912937778]
        lats (11136, 5568)  : [-81.13614609568566, ... 81.13614609568566]
        values (11136, 5568): [0.0, ... 26.103036880493164]
        values* (768, 768)  : [0.5235565900802612, ... 13.425487518310547]
        pixelWidth          : 3000.0
        pixelHeight         : -3000.0
        originX             : -855100.436345
        originY             : -2638000.0
    """

    if to_merc:
        ##
        # Make/Re-project to Mercator
        MERC_EPSG, merc_area_def = merc_vals

        ##
        # Join our filename based on the input file's basename
        outname = outname.replace(".tif", '_merc.tif')

        if verbose:
            print("\nCCAST Mercator")
            tmp = values_orig.flatten()
            tmp = tmp[~np.isnan(tmp)]
            print(f"\tvalues_orig {values_orig.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            del tmp

        ##
        # Resample our data to the area of interest
        #   values* (768, 768): [0.5609534978866577, ... 14.248218536376953]
        values = pr.kd_tree.resample_nearest(swath_def, values_orig,
                                                 merc_area_def,
                                                 radius_of_influence=radius,
                                                 epsilon=epsilon,
                                                 fill_value=False)
        if verbose:
            tmp = values.flatten()
            tmp = tmp[~np.isnan(tmp)]
            print(f"\tvalues* {values.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            del tmp

        ##
        # Define some metadata
        cols = values.shape[1]
        rows = values.shape[0]
        pixelWidth = (merc_area_def.area_extent[2] - merc_area_def.area_extent[0]) / cols
        pixelHeight = (merc_area_def.area_extent[1] - merc_area_def.area_extent[3]) / rows
        originX = merc_area_def.area_extent[0]
        originY = merc_area_def.area_extent[3]
        if verbose:
            print(f"pixelWidth : {pixelWidth}")
            print(f"pixelHeight: {pixelHeight}")
            print(f"originX    : {originX}")
            print(f"originY    : {originY}")

        ##
        # Create the file
        driver = gdal.GetDriverByName("GTiff")
        outRaster = driver.Create(outname, cols, rows, 1)

        ##
        # Save the metadata
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

        ##
        # Create a new band and write the data
        #   <class 'osgeo.gdal.Dataset'>
        outband = outRaster.GetRasterBand(1)
        outband.SetNoDataValue(nodata) #specified no data value by user
        outband.WriteArray(np.array(values)) # writing the values
        outRasterSRS = osr.SpatialReference() # create CRS instance
        outRasterSRS.ImportFromEPSG(MERC_EPSG) # get info for EPSG 4326
        outRaster.SetProjection(outRasterSRS.ExportToWkt()) # set CRS as WKT

        ##
        # Clean up
        outband.FlushCache()
        outband = None
        outRaster = None

        """
        CCAST Mercator
            values_orig (11136, 5568): [0.0, ... 26.103036880493164]
            values* (768, 768): [0.5609534978866577, ... 14.248218536376953]
            pixelWidth : 5594.088533248166
            pixelHeight: -4241.094870654907
            originX    : -536174.1912289965
            originY    : 8397504.685448818
        """

    return
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
