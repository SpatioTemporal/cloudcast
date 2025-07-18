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
import typing
from typing import Optional, Union, TypeAlias

# Third-Party Imports
import numpy as np
from osgeo import gdal
from osgeo import osr
import pyresample as pr
from satpy import Scene
import cartopy
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

##
# Type Aliases
SetOut: TypeAlias = tuple[int, int, int, typing.TypeVar('cartopy.crs'),
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')]]

###############################################################################
# PUBLIC nat2tif()
# ----------------
def nat2tif(fname: str, fvar: str, reader: str, outdir: str, label: str, atag: str, geo_dat: SetOut, to_full:bool, to_ccast: bool, to_euro: bool, to_merc: bool, to_lcc: bool) -> None:
    """
    Read/Convert Meteosat Second Generation (MSG) Native Archive Format (.nat) file to GTiff.

    Parameters
    ----------
    fname : str
        Full path to the input file
    fvar : str
        The variable to extract from the file
    reader : str
        The reader to use to open the file
    outdir : str
        The directory to save the output file
    label : str
        The label to add to the output file
    atag : str
        The tag to add to the output file
    geo_dat : SetOut
        The geospatial data to use for the conversion
    to_full : bool
        Return MSG full frame data (no reprojection, subsetting or interpolation)
    to_ccast : bool
        Return the CloudCast resolution and domain (reprojection, subsetting and interpolation)
    to_euro : bool
        Return the CloudCast european resolution and domain (reprojection, subsetting and interpolation)
    to_merc : bool
        Return in Mercator projection (reprojection, subsetting and interpolation)
    to_lcc : bool
        Return in Lambert Conformal projection (reprojection, subsetting and interpolation)
    """
    verbose = [False, True][1]

    # to_euro
    euro_nrows = 928
    euro_ncols = 1530

    ##
    # Unpack setup setup data
    (MSG_EPSG, MERC_EPSG, LCC_EPSG, geod_crs, ccast_area_def, ccast_crs, ccast_merc_area_def,
     ccast_merc_crs, ccast_lcc_area_def, ccast_lcc_crs, msg_area_def, msg_crs,
     msg_merc_area_def, msg_merc_crs, raw_area_def, raw_crs, raw_merc_area_def,
     raw_merc_crs) = geo_dat

    use_lcc_crs = None
    use_lcc_area_def = None
    if to_full:
        use_crs = msg_crs
        use_area_def = msg_area_def
        use_merc_crs = msg_merc_crs
        use_merc_area_def = msg_merc_area_def
    elif to_euro:
        use_crs = raw_crs
        use_area_def = raw_area_def
        use_merc_crs = raw_merc_crs
        use_merc_area_def = raw_merc_area_def
    elif to_ccast:
        use_crs = ccast_crs
        use_area_def = ccast_area_def
        use_merc_crs = ccast_merc_crs
        use_merc_area_def = ccast_merc_area_def
        use_lcc_crs = ccast_lcc_crs
        use_lcc_area_def = ccast_lcc_area_def
    dtype = "float32"
    radius = 16000
    epsilon = 0.5
    nodata = -3.4E+38

    ##
    # Open the file w/ satpy, which uses Xarray
    #   <class 'satpy.scene.Scene'>
    scn = Scene(filenames = {reader: [fname]})

    ##
    # Check the specified data set is actually available
    scn_names = scn.all_dataset_names()
    if fvar not in scn_names:
        raise Exception(f"Specified variable {fvar} is not available.")

    ##
    # Get calibration
    # if fvar in ['VIS006', "HRV"]:
    #     calibration = "reflectance"
    # else:
    #     calibration = "brightness_temperature"
    # calibration = "radiance"
    """
    IR_039 calibration = "radiance"
        data_vals (928, 1530): [0.003658667206764221, ... 1.4963948726654053]

    IR_039 calibration = "", which I think defaults to "brightness_temperature"
        data_vals (928, 1530): [204.75961303710938, ... 310.7130126953125]
    """
    scn.load([fvar])
    fvar__ = scn[fvar]
    fvar__atts = fvar__.attrs
    fvar_units = fvar__atts['units']
    fvar_name = fvar__atts['standard_name']
    use_calibration = fvar__atts['calibration']
    if verbose:
        print(f"\tcalibration '{use_calibration}' w/ units of {fvar_units}")
    del fvar__

    ##
    # Load the data, different calibration can be chosen
    scn.load([fvar], calibration=use_calibration)
    # scn.load([fvar])

    ##
    # Extract the longitude and latitude data
    """
    if fvar == HRV
        lons (11136, 5568)     : [-65.47019958496094, ... 81.21913146972656]
        lats (11136, 5568)     : [-81.13614654541016, ... 81.13614654541016]
    else:
        lons (3712, 3712)      : [-81.12566375732422, ... 81.12566375732422]
        lats (3712, 3712)      : [-81.0744857788086,  ... 81.0744857788086]
    """
    lons, lats = scn[fvar].area.get_lonlats()
    if to_euro:
        ##
        # MSG subset (928, 1530)
        lons = lons[-928:, 1092:2622]
        lats = lats[-928:, 1092:2622]
    if verbose:
        if to_full or to_euro:
            print(f"\nMSG LatLon {fvar}")
        else:
            print(f"\nCCAST Stereographic {fvar}")
        tmp = lons.flatten()
        tmp = tmp[np.abs(tmp) <= 180.0]
        print(f"\tlons {lons.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        tmp = lats.flatten()
        tmp = tmp[np.abs(tmp) <= 90.0]
        print(f"\tlats {lats.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        del tmp

    ##
    # Extract the data values
    """
    if fvar == HRV
        data_vals (11136, 5568): [0.0, ... 26.103036880493164]
    else:
        data_vals (3712, 3712) : [0.0, ... 3.5562241077423096]
    """
    data_vals = scn[fvar].values
    if to_euro:
        data_vals = data_vals[-928:, 1092:2622]
    if verbose:
        tmp = data_vals.flatten()
        tmp = tmp[~np.isnan(tmp)]
        print(f"\tdata_vals {data_vals.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        del tmp

    # os._exit(1)

    ##
    # Change the datatype of the arrays depending on the present data this can be changed
    lons = lons.astype(dtype)
    lats = lats.astype(dtype)
    data_vals = data_vals.astype(dtype)
    if to_merc or to_lcc:
        data_vals_orig = copy.deepcopy(data_vals)

    if to_ccast:
        ##
        # Apply a swath definition for our output raster
        #   <class 'pyresample.geometry.SwathDefinition'>
        swath_def = pr.geometry.SwathDefinition(lons=lons, lats=lats)

        ##
        # Resample our data to the area of interest
        #   data_vals* (768, 768): [0.5235565900802612, ... 13.425487518310547]
        data_vals = pr.kd_tree.resample_nearest(swath_def, data_vals,
                                                 use_area_def,
                                                 radius_of_influence=radius, # in meters
                                                 epsilon=epsilon,
                                                 fill_value=False)
        if verbose:
            tmp = data_vals.flatten()
            tmp = tmp[~np.isnan(tmp)]
            print(f"\tdata_vals* {data_vals.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            del tmp

        # data_vals (3712, 3712): [0.0, ... 3.5562241077423096]
        # data_vals* (768, 768): [0.02926933765411377, ... 0.8268587589263916]

    ##
    # Check if outdir exists, otherwise create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ##
    # Join our filename based on the input file's basename
    outname = os.path.join(outdir, os.path.basename(fname)[:-4] + "_" + str(label) + f"_{atag}.tif")

    ##
    # Define some metadata
    #   pixelWidth : 3000.0
    #   pixelHeight: -3000.0
    #   originX    : -855100.436345
    #   originY    : -2638000.0
    cols = data_vals.shape[1]
    rows = data_vals.shape[0]
    pixelWidth = (use_area_def.area_extent[2] - use_area_def.area_extent[0]) / cols
    pixelHeight = (use_area_def.area_extent[1] - use_area_def.area_extent[3]) / rows
    originX = use_area_def.area_extent[0]
    originY = use_area_def.area_extent[3]
    if verbose:
        print(f"\tpixelWidth : {pixelWidth}")
        print(f"\tpixelHeight: {pixelHeight}")
        print(f"\toriginX    : {originX}")
        print(f"\toriginY    : {originY}")

    ##
    # Create the file
    driver = gdal.GetDriverByName("GTiff")
    #   outRaster : <osgeo.gdal.var>
    outRaster = driver.Create(outname, cols, rows, 1)

    ##
    # Save the metadata
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))

    ##
    # Create a new band and write the data
    #   <class 'osgeo.gdal.var'>
    outband = outRaster.GetRasterBand(1)
    outband.SetNoDataValue(nodata) #specified no data value by user
    outband.WriteArray(np.array(data_vals)) # writing the values
    outRasterSRS = osr.SpatialReference() # create CRS instance
    outRasterSRS.ImportFromEPSG(MSG_EPSG) # get info for EPSG 4326
    outRaster.SetProjection(outRasterSRS.ExportToWkt()) # set CRS as WKT

    # Band 1 Block=1530x5 Type=Byte, ColorInterp=Gray
    #   NoData Value=-3.4e+38

    ##
    # Clean up
    outband.FlushCache()
    outband = None
    outRaster = None

    """
    MSG LatLon IR_039
        lons (3712, 3712): [-81.12566252856628, ... 81.12566252856628]
        lats (3712, 3712): [-81.07448485965803, ... 81.07448485965803]
        data_vals (3712, 3712): [0.0, ... 3.5562241077423096]
        pixelWidth : 1.3787170936321391e-05
        pixelHeight: 1.3787170936321391e-05
        originX    : -81.12566375732422
        originY    : 81.0744857788086

    MSG LatLon HRV
        lons (11136, 5568): [-65.47019634595965, ... 81.21912937778]
        lats (11136, 5568): [-81.13614609568566, ... 81.13614609568566]
        data_vals (11136, 5568): [0.0, ... 26.103036880493164]
        pixelWidth : -0.0028135680604255064
        pixelHeight: 7.451950818642802e-06
        originX    : -65.47019958496094
        originY    : 81.13614654541016

    CCAST Stereographic
        lons (11136, 5568)  : [-65.47019634595965, ... 81.21912937778]
        lats (11136, 5568)  : [-81.13614609568566, ... 81.13614609568566]
        data_vals (11136, 5568): [0.0, ... 26.103036880493164]
        data_vals* (768, 768)  : [0.5235565900802612, ... 13.425487518310547]
        pixelWidth          : 3000.0
        pixelHeight         : -3000.0
        originX             : -855100.436345
        originY             : -2638000.0
    """

    if to_merc:
        ##
        # Make/Re-project to Mercator

        ##
        # Join our filename based on the input file's basename
        outname = outname.replace(".tif", '_merc.tif')

        if verbose:
            print("\nMercator")
            tmp = data_vals_orig.flatten()
            tmp = tmp[~np.isnan(tmp)]
            print(f"\tdata_vals_orig {data_vals_orig.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            del tmp

        ##
        # Resample our data to the area of interest
        #   data_vals* (768, 768): [0.5609534978866577, ... 14.248218536376953]
        data_vals = pr.kd_tree.resample_nearest(swath_def, data_vals_orig,
                                                 use_merc_area_def,
                                                 radius_of_influence=radius,
                                                 epsilon=epsilon,
                                                 fill_value=False)
        if verbose:
            tmp = data_vals.flatten()
            tmp = tmp[~np.isnan(tmp)]
            print(f"\tdata_vals* {data_vals.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            del tmp

        ##
        # Define some metadata
        cols = data_vals.shape[1]
        rows = data_vals.shape[0]
        pixelWidth = (msg_merc_area_def.area_extent[2] - msg_merc_area_def.area_extent[0]) / cols
        pixelHeight = (msg_merc_area_def.area_extent[1] - msg_merc_area_def.area_extent[3]) / rows
        originX = msg_merc_area_def.area_extent[0]
        originY = msg_merc_area_def.area_extent[3]
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
        #   <class 'osgeo.gdal.var'>
        outband = outRaster.GetRasterBand(1)
        outband.SetNoDataValue(nodata) #specified no data value by user
        outband.WriteArray(np.array(data_vals)) # writing the values
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
            data_vals_orig (11136, 5568): [0.0, ... 26.103036880493164]
            data_vals* (768, 768): [0.5609534978866577, ... 14.248218536376953]
            pixelWidth : 5594.088533248166
            pixelHeight: -4241.094870654907
            originX    : -536174.1912289965
            originY    : 8397504.685448818
        """

    if to_lcc and to_ccast:
        ##
        # Make/Re-project to Lambert Conformal

        ##
        # Join our filename based on the input file's basename
        outname = outname.replace(".tif", '_lcc.tif')

        if verbose:
            print("\nLambert Conformal")
            tmp = data_vals_orig.flatten()
            tmp = tmp[~np.isnan(tmp)]
            print(f"\tdata_vals_orig {data_vals_orig.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            del tmp

        ##
        # Resample our data to the area of interest
        #   data_vals* (768, 768): [0.02561067044734955, ... 3.5562241077423096]
        data_vals = pr.kd_tree.resample_nearest(swath_def, data_vals_orig,
                                                 use_lcc_area_def,
                                                 radius_of_influence=radius,
                                                 epsilon=epsilon,
                                                 fill_value=False)
        if verbose:
            tmp = data_vals.flatten()
            tmp = tmp[~np.isnan(tmp)]
            print(f"\tdata_vals* {data_vals.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            del tmp

        ##
        # Define some metadata
        cols = data_vals.shape[1]
        rows = data_vals.shape[0]
        pixelWidth = (use_lcc_area_def.area_extent[2] - use_lcc_area_def.area_extent[0]) / cols
        pixelHeight = (use_lcc_area_def.area_extent[1] - use_lcc_area_def.area_extent[3]) / rows
        originX = use_lcc_area_def.area_extent[0]
        originY = use_lcc_area_def.area_extent[3]
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
        #   <class 'osgeo.gdal.var'>
        outband = outRaster.GetRasterBand(1)
        outband.SetNoDataValue(nodata) #specified no data value by user
        outband.WriteArray(np.array(data_vals)) # writing the values
        outRasterSRS = osr.SpatialReference() # create CRS instance
        outRasterSRS.ImportFromEPSG(LCC_EPSG)
        outRaster.SetProjection(outRasterSRS.ExportToWkt()) # set CRS as WKT

        ##
        # Clean up
        outband.FlushCache()
        outband = None
        outRaster = None

        """
        CCAST Lambert Conformal
            data_vals_orig (3712, 3712): [0.0, ... 3.5562241077423096]
            data_vals* (768, 768):       [0.02561067044734955, ... 3.5562241077423096]
        pixelWidth : 3412.1674980195094
        pixelHeight: -2796.642223159081
        originX    : -1216087.7976068966
        originY    : 1248209.190564286
        """

    return
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
