#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Read Meteosat Second Generation (MSG) Native Archive Format (.nat) file.

readnat.py
~~~~~~~~~~
"""
# Standard Imports
import os
import typing
from typing import Optional, Union

# Third-Party Imports
import numpy as np
import numpy.typing as npt
# from osgeo import gdal
# from osgeo import osr
# import pyresample as pr
#   conda install -c conda-forge satpy
from satpy import Scene

# STARE Imports

# Local Imports

##
# List of Public objects from this module.
__all__ = ['readnat']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

###############################################################################
# PUBLIC readnat()
# ----------------
def readnat(file: str, calibration: str, dataset: Union[str, list], reader: str, dtype: str) -> tuple[npt.ArrayLike, npt.ArrayLike, npt.ArrayLike]:
    ##
    # Open the file w/ satpy, which uses Xarray
    #   <class 'satpy.scene.Scene'>
    #   ['aggregate', 'all_composite_ids', 'all_composite_names', 'all_dataset_ids', 'all_dataset_names', 'all_modifier_names',
    #    'all_same_area', 'all_same_proj', 'attrs', 'available_composite_ids', 'available_composite_names', 'available_dataset_ids',
    #    'available_dataset_names', 'chunk', 'coarsest_area', 'compute', 'copy', 'crop', 'end_time', 'finest_area',
    #    'generate_possible_composites', 'get', 'images', 'iter_by_area', 'keys', 'load', 'max_area', 'min_area', 'missing_datasets',
    #    'persist', 'resample', 'save_dataset', 'save_datasets', 'sensor_names', 'show', 'slice', 'start_time', 'to_geoviews',
    #    'to_hvplot', 'to_xarray', 'to_xarray_dataset', 'unload', 'values', 'wishlist']
    scn = Scene(filenames = {reader: [file]})

    ##
    # Check the specified data set is actually available
    scn_names = scn.all_dataset_names()
    single_var = True if isinstance(dataset, str) else False
    if single_var:
        if dataset not in scn_names:
            raise Exception(f"Specified dataset {dataset} is not available.")
    else:
        for ds in dataset:
            if ds not in scn_names:
                raise Exception(f"Specified dataset {ds} is not available.")

    ##
    # Load the data, different calibration can be chosen
    scn.load([dataset] if single_var else dataset, calibration=calibration)

    ##
    # Extract the longitude and latitude data
    lons, lats = scn[dataset if single_var else dataset[0]].area.get_lonlats()

    ##
    # Extract the data values
    if single_var:
        values = scn[dataset].values
    else:
        # values = scn.to_xarray_dataset(datasets=dataset)
        values = []
        for ds in dataset:
            values.append(scn[ds].values)

    ##
    # Change the datatype of the arrays depending on the present data this can be changed
    #   lons   (11136, 5568)
    #   lats   (11136, 5568)
    #   values (11136, 5568) or (4, 3712, 3712)
    lons = lons.astype(dtype)
    lats = lats.astype(dtype)
    if single_var:
        values = values.astype(dtype)
    else:
        tmp_shp = values[0].shape
        new_values = np.zeros((len(dataset), tmp_shp[0], tmp_shp[1]), dtype=dtype)
        for didx in range(len(dataset)):
            new_values[didx, :, :] = values[didx][:, :].astype(dtype)
        values = new_values
        del new_values

    # print(f"lons {lons.shape}")
    # print(f"lats {lats.shape}")
    # print(f"values {values.shape}")

    return lons, lats, values
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
