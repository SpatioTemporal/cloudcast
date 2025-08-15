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
from pathlib import Path
import pickle

# Third-Party Imports
import numpy as np
import numpy.typing as npt
# from osgeo import gdal
# from osgeo import osr
import pyresample as pr
#   conda install -c conda-forge satpy
from satpy import Scene
import cartopy
import cartopy.crs as ccrs
from zipfile import ZipFile


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
# PUBLIC natread()
# ----------------
def natread(fname: str, fvar: str, reader: str, to_euro: bool, euro_lons: npt.ArrayLike, euro_lats: npt.ArrayLike, to_ccast:bool, geo_dat: SetOut, fromzip: bool) -> None:

    verbose = [False, True][1]

    ##
    # Unpack setup setup data
    (MSG_EPSG, MERC_EPSG, LCC_EPSG, geod_crs, ccast_area_def, ccast_crs, ccast_merc_area_def,
     ccast_merc_crs, ccast_lcc_area_def, ccast_lcc_crs, msg_area_def, msg_crs,
     msg_merc_area_def, msg_merc_crs, raw_area_def, raw_crs, raw_merc_area_def,
     raw_merc_crs) = geo_dat

    use_lcc_crs = None
    use_lcc_area_def = None
    if to_euro:
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
    if verbose:
        print(f"\nReading {fvar}")

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

    if fromzip:
        zpath = Path(fname)
        zfile = f"{zpath.parent}.zip"
        ffile = f"{zpath.name}"
        with ZipFile(zfile, 'r') as zipper:
            fname = zipper.extract(ffile)
    scn = Scene(filenames = {reader:[fname]})
    # print(scn.values())

    ##
    # Get calibration
    scn.load([fvar])
    fvar__ = scn[fvar]
    fvar__atts = fvar__.attrs
    fvar_units = fvar__atts['units']
    fvar_name = fvar__atts['standard_name']
    use_calibration = fvar__atts['calibration']

    # print(fvar__atts['orbital_parameters'])
    print(fvar__atts.keys())

    # print(f"{fvar__atts['platform_name'] = }")
    # print(f"{fvar__atts['sensor'] = }")
    # print(f"{fvar__atts['time_parameters'] = }")
    # os._exit(0)

    if verbose:
        print(f"\tcalibration '{use_calibration}' w/ units of {fvar_units}")
    del fvar__

    ## use_calibration = "counts"
    ## use_calibration = "radiance"

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
    lons, lats, data_vals = readnat(file=fname, calibration=use_calibration, dataset=fvar, reader=reader, dtype="float32")
    if verbose:
        tmp = lons.flatten()
        tmp_len = len(tmp)
        tmp = tmp[np.abs(tmp) <= 180.0]
        tmp1_len = len(tmp)
        len_frac = 100.0 * (tmp1_len / tmp_len)
        print(f"\n\tlons {lons.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        tmp = lats.flatten()
        tmp_len = len(tmp)
        tmp = tmp[np.abs(tmp) <= 90.0]
        tmp1_len = len(tmp)
        len_frac = 100.0 * (tmp1_len / tmp_len)
        print(f"\tlats {lats.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        tmp = data_vals.flatten()
        tmp_len = len(tmp)
        tmp = tmp[~np.isnan(tmp)]
        tmp = tmp[tmp > 0.0]
        if use_calibration == "reflectance":
            tmp = tmp[tmp <= 100.0]
        tmp1_len = len(tmp)
        len_frac = 100.0 * (tmp1_len / tmp_len)
        print(f"\tdata_vals {data_vals.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        del tmp

    if fromzip:
        os.remove(fname)

    if to_ccast:
        dtype = "float32"
        radius = 16000
        epsilon = 0.5
        nodata = -3.4E+38

        ##
        # Apply a swath definition for our output raster
        #   <class 'pyresample.geometry.SwathDefinition'>
        swath_def = pr.geometry.SwathDefinition(lons=lons, lats=lats)

        ##
        # Resample our data to the area of interest
        #   data_vals* (768, 768): [0.5235565900802612, ... 13.425487518310547]

        # if use_calibration == "reflectance":
        #     data_vals = np.nan_to_num(data_vals, nan=0.0)

        ccast_data_vals = pr.kd_tree.resample_nearest(swath_def, data_vals,
                                                      use_area_def,
                                                      radius_of_influence=radius, # in meters
                                                      epsilon=epsilon,
                                                      fill_value=False)
        ccast_lons = pr.kd_tree.resample_nearest(swath_def, lons,
                                                 use_area_def,
                                                 radius_of_influence=radius, # in meters
                                                 epsilon=epsilon,
                                                 fill_value=False)
        ccast_lats = pr.kd_tree.resample_nearest(swath_def, lats,
                                                 use_area_def,
                                                 radius_of_influence=radius, # in meters
                                                 epsilon=epsilon,
                                                 fill_value=False)

        ##
        # Save to File for working with ERA-5
        lonlatfile = "/Users/mbauer/tmp/CloudCast/ccast_lonlat.npy"
        with open(lonlatfile, 'wb') as f:
            np.save(f, ccast_lons)
            np.save(f, ccast_lats)
        if verbose:
            tmp = ccast_lons.flatten()
            tmp_len = len(tmp)
            tmp = tmp[np.abs(tmp) <= 180.0]
            tmp1_len = len(tmp)
            len_frac = 100.0 * (tmp1_len / tmp_len)
            print(f"\n\tccast_lons {ccast_lons.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")

            tmp = ccast_lats.flatten()
            tmp_len = len(tmp)
            tmp = tmp[np.abs(tmp) <= 90.0]
            tmp1_len = len(tmp)
            len_frac = 100.0 * (tmp1_len / tmp_len)
            print(f"\tccast_lats {ccast_lats.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")

            tmp = ccast_data_vals.flatten()
            tmp_len = len(tmp)
            tmp = tmp[~np.isnan(tmp)]
            tmp = tmp[tmp > 0.0]
            if use_calibration == "reflectance":
                tmp = tmp[tmp <= 100.0]
            tmp1_len = len(tmp)
            len_frac = 100.0 * (tmp1_len / tmp_len)
            print(f"\tccast_data_vals {ccast_data_vals.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")

            print(np.std(tmp))  # 13.81

            tmp = tmp.astype(np.uint16)
            print(f"\tccast_data_vals {ccast_data_vals.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
            print(np.std(tmp))  # 13.813447
            # ccast_data_vals (768, 768) 100.00%: [210, ... 292]
            # 13.83

            del tmp


        return


    if to_euro:
        #   euro_lons (928, 1530): [-69.2706298828125   ... 69.2706298828125]
        #   euro_lats (928, 1530): [ 26.67105484008789  ... 81.09877014160156]
        #
        #   lons (3712, 3712)    : [-81.12566375732422, ... 81.12566375732422]
        #   lats (3712, 3712)    : [-81.0744857788086,  ... 81.0744857788086]
        #
        #   fvar == HRV
        #   lons (11136, 5568)   : [-65.47019958496094, ... 81.21913146972656]
        #   lats (11136, 5568)   : [-81.13614654541016, ... 81.13614654541016]
        """
            euro_lats (928, 1530) vs lats (3712, 3712)
        """
        if fvar == 'HRV':
            pass
        else:
            ##
            # Look for i/j (col/row) alignment
            r"""
            euro_lats organized so that they start in North and move South, while euro_lons start in the West and move East

                euro_lats 0000: NAN
                ...
                euro_lats 0050: NAN
                euro_lats 0051: i_start = 0741   +81.0988 ... i_end = 0789   +81.0988
                ...
                euro_lats 0927: i_start = 0000   +27.1626 ... i_end = 1529   +27.1612

                euro_lons 0000: NAN
                ...
                euro_lons 0050: NAN
                euro_lons 0051: i_start = 0741    -4.8004 ... i_end = 0789    +4.8004
                ...
                euro_lons 0927: i_start = 0000   -24.7438 ... i_end = 1529   +24.7077

            lats organized so that they start in South and move North, while

                lats 0000: NAN
                ...
                lats 0051: i_start = 1808   -81.0745 ... i_end = 1903   -81.0745
                ...
                lats 1855: i_start = 0045    -0.0157 ... i_end = 3666    -0.0157
                ....
                lats 1856: i_start = 0045    +0.0157 ... i_end = 3666    +0.0157    Start N Hemi
                ...
                lats 2652: i_start = 0229   +26.0090 ... i_end = 3482   +26.0090    Approaching Southern edge of euro_lats
                ...
                lats 3660: i_start = 1808   +81.0745 ... i_end = 1903   +81.0745
                lats 3661: NAN
                ...
                lats 3711: NAN

                lons 0000: NAN
                ...
                lons 0050: NAN
                lons 0051: i_start = 1808    +9.5090 ... i_end = 1903    -9.5090
                ...
                lons 3660: i_start = 1808    +9.5090 ... i_end = 1903    -9.5090
                lons 3661: NAN
                ...
                lons 3711: NAN

            Need to reverse order of both euro_lats and euro_lons

                euro_lats 0000: i_start = 0000   +27.1612 ... i_end = 1529   +27.1626
                ...
                euro_lats 0876: i_start = 0740   +81.0988 ... i_end = 0788   +81.0988
                euro_lats 0877: NAN
                ...
                euro_lats 0927: NAN

                euro_lons 0000: i_start = 0000   +24.7077 ... i_end = 1529   -24.7438
                ...
                euro_lons 0876: i_start = 0740    +4.8004 ... i_end = 0788    -4.8004
                euro_lons 0877: NAN
                ...
                euro_lons 0927: NAN
            """
            euro_lats = euro_lats[::-1, ::-1]
            euro_lons = euro_lons[::-1, ::-1]
            # for jidx in range(euro_lons.shape[0]):
            #     msg = f"euro_lons {jidx:04d}:"
            #     tmp = euro_lons[jidx, :]

            #     # for euro_lats/euro_lons
            #     good_idx = np.nonzero(~np.isnan(tmp))[0]

            #     # for lats/lons
            #     # good_idx = np.nonzero(~np.isinf(tmp))[0]
            #     if len(good_idx) == 0:
            #         print(f"{msg} NAN")
            #         continue
            #     else:
            #         print(f"{msg} i_start = {good_idx[0]:04d} {tmp[good_idx[0]]:+10.4f} ... i_end = {good_idx[-1]:04d} {tmp[good_idx[-1]]:+10.4f}")
            # os._exit(1)

            ##
            # The minimum height is 928 rows of the starting 3712 rows
            #   so skip_to_bot = 2652 leaves enough (1060) for a fit.
            skip_to_bot = 2652

            #   check_hieght = 1060
            check_hieght = lats.shape[0] - 2652
            # print(f"{check_hieght = }")

            #   col_span = 1530
            col_span = euro_lats.shape[1]
            # print(f"{col_span = }")

            """
            raw_top_row 0927
            top_row     3712
            raw_bot_row 0000
            bot_row     2784

            test_lats (928, 3712)

            Min Lat Col 1092 w/ diff 40847800
                1092+1530 = 2622
                3712-2622 = 1090

            Starting with
                 raw_lats (928, 1530) vs ccast_lats (3712, 3712)
                 raw_lons (928, 1530) vs ccast_lons (3712, 3712)

            We want a subset of MSG_lats/MSG_lons that has the dimensions (928, 1530) and has the minimum absolute summed difference in lat and lon between
            the subset of MSG_lats/MSG_lons and raw_lats/raw_lons.

            To start I simply anchored the subset so that the top (highest latitude) MSG and raw rows match up (both ~81N).
            Then the subset can be extended down 928 rows to lower latitudes (both ~27N).

            This seems like a good start.

            Next, I found the summed absolute lat + lon difference between the subset and the target raw_lats/raw_lons.
            That is, I took a raw sized (928, 1530) slice from the MSG lon/lat (3712, 3712) and found the absolute lat + lon difference.
            I did this by sliding the subset from the left edge to the right edge of the MSG width (3712) which allows for a 1530 wide subset.
            Doing this I found the difference dropped consistently until column 1092, after which it climbed again.

            This is a happy solution as a 1530 wide subset starting at column 1092 happens to be almost exactly 1090 columns from each edge.

            Thus basically, it seems likely that the (928, 1530) from (3712, 3712) was done by taking the top 928 rows and then the middle 1530 columns.

            I'll clean this up and post the code (and check that better fits aren't to be found).

            occurs if we slide the subset
            over to very close to the middle (edges equally far from edges of the MSG grid 3712 wide) so

            1092 to 2622 is a slice of 1530 (correct) and very close to 1092 from the column edges (seems good)

                Columns 1092-2622:
                    Average abs lat difference is 0.0508 degrees
                    Average abs lon difference is 0.0588 degrees

            I found that we anchoring the subset so that the top (highest latitude) rows match up and then
            sliding a 1530 column wide and

            1092 to 2622 is a slice of 1530 (correct) and very close to 1092 from the column edges (seems good)

            """
            # raw_top_row = euro_lats.shape[0] - 1
            # top_row = lats.shape[0]
            # raw_bot_row = 0
            # bot_row = lats.shape[0] - euro_lats.shape[0]
            # print(f"raw_top_row {raw_top_row:04d}")
            # print(f"top_row     {top_row:04d}")
            # print(f"raw_bot_row {raw_bot_row:04d}")
            # print(f"bot_row     {bot_row:04d}")

            # test_lats = lats[bot_row:top_row, :]
            # test_lats = np.nan_to_num(test_lats, nan=0, posinf=0, neginf=0)
            # print(f"test_lats {test_lats.shape}")
            # test_euro_lats = np.nan_to_num(euro_lats, nan=0, posinf=0, neginf=0)

            # test_lons = lons[bot_row:top_row, :]
            # test_lons = np.nan_to_num(test_lons, nan=0, posinf=0, neginf=0)
            # print(f"test_lons {test_lons.shape}")
            # test_euro_lons = np.nan_to_num(euro_lons, nan=0, posinf=0, neginf=0)

            ##
            # MSG subset (928, 1530)
            lons_sub = lons[-928:, 1092:2622]
            lats_sub = lats[-928:, 1092:2622]
            data_vals_sub = data_vals[-928:, 1092:2622]
            if verbose:
                tmp = lons_sub.flatten()
                tmp_len = len(tmp)
                tmp = tmp[np.abs(tmp) <= 180.0]
                tmp1_len = len(tmp)
                len_frac = 100.0 * (tmp1_len / tmp_len)
                print(f"\n\tlons_euro {lons_sub.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
                tmp = lats_sub.flatten()
                tmp_len = len(tmp)
                tmp = tmp[np.abs(tmp) <= 90.0]
                tmp1_len = len(tmp)
                len_frac = 100.0 * (tmp1_len / tmp_len)
                print(f"\tlats_euro {lats_sub.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
                tmp = data_vals_sub.flatten()
                tmp_len = len(tmp)
                tmp = tmp[~np.isnan(tmp)]
                tmp = tmp[tmp > 0.0]
                if use_calibration == "reflectance":
                    tmp = tmp[tmp <= 100.0]
                tmp1_len = len(tmp)
                len_frac = 100.0 * (tmp1_len / tmp_len)
                print(f"\tdata_vals_euro {data_vals_sub.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
                del tmp
            return


    #         ##
    #         # Scan along the row (by column) looking for a col_span set of lat diffs
    #         ncols = lats.shape[1]
    #         min_diff_col = -1
    #         min_diff = 1e10
    #         for iidx in range(ncols):
    #             span_end = iidx + col_span
    #             if span_end > ncols - 1:
    #                 break
    #             test_span_lats = test_lats[:, iidx:span_end]
    #             test_span_lons = test_lons[:, iidx:span_end]
    #             lat_diff = np.sum(np.abs(np.subtract(test_euro_lats, test_span_lats)))
    #             lon_diff = np.sum(np.abs(np.subtract(test_euro_lons, test_span_lons)))
    #             total_diff = lat_diff + lon_diff
    #             if total_diff < min_diff:
    #                 min_diff = total_diff
    #                 min_diff_col = iidx
    #             # print(f"\t{iidx:04d}-{span_end}: {total_diff:10.1f} {np.mean(np.abs(np.subtract(test_euro_lats, test_span_lats)))} {np.mean(np.abs(np.subtract(test_euro_lons, test_span_lons)))}")
    #             # if iidx == 1090:
    #             #     print(f"\t{iidx:04d}-{span_end}: {total_diff:10.1f} {np.mean(np.abs(np.subtract(test_euro_lons[464, :], test_span_lons[464, :])))}")
    #             #     print("diff", np.abs(np.subtract(test_euro_lons[464, :], test_span_lons[464, :])))
    #             #     print(f"{test_euro_lons[464, :] = }")
    #             #     print(f"{test_span_lons[464, :] = }\n")
    #         print(f"Min Lat Col {min_diff_col} w/ diff {min_diff}")
    #         # os._exit(1)

    #         # # row_diff = np.zeros((check_hieght, lats.shape[0]))
    #         # row_diff = {}
    #         # for jidx in range(lats.shape[0]):
    #         #     if jidx < skip_to_bot:
    #         #         continue
    #         #     msg = f"\tlats {jidx:04d}:"
    #         #     tmp = lats[jidx, :]
    #         #     good_idx = np.nonzero(np.abs(tmp) <= 90.0)[0]
    #         #     if len(good_idx) == 0:
    #         #         print(f"{msg} NAN")
    #         #         continue
    #         #     else:
    #         #         print(f"{msg} i_start = {good_idx[0]:04d} {tmp[good_idx[0]]:+10.4f} ... i_end = {good_idx[-1]:04d} {tmp[good_idx[-1]]:+10.4f}")
    #         #         ##
    #         #         # Scan along the row (by column) looking for a col_span set of lat diffs
    #         #         for iidx in range(lats.shape[1]):
    #         #             test_span = lats[jidx, iidx:iidx + col_span]
    #         #             raw_span =
    #         #             # r_diff =
    #         #             os._exit(1)

    #     #   euro_lons (928, 1530): [-69.2706298828125   ... 69.2706298828125]
    #     #   euro_lats (928, 1530): [ 26.67105484008789  ... 81.09877014160156]
    #     #
    #     #   lons (3712, 3712)    : [-81.12566375732422, ... 81.12566375732422]
    #     #   lats (3712, 3712)    : [-81.0744857788086,  ... 81.0744857788086]

    #     # return
    #     # os._exit(1)

    #     # 3712 x 3712
    #     # 11136 x 5568

    #     # # lats (928, 3712) (1932430,): [26.656396865844727, ... 81.0744857788086]
    #     # # lats = lats[3712 - euro_nrows:, :]

    #     # # row (jj) starts at S and moves N
    #     # # col (ii) starts in E and moves W
    #     # # lons[jj, ii] where jj is row and ii is column
    #     # # 45W - 30E

    #     # # # jj = 2652 ii = 3146: (lat, lon) (23.81135368347168, -45.049129486083984)
    #     # # jj = 2652
    #     # # ii = 3146
    #     # # print(f"{jj = :4d} {ii = :4d}: (lat, lon) ({lats[jj, ii]}, {lons[jj, ii]})")

    #     # # # jj = 2652 ii =  910: (lat, lon) (23.114727020263672, 30.13174819946289)
    #     # # jj = 2652
    #     # # ii = 910
    #     # # print(f"{jj = :4d} {ii = :4d}: (lat, lon) ({lats[jj, ii]}, {lons[jj, ii]})")

    #     # # jj = 2652 ii =  910: (lat, lon) (23.114727020263672, 30.13174819946289)
    #     # jj = 2652
    #     # ii = 3146 - euro_ncols
    #     # print(f"{jj = :4d} {ii = :4d}: (lat, lon) ({lats[jj, ii]}, {lons[jj, ii]})")

    #     # 3146 - 1530 = 1616
    #     # # euro_nrows = 928
    #     # # euro_ncols = 1530
    #     # # 3146-910 = 2236
    #     # # 2236 - 1530 = 706
    #     # os._exit(1)

    #     # dones = 0
    #     # for jj in range(3712):
    #     #     for ii in range(3712):
    #     #         if np.abs(lons[jj, ii]) <= 180.0 and np.abs(lats[jj, ii] <= 90.0):
    #     #             if lats[jj, ii] >= 26.0 and lons[jj, ii] <= -45:
    #     #                 print(f"{jj = :4d} {ii = :4d}: (lat, lon) ({lats[jj, ii]}, {lons[jj, ii]})")
    #     #                 dones += 1
    #     #     # if dones > 5:
    #     #     #     os._exit(1)

    #     lats = lats[3712 - euro_nrows:, :]
    #     lons = lons[3712 - euro_nrows:, :]

    #     # 45W - 30E
    #     # 3712-1530
    #     if verbose:
    #         tmp = lons.flatten()
    #         tmp = tmp[np.abs(tmp) <= 180.0]
    #         print(f"\n\tlons {lons.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
    #         tmp = lats.flatten()
    #         tmp = tmp[np.abs(tmp) <= 90.0]
    #         print(f"\tlats {lats.shape} {tmp.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
    #         del tmp

    # # ny, nx = lons.shape
    # # # for yidx in range(53, 60):
    # # for yidx in range(ny):
    # #     left_good_idx = -1
    # #     right_good_idx = -1
    # #     ## print(f"Checking {yidx = :5d}")
    # #     for xidx in range(nx):
    # #         ## print(f"\tChecking {xidx = :5d} {data_vals[yidx, xidx]}")
    # #         if np.isnan(data_vals[yidx, xidx]):
    # #             # Invalid data value index
    # #             if left_good_idx >= 0:
    # #                 # Hit righthand end in the lat/row
    # #                 right_good_idx = xidx - 1
    # #                 print(f"* row {yidx:5d}: {left_good_idx:5d} {right_good_idx:5d}")
    # #                 row_lons = lons[yidx, left_good_idx:right_good_idx + 1]
    # #                 row_lats = lats[yidx, left_good_idx:right_good_idx + 1]
    # #                 print(f"\trow_lons ({len(row_lons)}): [{np.amin(row_lons):+9.3f} ... {np.amax(row_lons):+9.3f}]")
    # #                 print(f"\trow_lats ({len(row_lats)}): [{np.amin(row_lats):+9.3f} ... {np.amax(row_lats):+9.3f}]")
    # #                 break
    # #             continue
    # #         if left_good_idx < 0:
    # #             # First valid value in the lat/row on the lefthand side of the array
    # #             left_good_idx = xidx
    # #             # print(f"\t\tSet {left_good_idx = :5d}")
    # #         else:
    # #             # Post First value value in lat/row
    # #             if xidx == nx - 1:
    # #                 # Entire row valid
    # #                 right_good_idx = xidx
    # #                 # print(f"\t\tSet {right_good_idx = :5d}")
    # #                 print(f"+ row {yidx:5d}: {left_good_idx:5d} {right_good_idx:5d}")
    # #                 row_lons = lons[yidx, left_good_idx:right_good_idx + 1]
    # #                 row_lats = lats[yidx, left_good_idx:right_good_idx + 1]
    # #                 print(f"\trow_lons ({len(row_lons)}): [{np.amin(row_lons):+9.3f} ... {np.amax(row_lons):+9.3f}]")
    # #                 print(f"\trow_lats ({len(row_lats)}): [{np.amin(row_lats):+9.3f} ... {np.amax(row_lats):+9.3f}]")
    # #                 continue
    # #             right_good_idx = xidx
    # #             # print(f"\t\tSet ** {right_good_idx = :5d}")
    # #     if left_good_idx < 0:
    # #         # Entire row invalid
    # #         print(f"- row {yidx:5d}:")
    # # os._exit(1)

    # Find spatial extent of non-missing data
    tmp_vals = data_vals.flatten()
    tmp_lons = lons.flatten()
    tmp_lats = lats.flatten()
    good_lons = []
    good_lats = []
    good_vals = []
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
            good_vals.append(tmp_vals[ii])
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
        print(f"good_vals ({len(good_vals)}): [{np.amin(good_vals)} ... {np.amax(good_vals)}]")
        # print(f"ll_corner: {ll_corner}")
        # print(f"ur_corner: {ur_corner}")

    return

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
