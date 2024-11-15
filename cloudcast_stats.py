#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""

cloudcast_stats
~~~~~~~~~~~~~~~

$ python cloudcast_stats.py


curl --output  full_raw_cloud.zip https://vision.eng.au.dk/data/CloudDataset/full_raw_cloud.zip

"""
# Standard Imports
import os
import pickle
from calendar import monthrange
from pathlib import Path

# Third-Party Imports
import numpy as np
import numpy.ma as ma
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import BORDERS
import cartopy.feature as cfeature
from pyresample.geometry import AreaDefinition

# STARE Imports

# Local Imports
from find_grid_area import find_grid_area

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

# Define Global Constants and State Variables
# -------------------------------------------

use_raw = [False, True][1]
use_path = -1 if use_raw else 0

##
# Display region (not data values)
as_region = [False, True][0]

##
# Define path to folder
BASE_PATH = ["/Users/mbauer/tmp/", "/Volumes/saved/"][1]
FILE_PATH = [f"{BASE_PATH}data/CloudCast/full_cropped_cloud/", f"{BASE_PATH}data/CloudCast/small_cloud/", "/Volumes/saved/data/CloudCast/full_raw_cloud/"][use_path]

CCAST_YYYY = (2017, 2018)
CCAST_MM = tuple(range(1, 13))
if use_raw:
    # CCAST_NAMES = []
    # #   2017-01/CT/
    # #   S_NWC_CT_MSG3_MSG-N-VISIR_20170401T000000Z.nc
    # #   S_NWC_CT_MSG3_MSG-N-VISIR_20170401T001500Z.nc
    # for yyyy in CCAST_YYYY:
    #     for mm in CCAST_MM:
    #         dend = monthrange(yyyy, mm)[1]
    #         for dd in range(1, dend + 1):
    #             for hh in range(0, 1440, 15):
    #                 hour = hh  // 60
    #                 minute = hh % 60
    #                 CCAST_NAMES.append(f"{yyyy:04d}-{mm:02d}/CT/S_NWC_CT_MSG3_MSG-N-VISIR_{yyyy:04d}{mm:02d}{dd:02d}T{hour:02d}{minute:02d}00Z.nc")
    # CCAST_NAMES =  tuple(CCAST_NAMES)
    # # find '/Volumes/saved/data/CloudCast/full_raw_cloud' -type f | wc -l
    # # 70080
    # print(len(CCAST_NAMES))
    # # print(CCAST_NAMES[0], CCAST_NAMES[-1])

    # 70039
    CCAST_NAMES = list(Path('/Volumes/saved/data/CloudCast/full_raw_cloud').rglob("*.nc"))
    # CCAST_NAMES = sorted(list(set([os.path.split(_)[-1] for _ in CCAST_NAMES])))
    # print(len(CCAST_NAMES))
    # print(CCAST_NAMES[0], CCAST_NAMES[-1])
else:
    CCAST_NAMES = tuple((f"{yyyy:04d}M{mm:02d}.nc" for yyyy in CCAST_YYYY for mm in CCAST_MM))
# number of time steps in each monthly file, does not apply to use_raw
CCAST_TSTEPS_FILE = (2976, 2685, 2976, 2877, 2976, 2879, 2976, 2976, 2878, 2976, 2875, 2976,
                     2976, 2686, 2976, 2880, 2975, 2877, 2958, 2976, 2880, 2976, 2879, 2974)

if use_raw:
    # 70039
    CCAST_TSTEPS = len(CCAST_NAMES)
else:
    # 70039
    CCAST_TSTEPS = sum(CCAST_TSTEPS_FILE)

if use_raw:
    CCAST_CTYPES = {0: "No Cloud",
                    1: "Cloud-Free Land",
                    2: "Cloud-Free Sea",
                    3: "Snow Over Land",
                    4: "Sea Ice",
                    5: "Very Low Clouds",
                    6: "Low Clouds",
                    7: "Mid-Level Clouds",
                    8: "High Opaque Clouds",
                    9: "Very High Opaque Clouds",
                    10: "Fractional Clouds",
                    11: "High Semitransparent Thin Clouds",
                    12: "High Semitransparent Moderately Thick Clouds",
                    13: "High Semitransparent Thick Clouds",
                    14: "High Semitransparent Above Low or Medium Clouds",
                    15: "High Semitransparent Above Snow/Ice"}
    # Post correction
    CCAST_CTYPES = {0: "No Cloud",
                    1: "Very Low Cloud",
                    2: "Low Cloud",
                    3: "Mid-level Cloud",
                    4: "High Opaque Cloud",
                    5: "Very High Opaque Cloud",
                    6: "Fractional Cloud",
                    7: "High Semitransparent Thin Cloud",
                    8: "High Semitransparent Moderately Thick Cloud",
                    9: "High Semitransparent Thick Cloud",
                    10: "High Semitransparent above Low/Medium Cloud"}
else:
    CCAST_CTYPES = {0: "No Cloud",
                    1: "Very Low Cloud",
                    2: "Low Cloud",
                    3: "Mid-level Cloud",
                    4: "High Opaque Cloud",
                    5: "Very High Opaque Cloud",
                    6: "Fractional Cloud",
                    7: "High Semitransparent Thin Cloud",
                    8: "High Semitransparent Moderately Thick Cloud",
                    9: "High Semitransparent Thick Cloud",
                    10: "High Semitransparent above Low/Medium Cloud"}
N_CTYPES = len(CCAST_CTYPES.keys())

CCAST_CTYPES_Z = {"No and Fractional Cloud": (0, 6),
                  "Low Cloud": (1, 2),
                  "Mid Cloud": (3,),
                  "High Cloud": (4, 5, 7, 8, 9, 10)}
N_CTYPES_Z = 4

# WGS-84 Earth equatorial radius at sea level (meters)
GLOBE = ccrs.Globe(datum='WGS84', ellipse='WGS84')
# Geodetic:
#   A 3D/spherical CRS based on latitude and longitude where geographical distance and coordinates are measured in degrees.
GEOD_CRS = ccrs.Geodetic(globe=GLOBE)

##
# Define CloudCast spatial stereographic projection coordinates CRS
if use_raw:
    CCAST_HEIGHT = 928
    CCAST_WIDTH = 1530
    CCAST_1D_LEN = CCAST_HEIGHT * CCAST_WIDTH
    # lower_left_xy = [-855100.436345, -4942000.0]
    # upper_right_xy = [1448899.563655, -2638000.0]
    upper_left_xy = [-2296808.8, 5570249.0]
    lower_right_xy = [2293808.2, 2785874.8]
    lat_min_max = [26.671055, 81.09877]
    lon_min_max = [-69.27063, 69.27063]

    ##
    # Geostationary Projection (GEOS) EPSG
    #   https://proj4.org/en/9.2/operations/projections/geos.html
    #       +proj=geos +h=42164000.0 +R=6378000.0 +lon_0=0 +sweep=y
    #
    #   https://scitools.org.uk/cartopy/docs/latest/reference/projections.html#geostationary
    CCAST_CRS = ccrs.Geostationary(central_longitude=0.0, satellite_height=35785863, sweep_axis='y')
    #   gdal_projection = "+proj=geos +a=6378137.000000 +b=6356752.300000 +lon_0=0.000000 +h=35785863.000000 +sweep=y";

    area_id = f"GEOS"
    description = "Full Disk"
    proj_id = f"GEOS"
    area_def = AreaDefinition(area_id, description, proj_id,
                              {'lat_0': '0.00', 'lat_ts': '0.00', 'lon_0': '0.00', 'proj': 'geos', 'h': '35785863.0', 'sweep': 'y'},
                              CCAST_HEIGHT, CCAST_WIDTH,
                              (upper_left_xy[0], upper_left_xy[1], lower_right_xy[0], lower_right_xy[1]))

    ##
    # Form a cartopy CRS
    #   <class 'pyresample.utils.cartopy.Projection'>
    CCAST_CRS = area_def.to_cartopy_crs()
    # print(CCAST_CRS)

    #   (-2296808.8, 2293808.2, 5570249.0, 2785874.8)
    # print(CCAST_CRS.bounds)
else:
    CCAST_HEIGHT = 768
    CCAST_WIDTH = 768
    CCAST_1D_LEN = CCAST_HEIGHT * CCAST_WIDTH
    lower_left_xy = [-855100.436345, -4942000.0]
    upper_right_xy = [1448899.563655, -2638000.0]
    area_def = AreaDefinition('areaD', 'Europe', 'areaD',
                              {'lat_0': '90.00', 'lat_ts': '50.00',
                               'lon_0': '5', 'proj': 'stere', 'ellps': 'WGS84'},
                              CCAST_HEIGHT, CCAST_WIDTH,
                              (lower_left_xy[0], lower_left_xy[1],
                               upper_right_xy[0], upper_right_xy[1]))
    CCAST_CRS = area_def.to_cartopy_crs()

SAVE_DIR = f"{BASE_PATH}hidden/cloudcast/"
if use_raw:
    SAVE_FILE = "cloudcast_stats_raw.pkl"
else:
    SAVE_FILE = "cloudcast_stats.pkl"

###############################################################################
# PUBLIC read_ccast()
# -------------------
def read_ccast(fname: str, get_coords: bool, as_raw: bool):
    verbose = [False, True][0]
    x_coords = []
    y_coords = []
    tstamps = []

    # Load dataset using xarray
    # -------------------------------------------
    if verbose:
        print(f"\tread_ccast({fname})")

    r"""
    read_data =
        Variable
            xarray.DataArray(lat: 768, lon: 768, time: 2879) > Size: 2GB
        Coordinates:
          * lat      (768) float64          6kB -8.536e+05 -8.506e+05 ... 1.444e+06 1.447e+06
          * lon      (768) float64          6kB -2.64e+06 -2.642e+06  ... -4.938e+06 -4.94e+06
          * time     (2976) datetime64[ns] 23kB 2017-06-01T00:09:17   ... 2017-06-30T23...

    read_dataset =
        Dimensions:            (ny: 928, nx: 1530, pal01_colors: 256, pal_RGB: 3,
                                pal02_colors: 256, pal03_colors: 256)
        Coordinates:
            lat                (ny, nx) float32 6MB ...
            lon                (ny, nx) float32 6MB ...
          * ny                 (ny) float32 4kB 5.569e+06 5.566e+06 ... 2.787e+06
          * nx                 (nx) float32 6kB -2.295e+06 -2.292e+06 ... 2.292e+06
        Dimensions without coordinates: pal01_colors, pal_RGB, pal02_colors,
                                        pal03_colors
        Data variables:
            ct                 (ny, nx) float32 6MB ...
            ct_cumuliform      (ny, nx) float32 6MB ...
            ct_multilayer      (ny, nx) float32 6MB ...
            ct_status_flag     (ny, nx) float32 6MB ...
            ct_conditions      (ny, nx) float32 6MB ...
            ct_quality         (ny, nx) float32 6MB ...
            ct_pal             (pal01_colors, pal_RGB) uint8 768B ...
            ct_cumuliform_pal  (pal02_colors, pal_RGB) uint8 768B ...
            ct_multilayer_pal  (pal03_colors, pal_RGB) uint8 768B ...
        Attributes: (12/46)
            Conventions:                  CF-1.6
            title:                        NWC GEO Cloud Type Product
            history:                      2019-06-11T07:46:14Z (null) Product Created...
            institution:                  Aarhus University
            source:                       NWC/GEO version v2018
            comment:                      Copyright 2019, EUMETSAT, All Rights Reserved
            ...                           ...
            product_quality:              66.86354
            product_completeness:         99.46954
            geospatial_lat_max:           81.09877
            geospatial_lat_min:           26.671055
            geospatial_lon_max:           69.27063
            geospatial_lon_min:           -69.27063
    """
    if as_raw:
        read_dataset = xr.open_dataset(fname)
        # (ny: 928, nx: 1530)>
        read_data = read_dataset.variables['ct'][:]
    else:
        read_data = xr.open_dataarray(fname)

    # Pre-processing to match cloud types in paper
    # -------------------------------------------

    ##
    # Remove classes 1, 2, 3 and 4, which are cloud-free land, cloud-free sea, snow over land and sea ice.
    """
    <xarray.DataArray (lat: 768, lon: 768, time: 2976)> Size: 7GB
    array([[[nan, nan, 10., ...,  7.,  7.,  8.],
            [nan, nan,  5., ...,  7.,  7.,  8.],
            [10., 10., 10., ...,  7.,  7., 12.],
            ...,
            [10., 10.,  5., ..., nan, nan, nan],
            [10., 10.,  5., ..., nan, nan, nan],
            [10., 10., 10., ..., nan,  5., nan]],

           [[10., nan, 10., ...,  7., 12.,  8.],
            [nan, nan,  5., ...,  7.,  7.,  8.],
            [10., 10., 10., ...,  7.,  7., 12.],
            ...,
            [10., 10.,  5., ..., nan, nan, nan],
            [10., 10., 10., ..., nan, nan, nan],
            [10., 10., 10., ..., nan, nan, nan]],

           [[10., nan, 10., ...,  7., 12.,  8.],
            [nan, 10., 10., ...,  7., 12.,  8.],
            [ 5., 10., 10., ..., 10.,  7.,  8.],
            ...,
    ...
            ...,
            [nan, nan, nan, ..., nan, nan, nan],
            [nan, nan, nan, ..., nan, nan, nan],
            [ 6., 10., 10., ..., 10., 10., nan]],

           [[ 7.,  7.,  7., ..., 12.,  8., 13.],
            [ 7.,  7.,  7., ..., 13., 13., 13.],
            [ 7.,  7.,  7., ..., 13., 12., 13.],
            ...,
            [ 6., nan, 10., ..., nan, nan, nan],
            [ 6., 10., 10., ..., nan, nan, nan],
            [ 6., 10., 10., ..., nan, nan, nan]],

           [[ 7.,  7.,  7., ..., 12., 13., 13.],
            [ 7.,  7.,  7., ..., 13., 13., 13.],
            [ 7.,  7.,  7., ..., 13., 12., 13.],
            ...,
            [ 6., nan, 10., ..., nan, nan, nan],
            [ 6., nan, 10., ..., nan, nan, nan],
            [ 6., 10., 10., ..., nan, nan, nan]]], dtype=float32)
    """
    read_data = read_data.where(read_data > 4)

    ##
    # Subtract 4 to correspond to paper cloud types
    read_data = read_data - 4

    ##
    # Set nans to zero
    read_data = read_data.fillna(0)

    ##
    # Change to int to save memory (7Gb to 2Gb)
    read_data = read_data.astype(np.uint8)
    # print(read_data)
    # os._exit(1)

    ##
    # Extract numpy array (768, 768, 2976) of np.uint8
    ccast_dat = read_data.values

    if as_raw:
        tstamps = []
    else:
        tstamps = read_data.coords['time'].values

    if get_coords:
        ##
        # Pull spatial stereographic projection coordinates
        if as_raw:
            #   x_coords (928, 1530): [-69.3 ... 69.3]
            #   y_coords (928, 1530): [26.67105484008789, ... 81.09877014160156]
            x_coords = read_dataset.variables['lon'].values
            y_coords = read_dataset.coords['lat'].values
            if verbose:
                tmp = x_coords.flatten()
                tmp = tmp[np.abs(tmp) <= 180.0]
                print(f"x_coords {x_coords.shape}: [{np.amin(tmp):8.1f} ... {np.amax(tmp):8.1f}]")
                tmp = y_coords.flatten()
                tmp = tmp[np.abs(tmp) <= 90.0]
                print(f"y_coords {y_coords.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
        else:
            #   Note lon and lat seem flipped and y_coords needs to be reversed.
            #       y_coords (768): [-4940500.0 ... -2639500.0]
            #       x_coords (768): [-853600.4  ...  1447399.6]
            # It seems the y-axis needs to be reversed to match this corners
            #   lower_left_xy = [-855100.436345, -4942000.0]
            #   upper_right_xy = [1448899.563655, -2638000.0]
            # Giving
            y_coords = read_data.coords['lon'].values
            y_coords = y_coords[::-1]
            x_coords = read_data.coords['lat'].values
            if verbose:
                print(f"y_coords ({len(y_coords)}): [{y_coords[0]:8.1f} ... {y_coords[-1]:8.1f}]")
                print(f"x_coords ({len(x_coords)}): [{x_coords[0]:8.1f} ... {x_coords[-1]:8.1f}]")

        # ##
        # # Save x_coords, y_coords for lon/lat domain matching with MSG
        # with open(f"{SAVE_DIR}raw_coords.pkl", 'wb') as f:
        #     tmp = (x_coords, y_coords)
        #     pickle.dump(tmp, f)
        # os._exit(1)

    # ##
    # # To match flip and reverse of x_coords and y_coords
    # ccast_dat = np.swapaxes(ccast_dat, 0, 1)
    # ccast_dat = ccast_dat[:, ::-1, :]

    ##
    # Ensure memory released
    del read_data
    if as_raw:
        del read_dataset

    return ccast_dat, x_coords, y_coords, tstamps


###############################################################################
# PUBLIC main()
# -------------
def main():

    verbose = [False, True][1]
    just_one = [False, True][0]
    recall_stats = [False, True][0]

    ##
    # Define stat storage
    ctype_cnts = np.zeros((N_CTYPES, CCAST_TSTEPS), dtype=int)
    ctype_map_cnts = np.zeros((N_CTYPES, CCAST_1D_LEN), dtype=int)

    sfile = f"{SAVE_DIR}{SAVE_FILE}"

    if recall_stats:
        with open(sfile, 'rb') as f:
            tmp = pickle.load(f)
            ctype_cnts, ctype_map_cnts = tmp
            del tmp
    else:
        ##
        # Loop over all files/time
        for ridx, rfile in enumerate(CCAST_NAMES):
            get_file = str(rfile) if use_raw else f"{FILE_PATH}{rfile}"
            if verbose:
                print(f"Reading {ridx:2d} {get_file.split('/')[-1]}")

            # # TMP
            # if get_file.split('/')[-1].startswith("2017"):
            #     continue
            # if ridx <= 0:
            #     continue

            if just_one:
                # Need to save see TMP below
                with open(f"{SAVE_DIR}tmp.pkl", 'rb') as f:
                    tmp = pickle.load(f)
                ccast_dat, x_coords, y_coords, tstamps = tmp
            else:
                ##
                # Read data file
                ccast_dat, x_coords, y_coords, tstamps = read_ccast(get_file, 1 if ridx == 0 else 0, as_raw=use_raw)
                ## print(f"{get_file.split('/')[-1]} {len(tstamps)}")

                # # TMP save/recall for debugging w/ just_one
                # with open(f"{SAVE_DIR}tmp.pkl", 'wb') as f:
                #     tmp = (ccast_dat, x_coords, y_coords, tstamps)
                #     pickle.dump(tmp, f)
                # os._exit(1)

            if ridx == 0:
                ##
                # Convert from projection (Cartesian Coordinates) to geographic (Spherical Coordinates)
                if use_raw:
                    #   ccast_lons (928, 1530): [-69.2706298828125 ... 69.2706298828125]
                    #   ccast_lats (928, 1530): [26.67105484008789 ... 81.09877014160156]
                    ccast_lons = x_coords
                    ccast_lats = y_coords
                    if verbose:
                        tmp = ccast_lons.flatten()
                        tmp = tmp[np.abs(tmp) <= 180.0]
                        print(f"ccast_lons {ccast_lons.shape}: [{np.amin(tmp)} ... {np.amax(tmp)}]")
                        tmp = ccast_lats.flatten()
                        tmp = tmp[np.abs(tmp) <= 90.0]
                        print(f"ccast_lats {ccast_lats.shape}: [{np.amin(tmp)} ... {np.amax(tmp)}]")
                else:
                    #   ccast_lons (768): [-12.921 ... 21.329] W to E
                    #   ccast_lats (768): [ 62.403 ... 40.928] N to S
                    transformed = GEOD_CRS.transform_points(CCAST_CRS, x_coords, y_coords)
                    ccast_lons = transformed[..., 0].tolist()
                    ccast_lats = transformed[..., 1].tolist()
                    ccast_nlons = len(ccast_lons)
                    ccast_nlats = len(ccast_lats)
                    if verbose:
                        print(f"ccast_lons ({ccast_nlons}): [{ccast_lons[0]} ... {ccast_lons[-1]}]")
                        print(f"ccast_lats ({ccast_nlats}): [{ccast_lats[0]} ... {ccast_lats[-1]}]")

                if not use_raw:
                    ##
                    # Surface area of the data-grid (Latitude-Longitude Quadrangle) [km^2]
                    #   Total Surface Area: 6,548,627 km^2
                    ccast_grid_area = find_grid_area(ccast_lons, ccast_lats)
                    ccast_total_area = sum(ccast_grid_area[_] * ccast_nlons for _ in range(ccast_nlats))
                    if verbose:
                        print(f"\tTotal Surface Area: {int(ccast_total_area):,} km^2")

            if use_raw:
                ##
                # Collect some 1D stats
                ccast_dat_flat = ccast_dat.flatten()

                ##
                # Populate counting map
                for ctype in range(N_CTYPES):
                    ##
                    # Indices of this ctype
                    cidx = np.where(ccast_dat_flat == ctype)[0]
                    # print(cidx)
                    ##
                    # Increment map
                    ctype_map_cnts[ctype, cidx] += 1
                    # print(f"\t{ctype} {len(cidx)} {np.sum(ctype_map_cnts[ctype, :])}")

                unique, counts = np.unique(ccast_dat_flat, return_counts=True)
                # total_cnt = sum(counts)
                for cidx, ctype in enumerate(unique):
                    # print(f"{cidx} {ctype} {counts[cidx]}")
                    if ctype > 10:
                        print(f"\tFound ctype {ctype} {counts[cidx]}")
                        continue
                    ctype_cnts[ctype, ridx] = counts[cidx]
            else:
                ##
                # Collect some 1D stats
                #   ccast_dat [CCAST_HEIGHT, CCAST_WIDTH, CCAST_TSTEPS_FILE[ridx]]
                for tidx in range(CCAST_TSTEPS_FILE[ridx]):
                    ccast_dat_flat = ccast_dat[:, :, tidx].flatten()
                    # print(f"{tidx:04d} Time Inded {tstamps[tidx]}")
                    ##
                    # Populate counting map
                    for ctype in range(N_CTYPES):
                        ##
                        # Indices of this ctype
                        cidx = np.where(ccast_dat_flat == ctype)[0]

                        ##
                        # Increment map
                        ctype_map_cnts[ctype, cidx] += 1
                        # print(f"\t{ctype} {len(cidx)} {np.sum(ctype_map_cnts[ctype, :])}")

                    unique, counts = np.unique(ccast_dat_flat, return_counts=True)
                    # total_cnt = sum(counts)
                    for cidx, ctype in enumerate(unique):
                        # print(f"{cidx} {ctype} {counts[cidx]}")
                        if ctype > 10:
                            print(f"\tFound ctype {ctype} {counts[cidx]}")
                            continue
                        ctype_cnts[ctype, tidx] = counts[cidx]
                    #  os._exit(1)
            if just_one:
                break
        ##
        # Save for later
        with open(sfile, 'wb') as f:
            tmp = (ctype_cnts, ctype_map_cnts)
            pickle.dump(tmp, f)

    ##
    # Stat IO
    """
    Reading  0 2017M01.nc
        Total Surface Area: 6,548,627 km^2
    Reading  1 2017M02.nc
    Reading  2 2017M03.nc
        Found ctype 251 199
    Reading  3 2017M04.nc
    Reading  4 2017M05.nc
    Reading  5 2017M06.nc
    Reading  6 2017M07.nc
    Reading  7 2017M08.nc
        Found ctype 251 4
        Found ctype 251 183
    Reading  8 2017M09.nc
        Found ctype 251 346
        Found ctype 251 1595
        Found ctype 251 396
    Reading  9 2017M10.nc
    Reading 10 2017M11.nc
    Reading 11 2017M12.nc
    Reading 12 2018M01.nc
    Reading 13 2018M02.nc
        Found ctype 251 4194
    Reading 14 2018M03.nc
        Found ctype 251 8781
    Reading 15 2018M04.nc
        Found ctype 251 573
    Reading 16 2018M05.nc
        Found ctype 251 8709
        Found ctype 251 317
        Found ctype 251 444
        Found ctype 251 1736
        Found ctype 251 285684
    Reading 17 2018M06.nc
    Reading 18 2018M07.nc
    Reading 19 2018M08.nc
    Reading 20 2018M09.nc
        Found ctype 251 427
        Found ctype 251 1082
        Found ctype 251 3324
        Found ctype 251 246
    Reading 21 2018M10.nc
        Found ctype 251 83771
    Reading 22 2018M11.nc
    Reading 23 2018M12.nc

    Total Pixel Counts by CType
     0:       313508079   17.70448%     No clouds or missing data
     1:       263946051   14.90560%     Very low clouds
     2:       159297076    8.99585%     Low clouds
     3:       197203809   11.13652%     Mid-level clouds
     4:       184449485   10.41626%     High opaque clouds
     5:        21787279    1.23037%     Very high opaque clouds
     6:       142539617    8.04952%     Fractional clouds
     7:        65357760    3.69089%     High semitransparent thin clouds
     8:       174101311    9.83188%     High semitransparent moderately thick clouds
     9:       208207821   11.75794%     High semitransparent thick clouds
    10:        40386036    2.28069%     High semitransparent above low or medium clouds
    """
    # print("\nTotal Pixel Counts by CType")
    # tsum = np.sum(ctype_cnts)
    # for cidx, ctype in enumerate(CCAST_CTYPES.keys()):
    #     csum = np.sum(ctype_cnts[cidx, :])
    #     print(f"{ctype:2d}: {csum:15d} {100.0 * (csum / tsum):10.5f}%\t{CCAST_CTYPES[ctype]}")
    # # os._exit(1)

    ##
    # Work with counting map
    max_pixels = CCAST_1D_LEN * CCAST_TSTEPS
    ctype_map_occurence = np.zeros((N_CTYPES, CCAST_1D_LEN), dtype=float)
    ctype_z_map_occurence = np.zeros((N_CTYPES_Z, CCAST_1D_LEN), dtype=float)
    ##
    # Set to False if you do not want background image
    use_nasa_background = [False, True][0]
    # freq_map_cmap = "plasma"
    freq_map_cmap = "brg"

    """
    Pixel Time Occurrence/Frequency by CType
     0: Min    2.51146% Mean   27.04586% Max   61.10453%    No cloud
     1: Min    0.30412% Mean   13.12638% Max   29.21515%    Very low cloud
     2: Min    3.26247% Mean   12.25170% Max   23.14282%    Low cloud
     3: Min    3.50805% Mean   11.39559% Max   31.93364%    Mid-level cloud
     4: Min    4.61172% Mean    9.71385% Max   18.65818%    High opaque cloud
     5: Min    0.28698% Mean    0.79755% Max    1.96605%    Very high opaque cloud
     6: Min    0.02427% Mean    7.93259% Max   25.19453%    Fractional cloud
     7: Min    0.20560% Mean    3.48593% Max    6.45783%    High semitransparent thin cloud
     8: Min    4.44895% Mean    6.44301% Max   10.85966%    High semitransparent moderately thick cloud
     9: Min    2.28444% Mean    5.30360% Max   13.38683%    High semitransparent thick cloud
    10: Min    0.89664% Mean    2.50297% Max    6.40358%    High semitransparent above low/medium cloud
    """
    # print("\nPixel Time Occurrence/Frequency by CType")
    # for ctype in range(N_CTYPES):
    #     ##
    #     # Fractional time coverage
    #     if just_one:
    #         ctype_map_occurence[ctype, :] = np.divide(ctype_map_cnts[ctype, :], 2976) * 100.0
    #     else:
    #         ctype_map_occurence[ctype, :] = np.divide(ctype_map_cnts[ctype, :], CCAST_TSTEPS) * 100.0
    #     # ctype_map_occurence[ctype, :] = ctype_map_cnts[ctype, :]
    #     print(f"{ctype:2d}: Min {np.amin(ctype_map_occurence[ctype, :]):10.5f}% Mean {np.mean(ctype_map_occurence[ctype, :]):10.5f}% Max {np.amax(ctype_map_occurence[ctype, :]):10.5f}%\t{CCAST_CTYPES[ctype]}")

    #     ##
    #     # For plotting purposes, we remove 0 (so we can actually see land)
    #     mapdat = np.copy(ctype_map_occurence[ctype, :])
    #     # mapdat = ma.masked_where(mapdat > 0, mapdat, copy=True)
    #     mapdat = np.reshape(mapdat, (CCAST_HEIGHT, CCAST_WIDTH))

    #     pname = f"{SAVE_DIR}ccast_freq_map_{ctype:02d}.png"
    #     fig = plt.figure(figsize=(10, 8))
    #     ax = plt.axes(projection=CCAST_CRS)
    #     if use_nasa_background:
    #         ax.background_img(name='BM', resolution='low')
    #     else:
    #         ax.add_feature(cfeature.COASTLINE, alpha=0.5)
    #     a_image = plt.imshow(mapdat, cmap=freq_map_cmap, transform=CCAST_CRS, extent=CCAST_CRS.bounds, origin='upper')
    #     fig.colorbar(a_image, ax=ax, label=f"Time Frequency of Occurrence for CloudType {ctype}\n{CCAST_CTYPES[ctype]} [%]")
    #     fig.savefig(pname, facecolor='w', edgecolor='w', orientation='landscape', dpi=300)
    #     plt.clf()
    #     plt.close('all')

    """
    Pixel Time Occurrence/Frequency by Layer CType
     0: Min    9.34051% Mean   34.97845% Max   67.87932%    No and Fractional Cloud
     1: Min    5.53692% Mean   25.37808% Max   45.81590%    Low Cloud
     2: Min    3.50805% Mean   11.39559% Max   31.93364%    Mid Cloud
     3: Min   20.09595% Mean   28.24690% Max   41.37695%    High Cloud

    Raw
    Pixel Time Occurrence/Frequency by Layer CType
     0: Min    1.92036% Mean   51.75511% Max  100.00000%    No and Fractional Cloud
     1: Min    0.00000% Mean   21.81110% Max   55.05790%    Low Cloud
     2: Min    0.00000% Mean    6.97824% Max   35.27035%    Mid Cloud
     3: Min    0.00000% Mean   19.45555% Max   79.03168%    High Cloud
    """
    if as_region:
        pass
    else:
        if use_raw:
            print("\nRAW Pixel Time Occurrence/Frequency by Layer CType")
        else:
            print("\nPixel Time Occurrence/Frequency by Layer CType")
        for cidx, ctype in enumerate(CCAST_CTYPES_Z.keys()):
            ##
            # Fractional time coverage
            if just_one:
                ctype_z_map_occurence[cidx, :] = np.divide(ctype_map_cnts[ctype, :], 2976) * 100.0
            else:
                tmp = np.zeros((CCAST_1D_LEN,), dtype=int)
                for didx in CCAST_CTYPES_Z[ctype]:
                    tmp += ctype_map_cnts[didx, :]
                ctype_z_map_occurence[cidx, :] = np.divide(tmp, CCAST_TSTEPS) * 100.0
            print(f"{cidx:2d}: Min {np.amin(ctype_z_map_occurence[cidx, :]):10.5f}% Mean {np.mean(ctype_z_map_occurence[cidx, :]):10.5f}% Max {np.amax(ctype_z_map_occurence[cidx, :]):10.5f}%\t{ctype}")

            ##
            # For plotting purposes, we remove 0 (so we can actually see land)
            mapdat = np.copy(ctype_z_map_occurence[cidx, :])
            mapdat = np.reshape(mapdat, (CCAST_HEIGHT, CCAST_WIDTH))

            if use_raw:
                pname = f"{SAVE_DIR}ccast_raw_freq_map_z_{cidx:02d}.png"
            else:
                pname = f"{SAVE_DIR}ccast_freq_map_z_{cidx:02d}.png"
            fig = plt.figure(figsize=(10, 8))
            ax = plt.axes(projection=CCAST_CRS)
            if use_nasa_background:
                ax.background_img(name='BM', resolution='low')
            else:
                ax.add_feature(cfeature.COASTLINE, alpha=0.5)

            # norm=None
            if use_raw:
                a_image = plt.imshow(mapdat, cmap=freq_map_cmap, transform=CCAST_CRS, extent=CCAST_CRS.bounds, vmin=0, vmax=100, origin='lower')
            else:
                a_image = plt.imshow(mapdat, cmap=freq_map_cmap, transform=CCAST_CRS, extent=CCAST_CRS.bounds, vmin=0, vmax=100, origin='upper')
            fig.colorbar(a_image, ax=ax, label=f"Time Frequency of Occurrence\n{ctype} [%]")
            fig.savefig(pname, facecolor='w', edgecolor='w', orientation='landscape', dpi=300)
            plt.clf()
            plt.close('all')

    ##
    # Make RGB Composite
    print("\n\nRGB Composite")
    def scale_min_max(in_ndarray):
        return (in_ndarray - np.nanmin(in_ndarray)) / (np.nanmax(in_ndarray) - np.nanmin(in_ndarray))

    if as_region:
        rgb_normalized = np.ones((CCAST_HEIGHT, CCAST_WIDTH))
    else:
        for cidx, ctype in enumerate(CCAST_CTYPES_Z.keys()):
            if cidx == 0:
                continue
            elif cidx == 1:
                # Low Cloud
                #   red_band (768, 768) = array([[34.93339425, ..., 28.85820757], ... 19.94174674]])
                red_band = np.copy(ctype_z_map_occurence[cidx, :])
                red_band = np.reshape(red_band, (CCAST_HEIGHT, CCAST_WIDTH))
            elif cidx == 2:
                # Mid Cloud
                green_band = np.copy(ctype_z_map_occurence[cidx, :])
                green_band = np.reshape(green_band, (CCAST_HEIGHT, CCAST_WIDTH))
            elif cidx == 3:
                # High Cloud
                blue_band = np.copy(ctype_z_map_occurence[cidx, :])
                blue_band = np.reshape(blue_band, (CCAST_HEIGHT, CCAST_WIDTH))

        ##
        # Define color map
        # color_map = {1: np.array([255, 0, 0]), # red
        #              2: np.array([0, 255, 0]), # green
        #              3: np.array([0, 0, 255])} # blue

        #   red_normalized from 0 to 1
        red_normalized = scale_min_max(red_band)
        ## red_normalized = 0.0 * red_normalized
        #   red_normalized from 0 to 255
        ## red_normalized = np.round(red_normalized * 255.0, decimals=0).astype(int)

        green_normalized = scale_min_max(green_band)
        ## green_normalized = 0.0 * green_normalized
        ## green_normalized = np.round(green_normalized * 255.0, decimals=0).astype(int)

        blue_normalized = scale_min_max(blue_band)
        ## blue_normalized = 0.0 * blue_normalized
        ## blue_normalized = np.round(blue_normalized * 255.0, decimals=0).astype(int)

        #   rgb_normalized = [[[186 131 147] ... [187 131 147]]]
        # rgb_normalized = np.dstack((red_normalized, green_normalized, blue_normalized)).astype(int)
        rgb_normalized = np.dstack((red_normalized, green_normalized, blue_normalized))

        r"""
        RGB Composite
            Red  : Min     0 Mean   126 Max   255
            Green: Min     0 Mean    71 Max   255
            Blue : Min     0 Mean    98 Max   255
            RGB  : Min     0 Mean    98 Max   255
        """
        print(f"\tRed  : Min {np.amin(red_normalized):6.3f} Mean {np.mean(red_normalized):6.3f} Max {np.amax(red_normalized):6.3f}")
        print(f"\tGreen: Min {np.amin(green_normalized):6.3f} Mean {np.mean(green_normalized):6.3f} Max {np.amax(green_normalized):6.3f}")
        print(f"\tBlue : Min {np.amin(blue_normalized):6.3f} Mean {np.mean(blue_normalized):6.3f} Max {np.amax(blue_normalized):6.3f}")
        print(f"\tRGB  : Min {np.amin(rgb_normalized):6.3f} Mean {np.mean(rgb_normalized):6.3f} Max {np.amax(rgb_normalized):6.3f}")

    if use_raw:
        pname = f"{SAVE_DIR}ccast_raw_freq_map_z_RGB.png"
    else:
        pname = f"{SAVE_DIR}ccast_freq_map_z_RGB.png"
    if as_region:
        pname = pname.replace(".png", "_region.png")

    fig = plt.figure(figsize=(10, 8))
    ax = plt.axes(projection=CCAST_CRS)
    ax.add_feature(cfeature.COASTLINE, alpha=0.5)

    # from matplotlib.colors import LinearSegmentedColormap
    # cast_rgb = LinearSegmentedColormap('CAST_RGB', cdict1)

    # norm=None
    if use_raw:
        if as_region:
            a_image = plt.imshow(rgb_normalized, transform=CCAST_CRS, extent=CCAST_CRS.bounds, origin='lower', cmap="bone_r", vmin=0.0, vmax=1.0, alpha=0.25)
            ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=2, color='black', alpha=0.5, linestyle='--')
        else:
            a_image = plt.imshow(rgb_normalized, transform=CCAST_CRS, extent=CCAST_CRS.bounds, origin='lower')
    else:
        # a_image = plt.imshow(rgb_normalized, cmap=freq_map_cmap, transform=CCAST_CRS, extent=CCAST_CRS.bounds, vmin=0, vmax=1, origin='upper')
        a_image = plt.imshow(rgb_normalized, transform=CCAST_CRS, extent=CCAST_CRS.bounds, origin='upper')
        # a_image = plt.imshow(rgb_normalized, cmap=freq_map_cmap, transform=CCAST_CRS, extent=CCAST_CRS.bounds, origin='upper')

    # fig.colorbar(a_image, ax=ax, label=f"Time Frequency of Occurrence\n{ctype} [%]")

    fig.savefig(pname, facecolor='w', edgecolor='w', orientation='landscape', dpi=300)
    plt.clf()
    plt.close('all')

# ---Start of main code block.
if __name__ == '__main__':
    ##
    # Run the main routine
    main()

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<



