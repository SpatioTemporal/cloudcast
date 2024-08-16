#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""

cloudcast_stats
~~~~~~~~~~~~~~~

$ python cloudcast_stats.py
"""
# Standard Imports
import os
import pickle

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

##
# Define path to folder
FILE_PATH = ["/Volumes/val/data/CloudCast/full_cropped_cloud/", "/Volumes/val/data/CloudCast/small_cloud/"][0]

CCAST_YYYY = (2017, 2018)
CCAST_MM = tuple(range(1, 13))
CCAST_NAMES = tuple((f"{yyyy:04d}M{mm:02d}.nc" for yyyy in CCAST_YYYY for mm in CCAST_MM))
# number of time steps in each monthly file
CCAST_TSTEPS_FILE = (2976, 2685, 2976, 2877, 2976, 2879, 2976, 2976, 2878, 2976, 2875, 2976,
                     2976, 2686, 2976, 2880, 2975, 2877, 2958, 2976, 2880, 2976, 2879, 2974)

# 70039
CCAST_TSTEPS = sum(CCAST_TSTEPS_FILE)

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

SAVE_DIR = "/Volumes/val/hidden/cloudcast/"
SAVE_FILE = "cloudcast_stats.pkl"


###############################################################################
# PUBLIC read_ccast()
# -------------------
def read_ccast(fname: str, get_coords: bool):
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
          * lat      (lat) float64          6kB -8.536e+05 -8.506e+05 ... 1.444e+06 1.447e+06
          * lon      (lon) float64          6kB -2.64e+06 -2.642e+06  ... -4.938e+06 -4.94e+06
          * time     (time) datetime64[ns] 23kB 2017-06-01T00:09:17   ... 2017-06-30T23...

    """
    read_data = xr.open_dataarray(fname)
    # print(read_data)
    # os._exit(1)

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
    tstamps = read_data.coords['time'].values

    if get_coords:
        ##
        # Pull spatial stereographic projection coordinates
        x_coords = read_data.coords['lon'].values
        y_coords = read_data.coords['lat'].values

    ##
    # Ensure memory released
    del read_data

    return ccast_dat, x_coords, y_coords, tstamps


###############################################################################
# PUBLIC main()
# -------------
def main():

    verbose = [False, True][1]
    just_one = [False, True][0]
    recall_stats = [False, True][1]

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
            get_file = f"{FILE_PATH}{rfile}"
            if verbose:
                print(f"Reading {ridx:2d} {get_file.split('/')[-1]}")

            # # TMP
            # if get_file.split('/')[-1].startswith("2017"):
            #     continue
            # if ridx <= 0:
            #     continue

            if just_one:
                with open('tmp.pkl', 'rb') as f:
                    tmp = pickle.load(f)
                ccast_dat, x_coords, y_coords, tstamps = tmp
            else:
                ##
                # Read data file
                ccast_dat, x_coords, y_coords, tstamps = read_ccast(get_file, 1 if ridx == 0 else 0)
                ## print(f"{get_file.split('/')[-1]} {len(tstamps)}")

                # # TMP save/recall for debugging
                # # with open('tmp.pkl', 'wb') as f:
                # #     tmp = (ccast_dat, x_coords, y_coords, tstamps)
                # #     pickle.dump(tmp, f)

            if ridx == 0:
                ccast_x = x_coords
                ccast_y = y_coords
                ##
                # Convert from projection (Cartesian Coordinates) to geographic (Spherical Coordinates)
                transformed = GEOD_CRS.transform_points(CCAST_CRS, ccast_x, ccast_y)
                ccast_lons = transformed[..., 0].tolist()
                ccast_lats = transformed[..., 1].tolist()
                ccast_nlons = len(ccast_lons)
                ccast_nlats = len(ccast_lats)

                ##
                # Surface area of the data-grid (Latitude-Longitude Quadrangle) [km^2]
                #   Total Surface Area: 6,548,627 km^2
                ccast_grid_area = find_grid_area(ccast_lons, ccast_lats)
                ccast_total_area = sum(ccast_grid_area[_] * ccast_nlons for _ in range(ccast_nlats))
                if verbose:
                    print(f"\tTotal Surface Area: {int(ccast_total_area):,} km^2")

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
    freq_map_cmap = "plasma"

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
    """
    print("\nPixel Time Occurrence/Frequency by Layer CType")
    for cidx, ctype in enumerate(CCAST_CTYPES_Z.keys()):
        ##
        # Fractional time coverage
        if just_one:
            os.exit(1)
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

        pname = f"{SAVE_DIR}ccast_freq_map_z_{cidx:02d}.png"
        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes(projection=CCAST_CRS)
        if use_nasa_background:
            ax.background_img(name='BM', resolution='low')
        else:
            ax.add_feature(cfeature.COASTLINE, alpha=0.5)
        a_image = plt.imshow(mapdat, cmap=freq_map_cmap, transform=CCAST_CRS, extent=CCAST_CRS.bounds, origin='upper')
        fig.colorbar(a_image, ax=ax, label=f"Time Frequency of Occurrence\n{ctype} [%]")
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


