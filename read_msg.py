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
import sys
import pickle
from pathlib import Path

# Third-Party Imports
import numpy as np
import numpy.ma as ma
from osgeo import gdal
import pyresample as pr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from satpy import Scene
from zipfile import ZipFile
import netCDF4
from scipy.interpolate import griddata
from tqdm import tqdm
from multiprocessing import Pool

# STARE Imports

# Local Imports
from cloudcast.cfg.setup_msg import setup_msg
from cloudcast.util.natread import natread
from cloudcast.util.nat2tif import nat2tif
from cloudcast.util.calculate_solar_angles import calculate_solar_angles
from cloudcast.cfg.setup_era5 import setup_era5
from cloudcast.util.read_era5 import read_era5
from cloudcast.util.spatial_downscale import spatial_downscale
from cloudcast.plot.plot_era5_ccast import plot_era5_ccast

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

# Define Global Constants and State Variables
# -------------------------------------------

##
# Multiprocessing cores
N_CORES = 4

###############################################################################
# PUBLIC main()
# -------------
def main(make_tif: bool, read_nat: bool, make_fig: bool, use_nat: bool, zip2nat: bool, as_full: bool, as_euro: bool, as_ccast: bool, as_merc: bool, as_lcc: bool, as_region: bool, verbose: bool, mod_t_files: str, use_tag: str, use_dataset: str, freq_map_cmap: str, SUB_PATH: str, FNAME: str, TNAME: str, ZNAME: str, t_src_dir: str, lonlatfile: str, t_levs: list[str], t_years: list[int]) -> None:

    ##
    # Define reader (GDAL)
    #   https://satpy.readthedocs.io/en/stable/api/satpy.readers.seviri_l1b_native.html
    reader = "seviri_l1b_native"

    if as_euro:
        ##
        # Read raw_lons, raw_lats for CloudCast raw for lon/lat domain matching with MSG
        #   raw_lons (928, 1530): [-69.2706298828125  ... 69.2706298828125]
        #   raw_lats (928, 1530): [ 26.67105484008789 ... 81.09877014160156]
        with open(f"/Users/mbauer/tmp/CloudCast/raw_coords.pkl", 'rb') as f:
            tmp = pickle.load(f)
            raw_lons, raw_lats = tmp
            del tmp
        if verbose:
            tmp = raw_lons.flatten()
            tmp = tmp[np.abs(tmp) <= 180.0]
            print(f"\traw_lons {raw_lons.shape}: [{np.amin(tmp)} ... {np.amax(tmp)}]")
            tmp = raw_lats.flatten()
            tmp = tmp[np.abs(tmp) <= 90.0]
            print(f"\traw_lats {raw_lats.shape}: [{np.amin(tmp)} ... {np.amax(tmp)}]")
    else:
        raw_lons = np.zeros((1,1))
        raw_lats = np.zeros((1,1))

    ##
    # Run MSG setup
    geo_stuff = setup_msg(use_dataset, as_ccast, as_euro, as_merc, as_lcc)
    (MSG_EPSG, MERC_EPSG, LCC_EPSG, geod_crs,
     ccast_area_def, ccast_crs, ccast_merc_area_def, ccast_merc_crs, ccast_lcc_area_def, ccast_lcc_crs,
     msg_area_def, msg_crs, msg_merc_area_def, msg_merc_crs, raw_area_def, raw_crs, raw_merc_area_def, raw_merc_crs) = geo_stuff
    CCAST_HEIGHT = 768
    CCAST_WIDTH = 768
    # labels = ("MSG_EPSG", "MERC_EPSG", "geod_crs", "ccast_area_def", "ccast_crs", "ccast_merc_area_def", "ccast_lcc_area_def", "ccast_lcc_crs",
    #  "ccast_merc_crs", "msg_area_def", "msg_crs", "msg_merc_area_def", "msg_merc_crs")
    # for ii, ival in enumerate(tmp):
    #     print(f"\n{labels[ii]}: {ival}")
    # return

    ##
    # Interpolate ERA-5 Temperature Fields to match MSG
    if mod_t_files:
        read_spatial = [False, True][0]
        do_stats = [False, True][0]
        make_animation = [False, True][0]
        fuze_years = [False, True][1]
        if make_animation:
            read_spatial = False
            do_stats = False
            fuze_years = False

        ##
        # Run ERA5 setup
        do_dir = f"{t_src_dir}{t_years[0]:4d}/{t_levs[0]}/"
        do_file = sorted([f"{do_dir}{_}" for _ in os.listdir(do_dir) if _.endswith('.nc')])[0]
        output = setup_era5(do_file)
        (era5_lons, era5_lats, era5_lat_edges, era5_lon_edges) = output
        del output

        ##
        # The ERA-5 mesh to interpolate from
        grid_lons, grid_lats = np.meshgrid(era5_lons, era5_lats)
        sparse_points = np.stack([grid_lons.ravel(), grid_lats.ravel()], -1)  # shape (N, 2) in 2d

        ##
        # Read MSG lon and lats from file (saved in natread.py)
        #   fine mesh to interpolate into
        with open(lonlatfile, 'rb') as f:
            ccast_lons = np.load(f)
            ccast_lats = np.load(f)
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
            del tmp

        # done_this = [('T850', 2017), ('T850', 2018),
        #              ('T700', 2017), ('T700', 2018),
        #              ('T500', 2017), ('T500', 2018),
        #              ('T250', 2017), ('T250', 2018)]
        done_this = []
        # done_this = [('T850', 2017)]
        for do_lev in t_levs:
            for do_yr in t_years:
                check_done = (do_lev, do_yr)
                if check_done in done_this:
                    print(f"\tDid {check_done}")
                    continue
                r"""
                    ERA-5 Latitudes  161: [ +70.000, ...  +30.000]
                    ERA-5 Longitudes 241: [ -20.000, ...  +40.000]

                    Grid Spacing 0.25 x 0.25 deg
                                 ~30 x 30 km

                    Grid Centers
                        (70, -20)----------------(70, 40)
                        |                               |
                        |                               |
                        |                               |
                        |                               |
                        (30, -20)----------------(30, 40)

                    Grid Edges
                        (70.125, -20.125)------+------(70.125, -19.875)-----------------(70.125, 39.875)------+------(70.125, 40.125)
                                            (70, -20)                                                      (70, 40)
                        (69.875, -20.125)------+------(69.875, -19.875)                 (69.875, 39.875)------+------(69.875, 40.125)
                        |                                                                                                           |
                        |                                                                                                           |
                        |                                                                                                           |
                        |                                                                                                           |
                        |                                                                                                           |
                        (30.125, -20.125)------+------(30.125, -19.875)                 (30.125, 39.875)------+------(30.125, 40.125)
                                            (30, -20)                                                      (30, 40)
                        (29.875, -20.125)------+------(29.875, -19.875)-----------------(29.875, 39.875)------+------(29.875, 40.125)

                    For interpolation we will project the data to put it in a linear Cartesian framework.

                    2017 ERA-5 850 hPa Temperature (8760, 161, 241): min  247.143, mean  277.832, max  309.442 K

                """
                do_dir = f"{t_src_dir}{do_yr:4d}/{do_lev}/"
                do_file = sorted([f"{do_dir}{_}" for _ in os.listdir(do_dir) if _.endswith('.nc')])[0]
                print(f"\tReading {do_file}")

                spatial_file = do_file.replace(".nc", "_spatial.npz")
                final_file = do_file.replace(".nc", "_ccast.npz")

                if fuze_years:
                    if do_yr == 2018:
                        continue
                    first_step = 35036
                    first_file = spatial_file.replace(".npz", f"_{first_step:04d}.npz")
                    # print(f"Reading {first_file}")
                    b_ds = np.load(first_file)
                    first_temp = b_ds['arr_0']
                    del b_ds

                    last_file = spatial_file.replace(".npz", f"_0000.npz")
                    last_file = last_file.replace("2017", "2018")
                    # print(f"Reading {last_file}")
                    b_ds = np.load(last_file)
                    last_temp = b_ds['arr_0']
                    del b_ds

                    n15_file = first_file.replace(f"{first_step:04d}", f"{first_step + 1:04d}")
                    n30_file = first_file.replace(f"{first_step:04d}", f"{first_step + 2:04d}")
                    n45_file = first_file.replace(f"{first_step:04d}", f"{first_step + 3:04d}")
                    # print(f"{n15_file} -> {n30_file} -> {n45_file}")

                    temp_15 = np.zeros((CCAST_HEIGHT, CCAST_WIDTH), dtype=float)
                    temp_30 = np.zeros((CCAST_HEIGHT, CCAST_WIDTH), dtype=float)
                    temp_45 = np.zeros((CCAST_HEIGHT, CCAST_WIDTH), dtype=float)

                    ##
                    # Pixel by Pixel linear fit over time insert 15, 30, 45
                    for jj in range(CCAST_HEIGHT):
                        for ii in range(CCAST_WIDTH):
                            start_temp = first_temp[jj, ii]
                            end_temp = last_temp[jj, ii]
                            slope = end_temp - start_temp
                            # print(f"\t\t{jj = :3d} {ii = :3d} {start_temp:10.5f} -> {end_temp:10.5f}\tslope = {slope}")
                            ##
                            # Linear Fit
                            temp_15[jj, ii] = (slope * 0.25) + start_temp
                            temp_30[jj, ii] = (slope * 0.5) + start_temp
                            temp_45[jj, ii] = (slope * 0.75) + start_temp
                            # print(f"\t\t\t{start_temp:10.5f} {temp_15[jj, ii]:10.5f} {temp_30[jj, ii]:10.5f} {temp_45[jj, ii]:10.5f} {end_temp:10.5f}")
                            # break
                        # break
                    np.savez_compressed(n15_file, temp_15)
                    np.savez_compressed(n30_file, temp_30)
                    np.savez_compressed(n45_file, temp_45)

                    if do_lev == 'T250':
                        return
                    else:
                        continue

                ##
                # Read raw Netcdf file
                ncfile = netCDF4.Dataset(do_file, mode='r', format='NETCDF4_CLASSIC')
                era5_temp = np.squeeze(ncfile.variables["t"][:], axis=1)
                ncfile.close()
                tsteps = era5_temp.shape[0]
                ccast_tsteps = tsteps * 4
                nlats = len(era5_lats)
                nlons = len(era5_lons)
                if verbose or do_stats:
                    print(f"\n\tERA-5 Latitudes                           ({nlats}): [{era5_lats[-1]:+8.3f}, ... {era5_lats[0]:+8.3f}]")
                    print(f"\tERA-5 Longitudes                          ({len(era5_lons)}): [{era5_lons[0]:+8.3f}, ... {era5_lons[-1]:+8.3f}]")
                    print(f"\t{do_yr} ERA-5 {do_lev[1:]} hPa Temperature {era5_temp.shape}: min {np.amin(era5_temp):8.3f}, mean {np.mean(era5_temp):8.3f}, max {np.amax(era5_temp):8.3f} K")
                    print(f"\tTime Steps                                     : {tsteps}")
                ##
                # Remove Mask
                era5_temp = ma.compressed(era5_temp)
                oned_shape = era5_temp.shape
                if oned_shape == (nlons * nlats * tsteps,):
                    era5_temp = np.reshape(era5_temp, shape=(tsteps, nlats, nlons))
                else:
                    raise Exception(f"Shape error {oned_shape} != {(nlons * nlats * tsteps,)}")

                if do_stats:
                    ##
                    # Trim ERA-5 so roughly same domain as CCAST
                    #    40 < era5_lats < 63
                    #   -13 < era5_lons < 34
                    jtmp = [ij for ij, _ in enumerate(era5_lats) if 40.0 < _ < 63.0]
                    # print(jtmp)
                    # print(era5_lats[jtmp[0]], era5_lats[jtmp[-1]])
                    itmp = [ij for ij, _ in enumerate(era5_lons) if -13.0 < _ < 34.0]
                    # print(itmp)
                    # print(era5_lons[itmp[0]], era5_lons[itmp[-1]])
                    era5_temp = era5_temp[:, jtmp[0]:jtmp[-1] + 1, itmp[0]:itmp[-1] + 1]

                    era5_range = (float(np.amin(era5_temp)), float(np.amax(era5_temp)))
                    era5_q75, era5_q25 = np.percentile(era5_temp, [75 ,25])
                    era5_iqr = float(era5_q75 - era5_q25)
                    era5_mean = float(np.mean(era5_temp))
                    era5_median = float(np.median(era5_temp))
                    print(f"\tRange  : ({era5_range[0]:8.3f}, {era5_range[1]:8.3f} {era5_range[1] - era5_range[0]:8.3f} K")
                    print(f"\tQ75    : {era5_q75:8.3f} K")
                    print(f"\tQ25    : {era5_q25:8.3f} K")
                    print(f"\tIQR    : {era5_iqr:8.3f} K")
                    print(f"\tMean   : {era5_mean:8.3f} K")
                    print(f"\tMedian : {era5_median:8.3f} K")
                    print("\n Time Series")
                    for tidx in range(tsteps):
                        zvals = era5_temp[tidx, :, :]
                        zvals_range = (float(np.amin(zvals)), float(np.amax(zvals)))
                        zvals_q75, zvals_q25 = np.percentile(zvals, [75 ,25])
                        zvals_iqr = float(zvals_q75 - zvals_q25)
                        zvals_mean = float(np.mean(zvals))
                        zvals_median = float(np.median(zvals))
                        print(f"\t{tidx:4d}\t{zvals_range[0]:8.3f}\t{zvals_q25:8.3f}\t{zvals_mean:8.3f}\t{zvals_median:8.3f}\t{zvals_q75:8.3f}\t{zvals_range[1]:8.3f}")

                    tmp = ccast_lats[np.abs(ccast_lats) <= 90.0]
                    print(f"\n\n\tCCAST Latitudes                           ({len(ccast_lats)}): [{np.amin(tmp):+8.3f}, ... {np.amax(tmp):+8.3f}]")
                    tmp = ccast_lons[np.abs(ccast_lons) <= 180.0]
                    print(f"\tCCAST Longitudes                          ({len(ccast_lons)}): [{np.amin(tmp):+8.3f}, ... {np.amax(tmp):+8.3f}]")
                    print(f"\tTime Steps                                     : {ccast_tsteps}")
                    for tidx in range(ccast_tsteps):
                        s_file = spatial_file.replace(".npz", f"_{tidx:04d}.npz")
                        # print(f"Reading {s_file}")
                        b_ds = np.load(s_file)
                        ccast_temp = b_ds['arr_0']
                        del b_ds
                        ccast_range = (float(np.amin(ccast_temp)), float(np.amax(ccast_temp)))
                        ccast_q75, ccast_q25 = np.percentile(ccast_temp, [75 ,25])
                        ccast_iqr = float(ccast_q75 - ccast_q25)
                        ccast_mean = float(np.mean(ccast_temp))
                        ccast_median = float(np.median(ccast_temp))
                        # print(f"\t{tidx:4d}\t{ccast_range[0]:8.3f}\t{ccast_q25:8.3f}\t{ccast_mean:8.3f}\t{ccast_median:8.3f}\t{ccast_q75:8.3f}\t{ccast_range[1]:8.3f}")
                    return

                looper = list(range(ccast_tsteps))
                #     0    1    2    3    4    5    6   7   8   9
                #     R                   R                 R
                real_times = looper[::4]
                del looper
                if make_animation:
                    # ffmpeg -framerate 48 -pattern_type glob -i '*.png' -vcodec libx264 -pix_fmt yuv420p -s 1920x1080 -crf 0 movie.mp4
                    wide_domain = [False, True][0]
                    as_lcc = [False, True][0]
                    make_mesh = [False, True][0]

                    if check_done[0] == 'T850':
                        # From min/max of original ERA-5 data
                        minmax = (247, 305)

                    the_title = f"ERA-5 TSTEP of {do_yr} {do_lev[1:]} hPa as CCAST via CUBIC Interpolation"
                    # looper = range(ccast_tsteps) if verbose else tqdm(range(ccast_tsteps), total=len(list(range(ccast_tsteps))), desc=f"Time Animation...")
                    # for tidx in looper:
                    #     s_file = spatial_file.replace(".npz", f"_{tidx:04d}.npz")
                    #     b_ds = np.load(s_file)
                    #     temp_dat = b_ds['arr_0']
                    #     del b_ds
                    #     plot_era5_ccast(temp_dat, ccast_lons, ccast_lats, do_lev, do_yr, tidx, wide_domain,
                    #                     as_lcc, make_mesh, t_src_dir, the_title.replace("TSTEP", f"step {tidx:5d}"),
                    #                     minmax)

                    done_files = []
                    print("\tRun Animation...")
                    with Pool(N_CORES) as pool:
                        results = pool.starmap(plot_era5_ccast, ((spatial_file.replace(".npz", f"_{tidx:04d}.npz"), ccast_lons, ccast_lats, do_lev, do_yr, tidx, wide_domain, as_lcc, make_mesh, t_src_dir, the_title.replace("TSTEP", f"step {tidx:5d}"), minmax) for tidx in range(ccast_tsteps)), chunksize=None)
                        for res in results:
                            done_files.append(res)
                    del res, results
                    print(f"\tDone {len(done_files)} of {len(list(range(ccast_tsteps)))}")
                    print("Done Interpolation Animation")
                    return

                if do_stats or make_animation or fuze_years:
                    done_this.append(check_done)
                    continue

                ##
                # Interpolate ERA-5 to MSG CloudCast
                #   [tsteps, nlats, nlons] -> [ccast_tsteps, CCAST_HEIGHT, CCAST_WIDTH]
                #   [  8760,   161,   241] -> [       35040,          768,         768]
                if read_spatial:
                    del era5_temp
                    pass
                else:
                    ##
                    # Spatial Interpolate ERA-5 to MSG CloudCast
                    done_files = []
                    tmp_clons = ccast_lons.ravel()
                    tmp_clats = ccast_lats.ravel()
                    print("\tRun spatial_downscale...")
                    with Pool(N_CORES) as pool:
                        results = pool.starmap(spatial_downscale, ((era5_temp[ridx, :, :], tidx, sparse_points, tmp_clons, tmp_clats, spatial_file, CCAST_HEIGHT, CCAST_WIDTH) for ridx, tidx in enumerate(real_times)), chunksize=None)
                        for res in results:
                            done_files.append(res)
                    del res, results
                    print(f"\tDone {len(done_files)} of {len(real_times)}")

                    # looper = enumerate(real_times) if verbose else tqdm(enumerate(real_times), total=len(real_times), desc=f"Spatial Interpolation...")
                    # for ridx, tidx in looper:
                    #     grid_centers_z_vals = era5_temp[ridx, :, :]
                    #     if verbose:
                    #         print(f"\n\t\t{ridx = :5d} {tidx = :4d} ERA-5 Temperature {grid_centers_z_vals.shape}: min {np.amin(grid_centers_z_vals):8.3f}, mean {np.mean(grid_centers_z_vals):8.3f}, max {np.amax(grid_centers_z_vals):8.3f} K")
                    #     ##
                    #     # Spatial Interpolation
                    #     use_method = ('nearest', 'linear', 'cubic')[2]
                    #     fine_mesh_vals = griddata(sparse_points, grid_centers_z_vals.ravel(), (ccast_lons.ravel(), ccast_lats.ravel(), CCAST_HEIGHT, CCAST_WIDTH), method=use_method)
                    #     ##
                    #     # Reshape
                    #     fine_mesh_vals = np.reshape(fine_mesh_vals, (CCAST_HEIGHT, CCAST_WIDTH))
                    #     if verbose:
                    #         print(f"\t\tInterpolated ERA-5 Temperature {fine_mesh_vals.shape}: min {np.amin(fine_mesh_vals):8.3f}, mean {np.mean(fine_mesh_vals):8.3f}, max {np.amax(fine_mesh_vals):8.3f} K")
                    #     ##
                    #     # Save tmp file
                    #     s_file = spatial_file.replace(".npz", f"_{tidx:04d}.npz")
                    #     # print(f"\tSaving {s_file}")
                    #     np.savez_compressed(s_file, fine_mesh_vals)
                    #     # if ridx > 10:
                    #     #     break
                    # del fine_mesh_vals

                ##
                # Loop over time to insert 15 minute entries via linear fit.
                looper = real_times if verbose else tqdm(real_times, total=len(real_times), desc=f"Temporal Interpolation...")
                for tidx in looper:
                    s_file = spatial_file.replace(".npz", f"_{tidx:04d}.npz")
                    b_ds = np.load(s_file)
                    at_time_end = b_ds['arr_0']
                    del b_ds
                    if tidx == 0:
                        at_time_start = at_time_end
                        continue
                    oidx = tidx - 4
                    # print(f"\t{oidx = :4d} -> {tidx:4d}")

                    n15_file = s_file.replace(f"{tidx:04d}", f"{oidx + 1:04d}")
                    n30_file = s_file.replace(f"{tidx:04d}", f"{oidx + 2:04d}")
                    n45_file = s_file.replace(f"{tidx:04d}", f"{oidx + 3:04d}")
                    # print(f"{n15_file} -> {n30_file} -> {n45_file}")

                    temp_15 = np.zeros((CCAST_HEIGHT, CCAST_WIDTH), dtype=float)
                    temp_30 = np.zeros((CCAST_HEIGHT, CCAST_WIDTH), dtype=float)
                    temp_45 = np.zeros((CCAST_HEIGHT, CCAST_WIDTH), dtype=float)
                    ##
                    # Pixel by Pixel linear fit over time insert 15, 30, 45
                    for jj in range(CCAST_HEIGHT):
                        for ii in range(CCAST_WIDTH):
                            start_temp = at_time_start[jj, ii]
                            end_temp = at_time_end[jj, ii]
                            slope = end_temp - start_temp
                            # print(f"\t\t{jj = :3d} {ii = :3d} {start_temp:10.5f} -> {end_temp:10.5f}\tslope = {slope}")
                            ##
                            # Linear Fit
                            temp_15[jj, ii] = (slope * 0.25) + start_temp
                            temp_30[jj, ii] = (slope * 0.5) + start_temp
                            temp_45[jj, ii] = (slope * 0.75) + start_temp
                            # print(f"\t\t\t{start_temp:10.5f} {temp_15[jj, ii]:10.5f} {temp_30[jj, ii]:10.5f} {temp_45[jj, ii]:10.5f} {end_temp:10.5f}")
                            # break
                        # break
                    np.savez_compressed(n15_file, temp_15)
                    np.savez_compressed(n30_file, temp_30)
                    np.savez_compressed(n45_file, temp_45)
                    at_time_start = at_time_end

                done_this.append(check_done)
            #     break  # year
            # break # lev
        print("Done mod_t_files")
        return

    if read_nat:
        ##
        # Read the file
        natread(fname=FNAME, fvar=use_dataset, reader=reader, to_euro=as_euro, euro_lons=raw_lons, euro_lats=raw_lats, to_ccast=as_ccast, geo_dat=geo_stuff, fromzip=zip2nat)
        return

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
        """

        ##
        # Convert nat to GeoTIF
        nat2tif(fname=FNAME, fvar=use_dataset, reader=reader, outdir=SUB_PATH, label=use_dataset, atag=use_tag,
                geo_dat=geo_stuff, to_full=as_full, to_ccast=as_ccast, to_euro=as_euro, to_merc=as_merc, to_lcc=as_lcc)

    if make_fig:
        ##
        # Plot the tif made by make_tif or in case of as_merc modified with gdal
        #   $
        #   <class 'osgeo.gdal.Dataset'>
        #   ds.GetGeoTransform()  = (-855100.436345, 3000.0, 0.0, -2638000.0, 0.0, -3000.0)
        #   ds.GetProjection()    = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]'

        the_title = f"{use_dataset}"
        if use_nat:
            if zip2nat:
                zpath = Path(FNAME)
                zfile = f"{zpath.parent}.zip"
                ffile = f"{zpath.name}"
                with ZipFile(zfile, 'r') as zipper:
                    TNAME = zipper.extract(ffile)
            ##
            # Open the file w/ satpy, which uses Xarray
            #   <class 'satpy.scene.Scene'>
            scn = Scene(filenames = {reader: [TNAME]})

            ##
            # Get calibration
            scn.load([use_dataset])
            fvar__ = scn[use_dataset]
            fvar__atts = fvar__.attrs
            # print(fvar__atts)
            fvar_units = fvar__atts['units']
            fvar_name = fvar__atts['standard_name']
            use_calibration = fvar__atts['calibration']
            if verbose:
                print(f"\tcalibration '{use_calibration}' w/ units of {fvar_units}")
            del fvar__

            the_title = f"{use_dataset} as {use_calibration}"

            # calibration = "radiance"
            # dtype = "float32"
            # radius = 16000
            # epsilon = 0.5
            # nodata = -3.4E+38

            ##
            # Load the data, different calibration can be chosen
            scn.load([use_dataset], calibration=use_calibration)
            # scn.load([use_dataset])
            data_vals = scn[use_dataset].values
            if as_euro:
                data_vals = data_vals[-928:, 1092:2622]
            elif as_ccast:
                if as_merc:
                    use_crs = ccast_merc_crs
                    use_area_def = ccast_merc_area_def
                elif as_lcc:
                    use_crs = ccast_lcc_crs
                    use_area_def = ccast_lcc_area_def
                else:
                    use_crs = ccast_crs
                    use_area_def = ccast_area_def
                dtype = "float32"
                radius = 16000
                epsilon = 0.5
                nodata = -3.4E+38

                #   use_crs.bounds = (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)
                # print(f"{use_crs.bounds = }")

                r"""
                from setup_msg():
                    MSG_HEIGHT = 3712
                    MSG_WIDTH = 3712
                    lower_left_xy  = [-75.26545715332031, -78.95874786376953]
                    upper_right_xy = [75.56462097167969, 78.29975891113281]

                """
                #   lons (3712, 3712): -81.12548663687087 ... 81.12472104087496
                #   lats (3712, 3712): -81.07327794830798 ... 81.07426055288249
                lons, lats = scn[use_dataset].area.get_lonlats()
                print(f"lons {lons.shape}: {np.amin(lons)} ... {np.amax(lons[np.isfinite(lons)])}")
                print(f"lats {lats.shape}: {np.amin(lats)} ... {np.amax(lats[np.isfinite(lats)])}")

                first_lon = 0
                max_width = [0, 0]
                row_start = -1
                row_end = -1
                for jj in range(3712):
                    print(f"{jj = :3d}: {lats[0, jj]}")
                    col_cnt = 0
                    col_hit = 0
                    for ii in range(3712):
                        # if np.isfinite(lons[jj, ii]):
                        #     # print(f"\n{jj = }, {ii = }: {lons[jj, ii] = }")
                        #     col_cnt += 1
                        #     col_hit = 1
                        #     if row_start == -1:
                        #         row_start = jj
                        if np.isfinite(lons[ii, jj]):
                            # print(f"\t{jj = :3d}, {ii = :3d}: {float(lons[ii, jj])}")
                            col_cnt += 1
                            col_hit = 1
                            if row_start == -1:
                                row_start = jj
                    if not col_hit and row_end == -1 and row_start != -1:
                        row_end = jj
                    if col_cnt > max_width[1]:
                        max_width[0] = f"{jj = :3d}, {ii = :3d}"
                        max_width[1] = col_cnt
                    print(f"\t{col_cnt = }")

                    # if first_lon:
                    #     break
                print(f"\n{row_start = }, {row_end = }")
                print(f"{max_width = }")
                print("Done")
                return
                """
                [jj, ii]
                    row_start = 51, row_end = 3661              3712-3661=51
                    max_width = ['jj = 1807, ii = 3711', 3622]
                [ii, jj]
                    row_start = 45, row_end = 3667              3712-3667=45
                    max_width = ['jj = 1808, ii = 3711', 3610]
                """
                print(f"{lons[0, 0] = }")
                print(f"{lons[-1, 0] = }")
                print(f"{lons[0, -1] = }")
                print(f"{lons[-1, -1] = }")
                print()
                print(f"{lats[0, 0] = }")
                print(f"{lats[-1, 0] = }")
                print(f"{lats[0, -1] = }")
                print(f"{lats[-1, -1] = }")

                return

                #   xvals (3712, 3712): -5567248.017211914 ... 5567247.798461914
                #   yvals (3712, 3712): -5567247.798461914 ... 5567248.017211914
                xvals, yvals = scn[use_dataset].area.get_proj_coords()
                print(f"xvals {xvals.shape}: {np.amin(xvals)} ... {np.amax(xvals)}")
                print(f"yvals {yvals.shape}: {np.amin(yvals)} ... {np.amax(yvals)}")

                ##
                # Apply a swath definition for our output raster
                #   <class 'pyresample.geometry.SwathDefinition'>
                swath_def = pr.geometry.SwathDefinition(lons=lons, lats=lats)

                ##
                # Resample our data to the area of interest
                """
                data_vals  (768, 768): [211.08590698242188, ... 288.95416259765625]

                from setup_msg():
                    ccast_area_def.corners = [
                        (-12.920934886492649, 62.403066090517555), UL
                        (33.73865749382469,   60.15059617915855),  UR
                        (21.32880156090482,   40.92817004367345),  LR
                        (-4.802566482888071,  42.068097533886025)] LL

                    ccast_area_def.area_extent  (-855100.436345, -4942000.0, 1448899.563655, -2638000.0)
                    ccast_crs.bounds      (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)

                    lower_left_xy  = [-855100.436345, -4942000.0] => (-4.816534709314307, 42.053336570266744)
                    upper_right_xy = [1448899.563655, -2638000.0] => (33.77742545811928,  60.15622631725923)

                # Seems good
                ccast_lons (768, 768): [-12.939207254994848, ... 33.70159752849037]
                ccast_lats (768, 768): [40.93204238308001, ... 63.7329345658077]

                # Not quite the same as setup_msg()
                ccast_x (768, 768)   : [-604581.3292236328, ... 1666723.7994384766]
                ccast_y (768, 768)   : [3929027.9376220703, ... 5144191.18347168]


                """
                data_vals = pr.kd_tree.resample_nearest(swath_def, data_vals,
                                                        use_area_def,
                                                        radius_of_influence=radius, # in meters
                                                        epsilon=epsilon,
                                                        fill_value=False)

                tmp = data_vals.flatten()
                tmp = tmp[~np.isnan(tmp)]
                print(f"data_vals {data_vals.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
                del tmp

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
                tmp = ccast_lons.flatten()
                print(f"ccast_lons {ccast_lons.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")

                tmp = ccast_lats.flatten()
                print(f"ccast_lats {ccast_lats.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")

                ccast_x = pr.kd_tree.resample_nearest(swath_def, xvals,
                                                         use_area_def,
                                                         radius_of_influence=radius, # in meters
                                                         epsilon=epsilon,
                                                         fill_value=False)
                ccast_y = pr.kd_tree.resample_nearest(swath_def, yvals,
                                                         use_area_def,
                                                         radius_of_influence=radius, # in meters
                                                         epsilon=epsilon,
                                                         fill_value=False)
                tmp = ccast_x.flatten()
                print(f"ccast_x {ccast_x.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
                tmp = ccast_y.flatten()
                print(f"ccast_y {ccast_y.shape}: [{np.amin(tmp)}, ... {np.amax(tmp)}]")

                ##
                # Define some metadata
                """
                CCAST_HEIGHT = 768
                CCAST_WIDTH = 768
                lower_left_xy = [-855100.436345, -4942000.0]
                upper_right_xy = [1448899.563655, -2638000.0]

                ccast_x (768, 768)   : [-604581.3292236328, ... 1666723.7994384766]
                ccast_y (768, 768)   : [3929027.9376220703, ... 5144191.18347168]

                pixelWidth : 3000.0
                pixelHeight: -3000.0
                originX    : -855100.436345
                originY    : -2638000.0
                """
                cols = data_vals.shape[1]
                rows = data_vals.shape[0]
                pixelWidth = (use_area_def.area_extent[2] - use_area_def.area_extent[0]) / cols
                pixelHeight = (use_area_def.area_extent[1] - use_area_def.area_extent[3]) / rows
                originX = use_area_def.area_extent[0]
                originY = use_area_def.area_extent[3]
                print(f"pixelWidth : {pixelWidth}")
                print(f"pixelHeight: {pixelHeight}")
                print(f"originX    : {originX}")
                print(f"originY    : {originY}")


                # used_dtype = [int, float][1]
                # ccast_x = np.zeros((768, 768), dtype=used_dtype)
                # ccast_y = np.zeros((768, 768), dtype=used_dtype)
                # for jj in range(768):
                #     for ii in range(768):
                #         # [jj, ii]
                #         # newx, newy = scn[use_dataset].area.get_array_indices_from_lonlat(ccast_lons[jj, ii], ccast_lats[jj, ii])
                #         # Matches ccast_x, ccast_y above
                #         newx, newy = scn[use_dataset].area.get_projection_coordinates_from_lonlat(ccast_lons[jj, ii], ccast_lats[jj, ii])
                #         ccast_x[jj, ii] = newx
                #         ccast_y[jj, ii] = newy

                #         #  [ii, jj]
                #         # newx, newy = scn[use_dataset].area.get_array_indices_from_lonlat(ccast_lons[ii, jj], ccast_lats[ii, jj])
                #         # # Matches ccast_x, ccast_y above
                #         # newx, newy = scn[use_dataset].area.get_projection_coordinates_from_lonlat(ccast_lons[ii, jj], ccast_lats[ii, jj])
                #         # ccast_x[ii, jj] = newx
                #         # ccast_y[ii, jj] = newy

                # # [jj, ii]
                # lower_left_i = ccast_x[0, 0]
                # lower_left_j = ccast_y[0, 0]
                # lower_right_i = ccast_x[0, -1]
                # lower_right_j = ccast_y[0, -1]
                # upper_left_i = ccast_x[-1, 0]
                # upper_left_j = ccast_y[-1, 0]
                # upper_right_i = ccast_x[-1, -1]
                # upper_right_j = ccast_y[-1, -1]

                # # # [ii, jj]
                # # lower_left_i = ccast_x[0, 0]
                # # lower_left_j = ccast_y[0, 0]
                # # lower_right_i = ccast_x[-1, 0]
                # # lower_right_j = ccast_y[-1, 0]
                # # upper_left_i = ccast_x[0, -1]
                # # upper_left_j = ccast_y[0, -1]
                # # upper_right_i = ccast_x[-1, -1]
                # # upper_right_j = ccast_y[-1, -1]

                # print(f"upper_left ({upper_left_i}, {upper_left_j}) .... ({upper_right_i}, {upper_right_j}) upper_right")
                # print(f"lower_left ({lower_left_i}, {lower_left_j}) .... ({lower_right_i}, {lower_right_j}) lower_right")
                # print()
                # print(f"ccast_x: [{np.amin(ccast_x)}, ... {np.amax(ccast_x)}]")
                # print(f"ccast_y: [{np.amin(ccast_y)}, ... {np.amax(ccast_y)}]")

                # # [ii, jj]
                # #   upper_left (1323, 3506) .... (1300, 3165) upper_right
                # #
                # #   lower_left (2057, 3553) .... (1982, 3204) lower_right
                # #
                # #   ccast_x: [1300, ... 2057] 2057-1300=757
                # #   ccast_y: [3165, ... 3570] 3570-3165=405
                # #
                # #   upper_left (1597714, 4952165) .... (1666723, 3929027) upper_right
                # #   lower_left (-604581, 5093184) .... (-379551, 4046043) lower_right
                # #
                # #   ccast_x: [-604581, ... 1666723]
                # #   ccast_y: [3929027, ... 5144191]

                # # [jj, ii]
                # #   upper_left (1982, 3204) .... (1300, 3165) upper_right
                # #   lower_left (2057, 3553) .... (1323, 3506) lower_right
                # #
                # #   ccast_x: [1300, ... 2057] 2057-1300=757
                # #   ccast_y: [3165, ... 3570] 3570-3165=405
                # #
                # #   upper_left (-379551.09851074027, 4046043.657592753) .... (1666723.7994384738, 3929027.937622063) upper_right
                # #   lower_left (-604581.3292236318,  5093184.331176753) .... (1597714.5286865155, 4952165.386596654) lower_right
                # #
                # #   ccast_x: [-604581.329223633, ... 1666723.7994384798]
                # #   ccast_y: [3929027.9376220573, ... 5144191.18347169]

                return

                # data_vals = data_vals.astype(np.uint16)
                # print(data_vals)
                # the_title = f"{use_dataset} as {use_calibration} and int"

                # print(data_vals)
            if zip2nat:
                os.remove(TNAME)
        else:
            if as_full:
                if as_merc:
                    tif_name = TNAME.replace(".tif", "_merc.tif")
                else:
                    tif_name = TNAME
            elif as_euro:
                if as_merc:
                    tif_name = TNAME.replace(".tif", "_merc.tif")
                else:
                    tif_name = TNAME
            elif as_ccast:
                if as_merc:
                    tif_name = TNAME.replace(".tif", "_merc.tif")
                elif as_lcc:
                    tif_name = TNAME.replace(".tif", "_lcc.tif")
                else:
                    tif_name = TNAME

            ds = gdal.Open(tif_name)
            band = ds.GetRasterBand(1)
            data_vals = band.ReadAsArray()

        # print(use_dataset, np.amin(data_vals), np.amax(data_vals))
        # tmp = data_vals.flatten()
        # tmp = tmp[~np.isnan(tmp)]
        # print(np.amin(tmp), np.amax(tmp))
        # return

        if as_region:
            data_vals = np.where(data_vals >= 0, 1, 0)
        if as_full:
            if as_merc:
                use_crs = msg_merc_crs
                use_area_def = msg_merc_area_def
            else:
                use_crs = msg_crs
                use_area_def = msg_area_def
        elif as_euro:
            if as_merc:
                use_crs = raw_merc_crs
                use_area_def = raw_merc_area_def
            else:
                use_crs = raw_crs
                use_area_def = raw_area_def
        elif as_ccast:
            if as_merc:
                use_crs = ccast_merc_crs
                use_area_def = ccast_merc_area_def
            elif as_lcc:
                use_crs = ccast_lcc_crs
                use_area_def = ccast_lcc_area_def
            else:
                use_crs = ccast_crs
                use_area_def = ccast_area_def
        use_extent = use_crs.bounds

        fig = plt.figure(figsize=(10, 8))
        as_geos = False
        if as_full:
            # ax = plt.axes(projection=ccrs.Miller())
            # ax = plt.axes(projection=ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0))
            # ax = plt.axes(projection=ccrs.Mercator(central_longitude=0.0, min_latitude=-80.0, max_latitude=80.0))
            #ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0))
            ax = plt.axes(projection=ccrs.Geostationary(central_longitude=0.0)); as_geos = True

            # proj = ccrs.PlateCarree(central_longitude=0.0, globe=ccrs.Globe(datum='WGS84', ellipse='WGS84'))
            # ax = plt.axes(projection=proj)
            # ax = plt.axes(projection=use_crs)
        elif as_euro:
            # ax = plt.axes(projection=ccrs.Miller())
            # ax = plt.axes(projection=ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0))
            # ax = plt.axes(projection=ccrs.Mercator(central_longitude=0.0, min_latitude=-80.0, max_latitude=80.0))
            #ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0))
            # ax = plt.axes(projection=ccrs.Geostationary(central_longitude=0.0)); as_geos = True

            ax = plt.axes(projection=use_crs)
            # proj = ccrs.PlateCarree(central_longitude=0.0, globe=ccrs.Globe(datum='WGS84', ellipse='WGS84'))
            # ax = plt.axes(projection=proj)
            # ax = plt.axes(projection=use_crs)
        elif as_ccast:
            ax = plt.axes(projection=use_crs)
            # ax = plt.axes(projection=ccrs.LambertConformal(central_longitude=0.0, central_latitude=0.0, cutoff=-30))
        else:
            ax = plt.axes(projection=use_crs)

        if as_ccast:
            #   CCAST Sterographic
            #       values: [0.5235565900802612, ... 13.425487518310547]
            #   CCAST Mercator
            #       values: [0.5609534978866577, ... 14.248218536376953]
            the_alpha = 1.0
            # if use_dataset == "HRV":
            #     # bounds = np.linspace(0.0, 15.0, num=16, endpoint=True)
            #     bounds = np.linspace(0.0, 25.0, num=16, endpoint=True)
            #     the_min = 0.0
            #     the_max = 25.0
            # else:
            #     bounds = np.linspace(0.0, 4.0, num=16, endpoint=True)
            #     the_min = 0.0
            #     the_max = 4.0

            #     bounds = np.linspace(0.0, 60.0, num=16, endpoint=True)
            #     the_min = 0.0
            #     the_max = 60.0

            # if as_region:
            #     bounds = np.linspace(0.0, 1.0, num=2, endpoint=True)
            #     the_min = 0.0
            #     the_max = 1.0
            #     the_alpha = 0.25
            the_min = np.round(np.amin(data_vals), decimals=0)
            the_max = np.round(np.amax(data_vals), decimals=0)
            bounds = np.linspace(the_min, the_max, num=16, endpoint=True)

            if as_lcc:
                a_image = plt.imshow(data_vals, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=the_alpha)
                ax.set_extent(use_extent, crs=use_crs)
            else:
                a_image = plt.imshow(data_vals, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=the_alpha)
                ax.set_extent(use_extent, crs=use_crs)
        elif as_merc:
            # ax.set_extent(use_extent, crs=geod_crs)
            a_image = plt.imshow(data_vals, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper')
            ax.set_extent(use_extent, crs=use_crs)
        else:
            if use_dataset == "HRV":
                use_extent = [-2750000, 2750000, -5500000, 5500000]
            else:
                if as_euro:
                    use_extent = (-2296808.8, 2293808.2, 5570249.0, 2785874.8)
                else:
                    use_extent = [-5500000, 5500000, -5500000, 5500000]
            # use_extent = [-6500000, 6500000, -6500000, 6500000]
            # use_extent = [-6500000, 6500000, -6500000, 6500000]
            # print(use_extent)

            #   CCAST Sterographic
            #       values: [0.5235565900802612, ... 13.425487518310547]
            #   CCAST Mercator
            #       values: [0.5609534978866577, ... 14.248218536376953]
            the_alpha = 1.0
            if use_dataset == "HRV":
                # bounds = np.linspace(0.0, 15.0, num=16, endpoint=True)
                the_min = 0.0
                the_max = 25.0
                bounds = np.linspace(the_min, the_max, num=11, endpoint=True)
            elif use_dataset == "IR_108":
                the_min = 0.0
                the_max = 120.0
                bounds = np.linspace(the_min, the_max, num=20, endpoint=True)
                if use_nat:
                    # print(use_dataset, np.amin(data_vals), np.amax(data_vals))
                    # tmp = data_vals.flatten()
                    # tmp = tmp[~np.isnan(tmp)]
                    # print(use_dataset, np.amin(tmp), np.amax(tmp))
                    # IR_108 211.34976 305.07642
                    the_min = 200
                    the_max = 310
                    bounds = np.linspace(the_min, the_max, num=23, endpoint=True)
            elif use_dataset == "IR_120":
                the_min = 0.0
                the_max = 140.0
                bounds = np.linspace(the_min, the_max, num=20, endpoint=True)
                if use_nat:
                    # print(use_dataset, np.amin(data_vals), np.amax(data_vals))
                    # tmp = data_vals.flatten()
                    # tmp = tmp[~np.isnan(tmp)]
                    # print(use_dataset, np.amin(tmp), np.amax(tmp))
                    # # IR_120 202.00209 305.52527
                    the_min = 200
                    the_max = 310
                    bounds = np.linspace(the_min, the_max, num=23, endpoint=True)
            elif use_dataset == "VIS006":
                the_min = 0.0
                the_max = 10.0
                bounds = np.linspace(the_min, the_max, num=10, endpoint=True)
                if use_nat:
                    print(use_dataset, np.amin(data_vals), np.amax(data_vals))
                    tmp = data_vals.flatten()
                    tmp = tmp[~np.isnan(tmp)]
                    print(use_dataset, np.amin(tmp), np.amax(tmp))
                    the_min = 0
                    the_max = 64
                    bounds = np.linspace(the_min, the_max, num=17, endpoint=True)
            elif use_dataset == "IR_039":
                the_min = 0.0
                the_max = 1.0
                bounds = np.linspace(the_min, the_max, num=11, endpoint=True)
                if use_nat:
                    # print(use_dataset, np.amin(data_vals), np.amax(data_vals))
                    # tmp = data_vals.flatten()
                    # tmp = tmp[~np.isnan(tmp)]
                    # print(use_dataset, np.amin(tmp), np.amax(tmp))
                    # # IR_039 204.75961 310.713
                    the_min = 200
                    the_max = 310
                    bounds = np.linspace(the_min, the_max, num=23, endpoint=True)
            else:
                the_min = 0.0
                the_max = 4.0
                bounds = np.linspace(the_min, the_max, num=16, endpoint=True)
                # print(use_dataset, np.amin(data_vals), np.amax(data_vals))
            if as_region:
                bounds = np.linspace(0.0, 1.0, num=2, endpoint=True)
                the_min = 0.0
                the_max = 1.0
                the_alpha = 0.25
            a_image = plt.imshow(data_vals, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=the_alpha)
        ax.add_feature(cfeature.COASTLINE, alpha=1)
        ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=2, color='black', alpha=0.5, linestyle='--')

        plt.title(the_title)

        if not as_region:
            ##
            # Color bar
            use_shrink = 0.9
            if not as_geos:
                use_shrink = 0.5

            cbar = fig.colorbar(a_image, location="right", orientation='vertical', extend="neither", ticks=bounds, shrink=use_shrink,
                                spacing='uniform', fraction=0.15, pad=0.1, aspect=30, drawedges=False,
                                ax=ax)
            cbar.ax.tick_params(labelsize=8)
            cbar.solids.set(alpha=1)

        plt.show()
        # To trim region

    return

#---Start of main code block.
if __name__== '__main__':

    ##
    # Convert a nat file to geotif
    make_tif = [False, True][0]

    ##
    # Read the nat file only (not geotif or plots)
    read_nat = [False, True][0]

    ##
    # Make a figure
    make_fig = [False, True][0]

    ##
    # Read nat rather than tif
    use_nat = [False, True][0]

    ##
    # Read .nat from .zip (use with use_nat and read_nat)
    zip2nat = [False, True][1]

    ##
    # Select Domain (only one can be True)

    ##
    # Return MSG full frame data (no reprojection, subsetting or interpolation)
    as_full = [False, True][0]

    ##
    # Return the CloudCast european resolution and domain (reprojection, subsetting and interpolation)
    as_euro = [False, True][0]

    ##
    # Return the CloudCast resolution and domain (reprojection, subsetting and interpolation)
    as_ccast = [False, True][1]

    if sum([int(as_full), int(as_euro), int(as_ccast)]) > 1:
        raise Exception("Can only use one of as_full as_euro as_ccast.")
    if sum([int(as_full), int(as_euro), int(as_ccast)]) == 0:
        raise Exception("Must select one of as_full as_euro as_ccast.")

    ##
    # Return using a Mercator projection (vs a Sterographic projection) (reprojection and interpolation, no subsetting)
    as_merc = [False, True][0]

    #
    # Return using a Lambert Conformal projection (vs a Sterographic projection) (reprojection and interpolation, no subsetting)
    #   Doesn't work for full because of latitude span
    as_lcc = [False, True][0]
    if as_lcc and as_full:
        raise Exception("Due to latitude span as_full can't be used with as_lcc.")

    ##
    # Tag to id new files
    use_tag = 'full' if as_full else 'ccast' if as_ccast else 'euro' if as_euro else ''
    if as_merc:
        use_tag = f"{use_tag}_merc"
    if as_lcc:
        use_tag = f"{use_tag}_lcc"

    ##
    # Display region (not data values)
    as_region = [False, True][0]

    ##
    # Obvious
    verbose = [False, True][0]

    ##
    # From scn.all_dataset_names() below
    ds_names = ("HRV", "IR_016", "IR_039", "IR_087", "IR_097", "IR_108", "IR_120",
                "IR_134", "VIS006", "VIS008", "WV_062", "WV_073")

    ##
    # Name(s) of MSG variable to work with.
    #   These are the base cloudcast fields  ["VIS006", "IR_039", "IR_108", "IR_120"]

    ##
    # Just HRV  Window/Water Vapor Channel
    # use_dataset = ds_names[ds_names.index("HRV")]

    ##
    # Just IR_039 IR Window Channel
    # use_dataset = ds_names[ds_names.index("IR_039")]

    ##
    # Just IR_087 IR Window Channel
    use_dataset = ds_names[ds_names.index("IR_087")]

    ##
    # Just WV_073 IR Water vapor
    # use_dataset = ds_names[ds_names.index("WV_073")]

    ##
    # Just IR_108 IR Window Channel
    # use_dataset = ds_names[ds_names.index("IR_108")]

    ##
    # Just IR_120 IR Window Channel
    # use_dataset = ds_names[ds_names.index("IR_120")]

    ##
    # Just VIS006 Visible Window Channel
    # use_dataset = ds_names[ds_names.index("VIS006")]

    ##
    # Just VIS008 Visible Window Channel
    # use_dataset = ds_names[ds_names.index("VIS008")]

    ##
    # Define path to folder
    FILE_PATH = "/Users/mbauer/tmp/CloudCast/msg/"
    BASENAME = ["MSG3-SEVI-MSG15-0100-NA-20170601001240.835000000Z-NA",
                "MSG3-SEVI-MSG15-0100-NA-20170102002740.606000000Z-NA",
                "MSG3-SEVI-MSG15-0100-NA-20170102122740.989000000Z-NA",
                "MSG3-SEVI-MSG15-0100-NA-20170104005740.099000000Z-NA"][0]
    SUB_PATH = f"{FILE_PATH}{BASENAME}/"
    FNAME = f"{SUB_PATH}{BASENAME}.nat"
    TNAME = f"{SUB_PATH}{BASENAME}_{use_dataset}_{use_tag}.tif"
    ZNAME = f"{FILE_PATH}{BASENAME}.zip"
    if use_nat:
        TNAME = FNAME
    if verbose:
        print(f"{'Making' if make_tif else 'Reading'}: {TNAME}")

    ##
    # Color map for imaging
    freq_map_cmap = ["plasma", "gist_ncar_r", "bone_r", "nipy_spectral"][2]
    if as_region:
        freq_map_cmap = "bone_r"
    if use_nat:
        freq_map_cmap = "bone"

    mod_t_files = [False, True][1]
    t_levs = ["T850", "T700", "T500", "T250"]
    t_src_dir = "/Volumes/saved4/ERA5/"
    t_years = [2017, 2018]
    # fine mesh to interpolate into  saved in natread.py
    lonlatfile = "/Users/mbauer/tmp/CloudCast/ccast_lonlat.npy"

    sys.exit(main(make_tif, read_nat, make_fig, use_nat, zip2nat, as_full, as_euro,
                  as_ccast, as_merc, as_lcc, as_region, verbose, mod_t_files, use_tag,
                  use_dataset, freq_map_cmap, SUB_PATH, FNAME, TNAME, ZNAME, t_src_dir,
                  lonlatfile, t_levs, t_years))
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
