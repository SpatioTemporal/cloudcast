#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Rewrite ERA-5 to MSG data

restruct_era.py
~~~~~~~~~~~~~~~

$ python restruct_era.py

$ cd /Volumes/saved4/ERA5/MLU

$ tar -cvf era5_2017.tar -C /Volumes/saved4/ERA5/MLU/2017 .
$ tar -cvf era5_2018.tar -C /Volumes/saved4/ERA5/MLU/2019 .

$ tar -tvf era5_2017.tar
$ tar -xvf era5_2017.tar -C path_to_expand_to

"""
# Standard Imports
import os
import sys
from pathlib import Path
from calendar import monthrange

# Third-Party Imports
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool

# Local Imports

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
# PUBLIC collapse_files()
# -----------------------
def collapse_files(cnames: list[list[str]], fname:str, tidx: int, xydim: int, nlevs: int, verbose: bool) -> str:
    ##
    # Cast final array
    temp_out = np.zeros((nlevs, xydim, xydim), dtype=np.float32)

    ##
    # Read Each level
    for ilev in range(nlevs):
        if verbose:
            print(f"\tReading {cnames[ilev][tidx]}")
        b_ds = np.load(cnames[ilev][tidx])
        in_temp = b_ds['arr_0']
        del b_ds
        temp_out[ilev, :, :] = in_temp[:, :]

    ##
    # Save Result
    np.savez_compressed(fname, temp_out)

    return fname

###############################################################################
# PUBLIC main()
# -------------
def main() -> None:

    verbose = [False, True][0]
    CCAST_HEIGHT = 768
    CCAST_WIDTH = 768
    t_levs = ["T850", "T700", "T500", "T250"]
    t_src_dir = "/Volumes/saved4/ERA5/"
    t_dst_dir = "/Volumes/saved4/ERA5/MLU/"
    t_years = [2017, 2018]

    ##
    # Loop over years
    for do_yr in t_years:
        ##
        # Make Timestamps (for filenames)
        timestamps = []
        for do_mm in range(1, 13):
            num_days = monthrange(do_yr, do_mm)[1]
            for do_dd in range(1, num_days + 1):
                for do_hh in range(24):
                    for do_min in (0, 15, 30, 45):
                        timestamps.append(f"{do_yr:04d}{do_mm:02d}{do_dd:02d}{do_hh:02d}{do_min:02d}00")
        if verbose:
            print(f"Created {len(timestamps)} timestamps [{timestamps[0]}... {timestamps[-1]}]")

        ##
        # Fetch file names for each level
        ccast_names_850 = sorted(list(Path(f"{t_src_dir}{do_yr:4d}/{t_levs[0]}/").rglob("*.npz")))
        ccast_names_700 = sorted(list(Path(f"{t_src_dir}{do_yr:4d}/{t_levs[1]}/").rglob("*.npz")))
        ccast_names_500 = sorted(list(Path(f"{t_src_dir}{do_yr:4d}/{t_levs[2]}/").rglob("*.npz")))
        ccast_names_250 = sorted(list(Path(f"{t_src_dir}{do_yr:4d}/{t_levs[3]}/").rglob("*.npz")))
        ccast_names = [ccast_names_850, ccast_names_700, ccast_names_500, ccast_names_250]

        ##
        # For each timestep, read all 4 levels, put in common array, form filename and save to file.
        nlevs = len(t_levs)
        nfiles = len(ccast_names_850)
        done_files = []
        print(f"\tConverting {do_yr} ...")
        with Pool(N_CORES) as pool:
            results = pool.starmap(collapse_files, ((ccast_names, f"{t_dst_dir}ERA5_temperature_{timestamps[tidx]}.npz", tidx, CCAST_HEIGHT, nlevs, verbose) for tidx in range(nfiles)), chunksize=None)
            for res in results:
                done_files.append(res)
        del res, results
        print(f"\tDone {len(done_files)} of {nfiles}")


        # looper = range(len(ccast_names_850)) if verbose else tqdm(range(len(ccast_names_850)), total=len(ccast_names_850), desc=f"Converting {do_yr} ...")
        # for tidx in looper:
        #     if verbose:
        #         print(f"Doing {tidx}: {timestamps[tidx]}")

        #     temp_out = np.zeros((4, CCAST_HEIGHT, CCAST_WIDTH), dtype=np.float32)

        #     ##
        #     # Read Each level
        #     for ilev in range(len(t_levs)):
        #         if verbose:
        #             print(f"\tReading {ccast_names[ilev][tidx]}")
        #         b_ds = np.load(ccast_names[ilev][tidx])
        #         in_temp = b_ds['arr_0']
        #         del b_ds
        #         temp_out[ilev, :, :] = in_temp[:, :]

        #     ##
        #     # Save Result
        #     fname = f"{t_dst_dir}ERA5_temperature_{timestamps[tidx]}.npz"
        #     if verbose:
        #         print(f"\n\tSaving {fname}")
        #     np.savez_compressed(fname, temp_out)

        #     # Done tidx
        #     # return

        # Done year
        # return

    # Done All
    return

#---Start of main code block.
if __name__== '__main__':

    sys.exit(main())

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<

