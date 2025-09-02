#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""Rearrange lon/lat to work with matplotlib pcolormesh.

plot_era5_ccast.py
~~~~~~~~~~~~~~~~~~
"""

# Standard Imports
import os
import sys

# Third-Party Imports
import numpy as np
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pyresample as pr

# Local Imports
from cloudcast.plot.whitegreenblueyellow_colors import read_wgby_colors_from_json
from cloudcast.plot.fix_lon_lat_for_pcolormesh import fix_lon_lat_for_pcolormesh

##
# List of Public objects from this module.
__all__ = ["plot_era5_ccast"]

##
# Markup Language Specification (see Google Python Style Guide https://google.github.io/styleguide/pyguide.html)
__docformat__ = "Google en"
# ------------------------------------------------------------------------------

###############################################################################
# PUBLIC plot_era5_ccast()
# ------------------------
def plot_era5_ccast(temp_file, ccast_lons, ccast_lats, do_lev, do_yr, tidx, wide_domain, as_lcc, make_mesh, t_src_dir, the_title, minmax):

    b_ds = np.load(temp_file)
    temp_dat = b_ds['arr_0']
    del b_ds

    ##
    # MSG base projection WGS84 - World Geodetic System 1984
    #   https://epsg.io/4326
    #   +proj=longlat +datum=WGS84 +no_defs +type=crs
    msg_crs = ccrs.PlateCarree()

    ##
    # WGS-84 Earth equatorial radius at sea level (meters)
    #   STARE uses ccrs.Globe(datum='WGS84', ellipse='WGS84') https://epsg.io/4326
    #   Clarke 1866 ellipsoid https://epsg.io/7008-ellipsoid
    globe = ccrs.Globe(datum='WGS84', ellipse='WGS84')

    # Geodetic:
    #   A 3D/spherical CRS based on latitude and longitude where geographical distance and coordinates are measured in degrees.
    geod_crs = ccrs.Geodetic(globe=globe)

    if not 'whitegreenblueyellow_colors' in mpl.colormaps:
        whitegreenblueyellow_colors_data = read_wgby_colors_from_json()
        whitegreenblueyellow_colors = mcolors.LinearSegmentedColormap('whitegreenblueyellow_colors', segmentdata=whitegreenblueyellow_colors_data, N=256)
        mpl.colormaps.register(name='whitegreenblueyellow_colors', cmap=whitegreenblueyellow_colors)
        color_scheme = mpl.colormaps.get_cmap('whitegreenblueyellow_colors')
    else:
        color_scheme = mpl.colormaps.get_cmap('whitegreenblueyellow_colors')

    freq_map_cmap = ["plasma", "gist_rainbow", "jet", "brg", 'whitegreenblueyellow_colors'][4]

    if make_mesh:
        plon_mid, plat_mid = fix_lon_lat_for_pcolormesh(np.copy(ccast_lons), np.copy(ccast_lats))
        x_lon, y_lat = np.meshgrid(plon_mid, plat_mid)

    # pname = f"{t_src_dir}era5_2_msg_{do_lev}_cubic_{freq_map_cmap}.png"
    pname = f"{t_src_dir}era5_2_msg_{do_lev}_cubic_{tidx:05d}.png"
    if make_mesh:
        pname = pname.replace(".png", "_mesh.png")
    if as_lcc:
        pname = pname.replace(".png", "_lcc.png")
    if wide_domain:
        pname = pname.replace(".png", "_wide.png")

    ##
    # Create some information on the reference system
    CCAST_HEIGHT = 768
    CCAST_WIDTH = 768
    lower_left_xy = [-855100.436345, -4942000.0]
    upper_right_xy = [1448899.563655, -2638000.0]

    ##
    # Define the area
    ccast_area_def = pr.geometry.AreaDefinition('areaD', 'Europe', 'areaD',
                                                {'lat_0': '90.00', 'lat_ts': '50.00', 'lon_0': '5', 'proj': 'stere', 'ellps': 'WGS84'},
                                                CCAST_HEIGHT, CCAST_WIDTH,
                                                (lower_left_xy[0], lower_left_xy[1], upper_right_xy[0], upper_right_xy[1]))
    ##
    # Form a cartopy CRS
    ccast_crs = ccast_area_def.to_cartopy_crs()
    if as_lcc:
        ccast_lcc_crs = ccrs.LambertConformal(central_longitude=10.0, central_latitude=51.0, cutoff=0)
        area_id = "EPSG:9801"
        description = "Lambert"
        proj_id = "EPSG:9801"

        #   ll_xy: (-4.816534709314301, 42.053336570266744)
        #   ur_xy: (33.77742545811928,  60.15622631725923)
        ur_xy = ccast_lcc_crs.transform_point(upper_right_xy[0], upper_right_xy[1], ccast_crs)
        ll_xy = ccast_lcc_crs.transform_point(lower_left_xy[0], lower_left_xy[1], ccast_crs)
        projection = {'proj': 'lcc', 'lon_0': 10, 'lat_0': 51, 'lat_1': 40, 'lat_2': 60, 'type': 'crs'}
        ccast_lcc_area_def = pr.geometry.AreaDefinition(area_id, 'Europe', proj_id,
                                                        projection, CCAST_HEIGHT, CCAST_WIDTH,
                                                        (ll_xy[0], ll_xy[1], ur_xy[0], ur_xy[1]))
        ccast_lcc_crs = ccast_lcc_area_def.to_cartopy_crs()
        #              (min_lon,             max_lon,            min_lat,            max_lat)
        # use_extent = (-1216087.7976068966, 1404456.8408720866, -899612.0368218884, 1248209.190564286)

        # ll_lon: (-4.057268006224109, 38.90774971765646)
        # ur_lon: (36.732528675801106, 61.4087896613198)
        use_extent = ccast_lcc_crs.bounds
        ur_lon = geod_crs.transform_point(use_extent[1], use_extent[1], ccast_lcc_crs)
        ll_lon = geod_crs.transform_point(use_extent[0], use_extent[0], ccast_lcc_crs)
        # print(f"\tll_lon: ({float(ll_lon[0])}, {float(ll_lon[1])})")
        # print(f"\tur_lon: ({float(ur_lon[0])}, {float(ur_lon[1])})")
        use_crs = ccast_lcc_crs
    else:
        #   use_extent = (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)
        use_extent = ccast_crs.bounds
        use_crs = ccast_crs
    # print(f"\n\tuse_extent = ({float(use_extent[0])}, {float(use_extent[1])}, {float(use_extent[2])}, {float(use_extent[3])})")

    base_lcc_crs = ccrs.LambertConformal(central_longitude=10.0, central_latitude=51.0, cutoff=20)

    fig = plt.figure(figsize=(10, 8))
    if wide_domain:
        ax = plt.axes(projection=base_lcc_crs)
    else:
        ax = plt.axes(projection=ccast_crs)

    if wide_domain:
        # (min_lon, max_lon, min_lat, max_lat)
        if as_lcc:
            ax.set_extent([-10, 30, 35, 70], crs=msg_crs)
        else:
            ax.set_extent([-10, 30, 35, 70], crs=msg_crs)

    use_alpha = 1.0
    if do_lev == "T850":
        range_offset = 1
        npnts =  33
    elif do_lev == "T700":
        range_offset = 0
        npnts =  35
    elif do_lev == "T500":
        range_offset = 0
        npnts =  39
    elif do_lev == "T250":
        range_offset = 0
        npnts =  31
    else:
        range_offset = 0
        npnts =  33

    # the_min, the_max = minmax
    the_min = np.round(np.amin(temp_dat), decimals=0) + range_offset
    the_max = np.round(np.amax(temp_dat), decimals=0) - range_offset
    # print(the_min, the_max, range_offset)
    # for test_npnts in range(5, 40):
    #     bounds = np.linspace(the_min, the_max, num=test_npnts, endpoint=True)
    #     print(f"\nnpnts {test_npnts}: {bounds.tolist()}")
    # return

    bounds = np.linspace(the_min, the_max, num=npnts, endpoint=True)
    norm = mcolors.Normalize(vmin=bounds[0], vmax=bounds[-1])
    if make_mesh:
        # plon_mid, plat_mid = fix_lon_lat_for_pcolormesh(np.copy(lons), np.copy(lats_nh))
        # x_lon, y_lon = np.meshgrid(plon_mid, plat_mid)

        # param_dict = {"cmap": "Blues_r", "alpha": 1, "shading": 'nearest', "edgecolors": 'None', "norm": norm}
        param_dict = {"cmap": freq_map_cmap, "alpha": 1, "shading": 'nearest', "edgecolors": 'None', "norm": norm}
        aplt =  geo_axes.pcolormesh(x_lon, y_lat, ma.masked_equal(temp_dat, 0), transform=use_crs, zorder=1, **param_dict)
    else:
        # geo_axes.contourf(lons, lats_nh, ma.masked_equal(map_this, 0), 1, transform=ccrs.PlateCarree())
        i_ways = ['none', 'antialiased', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36',
                  'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel',
                  'mitchell', 'sinc', 'lanczos', 'blackman']
        iway = i_ways[1]
        # param_dict = {"extent": map_extent, "interpolation": iway, "cmap": color_scheme, "origin": 'lower', "norm": norm}
        # aplt = geo_axes.imshow(ma.masked_equal(map_this, 0), transform=map_crs, **param_dict)
        if as_lcc:
            a_image = plt.imshow(temp_dat, interpolation=iway, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=use_alpha)
        else:
            a_image = plt.imshow(temp_dat, interpolation=iway, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=use_alpha)

    # ax.set_extent(use_extent, crs=use_crs)
    # (min_lon, max_lon, min_lat, max_lat)
    # ax.set_extent([-10, 38, 35, 62], crs=msg_crs)

    ax.add_feature(cfeature.COASTLINE, alpha=1)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, linewidth=2, color='black', alpha=0.5, linestyle='--')
    plt.title(the_title)

    ##
    # Color bar
    use_shrink = 0.5
    cbar = fig.colorbar(a_image, location="right", orientation='vertical', extend="neither", ticks=bounds, shrink=use_shrink,
                        spacing='uniform', fraction=0.15, pad=0.1, aspect=30, drawedges=False,
                        ax=ax)
    cbar.ax.tick_params(labelsize=10)
    cbar.solids.set(alpha=1)
    cbar.set_ticks(bounds[0::4])
    cbar.set_label('Temperature [K]', rotation=270, labelpad=15)

    # plt.show()
    fig.savefig(pname, dpi=300, facecolor='w', edgecolor='w',
                orientation='landscape', bbox_inches='tight', pad_inches=0.02)
    plt.clf()
    plt.close('all')

    return

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
