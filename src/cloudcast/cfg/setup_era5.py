#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Setup geometry for Read ERA-5.

setup_era5.py
~~~~~~~~~~~~
"""
# Standard Imports
import os
import typing
from typing import Optional, Union, TypeAlias

# Third-Party Imports
import numpy as np
import pyresample as pr
import cartopy
import cartopy.crs as ccrs
import netCDF4

# STARE Imports

# Local Imports

##
# List of Public objects from this module.
__all__ = ['setup_era5']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

def LON_TO_180(x): return ((x + 180.0) % 360.0) - 180.0
def LON_TO_360(x): return (x + 360.0) % 360.0

# ##
# # Type Aliases
# SetOut: TypeAlias = tuple[int, int, int, typing.TypeVar('cartopy.crs'),
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
#                      Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')]]

###############################################################################
# PUBLIC setup_era5()
# ------------------
def setup_era5(a_file: str, do_year: str, do_lev: str): #  -> SetOut:
    """
    Setup geometry for Meteosat Second Generation (MSG).

    Parameters
    ----------
    a_file : str
        Path to a ERA-5 to gather data

    Returns
    -------
    SetOut
        Tuple of geometry information for the specified projection and variable.
    """
    ##
    # MSG base projection WGS84 - World Geodetic System 1984
    #   https://epsg.io/4326
    #   +proj=longlat +datum=WGS84 +no_defs +type=crs
    MSG_EPSG = 32662
    msg_crs = ccrs.PlateCarree()

    # Lambert Conic Conformal (1SP) EPSG
    LCC_EPSG = 9801

    ##
    # Lambert Conic Conformal (2SP) EPSG
    LCC_EPSG = 9802

    ##
    # WGS-84 Earth equatorial radius at sea level (meters)
    #   STARE uses ccrs.Globe(datum='WGS84', ellipse='WGS84') https://epsg.io/4326
    #   Clarke 1866 ellipsoid https://epsg.io/7008-ellipsoid
    globe = ccrs.Globe(datum='WGS84', ellipse='WGS84')

    # Geodetic:
    #   A 3D/spherical CRS based on latitude and longitude where geographical distance and coordinates are measured in degrees.
    geod_crs = ccrs.Geodetic(globe=globe)

    verbose = [False, True][1]
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

    Domain based on Grid Edges
        UL_edge(+70.125, -20.125)------------------UR_edge(+70.125, +40.125)
        |                                                                  |
        LL_edge(+29.875, -20.125)------------------LR_edge(+29.875, +40.125)

    Domain based on Grid Centers
        UL(+70.000, -20.000)------------------UR(+70.000, +40.000)
        |                                                        |
        LL(+30.000, -20.000)------------------LR(+30.000, +40.000)

    The MSG CloudCast Domain is (ccast_crs)
        UL(+62.403, -12.921)------------------UR(+60.151, +33.739)
        |
        LL(+42.068,  -4.803)------------------LR(+40.928, +21.329)

    So the MSG domain fits easily within the bounds of the ERA-5 source data.



    """
    ncfile = netCDF4.Dataset(a_file, mode='r', format='NETCDF4_CLASSIC')

    era5_lats = ncfile.variables["latitude"][:]
    era5_lons = ncfile.variables["longitude"][:]
    nlats = len(era5_lats)
    nlons = len(era5_lons)
    dlat = dlon = 0.25
    ddlat = ddlon = 0.25 * 0.5

    era5_lat_edges = dict([(round(float(i) + ddlat, 3), 1) for i in era5_lats])
    era5_lat_edges.update([(round(float(i) - ddlat, 3), 1) for i in era5_lats])
    era5_lat_edges = sorted(list(era5_lat_edges.keys()), reverse=True)

    era5_lon_edges = dict([(float(i) + ddlon, 1) for i in era5_lons])
    era5_lon_edges.update([(float(i) - ddlon, 1) for i in era5_lons])
    era5_lon_edges = sorted(list(era5_lon_edges.keys()))

    # Corners of bounding box
    ll_x_edge = era5_lon_edges[0]
    ll_y_edge = era5_lat_edges[-1]

    ul_x_edge = era5_lon_edges[0]
    ul_y_edge = era5_lat_edges[0]

    ur_x_edge = era5_lon_edges[-1]
    ur_y_edge = era5_lat_edges[0]

    lr_x_edge = era5_lon_edges[-1]
    lr_y_edge = era5_lat_edges[-1]

    ll_x = era5_lons[0]
    ll_y = era5_lats[-1]
    ul_x = era5_lons[0]
    ul_y = era5_lats[0]
    ur_x = era5_lons[-1]
    ur_y = era5_lats[0]
    lr_x = era5_lons[-1]
    lr_y = era5_lats[-1]

    if verbose:
        print(f"\tERA-5 Latitudes                           ({nlats}): [{era5_lats[0]:+8.3f}, ... {era5_lats[-1]:+8.3f}]")
        print(f"\tERA-5 Longitudes                          ({len(era5_lons)}): [{era5_lons[0]:+8.3f}, ... {era5_lons[-1]:+8.3f}]")
        print("\n")
        print(f"\tUL_edge({ul_y_edge:+7.3f}, {ul_x_edge:+7.3f})------------------UR_edge({ur_y_edge:+7.3f}, {ur_x_edge:+7.3f})")
        print("\t|                                                                  |")
        print(f"\tLL_edge({ll_y_edge:+7.3f}, {ll_x_edge:+7.3f})------------------LR_edge({lr_y_edge:+7.3f}, {lr_x_edge:+7.3f})")
        print("\n")
        print(f"\tUL({ul_y:+7.3f}, {ul_x:+7.3f})------------------UR({ur_y:+7.3f}, {ur_x:+7.3f})")
        print("\t|                                                        |")
        print(f"\tLL({ll_y:+7.3f}, {ll_x:+7.3f})------------------LR({lr_y:+7.3f}, {lr_x:+7.3f})")

    ##
    # Read MSG lon and lats from file (saved in natread.py)
    lonlatfile = "/Users/mbauer/tmp/CloudCast/ccast_lonlat.npy"
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

    from scipy.interpolate import griddata

    ##
    # The ERA-5 mesh to interpolate from (38801, 2)
    grid_lons, grid_lats = np.meshgrid(era5_lons, era5_lats)
    sparse_points = np.stack([grid_lons.ravel(), grid_lats.ravel()], -1)  # shape (N, 2) in 2d

    ##
    # The MSG fine mesh to interpolate into (589824, 2)
    # fine_points = np.stack([ccast_lons.ravel(), ccast_lats.ravel()], -1)  # shape (N, 2) in 2d
    ccast_nlons = len(ccast_lons)
    ccast_nlats = len(ccast_lats)

    ##
    # A simple ERA-5 Temperature Field
    era5_temp = np.squeeze(ncfile.variables["t"][:], axis=1)
    grid_centers_z_vals = era5_temp[0, :, :]
    del era5_temp
    print(f"\tERA-5 Temperature {grid_centers_z_vals.shape}: min {np.amin(grid_centers_z_vals):8.3f}, mean {np.mean(grid_centers_z_vals):8.3f}, max {np.amax(grid_centers_z_vals):8.3f} K")

    # z_dense_smooth_griddata = interp.griddata(sparse_points, z_sparse_smooth.ravel(), (x_dense, y_dense), method='cubic')

    use_method = ('nearest', 'linear', 'cubic')[2]
    fine_mesh_vals = griddata(sparse_points, grid_centers_z_vals.ravel(), (ccast_lons.ravel(), ccast_lats.ravel()), method=use_method)
    # Reshape
    fine_mesh_vals = np.reshape(fine_mesh_vals, (ccast_nlats, ccast_nlons))
    print(f"\tInterpolated ERA-5 Temperature {fine_mesh_vals.shape}: min {np.amin(fine_mesh_vals):8.3f}, mean {np.mean(fine_mesh_vals):8.3f}, max {np.amax(fine_mesh_vals):8.3f} K")
    # print(np.any(np.isnan(fine_mesh_vals)))


    ##
    # Plot
    import cartopy.feature as cfeature
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import matplotlib.colors as mcolors
    from cloudcast.plot.whitegreenblueyellow_colors import read_wgby_colors_from_json
    from cloudcast.plot.fix_lon_lat_for_pcolormesh import fix_lon_lat_for_pcolormesh

    wide_domain = [False, True][0]
    as_lcc = [False, True][0]

    the_title = f"ERA-5 {do_lev[1:]} hPa as MSG via {use_method.capitalize()} Interpolation"

    if not 'whitegreenblueyellow_colors' in mpl.colormaps:
        whitegreenblueyellow_colors_data = read_wgby_colors_from_json()
        whitegreenblueyellow_colors = mcolors.LinearSegmentedColormap('whitegreenblueyellow_colors', segmentdata=whitegreenblueyellow_colors_data, N=256)
        mpl.colormaps.register(name='whitegreenblueyellow_colors', cmap=whitegreenblueyellow_colors)
        color_scheme = mpl.colormaps.get_cmap('whitegreenblueyellow_colors')
    else:
        color_scheme = mpl.colormaps.get_cmap('whitegreenblueyellow_colors')

    freq_map_cmap = ["plasma", "gist_rainbow", "jet", "brg", 'whitegreenblueyellow_colors'][4]

    make_mesh = [False, True][0]
    if make_mesh:
        plon_mid, plat_mid = fix_lon_lat_for_pcolormesh(np.copy(ccast_lons), np.copy(ccast_lats))
        x_lon, y_lat = np.meshgrid(plon_mid, plat_mid)

    pname = f"era5_2_msg_{do_lev}_{use_method}_{freq_map_cmap}.png"
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

    #   ll_lon: (-4.816534709314307, 42.053336570266744)
    #   ur_lon: (33.77742545811928, 60.15622631725923)
    # use_extent = ccast_crs.bounds
    # ur_latlon = geod_crs.transform_point(use_extent[1], use_extent[3], ccast_crs)
    # ll_latlon = geod_crs.transform_point(use_extent[0], use_extent[2], ccast_crs)
    # print(f"\tll_lon: ({float(ll_latlon[0])}, {float(ll_latlon[1])})")
    # print(f"\tur_lon: ({float(ur_latlon[0])}, {float(ur_latlon[1])})")

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
        print(f"\tll_lon: ({float(ll_lon[0])}, {float(ll_lon[1])})")
        print(f"\tur_lon: ({float(ur_lon[0])}, {float(ur_lon[1])})")
        use_crs = ccast_lcc_crs
    else:
        #   use_extent = (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)
        use_extent = ccast_crs.bounds
        use_crs = ccast_crs
    print(f"\n\tuse_extent = ({float(use_extent[0])}, {float(use_extent[1])}, {float(use_extent[2])}, {float(use_extent[3])})")


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

    the_min = np.round(np.amin(fine_mesh_vals), decimals=0) + range_offset
    the_max = np.round(np.amax(fine_mesh_vals), decimals=0) - range_offset

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
        aplt =  geo_axes.pcolormesh(x_lon, y_lat, ma.masked_equal(fine_mesh_vals, 0), transform=use_crs, zorder=1, **param_dict)
    else:
        # geo_axes.contourf(lons, lats_nh, ma.masked_equal(map_this, 0), 1, transform=ccrs.PlateCarree())
        i_ways = ['none', 'antialiased', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36',
                  'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel',
                  'mitchell', 'sinc', 'lanczos', 'blackman']
        iway = i_ways[1]
        # param_dict = {"extent": map_extent, "interpolation": iway, "cmap": color_scheme, "origin": 'lower', "norm": norm}
        # aplt = geo_axes.imshow(ma.masked_equal(map_this, 0), transform=map_crs, **param_dict)
        if as_lcc:
            a_image = plt.imshow(fine_mesh_vals, interpolation=iway, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=use_alpha)
        else:
            a_image = plt.imshow(fine_mesh_vals, interpolation=iway, cmap=freq_map_cmap, transform=use_crs, extent=use_extent, origin='upper', vmin=the_min, vmax=the_max, alpha=use_alpha)

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

    ##
    # Save to NETCDF
    savefile = f"/Users/mbauer/tmp/CloudCast/test_era_2_msg_{do_year}_{do_lev}_{use_method}.nc"
    ncfile = netCDF4.Dataset(savefile, mode='w', format='NETCDF4_CLASSIC')
    lat_dim = ncfile.createDimension('ny', ccast_nlats) # latitude axis
    lon_dim = ncfile.createDimension('nx', ccast_nlons) # longitude axis

    # Note different that 1d version where create dim 'lat' and a 1d variable named 'lat'
    #   here we use dim name 'ny' and var name 'lat'
    r"""
        dimensions:
         x = 356;
         y = 236;
         height_2m = 1;
         time = UNLIMITED; // (365 currently
         nb2 = 2;

    variables:
         float lon(y=236, x=356);
             :standard_name = "longitude";
             :long_name = "longitude";
             :units = "degrees_east";
             :_CoordinateAxisType = "Lon";

         float lat(y=236, x=356);
             :standard_name = "latitude";
             :long_name = "latitude";
             :units = "degrees_north";
             :_CoordinateAxisType = "Lat";

         float height_2m(height_2m=1);
             :standard_name = "height";
             :long_name = "height above the surface";
             :units = "m";
             :positive = "up";
             :axis = "Z";

         double time(time=365);
             :standard_name = "time";
             :long_name = "time";
             :bounds = "time_bnds";
             :units = "seconds since 1979-01-01 00:00:00";
             :calendar = "proleptic_gregorian";

         double time_bnds(time=365, nb2=2);
             :units = "seconds since 1979-01-01 00:00:00";
             :calendar = "proleptic_gregorian";

         float TMAX_2M(time=365, height_2m=1, y=236, x=356);
             :standard_name = "air_temperature";
             :long_name = "2m maximum temperature";
             :units = "K";
             :coordinates = "lon lat";
             :cell_methods = "time: maximum";
             :_FillValue = -9.0E33f; // float
    """


    lat = ncfile.createVariable('lat', np.float32, ('ny', 'nx'), least_significant_digit=4, fill_value=-999)
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lat.axis = "Y"
    lat.CoordinateAxisType = 'Lat'

    lon = ncfile.createVariable('lon', np.float32, ('ny', 'nx'), least_significant_digit=4, fill_value=-999)
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    lon.axis = "X"
    lon.CoordinateAxisType = 'Lon'

    outdata = ncfile.createVariable('t', np.float32, ('ny', 'nx'), least_significant_digit=3, fill_value=0)
    outdata.units = 'K'
    outdata.long_name = f"ERA-5 {do_year} {do_lev} hPa Temperature"
    outdata.setncattr('coordinates', 'lat lon')

    lat[:] = ccast_lats
    lon[:] = ccast_lons
    outdata[:] = fine_mesh_vals

    # nc_crs = ncfile.createVariable('spatial_ref', 'i4')
    # nc_crs.spatial_ref='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

    ncfile.close()


    # ##
    # # Create some information on the reference system
    # CCAST_HEIGHT = 768
    # CCAST_WIDTH = 768
    # lower_left_xy = [-855100.436345, -4942000.0]
    # upper_right_xy = [1448899.563655, -2638000.0]

    # # Define the area
    # #   <class 'pyresample.geometry.AreaDefinition'>
    # ccast_area_def = pr.geometry.AreaDefinition('areaD', 'Europe', 'areaD',
    #                                             {'lat_0': '90.00', 'lat_ts': '50.00', 'lon_0': '5', 'proj': 'stere', 'ellps': 'WGS84'},
    #                                             CCAST_HEIGHT, CCAST_WIDTH,
    #                                             (lower_left_xy[0], lower_left_xy[1], upper_right_xy[0], upper_right_xy[1]))
    # #   +proj=stere +lat_0=90 +lat_ts=50 +lon_0=5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs +type=crs
    # # print(ccast_area_def.proj4_string)

    # ##
    # # Form a cartopy CRS
    # #   <class 'pyresample.utils.cartopy.Projection'>
    # ccast_crs = ccast_area_def.to_cartopy_crs()
    # # print(ccast_merc_crs.bounds)
    # #   ccast_crs.bounds: (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)


    # dtype = "float32"
    # radius = 16000
    # epsilon = 0.5
    # nodata = -3.4E+38
    # ##
    # # Apply a swath definition for our output raster
    # #   <class 'pyresample.geometry.SwathDefinition'>
    # exstract_def = pr.geometry.SwathDefinition(lons=era5_lons, lats=era5_lats)

    # ##
    # # Resample our data to the area of interest
    # ccast_data_vals = pr.kd_tree.resample_nearest(exstract_def, data_vals,
    #                                               ccast_area_def,
    #                                               radius_of_influence=radius, # in meters
    #                                               epsilon=epsilon,
    #                                               fill_value=False)
    # ccast_lons = pr.kd_tree.resample_nearest(exstract_def, lons,
    #                                          ccast_area_def,
    #                                          radius_of_influence=radius, # in meters
    #                                          epsilon=epsilon,
    #                                          fill_value=False)
    # ccast_lats = pr.kd_tree.resample_nearest(exstract_def, lats,
    #                                          ccast_area_def,
    #                                          radius_of_influence=radius, # in meters
    #                                          epsilon=epsilon,
    #                                          fill_value=False)


# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
