#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
find_grid_area
~~~~~~~~~~~~~~
"""
# Standard Imports
import math

# Third-Party Imports
import numpy.typing as npt

# STARE Imports

# Local Imports

##
# List of Public objects from this module.
__all__ = ['find_grid_area']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

###############################################################################
# PUBLIC find_grid_area()
# ---------------------
def find_grid_area(lons: npt.ArrayLike, lats: npt.ArrayLike) -> list[float]:
    """_Calculate data-grid (Latitude-Longitude Quadrangle) surface area as a function of latitude.

    Parameters
    ----------
    lons : npt.ArrayLike
        Longitude values in degrees
    lats : npt.ArrayLike
        Latitude values in degrees

    Returns
    -------
    list[float]
        Surface area of the data-grid (Latitude-Longitude Quadrangle) as a function of latitude [km^2].
    """
    earth_radius_km = 6371.2
    earth_radius_kmsq = earth_radius_km * earth_radius_km
    grid_area = []
    nlats = len(lats)
    lon_180 = lambda x: ((x + 180.0) % 360.0) - 180.0
    """
    For MERRA2, polar_method = "C" at Equator

                                       >            dlon = 0.625 deg or 0.01091 rad           <

    (-0.3125, +0.75)+------------------+------------------+(+0.3125, +0.75)+------------------+------------------+(+0.9375, +0.75)
             |                                                     |                                                      |
             +                  (0.0000, +0.500)                   +                   (+0.625, +0.500)                   +              |
             |                                                     |                                                      |
    (-0.3125, +0.25)+------------------+------------------+(+0.3125, +0.25)+------------------+------------------+(+0.9375, +0.25)      dlat = 0.5 deg
             |                                                     |                                                      |
             +                  (0.0000,  0.000)                   +                   (+0.625,  0.000)                   +              |
             |                                                     |                                                      |
    (-0.3125, -0.25)+------------------+------------------+(+0.3125, -0.25)+------------------+------------------+(+0.9375, -0.25)

    Surface Area of Latitude-Longitude Quadrangle:
        * Latitudes and Longitudes in radians
            - lat_1 > lat_0
            - lon_1 > lon_0, longitudes must be in +/-180 format (not 0/360)
            - R is the mean radius of the Earth

        (sin(lat_1) - sin(lat_0) * (lon_1 - lon_0) * R^2

        Check by using global size grid:
            lat_1 = radians(90.0)   =
            lat_0 = radians(-90.0)  =
            lon_1 = radians(180.0)  =
            lon_0 = radians(-180.0) =

        Total Surface Area: 510,096,496 km^2
        Should be close to  510,072,000 km^2 using mean radius of 6,371 km

    From diagnostic output below:
        180 | -0.250 --  +0.000 --  +0.250|  3,864.0 km^2

    From: https://www.engr.scu.edu/~emaurer/tools/calc_cell_area_cgi.pl

        Area = 3863.828 km2

        Grid Cell Box Dimensions:
            Top Length ( 0.25 deg lat)   : 69.573 km
            Bottom Length (-0.25 deg lat): 69.573 km
            Sides (centered 0.00 deg lat): 55.659 km

        A = L * W
          = 69.573 * 55.659
          = 3872.362 km^2
    """
    # (Comment out normally)
    # lat_1 = math.radians(90.0)
    # lat_0 = math.radians(-90.0)
    # lon_1 = math.radians(180.0)
    # lon_0 = math.radians(-180.0)
    # lat_delta = abs(math.sin(lat_1) - math.sin(lat_0))
    # lon_delta = abs(lon_1 - lon_0)
    # tarea = lat_delta * lon_delta * earth_radius_kmsq
    # print(f"Total Surface Area: {int(tarea):,} km^2")
    # raise Exception("Stop Here")

    # Longitude span (radians) of the Quadrangle (Lons must in +/-180 format)
    dlon = abs(math.radians(lon_180(lons[1])) - math.radians(lon_180(lons[0])))

    # Latitude half-span (degrees) of the Quadrangle
    dlat = 0.5 * abs(lats[10] - lats[11])

    # Check whole Earth (Comment out normally)
    """
        Total Surface Area: 510,096,496 km^2
        Should be close to  510,072,000 km^2 using mean radius of 6,371 km
    """
    # whole_earth = 1
    # dlon = math.radians(360.0)
    # dlat = 180.0

    multiplier = earth_radius_kmsq * dlon
    # print(f"dlon {dlon:.4G} dlat {dlat:.4G} multiplier {multiplier:.3f}")
    for jj in range(nlats):
        lata = lats[jj] - dlat
        latb = lats[jj] + dlat
        if lata < -90.0:
            lata = -90.0
        if latb > 90.0:
            latb = 90.0
        garea = multiplier * abs(math.sin(math.radians(latb)) - math.sin(math.radians(lata)))
        grid_area.append(garea)

        # if whole_earth:
        #     print(f"Total Surface Area: {int(garea):,} km^2")
        #     raise Exception("Stop Here")
        #     break
        # print(f"{jj:03d} |{lata:+7.3f} -- {lats[jj]:+7.3f} -- {latb:+7.3f}| {int(garea):7d}")

    # total_area = sum(grid_area[_] * len(lons) for _ in range(nlats))
    # print(f"Total Surface Area: {int(total_area):,} km^2")
    """
    Total Surface Area: 510,096,496 km^2
    Should be close to  510,072,000 km^2 using mean radius of 6,371 km
    """
    return grid_area

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
