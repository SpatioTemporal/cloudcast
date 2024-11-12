#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Setup geometry for Meteosat Second Generation (MSG).

setup_msg.py
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

# STARE Imports

# Local Imports

##
# List of Public objects from this module.
__all__ = ['setup_msg']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

def LON_TO_180(x): return ((x + 180.0) % 360.0) - 180.0
def LON_TO_360(x): return (x + 360.0) % 360.0

##
# Type Aliases
SetOut: TypeAlias = tuple[int, int, int, typing.TypeVar('cartopy.crs'),
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')],
                     Union[None, pr.geometry.AreaDefinition], Union[None, typing.TypeVar('cartopy.crs')]]

###############################################################################
# PUBLIC setup_msg()
# ------------------
def setup_msg(use_var: str, to_ccast:bool, to_euro:bool, to_merc:bool, to_lcc:bool) -> SetOut:
    """
    Setup geometry for Meteosat Second Generation (MSG).

    Parameters
    ----------
    use_var : str
        MSG variable/channel to use
    to_ccast : bool
        Return geometry for CloudCast
    to_euro : bool
        Return geometry for Europe subset
    to_merc : bool
        Return geometry for Mercator projection
    to_lcc : bool
        Return geometry for Lambert Conformal Conic projection

    Returns
    -------
    SetOut
        Tuple of geometry information for the specified projection and variable.
    """
    verbose = [False, True][1]

    euro_nrows = 928
    euro_ncols = 1530

    # Type of Mercator to use
    use_pseudo = [False, True][0]

    ##
    # MSG base projection WGS84 - World Geodetic System 1984
    #   https://epsg.io/4326
    #   +proj=longlat +datum=WGS84 +no_defs +type=crs
    #MSG_EPSG = 4326
    # msg_crs = ccrs.Geodetic()
    MSG_EPSG = 32662
    msg_crs = ccrs.PlateCarree()

    ##
    # Mercator EPSG
    MERC_EPSG = 3857 if use_pseudo else 3395

    ##
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

    # Default values
    ccast_area_def = None
    ccast_crs = None
    ccast_merc_area_def = None
    ccast_merc_crs = None
    ccast_lcc_area_def = None
    ccast_lcc_crs = None

    msg_area_def = None
    msg_merc_area_def = None
    msg_merc_crs = None

    if to_ccast:
        ##
        # Create some information on the reference system
        CCAST_HEIGHT = 768
        CCAST_WIDTH = 768
        lower_left_xy = [-855100.436345, -4942000.0]
        upper_right_xy = [1448899.563655, -2638000.0]

        # From cloudcast repo latlonraster() in cloudcast/base/plotutils.py
        # img_size = (768, 768)
        # x = np.linspace(-1065644.490, 1306855.478, 768)
        # y = np.linspace(9683729.573, 7011229.349, 768)

        # Define the area
        #   <class 'pyresample.geometry.AreaDefinition'>
        ccast_area_def = pr.geometry.AreaDefinition('areaD', 'Europe', 'areaD',
                                                    {'lat_0': '90.00', 'lat_ts': '50.00', 'lon_0': '5', 'proj': 'stere', 'ellps': 'WGS84'},
                                                    CCAST_HEIGHT, CCAST_WIDTH,
                                                    (lower_left_xy[0], lower_left_xy[1], upper_right_xy[0], upper_right_xy[1]))
        #   +proj=stere +lat_0=90 +lat_ts=50 +lon_0=5 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs +type=crs
        # print(ccast_area_def.proj4_string)

        ##
        # Form a cartopy CRS
        #   <class 'pyresample.utils.cartopy.Projection'>
        ccast_crs = ccast_area_def.to_cartopy_crs()
        # print(ccast_merc_crs.bounds)
        #   ccast_crs.bounds: (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)

        """
        ccast_crs
            ccast_area_def.corners = [
                (-12.920934886492649, 62.403066090517555), UL
                (33.73865749382469,   60.15059617915855),  UR
                (21.32880156090482,   40.92817004367345),  LR
                (-4.802566482888071,  42.068097533886025)] LL

            ccast_area_def.area_extent  (-855100.436345, -4942000.0, 1448899.563655, -2638000.0)
            ccast_crs.bounds      (-855100.436345, 1448899.563655, -4942000.0, -2638000.0)

            lower_left_xy  = [-855100.436345, -4942000.0] => (-4.816534709314307, 42.053336570266744)
            upper_right_xy = [1448899.563655, -2638000.0] => (33.77742545811928,  60.15622631725923)
        """
        # print(ccast_area_def.corners)
        # print(ccast_area_def.area_extent)
        # print(ccast_crs.bounds)
        # ur_xy = geod_crs.transform_point(upper_right_xy[0], upper_right_xy[1], ccast_crs)
        # ll_xy = geod_crs.transform_point(lower_left_xy[0], lower_left_xy[1], ccast_crs)
        # print(f"ll_xy: {ll_xy}")
        # print(f"ur_xy: {ur_xy}")
        if to_merc:
            ##
            # Form a Mercator CSR version of ccast_crs
            ccast_merc_crs = ccrs.epsg(MERC_EPSG)
            if use_pseudo:
                ##
                # Mercator version EPSG:3857 for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI
                #   https://epsg.io/3857
                area_id = "EPSG:3857"
                description = "Pseudo-Mercator"
                proj_id = "EPSG:3857"
            else:
                ##
                # World Mercator version EPSG:3395 for WGS 84
                #   https://epsg.io/3395
                area_id = "EPSG:3395"
                description = "Mercator"
                proj_id = "EPSG:3395"

            #   ll_xy: (-536174.1912289965, 5140343.824785849)
            #   ur_xy: (3760085.802305594, 8397504.685448818)
            ur_xy = ccast_merc_crs.transform_point(upper_right_xy[0], upper_right_xy[1], ccast_crs)
            ll_xy = ccast_merc_crs.transform_point(lower_left_xy[0], lower_left_xy[1], ccast_crs)
            # print(f"ll_xy: {ll_xy}")
            # print(f"ur_xy: {ur_xy}")

            ##
            # Specify projection parameters
            if use_pseudo:
                projection = {'proj': 'merc', 'lat_ts': 0, 'lon_0': 0, 'a': 6378137, 'b': 6378137, 'x_0': 0, 'y_0': 0, 'k': 1, 'units': 'm'}
            else:
                projection = {'proj': 'merc','x_0': 0, 'y_0': 0, 'k': 1, 'units': 'm', 'no_defs': None, 'datum': 'WGS84', 'type': 'crs'}

            ##
            # Define the area
            ccast_merc_area_def = pr.geometry.AreaDefinition(area_id, 'Europe', proj_id,
                                                             projection, CCAST_HEIGHT, CCAST_WIDTH,
                                                             (ll_xy[0], ll_xy[1], ur_xy[0], ur_xy[1]))
            #   proj4_string: +proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs
            # print(f"proj4_string: {ccast_merc_area_def.proj4_string}")

            ccast_merc_crs = ccast_merc_area_def.to_cartopy_crs()

            """
            ccast_merc_crs
                merc_area_corners = [(-4.7914084331636335, 60.14672953639852),  UL
                                     (33.75229918196862, 60.14672953639852),    UR
                                     (33.75229918196862, 42.06753197287029),    LR
                                     (-4.7914084331636335, 42.06753197287029)]  LL

                ccast_merc_area_def.area_extent (-536174.1912289965, 5140343.824785849, 3760085.802305594, 8397504.685448818)
                ccast_merc_crs.bounds           (-536174.1912289965, 5140343.824785849, 3760085.802305594, 8397504.685448818)

                ll_xy  = (-536174.1912289965, 5140343.824785849) => (-4.816534709314306, 42.05333657026676)
                ur_xy  = (3760085.802305594, 8397504.685448818)  => (33.77742545811928, 60.15622631725923)
            """
            # print(ccast_merc_area_def.corners)
            # print(ccast_merc_area_def.area_extent)
            # print(ccast_merc_crs.bounds)
            # ur_xya = geod_crs.transform_point(ur_xy[0], ur_xy[1], ccast_merc_crs)
            # ll_xya = geod_crs.transform_point(ll_xy[0], ll_xy[1], ccast_merc_crs)
            # print(f"ll_xy: {ll_xya}")
            # print(f"ur_xy: {ur_xya}")
            # os._exit(1)
        if to_lcc:
            ##
            # Form a Lambert Conformal CSR version of ccast_crs
            ccast_lcc_crs = ccrs.LambertConformal(central_longitude=10.0, central_latitude=51.0, cutoff=-30)
            # print((-12.920934886492649 + 33.73865749382469) * 0.5)
            # print((40.92817004367345 + 62.403066090517555) * 0.5)

            ##
            # Lambert Conic Conformal (1SP) EPSG:9801
            #   https://epsg.io/9801-method
            area_id = "EPSG:9801"
            description = "Lambert"
            proj_id = "EPSG:9801"

            #   ll_xy: (-1216087.7976068966, -899612.0368218884)
            #   ur_xy: (1404456.8408720866,  1248209.190564286)
            ur_xy = ccast_lcc_crs.transform_point(upper_right_xy[0], upper_right_xy[1], ccast_crs)
            ll_xy = ccast_lcc_crs.transform_point(lower_left_xy[0], lower_left_xy[1], ccast_crs)
            # print(f"ll_xy: {ll_xy}")
            # print(f"ur_xy: {ur_xy}")

            # #   ll_xy: (-4.816534709314301, 42.053336570266744)
            # #   ur_xy: (33.77742545811928,  60.15622631725923)
            # ur_xy = geod_crs.transform_point(ur_xy[0], ur_xy[1], ccast_lcc_crs)
            # ll_xy = geod_crs.transform_point(ll_xy[0], ll_xy[1], ccast_lcc_crs)
            # print(f"ll_xy: {ll_xy}")
            # print(f"ur_xy: {ur_xy}")

            ##
            # Specify projection parameters
            projection = {'proj': 'lcc', 'lon_0': 10, 'lat_0': 51, 'lat_1': 40, 'lat_2': 60, 'type': 'crs'}

            ##
            # Define the area
            ccast_lcc_area_def = pr.geometry.AreaDefinition(area_id, 'Europe', proj_id,
                                                             projection, CCAST_HEIGHT, CCAST_WIDTH,
                                                             (ll_xy[0], ll_xy[1], ur_xy[0], ur_xy[1]))
            #   proj4_string: +proj=lcc +lat_0=0 +lon_0=10 +lat_1=51 +lat_2=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs
            # print(f"proj4_string: {ccast_lcc_area_def.proj4_string}")

            ccast_lcc_crs = ccast_lcc_area_def.to_cartopy_crs()

            """
            ccast_lcc_crs
                lcc_area_corners = [(-12.47631979109102, 60.62706437795364),  UL
                                    (35.7080370218538, 60.08967574868918),    UR
                                    (26.973396040999436, 41.361735203234744), LR
                                    (-4.757808894947807, 41.72218823476732)]  LL

                ccast_lcc_area_def.area_extent (-1216087.7976068966, -899612.0368218884, 1404456.8408720866, 1248209.190564286)
                ccast_lcc_crs.bounds           (-1216087.7976068966, 1404456.8408720866, -899612.0368218884, 1248209.190564286)

                ll_xy  = (-1216087.7976068966, -899612.0368218884) =>  (-4.774671331638412, 41.70674664619926)
                ur_xy  = (1404456.8408720866, 1248209.190564286)   =>  (35.745388874313626, 60.0962915843566)
            # """
            # print(ccast_lcc_area_def.corners)
            # print(ccast_lcc_area_def.area_extent)
            # print(ccast_lcc_crs.bounds)
            # ur_xya = geod_crs.transform_point(ur_xy[0], ur_xy[1], ccast_lcc_crs)
            # ll_xya = geod_crs.transform_point(ll_xy[0], ll_xy[1], ccast_lcc_crs)
            # print(f"ll_xy: {ll_xya}")
            # print(f"ur_xy: {ur_xya}")
            # os._exit(1)
    elif to_euro:
        print("Fix")
        os._exit(1)
    else:
        ##
        # Geostationary Projection (GEOS) EPSG
        #   https://proj4.org/en/9.2/operations/projections/geos.html
        #       +proj=geos +h=42164000.0 +R=6378000.0 +lon_0=0 +sweep=y
        #
        #   https://scitools.org.uk/cartopy/docs/latest/reference/projections.html#geostationary
        msg_crs = ccrs.Geostationary(central_longitude=0.0, satellite_height=42164000, globe=globe, sweep_axis='y')
        # msg_crs = ccrs.Geostationary(central_longitude=0.0)

        # print(msg_crs.boundary)
        # LINEARRING (5434181.528949673 0, 5421058.433025361 -377615.14583773003, 5381748.153532824 -753396.287082276, 5316427.634887653 -1125518.1235829785, 5225391.5086746905 -1492172.7295488187, 5109051.653353975 -1851578.1544986803, 4967936.525591165 -2201986.9209442344, 4802690.213302837 -2541694.3846535054, 4614071.149595708 -2869046.9234366952, 4402950.418356084 -3182449.9203697643, 4170309.5768173663 -3480375.5071729054, 3917237.9184494824 -3761370.0330624497, 3644929.1013382724 -4024061.223803627, 3354677.073098844 -4267164.99495919, 3047871.233382138 -4489491.882551729, 2725990.7890863204 -4689953.053674234, 2390598.275157427 -4867565.859165926, 2043332.234811273 -5021458.8905083425, 1685899.0763594205 -5150876.503779991, 1320064.1486001024 -5255182.775005638, 947642.101802744 -5333864.853676236, 570486.6254336371 -5386535.68466741, 190479.67567897929 -5412936.073245876, -190479.67567897862 -5412936.073245876, -570486.6254336376 -5386535.68466741, -947642.1018027435 -5333864.853676236, -1320064.148600103 -5255182.775005638, -1685899.07635942 -5150876.503779991, -2043332.2348112732 -5021458.890508342, -2390598.275157427 -4867565.859165926, -2725990.7890863186 -4689953.053674235, -3047871.2333821375 -4489491.882551729, -3354677.0730988444 -4267164.9949591905, -3644929.1013382724 -4024061.223803627, -3917237.918449484 -3761370.0330624497, -4170309.576817365 -3480375.5071729063, -4402950.418356084 -3182449.920369765, -4614071.149595708 -2869046.9234366957, -4802690.213302837 -2541694.384653504, -4967936.525591164 -2201986.9209442353, -5109051.653353973 -1851578.1544986812, -5225391.5086746905 -1492172.729548819, -5316427.634887653 -1125518.1235829785, -5381748.153532824 -753396.2870822754, -5421058.433025361 -377615.1458377312, -5434181.528949673 -6.629406305492799e-10, -5421058.433025361 377615.1458377298, -5381748.153532824 753396.2870822764, -5316427.6348876525 1125518.1235829794, -5225391.5086746905 1492172.7295488177, -5109051.653353974 1851578.15449868, -4967936.525591165 2201986.9209442344, -4802690.213302836 2541694.384653505, -4614071.149595709 2869046.923436695, -4402950.418356084 3182449.9203697634, -4170309.5768173663 3480375.5071729054, -3917237.918449482 3761370.03306245, -3644929.101338272 4024061.2238036273, -3354677.0730988425 4267164.9949591905, -3047871.233382137 4489491.88255173, -2725990.7890863223 4689953.053674233, -2390598.2751574274 4867565.8591659255, -2043332.2348112746 5021458.890508342, -1685899.0763594212 5150876.503779991, -1320064.148600103 5255182.775005638, -947642.1018027436 5333864.853676236, -570486.6254336365 5386535.684667411, -190479.67567897754 5412936.073245876, 190479.6756789804 5412936.073245876, 570486.6254336345 5386535.684667411, 947642.1018027416 5333864.853676236, 1320064.148600101 5255182.775005638, 1685899.0763594192 5150876.503779991, 2043332.234811273 5021458.8905083425, 2390598.2751574265 4867565.859165928, 2725990.7890863204 4689953.053674234, 3047871.233382139 4489491.882551729, 3354677.073098846 4267164.99495919, 3644929.1013382706 4024061.2238036282, 3917237.918449481 3761370.033062451, 4170309.5768173644 3480375.507172907, 4402950.418356084 3182449.9203697653, 4614071.149595707 2869046.9234366952, 4802690.213302837 2541694.3846535054, 4967936.525591165 2201986.920944234, 5109051.653353974 1851578.1544986798, 5225391.5086746905 1492172.7295488175, 5316427.634887653 1125518.1235829815, 5381748.1535328245 753396.2870822785, 5421058.433025361 377615.14583773183, 5434181.528949673 1.3258812610985598e-9, 5434181.528949673 0)
        # tmp = msg_crs.boundary
        # tmpr = list(tmp.coords)
        # tmpr_x = [_[0] for _ in tmpr]
        # tmpr_y = [_[1] for _ in tmpr]
        # print(min(tmpr_x), max(tmpr_x))
        # print(min(tmpr_y), max(tmpr_y))
        # -5434177.815885395 5434177.815885395
        # -5412932.376718026 5412932.376718026

        # print(msg_crs.proj4_params)
        # {'datum': 'WGS84', 'ellps': 'WGS84', 'proj': 'geos', 'lon_0': 0.0, 'lat_0': 0.0, 'h': 35786000, 'x_0': 0, 'y_0': 0, 'units': 'm', 'sweep': 'y'}

        # print(msg_crs.to_proj4())
        # +proj=geos +datum=WGS84 +ellps=WGS84 +lon_0=0.0 +lat_0=0.0 +h=35786000 +x_0=0 +y_0=0 +units=m +sweep=y +no_defs +type=crs

        #   (-5434181.528949673, 5434181.528949673)
        #   (-5412936.073245876, 5412936.073245876)
        # print(msg_crs.x_limits)
        # print(msg_crs.y_limits)

        if use_var == "HRV":
            MSG_HEIGHT = 11136
            MSG_WIDTH = 5568
            # lower_left_xy = [-57.18206024169922, -78.96027374267578]
            # upper_right_xy = [79.4859390258789, 78.3008804321289]
            lower_left_xy = [-79.4859390258789, -78.96027374267578]
            upper_right_xy = [79.4859390258789, 78.3008804321289]
        else:
            if to_euro:
                MSG_HEIGHT = euro_nrows
                # MSG_WIDTH = euro_ncols
                MSG_WIDTH = 3712
                lower_left_xy = [-75.26545715332031, 26.65]
                upper_right_xy = [75.56462097167969, 78.29975891113281]
            else:
                MSG_HEIGHT = 3712
                MSG_WIDTH = 3712
                lower_left_xy = [-75.26545715332031, -78.95874786376953]
                upper_right_xy = [75.56462097167969, 78.29975891113281]

        # area_id = f"EPSG:{MSG_EPSG}"
        area_id = f"GEOS"
        description = "Partial Disk" if use_var == "HRV" else  "Full Disk"
        proj_id = f"GEOS"

        ##
        # Define the area
        #   <class 'pyresample.geometry.AreaDefinition'>
        # msg_area_def = pr.geometry.AreaDefinition(area_id, description, proj_id,
        #                                             {'lat_0': '0.00', 'lat_ts': '0.00', 'lon_0': '0.00', 'proj': 'longlat', 'ellps': 'WGS84'},
        #                                             MSG_HEIGHT, MSG_WIDTH,
        #                                             (lower_left_xy[0], lower_left_xy[1], upper_right_xy[0], upper_right_xy[1]))
        # #   +proj=longlat +lat_0=0.0 +lon_0=0.0 +ellps=WGS84 +type=crs
        # # print(msg_area_def.proj4_string)

        # msg_area_def = pr.geometry.AreaDefinition(area_id, description, proj_id,
        #                                           {'lat_0': '0.00', 'lat_ts': '0.00', 'lon_0': '0.00', 'proj': 'eqc', 'units': 'm', 'no_defs': None, 'datum': 'WGS84', 'type': 'crs'},
        #                                           MSG_HEIGHT, MSG_WIDTH,
        #                                           (lower_left_xy[0], lower_left_xy[1], upper_right_xy[0], upper_right_xy[1]))
        # #   +proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs
        # # print(msg_area_def.proj4_string)

        # msg_area_def = pr.geometry.AreaDefinition(area_id, description, proj_id,
        #                                           {'lat_0': '0.00', 'lat_ts': '0.00', 'lon_0': '0.00', 'proj': 'geos', 'h': '42164000.0', 'R': '6378000.0', 'sweep': 'y'},
        #                                           MSG_HEIGHT, MSG_WIDTH,
        #                                           (lower_left_xy[0], lower_left_xy[1], upper_right_xy[0], upper_right_xy[1]))
        # #   +proj=geos +lat_0=0.0 +lon_0=0.0 +h=42164000.0 +R=6378000.0 +sweep=y +type=crs
        # # print(msg_area_def.proj4_string)

        msg_area_def = pr.geometry.AreaDefinition(area_id, description, proj_id,
                                                  {'lat_0': '0.00', 'lat_ts': '0.00', 'lon_0': '0.00', 'proj': 'geos', 'h': '35786000.0', 'sweep': 'y'},
                                                  MSG_HEIGHT, MSG_WIDTH,
                                                  (lower_left_xy[0], lower_left_xy[1], upper_right_xy[0], upper_right_xy[1]))

        ##
        # Form a cartopy CRS
        #   <class 'pyresample.utils.cartopy.Projection'>
        msg_crs = msg_area_def.to_cartopy_crs()
        # print(msg_crs)

        #   (-57.18206024169922,  79.4859390258789, -78.96027374267578, 78.3008804321289)
        #   (-81.12566375732422, -81.0744857788086,  81.12566375732422, 81.0744857788086)
        # print(msg_crs.bounds)

        # # Kim
        # ##
        # # The issue is transform_point gives nan for values I assume are actually outside of field of view?
        # # msg_crs = ccrs.epsg(MERC_EPSG)
        # # msg_crs = ccrs.Geostationary()
        # lower_left_xy = [-47.18, -78.96]
        # upper_right_xy = [79.49, 78.31]
        # ur_xy = msg_crs.transform_point(upper_right_xy[0], upper_right_xy[1], geod_crs)
        # ll_xy = msg_crs.transform_point(lower_left_xy[0], lower_left_xy[1], geod_crs)
        # print(f"ll_xy: {ll_xy}")
        # print(f"ur_xy: {ur_xy}")

        # # ur_xya = geod_crs.transform_point(ur_xy[0], ur_xy[1], msg_crs)
        # # ll_xya = geod_crs.transform_point(ll_xy[0], ll_xy[1], msg_crs)
        # # print(f"ll_xy: {ll_xya}")
        # # print(f"ur_xy: {ur_xya}")
        # os._exit(1)
        if to_merc:
            ##
            # Form a Mercator CSR version of ccast_crs
            msg_merc_crs = ccrs.epsg(MERC_EPSG)
            if use_pseudo:
                ##
                # Mercator version EPSG:3857 for WGS 84 / Pseudo-Mercator -- Spherical Mercator, Google Maps, OpenStreetMap, Bing, ArcGIS, ESRI
                #   https://epsg.io/3857
                area_id = "EPSG:3857"
                description = "Pseudo-Mercator"
                proj_id = "EPSG:3857"
            else:
                ##
                # World Mercator version EPSG:3395 for WGS 84
                #   https://epsg.io/3395
                area_id = "EPSG:3395"
                description = "Mercator"
                proj_id = "EPSG:3395"

            #   ll_xy: (-9030867.579731662, 16261570.174893504)
            #   ur_xy: (-9025170.473223472, 16224751.37493146)
            ur_xy = msg_merc_crs.transform_point(upper_right_xy[0], upper_right_xy[1], geod_crs)
            ll_xy = msg_merc_crs.transform_point(lower_left_xy[0], lower_left_xy[1], geod_crs)
            # print(f"ll_xy: {ll_xy}")
            # print(f"ur_xy: {ur_xy}")

            ##
            # Specify projection parameters
            if use_pseudo:
                projection = {'proj': 'merc', 'lat_ts': 0, 'lon_0': 0, 'a': 6378137, 'b': 6378137, 'x_0': 0, 'y_0': 0, 'k': 1, 'units': 'm'}
            else:
                projection = {'proj': 'merc','x_0': 0, 'y_0': 0, 'k': 1, 'units': 'm', 'no_defs': None, 'datum': 'WGS84', 'type': 'crs'}

            ##
            # Define the area
            msg_merc_area_def = pr.geometry.AreaDefinition(area_id, 'Europe', proj_id,
                                                           projection, MSG_HEIGHT, MSG_WIDTH,
                                                           (ll_xy[0], ll_xy[1], ur_xy[0], ur_xy[1]))
            #   proj4_string: +proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs
            # print(f"proj4_string: {msg_merc_area_def.proj4_string}")
            msg_merc_crs = msg_merc_area_def.to_cartopy_crs()

            """
            msg_merc_crs
                msg_merc_area_corners  =[(-81.12565686373874, 81.07449269207677), UL
                                         (-81.07449267239406, 81.07449269207677), UR
                                         (-81.07449267239406, 81.12565688338495), LR
                                         (-81.12565686373874, 81.12565688338495)] LL

                msg_merc_area_def.area_extent (-9030867.579731662, 16261570.174893504, -9025170.473223472, 16224751.37493146)
                msg_merc_crs.bounds           (-9030867.579731662, -9025170.473223472, 16261570.174893504, 16224751.37493146)

                ll_xy  = (-9030867.579731662, 16261570.174893504) => (-81.12566375732422, 81.12566375732422)
                ur_xy  = (-9025170.473223472, 16224751.37493146)  => (-81.07449267239406, 81.07449269207677)
            """
            # print(msg_merc_area_def.corners)
            # print(msg_merc_area_def.area_extent)
            # print(msg_merc_crs.bounds)
            # ur_xya = geod_crs.transform_point(ur_xy[0], ur_xy[1], msg_merc_crs)
            # ll_xya = geod_crs.transform_point(ll_xy[0], ll_xy[1], msg_merc_crs)
            # print(f"ll_xy: {ll_xya}")
            # print(f"ur_xy: {ur_xya}")
            # os._exit(1)

    return (MSG_EPSG, MERC_EPSG, LCC_EPSG, geod_crs,
            ccast_area_def, ccast_crs, ccast_merc_area_def, ccast_merc_crs, ccast_lcc_area_def, ccast_lcc_crs,
            msg_area_def, msg_crs, msg_merc_area_def, msg_merc_crs)

# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
# >>>> END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE | END OF FILE <<<<
# >>>> ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: <<<<
