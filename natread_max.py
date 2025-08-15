#! /usr/bin/env python -tt
# -*- coding: utf-8; mode: python -*-
r"""
Read Meteosat Second Generation (MSG) Native Archive Format (Max .nat) file.

natread_max.py
~~~~~~~~~~~~~~
"""
# Standard Imports
import os


# Third-Party Imports
import numpy as np
from osgeo import gdal
from osgeo import osr
# import pyresample as pr
#   conda install -c conda-forge satpy
# from satpy import Scene
import cartopy
import cartopy.crs as ccrs
import pprint

# Local Imports

##
# List of Public objects from this module.
__all__ = ['natread_max']

##
# Markup Language Specification (see NumpyDoc Python Style Guide https://numpydoc.readthedocs.io/en/latest/format.html)
__docformat__ = "Numpydoc"
# ------------------------------------------------------------------------------

def get_pixel_coords(geotransform, lon, lat):
    x_origin = geotransform[0]
    y_origin = geotransform[3]
    pixel_width = geotransform[1]
    pixel_height = geotransform[5]
    x_pixel = round((lon - x_origin) / pixel_width)
    y_pixel = round((lat - y_origin) / pixel_height)
    return x_pixel, y_pixel


# Define Global Constants and State Variables
# -------------------------------------------

##
# From satpy.scene.Scene
#   ['HRV', 'IR_016', 'IR_039', 'IR_087', 'IR_097', 'IR_108', 'IR_120', 'IR_134', 'VIS006', 'VIS008', 'WV_062', 'WV_073']


# CloudCast: A Satellite-Based Dataset and Baseline for Forecasting Clouds
#   Nielsen and Iosifidis (2021)
#     Primary  : IR 8.7, IR 10.8, IR 12.0 um
#     Secondary: VIS 0.6, WV 7.3 um
#
#  MSG Instrument Channels: VIS 0.6  Band 1   * Used by Nielsen and Iosifidis (2021)
#                           VIS 0.8  Band 2
#                           IR 1.6   Band 3
#                           IR 3.9   Band 4
#                           WV 6.2   Band 5
#                           WV 7.3   Band 6   * Used by Nielsen and Iosifidis (2021)
#                           IR 8.7   Band 7   * Used by Nielsen and Iosifidis (2021)
#                           IR 9.7   Band 8
#                           IR 10.8  Band 9   * Used by Nielsen and Iosifidis (2021)
#                           IR 12.0  Band 10  * Used by Nielsen and Iosifidis (2021)
#                           IR 13.4  Band 11
#                           HRV      Band 12
band_ids = [1, 3, 4, 7, 9, 10]
band_def = ('VIS 0.6', 'IR 1.6', 'IR 3.9', 'IR 8.7', 'IR 12.0')

band_ids_nielsen = [1, 6, 7, 9, 10]
band_def_nielsen = ('VIS 0.6', 'WV 7.3', 'IR 8.7', 'IR 10.8', 'IR 12.0')
band_ids = band_ids_nielsen
band_def = band_def_nielsen

filename = "/Users/mbauer/tmp/CloudCast/msg/20170601001240.nat"
# filename = "/Users/mbauer/tmp/CloudCast/msg/MSG3-SEVI-MSG15-0100-NA-20170601001240.835000000Z-NA/MSG3-SEVI-MSG15-0100-NA-20170601001240.835000000Z-NA.nat"

gdal.UseExceptions()

##
# Info Dump
# with gdal.Open(filename) as ds:
#     info = gdal.Info(ds, format='json')
#     del info["stac"]  # to avoid cluttering below output
#     pprint.pprint(info, indent=2, width=100)
r"""
    { 'bands': [ { 'band': 1,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 01',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 2,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 02',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 3,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 03',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 4,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 04',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 5,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 05',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 6,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 06',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 7,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 07',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 8,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 08',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 9,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 09',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 10,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 10',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'},
                 { 'band': 11,
                   'block': [3712, 1],
                   'colorInterpretation': 'Undefined',
                   'description': 'band 11',
                   'max': 1023.0,
                   'metadata': {},
                   'min': 1.0,
                   'noDataValue': 0.0,
                   'type': 'UInt16'}],
      'coordinateSystem': { 'dataAxisToSRSAxisMapping': [1, 2],
                            'wkt': 'PROJCRS["Geostationary projection (MSG)",\n'
                                   '    BASEGEOGCRS["MSG Ellipsoid",\n'
                                   '        DATUM["MSG_DATUM",\n'
                                   '            ELLIPSOID["MSG_SPHEROID",6378137,298.257024882281,\n'
                                   '                LENGTHUNIT["metre",1,\n'
                                   '                    ID["EPSG",9001]]]],\n'
                                   '        PRIMEM["Greenwich",0,\n'
                                   '            ANGLEUNIT["degree",0.0174532925199433,\n'
                                   '                ID["EPSG",9122]]]],\n'
                                   '    CONVERSION["Geostationary Satellite (Sweep Y)",\n'
                                   '        METHOD["Geostationary Satellite (Sweep Y)"],\n'
                                   '        PARAMETER["Longitude of natural origin",0,\n'
                                   '            ANGLEUNIT["degree",0.0174532925199433],\n'
                                   '            ID["EPSG",8802]],\n'
                                   '        PARAMETER["Satellite Height",35785863,\n'
                                   '            LENGTHUNIT["metre",1,\n'
                                   '                ID["EPSG",9001]]],\n'
                                   '        PARAMETER["False easting",0,\n'
                                   '            LENGTHUNIT["metre",1],\n'
                                   '            ID["EPSG",8806]],\n'
                                   '        PARAMETER["False northing",0,\n'
                                   '            LENGTHUNIT["metre",1],\n'
                                   '            ID["EPSG",8807]]],\n'
                                   '    CS[Cartesian,2],\n'
                                   '        AXIS["(E)",east,\n'
                                   '            ORDER[1],\n'
                                   '            LENGTHUNIT["metre",1,\n'
                                   '                ID["EPSG",9001]]],\n'
                                   '        AXIS["(N)",north,\n'
                                   '            ORDER[2],\n'
                                   '            LENGTHUNIT["metre",1,\n'
                                   '                ID["EPSG",9001]]]]'},
      'cornerCoordinates': { 'center': [0.0, 1500.202],
                             'lowerLeft': [-5568748.276, -5567248.074],
                             'lowerRight': [5568748.276, -5567248.074],
                             'upperLeft': [-5568748.276, 5570248.477],
                             'upperRight': [5568748.276, 5570248.477]},
      'description': '/Users/mbauer/tmp/CloudCast/msg/20170601001240.nat',
      'driverLongName': 'EUMETSAT Archive native (.nat)',
      'driverShortName': 'MSGN',
      'files': ['/Users/mbauer/tmp/CloudCast/msg/20170601001240.nat'],
      'geoTransform': [ -5568748.275756836,
                        3000.4031658172607,
                        0.0,
                        5570248.477339745,
                        0.0,
                        -3000.4031658172607],
      'metadata': { '': { 'Date/Time': '20170601/00:12',
                          'Origin': '1 1',
                          'Radiometric parameters format': 'offset slope',
                          'ch01_cal': '-1.065267613158e+00 2.088760025799e-02',
                          'ch02_cal': '-1.421905456111e+00 2.788049913943e-02',
                          'ch03_cal': '-1.202993122861e+00 2.358810044825e-02',
                          'ch04_cal': '-1.865920103496e-01 3.658666869601e-03',
                          'ch05_cal': '-4.242236706827e-01 8.318111189856e-03',
                          'ch06_cal': '-1.969720376950e+00 3.862196817550e-02',
                          'ch07_cal': '-6.463960252982e+00 1.267443186859e-01',
                          'ch08_cal': '-5.302006445576e+00 1.039609106976e-01',
                          'ch09_cal': '-1.045681948659e+01 2.050356762077e-01',
                          'ch10_cal': '-1.133786848073e+01 2.223111466809e-01',
                          'ch11_cal': '-8.037951740768e+00 1.576068968778e-01'}},
      'size': [3712, 3712],
      'wgs84Extent': {'coordinates': [[]], 'type': 'Polygon'}}
"""

"""From https://svn.osgeo.org/gdal/tags/1.7.0/gdal/frmts/msgn/frmt_msgn.html
    MSGN -- Meteosat Second Generation (MSG) Native Archive Format (.nat)

    GDAL supports reading only of MSG native files. These files may have anything from 1 to 12 bands, all at 10-bit resolution.

    Includes support for the 12th band (HRV - High Resolution Visible).
    This is implemented as a subset, i.e., it is accessed by prefixing the filename with the tag "HRV:".

    Similarly, it is possible to obtain floating point radiance values in stead of the usual 10-bit digital numbers (DNs).
    This subset is accessed by prefixing the filename with the tag "RAD:".

    Georeferencing is currently supported, but the results may not be acceptable (accurate enough), depending on your requirements.
    The current workaround is to implement the CGMS Geostationary projection directly, using the code available from EUMETSAT.
"""


dataset = gdal.Open(f"RAD:{filename}")

print(f"\nReading: {dataset.GetDescription()}")
print(f"\tRasterCount: {dataset.RasterCount}")
# print(f"\tGetMetadata: {dataset.GetMetadata()}")

upper_left_xy = [-855100.436345, -4942000.0]
lower_right_xy = [1448899.563655, -2638000.0]

img = []
for bidx, band_id in enumerate(band_ids):
    print(f"\tReading Band {band_id}: {band_def[bidx]}")
    band = dataset.GetRasterBand(band_id)
    # print(band.GetDataset().GetDriver())
    # print(band.GetDataset().GetLayerCount())

    # print(band.GetHistogram())
    # break
    geotransform = dataset.GetGeoTransform()

    x_min, y_max = get_pixel_coords(geotransform, upper_left_xy[0], upper_left_xy[1])
    x_max, y_min = get_pixel_coords(geotransform, lower_right_xy[0], lower_right_xy[1])
    print(f"\t{x_min = }, {y_max = }, {x_max = }, {y_min = }")
    break
    # data_subset = band.ReadAsArray(x_min, y_min, x_max - x_min, y_max - y_min).tolist()
    #   ndarray (768, 768)
    data_subset = band.ReadAsArray(x_min, y_min, x_max - x_min, y_max - y_min)
    #img.append(data_subset)

    tmp = data_subset.flatten()
    tmp_len = len(tmp)
    tmp = tmp[~np.isnan(tmp)]
    tmp = tmp[tmp > 0.0]
    tmp1_len = len(tmp)
    len_frac = 100.0 * (tmp1_len / tmp_len)
    print(f"\t\tdata_subset {data_subset.shape} {len_frac:5.2f}%: [{np.amin(tmp)}, ... {np.amax(tmp)}]")
    del tmp

"""
First question:
    What 'bands' are we using? Based on Nielsen and Iosifidis (2021) I thought we were using ('VIS 0.6', 'WV 7.3', 'IR 8.7', 'IR 10.8', 'IR 12.0'),
    which if following the order in MSG_provenance.txt are bands (1, 6, 7, 9, 10). But the ones listed in Max's code are bands (1, 3, 4, 7, 9, 10).

    That could be just a misunderstanding on my part.

    CloudCast: A Satellite-Based Dataset and Baseline for Forecasting Clouds
      Nielsen and Iosifidis (2021)
        Primary  : IR 8.7, IR 10.8, IR 12.0 um
        Secondary: VIS 0.6, WV 7.3 um

    From MSG_provenance.txt
     MSG Instrument Channels: VIS 0.6  Band 1   * Used by Nielsen and Iosifidis (2021)
                              VIS 0.8  Band 2
                              IR 1.6   Band 3
                              IR 3.9   Band 4
                              WV 6.2   Band 5
                              WV 7.3   Band 6   * Used by Nielsen and Iosifidis (2021)
                              IR 8.7   Band 7   * Used by Nielsen and Iosifidis (2021)
                              IR 9.7   Band 8
                              IR 10.8  Band 9   * Used by Nielsen and Iosifidis (2021)
                              IR 12.0  Band 10  * Used by Nielsen and Iosifidis (2021)
                              IR 13.4  Band 11
                              HRV      Band 12

Second question:
    Even assuming the bands Max lists are correct the values I get are not, except perhaps for the visible channel (band 1).
    This is a hint because if you look at the .nat file specification there are 3 raster fields for each channel/band.
         DataID(name='VIS006', wavelength=WavelengthRange(min=0.56, central=0.635, max=0.71, unit='µm'), resolution=3000.403165817, calibration=<calibration.reflectance>, modifiers=()),
         DataID(name='VIS006', wavelength=WavelengthRange(min=0.56, central=0.635, max=0.71, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='VIS006', wavelength=WavelengthRange(min=0.56, central=0.635, max=0.71, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),

    In the base of the 'VIS006', the first raster is the one we want which is calibrated to 'reflectance' and has a range from
    0 to 100. And indeed the values in Max's file for band=1 are from 64 to 70. However, isn't this image at night over Europe?

    The source file 'MSG3-SEVI-MSG15-0100-NA-20170601001240.835000000Z-NA.zip' suggests YYYY = 2017, MM = 06, DD=01, HH=00 or just after 00 UTC.
    In that case I wouldn't expect more than a little reflectance due to maybe moonlight or city lights?
    Anyway, when I read that band I only get reflectance values around 1 to 2.

    My concern is the third/last raster 'counts'. My code reads this as ranging from  65 to 77, which is similar to what Max's file has.

    For the IR bands the first of the 3 rasters are as follows
         DataID(name='IR_087', wavelength=WavelengthRange(min=8.3, central=8.7, max=9.1, unit='µm'), resolution=3000.403165817, calibration=<calibration.brightness_temperature>, modifiers=()),
         DataID(name='IR_087', wavelength=WavelengthRange(min=8.3, central=8.7, max=9.1, unit='µm'), resolution=3000.403165817, calibration=<calibration.radiance>, modifiers=()),
         DataID(name='IR_087', wavelength=WavelengthRange(min=8.3, central=8.7, max=9.1, unit='µm'), resolution=3000.403165817, calibration=<calibration.counts>, modifiers=()),

    Here Max's file suggests values between 109 and 530

    This is what I get if I read the 'counts' calibration (between 108 and 521).
    If I read the brightness_temperature calibration I get values between 211 and 288, which seems correct of degrees Kelvin.

So there are three issues. The first is Max has read the 'count' rather than 'reflectance' or 'brightness_temperature' calibration.
Second, we need to be on board about the bands we are reading. The third and more minor issue is if we can/should do the reduction to integers and if that is done by truncation or some form of rounding?

Mike

Results
    Reading: /Users/mbauer/tmp/CloudCast/msg/20170601001240.nat
        Reading Band 1: VIS 0.6
            data_subset (768, 768) 100.00%: [64, ... 70]
        Reading Band 3: IR 1.6
            data_subset (768, 768) 100.00%: [87, ... 108]
        Reading Band 4: IR 3.9
            data_subset (768, 768) 100.00%: [51, ... 224]
        Reading Band 7: IR 8.7
            data_subset (768, 768) 100.00%: [109, ... 530]
        Reading Band 9: IR 12.0
            data_subset (768, 768) 100.00%: [134, ... 535]

Using the Nielsen and Iosifidis (2021) channel over from MSG_provenance.txt
    Reading: /Users/mbauer/tmp/CloudCast/msg/20170601001240.nat
        Reading Band 1: VIS 0.6
            data_subset (768, 768) 100.00%: [64, ... 70]
        Reading Band 6: WV 7.3
            data_subset (768, 768) 100.00%: [124, ... 582]
        Reading Band 7: IR 8.7
            data_subset (768, 768) 100.00%: [109, ... 530]
        Reading Band 9: IR 10.8
            data_subset (768, 768) 100.00%: [134, ... 535]
        Reading Band 10: IR 12.0
            data_subset (768, 768) 100.00%: [155, ... 558]

"""
