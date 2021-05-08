#!/usr/bin/python
# vim: set fileencoding=utf-8 fileformat=unix :
# -*- coding: utf-8 -*-
# vim: set ts=8 et sw=4 sts=4 sta :

########################################################################
# Copyright (C) 2017 by Carlos Veiga Rodrigues. All rights reserved.
# Author:  Carlos Veiga Rodrigues <cvrodrigues@gmail.com>
#
# This program can be redistribuited and modified
# under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation,
# either version 3 of the License or any later version.
# This program is distributed WITHOUT ANY WARRANTY,
# without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.
# For more details consult the GNU Lesser General Public License
# at http://www.gnu.org/licenses/lgpl.html.
#
# This script was initially based from:
# * phstl.py from Jim DeVona at http://github.com/anoved/phstl (MIT License)
#   from commit id 5717d88, 2016-11-09
# * function ComputeUTMProj4defs based on the LatLongUTMconversion.py
#   by Han Ul Yoon <hyoon24@uiuc.edu>, acessed 2017-05-12
#   http://robotics.ai.uiuc.edu/~hyoon24/LatLongUTMconversion.py
#
# ChangeLog (date - version - author):
# * May 2017 - 1.0 - Carlos Veiga Rodrigues
#
########################################################################


import sys
import argparse
import numpy as np
from collections import deque
from struct import pack, unpack
import pyproj
try:
    from osgeo import gdal
    from osgeo.gdalconst import *
    from osgeo import osr
    from osgeo import ogr
except ImportError:
    import gdal
    from gdalconst import *
    import osr
    import ogr

gdal.UseExceptions()
gdal.TermProgress = gdal.TermProgress_nocb

swap_GeoTransform = True


## parser
ap = argparse.ArgumentParser(\
    description="Convert a GDAL raster (like a GeoTIFF heightmap)"\
    +" to an STL terrain surface.")
ap.add_argument('-v', '--verbose', action='store_true', default=False,\
    help="Print log messages")
ap.add_argument('--band', action='store', default=1, type=int,\
    help="Raster data band, defaults to 1")
ap.add_argument('RASTER', help="Input GeoTIFF image")
ap.add_argument('STL', help="Output STL path")
args = ap.parse_args()


def fail(msg):
    print(msg, file=sys.stderr)
    exit(1)

def verbose(msg):
    if args.verbose:
        print(msg, file=sys.stderr)
    return

def ComputeUTMProj4defs (lon, lat):
    """
    Estimate UTM Zone and exports definitions for Proj4 input.
    Equations from USGS Bulletin 1532. Coordinates in decimal degrees.
    Based from:  http://robotics.ai.uiuc.edu/~hyoon24/LatLongUTMconversion.py
    """
    # Make sure the longitude is between -180.0 and 179.9
    lon = (lon + 180) - int((lon + 180) / 360) * 360 - 180
    Zone = int((lon + 180.)/6.) + 1 
    if lat >= 56 and lat < 64 and lon >= 3 and lon < 12:
        Zone = 32
    # Special zones for Svalbard
    if lat >= 72 and lat < 84:
        if  lon >= 0  and lon < 9:
            Zone = 31
        elif lon >= 9 and lon < 21:
            Zone = 33
        elif lon >= 21 and lon < 33:
            Zone = 35
        elif lon >= 33 and lon < 42:
            Zone = 37
    Proj4def = "+proj=utm +ellps=WGS84 +zone=%d" % Zone
    if lat < 0:
        Proj4def += " +south"
    return pyproj.Proj(Proj4def)

def NormalVector(t):
    """
    Calculate the normal vector of a triangle (unit vector perpendicular to
    triangle surface, pointing away from the "outer" face of the surface).
    Computed using 32-bit float operations for consistency.

    Parameters:  triangle vertices (nested x y z tuples)

    Returns:     normal vector (x y z tuple)
    """
    (ax, ay, az) = t[0]
    (bx, by, bz) = t[1]
    (cx, cy, cz) = t[2]
    # first edge
    e1x = np.float32(ax) - np.float32(bx)
    e1y = np.float32(ay) - np.float32(by)
    e1z = np.float32(az) - np.float32(bz)
    # second edge
    e2x = np.float32(bx) - np.float32(cx)
    e2y = np.float32(by) - np.float32(cy)
    e2z = np.float32(bz) - np.float32(cz)
    # cross product
    cpx = np.float32(e1y * e2z) - np.float32(e1z * e2y)
    cpy = np.float32(e1z * e2x) - np.float32(e1x * e2z)
    cpz = np.float32(e1x * e2y) - np.float32(e1y * e2x)
    # return cross product vector normalized to unit length
    mag = np.sqrt(np.power(cpx, 2) + np.power(cpy, 2) + np.power(cpz, 2))
    return (cpx/mag, cpy/mag, cpz/mag)

class stlwriter():
    """
    stlwriter is a simple class for writing binary STL meshes.
    Class instances are constructed with a predicted face count.
    The output file header is overwritten upon completion with
    the actual face count.
    """
    # path: output binary stl file path
    # facet_count: predicted number of facets
    def __init__(self, path, facet_count=0):
        self.f = open(path, 'wb')
        # track number of facets actually written
        self.written = 0
        # write binary stl header with predicted facet count
        self.f.write(('\0' * 80).encode())
        # (facet count is little endian 4 byte unsigned int)
        self.f.write(pack('<I', facet_count))
        return
    #
    # t: ((ax, ay, az), (bx, by, bz), (cx, cy, cz))
    def add_facet(self, t):
        # facet normals and vectors are little endian 4 byte float triplets
        # strictly speaking, we don't need to compute NormalVector,
        # as other tools could be used to update the output mesh.
        self.f.write(pack('<3f', *NormalVector(t)))
        for vertex in t:
            self.f.write(pack('<3f', *vertex))
        # facet records conclude with two null bytes (unused "attributes")
        self.f.write(('\0\0').encode())
        self.written += 1
        return
    #
    def done(self):
        # update final facet count in header before closing file
        self.f.seek(80)
        self.f.write(pack('<I', self.written))
        self.f.close()
        return
    #
    def __enter__(self):
        return self
    #
    def __exit__(self, exc_type, exc_value, traceback):
        self.done()
        return



## open dataset
try:
    ds = gdal.Open(args.RASTER)
except RuntimeError as e:
    fail(str(e).strip())

## dataset coordinate system
OSRds = osr.SpatialReference()
OSRds.ImportFromWkt(ds.GetProjection())    # ds.GetProjection() in Wkt
P4ds = pyproj.Proj(OSRds.ExportToProj4())  # get projection in Proj4
EPSGds = None
if 0 == OSRds.AutoIdentifyEPSG():
    EPSGds = OSRds.GetAuthorityCode(None)  # guess EPSG code


## raster dimensions
nx, ny = ds.RasterXSize, ds.RasterYSize
nxm1, nym1 = nx - 1, ny - 1
verbose("raster size %dx%d pixels" % (nx, ny))

## geo transforms
GT = ds.GetGeoTransform()
# i, j = nx/3, ny*2/3
# x, y = gdal.ApplyGeoTransform(GT, i+.5, j+.5)
## using arrays instead...
## i = np.arange(nx) ; j = np.arange(ny) ;
# x = GT[0] + GT[1] * (i[None,:] + .5) + GT[2] * (j[:,None] + .5)
# y = GT[3] + GT[4] * (i[None,:] + .5) + GT[5] * (j[:,None] + .5)
IGT = gdal.InvGeoTransform(GT)[1]
# i, j = np.round(np.array(gdal.ApplyGeoTransform(IGT, x, y)) - .5).astype(int)
## using arrays instead...
# i = np.round(IGT[0] + IGT[1] * x[:] + IGT[2] * y[:] - .5).astype(int)
# j = np.round(IGT[3] + IGT[4] * x[:] + IGT[5] * y[:] - .5).astype(int)

## get UTM projected coordinates
if 1 == OSRds.IsProjected():
    # i = np.arange(nx)
    # j = np.arange(ny)
    i = np.array([0, nxm1])
    j = np.array([0, nym1])
    x = GT[0] + GT[1] * (i[None,:] + .5) + GT[2] * (j[:,None] + .5)
    y = GT[3] + GT[4] * (i[None,:] + .5) + GT[5] * (j[:,None] + .5)
    P4xy = P4ds
    OSRxy = OSRds
    EPSGxy = EPSGds
    print("Projected CS,  EPSG:%s" % EPSGxy)
elif 1 == OSRds.IsGeographic():
    print("Geographic CS, EPSG:%s" % EPSGds)
    # i = np.arange(nx)
    # j = np.arange(ny) 
    i = np.array([0, nxm1])
    j = np.array([0, nym1])
    lon = GT[0] + GT[1] * (i[None,:] + .5) + GT[2] * (j[:,None] + .5)
    lat = GT[3] + GT[4] * (i[None,:] + .5) + GT[5] * (j[:,None] + .5)
    ## compute UTM zone from centre lon lat
    lonc, latc = gdal.ApplyGeoTransform(GT, nx/2., ny/2.)
    P4xy = ComputeUTMProj4defs (lonc, latc)
    OSRxy = osr.SpatialReference()
    OSRxy.ImportFromProj4(P4xy.srs)
    EPSGxy = None
    if 0 == OSRxy.AutoIdentifyEPSG():
        EPSGxy = OSRxy.GetAuthorityCode(None)  # guess EPSG code
    x, y = pyproj.transform(P4ds, P4xy, lon, lat)
    print("Projected CS,  EPSG:%s" % EPSGxy)
else:
    print("??  Error, cannot determine Coordinate Reference System")
    print("??  Wkt string: " + ds.GetProjection())
    print("??  Proj4 definitions: " + P4ds.srs)
    print("??  EPSG:%s" % EPSGds)
    fail("Aborting")

verbose("xmin, xmax = %g m, %g m" % (np.min(x), np.max(x)))
verbose("ymin, ymax = %g m, %g m" % (np.min(y), np.max(y)))


verbose("Working on band %d" % args.band)
band = ds.GetRasterBand(args.band)
nd = band.GetNoDataValue()

## map GDAL pixel data type to corresponding struct format character
typemap = { \
    gdal.GDT_Byte:    'B',\
    gdal.GDT_UInt16:  'H',\
    gdal.GDT_Int16:   'h',\
    gdal.GDT_UInt32:  'I',\
    gdal.GDT_Int32:   'i',\
    gdal.GDT_Float32: 'f',\
    gdal.GDT_Float64: 'd'}

typeName = gdal.GetDataTypeName(band.DataType)
if band.DataType not in typemap:
    fail("Unsupported data type: %s" % typeName)

## rowformat is used to unpack a row of raw image data to numeric form
rowformat = typemap.get(band.DataType) * nx
verbose("data type   = %s" % typeName)
verbose("type format = %s" % typemap.get(band.DataType))

## statistics on z field
zmin, zmax, zavg, zstd = band.GetStatistics(True, True)
verbose("min(z) = %g" % zmin)
verbose("max(z) = %g" % zmax)
verbose("avg(z) = %g" % zavg)
verbose("std(z) = %g" % zstd)


## rolling pixel buffer has space for two rows of image data
## old data is automatically discarded as new data is loaded
pixels = deque(maxlen = (2 * nx))
verbose("buffer size = %s" % str(pixels.maxlen))

## initialize pixel buffer with first row of data from the image window
pixels.extend(unpack(rowformat, \
    band.ReadRaster(0, 0, nx, 1, nx, 1, band.DataType)))

## precalculate output mesh size (STL is 50 bytes/facet + 84 byte header)
## actual facet count and file size may differ (be less) if pixels are
## skipped as nodata or out of range
facetcount = nxm1 * nym1 * 2
filesize = (facetcount * 50) + 84
verbose("predicted (max) facet count = %s" % str(facetcount))
verbose("predicted (max) STL file size = %s bytes" % str(filesize))

i0, j0 = 0, 0
with stlwriter(args.STL, facetcount) as mesh:
    for j in range(nym1):
        ## each row, extend pixel buffer with the next row of data
        ## from the image window
        pixels.extend(unpack(rowformat, \
            band.ReadRaster(i0, j0 + j + 1, nx, 1, nx, 1, band.DataType)))
        for i in range(nxm1):
            ## z values of this pixel (a) and its neighbors (b, c and d)
            av = pixels[i]
            bv = pixels[nx + i]
            cv = pixels[i + 1]
            dv = pixels[nx + i + 1]
            ## apply transforms to obtain output mesh coordinates of the
            ## four corners composed of raster points a (x, y), b, c,
            ## and d (x + 1, y + 1):
            ##
            ## a-c   a-c     c
            ## |/| = |/  +  /|
            ## b-d   b     b-d
            ##
            ## points b and c are required for both facets, so if either
            ## are not valid, we can skip this pixel altogether
            if nd == bv or nd == cv:
                ## point b or c invalid
                continue
            ## remember that we are looking for the vertexes, not the centres
            b = [ \
                GT[0] + (i0 + i) * GT[1] + (j0 + j + 1) * GT[2],\
                GT[3] + (i0 + i) * GT[4] + (j0 + j + 1) * GT[5],\
                float(bv) ]
            c = [ \
                GT[0] + (i0 + i + 1) * GT[1] + (j0 + j) * GT[2],\
                GT[3] + (i0 + i + 1) * GT[4] + (j0 + j) * GT[5],\
                float(cv) ]
            if 1 == OSRds.IsGeographic():
                b[0], b[1] = pyproj.transform(P4ds, P4xy, b[0], b[1])
                c[0], c[1] = pyproj.transform(P4ds, P4xy, c[0], c[1])
            if nd != av:
                ## point a is valid
                a = [ \
                    GT[0] + (i0 + i) * GT[1] + (j0 + j) * GT[2],\
                    GT[3] + (i0 + i) * GT[4] + (j0 + j) * GT[5],\
                    float(av) ]
                if 1 == OSRds.IsGeographic():
                    a[0], a[1] = pyproj.transform(P4ds, P4xy, a[0], a[1])
                mesh.add_facet((a, b, c))
            if nd != dv:
                ## point d is valid
                d = [ \
                    GT[0] + (i0 + i + 1) * GT[1] + (j0 + j + 1) * GT[2],\
                    GT[3] + (i0 + i + 1) * GT[4] + (j0 + j + 1) * GT[5],\
                    float(dv) ]
                if 1 == OSRds.IsGeographic():
                    d[0], d[1] = pyproj.transform(P4ds, P4xy, d[0], d[1])
                mesh.add_facet((d, c, b))
        ## update progress each row
        gdal.TermProgress(float(j + 1) / nym1)

verbose("actual facet count: %s" % str(mesh.written))

