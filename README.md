<!-- vim: set fileencoding=utf-8 fileformat=unix : -->
<!-- vim: set spell spelllang=en : -->
<!-- -*- coding: utf-8 -*- -->
<!-- vim: set ts=8 et sw=4 sts=4 sta : -->

# dem2stl
Converts a Digital Elevation Model (DEM) topography surface to the STereoLithography (STL) format using the [GDAL](http://www.gdal.org/) library.



## Summary

The [`dem2stl.py`](http://github.com/cvr/dem2stl) python script converts a [Digital Elevation Model (DEM)](http://en.wikipedia.org/wiki/Digital_elevation_model) (such as GeoTIFF) to a [STereoLithography (STL)]("http://en.wikipedia.org/wiki/STL_(file_format)) file.
The result is a topography surface that may be used to generate computational grids such as those required by [`OpenFOAM`](http://openfoam.org/) meshing programs.

The [`dem2stl.py`](http://github.com/cvr/dem2stl) script can cope with DEM raster files referenced to geographic coordinates (longitude and latitude), internally projecting into an UTM coordinate system using [Proj4](http://github.com/jswhit/pyproj).


## Example

The following example uses the GeoTIFF sample file [`GeogToWGS84GeoKey5.tif`](http://download.osgeo.org/geotiff/samples/GeogToWGS84GeoKey/GeogToWGS84GeoKey5.tif):

```sh
./dem2stl.py example/GeogToWGS84GeoKey5.tif example/GeogToWGS84GeoKey5.stl
```

The [`GeogToWGS84GeoKey5.stl`](example/GeogToWGS84GeoKey5.png) output file
may be viewed using [meshlab](http://www.meshlab.net/):

![`GeogToWGS84GeoKey5.stl`](example/GeogToWGS84GeoKey5.png)



## Usage

```sh
usage: dem2stl.py [-h] [-v] [--band BAND] RASTER STL

Convert a GDAL raster (like a GeoTIFF heightmap) to an STL terrain surface.

positional arguments:
  RASTER         Input GeoTIFF image
  STL            Output STL path

optional arguments:
  -h, --help     show this help message and exit
  -v, --verbose  Print log messages
  --band BAND    Raster data band, defaults to 1
```


## Requirements to run

The following python modules are required to run
[`dem2stl.py`](http://github.com/cvr/dem2stl):
  * [GDAL/OGR](http://trac.osgeo.org/gdal/wiki/GdalOgrInPython)
  * [Proj4](http://github.com/jswhit/pyproj)
  * [NumPy](http://www.numpy.org/)
  * `argparse`
  * `collections`
  * `struct`
  * `numpy`

On Debian-based systems these modules may be installed running:
```sh
apt-get install libpython2.7-stdlib libpython2.7-minimal \
    python-numpy python-gdal python-pyproj
```


## Acknowledgments

[`dem2stl.py`](http://github.com/cvr/dem2stl) was based from the following
python scripts:
* [`phstl.py`](http://github.com/anoved/phstl) by Jim DeVona (MIT License, commit id `5717d88`, 2016-11-09). 
* The computation of the UTM zone was based on [`LatLongUTMconversion.py`](http://robotics.ai.uiuc.edu/~hyoon24/LatLongUTMconversion.py) by [Han Ul Yoon](mailto:hyoon24@uiuc.edu) (acessed 2017-05-12).


