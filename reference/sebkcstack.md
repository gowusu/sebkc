# Stacking, gap-filling, cloud removal of Landsat image ET Model

This function prepares landsat 5, 7, and 8 data for
[`landsat578`](https://gowusu.github.io/sebkc/reference/landsat578.md).
First, it reads meta data and returns many input parameters for
[`landsat578`](https://gowusu.github.io/sebkc/reference/landsat578.md)
and [`sebal`](https://gowusu.github.io/sebkc/reference/sebal.md).
Second, it stacks 7 landsat layers: it stacks bands 1-7 for landsat 5
and 7; and bands 2-7,10 for landsat 8.. Third, it removes landsat clouds
based on Martinuzzi et al(2009) method. Fourth, it gap-fills landsat 7
strips due to SLC-off as well as gap-filling NAs in other landsat
sensors. It can also process any selected bands.

## Usage

``` r
sebkcstack(
  folder,
  meta_pattern = "MTL.txt",
  remove_cloud = FALSE,
  gap_fill = FALSE,
  gap_loop = 8,
  clip = NULL
)

# Default S3 method
sebkcstack(
  folder,
  meta_pattern = "MTL.txt",
  remove_cloud = FALSE,
  gap_fill = FALSE,
  gap_loop = 8,
  clip = NULL
)
```

## Arguments

- folder:

  The folder or directory of landsat files

- meta_pattern:

  The regular expression pattern to identify metadata file. The default
  is "MTL.txt". The function works on the new landsat metadata only.

- remove_cloud:

  TRUE FALSE, 1, or 2. To indicate the cloud cover removal from the
  image. Default is FALSE. If it is set to 1, the overlapping areas of
  band 1 DN values between 120 and 255 and band 6 DN values of 102 to
  128 are used. If it is set to 2 the union of the DN ranges are used;
  and clouds as well as urban, barren, quarries, rocks, and sand are
  removed; otherwise only clouds are removed. See Martinuzzi et al(2009)
  for more details. At the moment only landsat 5 and 7 clouds are
  removed.

- gap_fill:

  Logical (TRUE or FALSE).It replaces strips in landsat 7 or NAs in
  other data set. It uses
  [`focal`](https://rspatial.github.io/terra/reference/focal.html)
  neighbourhood approach. It is an automated Version of USGS
  "Gap-Filling Landsat 7 SLC-off Single Scenes Using Erdas Imagine TM".
  Like USGS, it also takes approximately one hour to complete. see
  details at http://landsat.usgs.gov/ERDAS_Approach.php

- gap_loop:

  numeric. The number of loops or repetitions that gap-fill should be
  run. The default is 8 for one landsat scene. If the landsat image is
  smaller, it can be reduced.

- clip:

  extent object or raster object or polygon from which an Extent object
  can be extracted. A polygon or raster will be reprojected to conform
  the data to be cropped. for example clip can take the form of c(xmin,
  xmax, ymin, ymax)

## Value

- data::

  raster stacked data. If crop is not NULL then crop Extent object is
  returned

- brescale::

  brescale of the stacked data

- grescale::

  grescale of the stacked data

- sunelev::

  sun elevation for stacked data

- date::

  Acquired date of the raster stacked data

- viewangle::

  view angle of landsat 8

- azimuth::

  azimuth angle of the image

- sensor::

  Type of satellite

- K1::

  K1 of landsat 8 thermal band 10

- K2::

  K2 of landsat 8 thermal band 10

- cc::

  cloud cover percentage

- quality::

  Image quality

- row::

  WRS_ROW of the image

- path::

  WRS_PATH of the image

## Details

The function can simulate any of the activities described in the
description section. It takes several minutes to replace strips (1 hour)
and remove clouds (30 minutes) on one landsat scene. But stacking of
bands and reading of metadata take few seconds. One can save the stacked
and processed data in the same "folder" with
[`writesebkc`](https://gowusu.github.io/sebkc/reference/writesebkc.md).
The output file name is "sebkcstack.tif".

## References

Martinuzzi, S., Gould, W.A., Ramos Gonzales, O.M. 2007. Creating
Cloud-Free Landsat ETM+ Data Sets in Tropical Landscapes: Cloud and
Cloud-Shadow Removal. USDA Forest Service General Technical Report
IITF-GTR-32.

USGS Using Landsat 7 Data http://landsat.usgs.gov/ERDAS_Approach.php

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
folder=system.file("extdata","stack",package="sebkc")
stack=sebkcstack(folder=folder)
#returns data and grescale
stack$data
stack$grescale
} # }
```
