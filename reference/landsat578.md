# Radiometric correction and processing of Landsat Images

`landsat578` performs radiometric collection of landsat 5,7,and 8
images. It returns albedo, NDVI, fractional vegetation cover, surface
temperature SAVI,radiance and reflectance of the image

## Usage

``` r
landsat578(
  data,
  welev,
  sensor = NULL,
  gain = "low",
  sunelev = NULL,
  date = NULL,
  weights = "auto",
  grescale = "auto",
  brescale = "auto",
  ESUNx = "auto",
  Mp = NULL,
  Ap = NULL,
  K1 = NULL,
  K2 = NULL,
  Rp = 0.91,
  Rsky = 1.32,
  tNB = 0.886,
  DOY = NULL,
  clip = NULL,
  folder = NULL
)

# Default S3 method
landsat578(
  data = NULL,
  welev,
  sensor = NULL,
  gain = "low",
  sunelev = NULL,
  date = NULL,
  weights = "auto",
  grescale = "auto",
  brescale = "auto",
  ESUNx = "auto",
  Mp = NULL,
  Ap = NULL,
  K1 = NULL,
  K2 = NULL,
  Rp = 0.91,
  Rsky = 1.32,
  tNB = 0.886,
  DOY = NULL,
  clip = NULL,
  folder = NULL
)
```

## Arguments

- data:

  An original folder containing Landsat data or Raster
  [`brick`](https://rdrr.io/pkg/raster/man/brick.html) or raster
  [`stack`](https://rdrr.io/pkg/raster/man/stack.html). Arrange the
  bands in ascending order: c(band1,band2,band3,band4,band5,band6,band7)
  for Landsat 5 and 7. And a similar order for Landsat 8. On folder see
  [`sebkcstack`](https://gowusu.github.io/sebkc/reference/sebkcstack.md)
  See also details.

- welev:

  Weather station elevation in meters

- sensor:

  Numeric. Type of landsat satellite; it takes 5,7, or 8

- gain:

  Band-specific sensor gain. It takes either "high" or "low"

- sunelev:

  Angle of Sun elevation in degrees, as found in the meta data of the
  satellite image

- date:

  DATE_ACQUIRED of the image, see metadata

- weights:

  albedo bands weights for each band. list as in order
  c(band1,band2,band3,band4,band5,band6,band7) for landsat 5 and 7. The
  default is "auto" where sebal weights are used.

- grescale:

  list as in order weights. Band-specific sensor values. The default for
  landsat 5 and 7 is from Chander et.al (2009). For Landsat 8 use ML
  Radiance scale in the metadata; i.e. Band specific multiplicative
  rescaling factor from the metadata (MTL file) (RADIANCE_MULT_BAND_x,
  where x is the band number).

- brescale:

  list. Band-specific sensor values. The default for landsat 5 and 7 is
  from Chander et.al (2009). for Landsat 8 use AL Radiance scale in the
  metadata; i.e. Band specific multiplicative rescaling factor from the
  metadata (MTL file) (RADIANCE_MULT_BAND_x, where x is the band
  number).

- ESUNx:

  list. The mean solar exo-atmospheric irradiance for each band. It is
  used for reflectance correction. The default is taken from SEBAL
  manual

- Mp:

  Reflectance scale for Landsat 8. Band specific multiplicative
  rescaling factor from the metadata (MTL file)
  (REFLECTANCE_MULT_BAND_x, where x is the band number).

- Ap:

  Reflectance scale for Landsat 8. Band specific additive rescaling
  factor from the metadata (MTL file) (REFLECTANCE_ADD_BAND_x, where x
  is the band number)

- K1:

  Temperature constants for each Landsat image. Default is taken from
  landsat handbook

- K2:

  Temperature constants for each Landsat image. Default is taken from
  landsat handbook

- Rp:

  path radiance in the 10.4-12.5 um band. MODTRAN runs can accurately
  estimate this. The default is taken from Allen et al(2007) to be 0.91,

- Rsky:

  narrow band downward thermal radiation. MODTRAN runs can accurately
  estimate this. The default is taken from Allen et al(2007) to be 1.32

- tNB:

  narrow band transmissivity of air (10.4-12.5) um range. MODTRAN runs
  can accurately estimate this. The default is taken from Allen et
  al(2007) to be 0.866

- DOY:

  Numeric or Date \[YYYY-mm-dd\]. Day of the Year. If you give data in
  the form of date \[YYYY-mm-dd\], it will be converted to DOY

- clip:

  extent object or raster object or polygon from which an Extent object
  can be extracted. A polygon or raster will be reprojected to conform
  the data to be cropped. for example clip can take the form of c(xmin,
  xmax, ymin, ymax)

- folder:

  An original directory of the images that contains all landsat 5,7 or 8
  bands and metadata. At the moment only Landsat folder is supported

## Value

- albedo::

  albedo values

- Ts::

  Surface Temperature

- NDVI::

  NDVI values

- SAVI::

  SAVI values

- fc::

  Fractional vegetation cover

- radiance::

  radiance values

- reflectance::

  TOA reflectance

## Details

If you set the data to landsat folder that contains landsat band and
meta data
[`sebkcstack`](https://gowusu.github.io/sebkc/reference/sebkcstack.md)
is used to compute the parameters. The user must specify welev and may
use default parameters. The data must always be a raster brick or stack.
To use the default values please order bind the bands
c(band1,band2,band3,band4,band5,band6,band7). The defaults values are
for landsat 5 and 7. For landsat 8, use the meta data.

## References

- Chander, G., Markham, B. L., & Helder, D. L. 2009. Summary of current
  radiometric calibration coefficients for Landsat MSS, TM, ETM+, and
  EO-1 ALI sensors. Remote Sensing of Environment, 113(5): 893-903.

- Waters, R. 2002. SEBAL: Surface Energy Balance Algorithms for Land.
  Idaho Implementation, Advanced Training and Users Manual. Kimberly,
  Idaho: University of Idaho.

- Survey, U. S. G. 2015. Landsat 8 (L8) data users handbook: 97p. from
  https://landsat.usgs.gov/documents/Landsat8DataUsersHandbook.pdf

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
#Automatic detection of parameters
folder=system.file("extdata","stack",package="sebkc")
modauto=landsat578(data=folder, welev=362)

#Manual estimation with input parameters
landsat7_6feb2004=brick(system.file("extdata","landsat7_6feb2004.grd",
package="sebkc"))
mod=landsat578(data=landsat7_6feb2004,welev=317.1,sensor=7,gain="high",
sunelev=50.71154048,weights="auto",grescale="auto",brescale="auto",
ESUN="auto",K1=666.09,K2=1282.71,date="2002-2-6")
#acces albedo value
albedo=mod$albedo
} # }
```
