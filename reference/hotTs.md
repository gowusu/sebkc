# Automatic computation of hot pixel Radiometric Surface Temperature (Ts)

This function computes temperature of a hot pixel based on NDVI, Surface
Temperature, DEM and albedo. The main difference with
[`hotTs2`](https://gowusu.github.io/sebkc/reference/hotTs2.md) is that
preference is given to temperature selection before considering NDVI and
albedo selection. Moreover, this is a simple end member selection where
no clustering is done.

## Usage

``` r
hotTs(
  Ts,
  NDVI,
  albedo = NULL,
  DEM = NULL,
  cluster = 8,
  extent = "interactive",
  upper = 0.8,
  lower = 0.2,
  plot = TRUE,
  layout = "portrait",
  draw = "poly",
  folder = NULL,
  welev = NULL,
  clip = NULL,
  iter.max = 500,
  seed = NULL
)

# Default S3 method
hotTs(
  Ts = NULL,
  NDVI,
  albedo = NULL,
  DEM = NULL,
  cluster = 8,
  extent = "interactive",
  upper = 0.8,
  lower = 0.2,
  plot = TRUE,
  layout = "portrait",
  draw = "poly",
  folder = NULL,
  welev = NULL,
  clip = NULL,
  iter.max = 500,
  seed = NULL
)
```

## Arguments

- Ts:

  A RasterLayer data that indicates radiometric surface temperature
  values from remote sensing image, preferably in Kelvin (K). You can
  also point to raster file on your computer.

- NDVI:

  A RasterLayer data that indicates NDVI (Normalized Difference
  Vegetation Index) values. You can also point to a file on your
  computer.

- albedo:

  A RasterLayer data that has albedo values. You can also provide file
  path or location on your computer

- DEM:

  A digital elevation model\[m\]

- cluster:

  numeric, a value indicating the number of clusters.

- extent:

  A subset of the raster file than include c(xmin, xmax; serastercond
  row: ymin, ymax). You can set it to "interactive" where you can
  digitise polygon. You can set to "auto" where it will be automatically
  generated. See layout parameter below for more details. see
  [`extent`](https://rdrr.io/pkg/raster/man/extent.html)

- upper:

  The quantile probability at which Ts should be included in the
  selection project. It ranges from 0-1. The default is 0.8

- lower:

  The quantile probability at which NDVI should be included in the
  selection project. The default is 0.2

- plot:

  logical:TRUE or FALSE. To indicate if histograms or triangle should be
  plotted

- layout:

  character. It takes "portrait" or "landscape". If extent is set to
  "auto" the middle part of the image will be used so that NAs can be
  avoided. You indicate whether the automatic selection should indicate
  west-east (landscape) or north-south direction (portrait)

- draw:

  interactive. The method of digitisation if extent is set "auto". It
  takes "poly" or "rect"

- folder:

  An original directory of the images that contains all landsat 5,7 or 8
  bands and metadata. At the moment only Landsat folder is supported

- welev:

  Weather station elevation in meters

- clip:

  extent object or raster object or polygon from which an Extent object
  can be extracted. A polygon or raster will be reprojected to conform
  the data to be cropped. for example clip can take the form of c(xmin,
  xmax, ymin, ymax)

- iter.max:

  maximum iterations of sensible heat calculation.

- seed:

  optional single number. Passed to
  [`set.seed`](https://rdrr.io/r/base/Random.html) before the k-means
  clustering so the pixel selection is reproducible. The default `NULL`
  leaves the random-number stream untouched.

## Value

- Ts::

  The selected or input Ts

- Tshot::

  The Temperature of hot pixel

- TX::

  The average temperature of the hot pixels

- x::

  The x coordinate of the hot pixel

- y::

  The y coordinate of the hot pixel

- xyhot: :

  combined x and y coordinates

- candidates::

  The similar candidates' pixels that can be used

## Details

The function first divides the selected area into a number of clusters.
The mean and standard deviations are computed. The cluster means with
the lowest standard deviations and qualifying for upper (Ts) and lower
(NDVI) limits (see above) are selected. The provided upper and lower
probabilities are also used to select the final Ts.

## References

Allen, R. G., Burnett, B., Kramber, W., Huntington, J., Kjaersgaard, J.,
Kilic, A., Kelly, C., & Trezza, R. 2013. Automated Calibration of the
METRIC-Landsat Evapotranspiration Process. JAWRA Journal of the American
Water Resources Association, 49(3): 563-576.

## See also

[`coldTs`](https://gowusu.github.io/sebkc/reference/coldTs.md)

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
 folder=system.file("extdata","stack",package="sebkc")
  modhot=hotTs(folder=folder,welev=170,extent="auto",cluster=3)
  
   #use object of  landsat578
folder=system.file("extdata","stack",package="sebkc")
data=landsat578(data=folder, welev=170)
modfhot=hotTs(folder=data,welev=170,extent="auto",cluster=3)

#Or Ts = data
modhotTs=hotTs(data,welev=170,extent="auto",cluster=3)
albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
Ts=raster(system.file("extdata","Ts.grd",package="sebkc"))
NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
#simulatte with default parameters without albedo
modhot=hotTs(Ts=Ts,NDVI=NDVI,extent="auto",cluster=3)
#simulate with albedo but set extent to "auto" and cluster to 7
modhot=hotTs(Ts=Ts,NDVI=NDVI,albedo=albedo,cluster=7,
extent="auto")
} # }
```
