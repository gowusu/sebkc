# Optional Automatic computation of cold pixel Radiometric Surface Temperature

This function computes temperature of a cold pixel based on NDVI,
Surface Temperature, DEM and albedo. The main difference with
[`coldTs`](https://gowusu.github.io/sebkc/reference/coldTs.md) is that
preference is given to temperature selection before considering NDVI and
albedo selection. Moreover, this is a simple end member selection where
no clustering is done.

## Usage

``` r
coldTs2(
  Ts,
  NDVI,
  albedo = NULL,
  sunangle = NULL,
  DEM = NULL,
  cluster = 10,
  extent = "interactive",
  upper = 0.2,
  lower = 0.1,
  NDVI.probs = 0.95,
  plot = TRUE,
  layout = "portrait",
  draw = "poly",
  folder = NULL,
  welev = NULL
)

# Default S3 method
coldTs2(
  Ts = NULL,
  NDVI,
  albedo = NULL,
  sunangle = NULL,
  DEM = NULL,
  cluster = 10,
  extent = "interactive",
  upper = 0.2,
  lower = 0.1,
  NDVI.probs = 0.6,
  plot = TRUE,
  layout = "portrait",
  draw = "poly",
  folder = NULL,
  welev = NULL
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

- sunangle:

  sun angle (degree) above the horizon.

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

  The quantile probability at which NDVI should be included in the
  selection project. It ranges from 0-1. The default is 0.95

- lower:

  The quantile probability at which Ts should be included in the
  selection project.The default is 0.2

- NDVI.probs:

  The probability of selecting quantile NDVI that is above the
  threshold. The value of 0.6 means only NDVI values that is above 60
  the analyses.

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

## Value

- Ts::

  The selected or input Ts \[K\]

- Tscold::

  The Temperature of cold pixel \[K\]

- TC::

  The average temperature of the cold pixels \[K\]

- x::

  The x coordinate of the cold pixel

- y::

  The y coordinate of the cold pixel

- xycold: :

  combined x and y coordinates

- candidates::

  The similar candidates' pixels that can be used

## Details

The function first divides the selected area into a number of clusters.
The mean and standard deviations are computed. The cluster means with
the lowest standard devaitions and qualifying for upper (NDVI) and lower
(Ts) limits (see above) are selected. The provided upper and lower
probabilities are also used to select the final Ts.

## References

Allen, R. G., Burnett, B., Kramber, W., Huntington, J., Kjaersgaard, J.,
Kilic, A., Kelly, C., & Trezza, R. 2013. Automated Calibration of the
METRIC-Landsat Evapotranspiration Process. JAWRA Journal of the American
Water Resources Association, 49(3): 563-576.

## See also

[`coldTs`](https://gowusu.github.io/sebkc/reference/coldTs.md),
[`hotTs`](https://gowusu.github.io/sebkc/reference/hotTs.md),
[`hotTs2`](https://gowusu.github.io/sebkc/reference/hotTs2.md)

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
#using landsat folder
folder=system.file("extdata","stack",package="sebkc")
modcold=coldTs2(folder=folder,welev=170,extent="auto")

#use object of  
folder=system.file("extdata","stack",package="sebkc")
data=landsat578(data=folder, welev=362)
modcold=coldTs2(folder=data,welev=170,extent="auto",cluster=3)

modcoldTs2=coldTs2(data,welev=170,extent="auto",cluster=3)

#using input different input data
albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
Ts=raster(system.file("extdata","Ts.grd",package="sebkc"))
NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
#simulate with default parameters without albedo
modcold=coldTs2(Ts=Ts,NDVI=NDVI,sunangle=50.7)
#simulate with albedo but set extent to "auto" and cluster to 7
modcold=coldTs2(Ts=Ts,NDVI=NDVI,albedo=albedo,sunangle=50,cluster=7,
extent="auto")
} # }
```
