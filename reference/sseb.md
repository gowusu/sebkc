# Simplified Surface Energy Balance

Senay etal (2011) Simplified Surface Energy Balance (SSEB) approach for
estimating landscape ET

## Usage

``` r
sseb(
  Ts,
  welev = NULL,
  TH = "auto",
  TC = "auto",
  sunelev = NULL,
  NDVI = NULL,
  DEM = NULL,
  albedo = NULL,
  NDVImax = 0.7,
  cc = 0.65,
  KL = 0.0065,
  ETo = NULL,
  wmo = NULL,
  airport = NULL,
  Krs = 0.16,
  surface = "grass",
  x = 1,
  Tmax = NULL,
  Tmin = NULL,
  zx = NULL,
  u = 2,
  DOY = NULL,
  latitude = NULL,
  n = NULL,
  RHmax = NULL,
  RHmin = NULL,
  clip = NULL,
  folder = NULL
)

# Default S3 method
sseb(
  Ts = NULL,
  welev = NULL,
  TH = "auto",
  TC = "auto",
  sunelev = NULL,
  NDVI = NULL,
  DEM = NULL,
  albedo = NULL,
  NDVImax = 0.7,
  cc = 0.65,
  KL = 0.0065,
  ETo = NULL,
  wmo = NULL,
  airport = NULL,
  Krs = 0.16,
  surface = "grass",
  x = 1.2,
  Tmax = NULL,
  Tmin = NULL,
  zx = NULL,
  u = 2,
  DOY = NULL,
  latitude = NULL,
  n = NULL,
  RHmax = NULL,
  RHmin = NULL,
  clip = NULL,
  folder = NULL
)
```

## Arguments

- Ts:

  A RasterLayer data that indicates radiometric surface temperature
  values from remote sensing image, preferably in Kelvin (K). You can
  also point to raster file on your computer.

- welev:

  Weather station elevation in meters

- TH:

  the average of the representative 3 hot pixels selected for hot "bare"
  areas. If it is set to "auto" TH will be estimated from
  [`hotTs`](https://gowusu.github.io/sebkc/reference/hotTs.md)

- TC:

  the average of representative 3 cold pixels selected from the
  irrigated fields. When it is set to "auto" TC will be estimated from
  [`coldTs`](https://gowusu.github.io/sebkc/reference/coldTs.md)

- sunelev:

  Angle of Sun elevation in degrees, as found in the meta data of the
  satellite image

- NDVI:

  A RasterLayer data that indicates NDVI (Normalized Difference
  Vegetation Index) values. You can also point to a file on your
  computer.

- DEM:

  A digital elevation model\[m\]

- albedo:

  A RasterLayer data that has albedo values. You can also provide file
  path or location on your computer

- NDVImax:

  numeric, NDVI value for heavily green-vegetated area in the image

- cc:

  numeric, the correction coefficient if the NDVI approaches 0.0

- KL:

  numeric, the assumed lapse rate of air moving along the terrain

- ETo:

  numeric, a standardized clipped grass reference ET. when it is set to
  NULL, it is automatically estimated from
  [`ETo`](https://gowusu.github.io/sebkc/reference/ETo.md)

- wmo:

  numeric. World Meteorological Organization weather station code. You
  can use the Worldwide Station List at <https://www.wunderground.com/>
  \[Accessed on 2016-5-6\] or <https://www.wetterzentrale.de/>
  \[Accessed on 2016-5-6\]

- airport:

  numeric. IATA or ICAO code. They are alphabetically listed at
  <https://en.wikipedia.org/wiki/List_of_airports_by_IATA_code:_A>
  \[Accessed on 2016-5-6\]

- Krs:

  Relationship between the fraction of extraterrestrial radiation that
  reaches the earth's surface, Rs/Ra, and the air temperature difference
  Tmax - Tmin for interior (Krs = 0.16) and coastal (Krs = 0.19) regions

- surface:

  character. It is either "grass" or "alfalfa"

- x:

  multiplier needed to estimate ET for tall, full cover crops such as
  alfalfa, corn and wheat. Default is 1.2.

- Tmax:

  Numeric. Maximum air Temperature in degree Celsius

- Tmin:

  Numeric. Minimum air Temperature in degree Celsius

- zx:

  The height above the weather station where the wind speed is measured
  \[m\]

- u:

  The satellites overpass time wind speed at the weather station \[m/s\]

- DOY:

  Numeric or Date \[YYYY-mm-dd\]. Day of the Year. If you give data in
  the form of date \[YYYY-mm-dd\], it will be converted to DOY

- latitude:

  geographical coordinates in decimal degrees. It should be negative for
  southern hemisphere

- n:

  Sunshine hours

- RHmax:

  Maximum Relative Humidity in percent

- RHmin:

  Minimum Relative Humidity in percent

- clip:

  extent object or raster object or polygon from which an Extent object
  can be extracted. A polygon or raster will be reprojected to conform
  the data to be cropped. for example clip can take the form of c(xmin,
  xmax, ymin, ymax)

- folder:

  An original directory of the images that contains all landsat 5,7 or 8
  bands and metadata. At the moment only Landsat folder is supported

## Value

- ET24::

  Actual 24 hour ET

- EF::

  ET fraction

- ETo::

  24 hour reference ETo

## References

Senay, G. B., Budde, M. E., & Verdin, J. P. 2011. Enhancing the
Simplified Surface Energy Balance (SSEB) approach for estimating
landscape ET: Validation with the METRIC model. Agricultural Water
Management, 98(4): 606-618.

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
Ts=raster(system.file("extdata","Ts.grd",package="sebkc"))
NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
#minimal data
modsseb=sseb(Ts=Ts,TH="full",TC="full",sunelev=50,NDVI=NDVI,DEM=NULL,
albedo=albedo,NDVImax=0.7,cc=0.65,KL=0.0065,
ETo="auto",x=1.2,Tmax=31,Tmin=28,zx=10,
u=2,DOY=37,latitude=5.6,n=NULL,RHmax=NULL,RHmin=NULL)

#uaing the folder parameters
folder=system.file("extdata","stack",package="sebkc")
ssebauto=sseb(folder=folder,welev=317,Tmax=31,Tmin=28, latitude=5.6)

#Another interactive example
#get Ts, Albedo, NDVI, SAVI, sunelev, DOY from landsat data 
welev=278
data=landsat578(rawdata,welev=welev)
#perform semi-auto simulation 
#Determine xyhot. Digitize polygon on the Ts map
modhot=hotTs(data,welev=300,extent="auto",cluster=2)
#determine the cold. Digitize polygon on the Ts map
modcold=coldTs(data,welev=275,extent="auto",cluster=2)
xyhot=modhot$Tshot
xycold=modcold$Tscold
#use object of \code{\link{coldTs}} or \code{\link{hotTs}} in different ways
modsseb=sseb(Ts=Ts,TH=xyhot,TC=modcold,sunelev=50,NDVI=NDVI,DEM=NULL,
albedo=albedo,NDVImax=0.7,cc=0.65,KL=0.0065,
ETo="auto",x=1.2,Tmax=31,Tmin=28,zx=10,
u=2,DOY=37,latitude=5.6,n=NULL,
RHmax=NULL,RHmin=NULL)
} # }
```
