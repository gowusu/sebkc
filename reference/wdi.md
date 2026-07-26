# Water Deficit Index (WDI)

Moran et. al (1994) Estimation of crop water deficit using surface-air
temperature and spectral vegetation index.

## Usage

``` r
wdi(
  Ts = NULL,
  Ta,
  NDVI,
  albedo = NULL,
  Tsmax = "auto",
  Tsmin = "auto",
  sunelev = NULL,
  folder = NULL,
  welev = NULL,
  Rn24 = NULL
)

wdi(
  Ts = NULL,
  Ta,
  NDVI,
  albedo = NULL,
  Tsmax = "auto",
  Tsmin = "auto",
  sunelev = NULL,
  folder = NULL,
  welev = NULL,
  Rn24 = NULL
)
```

## Arguments

- Ts:

  A RasterLayer data that indicates radiometric surface temperature
  values from remote sensing image, preferably in Kelvin (K). You can
  also point to raster file on your computer.

- Ta:

  numeric. Air temperature. If not available you can use

- NDVI:

  A RasterLayer data that indicates NDVI (Normalized Difference
  Vegetation Index) values. You can also point to a file on your
  computer.

- albedo:

  A RasterLayer data that has albedo values. You can also provide file
  path or location on your computer

- Tsmax:

  numeric. Maximum sensor surface temperature \[K\]

- Tsmin:

  numeric. Minimum sensor surface temperature \[K\]
  [`sebal`](https://gowusu.github.io/sebkc/reference/sebal.md) to get
  pixel based Ta.

- sunelev:

  Angle of Sun elevation in degrees, as found in the meta data of the
  satellite image

- folder:

  An original directory of the images that contains all landsat 5,7 or 8
  bands and metadata. At the moment only Landsat folder is supported

- welev:

  Weather station elevation in meters

- Rn24:

  A 24 hour net radiation \[W/m2\]. It is needed to estimate daily ET.
  It can be computed with
  [`ETo`](https://gowusu.github.io/sebkc/reference/ETo.md) with even a
  minimal data of Tmax and Tmin. You can also point to raster file on
  your computer

## Value

- EF::

  standardised wdi excluding negatives and one plus

- index::

  raw wdi including negatives and more than one

## Details

Note that surface temperature (Tsmin and Tsmax) are different from air
temperature (Tmin and Tmax) in
[`ETo`](https://gowusu.github.io/sebkc/reference/ETo.md),
[`sseb`](https://gowusu.github.io/sebkc/reference/sseb.md),
[`sebs`](https://gowusu.github.io/sebkc/reference/sebs.md)

You can get Ta from
[`sebal`](https://gowusu.github.io/sebkc/reference/sebal.md) Ta

## References

Moran, M. S., Clarke, T. R., Inoue, Y., & Vidal, A. 1994. Estimating
crop water deficit using the relation between surface-air temperature
and spectral vegetation index. Remote Sensing of Environment, 49(3):
246-263.

## See also

[`sseb`](https://gowusu.github.io/sebkc/reference/sseb.md) and
[`sebi`](https://gowusu.github.io/sebkc/reference/sebi.md)

## Author

George Owusu

## Examples

``` r
 if (FALSE) { # \dontrun{
albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
Ts=raster(system.file("extdata","Ts.grd",package="sebkc"))
NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
modWDI=wdi(Ts=Ts,Ta=299,NDVI=NDVI,albedo=NULL,Tsmax="auto",
Tsmin="auto",sunelev=50)
plot(modWDI$EF)

#Using folder parameter
folder=system.file("extdata","stack",package="sebkc")
wdiauto=wdi(folder=folder,welev=317,Tsmax=31,Tsmin=28,Ta=290)

#' #Another interactive example
#get Ts, Albedo, NDVI, SAVI, sunelev, DOY from landsat data 
welev=278
data=landsat578(rawdata,welev=welev)
#perform semi-auto simulation 
#Determine xyhot. Digitize polygon on the Ts map
modhot=hotTs(data,welev=300,extent="auto",cluster=2)
#determine the cold. Digitize polygon on the Ts map
modcold=coldTs(data,welev=275,extent="auto",cluster=2)
Tsmax=modhot$Tshot
Tsmin=modcold$Tscold
#use object of \code{\link{coldTs}} or \code{\link{hotTs}} in different ways
modWDI2=wdi(data,Ta=299,Tsmax=Tsmax,Tsmin=modcold,sunelev=50)
plot(modWDI2$EF)
} # }
```
