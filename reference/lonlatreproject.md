# Reproject WGS 84 longitude and latitude into a new raster coordinates and vice versa

Reproject WGS 84 longitude and latitude into a new raster coordinates
and vice versa

## Usage

``` r
lonlatreproject(x, y = NULL, var = NA, map, epsg = 4326)
```

## Arguments

- x:

  a longitude or raster or vector object

- y:

  a latitude or raster or vector object

- var:

  numeric variable. The length should be same as the length of latitude
  and longitude

- map:

  Raster\* object that will be used to define the new coordinates.

- epsg:

  coordinate reference number. The default is WGS 84: 4326

## Value

New x, y and map

## Details

If x is a list of longitudes then new horizontal coordinates based on
map coordinates. If x is a raster object with horizontal coordinates
then it will be projected to WGS 84.

## Author

George Owusu

## Examples

``` r
folder=system.file("extdata","stack",package="sebkc")
data=landsat578(data=folder, welev=362)
#> Computing Radiance: this may take several minutes........
#> Computing Reflectance
#> Computing Albedo
#> Computing NDVI
#> Computing SAVI
#> Computing LAI
#> Computing emissivities
#> Computing corrected thermal radiance
#> Computing Surface temperature
temp=data$Ts
latitude=seq(7.544,7.590,0.001)
longitude=seq(-1.211,-1.187,0.001)
var=seq(20,30,0.2)
latitude=latitude[1:length(longitude)]
var=var[1:length(longitude)]
newlonlat=lonlatreproject(longitude,latitude,var=var,map=temp)
#> Warning: GDAL Message 1: +init=epsg:XXXX syntax is deprecated. It might return a CRS with a non-EPSG compliant axis order. Further messages of this type will be suppressed.
longitude=newlonlat$longitude
map=newlonlat$map

```
