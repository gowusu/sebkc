# Inverse distance interpolation

This function performs inverse distance interpolation based on
projection and geometry of input spatial data 'ext'.

## Usage

``` r
invdist(longitude, latitude, var, ext, idp = 0.5)
```

## Arguments

- longitude:

  The longitude of the measurement site i.e. geographical coordinates in
  decimal degrees for the weather station. It should be negative for
  West and positive for East.

- latitude:

  geographical coordinates in decimal degrees. It should be negative for
  southern hemisphere

- var:

  numeric variable. The length should be same as the length of latitude
  and longitude

- ext:

  RasterLayer object.

- idp:

  Neighbourhood parameter

## Value

var: interpolated values

## Author

George Owusu

## Examples

``` r
 if (FALSE) { # \dontrun{
#Get ext data
folder=system.file("extdata","stack",package="sebkc")
data=landsat578(data=folder, welev=362)
temp=data$Ts
latitude=seq(7.544,7.590,0.001)
longitude=seq(-1.211,-1.187,0.001)
var=seq(20,30,0.2)
latitude=latitude[1:length(longitude)]
var=var[1:length(longitude)]
map=invdist(longitude,latitude,var,ext=temp)
} # }
```
