# Reprojecting and extracting Raster values

The function reprojects longitudes and latitudes into map coordinates
and extract raster values

## Usage

``` r
xyextract(map, longitude, latitude)
```

## Arguments

- map:

  a raster object.

- longitude:

  The longitude of the measurement site i.e. geographical coordinates in
  decimal degrees for the weather station. It should be negative for
  West and positive for East.

- latitude:

  geographical coordinates in decimal degrees. It should be negative for
  southern hemisphere

## Value

data

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
lonlat=read.table(system.file("extdata","sys","xyvalues.txt",package="sebkc"),header=TRUE)
longitude=lonlat$longitude
latitude=lonlat$latitude
folder=system.file("extdata","stack",package="sebkc")
data=landsat578(data=folder, welev=362)
Ts=data$Ts
Tsdata=xyextract(Ts,longitude,latitude)

} # }
```
