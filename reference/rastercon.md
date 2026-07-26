# Spatial Conditional function on raster object

It computes equivalent of logical if function in R

## Usage

``` r
rastercon(condition, trueValue, falseValue)
```

## Arguments

- condition:

  logical

- trueValue:

  logical

- falseValue:

  logical

## Value

Return True or False

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
LAI=raster(system.file("extdata","LAI.grd",package="sebkc"))
  eNB=rastercon( NDVI<0 & albedo<0.47,0.99,rastercon(LAI>=3,0.98,0.97+(LAI*0.0033)))
  } # }
```
