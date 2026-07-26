# Internal function for doing interpolation

Internal function for doing interpolation

## Usage

``` r
mapmatrix(xydata, i, map)
```

## Arguments

- xydata:

  x and y data

- i:

  index

- map:

  geographical map

## Value

map

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
folder=system.file("extdata","stack",package="sebkc")
sebiauto=sebi(folder=folder,welev=317,Tmax=31,Tmin=28)
points=sampleRandom(sebiauto$EF,100,sp=TRUE)
pt=cbind(points@coords,points@data)
longitude=pt$x
latitude=pt$y
value=pt$layer
mapmatrix=mapmatrix(xydata,2)
map=invdist(longitude,latitude,var=value,ext=sebiauto$EF)
} # }
```
