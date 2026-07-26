# Reprojecting and cropping images

This function returns the subset of x object as defined by y extent
object. See
[`crop`](https://rspatial.github.io/terra/reference/crop.html) and
[`extent`](https://rdrr.io/pkg/raster/man/extent.html) for more details.
The difference between
[`crop`](https://rspatial.github.io/terra/reference/crop.html) and
cropsebkc is the ability for the latter to reproject y object with x
object projection.

## Usage

``` r
cropsebkc(x, y)
```

## Arguments

- x:

  object to be extracted

- y:

  extent object or raster object or polygon from which an Extent object
  can be extracted. A polygon or raster will be reprojected to conform x

## Value

cropped x and reprojected y

## See also

[`crop`](https://rspatial.github.io/terra/reference/crop.html)

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
folder=system.file("extdata","stack",package="sebkc")
stack=sebkcstack(folder=folder)
data=cropsebkc(stack$data,stack$data)
} # }
```
