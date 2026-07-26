# Writing sebkc object to a folder

This function writes spatial output of the following functions
[`sebal`](https://gowusu.github.io/sebkc/reference/sebal.md),
[`sebi`](https://gowusu.github.io/sebkc/reference/sebi.md),[`ssebi`](https://gowusu.github.io/sebkc/reference/ssebi.md),
[`wdi`](https://gowusu.github.io/sebkc/reference/wdi.md),[`sebs`](https://gowusu.github.io/sebkc/reference/sebs.md),[`sseb`](https://gowusu.github.io/sebkc/reference/sseb.md),
[`tseb`](https://gowusu.github.io/sebkc/reference/tseb.md) and
[`landsat578`](https://gowusu.github.io/sebkc/reference/landsat578.md)

## Usage

``` r
writesebkc(object, folder = NULL, xy = NULL, overwrite = TRUE)
```

## Arguments

- object:

  sebkc object

- folder:

  Folder the files should be written to.If it is set to NULL, it is
  written to the input folder of sebkc object.

- xy:

  A dataframe of of xy coordinates in in decimal degrees or meters in
  the order of c(x,y). If it is not set NULL, the corresponding values
  are extracted and written to the folder

- overwrite:

  logical whether the file should be over written

## Value

Writes output to a folder

## Examples

``` r
if (FALSE) { # \dontrun{
folder=system.file("extdata","stack",package="sebkc")
stack=landsat578(data=folder, welev=362)
writesebkc(stack)
sebaloutput=sebal(folder = folder,welev = 317)
writesebkc(sebaloutput)

} # }
```
