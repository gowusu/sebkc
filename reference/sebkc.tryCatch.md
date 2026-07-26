# Internal function to Catch errors and warnings

Internal function to Catch errors and warnings

## Usage

``` r
sebkc.tryCatch(expr)
```

## Arguments

- expr:

  an R expression to evaluate

## Value

value a list with 'value' and 'warning', where value may be an error
caught

## Examples

``` r
thiserror=str( sebkc.tryCatch( log( "a" ) ) )$value
#> List of 2
#>  $ value  :List of 2
#>   ..$ message: chr "non-numeric argument to mathematical function"
#>   ..$ call   : language log("a")
#>   ..- attr(*, "class")= chr [1:3] "simpleError" "error" "condition"
#>  $ warning: NULL
```
