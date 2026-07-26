# Goodness of fit Test

Goodness of fit Test

## Usage

``` r
fitting(obs, est, P = 1)
```

## Arguments

- obs:

  Observed data

- est:

  Estimated data

- P:

  Number of parameters

## Value

- r2::

  R-square \[-\]

- r2_adj::

  An adjusted R-square \[-\]

- d::

  Index of agreement\[-\]

- p_value::

  p value

- sse::

  Sum of square error

- sst::

  Sum of square total

- AIC::

  AIC

- BIC::

  BIC

- RMSE::

  Root Mean Square Error

- RMSE_adj::

  adjusted Root Mean Square Error

- NRMSE:

  Normalised Root Mean Square Error

- NRMSE_adj::

  adjusted Normalised Root Mean Square Error

- E::

  Nash Sutcliffe Coefficient

- b::

  coefficients of regression

- AAE::

  average absolute error

- ARE::

  average relative error

- MBE::

  Mean Bias error

- MPE::

  Mean Percentage error

## Examples

``` r
file=system.file("extdata","sys","irrigation.txt",package="sebkc")
data=read.table(file,header=TRUE)  
obs=data$ETo
est=data$ETc
P=4
mod=fitting(obs,est,P)
r2=mod$r2
```
