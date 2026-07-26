# Day length calculation

Day length calculation

## Usage

``` r
dl(latitude, DOY, model = "CBM", Tmax = NULL, Tmin = NULL, p = 0.5)
```

## Arguments

- latitude:

  geographical coordinates in decimal degrees. It should be negative for
  southern hemisphere

- DOY:

  Numeric or Date \[YYYY-mm-dd\]. Day of the Year. If you give data in
  the form of date \[YYYY-mm-dd\], it will be converted to DOY

- model:

  character Type of model:"CBM" (Schoolfield, 1982),"BGC" (Running &
  Coughlan, 1988), "CERES" (Ritchie, 1991) or "FAO56"

- Tmax:

  Numeric. Maximum air Temperature in degree Celsius

- Tmin:

  Numeric. Minimum air Temperature in degree Celsius

- p:

  numeric. CMB parameter

## Value

DL

## References

- Schoolfield R. (1982). Expressing daylength as a function of latitude
  and Julian date.

- Running S. W. and Coughlan J. C. (1988). A general model of forest
  ecosystem processes for regional applications I. Hydrologic balance,
  canopy gas exchange and primary production processes. Ecological
  Modelling, 42(2), 125-154
  http://dx.doi.org/10.1016/0304-3800(88)90112-3.

- Ritchie J. T. (1991). Wheat Phasic Development. Modeling Plant and
  Soil Systems. Agronomy Monograph, 31.

- Allen R. G., Pereira L. S., Raes D. and Smith M. (1998). Crop
  evapotranspiration: Guidelines for computing crop water requirements.
  FAO Irrigation and Drainage Paper, 56, 300.

## Author

George Owusu

## Examples

``` r
DOY="2001-8-1"
latitude=0
model="CBM"
Tmax=31
Tmin=26
mod=dl(latitude,DOY,model,Tmax,Tmin=Tmin)
#> Computing ETr24 or ETo24
```
