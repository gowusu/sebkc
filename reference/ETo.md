# 24 Hour FAO 56 and ASCE-EWRI Reference Evapotranspiration ETo

24 Hour FAO 56 and ASCE-EWRI Reference Evapotranspiration ETo

## Usage

``` r
ETo(
  Tmax,
  Tmin,
  DOY,
  latitude,
  longitude = NULL,
  z = NULL,
  uz = 2,
  altitude = NULL,
  RHmax = NULL,
  RHmin = NULL,
  n = NULL,
  Krs = 0.16,
  albedo = 0.23,
  as = 0.25,
  bs = 0.5,
  Rs = NULL,
  Rn = NULL,
  G = 0,
  EF = NULL,
  wmo = NULL,
  airport = NULL,
  scale = "sine function",
  Kc = NULL,
  surface = "grass"
)

# Default S3 method
ETo(
  Tmax,
  Tmin,
  DOY,
  latitude,
  longitude = NULL,
  z = NULL,
  uz = 2,
  altitude = NULL,
  RHmax = NULL,
  RHmin = NULL,
  n = NULL,
  Krs = 0.16,
  albedo = 0.23,
  as = 0.25,
  bs = 0.5,
  Rs = NULL,
  Rn = NULL,
  G = 0,
  EF = NULL,
  wmo = NULL,
  airport = NULL,
  scale = "evaporative fraction",
  Kc = NULL,
  surface = "grass"
)
```

## Arguments

- Tmax:

  Numeric. Maximum air Temperature in degree Celsius

- Tmin:

  Numeric. Minimum air Temperature in degree Celsius

- DOY:

  Numeric or Date \[YYYY-mm-dd\]. Day of the Year. If you give data in
  the form of date \[YYYY-mm-dd\], it will be converted to DOY

- latitude:

  geographical coordinates in decimal degrees. It should be negative for
  southern hemisphere

- longitude:

  The longitude of the measurement site i.e. geographical coordinates in
  decimal degrees for the weather station. It should be negative for
  West and positive for East.

- z:

  The height above surface where wind speed was measured, in metres.

- uz:

  The wind speed in m/s

- altitude:

  The elevation of the weather station, in metres

- RHmax:

  Maximum Relative Humidity in percent

- RHmin:

  Minimum Relative Humidity in percent

- n:

  Sunshine hours

- Krs:

  Relationship between the fraction of extraterrestrial radiation that
  reaches the earth's surface, Rs/Ra, and the air temperature difference
  Tmax - Tmin for interior (Krs = 0.16) and coastal (Krs = 0.19) regions

- albedo:

  albedo values. For grass it is 0.23\[dimensionless\]

- as:

  Regression constant, intercept, expressing the fraction of
  extra-terrestrial radiation reaching the earth on overcast days (n =
  0); either a calibrated value or as = 0.25 can be used.

- bs:

  regression constant, slope, expressing the fraction of
  extraterrestrial radiation reaching the earth on overcast days (n =
  0); either a calibrated value or bs = 0.5 can be used.

- Rs:

  shortwave incoming radiation \[MJ/m2/day\]

- Rn:

  Net daily radiation in MJ/m2/day. If it is NULL, it will be computed

- G:

  Daily Soil heat flux in W/m2

- EF:

  Spatial. ET fraction \[-\] from
  [sebal](https://gowusu.github.io/sebkc/reference/sebal.md) or
  [tseb](https://gowusu.github.io/sebkc/reference/tseb.md)
  [sebs](https://gowusu.github.io/sebkc/reference/sebs.md). It is
  advisable to link it to any of the two functions.

- wmo:

  numeric. World Meteorological Organization weather station code. You
  can use the Worldwide Station List at <https://www.wunderground.com/>
  \[Accessed on 2016-5-6\] or <https://www.wetterzentrale.de/>
  \[Accessed on 2016-5-6\]

- airport:

  numeric. IATA or ICAO code. They are alphabetically listed at
  <https://en.wikipedia.org/wiki/List_of_airports_by_IATA_code:_A>
  \[Accessed on 2016-5-6\]

- scale:

  character The type of scale models to upscale instantaneous ET to
  daily ET. it takes "Evaporative Fraction" (Bastiaanssen et al 1998);
  "Sine Function" (Jackson et al 1983); "Solar Radiation" (French et al
  2005). Remember the blank spaces between the parameters. You can also
  respectively use "evaporative" or "solar" or "sine".

- Kc:

  Crop coefficient

- surface:

  character. It is either "grass" or "alfalfa"

## Value

- ETo::

  24 hour reference ETo \[mm/day\]

- ETa::

  24 hour actual ET. It is not NULL if EF has a numeric value

- ETc::

  24 hour crop ET. It is not NULL if Kc has a numeric value

- Rn::

  Net daily Radiation \[W/m2\]

- Rn_mj_m2_day::

  Net daily Radiation \[MJ/m2/day\]

- y::

  psychrometric constant \[kPa/degree Celsius\]

- slope::

  slope of saturation vapour pressure curve \[kPa/degree Celsius\]

- vpd::

  Vapour pressure deficit\[kPa\]

## References

ALLEN, R. G., PEREIRA, L. S., RAES, D., & SMITH, M. 1996. Crop
Evapotranspiration (guidelines for computing crop water requirements)
FAO Irrigation and Drainage Paper No. 56: FAO.

## See also

[`ETohr`](https://gowusu.github.io/sebkc/reference/ETohr.md)

## Author

George Owusu

## Examples

``` r
# \donttest{
#with minimal data
PET=ETo(Tmax=31,Tmin=28,latitude=5.6,DOY="6-2-2001")
#> Computing ETr24 or ETo24
PET$ETo
#> [1] 2.373057
#with enough data as in FAO56 EXAMPLE 18 on 6 July in Uccle (Brussels, Belgium) 
#located at 50 degrees 48 minutes N and at 100 m above sea level
PET24=ETo(Tmax=21.5,Tmin=12.3,z=10,uz=2.78,altitude=100,
RHmax=84,RHmin=63,n=9.25,DOY=187,latitude=50.8)
#> Computing ETr24 or ETo24
PET24$ETo
#> [1] 3.880216
# }
```
