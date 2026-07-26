# Biomass and yield estimation of crops and plants

This function uses Wageningen procedure (De Wit, 1968; Slabbers et al.,
1979) to compute the point or spatial possible maximum yield of a
standard crop based on a Rs data. This function can also perform
internal interpolation of data if longitude and latitude are provided.
See datails and examples below.

## Usage

``` r
biomass(
  Tday = NULL,
  T24,
  Rs,
  latitude,
  longitude = NULL,
  HI = NULL,
  plantdate,
  harvestdate,
  c30 = 0.0108,
  adaptability = 3,
  LAI = 5,
  crop = NULL
)
```

## Arguments

- Tday:

  numeric. Average Day time temperature in degrees Celsius for the
  growing period. A Remote sensing surface Temperature can also be used

- T24:

  numeric. Average 24-hour (including night hours) temperature in
  degrees Celsius for the growing period.

- Rs:

  A shortwave radiation data \[ca^-2 day^-1\]. RasterLayer data can be
  supplied.

- latitude:

  geographical coordinates in decimal degrees. It should be negative for
  southern hemisphere

- longitude:

  The longitude of the measurement site i.e. geographical coordinates in
  decimal degrees for the weather station. It should be negative for
  West and positive for East.

- HI:

  numeric. Harvest Index. It can be list or spatial

- plantdate:

  character. A date of planting in the form "YYYY-mm-dd"

- harvestdate:

  character. A date of harvest in the form "YYYY-mm-dd"

- c30:

  numeric. A maintenance respiration coefficient (c30) at 30 oC. it
  takes either a value of 0.0108 for non-legume crop or 0.0283 for
  legumes

- adaptability:

  numeric. The climate adaptability group number. It takes values of
  1,2,3,or 4. it depends on the temperature profile of study area. See
  references (eg. FAO (1981)) for more explanation.

- LAI:

  Leaf Area Index, dimensionless

- crop:

  character. A crop type. For example "maize"

## Value

yield and biomass in t/ha

## Details

- **Point Estimation**:

  `biomass` can compute biomass and yield for one location. Each input
  data must be given only one value

- **Point Estimation for several locations** :

  The function can compute biomass and yields for more than one
  location. Input data must be a list or array.

- **Spatial Estimation**:

  Input data such as Rs, Tday, T24, HI or LAI can be spatial.

- **Spatial Estimation with internal interpolation**:

  It takes Tday, for example from satellite image, and do interpolation
  for the rest of the data. Both longitude, latitude, and the variables
  must be of the same length. If longitude is not provided no
  interpolation will be done.

@references

- De Wit, C. T. (1968). Plant production Miscellaneous Papers: Landbouw
  Hogeschool, Wageningen.

- Slabbers, P. J., Herrendorf, V. S., & Stapper, M. (1979). Evaluation
  of simplified water-crop yield models. Agricultural Water Management,
  2(2), 95-129. doi: http://dx.doi.org/10.1016/0378-3774(79)90026-X

- FAO (1981). Report on the Agro-ecological Zones Project: Food and
  Agriculture Organization of the United Nations.

## Author

George Owusu

## Examples

``` r
if (FALSE) { # \dontrun{
################# Point Estimation ###################
# Identify the climatic input data on 
latitude=11.18 #in decimal degrees
Rs= 452 #Average shortwave radiation over the growing period (cal-2day-1)
Tday = 26.9 #  Average day time difference (oC) over the growing period
T24 =25.3 #Average 24 hour mean temperature over the growing period
# Gather crop Information
crop="maize"
plantdate="2015-05-1" #Start date
harvestdate="2015-08-29" #end date 120 days
HI=0.4 #Harvest index for maize
adaptability=3 #crop climate adaptability group III
LAI = 4 #Leaf Area Index
pointmod=biomass(Tday,T24,Rs,latitude,longitude=NULL,HI,
plantdate,harvestdate,c30=0.0108,adaptability=3,LAI=5,crop=NULL)

################# Point Estimation for several locations#################
latitude=c(10,11.18,12,9)
longitude=c(-1.4,-2,-1,-2.6)
Rs= c(450,452,400,500)
Tday = c(25,26.9,24,20,30)
T24=c(24,25.3,23,28)
HI=c(0.4,0.4,0.5,0.25)
plantdate=c("1-03-2015","1-04-2015","1-05-2015")
harvestdate=c("29-08-2015","29-08-2015","29-08-2015")

# computing several points
pointSmod=biomass(Tday,T24,Rs,latitude,longitude=longitude,HI,
plantdate,harvestdate,c30=0.0108,adaptability=3,LAI=5,crop=NULL)
yield=pointSmod$yield #tabular data
biomass=pointSmod$output #tabular data
#access on biomass
pointSmod$output$yield
} # }
if (FALSE) { # \dontrun{
############## Spatial Estimation with internal interpolation ###########################
#generate spatial data
folder=system.file("extdata","stack",package="sebkc")
modauto=landsat578(data=folder, welev=362)
Tday=modauto$Ts-273.15 #temperature in degree celcius
LAI=modauto$LAI
latitude=seq(7.544,7.590,0.001) #for interpolation
longitude=seq(-1.211,-1.187,0.001) #for interpolation
latitude=latitude[1:length(longitude)] #to make sure they are of the same length
spatialmodint=biomass(Tday,T24,Rs,latitude,longitude=longitude,HI,
plantdate,harvestdate,c30=0.0108,adaptability=3,LAI=LAI,crop=NULL)
# Rs, T24, and HI are interpolated
plot(spatialmodint$yield)
############## Spatial Estimation  ###########################
Rs= 452
T24 =25.3
HI=0.4 
spatialmod=biomass(Tday,T24,Rs,latitude,longitude=longitude,HI,
plantdate,harvestdate,c30=0.0108,adaptability=3,
LAI=LAI,crop=NULL)
plot(spatialmod$yield)
} # }
```
