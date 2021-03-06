% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kc.R
\name{kc}
\alias{fit.kc}
\alias{kc}
\alias{kc.default}
\alias{plot.kc}
\title{Estimating FAO56 Single and Dual Crop coefficients with a water balance model in R}
\usage{
kc(ETo, P, RHmin, soil = "Sandy Loam", crop, I = 0, CNII = 67, u2 = 2,
  FC = NULL, WP = NULL, h = NULL, Zr = NULL, Ze = 0.1, kc = NULL,
  p = NULL, lengths = NULL, fw = 0.8, fc = NULL, kctype = "dual",
  region = NULL, regexpr = "sweet", initgw = 50, LAI = 3, rc = NULL,
  Kb = 0, Rn = NULL, hmodel = "gaussian", Zrmodel = "linear",
  xydata = NULL, REW = NULL, a1 = 360, b1 = -0.17, a2 = 240,
  b2 = -0.27, a3 = -1.3, b3 = 6.6, a4 = 4.6, b4 = -0.65,
  report = "stages", time = "06:00", HWR = 1, latitude = NULL,
  density = "full", DOY = NULL, slope = NULL, y = NULL, r1 = 100,
  nongrowing = "weeds", theta = NULL, profile = "variable",
  Kcmin = 0.15, Ky = NULL)

\method{kc}{default}(ETo, P, RHmin, soil = "Sandy Loam", crop, I = 0,
  CNII = 67, u2 = 2, FC = NULL, WP = NULL, h = NULL, Zr = NULL,
  Ze = 0.1, kc = NULL, p = NULL, lengths = NULL, fw = 1, fc = NULL,
  kctype = "dual", region = NULL, regexpr = "sweet", initgw = 50,
  LAI = 3, rc = NULL, Kb = 0, Rn = NULL, hmodel = "gaussian",
  Zrmodel = "linear", xydata = NULL, REW = NULL, a1 = 360, b1 = -0.17,
  a2 = 240, b2 = -0.27, a3 = -1.3, b3 = 6.6, a4 = 4.6, b4 = -0.65,
  report = "stages", time = "06:00", HWR = 1, latitude = NULL,
  density = "full", DOY = NULL, slope = NULL, y = NULL, r1 = 100,
  nongrowing = "weeds", theta = NULL, profile = "variable",
  Kcmin = 0.15, Ky = NULL)

plot.kc(x, output = "AW_measured", main = NULL, cex = 1, legend = TRUE,
  pos = "top", horiz = FALSE, ...)

fit.kc(x, output = "AW_measured")
}
\arguments{
\item{ETo}{Numeric. Reference evapotranspiration [mm]. This can be estimated with \code{\link{ETo}}}

\item{P}{numeric. Precipitation [mm]}

\item{RHmin}{numeric. Minimum Relative Humidity [per cent]}

\item{soil}{character. Name of soil texture as written in FAO56 TABLE 19}

\item{crop}{character. Name of the crop as written on 
FAO56 Table 11, TABLE 12 and TABLE 17}

\item{I}{numeric. Irrigation [mm]}

\item{CNII}{numeric. The initial Curve Number of the study site. 
This can be calibrated with \code{cal.kc}.}

\item{u2}{numeric. Wind speed of the nearby weather station [m/s] measured at 2m}

\item{FC}{numeric. soil water content at field capacity [-][m3 m-3]. 
This can be calibrated with \code{cal.kc}.}

\item{WP}{numeric. soil water content at wilting point [-][m3 m-3]. 
This can be calibrated with \code{cal.kc}.}

\item{h}{numeric. Maximum crop height [m]. See details below.
This can be calibrated with \code{cal.kc}}

\item{Zr}{numeric. Rooting depth.  See details below.
This can be calibrated with \code{cal.kc}}

\item{Ze}{numeric. The depth of the surface soil layer that is subject to drying 
by way of evaporation [0.10-0.15 m]. This can be calibrated with \code{\link{cal.kc}}}

\item{kc}{list. The tabulated Kc values for single crop c(kcini,kcmid,kcend)
or  dual crop model c(kcbini,kcbmid,kcbend). To use with \code{\link{sebal}} or 
any of the Surface Energy Balance models, it should be linked as list(model1,model2,model3) or
 c(model1$EF,model2$EF,model3$EF) but not a mixture of both. 
This can be calibrated with \code{\link{cal.kc}}}

\item{p}{The depletion coefficient. This can be calibrated with \code{\link{cal.kc}}}

\item{lengths}{list. The lengths of the crop growth stages. It can take dates in the form of 
c(start of planting ('YYYY-mm-dd), start of rapid growth, 
start of mid season, start of maturity, harvesting) 
or in the form number of days: c(initial days, dev days, mid days, late days). 
See FAO56 TABLE 11. If performing time series simulation, different dates or 
one season dates can be set; for instance two seasons simulation can be of the form 
c(initial days, dev days, mid days, late days, initial days, dev days, mid days, late days) 
or just c(initial days, dev days, mid days, late days)}

\item{fw}{list or numeric or character. This can take different fw values
throughout the growing season. It can also take one numeric value and this will
be transform depending on the type of wetting. If the wetting event is only Precipitation,
it is set to 1; if it is only Irrigation it is set to the provided fw value. If it is
both Precipitation and Irrigation weights are assigned for calculation of fw. 
fw can also take a character such as "Drip", "Sprinkler", "Basin", "Border", 
"alternated furrows", "Trickle", "narrow bed" and "wide bed". 
In this case fw values from FAO56 are assigned.}

\item{fc}{numeric. fraction of ground covered by tree canopy}

\item{kctype}{character. The type of model either "single" or "double". In spatial modelling involving
Evaporation Fraction, kctype can be set to "METRIC" or "SSEB" so that 
Actual Evapotranspiration will be estimated as ETo*EF else it will be
estimated with Rn as ETa=(EF*Rn)/2.45. Note that Rn is in MJm-2 day-1.}

\item{region}{character. The region of the crop as on FAO56 TABLE 11}

\item{regexpr}{character. Extra term that can be used to search crop names on FAO56 Tables}

\item{initgw}{numeric. Initial ground water level [m]}

\item{LAI}{numeric. Leaf Area Index}

\item{rc}{numeric Runoff coefficient [0-1]. 
This must be used if there is no information on Curve Number.
This can be calibrated with \code{\link{cal.kc}}}

\item{Kb}{ground water (base flow) recession parameter. 
The default is zero where no baseflow occurs}

\item{Rn}{Net daily radiation [w/m2]. Estimated this with  \code{\link{ETo}}}

\item{hmodel}{character. The type of crop growth model during the simulation. 
It takes "gaussian","linear", "exponential"}

\item{Zrmodel}{character. The type of root growth model during the simulation. 
It takes "gaussian","linear", "exponential"}

\item{xydata}{list c(longitude, latitude). The observation points of the model. 
This should be set only if the model is spatial.}

\item{REW}{numeric Readily Evaporable Water. It can be calibrated with \code{\link{cal.kc}}}

\item{a1}{soil water storage to maximum root depth (Zr) at field capacity}

\item{b1}{Liu et al. (2006) parameter for computing capillary rise [-0.17]. It can calibrated with
\code{\link{cal.kc}}}

\item{a2}{Liu et al. (2006) parameter for computing capillary rise [240]. It can be calibrated with
\code{\link{cal.kc}}}

\item{b2}{Liu et al. (2006) parameter for computing capillary rise [-0.27]. It can be calibrated with
\code{\link{cal.kc}}}

\item{a3}{Liu et al. (2006) parameter for computing capillary rise [-1.3]. It can be calibrated with
\code{\link{cal.kc}}}

\item{b3}{Liu et al. (2006) parameter for computing capillary rise [6.6]. It can be calibrated with
\code{\link{cal.kc}}}

\item{a4}{Liu et al. (2006) parameter for computing capillary rise [4.6]. It can be calibrated with
\code{\link{cal.kc}}}

\item{b4}{Liu et al. (2006) parameter for computing capillary rise [-0.65]. It can be calibrated with
\code{\link{cal.kc}}}

\item{report}{list. In case of spatial modelling, which days' map should be reported. 
The default is "stages" where the periods of stages maps are reported as values. 
In an ascending order it can take values such as 1:12 days or c(1,3,6:12) days.}

\item{time}{character. Time of the occurrence Irrigation or Precipitation.}

\item{HWR}{numeric. Default is 1. Height to width ratio of individual plants or groups of 
plants when viewed from the east or from the west [-].}

\item{latitude}{geographical coordinates in decimal degrees. It should be negative
for southern hemisphere}

\item{density}{character. density factor for the trees, 
vines or shrubs on how they are arranged.
It can take "full", "row" or "random" or "round".}

\item{DOY}{Numeric or Date [YYYY-mm-dd]. Day of the Year. 
If you give data in the form of date [YYYY-mm-dd], it will be converted to DOY}

\item{slope}{numeric. slope of saturation vapour pressure curve at air temperature T [kPa oC-1].}

\item{y}{numeric. psychometric constant  [kPa oC-1]}

\item{r1}{leaf resistance FOR STOMATAL CONTROL [s/m]}

\item{nongrowing}{character. The type of surface during non-growing season. 
It can take "weeds", "bare", or "mulch"}

\item{theta}{numeric. The measured soil water content [m3/m3] in the 
profile. Leave those days without measured data blank.}

\item{profile}{numeric or character. The profile depth to perform water balance.
It takes "variable" when updating the depth with root growth (default). Or taking "fixed" when
the maximum root zone is used. It can also take numeric value and this can be longer than
maximum root length.}

\item{Kcmin}{Minimum  Kc value}

\item{Ky}{list or array. yield response factor [-] in the form of 
c(initial Ky,crop dev't Ky, mid season Ky, late season Ky, total season Ky). 
For instance Maize Ky is c(0.4,0.4,1.3,0.5,1.25) and it is available at
 \url{http://www.fao.org/nr/water/cropinfo_maize.html}.
Other crops Ky can be found at \url{http://www.fao.org/nr/water/cropinfo.html}
If only one Ky value is supplied it will be repeated for all the seasons.}

\item{x}{a return object of the function.}

\item{output}{The kind of plot that should be displayed. It takes 
'AW_measured' for plotting Available and modelled soil moisture;
'TS_measured' for Total Available and modelled soil moisture;
'Kc' for displaying Kcb and Kc
'Evaporation' for displaying of Evaporation and Transpiration
'balance' for displaying Irrigation, precipitation, ETa and critical water deficit}

\item{main}{Title of the plot}

\item{cex}{the relative font size of the texts. 
Vlaues: "horiz";    logical: if TRUE, set the legend horizontally rather than vertically}

\item{legend}{logical. TRUE or FALSE. Whether a legend should be displayed or not.}

\item{pos}{the position of legend. It takes one of the following 'top', 'topright',
'topleft', 'bottom', 'bottomright' 'bottomleft', 'center'. This must be set so that
the legend should not obscure the main plot.}
}
\value{
The output of the function depends of the kctype and kc data type.
For single and dual crop coefficient, an "output" dataframe contains variables such as
TAW,RAW,Runoff,Ks,Kcb,h,Kc,ETc,DP,CR,Drend,Kr,De,AW,E,Ke,CN. In case kc is set to EF the function returns
Kc,Kcxy,ETc,ETcxy,Drend,Drendxy,TAW,TAWxy,RAW,RAWxy,AW,Awcxy,Ineed,Ineedxy
 \itemize{
\item{output:} {  dataframe of the return values}
\item{ETc_adj:} { The soil water stress adjusted ETc  [mm]}
\item{Ya_by_Ym:} {actual yield by potential yield}
\item{yield_decrease:} {The relative yield decrease: 1-Ya/Ym}
\item{ETdecrease:} {The relative evapotranspiration deficit :1- ETc_adj/ETc}
\item{yield_by_ET_decrease:} {he relative yield decrease: 1-Ya/Ym and 
the relative evapotranspiration deficit in one dataframe}
\item{Awcxy:} { Available Soil water for each location at each time step [mm]}
\item{CR:} { capillary rise [mm/day]}
\item{GW:} { Updated ground water depth}
\item{CN:} { Curve Number}
\item{De:} { cumulative depth of evaporation (depletion) 
from the soil surface layer [mm] at the end of the day [mm/day]}
\item{DP:} { Deep Percolation [mm/day]}
\item{Drend:} { cumulative depth of evapotranspiration (depletion) 
from the root zone at the end of the day [mm]. 
This is equal to the amount of water needed to bring soil water to field capacity}
\item{Drendxy:} { cumulative depth of evapotranspiration (depletion) 
from the root zone at the end of the day [mm] for each xydata location at each time step. 
This is equal to the amount of water needed to bring soil water to field capacity}
\item{E:} { evaporation [mm/day]}
\item{ETc:} { crop evapotranspiration [mm/day]. 
If water stress coefficient (Ks) is less than 1, 
ETc is crop evapotranspiration under non-standard conditions else ETc is
crop evapotranspiration under standard conditions}
\item{ETcxy:} { crop evapotranspiration [mm/day] for each location at each time step.
crop evapotranspiration under standard conditions}
\item{h:} { crop height [m]}
\item{Ineed:} { Amount of soil water that is below maximum allowed depletion [mm]}
\item{Ineedxy:} { Amount of soil water that is below maximum allowed depletion 
for each xydata location at each time step [mm] }
\item{Kc:} { standard crop coefficient [-]}
\item{Kc_adj:} { adjusted crop coefficient [-]}
\item{EF:} { Evaporation Fraction [-]}
\item{EFxy:} { Evaporation Fraction for each monitored location [-]}
\item{kcb:} { basal crop coefficient for the growing season [-]}
\item{Kcxy:} { crop coefficient for each location at each time step [-]}
\item{Ke:} { soil evaporation coefficient}
\item{Kr: } { soil evaporation reduction coefficient [-]}
\item{Ks:} { Water stress coefficient [-]}
\item{RAW: } { Readily Available Water [mm]}
\item{RAwcxy:} { Readily Available Soil water for each xydata location at each time step [mm]}
\item{Runoff:} { Runoff [mm/day]}
\item{TAW:} { Total Available Water [mm] }
\item{TAwcxy:} { Total Available Soil water for each location at each time step [mm]}
}
}
\description{
This function computes FAO56 crop coefficients and water balance. 
The output of the model includes crop evapotranspiration, water stress coefficient,
irrigation requirements, soil available water, etc. In addition to point estimation,
spatial modelling is integrated using Evaporation fractions from 
\code{\link{sebal}}, \code{\link{sebi}}, \code{\link{ETo}},  \code{\link{sseb}},  
\code{\link{sebs}}. \code{cal.kc} can be used to calibrate  the model.
}
\details{
\describe{
 \item{\strong{Using Existing Tables}}{The function include backend database of FAO Tables.
 The user can use the parameter values such as Zr, h, p, FC, WP, etc from FAO56 tables. 
 Most of the parameters with NULL values can take default values from FAO56 Tables. When a crop name results in 3 crops,
  the user can differentiate  that with  regexpr parameter. For example, 
  a crop parameter with a name "maize" can match
  two crops on FAO Table 12. Including  regexpr parameter with "sweet" with match only one. 
  }
 \item{\strong{Modelling crop height and Rooting Depth}}{ 
 At each time step a crop height or root depth is estimated. 
 The crop  or rooting depth growth model can take  exponential model of the form
 hday=0+0.75*max(h)*(1-exp(-day/0.75*max(days))) or    Zrday=0+0.75*Zr*(1-exp(-day/0.75*max(days)))
  The gaussian model is of the form hday=0+0.75*max(h)*(1-exp(-day^2/(0.75*max(days)^2))) or Zrday=0+0.75*(1-exp(-day^2/(0.75*max(day))^2))
  The linear model is of the form hday=(day*max(h))/max(days), Zrday=(day*Zr)/max(days)
    }
 \item{\strong{Point modeling with Single and dual crop coefficients}}{ 
 The function can estimate kc curve for the growing period by using the traditional
 kcini, kcmid and kcend and crop lengths growth stages, in case of single crop coefficients.
 Just like single crop coefficients, the dual crop coefficients kcbini, kcbmid, kcbend are
 embedded in kc parameter. The "kc" parameter must take 3 
 values as c(initial, mid, end). The values of "lengths" parameter must also be specified as
 c(initial, dev, mid, late). See example section for more details.}
  \item{\strong{Spatial modeling With Evaporation Fraction}}{ 
  Instead of supplying kcs values as described above, Evaporation Fraction (EF) from
  can be used. The corresponding EF at initital, mid amd end seasons can be supplied.
  The EF curve will be developed for each location at each time step and ETa derived 
  from (EF*Rn)/2.45.
  }
  \item{\strong{Interpolation of Input data}}{ Once EF is used instead of Kcs, one can
  provide  longitudes and latitudes in xydata parameter data frame. A 2D column data of P,
  I,ETo,RHmin, u2, Rn for each time step for each location can be provided. At each time step 
  an inverse distance interpolation will be perform to get continuous surface. In case of initial
  and static input column data such as CNII, FC, WP, initgw, Zr, h can also be interpolated.
  }
  \item{\strong{Calibration}}{A calibration of unmeasured variables can be done with \code{cal.kc}}
 }
}
\examples{
\dontrun{
#FAO56 Example 32: Calculation of the crop coefficient (Kcb + Ke) under sprinkler irrigation
FAO56Example32=kc(ETo=7,P=0,RHmin=20,soil="sandy loam",crop="Broccoli",I=10,
kctype = "dual",u2=3,kc=c(0.9,0,0),h=1,Zr=0.1)
FAO56Example32$output$fc
FAO56Example32$output$Kc_adj
#FAO56 Example 33: Calculation of the crop coefficient (Kcb + Ke) under furrow irrigation
FAO56Example33=kc(ETo=7,P=0,RHmin=20,soil="sandy loam",crop="Broccoli",I=10,
         kctype = "dual",u2=3,kc=c(0.9,0,0),h=1,Zr=0.1,fw=0.3)
#FAO56 Example 34: Calculation of the crop coefficient (Kcb + Ke) under drip irrigation
FAO56Example34=kc(ETo=7,P=0,RHmin=20,soil="sandy loam",crop="Broccoli",I=10,
           kctype = "dual",u2=3,kc=c(0.9,0,0),h=1,Zr=0.1,fw="drip")
         
file=system.file("extdata","sys","irrigation.txt",package="sebkc")
data=read.table(file,header=TRUE)  
P=data$P
rc=0
I=data$I
ETo=data$ETo
Zr=data$Zr
p=0.6
FC=0.23
WP=0.10
u2=1.6  
RHmin=35
soil="sandy loam"
crop="Broccoli"
############Single crop coefficient model###########
modsingle=kc(ETo,P=P,RHmin,soil,crop,I)
#dataframe of the return values
modsingle$output 
modsingle$output$TAW #Total Avalaible water [mm]
#Single with runoff coefficient
############Dual crop coefficient model###########
moddual=kc(ETo,P=P,RHmin,soil,crop,I,kctype = "dual",FC=FC,p=p,WP=WP,u2=u2,kc=c(0.3,0,0),Zr=Zr)
#Reproducing FAO56 Example 36
fc=(0.00667*1:12)+0.07333
h=1:12/1:12*0.3
 moddual=kc(ETo,P=P,RHmin,soil,crop,I,kctype = "dual",FC=FC,p=p,WP=WP,u2=u2,rc=0,
 kc=data$Kcb,Zr=Zr,lengths = c(4,4,2,2),h=h,hmodel="linear",fc=fc,fw=0.8)
 plot(moddual$output$Day,moddual$output$Kc_adj,type="l",ylim=c(0,1.2))
 lines(moddual$output$Day,moddual$output$Kcb,type="l", lty=2)
###########Spatial model I: Evaporation Fraction (EF) Model###################
#landsat folder with original files but a subset
folder=system.file("extdata","stack",package="sebkc") 

sebiauto=sebi(folder=folder,welev=317,Tmax=31,Tmin=28) #sebi model
kc2=stack(sebiauto$EF/2,sebiauto$EF,sebiauto$EF/3) #assign EF to kc
modEF=kc(ETo,P=P,RHmin,soil,crop,I,kc=kc2,Rn=8,lengths = c(4,4,2,2))
spplot(modEF$EF)

#############Spatial model II: interpolating precipitation ######################
h=1:12
Zr=1:12/12
xydata=data.frame(cbind(longitude=data$longitude,latitude=data$latitude))
print(xydata) #take a look at xydata
P2=data.frame(matrix(rep(P,12),ncol=12))
print(P2) #take a look at the format of precipitation data
modEFh=kc(ETo,P=P2,RHmin,soil,crop,I,kc=kc2,Rn=8,lengths = c(4,4,2,2),h=h,Zr=Zr,xydata=xydata)
#Amount of irrigation water needed 
spplot(modEFh$Drend)
#Amount of water Irrigation that is critically needed to reach maximum allowed depletion level
modEFh$Ineedxy
ETcxy=modEFh$ETcxy ##get ETc for each location at each time step
TAW=modEFh$TAW ## 


############FAO56 Example 41: Estimation of mid-season crop coefficient for Trees or weeds ########
FAOExample41 =kc(fc=0.3,DOY=200,latitude=40,u2=1.5,h=2,density="row",lengths=c(0,0,1,0),
ETo=0,P=0,RHmin=55,soil="sand",crop="Broccoli",kctype="dual")
FAOExample41$output$Kcb

#FAO56 Example 43: Estimation of Kcb mid from ground cover with reduction for stomatal control
FAOExample43 =kc(fc=0.196,DOY=180,latitude=30,u2=2,h=5,density="round",y=0.0676,slope=0.189,
r1=420,RHmin=25,lengths=c(0,0,1,0),ETo=0,P=0,soil="sand",crop="Broccoli",kctype="dual")
FAOExample43$output$Kcb
}
}
\author{
George Owusu
}
\references{
\itemize{
\item{}{Allen, R. G., PEREIRA, L. S., RAES, D., & SMITH, M. (1998). 
Crop Evapotranspiration (guidelines for computing crop water requirements) 
FAO Irrigation and Drainage Paper No. 56: FAO.}
\item{}{Allen, R.G., Wright, J.L., Pruitt, W.O., Pereira, L.S., Jensen, M.E., 2007. Water requirements.
In: Hoffman, G.J., Evans, R.G., Jensen, M.E., Martin, D.L., Elliot, R.L. (Eds.),
Design and Operation of Farm Irrigation Systems. , 2nd edition. ASABE, St. Joseph,
MI, pp. 208-288.}
\item{}{Liu, Y., Pereira, L.S., Fernando, R.M., 2006. Fluxes through the bottom boundary of 
the root zone in silty soils: parametric approaches to estimate groundwater 
contribution and percolation. Agric. Water Manage. 84, 27-40.}
 
}
}
\seealso{
\code{\link{cal.kc}}
}

