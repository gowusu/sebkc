###########################################
#rooting depth
#' @title Estimating FAO56 Single and Dual Crop coefficients with a water balance model in R
#' 
#' @description This function computes FAO56 crop coefficients and water balance. 
#' The output of the model includes crop evapotranspiration, water stress coefficient,
#' irrigation requirements, soil available water, etc. In addition to point estimation,
#' spatial modelling is integrated using Evaporation fractions from 
#' \code{\link{sebal}}, \code{\link{sebi}}, \code{\link{ETo}},  \code{\link{sseb}},  
#' \code{\link{sebs}}. \code{cal.kc} can be used to calibrate  the model.
#'
#' @param ETo Numeric. Reference evapotranspiration [mm]. This can be estimated with \code{\link{ETo}}
#' @param P numeric. Precipitation [mm] 
#' @param RHmin numeric. Minimum Relative Humidity [per cent] 
#' @param soil character. Name of soil texture as written in FAO56 TABLE 19
#' @param crop character. Name of the crop as written on 
#' FAO56 Table 11, TABLE 12 and TABLE 17 
#' @param I numeric. Irrigation [mm]
#' @param CNII numeric. The initial Curve Number of the study site. 
#' This can be calibrated with \code{cal.kc}.
#' @param u2 numeric. Wind speed of the nearby weather station [m/s] measured at 2m
#' @param FC numeric. soil water content at field capacity [-][m3 m-3]. 
#' This can be calibrated with \code{cal.kc}.
#' @param WP numeric. soil water content at wilting point [-][m3 m-3]. 
#' This can be calibrated with \code{cal.kc}.
#' @param h numeric. Maximum crop height [m]. See details below.
#' This can be calibrated with \code{cal.kc}
#' @param Zr numeric. Rooting depth.  See details below.  
#' This can be calibrated with \code{cal.kc}
#' @param Ze numeric. The depth of the surface soil layer that is subject to drying 
#' by way of evaporation [0.10-0.15 m]. This can be calibrated with \code{\link{cal.kc}}
#' @param kc list. The tabulated Kc values for single crop c(kcini,kcmid,kcend)
#' or  dual crop model c(kcbini,kcbmid,kcbend). To use with \code{\link{sebal}} or 
#' any of the Surface Energy Balance models, it should be linked as list(model1,model2,model3) or
#'  c(model1$EF,model2$EF,model3$EF) but not a mixture of both. 
#' This can be calibrated with \code{\link{cal.kc}}
#' @param p The depletion coefficient. This can be calibrated with \code{\link{cal.kc}}
#' @param lengths list. The lengths of the crop growth stages. It can take dates in the form of 
#' c(start of planting ('YYYY-mm-dd), start of rapid growth, 
#' start of mid season, start of maturity, harvesting) 
#' or in the form number of days: c(initial days, dev days, mid days, late days). 
#' See FAO56 TABLE 11. If performing time series simulation, different dates or 
#' one season dates can be set; for instance two seasons simulation can be of the form 
#' c(initial days, dev days, mid days, late days, initial days, dev days, mid days, late days) 
#' or just c(initial days, dev days, mid days, late days)
#' @param fw list or numeric or character. This can take different fw values
#' throughout the growing season. It can also take one numeric value and this will
#' be transform depending on the type of wetting. If the wetting event is only Precipitation,
#' it is set to 1; if it is only Irrigation it is set to the provided fw value. If it is
#' both Precipitation and Irrigation weights are assigned for calculation of fw. 
#' fw can also take a character such as "Drip", "Sprinkler", "Basin", "Border", 
#' "alternated furrows", "Trickle", "narrow bed" and "wide bed". 
#' In this case fw values from FAO56 are assigned. 
#' @param fc numeric. fraction of ground covered by tree canopy
#' @param kctype character. The type of model either "single" or "double". In spatial modelling involving
#' Evaporation Fraction, kctype can be set to "METRIC" or "SSEB" so that 
#' Actual Evapotranspiration will be estimated as ETo*EF else it will be
#' estimated with Rn as ETa=(EF*Rn)/2.45. Note that Rn is in MJm-2 day-1.
#' @param region character. The region of the crop as on FAO56 TABLE 11
#' @param regexpr character. Extra term that can be used to search crop names on FAO56 Tables
#' @param initgw numeric. Initial ground water level [m]
#' @param LAI numeric. Leaf Area Index
#' @param rc numeric Runoff coefficient [0-1]. 
#' This must be used if there is no information on Curve Number.
#' This can be calibrated with \code{\link{cal.kc}}
#' @param Kb ground water (base flow) recession parameter. 
#' The default is zero where no baseflow occurs
#' @param Rn Net daily radiation [w/m2]. Estimated this with  \code{\link{ETo}}
#' @param hmodel character. The type of crop growth model during the simulation. 
#' It takes "gaussian","linear", "exponential"
#' @param Zrmodel character. The type of root growth model during the simulation. 
#' It takes "gaussian","linear", "exponential"
#' @param xydata list c(longitude, latitude). The observation points of the model. 
#' This should be set only if the model is spatial.
#' @param REW numeric Readily Evaporable Water. It can be calibrated with \code{\link{cal.kc}}
#' @param a1 soil water storage to maximum root depth (Zr) at field capacity
#' @param b1 Liu et al. (2006) parameter for computing capillary rise [-0.17]. It can calibrated with
#' \code{\link{cal.kc}}
#' @param a2 Liu et al. (2006) parameter for computing capillary rise [240]. It can be calibrated with
#' \code{\link{cal.kc}}
#' @param b2 Liu et al. (2006) parameter for computing capillary rise [-0.27]. It can be calibrated with
#' \code{\link{cal.kc}}
#' @param a3 Liu et al. (2006) parameter for computing capillary rise [-1.3]. It can be calibrated with
#' \code{\link{cal.kc}}
#' @param b3 Liu et al. (2006) parameter for computing capillary rise [6.6]. It can be calibrated with
#' \code{\link{cal.kc}}
#' @param a4 Liu et al. (2006) parameter for computing capillary rise [4.6]. It can be calibrated with
#' \code{\link{cal.kc}}
#' @param b4 Liu et al. (2006) parameter for computing capillary rise [-0.65]. It can be calibrated with
#' \code{\link{cal.kc}}
#' @param report list. In case of spatial modelling, which days' map should be reported. 
#' The default is "stages" where the periods of stages maps are reported as values. 
#' In an ascending order it can take values such as 1:12 days or c(1,3,6:12) days.
#' @param time character. Time of the occurrence Irrigation or Precipitation. 
#' @param HWR numeric. Default is 1. Height to width ratio of individual plants or groups of 
#' plants when viewed from the east or from the west [-]. 
#' @param slope numeric. slope of saturation vapour pressure curve at air temperature T [kPa oC-1].
#' @param y numeric. psychometric constant  [kPa oC-1]
#' @param nongrowing character. The type of surface during non-growing season. 
#' It can take "weeds", "bare", or "mulch"
#' @param profile numeric or character. The profile depth to perform water balance.
#' It takes "variable" when updating the depth with root growth (default). Or taking "fixed" when
#' the maximum root zone is used. It can also take numeric value and this can be longer than
#' maximum root length. 
#' @param theta numeric. The measured soil water content [m3/m3] in the 
#' profile. Leave those days without measured data blank.
#' @inheritParams ETo
#' @param density character. density factor for the trees, 
#' vines or shrubs on how they are arranged.
#' It can take "full", "row" or "random" or "round". 
#' @param r1 leaf resistance FOR STOMATAL CONTROL [s/m]
#' @param Kcmin Minimum  Kc value 
#' @param Ky list or array. yield response factor [-] in the form of 
#' c(initial Ky,crop dev't Ky, mid season Ky, late season Ky, total season Ky). 
#' For instance Maize Ky is c(0.4,0.4,1.3,0.5,1.25) and it is available at
#'  \url{http://www.fao.org/nr/water/cropinfo_maize.html}.
#' Other crops Ky can be found at \url{http://www.fao.org/nr/water/cropinfo.html}
#' If only one Ky value is supplied it will be repeated for all the seasons.
#' @param x a return object of the function.
#' @param main Title of the plot
#' @param output The kind of plot that should be displayed. It takes 
#' 'AW_measured' for plotting Available and modelled soil moisture;
#' 'TS_measured' for Total Available and modelled soil moisture;
#' 'Kc' for displaying Kcb and Kc
#' 'Evaporation' for displaying of Evaporation and Transpiration
#' 'balance' for displaying Irrigation, precipitation, ETa and critical water deficit
#' @param legend logical. TRUE or FALSE. Whether a legend should be displayed or not.
#' @param pos the position of legend. It takes one of the following 'top', 'topright',
#' 'topleft', 'bottom', 'bottomright' 'bottomleft', 'center'. This must be set so that
#' the legend should not obscure the main plot.
#' @param cex the relative font size of the texts. 
#' Vlaues: "horiz";	logical: if TRUE, set the legend horizontally rather than vertically
#' 
#' @seealso \code{\link{cal.kc}}
#' @details 
#'  \describe{
#'  \item{\strong{Using Existing Tables}}{The function include backend database of FAO Tables.
#'  The user can use the parameter values such as Zr, h, p, FC, WP, etc from FAO56 tables. 
#'  Most of the parameters with NULL values can take default values from FAO56 Tables. When a crop name results in 3 crops,
#'   the user can differentiate  that with  regexpr parameter. For example, 
#'   a crop parameter with a name "maize" can match
#'   two crops on FAO Table 12. Including  regexpr parameter with "sweet" with match only one. 
#'   }
#'  \item{\strong{Modelling crop height and Rooting Depth}}{ 
#'  At each time step a crop height or root depth is estimated. 
#'  The crop  or rooting depth growth model can take  exponential model of the form
#'  hday=0+0.75*max(h)*(1-exp(-day/0.75*max(days))) or    Zrday=0+0.75*Zr*(1-exp(-day/0.75*max(days)))
#'   The gaussian model is of the form hday=0+0.75*max(h)*(1-exp(-day^2/(0.75*max(days)^2))) or Zrday=0+0.75*(1-exp(-day^2/(0.75*max(day))^2))
#'   The linear model is of the form hday=(day*max(h))/max(days), Zrday=(day*Zr)/max(days)
#'     }
#'  \item{\strong{Point modeling with Single and dual crop coefficients}}{ 
#'  The function can estimate kc curve for the growing period by using the traditional
#'  kcini, kcmid and kcend and crop lengths growth stages, in case of single crop coefficients.
#'  Just like single crop coefficients, the dual crop coefficients kcbini, kcbmid, kcbend are
#'  embedded in kc parameter. The "kc" parameter must take 3 
#'  values as c(initial, mid, end). The values of "lengths" parameter must also be specified as
#'  c(initial, dev, mid, late). See example section for more details.}
#'   \item{\strong{Spatial modeling With Evaporation Fraction}}{ 
#'   Instead of supplying kcs values as described above, Evaporation Fraction (EF) from
#'   can be used. The corresponding EF at initital, mid amd end seasons can be supplied.
#'   The EF curve will be developed for each location at each time step and ETa derived 
#'   from (EF*Rn)/2.45.
#'   }
#'   \item{\strong{Interpolation of Input data}}{ Once EF is used instead of Kcs, one can
#'   provide  longitudes and latitudes in xydata parameter data frame. A 2D column data of P,
#'   I,ETo,RHmin, u2, Rn for each time step for each location can be provided. At each time step 
#'   an inverse distance interpolation will be perform to get continuous surface. In case of initial
#'   and static input column data such as CNII, FC, WP, initgw, Zr, h can also be interpolated.
#'   }
#'   \item{\strong{Calibration}}{A calibration of unmeasured variables can be done with \code{cal.kc}}
#'  }
#' @references 
#' \itemize{
#' \item{}{Allen, R. G., PEREIRA, L. S., RAES, D., & SMITH, M. (1998). 
#' Crop Evapotranspiration (guidelines for computing crop water requirements) 
#' FAO Irrigation and Drainage Paper No. 56: FAO.}
#' \item{}{Allen, R.G., Wright, J.L., Pruitt, W.O., Pereira, L.S., Jensen, M.E., 2007. Water requirements.
#' In: Hoffman, G.J., Evans, R.G., Jensen, M.E., Martin, D.L., Elliot, R.L. (Eds.),
#' Design and Operation of Farm Irrigation Systems. , 2nd edition. ASABE, St. Joseph,
#' MI, pp. 208-288.}
#' \item{}{Liu, Y., Pereira, L.S., Fernando, R.M., 2006. Fluxes through the bottom boundary of 
#' the root zone in silty soils: parametric approaches to estimate groundwater 
#' contribution and percolation. Agric. Water Manage. 84, 27-40.}
#'  
#' }
#'  
#'
#' @return The output of the function depends of the kctype and kc data type.
#' For single and dual crop coefficient, an "output" dataframe contains variables such as
#' TAW,RAW,Runoff,Ks,Kcb,h,Kc,ETc,DP,CR,Drend,Kr,De,AW,E,Ke,CN. In case kc is set to EF the function returns
#' Kc,Kcxy,ETc,ETcxy,Drend,Drendxy,TAW,TAWxy,RAW,RAWxy,AW,Awcxy,Ineed,Ineedxy
#'  \itemize{
#' \item{output:} {  dataframe of the return values}
#' \item{ETc_adj:} { The soil water stress adjusted ETc  [mm]}
#' \item{Ya_by_Ym:} {actual yield by potential yield}
#' \item{yield_decrease:} {The relative yield decrease: 1-Ya/Ym}
#' \item{ETdecrease:} {The relative evapotranspiration deficit :1- ETc_adj/ETc}
#' \item{yield_by_ET_decrease:} {he relative yield decrease: 1-Ya/Ym and 
#' the relative evapotranspiration deficit in one dataframe}
#' \item{Awcxy:} { Available Soil water for each location at each time step [mm]}
#' \item{CR:} { capillary rise [mm/day]}
#' \item{GW:} { Updated ground water depth}
#' \item{CN:} { Curve Number}
#' \item{De:} { cumulative depth of evaporation (depletion) 
#' from the soil surface layer [mm] at the end of the day [mm/day]}
#' \item{DP:} { Deep Percolation [mm/day]}
#' \item{Drend:} { cumulative depth of evapotranspiration (depletion) 
#' from the root zone at the end of the day [mm]. 
#' This is equal to the amount of water needed to bring soil water to field capacity}
#' \item{Drendxy:} { cumulative depth of evapotranspiration (depletion) 
#' from the root zone at the end of the day [mm] for each xydata location at each time step. 
#' This is equal to the amount of water needed to bring soil water to field capacity}
#' \item{E:} { evaporation [mm/day]}
#' \item{ETc:} { crop evapotranspiration [mm/day]. 
#' If water stress coefficient (Ks) is less than 1, 
#' ETc is crop evapotranspiration under non-standard conditions else ETc is
#' crop evapotranspiration under standard conditions}
#' \item{ETcxy:} { crop evapotranspiration [mm/day] for each location at each time step.
#' crop evapotranspiration under standard conditions}
#' \item{h:} { crop height [m]}
#' \item{Ineed:} { Amount of soil water that is below maximum allowed depletion [mm]}
#' \item{Ineedxy:} { Amount of soil water that is below maximum allowed depletion 
#' for each xydata location at each time step [mm] }
#' \item{Kc:} { standard crop coefficient [-]}
#' \item{Kc_adj:} { adjusted crop coefficient [-]}
#' \item{EF:} { Evaporation Fraction [-]}
#' \item{EFxy:} { Evaporation Fraction for each monitored location [-]}
#' \item{kcb:} { basal crop coefficient for the growing season [-]}
#' \item{Kcxy:} { crop coefficient for each location at each time step [-]}
#' \item{Ke:} { soil evaporation coefficient}
#' \item{Kr: } { soil evaporation reduction coefficient [-]}
#' \item{Ks:} { Water stress coefficient [-]}
#' \item{RAW: } { Readily Available Water [mm]}
#' \item{RAwcxy:} { Readily Available Soil water for each xydata location at each time step [mm]}
#' \item{Runoff:} { Runoff [mm/day]}
#' \item{TAW:} { Total Available Water [mm] }
#' \item{TAwcxy:} { Total Available Soil water for each location at each time step [mm]}
#' }
#' 
#' 
#' @author George Owusu
#' @export
#'
#' @examples
#' \dontrun{
#' #FAO56 Example 32: Calculation of the crop coefficient (Kcb + Ke) under sprinkler irrigation
#' FAO56Example32=kc(ETo=7,P=0,RHmin=20,soil="sandy loam",crop="Broccoli",I=10,
#' kctype = "dual",u2=3,kc=c(0.9,0,0),h=1,Zr=0.1)
#' FAO56Example32$output$fc
#' FAO56Example32$output$Kc_adj
#' #FAO56 Example 33: Calculation of the crop coefficient (Kcb + Ke) under furrow irrigation
#' FAO56Example33=kc(ETo=7,P=0,RHmin=20,soil="sandy loam",crop="Broccoli",I=10,
#'          kctype = "dual",u2=3,kc=c(0.9,0,0),h=1,Zr=0.1,fw=0.3)
#' #FAO56 Example 34: Calculation of the crop coefficient (Kcb + Ke) under drip irrigation
#' FAO56Example34=kc(ETo=7,P=0,RHmin=20,soil="sandy loam",crop="Broccoli",I=10,
#'            kctype = "dual",u2=3,kc=c(0.9,0,0),h=1,Zr=0.1,fw="drip")
#'          
#' file=system.file("extdata","sys","irrigation.txt",package="sebkc")
#' data=read.table(file,header=TRUE)  
#' P=data$P
#' rc=0
#' I=data$I
#' ETo=data$ETo
#' Zr=data$Zr
#' p=0.6
#' FC=0.23
#' WP=0.10
#' u2=1.6  
#' RHmin=35
#' soil="sandy loam"
#' crop="Broccoli"
#' ############Single crop coefficient model###########
#' modsingle=kc(ETo,P=P,RHmin,soil,crop,I)
#' #dataframe of the return values
#' modsingle$output 
#' modsingle$output$TAW #Total Avalaible water [mm]
#' #Single with runoff coefficient
#  modsinglerc=kc(ETo,P=P,RHmin,soil,crop,I,rc=0.2)

#' ############Dual crop coefficient model###########
#' moddual=kc(ETo,P=P,RHmin,soil,crop,I,kctype = "dual",FC=FC,p=p,WP=WP,u2=u2,kc=c(0.3,0,0),Zr=Zr)
#' #Reproducing FAO56 Example 36
#' fc=(0.00667*1:12)+0.07333
#' h=1:12/1:12*0.3
#'  moddual=kc(ETo,P=P,RHmin,soil,crop,I,kctype = "dual",FC=FC,p=p,WP=WP,u2=u2,rc=0,
#'  kc=data$Kcb,Zr=Zr,lengths = c(4,4,2,2),h=h,hmodel="linear",fc=fc,fw=0.8)
#'  plot(moddual$output$Day,moddual$output$Kc_adj,type="l",ylim=c(0,1.2))
#'  lines(moddual$output$Day,moddual$output$Kcb,type="l", lty=2)
#' ###########Spatial model I: Evaporation Fraction (EF) Model###################
#' #landsat folder with original files but a subset
#' folder=system.file("extdata","stack",package="sebkc") 
#' 
#' sebiauto=sebi(folder=folder,welev=317,Tmax=31,Tmin=28) #sebi model
#' kc2=stack(sebiauto$EF/2,sebiauto$EF,sebiauto$EF/3) #assign EF to kc
#' modEF=kc(ETo,P=P,RHmin,soil,crop,I,kc=kc2,Rn=8,lengths = c(4,4,2,2))
#' spplot(modEF$EF)
#'
#'#############Spatial model II: interpolating precipitation ######################
#' h=1:12
#' Zr=1:12/12
#' xydata=data.frame(cbind(longitude=data$longitude,latitude=data$latitude))
#' print(xydata) #take a look at xydata
#' P2=data.frame(matrix(rep(P,12),ncol=12))
#' print(P2) #take a look at the format of precipitation data
#' modEFh=kc(ETo,P=P2,RHmin,soil,crop,I,kc=kc2,Rn=8,lengths = c(4,4,2,2),h=h,Zr=Zr,xydata=xydata)
#' #Amount of irrigation water needed 
#' spplot(modEFh$Drend)
#' #Amount of water Irrigation that is critically needed to reach maximum allowed depletion level
#' modEFh$Ineedxy
#' ETcxy=modEFh$ETcxy ##get ETc for each location at each time step
#' TAW=modEFh$TAW ## 
#' 
#' 
#' ############FAO56 Example 41: Estimation of mid-season crop coefficient for Trees or weeds ########
#' FAOExample41 =kc(fc=0.3,DOY=200,latitude=40,u2=1.5,h=2,density="row",lengths=c(0,0,1,0),
#' ETo=0,P=0,RHmin=55,soil="sand",crop="Broccoli",kctype="dual")
#' FAOExample41$output$Kcb
#' 
#' #FAO56 Example 43: Estimation of Kcb mid from ground cover with reduction for stomatal control
#' FAOExample43 =kc(fc=0.196,DOY=180,latitude=30,u2=2,h=5,density="round",y=0.0676,slope=0.189,
#' r1=420,RHmin=25,lengths=c(0,0,1,0),ETo=0,P=0,soil="sand",crop="Broccoli",kctype="dual")
#' FAOExample43$output$Kcb
#' }
# generic
kc<-function(ETo,P,RHmin,soil="Sandy Loam",crop,I=0,CNII=67,u2=2,FC=NULL,WP=NULL,h=NULL,
             Zr=NULL,Ze=0.1,kc=NULL,p=NULL,lengths=NULL,fw=0.8,fc=NULL,
             kctype="dual",region=NULL,regexpr="sweet",initgw=50,
             LAI=3,rc=NULL,Kb=0,Rn=NULL,hmodel="gaussian",Zrmodel="linear",
             xydata=NULL,REW=NULL,a1=360,b1=-0.17,a2=240,b2=-0.27,
             a3=-1.3,b3=6.6,a4=4.6,b4=-0.65,report="stages",time="06:00",
             HWR=1,latitude=NULL,density="full",DOY=NULL,slope=NULL,
             y=NULL,r1=100,nongrowing="weeds",theta=NULL,profile="variable",
        Kcmin=0.15,Ky=NULL
             ) UseMethod ("kc")

#' @rdname kc
#' @export
#default function
kc.default<-function(ETo,P,RHmin,soil="Sandy Loam",crop,I=0,CNII=67,u2=2,FC=NULL,WP=NULL,h=NULL,
                     Zr=NULL,Ze=0.1,kc=NULL,p=NULL,lengths=NULL,fw=1,fc=NULL,
                     kctype="dual",region=NULL,regexpr="sweet",initgw=50,
                     LAI=3,rc=NULL,Kb=0,Rn=NULL,hmodel="gaussian",Zrmodel="linear",
                     xydata=NULL,REW=NULL,a1=360,b1=-0.17,a2=240,b2=-0.27,
                     a3=-1.3,b3=6.6,a4=4.6,b4=-0.65,report="stages",time="06:00",
                     HWR=1,latitude=NULL,density="full",DOY=NULL,slope=NULL,
                     y=NULL,r1=100,nongrowing="weeds",theta=NULL,profile="variable",
        Kcmin=0.15,Ky=NULL)
{
  options(warn=-1)
  #get weather data
  
  if(class(kc[1])!="numeric"){
  if(length(kc)==3){
  if(!is.null(kc[[1]]$EF)){
    kctype=class(kc[[1]])
    kc[[1]]=kc[[1]]$EF
  }
  
  if(!is.null(kc[[2]]$EF)){
    kctype=class(kc[[2]])
    kc[[2]]=kc[[2]]$EF
  }
  
  if(!is.null(kc[[3]]$EF)){
    kctype=class(kc[[3]])
    if(!is.null(kc[[3]]$ETr)){
      kctype="METRIC" 
    }
    kc[[3]]=kc[[3]]$EF
  }
    kc=stack(kc[[1]],kc[[2]],kc[[3]])
  }
  }
  
    #print(kctype)
  lengths2=lengths
  #lengths5=lengths
  #lengths=lengths2
  if(!is.null(lengths)){
    if(!is.numeric(lengths)){
      init=as.numeric(as.Date(lengths[2], "%Y%m/%d")-as.Date(lengths[1], "%Y%m/%d"))
      if(is.na(init)){
        init=as.numeric(as.Date(lengths[2], "%Y-%m-%d")-as.Date(lengths[1], "%Y-%m-%d"))
      }
      
      dev=as.numeric(as.Date(lengths[3], "%Y%m/%d")-as.Date(lengths[2], "%Y%m/%d"))
      if(is.na(dev)){
        dev=as.numeric(as.Date(lengths[3], "%Y-%m-%d")-as.Date(lengths[2], "%Y-%m-%d"))
      }
      
      mid=as.numeric(as.Date(lengths[4], "%Y%m/%d")-as.Date(lengths[3], "%Y%m/%d"))
      if(is.na(mid)){
        mid=as.numeric(as.Date(lengths[4], "%Y-%m-%d")-as.Date(lengths[3], "%Y-%m-%d"))
      }
      
      late=as.numeric(as.Date(lengths[5], "%Y%m/%d")-as.Date(lengths[5], "%Y%m/%d"))
      if(is.na(late)){
        late=as.numeric(as.Date(lengths[5], "%Y-%m-%d")-as.Date(lengths[4], "%Y-%m-%d"))
      }
      
      lengths=c(init,dev,mid,late)
      lengthscum=cumsum(lengths)
    }
  
  if(!is.null(theta)&&length(theta)==2&&is.character(lengths2)&&!is.numeric(theta[1])){
    thisinit=(as.Date(lengths2[1], "%Y-%m-%d"))
   # print("yes")
    
    if(is.na(thisinit)){
      thisinit=(as.Date(lengths2[1], "%Y/%m/%d"))
    }
    thisend=(as.Date(lengths2[5], "%Y-%m-%d"))
    
    if(is.na(thisend)){
      thisend=(as.Date(lengths2[5], "%Y-%m-%d"))
    }
    
    #print("yes")
        
    modlen=as.data.frame(seq(thisinit,length=as.numeric(thisend-thisinit)+1,by=1))
    names(modlen)="date"
    
    #theta=data.frame(date=c("1960-04-01", "1960-04-06","1960-4-10","1960-5-4"),theta=c(0.2,0.3,0.4,0.5))
    names(theta)=c("date","theta")
    thisdata2 <- merge(modlen, theta,by=c("date"),all=TRUE)
    theta=NULL
   # date=NULL
    for(i in 1:nrow(thisdata2)){
      if(thisdata2$date[i]!=thisdata2$date[i+1]&&i<nrow(thisdata2)){
        theta=rbind(theta,thisdata2$theta[i])
       # date=rbind(date,thisdata2$date[i])
        
        #print(c(thisdata2$date[i]))
      }
      
    }
    theta=as.data.frame(rbind(theta,thisdata2$theta[i]))
    theta=theta[[1]]
    #print (thisdata2)
    #date=as.data.frame(rbind(date,thisdata2$date[i]))
    #print(theta)
    #names(theta)="date"
  }
  }
  
  if(!is.null(theta)&&length(theta)==2&&is.numeric(lengths2)){
    theta2=NULL
  for (i in 1:max(cumsum(lengths2))) {
    testtheta=theta[which(theta[1]==i),]
    if(nrow(testtheta)>0){
      theta2=rbind(theta2,testtheta[[2]])
      
    }else{
      theta2=rbind(theta2,NA)
      
    }
  } 
    theta=theta2
    
  }
  #print(lengths)
  if(class(ETo)=="ETo"||class(ETo)=="ETo2"){
    if(!is.null(ETo$ETo)){
      latitude=ETo$latitude
      u2=ETo$u2
      if(!is.null(y)){
        y=ETo$y
      }
      #y=ETo$y
      DOY=ETo$DOY
      if(!is.null(ETo$data)){
        thisETo=ETo$data
      }else{
        thisETo=ETo
      }
      
      if(!is.null(thisETo$u2)){
        u2=thisETo$u2
      }
      if(!is.null(thisETo$theta)){
        theta=thisETo$theta
      }
      RHmin=thisETo$RHmin
      
      if(is.null(slope)){
        thisETo$slope=NULL
      }
      Rn=thisETo$Rn
      if(!is.null(thisETo$DOY)){
        DOY=thisETo$DOY
      }
      ETo=thisETo$ETo
      
    }
    
  }
  #print(thisETo)
  daylenth=nrow(data.frame(ETo))-1
  day=1:daylenth
  
  #day=1:1522
  if(!is.null(lengths)&&!is.numeric(lengths2)&&day!=lengthscum[4]){
    #DOY=as.Date(DOY)
    #loopdates=as.Date(DOY)-100
    #thisdata=datasub(startdate=DOY[1],thisETo,substart=lengths2[1],subend=lengths2[5])$moddata
    #class(thisdata)="ETo2"
    years=unique( format(as.Date(DOY),'%Y'))
    lengths3=as.Date(lengths2)
    thisYear=format(lengths3,'%Y')
    thisMonth=( format(lengths3,'%m'))
    thisday=( format(lengths3,'%d'))
    startdate=DOY[1]
    enddate=DOY[length(DOY)]

    output=data.frame()
    yield_decrease=data.frame()
    water.balance=data.frame()
    DL=data.frame(Tday=numeric(),T24=numeric(),Tn=numeric(),Rs=numeric(),DL=numeric(),
                  plantdate=character(),harvestdate=character())
    lenmin=0
    lenmax=0
    if((length(lengths2)/5)==length(years)){
      lenmin=1
      lenmax=5
      #print("lenmin")
    }
    for(i in 1:length(years)){
      thislengths1=paste(years[i],thisMonth[1],thisday[1],sep="-")
      thislengths2=paste(years[i],thisMonth[2],thisday[2],sep="-")
      thislengths3=paste(years[i],thisMonth[3],thisday[3],sep="-")
      thislengths4=paste(years[i],thisMonth[4],thisday[4],sep="-")
      thislengths5=paste(years[i],thisMonth[5],thisday[5],sep="-")
      thislengths=c(thislengths1,thislengths2,thislengths3,
                        thislengths4,thislengths5)
      
      substart=thislengths[1]
      subend=thislengths[5]
      substart=as.Date(substart)
      subend=as.Date(subend)
      
      startdate=as.Date(startdate)
      enddate=as.Date(enddate)
      
      modlen=seq(startdate,length=as.numeric(enddate-startdate)+1,by=1)
      #as.numeric(format(modlen,"%m"))
      #print(thislengths)
      gcmdata=data.frame(cbind(modlen,thisETo))
      thisdata <- gcmdata[ which(gcmdata$modlen>=substart&gcmdata$modlen<= subend), ]
      #thisdata=datasub(s,thisETo,substart=lengths2[1],subend=lengths2[5])$moddata
      #class(thisdata)="ETo2"
      #print(lengths)
      #thislengthsminmax=lengths
      if(lenmin>0){
      thislengthsminmax=lengths2[lenmin:lenmax]

      init=as.numeric(as.Date(thislengthsminmax[2], "%Y%m/%d")-as.Date(thislengthsminmax[1], "%Y%m/%d"))
      if(is.na(init)){
        init=as.numeric(as.Date(thislengthsminmax[2], "%Y-%m-%d")-as.Date(thislengthsminmax[1], "%Y-%m-%d"))
      }
      
      dev=as.numeric(as.Date(thislengthsminmax[3], "%Y%m/%d")-as.Date(thislengthsminmax[2], "%Y%m/%d"))
      if(is.na(dev)){
        dev=as.numeric(as.Date(thislengthsminmax[3], "%Y-%m-%d")-as.Date(thislengthsminmax[2], "%Y-%m-%d"))
      }
      
      mid=as.numeric(as.Date(thislengthsminmax[4], "%Y%m/%d")-as.Date(thislengthsminmax[3], "%Y%m/%d"))
      if(is.na(mid)){
        mid=as.numeric(as.Date(thislengthsminmax[4], "%Y-%m-%d")-as.Date(thislengthsminmax[3], "%Y-%m-%d"))
      }
      
      late=as.numeric(as.Date(thislengthsminmax[5], "%Y%m/%d")-as.Date(thislengthsminmax[5], "%Y%m/%d"))
      if(is.na(late)){
        late=as.numeric(as.Date(thislengthsminmax[5], "%Y-%m-%d")-as.Date(thislengthsminmax[4], "%Y-%m-%d"))
      }
        
      lengths=c(init,dev,mid,late)

      
      lenmin=lenmin+5
      lenmax=lenmax+5
      }
      #print(theta)
      model=kc(ETo=thisdata$ETo,P=thisdata$P,RHmin=thisdata$RHmin,soil=soil,crop=crop,
               I=I,CNII=CNII,u2=thisdata$u2,
             FC=FC,WP=WP,h=h,Zr=Zr,Ze=Ze,kc=kc,p=p,fw=fw,fc=fc,lengths=lengths,
             kctype=kctype,region=region,regexpr=regexpr,initgw=initgw,LAI=LAI,rc=rc,
             Rn=thisdata$Rn,hmodel=hmodel,Zrmodel=Zrmodel,a1=a1,b1=b1,a2=a2,b2=b2,
             a3=a3,b3=b3,a4=a4,b4=b4,xydata=xydata,REW=REW
             ,report=report,time=time,HWR=HWR,latitude=latitude,
             density=density,DOY=thisdata$DOY,slope=thisdata$slope,y=y,r1=r1,
             nongrowing=nongrowing,theta=theta,profile=profile,Kcmin=Kcmin,Ky=Ky)
      #output2=rbind(output2,model$output)
      #print(length(model$output$P))
      
      thisyear=as.numeric(years[i])
      outputs=cbind(Year=years[i],model$output)
      output=rbind(output,outputs)
      Ya_by_Yms=cbind(Year=thisyear,model$yield_by_ET_decrease)
      yield_decrease=rbind(yield_decrease,Ya_by_Yms)
      thisDL=dl(latitude,DOY=thisdata$DOY,model= "CBM",Tmax=thisdata$Tmax,Tmin=thisdata$Tmin)
      DL=rbind(DL,data.frame(Tday=mean(thisDL$Tday),
      T24=mean(thisDL$T24, na.rm = TRUE),Tn=mean(thisDL$Tn, na.rm = TRUE),Rs=mean(thisDL$Rs, na.rm = TRUE),latitude=latitude,
      DL=mean(thisDL$DL, na.rm = TRUE),
      plantdate=substart,harvestdate=subend))
      water.balances=cbind(Year=years[i],model$water.balance.stages)
      water.balance=rbind(water.balance,water.balances)
     # print(substart)
      #print(average(thisDL))
      #print(years[i])
      #print(cbind(ETo=thisdata$ETo,P=thisdata$P,RHmin=thisdata$RHmin,u2=thisdata$u2,
           #  Rn=thisdata$Rn,DOY=thisdata$DOY,slope=thisdata$slope))
      #print(model$output)
      
#print (years[i])    
    }
    if(class(kc[[1]])=="RasterLayer"||class(kc[[2]])=="RasterLayer"
       ||class(kc[[3]])=="RasterLayer"){
      factor=model
    }else{
      factor=list(output=output,yield_by_ET_decrease=yield_decrease,DL=DL,water.balance.stages=water.balance,kctype=kctype)
    }#print(thisdata$Tmax)
    
  }else{
  ##The main model starts here
  
  parameters=c(ETo=ETo,P=P,RHmin=RHmin,soil=soil,crop=crop,I=I,CNII=CNII,u2=u2,
               FC=FC,WP=WP,h=h,Zr=Zr,Ze=Ze,kc=kc,p=p,lengths=lengths,fc=fc,
               kctype=kctype,region=region,regexpr=regexpr,initgw=initgw,LAI=LAI,rc=rc,
               Rn=Rn,hmodel=hmodel,Zrmodel=Zrmodel,a1=a1,b1=b1,a2=a2,b2=b2,
               a3=a3,b3=b3,a4=a4,b4=b4,xydata=xydata,REW=REW
               ,report=report,time=time,HWR=HWR,latitude=latitude,
               density=density,DOY=DOY,slope=slope,y=y,r1=r1,
               nongrowing=nongrowing,theta=theta,profile=profile,Ky=Ky)
 
   
  
  latitude2=latitude
  density2=density
  DOY2=DOY
  slope2=slope
  y2=y
  r12=r1
  
  
  if(!is.null(DOY)){
  if(!is.numeric(DOY)){
    DOY2=DOY
    DOY=strptime(DOY,"%Y-%m-%d")$yday+1
    if(is.na(DOY[1])){
      DOY=strptime(DOY2,"%Y%m/%d")$yday+1
    }
  }
  }
  Zr3=Zr
  h3=h
  kc3=kc
  
  EF=NULL
  allEF=NULL
  if(class(kc[[1]])=="RasterLayer"||class(kc[[2]])=="RasterLayer"
     ||class(kc[[3]])=="RasterLayer"){
    if(class(kc[[1]])=="RasterLayer"){
      map=kc[[1]]
    }else{
      if(class(kc[[2]])=="RasterLayer"){
        map=kc[[2]]
      }else{
        map=kc[[3]]
      }
    }
    EF=TRUE 
  }
  
  
  if(!is.null(EF)&&is.null(Rn)){
    return(print("Daily net radiation (Rn) is needed"))
  }
  ##preparing interpolation files
  Pinter=NULL
  if(!is.null(EF)&&!is.null(xydata)){
  if(!is.null(nrow(P))){
    Pinter=TRUE
  }
  }
  Rninter=NULL
  if(!is.null(EF)&&!is.null(xydata)){
    if(!is.null(nrow(Rn))){
      Rninter=TRUE
    }
  }
  RHinter=NULL
  if(!is.null(EF)&&!is.null(xydata)){
    if(!is.null(nrow(RHmin))){
      RHinter=TRUE
    }
  }
  
  u2inter=NULL
  if(!is.null(EF)&&!is.null(xydata)){
    if(!is.null(nrow(u2))){
      u2inter=TRUE
    }
  }
  ETointer=NULL
  if(!is.null(EF)&&!is.null(xydata)){
    if(!is.null(nrow(ETo))){
      ETointer=TRUE
    }
  }
  Iinter=NULL
  if(!is.null(EF)&&!is.null(xydata)){
    if(!is.null(nrow(I))){
      Iinter=TRUE
    }
  }
  thetainter=NULL
  if(!is.null(EF)&&!is.null(xydata)){
    if(!is.null(nrow(theta))){
      thetainter=TRUE
    }
  }
  folder=system.file("extdata","sys",package="sebkc")
  
 
  
  soilwater=read.csv(system.file("extdata","sys","soilwater2.csv",package="sebkc"),header=TRUE,stringsAsFactors=FALSE)
  stages=(read.csv(system.file("extdata","sys","stages.csv",package="sebkc"),stringsAsFactors=FALSE,header=TRUE))
  crop_H_Zr=na.omit(read.csv(system.file("extdata","sys","crop_H_Zr.csv",package="sebkc"),stringsAsFactors=FALSE,header=TRUE))
  singlekc=na.omit(read.csv(system.file("extdata","sys","singlekc.csv",package="sebkc"),stringsAsFactors=FALSE,header=TRUE))
  dualkc=na.omit(read.csv(system.file("extdata","sys","dualkc.csv",package="sebkc"),stringsAsFactors=FALSE,header=TRUE))
  ky=na.omit(read.csv(system.file("extdata","sys","ky.csv",package="sebkc"),stringsAsFactors=FALSE,header=TRUE))
  wetting=na.omit(read.csv(system.file("extdata","sys","fw.csv",package="sebkc"),stringsAsFactors=FALSE,header=TRUE))
   #for crop heigh and roots
  soilwater$ID=1:9
  range=0.75
  soil2=soil
  crop2=crop
  CN=CNII
  report2=NULL
  daylenths=nrow(data.frame(P))
  day=1:daylenths
  #print(day)
  ##all kc is EF
  if(!is.null(EF)){
    if(nlayers(kc)==max(day)){
    allEF=TRUE
    }
  }
  ETcxy=NULL
  ones=day/day
  if(I==0&&length(I)==1){
    I=day*0
  }else{
    if(length(I)==1){
      I=ones*I
    }
  }
  
  if(length(Rn)==1){
    Rn=ones*Rn
  }
  
  if(length(u2)==1){
    u2=u2*ones
  }
  if(length(RHmin)==1){
    RHmin=RHmin*ones
  }
  if(length(DOY)==1){
    DOY=ones*DOY
  }
  if(length(HWR)==1){
    HWR=ones*HWR
  }
  if(length(latitude)==1){
    latitude=ones*latitude
  }
  
  if(length(slope)==1){
    slope=ones*slope
  }
  
  if(length(y)==1){
    y=ones*y
  }
  
  if(length(r1)==1){
    r1=ones*r1
  }
  
  fc2=fc
  if(length(fc2)==1){
    fc2=ones*fc2
  }
  data=data.frame(cbind(ETo,P,I))
  
  #WP and FC interpolation
  if(!is.null(xydata)&&length(WP)>1&&!is.null(EF)){
    wpdata=data.frame(cbind(xydata,WP))
    longitude=wpdata[1]
    latitude=wpdata[2]
    var=wpdata[3]
    #print(wpdata[3])
    WP= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
    
  }
  
  if(!is.null(xydata)&&length(FC)>1&&!is.null(EF)){
    FCdata=data.frame(cbind(xydata,FC))
    longitude=FCdata[1]
    latitude=FCdata[2]
    var=FCdata[3]
    FC= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
  }
  if(!is.null(xydata)&&length(initgw)>1&&!is.null(EF)){
    wpdata=data.frame(cbind(xydata,initgw))
    longitude=wpdata[1]
    latitude=wpdata[2]
    var=wpdata[3]
    initgw= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
  }
  
  if(!is.null(xydata)&&length(CNII)>1&&!is.null(EF)){
    wpdata=data.frame(cbind(xydata,CNII))
    longitude=wpdata[1]
    latitude=wpdata[2]
    var=wpdata[3]
    CNII= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
    print(CNII)
  }
  if(is.null(Zr)||is.null(p)){
    crop= crop_H_Zr[(grep(crop2,crop_H_Zr$crop,ignore.case=TRUE)),] 
    
    if(nrow(crop)>1){
      crop= crop[(grep((regexpr),(crop),ignore.case=TRUE)),] 
    }
    if(nrow(crop)==0){
      return (print(paste(crop2," not found [for crop type name]"))) 
    }
    if(is.null(Zr)){
      Zr=mean(c(as.numeric(strsplit(crop$Zr, "-")[[1]][1]),as.numeric(strsplit(crop$Zr, "-")[[1]][2])))
    }
    
    if(is.null(p)){
      p=crop$p[1]
    }
  }
  
  Zrinter=NULL
  co=0
  c1=range*Zr
  a=range*max(day)
  if(length(Zr)==1){
    if(Zrmodel=="exponential"||Zrmodel=="exponent"||Zrmodel=="exp"){
      Zr=co+c1*(1-exp(-day/a))
      # Zr=(co+c1)%o%(1-exp(-day/a))
    }else{
      
      if(Zrmodel=="linear"||Zrmodel=="lin"||Zrmodel=="line"){
        
        Zr=(day*Zr)/max(day)
        
      }else{
        Zr=co+c1*(1-exp(-day^2/a^2))
      }
    }
  }else{
    Zrinter=TRUE
    if(Zrmodel=="exponential"||Zrmodel=="exponent"||Zrmodel=="exp"){
      
      Zr=(co+c1)%o%(1-exp(-day/a))
      Zr=t(Zr)
    }else{
      
      if(Zrmodel=="linear"||Zrmodel=="lin"||Zrmodel=="line"){
        
        Zr=(day%o%Zr)/max(day)
        Zr=t(Zr)
      }else{
        Zr=co+c1%o%(1-exp(-day^2/a^2))
        Zr=t(Zr)
      }
    } 
  }
  
  #soil characteristics
  soil=toupper(soil)
  soil=soilwater[(grep(soil2,soilwater$type,ignore.case=TRUE)),]
  if(nrow(soil)==0){
    return (print(paste(soil2," not found: check soil"))) 
  }
  
  #lengths of crop growth stages
  if(is.null(lengths)){
    #stages2= stages[(grep(toupper(crop2),toupper(stages$crop))),] 
    stages2= stages[(grep((crop2),(stages$crop),ignore.case=TRUE)),] 
    if(nrow(stages2)>1){
      stages2= stages2[(grep((region),(stages2$region),ignore.case=TRUE)),] 
    }
    
    if(nrow(stages2)==0){
      return (print(paste(crop2," not found crop growth stages"))) 
    }
    init=mean(as.numeric(stages2$init))
    dev=mean(as.numeric(stages2$dev))
    mid=mean(as.numeric(stages2$mid))
    late=mean(as.numeric(stages2$late))
    
  }else{
    init=lengths[1]
    dev=lengths[2]
    mid=lengths[3]
    late=lengths[4]
  }
  stages2=c(init,dev,mid,late)
  stages2cum=cumsum(stages2)
  ##dual crop coefficient
  if(is.null(kc)){
    kcs= dualkc[(grep((crop2),(dualkc$crop),ignore.case=TRUE)),] 
    if(nrow(kcs)>1){
      kcs= kcs[(grep((regexpr),(kcs$crop),ignore.case=TRUE)),] 
    }
    
    if(nrow(kcs)==0){
      return (print(paste(" Kc not found for", crop2))) 
    }
    kcbini=mean(as.numeric(kcs$Kcbini))
    kcbmid=mean(as.numeric(kcs$Kcbmid))
    kcbend=mean(as.numeric(kcs$Kcbend))
    if(is.na(kcbend)){
      kcbend=mean(as.numeric(c(strsplit(kcs$Kcbend, "-")[[1]][1],strsplit(kcs$Kcbend, "-")[[1]][2])))
      if(is.na(kcbend)){
        kcbend=mean(as.numeric(c(strsplit(kcs$Kcbend, ",")[[1]][1],strsplit(kcs$Kcbend, ",")[[1]][2])))
        
      }
    }
    
    if(is.na(kcbmid)){
      kcbmid=mean(as.numeric(c(strsplit(kcs$Kcmid, "-")[[1]][1],strsplit(kcs$Kcmid, "-")[[1]][2])))
      if(is.na(kcbmid)){
        kcbmid=mean(as.numeric(c(strsplit(kcs$Kcmid, ",")[[1]][1],strsplit(kcs$Kcmid, ",")[[1]][2])))
        
      }
    }
    
    if(is.na(kcbini)){
      kcbini=mean(as.numeric(c(strsplit(kcs$Kcbini, "-")[[1]][1],strsplit(kcs$Kcbini, "-")[[1]][2])))
      if(is.na(kcbini)){
        kcbini=mean(as.numeric(c(strsplit(kcs$Kcbini, ",")[[1]][1],strsplit(kcs$Kcbini, ",")[[1]][2])))
        
      }
    }
    
  }else{
    kcbini=kc[1]
    kcbmid=kc[2]
    kcbend=kc[3]  
  }
  
  #kcb=(day*kcb)/max(day)
  
  if(is.null(h)){
    h= singlekc[(grep((crop2),(singlekc$crop),ignore.case=TRUE)),] 
    if(nrow(h)==0){
      return (print(paste(" h not found for", crop2))) 
    }
    if(nrow(h)>1){
      h= h[(grep((regexpr),(h$crop),ignore.case=TRUE)),] 
    }
    h=mean(as.numeric(h$h))
  }
  ksingle=FALSE
  if(kctype=="single"){
    ksingle=TRUE
  }
  
  if(is.null(EF)){
    if(ksingle==TRUE){
      kcss= singlekc[(grep((crop2),(singlekc$crop),ignore.case=TRUE)),] 
      if(nrow(kcss)==0){
        #return (print(paste(" Kc not found for", crop2))) 
      }
      if(nrow(kcss)>1){
        kcss= kcss[(grep((regexpr),(kcss$crop),ignore.case=TRUE)),] 
      }
      
      if(is.null(kc)){
        
        kcini=mean(as.numeric(kcss$Kcini))
        kcmid=mean(as.numeric(kcss$Kcmid))
        kcend=mean(as.numeric(kcss$Kcend))
        if(is.na(kcend)){
          kcend=mean(as.numeric(c(strsplit(kcss$Kcend, "-")[[1]][1],strsplit(kcss$Kcend, "-")[[1]][2])))
          if(is.na(kcend)){
            kcend=mean(as.numeric(c(strsplit(kcss$Kcend, ",")[[1]][1],strsplit(kcss$Kcend, ",")[[1]][2])))
            
          }
        }
        
        if(is.na(kcmid)){
          kcmid=mean(as.numeric(c(strsplit(kcss$Kcmid, "-")[[1]][1],strsplit(kcss$Kcmid, "-")[[1]][2])))
          if(is.na(kcmid)){
            kcmid=mean(as.numeric(c(strsplit(kcss$Kcmid, ",")[[1]][1],strsplit(kcss$Kcmid, ",")[[1]][2])))
            
          }
        }
        
        if(is.na(kcini)){
          kcini=mean(as.numeric(c(strsplit(kcss$Kcini, "-")[[1]][1],strsplit(kcss$Kcini, "-")[[1]][2])))
          if(is.na(kcini)){
            kcini=mean(as.numeric(c(strsplit(kcss$Kcini, ",")[[1]][1],strsplit(kcss$Kcini, ",")[[1]][2])))
            
          }
        }
        #if(kcmid>0.45){
          #kcmid = kcmid + (0.04*(u2-2) - 0.004*(RHmin-45))* (h/3)^0.3
          #}
        
        #if(kcend>0.45){
          # kcend =kcend + (0.04*(u2-2) - 0.004*(RHmin-45))* (h/3)^0.3
          # }
        kcss=c(kcini,kcmid,kcend)
      }else{
        kcini=kc[1]
        kcmid=kc[2]
        kcend=kc[3]
        kcss=c(kcini,kcmid,kcend)
      }
      #kcsscum=cumsum(kcss)
    }
  }else{
    kcss=kc
  }

  if(!is.null(kcbmid)&&!is.null(kcbend)&&is.null(EF)){
    
    #if(kcbmid>0.45){
    # kcbmid = kcbmid + (0.04*(u2-2) - 0.004*(RHmin-45))* (h/3)^0.3
    # }
    
    #if(kcbend>0.45){
    # kcbend =kcbend + (0.04*(u2-2) - 0.004*(RHmin-45))* (h/3)^0.3
    # }
    kcss=c(kcbini,kcbmid,kcbend)
  }
  
  if(ksingle==FALSE&&is.null(EF)){
    
    kcss=c(kcbini,kcbmid,kcbend)
  }
  #print(kcss)
  
  #ky
  if(is.null(Ky)){ 
   ky1= ky[(grep((crop2),(ky$crop),ignore.case=TRUE)),]
    if(nrow(ky1)>0){
      ky1=ky1$ky
          Ky=mean(as.numeric(ky1))
      if(is.na(Ky)){
        Ky=mean(as.numeric(c(strsplit(ky1, "-")[[1]][1],strsplit(ky1, "-")[[1]][2])))

      }
          
    }    
  }
    #kcsscum=cumsum(kcss)
  
  ############Dual Kc#############
  Kcmin=Kcmin
  if(is.null(REW)){
  REW=mean(c(as.numeric(strsplit(soil$REW, "-")[[1]][1]),as.numeric(strsplit(soil$REW, "-")[[1]][2])))
  }
  #Maximum depth of water that can be evaporated from the topsoil change to 0.125
  if(length(Ze)==1){
    Ze=(1:max(day)/1:max(day))*Ze
  }
  
  
  if(is.null(FC)||is.null(WP)){
    if(is.null(FC)){
      FC=mean(c(as.numeric(strsplit(soil$FC, "-")[[1]][1]),as.numeric(strsplit(soil$FC, "-")[[1]][2])))
    }
    if(is.null(WP)){
      WP=mean(c(as.numeric(strsplit(soil$WP, "-")[[1]][1]),as.numeric(strsplit(soil$WP, "-")[[1]][2])))
    }
    
  }
  #print(kcss)
  hinter=NULL
  co=0
  c1=range*h
  a=range*max(day)
  
  if(length(h)==1){
  if(hmodel=="exponential"||hmodel=="exponent"||hmodel=="exp"){
    h=co+c1*(1-exp(-day/a))
    # h=(co+c1)%o%(1-exp(-day/a))
  }else{
  
  if(hmodel=="linear"||hmodel=="lin"||hmodel=="line"){
    
  h=(day*h)/max(day)
  
  }else{
    h=co+c1*(1-exp(-day^2/a^2))
  }
  }
  }else{
    hinter=TRUE
    if(hmodel=="exponential"||hmodel=="exponent"||hmodel=="exp"){
      
    h=(co+c1)%o%(1-exp(-day/a))
    h=t(h)
    }else{
      
      if(hmodel=="linear"||hmodel=="lin"||hmodel=="line"){
        
        h=(day%o%h)/max(day)
        h=t(h)
      }else{
        h=co+c1%o%(1-exp(-day^2/a^2))
        h=t(h)
      }
    } 
  }
  
  #FC=0.23
  #WP=0.10
  #TAW=1000*(FC-WP)*Zr
  #RAW=p*TAW

  outputdata=data.frame()
  
  outputdual=data.frame()
  
  
  REWcor=2.5
  TEWcor=10
  if(report=="stages"){
    report=stages2cum
  }
  if(report=="all"||report=="All"){
    report=day
  }
    report3=report
    Awcxy=NULL
    ETcxy=NULL
    TAWxy=NULL
    RAWxy=NULL
    Ineedxy=NULL
    time2=time
    if(length(time)!=max(day)){
      time2=rep(time[1],max(day))
    }
    fwproj=NULL
    fw2=fw
    fwp=NULL
    if(length(fw)==1){
      fw2=rep(fw,max(day))
      fwproj=TRUE
    }
    
    ETc_adjsum=0
    ETc_sum=0
    TSW_measured=NA
    AW_measured2=NA
    if(is.null(Zr)||is.null(p)||is.null(h)){
      return(print("These input parameters are needed Zr, p, and h "))
      
    }
    #check kcini
    thiskcini=NULL
  i=1
  #Adjust Kc for weather
  if(is.null(EF)){
  u2mid=mean(u2[stages2cum[2]:stages2cum[3]])
  RHminmid=mean(RHmin[stages2cum[2]:stages2cum[3]])
  hmid=mean(h[stages2cum[2]:stages2cum[3]])
  kcss[[2]]=kcss[[2]] + (0.04*(u2mid-2) - 0.004*(RHminmid-45))* (hmid/3)^0.3
  u2end=mean(u2[stages2cum[3]:stages2cum[4]])
  RHminend=mean(RHmin[stages2cum[3]:stages2cum[4]])
  hend=mean(h[stages2cum[3]:stages2cum[4]])
  kcss[[3]]=rastercon(kcss[[3]]>0.45,kcss[[3]] + (0.04*(u2end-2) - 0.004*(RHminend-45))* (hend/3)^0.3,kcss[[3]])
  #print(kcss)
  }
  #hmid=mean(h[stages2cum[2]:stages2cum[3]])
  #kcss[[2]]=kcss[[2]] + (0.04*(u2mid-2) - 0.004*(RHminmid-45))* (hmid/3)^0.3
  
 
  DEEP=FALSE
  timeDP=0
  #kcb2=data$Kcb
  if(is.null(EF)){
  if(!is.null(Zr3)){
    if(length(Zr3)==length(day)){
      Zr=Zr3
    }
  }
  
  if(!is.null(h3)){
    if(length(h3)==length(day)){
      h=h3
    }
  }
  }
  if(is.numeric(profile)){
    Zri=profile
  }else{
    if(profile=="fixed"){
      Zri=max(Zr)
    }else{
      Zri=Zr[1]
    }
  }
  
  if(!is.null(theta)){
    if(length(theta)==1){
     theta=c(theta,rep(NA,max(day)-1)) 
    }
    inittheta=theta[[1]]
  }else{
    inittheta=WP
  }
  if(!is.null(theta)){
  if(is.na(inittheta)){
    inittheta=WP 
  }
  }
  Drend=(1000*(FC-inittheta)*Zri*p[[1]])
  
  CR=0
  #fw=0.0001
 
  Dei=(1000*(FC-(0.5*WP))*Ze[1])
  AW_measured=NA
  AW_measuredxy=NULL
  pi=p
  #print(P)
  while(i<=length(day)){
    #if(!is.null(EF)){
    if(is.character(fw2[i])){
      fw= wetting[(grep((fw2[i]),(wetting$Wetting),ignore.case=TRUE)),]$fw[1]
    }else{
      fw=fw2[i]
    }
      if(!is.null(fwproj)){
        
      if(P[i]>0&&I[i]==0){
        fw=1 
        fwp=fw
      } 
        if(P[[i]]>0&&I[i]>0){
          total=P[[i]]+I[[i]]
          fw=((P[[i]]/total)*1)+((I[[i]]/total)*fw)
          fwp=fw
          
        }
      
      if(P[[i]]==0&&I[i]>0){
       fwp=NULL 
      }
      }
      fw=ifelse(!is.null(fwp),fwp,fw)
      #print(fw)
      
    #}
    #print(i)
    
    TEW=(1000*(FC-(0.5*WP))*Ze[i])
    #limits of h
    if(is.null(EF)){
    if(h[[i]]<0.1||h[[i]]>10){
      if(h[[i]]<0.1){
        h[[i]]=0.1 
      }
      if(h[[i]]>10){
        h[[i]]=10
      }
    }
    #limits of u2
    if(u2[[i]]<1||u2[[i]]>6){
      if(u2[[i]]<1){
        u2[[i]]=1 
      }
      if(u2[[i]]>6){
        u2[[i]]=6
      }
    }
    
    #limits of RHmin
    if((RHmin[[i]]<20||RHmin[[i]]>80)){
      if(RHmin[[i]]<20){
        RHmin[[i]]=20 
      }
      if(RHmin[[i]]>80){
        RHmin[[i]]=80
      }
    }
    }
    if(!is.null(hinter)&&!is.null(EF)&&!is.null(xydata)){
      hdata=cbind(xydata,h[i,])
      longitude=hdata[,1]
      latitude=hdata[,2]
      var=hdata[,3]
    hi= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
    #plot(hi)
    hi=rastercon(hi<0.1,0.1,rastercon(hi>10,10,hi))
    }else{
      hi=h[[i]]
    }
    
    if(!is.null(Zrinter)&&!is.null(EF)&&!is.null(xydata)){
      Zrdata=cbind(xydata,Zr[i,])
      longitude=Zrdata[,1]
      latitude=Zrdata[,2]
      var=Zrdata[,3]
      Zri= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
      #Zri=rastercon(Zri<0.1,0.1,rastercon(Zri>10,10,Zri))
    }else{
      if(is.numeric(profile)){
        Zri=profile
      }else{
        if(profile=="fixed"){
          Zri=max(Zr)
        }else{
          Zri=Zr[[i]]
        }
      }
     # Zri=Zr[[i]]
    }
    #do interpolation
    if(!is.null(EF)&&!is.null(RHinter)){
      if(nrow(xydata)==length(RHmin[i,])){
        Pdata=cbind(xydata,RHmin[i,])
        longitude=Pdata[,1]
        latitude=Pdata[,2]
        var=Pdata[,3]  
        RHmini= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
      }else{
        RHmini=mapmatrix(RHmin,i,map) 
      }
      
    }else{
      RHmini=RHmin[[i]]
    }
    
    if(!is.null(EF)&&!is.null(u2inter)){
      if(nrow(xydata)==length(u2[i,])){
        Pdata=cbind(xydata,u2[i,])
        longitude=Pdata[,1]
        latitude=Pdata[,2]
        var=Pdata[,3]  
        u2i= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
      }else{
        u2i=mapmatrix(u2,i,map)  
      }
    }else{
      u2i=u2[[i]]
    }
    
    if(!is.null(EF)&&!is.null(ETointer)){
      if(nrow(xydata)==length(ETo[i,])){
        Pdata=cbind(xydata,ETo[i,])
        longitude=Pdata[,1]
        latitude=Pdata[,2]
        var=Pdata[,3]  
        EToi= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
      }else{
        EToi=mapmatrix(ETo,i,map)  
      }
    }else{
      EToi=ETo[[i]]
    }
    if(!is.null(EF)&&!is.null(thetainter)){
      if(nrow(xydata)==length(theta[i,])){
        Pdata=cbind(xydata,theta[i,])
        longitude=Pdata[,1]
        latitude=Pdata[,2]
        var=Pdata[,3]  
        thetai= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
      }else{
        thetai=mapmatrix(theta,i,map)  
      }
    }else{
      thetai=theta[[i]]
    }
    if(!is.null(EF)&&!is.null(Iinter)){
      if(nrow(xydata)==length(I[i,])){
        Pdata=cbind(xydata,I[i,])
        longitude=Pdata[,1]
        latitude=Pdata[,2]
        var=Pdata[,3]  
        I= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
      }else{
        Ii=mapmatrix(I,i,map)  
      }
    }else{
      Ii=I[[i]]
    }
    #this=as.numeric(lapply(lists,function(x)if(i<x){x*2}else{0}))
    if(!is.null(allEF)){
      Kci=kc[[i]]
    }else{
    if(i<=stages2cum[1]){
      #kcini
      #print(stages2[1])
      
      if(ksingle==TRUE&&is.null(EF)){
        if(nrow(data)>stages2[1]){
          initdata=data[1:stages2[1],]
        }else{
          initdata=data
         # print(length(data))
          
        }
       # print(initdata)
        
        initdataI=initdata[initdata$I>0.2*ETo[[i]],]
        if(is.null(initdata[["P"]])){
          initdataP=initdata[initdata$P.RO>0.2*ETo[[i]],]
          rain="P.RO"
        }else{
          initdataP=initdata[initdata[["P"]]>0.2*ETo[[i]],]
          rain="P"
        }
        nw=nrow(initdataP)+nrow(initdataI)
        #print(nw)
        thisETo=mean(initdata$ETo)
        Eso=1.15*thisETo
        Pmean=mean(initdataP[[rain]])
        Imean=mean(initdataI$I)
        infiltration=Pmean+Imean
        if(is.na(Imean)&&!is.na(Pmean)){
          infiltration=Pmean
        }
        if(is.na(Pmean)&&!is.na(Imean)){
          infiltration=Imean
        }
        if(!is.na(Pmean)&&!is.na(Imean)){
          infiltration=mean(initdata$I)+mean(initdata$P)
          #print("YES")
        }
       if(is.na(infiltration)){
         infiltration=9
       }
        if(infiltration<=10){
          TEWcor=10
          REWcor=min(max(2.5,6/(thisETo^0.5)),7)
        }else{
          soil3=toupper(soil2)
          if(infiltration>=40&&(soil3=="SAND"||soil3=="LOAMY SAND")){
            TEWcor = min(15, 7 *(thisETo^0.5))
            REWcor = min(6, TEWcor - 0.01)
          }else{
            if(infiltration>=40&&(soil3=="LOAM"||soil3=="SILT"||soil3=="SANDY LOAM"
              ||soil3=="SILT LOAM"||soil3=="CLAY LOAM"||soil3=="SILTY CLAY"||soil3=="CLAY")){
              TEWcor = min(28, 13 *(thisETo^0.5))
              REWcor = min(9, TEWcor - 0.01)
            }
            
            
          }
        }
        
        if(infiltration>10&&infiltration<40){
          TEWcor=10
          REWcor=min(max(2.5,6/(thisETo^0.5)),7)
          t1 = REWcor / Eso
          tw=stages2[1]/(nw+0.5)
          if(tw>t1){
            kcss10=(TEWcor-(TEWcor-REWcor)*exp((-(tw-t1)*Eso*(1+((REWcor)/(TEWcor-REWcor))))/(TEWcor)))/(tw*thisETo)
            #kcss[1]=(TEWcor-(TEW-REW)*exp((-(tw-t1)*Eso*(1+((REW)/(TEW-REW))))/(TEW)))/(tw*ETo[[i]])
          }else{
            kcss10=Eso/thisETo
            
          }
          if((infiltration>10&&infiltration<40)&&(soil3=="SAND"||soil3=="LOAMY SAND")){
            TEWcor = min(15, 7 *(thisETo^0.5))
            REWcor = min(6, TEWcor - 0.01)
          }
          if((infiltration>10&&infiltration<40)&&(soil3=="LOAM"||soil3=="SILT"||soil3=="SANDY LOAM"
          ||soil3=="SILT LOAM"||soil3=="CLAY LOAM"||soil3=="SILTY CLAY"||soil3=="CLAY")){
            TEWcor = min(28, 13 *(thisETo^0.5))
            REWcor = min(9, TEWcor - 0.01)
          }
          t1 = REWcor / Eso
          if(tw>t1){
            kcss40=(TEWcor-(TEWcor-REWcor)*exp((-(tw-t1)*Eso*(1+((REWcor)/(TEWcor-REWcor))))/(TEWcor)))/(tw*thisETo)
            #kcss[1]=(TEWcor-(TEW-REW)*exp((-(tw-t1)*Eso*(1+((REW)/(TEW-REW))))/(TEW)))/(tw*ETo[[i]])
            
          }else{
            kcss40=Eso/thisETo
          }
          
          kcss[1]=kcss10+((infiltration-10)/(40-10))*(kcss40-kcss10)
          #print(cbind(day=i,kcss10=kcss10,kcss40=kcss40,kcss=kcss[1]))
          
        }else{
        t1 = REWcor / Eso
        tw=stages2[1]/(nw+0.5)
        if(tw>t1){
          kcss[1]=(TEWcor-(TEWcor-REWcor)*exp((-(tw-t1)*Eso*(1+((REWcor)/(TEWcor-REWcor))))/(TEWcor)))/(tw*thisETo)
          #kcss[1]=(TEWcor-(TEW-REW)*exp((-(tw-t1)*Eso*(1+((REW)/(TEW-REW))))/(TEW)))/(tw*ETo[[i]])
          
        }else{
          kcss[1]=Eso/thisETo
          
          }
        }
        
      }
     # print(cbind(REWcor=REWcor,TEWcor=TEWcor,t1=t1,tw=tw,Eso=Eso,Kc=kcss[[1]]))
      Kci=  kcss[[1]]*fw
      
     # print(fw)
      
      #print(c(nw=nw,t1=t1,tw=tw,REWcor=REWcor,TEWcor=TEWcor,REW=REW,Kc=kcss[1]))
      
    }else{
      if(i>stages2cum[1]&&i<=stages2cum[2]){
        Kci=kcss[[1]]+((i-stages2[1])/stages2[2])*(kcss[[2]]-kcss[[1]])
      }else{
        if(i>stages2cum[2]&&i<=stages2cum[3]){
          #mid
          #Kci=kcss[2]+((i-stages2cum[2])/stages2[3])*(kcss[4]-kcss[2])
          #if(is.null(EF)){
            #kcss[[2]]=kcss[[2]] + (0.04*(u2mid-2) - 0.004*(RHminmid-45))* (hi/3)^0.3
           # }
          Kci=kcss[[2]]
          #Kci=rastercon(Kci>0.45,Kci + (0.04*(u2i-2) - 0.004*(RHmini-45))* (hi/3)^0.3,Kci)
        }else{
          #late
          #if(is.null(EF)){
          #kcss[[3]]=rastercon(kcss[[3]]>0.45,kcss[[3]] + (0.04*(u2end-2) - 0.004*(RHminend-45))* (hi/3)^0.3,kcss[[3]])
          #}
          Kci=kcss[[2]]+((i-stages2cum[3])/stages2[4])*(kcss[[3]]-kcss[[2]])
          #print("yes")
        }
      }
    }
    }
    #print(Kci)
    if(!is.null(kc3)&&is.null(EF)){
      if(length(kc3)==length(day)){
        Kci=kc3[i]
      }
    }
    
    if(is.null(EF)){
      kcb=max(Kcmin,Kci)
      #kcb=Kci
      Kcmax = max(1.2 + (0.04*(u2[i]-2) - 0.004*(RHmini-45))* (h[[i]]/3)^0.3,(kcb+0.05))
      
      if(is.null(fc2)){
        fc=((kcb-Kcmin)/(max(Kcmax-Kcmin,0.01)))^(1+(0.5*h[[i]]) )
        #print(fc)
      }else{
        fc=fc2[[i]]
      }
      if(!is.null(latitude2)&&density2!="full"&&!is.null(DOY2)){
        if((is.null(latitude2)||is.null(DOY2))&&(nongrowing=="weeds"||nongrowing=="grass")){
          return (print("latitude and DOY needed"))
        }
        
        
        
        
       
          
    if((i>=stages2cum[2]&&i<=stages2cum[4])||(nongrowing=="weeds"||nongrowing=="grass")){
      
      if(density!="full"){
        kcbh=min(1.0+(0.1*max(h)),1.20)
        kcbfull=kcbh + (0.04*(u2i-2) - 0.004*(RHmini-45))* (max(h)/3)^0.3
        latituderad=latitude[i]*(pi/180)
        declination=0.409*sin((((2*pi*DOY[i])/(365)))-1.39)
        sinethisBeta=sin(latituderad)*sin(declination)+cos(latituderad)*cos(declination)
        thisBeta=asin(sinethisBeta)
        if(density=="row"||density=="Row"||density=="ROW"){
          fceff=min(fc*(1+(HWR[i]/tan(thisBeta))),1)
        }
        if(density=="random"||density=="round"||density=="Random"||density=="RANDOM"){
          fceff=min(fc/sin(thisBeta),1)
        }
        Kd=min(1,2*fc,fceff^(1/(1+(max(h)))))
        #Kcmin=kcss[1]
        Kci=Kcmin+Kd*(kcbfull-Kcmin)
        #print(Kci)
        }
      
    }
      }
    }
    if(!is.null(slope2)&&!is.null(r12)&&!is.null(y2)){
      Fr=(slope[i]+y[i]*(1+0.34*(u2i)))/(slope[i]+y[i]*(1+0.34*(u2i)*(r1[i]/100)))
    Kci=Kci*Fr
    }
    if(!is.null(nongrowing)){
    if(nongrowing=="Bare"||nongrowing=="bare"||nongrowing=="Bare Soil"||nongrowing=="bare soil"||
       nongrowing=="mulch"||nongrowing=="Mulch"||nongrowing=="MULCH"){
    Kci=kcss[1]  
    if(nongrowing=="mulch"||nongrowing=="Mulch"||nongrowing=="MULCH"){
      Kci=Kci*0.95 
      #
    }
    
    }
    }
    if(!is.null(nongrowing)){
      if(kctype!="single"&&(nongrowing=="Bare"||nongrowing=="bare"||nongrowing=="Bare Soil"||nongrowing=="bare soil"
                          ||nongrowing=="mulch"||nongrowing=="Mulch"||nongrowing=="MULCH")){
        
       if(nongrowing=="mulch"||nongrowing=="Mulch"||nongrowing=="MULCH"){
         fc=0
         fw=1
         TEW=TEW*0.95
         
       }else{
         Kci=0
        Kcmin=0
        
       }
      }
      
    }
    kcb=Kci
    
    #print(kcb)
    
    
    if(!is.null(EF)&&!is.null(Pinter)){
      if(nrow(xydata)==length(P[i,])){
        #print("YES")
        Pdata=cbind(xydata,P[i,])
        longitude=Pdata[,1]
        latitude=Pdata[,2]
        var=Pdata[,3]  
        Prec= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
        #print(i)
      }else{
        Prec=mapmatrix(P,i,map)  
      }
    }else{
      Prec=P[[i]]
    }
    if(is.numeric(rc)){
      Runoff=rc*Prec
    }else{
      CNI=CNII/(2.281-(0.01281*CNII))
      #CNI=(4.2*CNII)/(10-(0.058*CNII))
      CNIII=CNII/(0.427+(0.00573*CNII))
      #CNIII=(23*CNII)/(10+(0.13*CNII))
      DeAWCIII=0.5*REW
      DeAWCI=(0.7*REW)+(0.3*TEW)
      if(is.null(EF)){
        if(Dei<=DeAWCIII){
          CN= CNIII 
        }else{
          if(Dei<=DeAWCI){
            CN= CNI
          }else{
            if(Dei<DeAWCI&&Dei>DeAWCIII){
              CN=(((Dei-(0.5*REW))*CNI)+(((0.7*REW)+(0.3*TEW)-Dei)*CNIII))/((0.2*REW)+(0.3*TEW))  
            }else{
              CN=CNII
            }  
          }  
        }
      }else{
        CN=rastercon(Dei<=DeAWCIII,CNIII,-1000)
        #CN=rastercon(Dei<DeAWCI&&Dei>DeAWCIII,(((Dei-(0.5*REW))*CNI)+(((0.7*REW)+(0.3*TEW)-Dei)*CNIII))/((0.2*REW)+(0.3*TEW)),CN)
        CN=rastercon(Dei<DeAWCI,(((Dei-(0.5*REW))*CNI)+(((0.7*REW)+(0.3*TEW)-Dei)*CNIII))/((0.2*REW)+(0.3*TEW)),CN)
        CN=rastercon(Dei>DeAWCIII,(((Dei-(0.5*REW))*CNI)+(((0.7*REW)+(0.3*TEW)-Dei)*CNIII))/((0.2*REW)+(0.3*TEW)),CN)
        
        CN=rastercon(CN==-1000,CNII,CN)
        
        #CN=rastercon(Dei<=DeAWCIII,CNIII,rastercon(Dei<=DeAWCI,CNI,rastercon(Dei<DeAWCI&&Dei>DeAWCIII,(((Dei-(0.5*REW))*CNI)+(((0.7*REW)+(0.3*TEW)-Dei)*CNIII))/((0.2*REW)+(0.3*TEW)))),CNII)
      }
      #}
      #print(CN)
      
      S=(1000/CN) -10
      #S05=1.33*S^1.15
      Ia=0.05*S
      #Ia=0.2*S
      Pin=0.0393701*Prec
      #Pin=0.0393701*Pin
      #Runoff=rastercon(Pin<Ia,0,((Pin-Ia)^2)/((Pin+(0.95*S05)))
      #Runoff=rastercon(Pin<Ia,0,((Pin-Ia))^2/((Pin+0.95*S05)))*25.4
      Runoff=rastercon(Pin<Ia,0,((Pin-Ia)^2)/((Pin-Ia)+S))*25.4
      CNII=CN
    }
    #kc###################
    
    
    #print((I[i]/fw) - (P[i]-Runoff))
    if(is.null(EF)){
      #print(fw)
      Deistart=ifelse(strptime(time2[i], "%H:%M") <= strptime("12:00", "%H:%M"), max(Dei - (I[i]/fw) - (P[i]-Runoff),  0), Dei)
      #thisDr=ifelse(strptime(time2[i], "%H:%M") <= strptime("12:00", "%H:%M"), Drend , Drstart)
    }
    
    
   # Deistart=ifelse(strptime(time, "%H:%M") <= strptime("12:00", "%H:%M"), max(Dei - I[i]/fw - (P[i]-Runoff),  0), Dei)
    
    #Deistart=data$Dstart[[i]]
    #Dei=Dei-(P-Run)
    
      #print(c(Kcmax,kcb,fc))
      #if(is.na(fc)&&kcb<=0.15){
       # kcb=0.2
       # fc=((kcb-Kcmin)/(Kcmax-Kcmin))^(1+(0.5*h[[i]]) ) 
        
        # return(print(paste("fc cannot be calculated. Please choose a higher kcbini")))
    #}#
    
    if(is.null(EF)){
      
      few=min(1-fc,fw)

      if(fw2[i]=="drip"||fw2[i]=="DRIP"||fw2[i]=="Drip"){
        few=min(1-fc,(1-((2/3)*fc))*fw)
      }
      if(Deistart<=REW){
        Kr=1
      }else{
        Kr=(TEW-Deistart)/(TEW-REW)
      }
      Ke=max(min(Kr*(Kcmax-kcb),few*Kcmax),0)
      
      E=max(Ke*ETo[[i]],0)
      
    }
    #Kr=max((TEW-Deistart)/(TEW-REW),0)
    
    #runoff
    #Runoff=(P[[i]]*rc)
    #if(P[[i]]>0){
    
    
    #water loss out of the root zone by deep percolation on day i
    #print(c(Kr,Kcmin,kcb,few,Kcmax,Ke))
    
    
    #Dei=max()
    
    #Total available soil water
    TAW=1000*(FC-WP)*Zri
    #print(Zri)
    
    #readily available soil water
    RAW=pi*TAW
    
    #DP=max(P.RO[[i]]+I[[i]]-Drend,0)
    #root zone depletion at the end of day i
    
    
    
    if(is.null(EF)){
      Drstart=ifelse(strptime(time2[i], "%H:%M") <= strptime("12:00", "%H:%M"), max(Drend - I[i] - (P[i]-Runoff),  0), Drend)
      #thisDr=ifelse(strptime(time2[i], "%H:%M") <= strptime("12:00", "%H:%M"), Drend , Drstart)
    }else{
      
      #Deistart=ifelse(strptime(time, "%H:%M") <= strptime("12:00", "%H:%M"), max(Dei - I[[i]]/fw - (P[[i]]-Runoff),  0), Dei)
      thisDrend=Drend - I[i] - (Prec-Runoff)
      thisDrend[thisDrend<0]=0
      #print(P[i])
      
      #print(thisDrend)
      if(i>1){
       # Drstart=ifelse(strptime(time2[i], "%H:%M") <= strptime("12:00", "%H:%M"), thisDrend, Drend)
      if(strptime(time2[i], "%H:%M") <= strptime("12:00", "%H:%M")){
        Drstart=thisDrend
      }else{
        Drstart=Drend 
      }
        }else{
        Drstart=Drend
        
      }
      
      
    }
    
   # print(cbind(Drstart,Drend))

    
    if(is.null(EF)){
      if(Drstart<RAW){
        Ks=1 
      }else{
        Ks=max((TAW-Drstart)/((1-pi)*TAW),0)
      }
    }
    
    if(!is.null(EF)&&!is.null(Rninter)){
      if(nrow(xydata)==length(Rn[i,])){
        Pdata=cbind(xydata,Rn[i,])
        longitude=Pdata[,1]
        latitude=Pdata[,2]
        var=Pdata[,3]  
        Rni= invdist(longitude=longitude,latitude=latitude,var=var,ext=map) 
        #print(i)
      }else{
        Rni=mapmatrix(Rn,i,map)  
      }
    }else{
      Rni=Rn[[i]]
    }
    #print(Ks)
    #Ks=1
    if(is.null(EF)){
      if(ksingle==FALSE){
        kcb=Kci
        Kc=((Ks*kcb)+Ke)
        Kci1=(kcb+Ke)
      }else{
        Kci1=Kci
        kcb=Kci
        Kc=kcb*Ks  
      }
      
      Kc=max(Kcmin,Kc)
      #print(kcss)
      ETc=Kci1*EToi
      ETc2=Kc*EToi
      #ETc_adjust=ETc2
      pi=p+0.04*(5-ETc)
      
    }else{
      #print(Kci)
      if(kctype=="METRIC"||kctype=="metric"||kctype=="SSEB"||kctype=="sseb"){
        ETc2=Kci*ETo[[i]]
      }else{
        ETc2=(Kci*Rni)*0.035
          }
      pi=p+0.04*(5-ETo[i])
    }
    
    if(is.null(EF)){
      pi=max(0.1,pi)
      pi=min(pi,0.8)
    }else{
     # pi=maxValue(0.1,pi)
      pi[pi<0.1]=0.1
      pi[pi>0.8]=0.8
      
      #pi=minValue(pi,0.8)
    }
    #pi=max(0.1,pi)
    #pi=min(pi,0.8)
    
    #pi=0.6
    thisDr=rastercon(strptime(time2[i], "%H:%M") <= strptime("12:00", "%H:%M"), Drend , Drstart)
    
    
    DP=max((Prec-Runoff)+Ii-ETc2-thisDr,0)
    #plot(ETc2)
    #Dei=data$Dei[[i]]
    #Dei=max(Dei-(P[[i]]-Runoff)-(I[[i]]/fw)+(E/few)+Tewi+DP,0)
    ##capillary rise
    #initgw=(DP/1000)+initgw
    if(is.null(a1)){
    a1=FC*Zri*1000
    }
    #b1=-0.17
    if(is.null(a2)){
    a2=1.1*((FC+WP)/2)*Zri
    }
    #b2=-0.27
    #a3=-1.3
    soil4=toupper(soil2)
    if(soil4=="CLAY"||soil4=="SILTY CLAY"||soil4=="LOAM"){
      if(is.null(b3)){
      b3=6.7
      }
    }else{
      if(is.null(b3)){
      b3=6.2
      }
    }
    if(soil4=="SILTY LOAM"||soil4=="SILTY CLAY"||soil4=="LOAM"){
      if(is.null(a4)){
      a4=4.6
      }
    }else{
      if(is.null(a4)){
        a4=4.6
      }
    }
    
    if(soil4=="SILTY"||soil4=="LOAM"){
      if(is.null(b4)){
        b4=-0.65
      }
    }else{
      if(is.null(b4)){
        b4=-2.5
      }
    }
    Wc=a1*initgw^b1
    if(is.null(EF)){
      if(initgw>=3){
        Ws= a2*initgw^b2
      }else{
        Ws=240  #mm
      }
    }else{
      Ws= rastercon(initgw>=3,a2*initgw^b2,240)
    }
    ETm=ETc2
    if(is.null(EF)){
      if(ETm<=4){
        Dwc=(a3*ETm)+b3 
        k=1-exp(-0.61*LAI)
      }else{
        Dwc=1.4 #m
        k=3.8*ETm
      }
    }else{
      Dwc=rastercon(ETm<=4,(a3*ETm)+b3,1.4)  
      k=rastercon(ETm<=4,1-exp(-0.61*LAI),3.8/ETm) 
    }
    if(is.null(EF)){
      if(initgw<=Dwc){
        CRmax=k*ETm
      }else{
        CRmax=a4*initgw^b4
      }
    }else{
      CRmax=rastercon(initgw<=Dwc,k*ETm,a4*initgw^b4)
    }
    a=((Ws+FC)/2)*Zri*1000
    if(soil4=="SAND"||soil4=="SANDY LOAM"||soil4=="LOAMY SAND")
    {
      b=-0.0174
    }else{
      b=-0.0172
    }
    if(is.null(EF)){
      if(DP>0){
        timeDP=1
        DEEP=TRUE
      }else{
        if(DEEP==TRUE){
          timeDP=timeDP+1
        }
      }
    }else{
      timeDP=rastercon(DP>0,1,timeDP+1) 
      DEEP=rastercon(DP>0,TRUE,FALSE)  
      
    }
    #if(timeDP>0){
    Wa=a*timeDP^b
    #print(timeDP)
    #}
    if(is.null(EF)){
      if(Wa<Ws){
        CR=CRmax
      }else{
        if(Wa>=Ws&&Wa<=Wc){
          CR=CRmax*((Wc-Wa)/(Wc-Ws)) 
        }else{
          CR=0
        }
      }
    }else{
      CR=rastercon(Wa<Ws,CRmax,CRmax*((Wc-Wa)/(Wc-Ws)))
      CR=rastercon(Wa<Wc,CRmax*((Wc-Wa)/(Wc-Ws)),CR)
      
      #CR=rastercon(Wa>=Ws&Wa<=Wc,CRmax*((Wc-Wa)/(Wc-Ws)),0)
      #CR[Wa>=Ws&Wa<=Wc,]=CRmax*((Wc-Wa)/(Wc-Ws))
    }
    CR=max(CR,0)
    DPei=(Prec-Runoff)+(Ii/fw)-Dei
    if(is.null(EF)){
      DPei=max(DPei,0)
      kti=((1-(Dei/TEW))/(1-(Drend/TAW)))*(Ze[[i]]/Zr[[i]])^0.6
      Tewi=kti*kcb*Ks*ETo[[i]]
      Dei=max(Dei-((Prec-Runoff))-(Ii/fw)+(E/few)+(Tewi)+DPei,0)
    }else{
      DPei[DPei<0]=0
      Dei=max(Dei-((Prec-Runoff))-(Ii/fw)+ETc2+DPei,0)
    }
    
    #Drend=Drend-(Prec-Runoff)-Ii-CR+ETc2+DP
    if(!is.null(EF)){
      Drend=Drend-(Prec-Runoff)-Ii-CR+ETc2+DP
      Drend[Drend<0,]=0
      Drend[Drend>TAW,]=TAW
    }else{
      Drend=max(Drend-(Prec-Runoff)-Ii-CR+ETc2+DP,0)
      Drend=min(Drend,TAW)
      
    }
    
    #print(c(P[i],Runoff,I[i]))
    #print(paste("P",Prec))
    #print(paste("Runoff",Runoff))
    #print(paste("I",Ii))
    
    
    
    #print(CR)
    #print(ETc2)
    #print(DP)
    
    initgw=initgw-(DP/1000)+(CR/1000)+((DP*Kb)/1000)
    #AW=TAW-Drend
    mod_theta=FC-(Drstart/(Zri*1000))
    if(!is.null(EF)){
      mod_theta=FC-(Drend/(Zri*1000))
    }
    if(!is.null(EF)&&!is.null(theta)){
    mod_theta[mod_theta<0,]=0
    mod_theta[mod_theta>1,]=1
    }
    
    TSW=1000*(mod_theta)*Zri
    #print(mod_theta)
    AW=1000*(mod_theta-WP)*Zri
    #print(AW)
    
    if(!is.null(theta)){
    AW_measured=1000*(thetai-WP)*Zri
    TSW_measured=1000*(thetai)*Zri
    
    }
    #print(TSW)
    Ineed=rastercon(Drstart>=RAW,Drstart,0)
    if(!is.null(EF)){
      Ineed=rastercon(Drend>=RAW,Drend,0)
    }
    ETc_adjsum=ETc_adjsum+ETc2
    ETc_sum=ETc_sum+ETo[i]
    #print(EF)
    
    #print(RAW)
    #print(AW_measured)
    print (paste("Day", i,"of", max(day)))
    Ineed[Ineed<=0]=0
    if(!is.null(EF)){
     # print(day)
      
      #i2= format(i, nsmall = 200000000)
      #reportF=format(report, nsmall = 200000000)
      if(length(report[report==i])>0&&is.null(report2)){
        thisETc2=NULL
        ETa1=ETc2
        if(class(ETa1)!="RasterLayer"){
         thisETc2=kcss[[1]] 
         thisETc2[thisETc2>-1000,]=ETa1
         ETa1=thisETc2
         #rm(thisETc2)
        }
        Kc1=Kci
        Drend2= Drend 
        if(class(Drend2)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=Drend2
          Drend2=thisETc2
        }
        AW2=AW
        if(class(AW2)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=AW2
          AW2=thisETc2
          }
        TAW2=TAW
        if(class(TAW2)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=TAW2
          TAW2=thisETc2
        }
        RAW2=RAW
        if(class(RAW2)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=RAW2
          RAW2=thisETc2
        }
        Ineed2=Ineed
        if(class(Ineed2)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=Ineed2
          Ineed2=thisETc2
        }
        Runoff2=Runoff
        if(class(Runoff2)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=Runoff2
          Runoff2=thisETc2
        }
        CR2=CR
        if(class(CR2)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=CR2
          CR2=thisETc2
        }
        DP2=DP
        if(class(DP2)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=DP2
          DP2=thisETc2
        }
        gw=initgw
        if(class(gw)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=gw
          gw=thisETc2
        }
        report2=TRUE
        #print(Runoff2)
        
        #report[i]=NA
        if(!is.null(theta)){
        AW_measured2=AW_measured
        if(class(AW_measured2)!="RasterLayer"){
          thisETc2=Kci 
          thisETc2[thisETc2>-1000,]=AW_measured2
          AW_measured2=thisETc2
        }
        }
        
        rm(thisETc2)
      }else{
      if(length(report[report==i])>0&&!is.null(report2)){
        #print(report3)
        #print(Ineed2)
        
        ETa1=stack(ETa1,ETc2)
        Kc1=stack(Kc1,Kci)
        Drend2= stack(Drend2,Drend)
        AW2=stack(AW2,AW)
        
        #TAW2=stack(TAW2,TAW)
        #RAW2=stack(RAW2,RAW)
        
        Ineed2=stack(Ineed2,Ineed)
        
        if(class(Runoff)=="RasterLayer"){
          Runoff2=stack(Runoff2,Runoff)
        }else{
          Runoff2=rbind(Runoff2,Runoff)
        }
        
        if(class(CR)=="RasterLayer"){
          CR2=stack(CR2,CR)
        }else{
          CR2=rbind(CR2,CR)
        }
        
        if(class(DP)=="RasterLayer"){
          DP2=stack(DP2,DP)
        }else{
          DP2=rbind(DP2,DP)
        }
        #print(report[i])
        
        if(class(initgw)=="RasterLayer"){
          #if(class(gw)!="RasterLayer"){
            #gw=initgw
          #}
          #print(initgw)
          gw=stack(gw,initgw)
        }else{
          gw=rbind(gw,initgw)
        }
        #print(AW_measured)
        if(!is.null(theta)){
          
        if(class(AW_measured)=="RasterLayer"){
          AW_measured2=stack(AW_measured2,AW_measured)
        
        }else{
          #print("yess")
          
          AW_measured2=rbind(AW_measured2,AW_measured)
        }
        }
        if(!is.null(theta)&&length(report3)==daylenths){
          report[i]=NA
          
        }
        
      }
        #print(i)
       # print(report[i])
        
      }
      #print("KKKKKKKKKKKKKKKKKKKKK")
      runofftest=NULL
    if(i==1){
      if(!is.null(xydata)){
        #print("KKKKKKKKKKKKKKKKKKKKK")
        
          ETcxy=xyextract(ETc2,xydata[,1],xydata[,2])
          #print(Kcxy)
          
          colnames(ETcxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          Kcxy=xyextract(Kci,xydata[,1],xydata[,2])
          
          #print("kci......")
          colnames(Kcxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          Drendxy=xyextract(Drend,xydata[,1],xydata[,2])
          #print("Drend......")
          colnames(Drendxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          
          if(class(AW)=="RasterLayer"){
          Awcxy=xyextract(AW,xydata[,1],xydata[,2])
          colnames(Awcxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          }else{
            Awcxy= cbind(xydata[,1],xydata[,2],AW)
            colnames(Awcxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          }
          
          #TAWxy=xyextract(TAW,xydata[,1],xydata[,2])
          #colnames(TAWcxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          #RAWxy=xyextract(RAW,xydata[,1],xydata[,2])
          #colnames(RAWcxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          if(class(Ineed)=="RasterLayer"){
            Ineedxy=xyextract(Ineed,xydata[,1],xydata[,2])
          colnames(Ineedxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          }else{
            Ineedxy= cbind(xydata[,1],xydata[,2],Ineed)
            colnames(Ineedxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          }
          if(class(Runoff)=="RasterLayer"){
            Runoff2xy=xyextract(Runoff,xydata[,1],xydata[,2])
            colnames(Runoff2xy)=c("longitude", "latitude",paste("day",day[i],sep=""))}else{
            Runoff2xy= cbind(xydata[,1],xydata[,2],Runoff)
            colnames(Runoff2xy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          }
          if(class(CR)=="RasterLayer"){
            CR2xy=xyextract(CR,xydata[,1],xydata[,2])
            colnames(CR2xy)=c("longitude", "latitude",paste("day",day[i],sep=""))}else{
              CR2xy= cbind(xydata[,1],xydata[,2],CR)
              colnames(CR2xy)=c("longitude", "latitude",paste("day",day[i],sep=""))
            }
          
          if(class(DP)=="RasterLayer"){
            DP2xy=xyextract(DP,xydata[,1],xydata[,2])
            colnames(DP2xy)=c("longitude", "latitude",paste("day",day[i],sep=""))}else{
              DP2xy= cbind(xydata[,1],xydata[,2],DP)
              colnames(DP2xy)=c("longitude", "latitude",paste("day",day[i],sep=""))
            }
          
          if(class(initgw)=="RasterLayer"){
            gwxy=xyextract(initgw,xydata[,1],xydata[,2])
            colnames(gwxy)=c("longitude", "latitude",paste("day",day[i],sep=""))}else{
              gwxy= cbind(xydata[,1],xydata[,2],initgw)
              colnames(gwxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
            }
          
          if(!is.null(theta)){
            
          if(class(AW_measured)=="RasterLayer"){
            AW_measuredxy=xyextract(AW_measured,xydata[,1],xydata[,2])
            colnames(AW_measuredxy)=c("longitude", "latitude",paste("day",day[i],sep=""))}else{
              AW_measuredxy= cbind(xydata[,1],xydata[,2],AW_measured)
              colnames(AW_measuredxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
            }
          }
          
      }
      
      #RO=Runoff
    }else{
      
      if(!is.null(xydata)){
        thisextract=xyextract(ETc2,xydata[,1],xydata[,2])
        thisextract=thisextract[3]
        colnames(thisextract)=paste("day",day[i],sep="")
        ETcxy=cbind(ETcxy,thisextract)
        thisextract=xyextract(AW,xydata[,1],xydata[,2])
        thisextract=thisextract[3]
        colnames(thisextract)=paste("day",day[i],sep="")
        Awcxy=cbind(Awcxy,thisextract)
        #print(Awcxy)
        thisextract=xyextract(Kci,xydata[,1],xydata[,2])
        thisextract=thisextract[3]
        colnames(thisextract)=paste("day",day[i],sep="")
        Kcxy=cbind(Kcxy,thisextract)
        
        thisextract=xyextract(Drend,xydata[,1],xydata[,2])
        thisextract=thisextract[3]
        colnames(thisextract)=paste("day",day[i],sep="")
        Drendxy=cbind(Drendxy,thisextract)
        
        
        thisextract=xyextract(Ineed,xydata[,1],xydata[,2])
        thisextract=thisextract[3]
        colnames(thisextract)=paste("day",day[i],sep="")
        Ineedxy=cbind(Ineedxy,thisextract)
        #print("okkkkkkkkkkkk2")
        
        
        
        if(class(initgw)=="RasterLayer"){
          thisextract=xyextract(initgw,xydata[,1],xydata[,2])
          thisextract=thisextract[3]
          colnames(thisextract)=paste("day",day[i],sep="")
          gwxy=cbind(gwxy,thisextract)}else{
            gwxy=cbind(gwxy,initgw)
            #colnames(gwxy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          }
        
        if(class(Runoff)=="RasterLayer"){
          thisextract=xyextract(Runoff,xydata[,1],xydata[,2])
          thisextract=thisextract[3]
          colnames(thisextract)=paste("day",day[i],sep="")
          runofftest=NULL
          Runoff2xy=cbind(Runoff2xy,thisextract)}else{
            runofftest=TRUE
            Runoff2xy=cbind(Runoff2xy,Runoff)
            #colnames(Runoff2xy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          }
        if(class(CR)=="RasterLayer"){
          thisextract=xyextract(CR,xydata[,1],xydata[,2])
          thisextract=thisextract[3]
          colnames(thisextract)=paste("day",day[i],sep="")
          CR2xy=cbind(CR2xy,thisextract)}else{
            CR2xy=cbind(CR2xy,CR)            #colnames(CR2xy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          }
        
        if(class(DP)=="RasterLayer"){
          thisextract=xyextract(DP,xydata[,1],xydata[,2])
          thisextract=thisextract[3]
          colnames(thisextract)=paste("day",day[i],sep="")
          DP2xy=cbind(DP2xy,thisextract) }else{
            DP2xy=cbind(DP2xy,DP)            #colnames(DP2xy)=c("longitude", "latitude",paste("day",day[i],sep=""))
          }
        if(!is.null(theta)){
        if(class(AW_measured)=="RasterLayer"){
          thisextract=xyextract(AW_measured,xydata[,1],xydata[,2])
          thisextract=thisextract[3]
          colnames(thisextract)=paste("day",day[i],sep="")
          AW_measuredxy=cbind(AW_measuredxy,thisextract) 
          }else{
            AW_measuredxy=cbind(AW_measuredxy,AW_measured) 
          }
        }
        #print("yess")
      
      }
      
    }
    }
    #print(kcb)
    
    #plot(Drend)
    if(!is.null(DOY2)){
      today=DOY2[i]
    }else{
      today=NA
    }
    if(is.null(thetai)){
      thetai=NA
    }
    if(is.null(EF)){ 
      
      if(kctype!="single"){
        
        outputdual=rbind(outputdual,data.frame(Date=today,Day=i,ETo=ETo[[i]],Ineed=max(0,Ineed),TAW=TAW,
                                               RAW=RAW,Drstart=Drstart,Runoff=Runoff,
                                             Ks=Ks,Zr=Zr[i],h=h[i],
                                             Kc=kcb,Kcb=kcb,De=Dei,Drend=Drend,
                                             Kr=Kr,Ke=Ke,p=pi,
                                             DP=DP,CR=CR,GW=initgw,Evaporation=E,Transpiration=ETc2-E,Kc_adj=Kc,
                                             ETc_adj=ETc2,ETc=ETc,CN=CN,fw=fw,few=few,fc=fc,P=P[i],I=I[i],
                                             AW_modelled=max(AW,0),
                                             AW_measured=AW_measured,TSW_measured=TSW_measured,TSW_modelled=TSW,
                                             theta_measured=thetai,theta_modelled=mod_theta,crop=crop2,kctype=kctype))
                                             
      }else{
      
        outputdual=rbind(outputdual,data.frame(Date=today,Day=i,ETo=ETo[[i]],Zr=Zr[[i]],
                  TAW=TAW,RAW=RAW,Drstart=Drstart,
                                             Ineed=max(0,Ineed),Runoff=Runoff,p=pi,
                                             Ks=Ks,Kcb=kcb,h=h[[i]],Kc_adj=Kc,
                                             ETc_adj=ETc2,ETc=ETc,DP=DP,CR=CR,GW=initgw,
                                             Drend=Drend,CN=CN,fc=fc,P=P[i],I=I[i]
                                             ,AW_modelled=max(AW,0),AW_measured=AW_measured,
                                             TSW_measured=TSW_measured,TSW_modelled=TSW,
                                             theta_measured=thetai,theta_modelled=mod_theta,crop=crop2,kctype=kctype))
      #print(c(today=today,Day=i,ETo=ETo[[i]],Zr[i],TAW=TAW,RAW=RAW,Drstart=Drstart, Ks=Ks,fc=fc,lengths=lengths))
      }
    }
    
    #print(ETc2)
    #print(ETa1)
    #print(c(RAW,DP,Drstart,Drend))
    #print(cbind(Date=today,Day=i,ETo=ETo[[i]],Zr=Zr[[i]],TAW=TAW,Drstart=Drstart))
    i=i+1
  }
  #print(output)
  if(is.null(Ky)){ 
    Ky=1
  }
  if(is.null(EF)){
   
    if(length(Ky)==1){
      Ky=c(Ky,Ky,Ky,Ky,Ky)
    }
    if(length(Ky)>1){
      initoutput<- outputdual[ which(outputdual$Day <=stages2cum[1]),]
      ETdecrease_init=(1-(sum(initoutput$ETc_adj)/sum(initoutput$ETc)))
      Ya_by_Ym_init=(1-(Ky[1]*ETdecrease_init))
      yield_decrease_init=1-Ya_by_Ym_init
      cropoutput<- outputdual[ which(outputdual$Day>stages2cum[1]&outputdual$Day <=stages2cum[2]),]
      ETdecrease_crop=(1-(sum(cropoutput$ETc_adj)/sum(cropoutput$ETc)))
      Ya_by_Ym_crop=(1-(Ky[2]*ETdecrease_crop))
      yield_decrease_crop=1-Ya_by_Ym_crop
      midoutput<- outputdual[ which(outputdual$Day>stages2cum[2]&outputdual$Day <=stages2cum[3]),]
      ETdecrease_mid=(1-(sum(midoutput$ETc_adj)/sum(midoutput$ETc)))
      Ya_by_Ym_mid=(1-(Ky[3]*ETdecrease_mid))
      yield_decrease_mid=1-Ya_by_Ym_mid
      lateoutput<- outputdual[ which(outputdual$Day>stages2cum[3]&outputdual$Day <=stages2cum[4]),]
      ETdecrease_late=(1-(sum(lateoutput$ETc_adj)/sum(lateoutput$ETc)))
      Ya_by_Ym_late=(1-(Ky[4]*ETdecrease_late))
      yield_decrease_late=1-Ya_by_Ym_late
     # print(lateoutput)
      sumvariables=c("ETo","TAW","RAW","Drstart","Ineed","Runoff", "ETc_adj",
                     "ETc","DP","CR","GW","Drend","CN",          
                     "P" ,"I","AW_modelled","AW_measured","TSW_measured" ,"TSW_modelled")
      meanvariables=c("p","Ks","Kcb","Kc_adj","h","Zr","fc","theta_measured", "theta_modelled")
      initsum=(colSums(initoutput[sumvariables]))
      initmean=(colMeans(initoutput[meanvariables]))
      thisinit=c(initsum,initmean)
      cropsum=(colSums(cropoutput[sumvariables]))
      cropmean=(colMeans(cropoutput[meanvariables]))
      thiscrop=c(cropsum,cropmean)
      midsum=(colSums(midoutput[sumvariables]))
      midmean=(colMeans(midoutput[meanvariables]))
      thismid=c(midsum,midmean)
      latesum=(colSums(lateoutput[sumvariables]))
      latemean=(colMeans(lateoutput[meanvariables]))
      thislate=c(latesum,latemean)
      totalsum=(colSums(outputdual[sumvariables]))
      totalmean=(colMeans(outputdual[meanvariables]))
      thistotal=c(totalsum,totalmean)
    #names(thislate)=paste("late",names(thislate),sep="_")
    #print(names(thislate))
    }
    water.balance=data.frame(total=thistotal, initial=thisinit,crop=thiscrop,
                          mid=thismid,late=thislate)
    Ky=Ky[5]
    ETdecrease=(1-(sum(outputdual$ETc_adj)/sum(outputdual$ETc)))
    Ya_by_Ym=(1-(Ky*ETdecrease))
    yield_decrease=1-Ya_by_Ym
    yield_decrease=data.frame(total=yield_decrease, initial=yield_decrease_init,crop=yield_decrease_crop,
                              mid=yield_decrease_mid,late=yield_decrease_late)
    ETdecrease=data.frame(total=ETdecrease, initial=ETdecrease_init,crop=ETdecrease_crop,
                          mid=ETdecrease_mid,late=ETdecrease_late)
    
    Ya_by_Ym=data.frame(Ya_by_Ym.total=Ya_by_Ym, Ya_by_Ym.initial=Ya_by_Ym_init,Ya_by_Ym.crop=Ya_by_Ym_crop,
                        Ya_by_Ym.mid=Ya_by_Ym_mid,Ya_by_Ym.late=Ya_by_Ym_late)
    
    Ya_by_Yms=cbind(ETdecrease=ETdecrease,yield_decrease=yield_decrease,Ya_by_Ym)
    #Kc=c(Kcini=kcss[1],Kcmid=kcss[2],Kcend=kcss[2])

    factor=list(output=outputdual,yield_by_ET_decrease=Ya_by_Yms,
         parameters=parameters,water.balance.stages=water.balance,kctype=kctype)
    
  }else{
    
      names(AW2) =paste("Day",report3,sep="")
      names(Kc1)=paste("Day",report3,sep="")
      names(ETa1)=paste("Day",report3,sep="")
      names(Drend2)=paste("Day",report3,sep="")
     names(Ineed2)=paste("Day",report3,sep="")
     names(DP2)=paste("Day",report3,sep="")
     names(Runoff2)=paste("Day",report3,sep="")
     names(gw)=paste("Day",report3,sep="")
     names(CR2)=paste("Day",report3,sep="")
     if(!is.null(theta)){
       #print("yes")
     names(AW_measured2)=paste("Day",report3,sep="")
     #print("yes")
     
     }
     
     TAW=AW2+Drend2
     
     RAW=Drend2-Ineed2
     if(!is.null(xydata)){
     TAWxy=Awcxy+Drendxy
     TAWxy$longitude=Awcxy$longitude
     TAWxy$latitude=Awcxy$latitude
     RAWxy=Drendxy-Ineedxy
     RAWxy$longitude=Ineedxy$longitude
     RAWxy$latitude=Ineedxy$latitude
     
     }else{
       Drendxy=NULL
       TAWxy=NULL
       RAWxy=NULL
       Kcxy=NULL
       ETcxy=NULL
       Ineedxy=NULL
       AWxy=NULL
       Runoff2xy=NULL
       DP2xy=NULL
       CR2xy=NULL
       gwxy=NULL
       AW_measuredxy=NULL
     }
     if(!is.null(runofftest)&&!is.null(xydata)){
       names(Runoffxy)=paste("Day",day,sep="")
       Runoffxy$longitude=xydata$longitude
       Runoffxy$latitude=xydata$latitude
     }
      #names(TAW2)=paste("Day",report3,sep="")
      #names(RAW2)=paste("Day",report3,sep="")
    #}
      #TAW=TAW2,TAWxy=TAWxy,RAW=RAW2,RAWxy=RAWxy,Ineed=Ineed2,Ineedxy=Ineedxy,
      #print("okkkkkkkkkkkk2")
     ETdecrease=(1-((ETc_adjsum)/(ETc_sum)))
     Ya_by_Ym=(1-(Ky*ETdecrease))
     
     ETdecreasexy=NULL
     yield_decreasexy=NULL
     if(!is.null(xydata)){
     ETdecreasexy=xyextract(ETdecrease,xydata[,1],xydata[,2])
     colnames(ETdecreasexy)=c("longitude", "latitude","ETdecrease")
     yield_decreasexy=xyextract(1-Ya_by_Ym,xydata[,1],xydata[,2])
     colnames(yield_decreasexy)=c("longitude", "latitude","yield_decrease")
     }
     factor= list(EF=Kc1,EFxy=Kcxy,ETc_adj=ETa1,ETc_adj_xy=ETcxy,Drend=Drend2,
         Drendxy=Drendxy,TAW=TAW,TAWxy=TAWxy,RAW=RAW,RAWxy=RAWxy,
         AW_modelledxy=Awcxy,AW_modelled=AW2,AW_measured=AW_measured2,AW_measuredxy=AW_measuredxy,
         Ineed=Ineed2,Ineedxy=Ineedxy,Runoff=Runoff2,Runoffxy=Runoff2xy,
         DP=DP2,DPxy=DP2xy,CR=CR2,CRxy=CR2xy,gw=gw, gwxy=gwxy,xydata=xydata,map=map,
         ETdecrease=ETdecrease,ETdecreasexy=ETdecreasexy,
         yield_decrease=1-Ya_by_Ym,yield_decreasexy=yield_decreasexy)
  }
}
  factor$call<-match.call()
  
  class(factor)<-"kc"
  factor
}

#' @export
#' @rdname kc
plot.kc<-function(x,output="AW_measured",main=NULL,cex=1,legend=TRUE,pos="top",horiz=FALSE,...)
{
  #print(x$output$)
  if(is.null(x$EFxy)){
    #par(xpd = T, mar = par()$mar + c(0,0,0,7))
    refine="black"
    dots="p"
    label1="Measured"
    label2="modelled"
    
thismeasure=output
thismodel="AW_modelled"
thisylab="Available Soil Moisture[mm]"
  if(output=="TSW_measured"||output=="TSW_modelled"||output=="total"){
    thismeasure="TSW_measured"
    thismodel="TSW_modelled"
    thisylab="Total Available Soil Moisture[mm]"
    
  }
if(output=="theta"||output=="theta_modelled"||output=="theta_measured"){
  thismeasure="theta_measured"
  thismodel="theta_modelled"
  thisylab=" Soil Moisture Content[-]"
}

if(output=="ETc"||output=="ETc_adj"){
  thisylab="Evaporation and Transpiration[mm]"
  
  thismodel="ETc"
  thismeasure="ETc_adj"
  label1="ETc"
  label2="ETc_adj"  
}
if(output=="Evaporation"||output=="Transpiration"||output=="T"||output=="E"){
   thismodel="Evaporation"
   thismeasure="Transpiration"
   
  thisylab="Evaporation and Transpiration[mm]"
  refine="black"
  label1="E"
  label2="T"
  if(x$kctype=="single"){
    thismodel="ETc"
    thismeasure="ETc_adj"
    label1="ETc"
    label2="ETc_adj"
  }
 }
if(output=="Kcb"||output=="kcb"||output=="Kc"||output=="KC"){
  thismodel="Kc"
  thismeasure="Kc_adj"
  thisylab=" Crop Coefficient[-]"
  refine="black"
  label1="Kc_adj"
  label2="Kcb"
  if(x$kctype=="single"){
    thismodel="Kcb"
    label2="Kc"
  }
  }
all=NULL
if(output=="balance"){
  all=TRUE
  thismodel="ETc_adj"
  thismeasure="DP"
  
  thismeasure="Transpiration"
  dots="l"
  thisylab="P, DP [mm]"
  
  barplot(x$output$P+x$output$I,main=main,...)
  par(new = TRUE)
  #lines(x$output$Drend,type="l")
  lines(x$output$ETc_adj,col="green",lty=2,lwd=2)
  lines(x$output$Ineed, col="red",lty=3,lwd=2)
  #lines(x$output$Runoff, col="red",lty=4,lwd=2)
  if(length(pos)>1){
    legend(pos[1], pos[2],(c("P+I","ETc_adj","Deficit")), cex=cex, col=c("gray","green","red"),lty=1:3, lwd=2, bty="n",horiz=horiz)
  }
  else{
  legend(pos, (c("P+I","ETc_adj","Deficit")), cex=cex, col=c("gray","green","red"),lty=1:3, lwd=2, bty="n",horiz=horiz)
  }
  
  thisdate=as.character(x$output$Date)
  thisdate<-strftime(as.Date(thisdate),"%d-%m-%Y")
  labelss=as.character(thisdate)
  thisday=max(x$output$Day)
  
  thislist=round(c(1,thisday*0.1,thisday*0.2,thisday*0.30,
                   thisday*0.40,thisday*0.50,
                   thisday*0.6,thisday*0.7,
                   thisday*0.8,thisday*0.90,1*thisday))
  thislabel=labelss[thislist]
  if(is.na(x$output$Date[1])){
    thislabel=x$output$Day[[thislist]]
  }
  text(thislist, par("usr")[3], labels = thislabel, srt = 45, pos = 1,offset=2, xpd = TRUE,cex=cex)
  
  
 
}else{
  #if(type=="all"||type=="available"||type=="soil water"||type=="TSW_measured"
    # ||type=="theta"||type=="AW_measured"||type=="AW_modelled"
     #||type=="TSW_modelled"||type=="theta_modelled"){
  day=x$output$Day[which(!is.na(x$output[thismeasure]))]
  measured=x$output[[thismeasure]][which(!is.na(x$output[[thismeasure]]))]
  date=x$output$Date[which(!is.na(x$output[[thismeasure]]))]
  if(is.na(x$output$Date[1])){
    plot(x$output$Day,x$output[[thismodel]],type="l",ylab=thisylab,main=main,cex=cex,...)
  }else{
    plot(x$output$Day,x$output[[thismodel]],type="l",xaxt = "n", xlab="",ylab=thisylab,main=main,cex=cex,...)
    thisdate=as.character(x$output$Date)
    
    thisdate<-strftime(as.Date(thisdate),"%d-%m-%Y")
    labelss=as.character(thisdate)
    thisday=max(x$output$Day)
    thislist=round(c(1,thisday*0.1,thisday*0.15,thisday*0.2,thisday*0.25,thisday*0.30,
                     thisday*0.35,thisday*0.40,thisday*0.45,thisday*0.50,
                     thisday*0.55,thisday*0.6,thisday*0.65,thisday*0.7,thisday*0.75,
                     thisday*0.8,thisday*0.85,thisday*0.90,thisday*0.95,1*thisday))
    thislist=round(c(1,thisday*0.1,thisday*0.2,thisday*0.30,
                     thisday*0.40,thisday*0.50,
                     thisday*0.6,thisday*0.7,
                     thisday*0.8,thisday*0.90,1*thisday))
    #at=c(1,50,100,120)
    axis(1, at=thislist, labels=FALSE)
    #axis(1, at=seq(1, 10, by=1), labels = FALSE)
    text(thislist, par("usr")[3], labels = labelss[thislist], srt = 45, pos = 1,offset=2, xpd = TRUE,cex=cex)
    lines(day,measured ,xlab="Date",type=dots,col=refine,lty=1)
    if(length(pos)>1){
      legend(pos[1],pos[2],(c(label1,label2)), cex=cex,lty=0:1, bty="n",horiz=horiz,pch=c("o",""))
      
    }else{
      legend(pos,(c(label1,label2)), cex=cex,lty=0:1, bty="n",horiz=horiz,pch=c("o",""))
    }
    
  }
  }
  
}
}
#' @export
#' @rdname kc
fit.kc<-function(x,output="AW_measured"){
  thismeasure=output
  thismodel="AW_modelled"
  thisylab="Available Soil Moisture[mm]"
  mod=NULL
  if(is.null(x$EFxy)){
  obs=x$output$AW_measured
  est=(x$output$AW_modelled)
  if(output=="TSW_measured"||output=="TSW_modelled"||output=="total"){
    obs=x$output$TSW_measured
    est=x$output$TSW_measured
  }
  if(output=="theta"||output=="theta_modelled"||output=="theta_measured"){
    obs=x$output$theta_measured
    est=x$output$theta_modelled
  }
  mod=cbind(crop=as.character(x$output$crop[1]),fitting(obs,est))
  }else{
    #print("YESSSS")
    obs= data.frame(t((x$AW_measuredxy)))
    est= data.frame(t((x$AW_modelledxy)))
    mod=NULL
    i=1
    while(i<=length(est)){
      
      mod=rbind(mod,cbind(longitude=est[1,][[i]],latidute=est[2,][[i]],fitting(obs[3:nrow(obs),][[i]],est[3:nrow(est),][[i]])))
      
      i=i+1
      
    }

  }
  mod
}
