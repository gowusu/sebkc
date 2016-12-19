#' Hourly FAO-56 and ASCE-EWRI Reference Evapotranspiration
#' 
#' The function calculates hourly or less than hour FAO-56 (grass) and 
#' ASCE-EWRI Reference (alfalfa) Evapotranspiration (ETo ETr). 
#' The "surface" parameter determines #' the type of ET. 
#' The function can also be used to scale hourly or less than hour actual evapotranspiration.
#' The daily ET is implemented in \code{\link{ETo}}.
#' @param Tmean Numeric. Hourly Mean Air Temperature in degree Celsius
#' @param RH Hourly Mean Relative Humidity in percent
#' @param Lz longitude of the centre of the local time zone [degrees west of Greenwich].
#' MUST BE POSITIVE
#' Lz = 0° for Greenwich.
#' @param Lm longitude of the measurement site [degrees west of Greenwich]. MUST BE POSITIVE
#' @param t1 The length of the calculation period in hour; 1 for hour, 0.5 for 30 minutes
#' 0.25 for 15 minutes. 
#' @param time The midpoint of the time of measurement[hour]; for example
#'  time is 12.5 for a period between 12:00 and 13:00
#' @param ETins model. It takes ET model from \code{\link{sebal}} or \code{\link{tseb}} and
#' retrieve instantaneous ET estimates. This helps to estimate ET24, T24, E24, EF, Kc, EFs, EFc
#' @param Rn net radiation [MJ/m2/hr]
#' @param Rs solar radiation [MJ/m2/hr]
#' @param map spatial raster data for inverse distance interpolation
# @param scaling numeric. The type of scaling models to upscale instantanous ET to 
# daily ET. it takes 1 is Evaporative Fraction (Bastiaanssen et al 1998); 
# 2 for Sine Function  (Jackson et al 1983); 3 for Solar Radiation (French et al 2005),
# 4 for time scaling of alfalfa (Allen et al 2007); 5 for time scaling of grass ()
#' @inheritParams sebal
#' @inheritParams ETo 
#' @inheritParams tseb
#' @inheritParams weather
#' @seealso \code{\link{ETo}}
#' @details 
#' \code{\link{sebal}}, \code{\link{sebi}},\code{\link{sebs}},\code{\link{sseb}},
#' \code{\link{tseb}} 
#' @author George Owusu
#' @return 
#' \itemize{
#' \item{ETo:} {hour reference ETo [mm/hr]}
#' \item{ETa:} {hour actual ET. It is not NULL if EF has a numeric value}
#' \item{ETc:} {hour crop ET. It is not NULL if Kc has a numeric value}
#' \item{Rn:} { Net  Radiation  [W/m2] }
#' \item{Rnwm_2:} { Net Radiation [MJ/m2/hr] }
#' \item{y:} { psychrometric constant [kPa/degree Celsius]}
#' \item{slope:} { slope of saturation vapour pressure 
#' curve [kPa/degree Celsius]}
#' \item{vpd:} { Vapour pressure deficit[kPa]}
#' }
#' @examples
#' \dontrun{
#' #FAO56 Example 19 14:00-15:00 h using Cd = 0.24 instead 0.34 in FAO56
#' ET1hr=ETohr(Tmean=38,RH=52,DOY=274,Lz=15,t1=1,time=14.5,Lm=16.25,latitude=16.22,uz=3.3,Rs=2.450)
#' ET1hr$ETo
#' }
#' @references
#' \itemize{ 
#' \item{}{ALLEN, R. G., PEREIRA, L. S., RAES, D., & SMITH, M. 1998. 
#' Crop Evapotranspiration (guidelines for computing crop water requirements) 
#' FAO Irrigation and Drainage Paper No. 56: FAO.}
#' 
#' \item{}{Jackson, R. D., Hatfield, J. L., Reginato, R. J., Idso, S. B., & Pinter Jr, P. J. (1983). 
#' Estimation of daily evapotranspiration from one time-of-day measurements. 
#' Agricultural Water Management, 7(1-3), 351-362. doi: http://dx.doi.org/10.1016/0378-3774(83)90095-1}
#' 
#' \item{}{French, A. N., Fitzgerald, G., Hunsaker, D., Barnes, E., Clarke, T., Lesch, S., . . .
#'  Pinter, P. (2005). 
#' Estimating spatially distributed cotton water use from thermal infrared aerial imagery. 
#' Paper presented at the World Water Congress 2005: 
#' Impacts of Global Climate Change - Proceedings of the 2005 World Water and 
#' Environmental Resources Congress, Reston, Va.}
#' 
#' \item{}{Colaizzi, P. D., Kustas, W. P., Anderson, M. C., Agam, N., Tolk, J. A., 
#' Evett, S. R., . . . O'Shaughnessy, S. A. (2012). Two-source energy balance model
#'  estimates of evapotranspiration using component and composite surface temperatures. 
#'  Advances in Water Resources, 50, 134-151. 
#'  doi: http://dx.doi.org/10.1016/j.advwatres.2012.06.004}
#'  }
#' @export
#' @rdname ETohr 
ETohr=function(Tmean,RH,DOY,Lz,t1,time,Lm,latitude,longitude=NULL,n=NULL,uz=2,Rs=NULL,
Rn=NULL,G=NULL,altitude=NULL,z=NULL, Krs=0.16,Tmax=NULL,Tmin=NULL,albedo=0.23,as=0.25,bs=0.5,
ETins=NULL,ETr24=NULL,wmo=NULL, airport=NULL,map=NULL,surface="grass",period="daytime")UseMethod ("ETohr")
#' @importFrom gstat gstat 
#' @import sp 
#' @export
#' @rdname ETohr 
ETohr.default=function(Tmean,RH,DOY,Lz,t1=1,time,Lm,latitude,longitude=NULL,n=NULL,uz=2,Rs=NULL,
  Rn=NULL,G=NULL,altitude=NULL,z=NULL, Krs=0.16,Tmax=NULL,Tmin=NULL,albedo=0.23,as=0.25,bs=0.5,
  ETins=NULL,ETr24=NULL,wmo=NULL, airport=NULL,map=NULL,surface="grass",period="daytime"){
  filedata=NULL
  timestep="hour"
  print("Computing ETr")
    if(!is.null(wmo)||!is.null(airport)){
   Tmean=sebkc.tryCatch(Tmean)$value
   
   if(class(Tmean)[1]=="simpleError"){
    Tmean=NULL 
   }
   if(is.null(Tmean$time)){
    time=time 
   }else{
     time=NULL
   }
   Tmean=sebkc.tryCatch(weather(data=Tmean,airport = airport,wmo=wmo, date=DOY,time=time))$value
   if(class(Tmean)[1]=="simpleError"){
    return (print(paste("Invalid DOY [YYYY-mm-dd] or airport or wmo",Tmean))) 
   }
   Tmean= Tmean$WMO$hour
   
     
  }
  
 
  
if(is.character(Tmean)||!is.null(names(Tmean))){
  
  if(!is.null(names(Tmean))){
    if(!is.null(Tmean$WMO$hour)){
      Tmean=Tmean$WMO$hour
    }
    filedata=Tmean 
  }else{
    filedata= read.csv(Tmean)
  }
   #print(names(filedata))
  if(!is.null(filedata$Tmean)){
    Tmean=filedata$Tmean
  }
   if(!is.null(filedata$RH)){
     RH=filedata$RH
   }
   if(!is.null(filedata$DOY)){
     DOY=filedata$DOY
   }
   
   if(!is.null(filedata$Lz)){
     Lz=filedata$Lz
   }
  if(!is.null(filedata$t1)){
    t1=filedata$t1
  } 
  if(!is.null(filedata$time)){
    time=filedata$time
  }  
   if(!is.null(filedata$uz)){
     uz=filedata$uz
   } 
   if(!is.null(filedata$latitude)){
     latitude=filedata$latitude
   } 
   if(!is.null(filedata$Lm)){
     Lm=filedata$Lm
   } 
   if(!is.null(filedata$z)){
     z=filedata$z
   } 
   if(!is.null(filedata$altitude)){
     altitude=filedata$altitude
   } 
   if(!is.null(filedata$Rs)){
     Rs=filedata$Rs
   }
   
   if(!is.null(filedata$Rn)){
     Rn=filedata$Rn
   }
   
   if(!is.null(filedata$n)){
     n=filedata$n
   }
   if(!is.null(filedata$Krs)){
     Krs=filedata$Krs
   }
  if(!is.null(filedata$Tmax)){
    Tmax=filedata$Tmax
  }
  if(!is.null(filedata$Tmin)){
    Tmin=filedata$Tmin
  }
   if(!is.null(filedata$albedo)){
     albedo=filedata$albedo
   }
   if(!is.null(filedata$as)){
     as=filedata$as
   }
   if(!is.null(filedata$bs)){
     bs=filedata$bs
   }
   if(!is.null(filedata$EF)){
     EF=filedata$EF
   }
   if(!is.null(filedata$Kc)){
     Kc=filedata$Kc
   }
  
  if(!is.null(filedata$timestep)){
    timestep=filedata$timestep
  }
  
  if(!is.null(filedata$surface)){
    surface=filedata$surface
  }
  
  if(!is.null(filedata$period)){
    period=filedata$period
  } 
  
  }  

  
#print(cbind(Tamax,Tamin))
#check hr or time
Cn=37
Cd=0.34
if(surface=="grass"){
  if(timestep=="hourly"||timestep=="hours"||timestep=="Hourly"||timestep=="hour"){
    Cn=37
    if(period=="nighttime"||period=="night"){
      Cd=0.96
    }else{
      Cd=0.24
    }
  }
}

if(surface=="ETr"||surface=="alfalfa"||surface=="Alfalfa"){
  if(timestep=="hourly"||timestep=="hours"||timestep=="Hourly"||timestep=="hour"){
    Cn=66
    
    if(period=="nighttime"||period=="night"){
      Cd=1.7
    }else{
      Cd=0.25
    }
  }else{
    Cn=1600
    Cd=0.38
  }
  
 
  
}

DOY2=DOY
if(!is.numeric(DOY)){
DOY=strptime(DOY,"%d-%m-%Y")$yday+1
if(is.na(DOY[1])){
  DOY=strptime(DOY2,"%Y-%m-%d")$yday+1
}
if(is.na(DOY[1])){
  DOY=strptime(DOY2,"%m/%d/%Y")$yday+1
}
if(is.na(DOY[1])){
  DOY=strptime(DOY2,"%m-%d-%Y")$yday+1
}
}

#print(DOY)
if(z!=2&&!is.null(z)){
u2=uz*(4.87/log((67.8*z)-5.42))
}else{
u2=uz
}
if(is.null(altitude)){
P=101.3
}else{
P=101.3*((293-(0.0065*altitude))/293)^5.26
}


slope=(4098*(0.6108*exp((17.27*Tmean)/(Tmean+237.3))))/(Tmean+237.3)^2
cp=1.013*10^-3
ratio=0.622
vaporization=2.45
y=(cp*P)/(ratio*vaporization)
es=0.6108*exp((17.27*Tmean)/(Tmean+237.3))
ea=es*(RH/100)
vpd=es-ea

#radiation
if(is.null(Rn)||is.null(Rs)){
lat=(pi/180)*latitude
Gsc=0.0820
dr=1+0.033*cos(((2*pi)/365)*DOY)
declination=0.409*sin((((2*pi*DOY)/(365)))-1.39)
ws=acos(-tan(lat)*tan(declination))
b=(2*pi*(DOY-81))/364
Sc=0.1645*sin(2*b)-0.1255*cos(b)-0.025*sin(b)
#w=(pi/12)*((time+(0.06667*(Lz-Lm)+Sc))-12)
w=(pi/12)*((time+0.06667*(Lz-Lm)+Sc)-12)
ws1=w-((pi*t1)/24)
ws2=w+((pi*t1)/24)
Ra=((12*60)/pi)*Gsc*dr*((ws2-ws1)*sin(lat)*sin(declination)+cos(lat)*cos(declination)*(sin(ws2)-sin(ws1)))
N=(24/pi)*ws
if(is.null(Rs)){
if(is.null(n)){
if(is.null(Tmax)||is.null(Tmin)||is.null(Krs)){
 print("daily Tmax, Tmin and Krs needed, please") 
}  
  Rs=Krs*sqrt(Tmax-Tmin)*Ra  
}else{  
  
  nN=n/N
#Rs=(as+(bs*nN))*Ra
#print(Rs)
#rc2=0.4560+(0.3566*nN)+(0.1874*((nN)^2))
rc=0.3139+(0.2480*nN)+(0.1276*((nN)^2))
#rc=0.1715+(0.8419*nN)+(-0.3206*((nN)^2))

if(n==0){
 rc=0.1811 
}
Rs=rc*Ra
}
#print(Rs)
#print(rc2*Ra)

}
#Rso=(0.75+(2*10^-5)*altitude)*Ra
Rso=(as+bs)*Ra
RsRso=Rs/Rso
Rns=(1-albedo)*Rs
stefan=2.043*10^-10
Rnl=(stefan*((Tmean+273.16)^4))*(0.34-(0.14*sqrt(ea)))*((1.35*RsRso)-0.35)
if(is.null(Rn)){
  Rn=Rns-Rnl
}
}
if(is.null(G)){
G=0.1*Rn
if(period=="night"||period=="nighttime"){
G=0.5*Rn  
}
}

wind=(1+(Cd*u2))
sloping=slope/(slope+(y*wind))
y2=y/(slope+y*(wind))
ES=0.408*(Rn-G)*sloping
yvpd=Cn/(Tmean+273.16)*u2*(vpd)*y2
ETo=ES+yvpd
ETc=NULL
#print(list(ES=ES,y2=y2,slope=slope,yvpd=yvpd,sloping=sloping,Cd=Cd,Cn=Cn))

ETa=NULL
E24=T24=EFc=EFs=EF=ET24=NULL
if(!is.null(ETins)){
  ETins.class=tolower(class(ETins))
  vaporisation=2430000
  if(ETins.class=="sebal"||ETins.class=="tseb"
     ||ETins.class=="metric"){
    if(!is.null(ETins$vaporisation)){
      #print(ETins)
      vaporisation=ETins$vaporisation
    } 
    
    ETint=(t1*3600*ETins$LE)/vaporisation
    
    EF=ETint/ETo

    if(ETins.class=="tseb"){
      Tint=(t1*3600*ETins$LEc)/vaporisation
      Sint=(t1*3600*ETins$LEs)/vaporisation
      #print(Tint)
      EFc=Tint/ETo
      EFs=Sint/ETo
      if(!is.null(ETr24)){
        E24=EFs*ETr24 
        T24=EFc*ETr24
      }
      }
  
  }else{
   EF=ETins/ETo 
  }
  if(!is.null(ETr24)){
    ET24=EF*ETr24 
  }
}
#print(paste("Cn",Cn,"Cd",Cd,"ETo",ETo))

if(!is.null(longitude)||!is.null(map)){
  if(length(longitude)>1){
    map=invdist(longitude,latitude,Rn,ext=map)
    Rn <- map$var
    map=invdist(longitude,latitude,ETo,ext=map)
    ETo <- map$var
  }
}
if(!is.null(filedata)){
  filedata$ETo=ETo  
  filedata$Rn=Rn  
  filedata$slope=slope  
  filedata$vpd=vpd 
  filedata$RH=RH  
  filedata$u2=u2
  filedata$Tmean=Tmean
  
}
factor<-list(ETo=ETo,Rn=11.6*Rn,Rs=Rs,G=G,Rn_mj_m2_hr=Rn,y=y,slope=slope,vpd=vpd,ETa=ET24,E24=E24,T24=T24,EFc=EFc,EFs=EFs,
             EF=EF,Kc=EF,xy=data,ETc=ETc,u2=u2,DOY=DOY2,
             latitude=latitude,Lm=Lm,data=filedata)
factor$call<-match.call()
class(factor)<-"ETohr"
factor
}

#file="C:/Users/GeoKings/Documents/George Owusu/UG PhD/literature/sebal/NASAweather1983-20052.csv"
#data=read.csv("C:/Users/GeoKings/Documents/George Owusu/UG PhD/literature/sebal/NASAweather1983-20052.csv")

#mod=ETo(file)
#mod=ETo(data)