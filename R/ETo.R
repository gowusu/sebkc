#' 24 Hour FAO 56 and ASCE-EWRI Reference Evapotranspiration ETo
#'
#' @param Tmax Numeric. Maximum air Temperature in degree Celsius
#' @param Tmin Numeric. Minimum air Temperature in degree Celsius
#' @param latitude geographical coordinates in decimal degrees. It should be negative
#' for southern hemisphere
#' @param longitude The longitude of the measurement site 
#' i.e. geographical coordinates in decimal degrees for the weather station.
#'  It should be negative for West and positive for East.
#' @param Rs shortwave incoming radiation [MJ/m2/day]
#' @param z The height above surface where wind speed was measured, in metres. 
#' @param uz The wind speed in m/s
#' @param altitude The elevation of the weather station, in metres
#' @param RHmax Maximum Relative Humidity in percent
#' @param RHmin Minimum Relative Humidity in percent
#' @param n Sunshine hours 
#' @param Krs Relationship between the fraction of extraterrestrial radiation that reaches the 
#' earth's surface, Rs/Ra, and the air temperature difference Tmax - Tmin for
#' interior (Krs = 0.16) and coastal (Krs = 0.19) regions
#' @param albedo albedo values. For grass it is 0.23[dimensionless] 
#' @param as Regression constant, intercept, expressing the fraction of extra-terrestrial 
#' radiation reaching the earth on overcast days (n = 0); 
#' either a calibrated value or as = 0.25 can be used.
#' @param bs regression constant, slope, expressing 
#' the fraction of extraterrestrial radiation
#' reaching the earth on overcast days (n = 0); either a calibrated value or bs = 0.5 can be used.
#' @param Rn Net daily radiation in MJ/m2/day. If it is NULL, 
#' it will be computed
#' @param scale character The type of scale models to upscale instantaneous ET to 
#' daily ET. it takes "Evaporative Fraction" (Bastiaanssen et al 1998); 
#' "Sine Function"  (Jackson et al 1983); "Solar Radiation" (French et al 2005). 
#' Remember the blank spaces between the parameters. 
#' You can also respectively use "evaporative" or "solar" or "sine".
# ,4 for time scale of alfalfa (Allen et al 2007); 5 for time scale of grass ()
#' @param surface character. It is either "grass" or "alfalfa"
#' @param Kc Crop coefficient
#' @param G Daily Soil heat flux in W/m2 
#' @param EF Spatial. ET fraction [-] from \link{sebal} or \link{tseb} \link{sebs}.
#'  It is advisable to link it to any of the two functions. 
#' @inheritParams sebal
#' @inheritParams weather
#' @inheritParams ETohr
#' @seealso 
#' \code{\link{ETohr}}
#' @author George Owusu
#' @return 
#' \itemize{
#' \item{ETo:} { 24 hour reference ETo [mm/day]}
#' \item{ETa:} { 24 hour actual ET. It is not NULL if EF has a numeric value}
#' \item{ETc:} { 24 hour crop ET. It is not NULL if Kc has a numeric value}
#' \item{Rn:} { Net daily Radiation  [W/m2] }
#' \item{Rn_mj_m2_day:} { Net daily Radiation [MJ/m2/day]}
#' \item{y:} { psychrometric constant [kPa/degree Celsius]}
#' \item{slope:} { slope of saturation vapour pressure 
#' curve [kPa/degree Celsius]}
#' \item{vpd:} { Vapour pressure deficit[kPa]}
#' }
#' @examples
#' \dontrun{
#' #with minimal data
#' PET=ETo(Tmax=31,Tmin=28,latitude=5.6,DOY="6-2-2001")
#' PET$ETo
#' #with enough data as in FAO56 EXAMPLE 18 on 6 July in Uccle (Brussels, Belgium) 
#' #located at 50 degrees 48 minutes N and at 100 m above sea level
#' PET24=ETo(Tmax=21.5,Tmin=12.3,z=10,uz=2.78,altitude=100,
#' RHmax=84,RHmin=63,n=9.25,DOY=187,latitude=50.8)
#' PET24$ETo
#' }
#' @references ALLEN, R. G., PEREIRA, L. S., RAES, D., & SMITH, M. 1996. 
#' Crop Evapotranspiration (guidelines for computing crop water requirements) 
#' FAO Irrigation and Drainage Paper No. 56: FAO.
#' @export
#' @rdname ETo 
ETo=function(Tmax,Tmin,DOY,latitude,longitude=NULL,z=NULL,uz=2,
  altitude=NULL,RHmax=NULL,RHmin=NULL,n=NULL,Krs=0.16,albedo=0.23,
      as=0.25,bs=0.5,Rs=NULL,Rn=NULL,G=0,EF=NULL,wmo=NULL, airport=NULL,scale="sine function",Kc=NULL,
  surface="grass") UseMethod ("ETo")
#' @importFrom gstat gstat 
#' @import sp 
#' @export
#' @rdname ETo 
ETo.default=function(Tmax,Tmin,DOY,latitude,longitude=NULL,z=NULL,uz=2,altitude=NULL,
RHmax=NULL,RHmin=NULL,n=NULL,Krs=0.16,albedo=0.23,
as=0.25,bs=0.5,Rs=NULL,Rn=NULL,G=0,EF=NULL,wmo=NULL, airport=NULL,scale="evaporative fraction",Kc=NULL,
surface="grass"){
  print("Computing ETr24 or ETo24")
  thislat=sebkc.tryCatch(latitude)$value
  if(class(thislat)[1]=="simpleError"){
    return(print("latitude is needed"))
    
  }
  filedata=NULL
  timestep="24-h"
  period="daytime"
  
  if(!is.null(wmo)||!is.null(airport)){
    Tmax=sebkc.tryCatch(Tmax)$value
    
    if(class(Tmax)[1]=="simpleError"){
      Tmax=NULL 
    }
   #print(Tmax)
    Tmax=sebkc.tryCatch(weather(data=Tmax,airport = airport,wmo=wmo, date=DOY))$value
    if(class(Tmax)[1]=="simpleError"){
      return (print(Tmax)) 
    }
    Tmax= Tmax$WMO$day
    
  }
  #altitude=welev
if(is.character(Tmax)||!is.null(names(Tmax))){
  
  if(!is.null(names(Tmax))){
    if(!is.null(Tmax$WMO$hour)){
      Tmax=Tmax$WMO$day
    }
    filedata=Tmax 
  }else{
    filedata= read.csv(Tmax)
  }
   #print(names(filedata))
  Tmax=filedata$Tmax
   if(!is.null(filedata$Tmin)){
     Tmin=filedata$Tmin
   }
  
   
   if(!is.null(filedata$DOY)){
     DOY=filedata$DOY
   }
   
   
   if(!is.null(filedata$uz)){
     uz=filedata$uz
   } 
   if(!is.null(filedata$latitude)){
     latitude=filedata$latitude
   } 
   if(!is.null(filedata$longitude)){
     longitude=filedata$longitude
   } 
   if(!is.null(filedata$z)){
     z=filedata$z
   } 
   if(!is.null(filedata$altitude)){
     altitude=filedata$altitude
   } 
   if(!is.null(filedata$RHmin)){
     RHmin=filedata$RHmin
   }
   if(!is.null(filedata$RHmax)){
     RHmax=filedata$RHmax
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

Tamax=Tmax
Tamin=Tmin
#print(cbind(Tamax,Tamin))
#check day or time
Cn=900
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

Tmean=(Tamax+Tamin)/2

slope=(4098*(0.6108*exp((17.27*Tmean)/(Tmean+237.3))))/(Tmean+237.3)^2
cp=1.013*10^-3
ratio=0.622
vaporization=2.45
y=(cp*P)/(ratio*vaporization)
etmax=0.6108*exp((17.27*Tamax)/(Tamax+237.3))
etmin=0.6108*exp((17.27*Tamin)/(Tamin+237.3))
es=(etmax+etmin)/2
if(is.null(RHmax)||is.null(RHmin)){
ea=etmin
}else{
ea=((etmin*(RHmax/100))+(etmax*(RHmin/100)))/2
}


vpd=es-ea

#radiation
if(is.null(Rn)){
lat=(pi/180)*latitude
Gsc=0.0820
dr=1+0.033*cos(((2*pi)/365)*DOY)
declination=0.409*sin((((2*pi*DOY)/(365)))-1.39)
ws=acos(-tan(lat)*tan(declination))
Ra=((24*60)/pi)*Gsc*dr*((ws*sin(lat)*sin(declination))+(cos(lat)*cos(declination)*sin(ws)))
N=(24/pi)*ws
nN=n/N
if(is.null(Rs)){
if(is.null(n)){
  if(is.null(Tmax)||is.null(Tmin)||is.null(Krs)){
    print("daily Tmax, Tmin and Krs needed, please") 
  } 
Rs=Krs*sqrt(Tamax-Tamin)*Ra
}else{
Rs=(as+(bs*nN))*Ra
}
}
if(as==0.25&&bs==0.5&&!is.null(altitude)){
  Rso=(0.75+(2*10^-5)*altitude)*Ra  
}else{
  Rso=(as+bs)*Ra
}
#
RsRso=Rs/Rso
Rns=(1-albedo)*Rs
tmaxk=Tamax+273.16
tmink=Tamin+273.16
stefan=4.903*10^-9
Rnl=(stefan*(((tmaxk^4)+(tmink^4))/2))*(0.34-(0.14*sqrt(ea)))*((1.35*RsRso)-0.35)
Rn=Rns-Rnl
}
wind=(1+(Cd*u2))
sloping=slope/(slope+(y*wind))
y2=y/(slope+y*(wind))
ES=0.408*(Rn-G)*sloping
yvpd=Cn/(Tmean+273.16)*u2*(vpd)*y2
ETo=ES+yvpd
ETc=NULL
if(!is.null(Kc)){
 ETc=ETo*Kc 
}
ETa=NULL
if(!is.null(EF)){
  model=NULL
  EF.class=tolower(class(EF))
  if(!is.numeric(scale)){
   scale= tolower(scale)
   if(scale=="evaporative"||scale=="evaporative fraction"||scale=="fraction"){
  scale=1
   }
   if(scale=="sine"||scale=="sine function"||scale=="function"){
     scale=2
   }
   if(scale=="solar"||scale=="radiation"||scale=="solar radiation"){
     scale=3
   }
     }
  if(EF.class=="sebal"||EF.class=="tseb"
     ||EF.class=="metric"){
    model=EF
    if(!is.null(EF$EF)){
      EF=EF$EF
    }else{
      EF=1.1*(EF$LE/(EF$Rn-EF$G))
      } 
    
  }
  
  if(!is.null(longitude)||!is.null(EF)){
    if(length(longitude)>1){
    map=invdist(longitude,latitude,Rn,ext=EF)
    Rn <- map$var
    map=invdist(longitude,latitude,ETo,ext=EF)
    ETo <- map$var
    }
  }
  #Evaporative Fraction
if(scale==1){
  ETa=(EF*Rn)/2.45
}
  #Sine Function
  if(!is.null(model)){
    #print(model$LE)
  if(scale==2){
    ETa=((model$LE/model$Rs)*Rs)/2.45
  }
  #Solar Radiation
  if(scale==3){
    ETa=(EF*Rs*(model$Rn/model$Rs))/2.45
      }
  }
}
#print(paste("Cn",Cn,"Cd",Cd,"ETo",ETo))
if(is.null(RHmin)){
  RHmin=(etmin/etmax)*100
}

if(!is.null(filedata)){
  filedata$ETo=ETo  
  filedata$Rn=11.6*Rn  
  filedata$slope=slope  
  filedata$vpd=vpd 
  filedata$RHmin=RHmin  
  filedata$u2=u2
  
}
factor<-list(ETo=ETo,Rn=11.6*Rn,Rs=Rs,Rn_mj_m2_day=Rn,y=y,slope=slope,vpd=vpd,ETa=ETa,
             EF=EF,xy=data,ETc=ETc,u2=u2,RHmin=RHmin,DOY=DOY2,
             latitude=latitude,longitude=longitude,data=filedata,N=N,nN=nN,Ra=Ra,Rns=Rns,Rnl=Rnl)
factor$call<-match.call()
class(factor)<-"ETo"
factor
}

#file="C:/Users/GeoKings/Documents/George Owusu/UG PhD/literature/sebal/NASAweather1983-20052.csv"
#data=read.csv("C:/Users/GeoKings/Documents/George Owusu/UG PhD/literature/sebal/NASAweather1983-20052.csv")

#mod=ETo(file)
#mod=ETo(data)