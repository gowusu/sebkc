#DOY="1-8-2001"
#latitude=0
#model="CBM"
#Tmax=31
#Tmin=26

#' Day length calculation
#' @inheritParams ETo
#'
#' @param model character Type of model:"CBM" (Schoolfield, 1982),"BGC" (Running & Coughlan, 1988),
#'  "CERES" (Ritchie, 1991) or "FAO56"
#' @param p numeric. CMB parameter
#'
#' @return DL
#' @author George Owusu
#' @export
#' 
#' @references 
#' \itemize{
#' \item{}{Schoolfield R. (1982). Expressing daylength as a function of latitude and Julian date.}
#' \item{}{Running S. W. and Coughlan J. C. (1988). A general model of forest ecosystem processes 
#' for regional applications I. Hydrologic balance, canopy gas exchange and primary production 
#' processes. Ecological Modelling, 42(2), 125-154 http://dx.doi.org/10.1016/0304-3800(88)90112-3.}
#' \item{}{Ritchie J. T. (1991). Wheat Phasic Development. Modeling Plant and Soil Systems. 
#'  Agronomy Monograph, 31.}
#' \item{}{Allen R. G., Pereira L. S., Raes D. and Smith M. (1998). Crop evapotranspiration: 
#' Guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper, 56, 300.}
#' }
#' 
#'
#' @examples
#' DOY="2001-8-1"
#' latitude=0
#' model="CBM"
#' Tmax=31
#' Tmin=26
#' mod=dl(latitude,DOY,model,Tmax,Tmin=Tmin)

dl=function(latitude,DOY,model="CBM",Tmax=NULL,Tmin=NULL,p=0.5){
  
  
if(!is.numeric(DOY)){
  DOY2=DOY
  DOY=strptime(DOY,"%Y-%m-%d")$yday+1
  if(is.na(DOY)){
    DOY=strptime(DOY2,"%d-%m-%Y")$yday+1
  }
}


#better
if(model=="CBM"||model=="cbm"){

phi=0.2163108+2*atan(0.9671396*tan(0.00860)*(DOY-186))
lat2=asin(0.39795*cos(phi))
DL=24-(24/pi)*(acos(((sin((p*pi)/180))+(sin((latitude*pi/180)*sin(lat2)))/(cos((latitude/180)*cos(phi))))))
}

if(model=="Brock"||model=="brock"||model=="BROCK"){
  phi=23.45*sin(360*((283+DOY)/(365)))
  hourangle=acos((-1*tan(latitude)*tan(phi)))
  DL=2*(hourangle/15)
  
}

if(model=="BGC"||model=="bgc"||model=="Bgc"){
  ampl=exp((7.42+0.045*latitude)/3600)
  DL=ampl*sin((DOY-79)*0.01721)+12
  
}

if(model=="CERES"||model=="ceres"||model=="Ceres"){
 phi=0.4093*sin(0.0172*(DOY-82.2))
#DL=7.639*acos(max(-0.87,((-sin((latitude*pi)/180)*sin(phi)-0.1047))/cos((latitude*pi)/180)*cos(phi)))
#DL=7.639*acos(max(-0.87,(-sin((latitude*pi)/180)*sin(phi)-0.1047)/(cos((latitude*pi)/180)*cos(phi))))
DL=7.639*acos(max(-0.87,((-sin((latitude*pi)/180)*sin(phi)-(0.1047))/(cos((latitude*pi)/180)*(cos(phi))))))
}
if(model=="FAO56"){
  declination=0.409*sin((((2*pi*DOY)/(365)))-1.39)
  ws=acos(-tan(latitude)*tan(declination))
  DL=(24/pi)*ws  
}  
  
  
  
  T24=NULL
  Tn=NULL
  Rs=NULL
  Tday=NULL
  if(!is.null(Tmin)&&!is.null(Tmax)){
Ta=0.5*(Tmax+Tmin)
h=12-(0.5*DL)
Tday=Ta+((Tmax-Tmin)/(4*pi))*((11+h)/(12-h))*sin(pi*((11-h)/(11+h)))
Tn=Ta+((Tmax-Tmin)/(4*pi))*((11+h)/(h))*sin(pi*((11-h)/(11+h)))
#T24=(Tn+Tday)/2

Ta2=0.5*(Tmax-Tmin)
T24=Ta+Ta2*((24-Tmax)/(pi*(DL-Tmax)))*sin((pi*(DL-Tmax))/(24-Tmax))
if(is.na(T24)){
 # print("yes")
  T24=(Tday+Tn)*0.5  
}
PET=ETo(Tmax=Tmax,Tmin=Tmin,latitude=latitude,DOY=DOY)
Rs= PET$Rs*23.9
  }
  

list(Tday=Tday, T24=T24, Tn=Tn,Rs=Rs,DL=DL)
}
