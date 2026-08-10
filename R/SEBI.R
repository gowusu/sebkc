
#' @title SEBI Evapotranspiration Fraction Estimation in R
#'
#' @description \code{sebi} computes Evapotranspiration Fraction Estimation
#' based on Menenti& Choudhury (1993) parameterization of 
#' land surface evaporation by means of location dependent potential 
#' evaporation and surface temperature range. 
#' @inheritParams sebal
#' @inheritParams ETo
#' 
#' @details Function \code{\link{ETo}} will be automatically used to 
#' estimate  parameters y (psychrometric constant), 
#' slope (slope of saturation vapour pressure curve) and 
#' vpd (Vapour pressure deficit).
#' @references 
#' Menenti, M., & Choudhury, B. 1993. Parameterization of land surface 
#' evaporation by means of location dependent potential evaporation and 
#' surface temperature range. Paper presented at the Proceedings of 
#' IAHS conference on Land Surface Processes.
#' @return EF Evapotranspiration Fraction
#' @author George Owusu
#'
#' @examples
#' \dontrun{
#' #using folder parameter 
#' folder=system.file("extdata","stack",package="sebkc")
#' sebiauto=sebi(folder=folder,welev=317,Tmax=31,Tmin=28)
#' 
#' #manual input of parameters
#' #data
#' albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="sebkc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
#' LAI=raster(system.file("extdata","LAI.grd",package="sebkc"))
#' #model
#' modsebi=sebi(albedo=albedo,Ts=Ts,Tmax=31,Tmin=28,RHmax=84,
#' RHmin=63,NDVI=NDVI,SAVI=NULL,iter.max=7,xyhot="full",
#' xycold="full",DOY=37,sunelev=50.71154048,welev=317.1,zx=10,
#' u=2,zomw=2,zom=NULL,LAI=LAI,Rn24=NULL,model="sebi")
#' 
#' 
#' #Another interactive example
#' #get Ts, Albedo, NDVI, SAVI, sunelev, DOY from landsat data 
#' welev=278
#' data=landsat578(rawdata,welev=welev)
#'  #perform semi-auto simulation 
#' #Determine xyhot. Digitize polygon on the Ts map
#' modhot=hotTs(data,welev=300,extent="auto",cluster=2)
#' #determine the cold. Digitize polygon on the Ts map
#' modcold=coldTs(data,welev=275,extent="auto",cluster=2)
#' xyhot=modhot$xyhot
#' xycold=modcold$xycold
#' #use object of \code{\link{coldTs}} or \code{\link{hotTs}} in different ways
#' modsebi=sebi(data,xyhot=modhot,xycold=xycold,welev=welev,
#' Tmax=Tmax,Tmin=Tmin,RHmax=RHmax,RHmin=RHmin)
#' }
#' @seealso \code{\link{sebal}} , \code{\link{ETo}}, \code{\link{sebs}}  
#' @export
#' @rdname sebi 

sebi=function(albedo,Ts,Tmax,Tmin,RHmax=NULL,RHmin=NULL,latitude=NULL,
      NDVI,SAVI,iter.max=7,xyhot="auto",xycold="auto",DOY,
      sunelev,welev,zx=10,u=2,zomw=2,zom=NULL,LAI=NULL,DEM=NULL,
lapse=0.0065,Rn24=NULL,wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",
model="sebi",clip=NULL,folder=NULL) UseMethod ("sebi")
#' @export
#' @rdname sebi 
sebi.default=function(albedo=NULL,Ts,Tmax,Tmin,RHmax=NULL,RHmin=NULL,
latitude=NULL,NDVI,SAVI,iter.max=7,xyhot="auto",xycold="auto",DOY,
sunelev,welev,zx=10,u=2,zomw=2,zom=NULL,
LAI=NULL,DEM=NULL,lapse=0.0065,Rn24=NULL,wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",
model=model,clip=NULL,folder=NULL){

  if(is.null(folder)){
    folder="/nothing454782ghf7poi8r.hope"
  }
  
  if(inherits(folder, "landsat578")||inherits(albedo, "landsat578")){
    file.info2=TRUE
  }else{
    file.info2=file.info(folder)[["isdir"]]
  }
  if(file.info2==TRUE&&!is.na(file.info2)){
    if(is.null(welev)){
      return(message("Please provide the parameter welev: 
                   the elavation of th weather station"))  
    }
    
    if(inherits(folder, "landsat578")||inherits(albedo, "landsat578")){
      if(inherits(albedo, "landsat578")){
        mod=albedo
      }else{
        mod=folder
      }
    }else{
      mod=landsat578(folder=folder, welev=welev,clip=clip)
      folder=mod
    }
    #mod=landsat578(data=folder, welev=welev)
    albedo=mod$albedo
    Ts=mod$Ts
    NDVI=mod$NDVI
    DOY=mod$DOY
    sunelev=mod$sunelev
    SAVI=mod$SAVI
    }
  PET=NULL
  
  if(!is.null(wmo)||!is.null(airport)){
    mod=sebkc.tryCatch(mod)$value
    thisDOY=sebkc.tryCatch(mod$date)$value
    if(inherits(thisDOY, "simpleError")){
      thisDOY=DOY 
    }
    ETr24=sebkc.tryCatch(ETo(DOY=thisDOY,airport = airport,wmo=wmo,latitude=latitude,
                             z=zx,Krs =Krs,altitude=welev))$value
    if(inherits(ETr24, "simpleError")){
      return (message(paste("Invalid DOY [YYYY-mm-dd] or airport or wmo",ETr24))) 
    }
    
    Tmin=ETr24$data$Tmin
    Tmax=ETr24$data$Tmax
    Rn24=ETr24$Rn
    PET=ETr24
    
  }
uw=u
waterdensity=1000
p=1.03547033
cp=1004.14
g=9.81
k=0.41 #von Karman's constant
z1=0.1 
z2=2
vaporisation=2430000

mod=sebal(albedo=albedo,Ts=Ts,NDVI=NDVI,SAVI=SAVI,
iter.max=iter.max,xyhot=xyhot,xycold=xycold,
DOY=DOY,sunelev=sunelev,welev=welev,
zx=zx,u=uw,zomw=zomw,zom=zom,LAI=LAI,DEM=DEM,
lapse=lapse,Rn24=Rn24,model="SEBAL",folder=folder)

if(is.null(PET)){
PET=ETo(Tmax=Tmax,Tmin=Tmin,z=zx,uz=u,altitude=welev,
RHmax=RHmax,RHmin=RHmin,n=9.2,DOY=mod$DOY,latitude=latitude)
}
y=PET$y
vpd=PET$vpd
slope=PET$slope
# Per-pixel, location-dependent dry/wet limits (Menenti & Choudhury, 1993). This
# replaces the earlier formulation, which derived the limits from the single hot
# and cold anchor pixels (mod$iter.coef) - contradicting SEBI's location-dependent
# premise - and which used Rn + G for the available energy where it should be
# Rn - G (corrected below, consistent with sebal()). SEBI takes the relative
# evaporation between the dry and wet surface-temperature limits directly as the
# evaporative fraction. Rn, G, H, rah and pair are per-pixel rasters from sebal().
avail=mod$Rn-mod$G
Hdry=avail                                                       #dry limit: maximum surface temperature range
Hwet=(avail-((mod$pair*cp)/mod$rah)*(vpd/y))/(1+(slope/y))       #wet limit
message("computing EF" )
EF=1-((mod$H-Hwet)/(Hdry-Hwet))                                  #relative evaporation from the Ts range
EF[avail<=0,]=0
EF[EF<0,]=0
EF[EF>1.1,]=1.1
ET24=NULL
if(!is.null(Rn24)){
  message("computing ETa" )
  
  ET24=EF*Rn24*0.035
}

factor<-list(EF=EF,ET24=ET24,folder=folder)
factor$call<-match.call()
class(factor)<-"sebi"
factor
}