
#' @title SEBS Estimation of Evapotranspiration Fraction
#'
#' @description \code{sebs} computes Su (2002) 
#' Surface Energy Balance System (SEBS) for estimation of turbulent heat fluxes
#' @inheritParams sebal
#' @inheritParams ETo
#' @inheritParams ETohr
#' @inheritParams weather
#' @inheritParams sebkcstack
#' @details Function \code{\link{ETo}} will be automatically used to 
#' estimate  parameters y (psychrometric constant), 
#' slope (slope of saturation vapour pressure curve) and 
#' vpd (Vapour pressure deficit). For information on data input procedure
#' see 'details' at \code{\link{sebal}}
#' @seealso \code{\link{sebal}} , \code{\link{ETo}}, \code{\link{sebi}}  
#' @references 
#' Su, Z. 2002. The Surface Energy Balance System (SEBS) for estimation 
#' of turbulent heat fluxes. Hydrol. Earth Syst. Sci., 6(1): 85-100.
#' @return EF Evapotranspiration Fraction
#' @author George Owusu
#'
#' @examples
#' \dontrun{
#' #define Input data
#' albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="sebkc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
#' LAI=raster(system.file("extdata","LAI.grd",package="sebkc"))
#' Tmax=31
#' Tmin=28
#' RHmax=84,
#' RHmin=63
#' #perform simulation
#' modsebs=sebs(albedo=albedo,Ts=Ts,Tmax=Tmax,Tmin=Tmin,RHmax=RHmax,
#' RHmin=RHmin,NDVI=NDVI,SAVI=NULL,iter.max=7,xyhot="full",
#' xycold="full",DOY=37,sunelev=50.71154048,welev=317.1,zx=10,
#' u=2,zomw=2,zom=NULL,LAI=LAI,Rn24=NULL,model="SEBS")
#'
#' 
#' #use original landsat 7 data by specifying folder path
#' folder=system.file("extdata","stack",package="sebkc")
#' modauto=sebs(folder = folder,welev = 380,Tmax=Tmax,
#' Tmin=Tmin,RHmax=RHmax,RHmin=RHmin)
#' #plot ET fraction
#' plot(modauto$EF)
#' 
#' #get Ts, Albedo, NDVI, SAVI, sunelev, DOY from landsat data 
#' welev=278
#' data=landsat578(rawdata,welev=welev)
#'  #perform semi-auto simulation 
#' #Determine xyhot. Digitize polygon on the Ts map
#' modhot=hotTs(data,welev=300,extent="digitize",cluster=2)
#' #determine the cold. Digitize polygon on the Ts map
#' modcold=coldTs(data,welev=275,extent="digitize",cluster=2)
#' xyhot=modhot$xyhot
#' xycold=modcold$xycold
#' #use object of \code{\link{coldTs}} or \code{\link{hotTs}} in different ways
#' modsebs=sebs(data,xyhot=modhot,xycold=xycold,welev=welev,
#' Tmax=Tmax,Tmin=Tmin,RHmax=RHmax,RHmin=RHmin)
#'  }
#' @export
#' @rdname sebs 

sebs=function(albedo,Ts,Tmax,Tmin,RHmax=NULL,RHmin=NULL,n=NULL,
      NDVI,SAVI,iter.max=7,xyhot,xycold,DOY,
      sunelev,welev,zx,u,zomw,zom=NULL,LAI=NULL,DEM=NULL,
      lapse=0.0065,Rn24=NULL,wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",
      latitude=NULL,model="SEBS",clip=NULL,folder=NULL) UseMethod ("sebs")
#' @export
#' @rdname sebs 
sebs.default=function(albedo=NULL,Ts,Tmax,Tmin,RHmax=NULL,RHmin=NULL,n=NULL,
NDVI,SAVI,iter.max=7,xyhot="auto",xycold="auto",DOY=NULL,
sunelev,welev,zx=10,u=2,zomw=2,zom=NULL,LAI=NULL,DEM=NULL,lapse=0.0065,Rn24=NULL,
wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",latitude=NULL,
model="SEBS",clip=NULL,folder=NULL){
  
  if(is.null(folder)){
    folder="/nothing454782ghf7poi8r.hope"
  }
  
  if(class(folder)=="landsat578"||class(albedo)=="landsat578"){
    file.info2=TRUE
  }else{
    file.info2=file.info(folder)[["isdir"]]
  }
  if(file.info2==TRUE&&!is.na(file.info2)){
    if(is.null(welev)){
      return(print("Please provide the parameter welev: 
                   the elavation of th weather station"))  
    }
    
    if(class(folder)=="landsat578"||class(albedo)=="landsat578"){
      if(class(albedo)=="landsat578"){
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
      if(class(thisDOY)[1]=="simpleError"){
        thisDOY=DOY 
      }
      ETr24=sebkc.tryCatch(ETo(DOY=thisDOY,airport = airport,wmo=wmo,latitude=latitude,
                               z=zx,Krs =Krs,altitude=welev))$value
      if(class(ETr24)[1]=="simpleError"){
        return (print(paste("Invalid DOY [YYYY-mm-dd] or airport or wmo",ETr24))) 
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

if(!is.numeric(DOY)&&!is.null(DOY)){
  DOY2=DOY
  DOY=strptime(DOY,"%d-%m-%Y")$yday+1
  if(is.na(DOY)){
    DOY=strptime(DOY2,"%Y-%m-%d")$yday+1
  }
}

mod=sebal(albedo=albedo,Ts=Ts,NDVI=NDVI,SAVI=SAVI,
iter.max=iter.max,xyhot=xyhot,xycold=xycold,
DOY=DOY,sunelev=sunelev,welev=welev,
zx=zx,u=uw,zomw=zomw,zom=zom,LAI=LAI,DEM=DEM,
lapse=lapse,Rn24=Rn24,model="SEBAL",folder=folder)

if(is.null(PET)){
PET=ETo(Tmax=Tmax,Tmin=Tmin,z=zx,uz=u,altitude=welev,
RHmax=RHmax,RHmin=RHmin,n=n,DOY=mod$DOY,latitude=1)
}
data=mod$iter.coef[nrow(mod$iter.coef),]
y=PET$y
vpd=PET$vpd
slope=PET$slope
Hdry=data$Rn_hot-data$G_hot
Hwet=((data$Rn_cold-data$G_cold)-((data$paircold*cp)/data$rah_cold)*(vpd/y))/((1+slope/y))
print ("computing EF" )
EFr=1-((mod$H-Hwet)/(Hdry-Hwet))
EF=1.1*((EFr*data$LEcold)/(mod$Rn-mod$G))
EF[EF<0,]=0
EF[EF>1.1,]=1.1

ET24=NULL
#print(Rn24)
if(!is.null(Rn24)){
  print ("computing ETa" )
  ET24=EF*Rn24*0.035
}
factor<-list(EF=EF,ET24=ET24,folder=folder)
factor$call<-match.call()
class(factor)<-"sebs"
factor
}




