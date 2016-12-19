#' @title Simplified Surface Energy Balance
#' 
#' @description Senay etal (2011) Simplified Surface Energy Balance (SSEB) 
#' approach for estimating landscape ET
#'
#' @param TH the average of the representative 
#' 3 hot pixels selected for hot "bare" areas. 
#' If it is set to "auto" TH will be estimated from \code{\link{hotTs}} 
#' @param TC the average of representative 3 cold 
#' pixels selected from the irrigated fields. 
#' When it is set to "auto" TC will be estimated from \code{\link{coldTs}} 
#' @inheritParams sebal
#' @inheritParams ETo 
#' @inheritParams coldTs 
#' @inheritParams sebkcstack
#' @inheritParams ETohr
#' @inheritParams weather
#' @param NDVImax numeric, NDVI value for heavily green-vegetated area in the image
#' @param cc numeric, the correction coefficient if the NDVI approaches 0.0
#' @param KL numeric, the assumed lapse rate of air moving along the terrain
#' @param ETo numeric, a standardized clipped grass reference ET. when it is set to NULL,
#' it is automatically estimated from \code{\link{ETo}}
#' @param x multiplier needed  to estimate ET for tall, full cover crops 
#' such as alfalfa, corn and wheat. Default is 1.2.
#' @author George Owusu
#' @return
#'  \itemize{
#' \item{ET24:} { Actual 24 hour ET}
#' \item{EF:} { ET fraction}
#' \item{ETo:} { 24 hour reference ETo}
#' }
#'
#' @examples
#' \dontrun{
#' albedo=raster(system.file("extdata","albedo.grd",package="SEBKc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="SEBKc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="SEBKc"))
#' #minimal data
#'modsseb=sseb(Ts=Ts,TH="full",TC="full",sunelev=50,NDVI=NDVI,DEM=NULL,
#'albedo=albedo,NDVImax=0.7,cc=0.65,KL=0.0065,
#'ETo="auto",x=1.2,Tmax=31,Tmin=28,zx=10,
#'u=2,DOY=37,latitude=5.6,n=NULL,RHmax=NULL,RHmin=NULL)
#' 
#' #uaing the folder parameters
#' folder=system.file("extdata","stack",package="SEBKc")
#' ssebauto=sseb(folder=folder,welev=317,Tmax=31,Tmin=28, latitude=5.6)
#' 
#' #Another interactive example
#' #get Ts, Albedo, NDVI, SAVI, sunelev, DOY from landsat data 
#' welev=278
#' data=landsat578(rawdata,welev=welev)
#' #perform semi-auto simulation 
#' #Determine xyhot. Digitize polygon on the Ts map
#' modhot=hotTs(data,welev=300,extent="digitize",cluster=2)
#' #determine the cold. Digitize polygon on the Ts map
#' modcold=coldTs(data,welev=275,extent="digitize",cluster=2)
#' xyhot=modhot$Tshot
#' xycold=modcold$Tscold
#' #use object of \code{\link{coldTs}} or \code{\link{hotTs}} in different ways
#' modsseb=sseb(Ts=Ts,TH=xyhot,TC=modcold,sunelev=50,NDVI=NDVI,DEM=NULL,
#' albedo=albedo,NDVImax=0.7,cc=0.65,KL=0.0065,
#' ETo="auto",x=1.2,Tmax=31,Tmin=28,zx=10,
#' u=2,DOY=37,latitude=5.6,n=NULL,
#' RHmax=NULL,RHmin=NULL)
#' }
#' @references 
#' Senay, G. B., Budde, M. E., & Verdin, J. P. 2011. 
#' Enhancing the Simplified Surface Energy Balance (SSEB) approach 
#' for estimating landscape ET: Validation with the METRIC model. 
#' Agricultural Water Management, 98(4): 606-618.
#' @export
#' @rdname sseb 

sseb=function(Ts,welev=NULL,TH="auto",TC="auto",sunelev=NULL,NDVI=NULL,DEM=NULL,
              albedo=NULL,NDVImax=0.7,cc=0.65,KL=0.0065,
              ETo=NULL,wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",
              x=1.0,Tmax=NULL,Tmin=NULL,zx=NULL,
              u=2,DOY=NULL,latitude=NULL,n=NULL,
              RHmax=NULL,RHmin=NULL,clip=NULL,folder=NULL) UseMethod ("sseb")
#' @export
#' @rdname sseb
sseb.default=function(Ts=NULL,welev=NULL,TH="auto",TC="auto",sunelev=NULL,
  NDVI=NULL,DEM=NULL,albedo=NULL,NDVImax=0.7,cc=0.65,KL=0.0065,
ETo=NULL,wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",
x=1.2,Tmax=NULL,Tmin=NULL,zx=NULL,
u=2,DOY=NULL,latitude=NULL,n=NULL,
RHmax=NULL,RHmin=NULL,clip=NULL,folder=NULL){
  
  if(is.null(folder)){
    folder="/nothing454782ghf7poi8r.hope"
  }
  
  if(class(folder)=="landsat578"||class(Ts)=="landsat578"){
    file.info2=TRUE
  }else{
    file.info2=file.info(folder)[["isdir"]]
  }
  if(file.info2==TRUE&&!is.na(file.info2)){
    if(is.null(welev)){
      return(print("Please provide the parameter welev: 
                     the elavation of th weather station"))  
    }
    
    if(class(folder)=="landsat578"||class(Ts)=="landsat578"){
      if(class(Ts)=="landsat578"){
        mod=Ts
      }else{
        mod=folder
      }
    }else{
      mod=landsat578(folder=folder, welev=welev,clip=clip)
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
  if(is.null(ETo)){
    
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
      ETo=ETr24$ETo
    
    }
  }
  
Tx=Ts    
multiplier=x   
threshold=NDVImax  
  

if(!is.numeric(DOY)&&!is.null(DOY)){
  DOY2=DOY
  DOY=strptime(DOY,"%d-%m-%Y")$yday+1
  if(is.na(DOY)){
    DOY=strptime(DOY2,"%Y-%m-%d")$yday+1
  }
}
NDVI[NDVI<0,]=0
NDVI[NDVI>1,]=1

if(!is.null(DEM)){
Tx=Tx+(KL*DEM)
}

if(class(TH)=="hotTs"){
  TH=TH$Tshot
}

if(class(TC)=="coldTs"){
  TC=TC$Tscold
}

testcold=NULL

if((TC=="auto"||TC=="Auto"||TC=="AUTO"||TC=="full"||TC=="Full"||TC=="FULL")&&!is.numeric(TC)){
  TC=tolower(TC)
  testcold=TRUE
}
else{
  if(length(TC)>2){
    testcold=TRUE 
  }else{
    TC=TC
  }
}
testhot=NULL
if((TH=="auto"||TH=="Auto"||TH=="AUTO"||TH=="FULL"||TH=="full"||TH=="Full")&&(!is.numeric(TH))){
  TH=tolower(TH)
  testhot=TRUE
}else{
  if(length(TH)>2){
    testhot=TRUE 
  }else{
    TH=TH
  }
}

if(!is.null(testhot)){
  modhot=hotTs(Ts=Tx,NDVI=NDVI,albedo=albedo,DEM=DEM,
               extent=TH,upper=0.80,lower=0.1)
  TH=modhot$TH
}
if(!is.null(testcold)){
  modcold=coldTs(Ts=Tx,sunangle=sunelev,NDVI=NDVI,albedo=albedo,
                 DEM=DEM,extent=TC,upper=0.95,lower=0.1)
  TC=modcold$TC
}
print ("computing ETf" )
ETf=(TH-Tx)/(TH-TC)

if(!is.null(NDVI)){
ETf=(((1-cc)*(NDVI/threshold)+(cc)))*ETf
}
if(is.null(ETo)||ETo=="auto"){
if(is.null(Tmax)||is.null(Tmin)||is.null(DOY)
||is.null(latitude)){
return(print("The following variables are needed:
Tmax,Tmin,DOY,latitude for estimation of ETo"))

}
  
PET=ETo(Tmax=Tmax,Tmin=Tmin,z=zx,uz=u,
altitude=welev,RHmax=RHmax,RHmin=RHmin,
n=n,DOY=DOY,latitude=latitude)

ETo=PET$ETo

}

ETf[ETf>1.1,]=1.1
ETf[ETf<0,]=0

print ("computing ETm" )

ETm=multiplier*ETo
print ("computing ETa" )

ETa=ETf*ETm
factor<-list(ET24=ETa, EF=ETf,ETo=ETo,folder=folder)
factor$call<-match.call()
class(factor)<-"ssebi"
factor
}

