#' @title Water Deficit Index (WDI) 
#' @description Moran et. al (1994)  
#' Estimation of crop water deficit using surface-air 
#' temperature and spectral vegetation index.
#' @param Ta numeric. Air temperature. If not available you can use 
#' @param Tsmax numeric. Maximum sensor surface temperature [K]
#' @param Tsmin numeric. Minimum sensor surface temperature [K]
#' \code{\link{sebal}} to get pixel based Ta.
#' @inheritParams sebal
#' @inheritParams ETo 
#' @inheritParams coldTs
#' @author George Owusu
#' @details Note that surface temperature (Tsmin and Tsmax)
#'  are different from air temperature (Tmin and Tmax) in
#' \code{\link{ETo}}, \code{\link{sseb}}, \code{\link{sebs}}
#' 
#' You can get Ta from  \code{\link{sebal}} Ta
#' @seealso \code{\link{sseb}} and \code{\link{sebi}}
#' @references Moran, M. S., Clarke, T. R., Inoue, Y., & Vidal, A. 1994.
#' Estimating crop water deficit using the relation between surface-air 
#' temperature and spectral vegetation index. Remote Sensing of Environment,
#'  49(3): 246-263.
#' @return
#'  \itemize{
#' \item{EF:} { standardised wdi excluding negatives and one plus}
#' \item{index:} { raw wdi including negatives and more than one}
#' }
#' @examples 
#'  \dontrun{
#' albedo=raster(system.file("extdata","albedo.grd",package="SEBKc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="SEBKc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="SEBKc"))
#' modWDI=wdi(Ts=Ts,Ta=299,NDVI=NDVI,albedo=NULL,Tsmax="auto",
#' Tsmin="auto",sunelev=50)
#' plot(modWDI$EF)
#' 
#' #Using folder parameter
#' folder=system.file("extdata","stack",package="SEBKc")
#' wdiauto=wdi(folder=folder,welev=317,Tsmax=31,Tsmin=28,Ta=290)
#' 
#' #' #Another interactive example
#' #get Ts, Albedo, NDVI, SAVI, sunelev, DOY from landsat data 
#' welev=278
#' data=landsat578(rawdata,welev=welev)
#' #perform semi-auto simulation 
#' #Determine xyhot. Digitize polygon on the Ts map
#' modhot=hotTs(data,welev=300,extent="digitize",cluster=2)
#' #determine the cold. Digitize polygon on the Ts map
#' modcold=coldTs(data,welev=275,extent="digitize",cluster=2)
#' Tsmax=modhot$Tshot
#' Tsmin=modcold$Tscold
#' #use object of \code{\link{coldTs}} or \code{\link{hotTs}} in different ways
#' modWDI2=wdi(data,Ta=299,Tsmax=Tsmax,Tsmin=modcold,sunelev=50)
#' plot(modWDI2$EF)
#' }
#' @export
#' @rdname wdi
wdi=function(Ts,Ta,NDVI,albedo=NULL,Tsmax="auto",Tsmin="auto",
             sunelev=NULL,folder=NULL,welev=NULL,Rn24=NULL)
  UseMethod ("wdi")
#' @export
#' @rdname wdi
wdi=function(Ts=NULL,Ta,NDVI,albedo=NULL,Tsmax="auto",Tsmin="auto",
             sunelev=NULL,folder=NULL,welev=NULL,Rn24=NULL){
  if(is.null(folder)){
    folder="/nothing454782ghf7poi8r.hope"
  }
  
  if(class(folder)=="landsat578"||class(Ts)=="landsat578"){
    file.info2=TRUE
  }else{
    file.info2=file.info(folder)[["isdir"]]
  }
  if(file.info2==TRUE&&!is.na(file.info2)){
    if(is.null(welev)&&class(Ts)!="landsat578"){
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
      mod=landsat578(folder=folder, welev=welev)
    }
    #mod=landsat578(data=folder, welev=welev)
    albedo=mod$albedo
    Ts=mod$Ts
    NDVI=mod$NDVI
    DOY=mod$DOY
    sunelev=mod$sunelev
    SAVI=mod$SAVI
    welev=mod$welev
    }
  NDVI[NDVI<0,]=0
  NDVI[NDVI>1,]=1
  
  #maximun Ts
  
  if(class(Tsmax)=="hotTs"){
    Tsmax=Tsmax$Tshot
  }
  
  if(class(Tsmin)=="coldTs"){
    Tsmin=Tsmin$Tscold
  }
  
  if(Tsmax=="auto"&&!is.numeric(Tsmax)){
    modhot=hotTs(Ts=Ts,NDVI=NDVI,albedo=albedo,DEM=NULL,cluster=10,
                 extent="auto")
    Tsmax=modhot$Tshot 
  }
  
  if(Tsmin=="auto"&&!is.numeric(Tsmin)){
    modcold=coldTs(Ts=Ts,NDVI=NDVI,sunangle=sunelev,albedo=albedo,
                   DEM=NULL,cluster=10, extent="auto",upper=0.95,lower=0.1)
    Tsmin=modcold$Tscold 
  }
  #minimun Ts
  
  Tacold=Tahot=Ta
  if(class(Ta)=="RasterLayer"){
    Tacold=getValues(Ta)[cellFromXY(Ta,c(modcold$x,modcold$y))]
    Tahot=getValues(Ta)[cellFromXY(Ta,c(modhot$x,modhot$y))]
  }
  
  
  wdi=((Tsmin-Tacold)-(Ts-Ta))/((Tsmin-Tacold)-(Tsmax-Tahot))
  index=wdi
  wdi[wdi>1,]=1
  wdi[wdi<0,]=0
  
  ET24=NULL
  if(!is.null(Rn24)){
    ET24=(wdi*Rn24)/2.45
  }
  factor<-list(EF=wdi, index=index,ET24=ET24,folder=folder)
  factor$call<-match.call()
  class(factor)<-"wdi"
  factor
}
