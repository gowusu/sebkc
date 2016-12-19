#' @title  Automatic computation of hot pixel Radiometric Surface Temperature (Ts)

#' @description This function computes temperature of a hot pixel based on NDVI, 
#' Surface Temperature, DEM and albedo. It plots histograms indicating
#' the selected values
#' @param upper The quantile probability at which Ts should be included 
#' in the selection project. It ranges from 0-1. The default is 0.8
#' @param lower The quantile probability at which NDVI should be included 
#' in the slection project.The default is 0.2
#' @inheritParams sebal
#' @inheritParams coldTs2
#' @inheritParams coldTs
#' @details The function first divides the selected area into a 
#' number of clusters. The mean and standard deviations are
#' computed. The cluster means with the lowest standard deviations
#' and qualifying for upper (Ts) and lower (NDVI) limits 
#' (see above) are selected. The provided upper and lower probabilities
#'  are also used to select the final Ts.
#' @seealso \code{\link{coldTs}},  \code{\link{hotTs}}, \code{\link{coldT2}} 
#' @return 
#' \itemize{
#' \item{Ts:} { The selected or input Ts}
#' \item{Tshot:} { The Temperature of hot pixel}
#' \item{TX:} { The average temperature of the hot pixels}
#' \item{x:} { The x cordinate of the hot pixel}
#' \item{y:} { The y cordinate of the hot pixel}
#' \item{xyhot: } {combined x and y coordinates}
#' \item{candidates:}{ The similar candidates' pixels that can be used}
#' }
#' @examples
#' \dontrun{
#'  folder=system.file("extdata","stack",package="SEBKc")
#'   modhot=hotTs2(folder=folder,welev=170,extent="auto",cluster=3)
#'   
#'    #use object of  landsat578
#' folder=system.file("extdata","stack",package="SEBKc")
#' data=landsat578(data=folder, welev=170)
#' modfhot=hotTs2(folder=data,welev=170,extent="auto",cluster=3)
#' 
#' #Or Ts = data
#' modhotTs2=hotTs2(data,welev=170,extent="auto",cluster=3)

#' albedo=raster(system.file("extdata","albedo.grd",package="SEBKc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="SEBKc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="SEBKc"))
#' #simulate  with default parameters without albedo
#' modhot=hotTs2(Ts=Ts,NDVI=NDVI,extent="auto",cluster=3)
#' #simulate with albedo but set extent to "auto" and cluster to 7
#' modhot=hotTs2(Ts=Ts,NDVI=NDVI,albedo=albedo,cluster=7,
#' extent="auto")
#' }
#' @author George Owusu
#' @references Allen, R. G., Burnett, B., Kramber, W., Huntington, J., 
#' Kjaersgaard, J., Kilic, A., Kelly, C., & Trezza, R. 2013. Automated 
#' Calibration of the METRIC-Landsat Evapotranspiration Process. 
#' JAWRA Journal of the American Water Resources Association, 49(3): 563-576.
#' @export
#' @rdname hotTs2
#'
 
hotTs2<-function(Ts,NDVI,albedo=NULL,sunangle=NULL,DEM=NULL,
                extent="interactive",upper=0.95,
lower=0.8,NDVI.probs=0.1,plot=TRUE,layout="portrait",draw="poly",
folder=NULL,welev=NULL) UseMethod ("hotTs2")
#' @export
#' @rdname hotTs2
hotTs2.default=function(Ts=NULL,NDVI,albedo=NULL,sunangle=NULL,DEM=NULL,
cluster=10,extent="interactive",upper=0.95,
lower=0.8,NDVI.probs=0.2,plot=TRUE,layout="portrait",draw="poly",folder=NULL,welev=NULL){
  
  if(is.null(folder)){
    folder="nothing454782ghf7poi8r.hope"
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
      mod=landsat578(data=folder, welev=welev)
    }
    
    albedo=mod$albedo
    Ts=mod$Ts
    NDVI=mod$NDVI
    sunangle=mod$sunelev
    }
  
  if(is.null(sunangle)){
    return (print("sunangle is needed"))
  }
  
ext=extent
albedothre=0.02
if(class(NDVI)!="RasterLayer"){
NDVI=raster(NDVI)
}
if(class(Ts)!="RasterLayer"){
Ts=raster(Ts)
}
if(class(albedo)!="RasterLayer"&&!is.null(albedo)){
albedo=raster(albedo)
}

if(is.numeric(ext[1])){
Ts <-crop(Ts, ext)
}else{
if(ext=="auto"||ext=="full"||extent=="Full"||extent=="FULL"){
  if(ext=="full"||extent=="Full"||extent=="FULL"){
    ext=extent(Ts)
  }else{
extents=extent(Ts)
xrange=extents[2]-extents[1]
yrange=extents[4]-extents[3]
H=0.2
V=0.4
if(layout=="H"||layout=="horinzontal"||layout=="landscape"){
H=V
V=H
}
x40=V*xrange
xmin=extents[1]+x40
xmax=extents[2]-x40
y20=H*yrange
ymin=extents[3]+y20
ymax=extents[4]-y20
extents2=c(xmin,xmax,ymin,ymax)
Ts=crop(Ts,extents2)
ext=extent(Ts)
}
}else{
if(draw=="poly"||draw=="line"){
plot(Ts,main="click and draw a line or polygon on a plotand right-click and
     select 'stop'. If you are using RStudio click on Finish to stop")
ext=drawPoly()
}else{
plot(Ts,main="click on two places on this map to select Area of Interest.
     If you are using RStudio click on Finish to stop")
ext <- drawExtent() #draw a box by clicking
}
Ts <-crop(Ts, ext)
}
}
if(!is.null(DEM)){
DEM = crop(DEM, ext)
Ts=Ts+(0.0065*DEM)
}
NDVI <- crop(NDVI, ext)
NDVI2=NDVI
if(!is.null(albedo)){
albedo2=crop(albedo, ext)
}
hotTs2=Ts
upper1=quantile(hotTs2,probs=upper)
lower1=quantile(hotTs2,probs=lower)
hotTs2[hotTs2>upper1[[1]],]=NA
hotTs2[hotTs2<lower1[[1]],]=NA
#NDVI2[is.na(hotTs2),]=NA
hotTs2[NDVI2<quantile(NDVI2,probs=NDVI.probs)[[1]],]=NA
p=data.frame(rasterToPoints(hotTs2))
hotTs2xy<-p[p$layer>0,]
#hotTs2xy=p[,c("x","y")]
p3=p
if(nrow(p)>1){
p2= p[order(p$layer),]
p3=p2[round(nrow(p2)/2),] 
}
x=p3$x
y=p3$y
ndvi=getValues(NDVI2)[cellFromXY(NDVI2,c(x,y))]

if(plot==TRUE||plot==T){
par(mfrow=c(1,2))
hist(Ts,main="Ts[Hot/Dry]",xlab="Temperature")
abline(v=p3,col="red")
hist(NDVI,main="NDVI",xlab="NDVI")
abline(v=ndvi,col="red")
}
hotTs2=p3$layer

factor<-list(x=x,y=y,xyhot=c(x,y),Ts=Ts,NDVI=ndvi,Tshot=hotTs2,
             p=p,candidates=p,TC=mean(p$layer))

factor$call<-match.call()

class(factor)<-"hotTs"
factor
}

