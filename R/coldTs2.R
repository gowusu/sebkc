#' @title  Optional Automatic computation of cold pixel Radiometric Surface Temperature
#'
#' @description This function computes temperature of a cold pixel based on NDVI, 
#' Surface Temperature, DEM and albedo. The main difference with \code{\link{coldTs}} is that preference 
#' is given to temperature selection before considering NDVI and albedo selection. Moreover, this
#' is a simple end member selection where no clustering is done.  
#' @param sunangle sun angle (degree) above the horizon.
#' @param NDVI.probs The probability of selecting quantile NDVI that is above the threshold. 
#' The value of 0.6 means only NDVI values that is above 60% quantile will be included in
#' the analyses.
#' @param extent A subset of the raster file than include 
#' c(xmin, xmax; serastercond row: ymin, ymax). You can set it to "interactive"
#' where you can digitise polygon. You can set to "auto" where it will be
#' automatically generated. See layout parameter below for more details.
#' see \code{\link[raster]{extent}} 
#' @param upper The quantile probability at which NDVI should be included 
#' in the selection project. It ranges from 0-1. The default is 0.95
#' @param lower The quantile probability at which Ts should be included 
#' in the selection project.The default is 0.2
#' @param plot logical:TRUE or FALSE. To indicate if histograms or triangle 
#' should be plotted
#' @param layout character. It takes "portrait" or "landscape". If extent is set to 
#' "auto" the middle part of the image will be used so that NAs can be 
#' avoided. You indicate whether the automatic selection should indicate
#' west-east (landscape) or north-south direction (portrait)
#' @param draw interactive. The method of digitisation if extent is set "auto". 
#' It takes "poly" or "rect"
#' @seealso \code{\link{coldTs}},  \code{\link{hotTs}}, \code{\link{hotTs2}}  
#' @inheritParams sebal
#' @inheritParams coldTs
#' @details The function first divides the selected area into a  
#' number of clusters. The mean and standard deviations are
#' computed. The cluster means with the lowest standard devaitions 
#' and qualifying for upper (NDVI) and lower (Ts) limits 
#' (see above) are selected. The provided upper and lower probabilities
#'  are also used to select the final Ts.
#' @return 
#' \itemize{
#' \item{Ts:} { The selected or input Ts [K]}
#' \item{Tscold:} { The Temperature of cold pixel [K] }
#' \item{TC:} { The average temperature of the cold pixels [K]}
#' \item{x:} { The x coordinate  of the cold pixel}
#' \item{y:} { The y coordinate  of the cold pixel}
#' \item{xycold: } {combined x and y coordinates}
#' \item{candidates:}{ The similar candidates' pixels that can be used}
#' }
#' @examples
#' \dontrun{
#' #using landsat folder
#' folder=system.file("extdata","stack",package="sebkc")
#' modcold=coldTs2(folder=folder,welev=170,extent="auto")
#'
#' #use object of  
#' folder=system.file("extdata","stack",package="sebkc")
#' data=landsat578(data=folder, welev=362)
#' modcold=coldTs2(folder=data,welev=170,extent="auto",cluster=3)
#' 
#' modcoldTs2=coldTs2(data,welev=170,extent="auto",cluster=3)
#' 
#' #using input different input data
#' albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="sebkc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
#' #simulate with default parameters without albedo
#' modcold=coldTs2(Ts=Ts,NDVI=NDVI,sunangle=50.7)
#' #simulate with albedo but set extent to "auto" and cluster to 7
#' modcold=coldTs2(Ts=Ts,NDVI=NDVI,albedo=albedo,sunangle=50,cluster=7,
#' extent="auto")
#' }
#' @author George Owusu
#' @references Allen, R. G., Burnett, B., Kramber, W., Huntington, J., 
#' Kjaersgaard, J., Kilic, A., Kelly, C., & Trezza, R. 2013. Automated 
#' Calibration of the METRIC-Landsat Evapotranspiration Process. 
#' JAWRA Journal of the American Water Resources Association, 49(3): 563-576.
#' @export
#' @rdname coldTs2
 
coldTs2<-function(Ts,NDVI,albedo=NULL,sunangle=NULL,DEM=NULL,
                extent="interactive",upper=0.2,
lower=0.1,NDVI.probs=0.95,plot=TRUE,layout="portrait",draw="poly",
folder=NULL,welev=NULL) UseMethod ("coldTs2")
#' @export
#' @rdname coldTs2
coldTs2.default=function(Ts=NULL,NDVI,albedo=NULL,sunangle=NULL,DEM=NULL,
cluster=10,extent="interactive",upper=0.2,
lower=0.1,NDVI.probs=0.6,plot=TRUE,layout="portrait",draw="poly",folder=NULL,welev=NULL){
  
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
coldTs2=Ts

upper1=quantile(coldTs2,probs=upper)
lower1=quantile(coldTs2,probs=lower)
coldTs2[coldTs2>upper1[[1]],]=NA
coldTs2[coldTs2<lower1[[1]],]=NA
abeldo_thresh=(0.001343*sunangle) +(0.3281* exp(-0.0188*sunangle))-albedothre
coldTs2[NDVI2>quantile(NDVI2,probs=NDVI.probs)[[1]],]=NA

if(!is.null(albedo)&&cellStats(albedo2, stat='max', na.rm=TRUE, asSample=TRUE)>abeldo_thresh){
albedo2[is.na(coldTs2),]=NA
coldTs2[albedo2<=abeldo_thresh,]=NA
}
#NDVI2[is.na(coldTs2),]=NA

p=data.frame(rasterToPoints(coldTs2))
coldTs2xy<-p[p$layer>0,]
#coldTs2xy=p[,c("x","y")]
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
hist(Ts,main="Ts[Cold/Wet]",xlab="Temperature")
abline(v=p3,col="red")
hist(NDVI,main="NDVI",xlab="NDVI")
abline(v=ndvi,col="red")
}
coldTs2=p3$layer

factor<-list(x=x,y=y,xycold=c(x,y),Ts=Ts,NDVI=ndvi,Tscold=coldTs2,
             p=p,candidates=p,TC=mean(p$layer))

factor$call<-match.call()

class(factor)<-"coldTs"
factor
}

