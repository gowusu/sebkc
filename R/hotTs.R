#' @title Automatic computation of hot pixel Radiometric Surface Temperature (Ts)

#' @description This function computes temperature of a hot pixel based on NDVI, 
#' Surface Temperature, DEM and albedo. The main difference with \code{\link{hotTs2}} is that preference 
#' is given to temperature selection before considering NDVI and albedo selection. Moreover, this
#' is a simple end member selection where no clustering is done.
#' @param upper The quantile probability at which Ts should be included 
#' in the selection project. It ranges from 0-1. The default is 0.8
#' @param lower The quantile probability at which NDVI should be included 
#' in the selection project. The default is 0.2
#' @inheritParams sebal
#' @inheritParams coldTs2
#' @inheritParams ETohr
#' @inheritParams weather
#' @inheritParams sebkcstack
#' @details The function first divides the selected area into a 
#' number of clusters. The mean and standard deviations are
#' computed. The cluster means with the lowest standard deviations
#' and qualifying for upper (Ts) and lower (NDVI) limits 
#' (see above) are selected. The provided upper and lower probabilities
#'  are also used to select the final Ts.
#' @seealso \code{\link{coldTs}} 
#' @return 
#' \itemize{
#' \item{Ts:} { The selected or input Ts}
#' \item{Tshot:} { The Temperature of hot pixel}
#' \item{TX:} { The average temperature of the hot pixels}
#' \item{x:} { The x coordinate  of the hot pixel}
#' \item{y:} { The y coordinate  of the hot pixel}
#' \item{xyhot: } {combined x and y coordinates}
#' \item{candidates:}{ The similar candidates' pixels that can be used}
#' }
#' @import raster
#' @examples
#' \dontrun{
#'  folder=system.file("extdata","stack",package="sebkc")
#'   modhot=hotTs(folder=folder,welev=170,extent="auto",cluster=3)
#'   
#'    #use object of  landsat578
#' folder=system.file("extdata","stack",package="sebkc")
#' data=landsat578(data=folder, welev=170)
#' modfhot=hotTs(folder=data,welev=170,extent="auto",cluster=3)
#' 
#' #Or Ts = data
#' modhotTs=hotTs(data,welev=170,extent="auto",cluster=3)

#' albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="sebkc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
#' #simulatte with default parameters without albedo
#' modhot=hotTs(Ts=Ts,NDVI=NDVI,extent="auto",cluster=3)
#' #simulate with albedo but set extent to "auto" and cluster to 7
#' modhot=hotTs(Ts=Ts,NDVI=NDVI,albedo=albedo,cluster=7,
#' extent="auto")
#' }
#' @author George Owusu
#' @references Allen, R. G., Burnett, B., Kramber, W., Huntington, J., 
#' Kjaersgaard, J., Kilic, A., Kelly, C., & Trezza, R. 2013. Automated 
#' Calibration of the METRIC-Landsat Evapotranspiration Process. 
#' JAWRA Journal of the American Water Resources Association, 49(3): 563-576.
#' @export
#' @rdname hotTs
#'
hotTs=function(Ts,NDVI,albedo=NULL,DEM=NULL,
      cluster=8,extent="interactive",upper=0.80,
    lower=0.2,plot=TRUE,layout="portrait",draw="poly",folder=NULL,welev=NULL,clip=NULL) 
  UseMethod ("hotTs")
#' @export
#' @rdname hotTs
hotTs.default=function(Ts=NULL,NDVI,albedo=NULL,DEM=NULL,
cluster=8,extent="interactive",upper=0.80,
lower=0.2,plot=TRUE,layout="portrait",draw="poly",folder=NULL,welev=NULL,clip=NULL){
  print("Selecting hot pixel")
  if(is.null(Ts)&&!is.null(folder)){
    #Ts=folder 
  }
  
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
  
albedothre=0.02
ext=extent

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
      #print("yes")
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
plot(Ts,main="Click and draw a line or polygon on a plot and right-click and 
     select 'stop' If you are using RStudio click on Finish to stop")
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
#Ts[Ts<0]=NA
}
NDVI2 <- crop(NDVI, ext)
if(!is.null(albedo)){
albedo2=crop(albedo, ext)
}
Ts.kmeans <- kmeans(na.omit(Ts[]), cluster, iter.max = 100, nstart = 3)
kmeansraster<-raster(Ts)
kmeansraster[]<-Ts.kmeans$cluster
Tsmean=zonal(Ts, kmeansraster, 'mean') 
TsSD=zonal(Ts, kmeansraster, 'sd') 
#TsCV=zonal(Ts, kmeansraster, fun='cv') 
#TsCV15=TsCV
#TsCV15[,2]=TsCV15[,2]*0.15

NDVImean=zonal(NDVI2, kmeansraster, 'mean') 
NDVISD=zonal(NDVI2, kmeansraster, 'sd') 
#NDVICV=zonal(NDVI2, kmeansraster, fun='cv') 
#NDVICV15=NDVICV
#NDVICV15[,2]=NDVICV15[,2]*0.15

#substract sd from the mean
NDVImean_sd=NDVImean
NDVImean_sd[,2]=NDVImean[,2]-NDVISD[,2]

#hot pixels AOI selection
NDVImax=NDVImean_sd[NDVImean_sd[,2]==max(NDVImean_sd[,2]),]
NDVImin=NDVImean_sd[NDVImean_sd[,2]==min(NDVImean_sd[,2]),]
NDVI95=quantile(NDVImean_sd[,2],  probs = lower)
NDVI95b=NDVImean_sd[NDVImean_sd[,2]<NDVI95,]
if(length(NDVI95b[1])==1){
Tsmean95=Tsmean[Tsmean[,1] %in%(NDVI95b[1]),]
Tsmean95SD=TsSD[TsSD[,1] %in%(NDVI95b[1]),]
Tsmean95_sd=Tsmean95
Tsmean95_sd[2]=Tsmean95[2]-Tsmean95SD[2]
Tsmean2_sda=quantile(Tsmean95_sd[2],  probs = upper)
Tsmean2b=Tsmean95_sd

}else
{
Tsmean95=Tsmean[Tsmean[,1] %in%(NDVI95b[,1]),]
Tsmean95SD=TsSD[TsSD[,1] %in%(NDVI95b[,1]),]
Tsmean95_sd=Tsmean95
Tsmean95_sd[,2]=Tsmean95[,2]-Tsmean95SD[,2]
Tsmean2_sda=quantile(Tsmean95_sd[,2],  probs = upper)
Tsmean2b=Tsmean95_sd[Tsmean95_sd[,2]>Tsmean2_sda[[1]],]
}


#cellStats(hotAOI,"mean")
##actual hot pixel calculation
hotAOIbin=kmeansraster==Tsmean2b[[1]]
###remove kmeansraster
hotAOINDVI=rastercon(hotAOIbin>0,NDVI2,0)
hotAOITs=rastercon(hotAOIbin>0,Ts,0)
hotAOINDVI[hotAOINDVI<=0]=NA
hotAOITs[hotAOITs<=0]=NA
hotAOITs0=quantile(hotAOITs,probs=upper)
hotAOINDVI95=quantile(hotAOINDVI,probs=lower)
hotAOINDVI95b=rastercon(hotAOINDVI<=hotAOINDVI95[[1]],hotAOINDVI,0)
hotAOINDVI95c=rastercon(hotAOITs>=hotAOITs0[[1]],hotAOINDVI95b,0)
hotAOINDVI95c[hotAOINDVI95c<=0]=NA
if(is.na(minValue(hotAOINDVI95c))&&is.na(maxValue(hotAOINDVI95c))){
hotAOINDVI95c=hotAOINDVI95b
hotAOINDVI95c[hotAOINDVI95c<=0]=NA
}
#correseponding Temperature
hotAOINDVI95cTs=rastercon(hotAOINDVI95c>=minValue(hotAOINDVI95c),hotAOITs,0)
hotAOINDVI95cTs[hotAOINDVI95cTs<=0]=NA
hotTs=hotAOINDVI95cTs
hotTs[hotTs<=0]=NA
p=data.frame(rasterToPoints(hotTs))
p3=p
if(nrow(p)>1){
p2= p[order(-p$layer),]
p3=p2[round(nrow(p2)/2),] 
}
x=p3$x
y=p3$y
ndvi=getValues(NDVI2)[cellFromXY(NDVI2,c(x,y))]

if(plot==TRUE||plot==T){
par(mfrow=c(1,2))
hist(Ts,main="Ts[Hot]",xlab="Temperature")
abline(v=p3,col="red")
hist(NDVI2,main="NDVI",xlab="NDVI")
abline(v=ndvi,col="red")
}
hotTs=p3$layer

factor<-list(x=x,y=y,xyhot=c(x,y),Tshot=hotTs,NDVI=ndvi,Ts=Ts,p=p,
             candidates=p,TH=mean(p$layer),cluster=kmeansraster)

factor$call<-match.call()

class(factor)<-"hotTs"
factor
}