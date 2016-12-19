#generic function
#default function
#' @title  Automatic computation of cold pixel Radiometric Surface Temperature
#'
#' @description This function computes temperature of a cold pixel based on NDVI, 
#' Surface Temperature, DEM and albedo. It plots histograms indicating
#' the selected values. The main difference with \code{\link{coldTs2}} is that preference 
#' is given to temperature selection before considering NDVI and albedo selection.
#' @param sunangle sun angle (degree) above the horizon.
#' @param cluster numeric, a value indicating the number of 
#' clusters.
#' @param extent A subset of the raster file than include 
#' c(xmin, xmax, ymin, ymax). You can set it to "interactive"
#' where you can digitise polygon. You can set to "auto" where it will be
#' automatically generated. You can set it to "full" where all the image dimension
#' will be used. See layout parameter below for more details.
#' see \code{\link[raster]{extent}} 
#' @param upper The quantile probability at which NDVI should be included 
#' in the selection project. It ranges from 0-1. The default is 0.95
#' @param lower The quantile probability at which Ts should be included 
#' in the slection project.The default is 0.2
#' @param plot logical:TRUE or FALSE. To indicate if histograms or triangle 
#' should be plotted
#' @param layout character. It takes "portrait" or "landscape". If extent is set to 
#' "auto" the middle part of the image will be used so that NAs can be 
#' avoided. You indicate whether the automatic slection should indicate
#' west-east (landscape) or north-south direction (portraite)
#' @param draw interactive. The method of digitisation if extent is set "auto". 
#' It takes "poly" or "rect"
#' @param iter.max The number iterations that kmeans cluster should use. 
#' @seealso \code{\link{hotTs}} 
#' @inheritParams sebal
#' @inheritParams ETohr
#' @inheritParams weather
#' @inheritParams sebkcstack
#' @details The function first divides the selected area into a  
#' number of clusters. The mean and standard devaitions are
#' computed. The cluster means with the lowest standard devaitions
#' and qualifying for upper (NDVI) and lower (Ts) limits 
#' (see above) are slected. The provided upper and lower probabilities
#'  are also used to select the final Ts.
#' @return 
#' \itemize{
#' \item{Ts:} { The selected or input Ts [K]}
#' \item{Tscold:} { The Temperature of cold pixe [K] }
#' \item{TC:} { The average temperature of the cold pixels [K]}
#' \item{x:} { The x cordinate of the cold pixel}
#' \item{y:} { The y cordinate of the cold pixel}
#' \item{xycold: } {combined x and y coordinates}
#' \item{candidates:}{ The similiar candidates' pixels that can be used}
#' }
#' @examples
#' \dontrun{
#' #using landsat folder
#' folder=system.file("extdata","stack",package="SEBKc")
#' modcold=coldTs(folder=folder,welev=170,extent="auto")
#'
#' #use object of  
#' folder=system.file("extdata","stack",package="SEBKc")
#' data=landsat578(data=folder, welev=362)
#' modcold=coldTs(folder=data,welev=170,extent="auto",cluster=3)
#' 
#' modcoldTs=coldTs(data,welev=170,extent="auto",cluster=3)
#' 
#' #using input different input data
#' albedo=raster(system.file("extdata","albedo.grd",package="SEBKc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="SEBKc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="SEBKc"))
#' #simulate with default parameters without albedo
#' modcold=coldTs(Ts=Ts,NDVI=NDVI,sunangle=50.7)
#' #simulate with albedo but set extent to "auto" and cluster to 7
#' modcold=coldTs(Ts=Ts,NDVI=NDVI,albedo=albedo,sunangle=50,cluster=7,
#' extent="auto")
#' }
#' @author George Owusu
#' @references Allen, R. G., Burnett, B., Kramber, W., Huntington, J., 
#' Kjaersgaard, J., Kilic, A., Kelly, C., & Trezza, R. 2013. Automated 
#' Calibration of the METRIC-Landsat Evapotranspiration Process. 
#' JAWRA Journal of the American Water Resources Association, 49(3): 563-576.
#' @export
#' @rdname coldTs
 
coldTs<-function(Ts,NDVI,albedo=NULL,sunangle=NULL,DEM=NULL,
                cluster=8,extent="interactive",upper=0.95,
lower=0.2,plot=TRUE,layout="portrait",draw="poly",
folder=NULL,welev=NULL,clip=NULL,
iter.max=100) UseMethod ("coldTs")
#' @export
#' @rdname coldTs
coldTs.default=function(Ts=NULL,NDVI,albedo=NULL,sunangle=NULL,DEM=NULL,
cluster=8,extent="interactive",upper=0.95,
lower=0.2,plot=TRUE,layout="portrait",draw="poly",folder=NULL,welev=NULL
,clip=NULL,
iter.max=100){
  print("Selecting cold pixel")
  
  
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
NDVI2 <- crop(NDVI, ext)
if(!is.null(albedo)){
albedo2=crop(albedo, ext)
}


upper1=quantile(Ts,probs=0.3)
lower1=quantile(Ts,probs=0.05)
#Ts[Ts>upper1[[1]],]=NA
#Ts[Ts<lower1[[1]],]=NA
Ts.kmeans <- kmeans(na.omit(Ts[]), cluster, iter.max = iter.max, nstart = 3)
kmeansraster<-raster(Ts)
kmeansraster[]<-Ts.kmeans$cluster
Tsmean=zonal(Ts, kmeansraster, 'mean') 
TsSD=zonal(Ts, kmeansraster, 'sd') 


NDVImean=zonal(NDVI2, kmeansraster, 'mean') 
NDVISD=zonal(NDVI2, kmeansraster, 'sd') 


#substract sd from the mean
NDVImean_sd=NDVImean
NDVImean_sd[,2]=NDVImean[,2]-NDVISD[,2]

#cold pixels AOI selection
NDVImax=NDVImean_sd[NDVImean_sd[,2]==max(NDVImean_sd[,2]),]
NDVImin=NDVImean_sd[NDVImean_sd[,2]==min(NDVImean_sd[,2]),]
NDVI95=quantile(NDVImean_sd[,2],  probs = upper,na.rm=TRUE)
NDVI95b=NDVImean_sd[NDVImean_sd[,2]>NDVI95,]
if(length(NDVI95b[1])==1){
Tsmean95=Tsmean[Tsmean[,1] %in%(NDVI95b[1]),]
Tsmean95SD=TsSD[TsSD[,1] %in%(NDVI95b[1]),]
Tsmean95_sd=Tsmean95
Tsmean95_sd[2]=Tsmean95[2]-Tsmean95SD[2]
Tsmean2_sda=quantile(Tsmean95_sd[2],  probs = lower)
Tsmean2b=Tsmean95_sd

}else
{
Tsmean95=Tsmean[Tsmean[,1] %in%(NDVI95b[,1]),]
Tsmean95SD=TsSD[TsSD[,1] %in%(NDVI95b[,1]),]
Tsmean95_sd=Tsmean95
Tsmean95_sd[,2]=Tsmean95[,2]-Tsmean95SD[,2]
Tsmean2_sda=quantile(Tsmean95_sd[,2],  probs = lower)
Tsmean2b=Tsmean95_sd[Tsmean95_sd[,2]<Tsmean2_sda[[1]],]
}


#cellStats(coldAOI,"mean")
##actual cold pixel calculation
coldAOIbin=kmeansraster==Tsmean2b[[1]]
###remove kmeansraster
coldAOINDVI=rastercon(coldAOIbin>0,NDVI2,0)
coldAOITs=rastercon(coldAOIbin>0,Ts,0)
coldAOINDVI[coldAOINDVI<=0]=NA
coldAOITs[coldAOITs<=0]=NA
coldAOITs0=quantile(coldAOITs,probs=lower)
coldAOINDVI95=quantile(coldAOINDVI,probs=upper)
coldAOINDVI95b=rastercon(coldAOINDVI>coldAOINDVI95[[1]],coldAOINDVI,0)
coldAOINDVI95c=rastercon(coldAOITs<coldAOITs0[[1]],coldAOINDVI95b,0)
coldAOINDVI95c[coldAOINDVI95c<=0]=NA
if(is.na(minValue(coldAOINDVI95c))&&is.na(maxValue(coldAOINDVI95c))){
coldAOINDVI95c=coldAOINDVI95b
coldAOINDVI95c[coldAOINDVI95c<=0]=NA
}
#correseponding Temperature
coldAOINDVI95cTs=rastercon(coldAOINDVI95c>=minValue(coldAOINDVI95c),coldAOITs,0)
coldAOINDVI95cTs[coldAOINDVI95cTs<=0]=NA
if(!is.null(albedo)){
abeldo_thresh=(0.001343*sunangle) +(0.3281* exp(-0.0188*sunangle))
albedo3=rastercon(Ts==coldAOINDVI95cTs,albedo2,0)
albedo3[albedo3<=0]=NA
#rm(albedo)
coldAOINDVI95cTsalbedo=rastercon(albedo3>=(abeldo_thresh-albedothre),coldAOINDVI95c,0)
coldAOINDVI95cTsalbedo2=rastercon(albedo3<=(abeldo_thresh+albedothre),coldAOINDVI95cTsalbedo,0)
coldAOINDVI95cTsalbedo2[coldAOINDVI95cTsalbedo2<=0]=NA
coldTs=rastercon(coldAOINDVI95cTsalbedo2>=minValue(coldAOINDVI95cTsalbedo2),coldAOINDVI95cTs,0)
coldTs[coldTs<=0]=NA
if(is.na(minValue(coldTs))&&is.na(maxValue(coldTs))){
coldTs=rastercon(coldAOINDVI95cTsalbedo>0,coldAOINDVI95cTs,0)
coldTs[coldTs<=0]=NA
}
if(is.na(minValue(coldTs))&&is.na(maxValue(coldTs))){
#print("yes")
coldTs=coldAOINDVI95cTs
}
coldTs[coldTs<=0]=NA
}else{
coldTs=coldAOINDVI95cTs
coldTs[coldTs<=0]=NA
}
p=data.frame(rasterToPoints(coldTs))
coldTsxy<-p[p$layer>0,]
#coldTsxy=p[,c("x","y")]
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
hist(NDVI2,main="NDVI",xlab="NDVI")
abline(v=ndvi,col="red")
}
coldTs=p3$layer

factor<-list(x=x,y=y,xycold=c(x,y),Ts=Ts,NDVI=ndvi,Tscold=coldTs,
             p=p,candidates=p,TC=mean(p$layer),cluster=kmeansraster)

factor$call<-match.call()

class(factor)<-"coldTs"
factor
}

