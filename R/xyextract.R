#' @title Reprojecting and extracting Raster values 
#' 
#' @description The function reprojects longitudes and latitudes into map
#' coordinates and extract raster values
#'
#' @param map a raster object.
#' @inheritParams ETo
#' @return data
#' @author George Owusu
#' @export 
#'
#' @examples
#' \dontrun{
#' lonlat=read.table(system.file("extdata","sys","xyvalues.txt",package="sebkc"),header=TRUE)
#' longitude=lonlat$longitude
#' latitude=lonlat$latitude
#' folder=system.file("extdata","stack",package="sebkc")
#' data=landsat578(data=folder, welev=362)
#' Ts=data$Ts
#' Tsdata=xyextract(Ts,longitude,latitude)
#' 
#' }

xyextract=function(map,longitude,latitude){
  
  #if(!is.null(longitude)){
  
  # if(!is.null(map$xydata)){
  #  longitude=map$xydata$longitude
  #  latitude=map$xydata$latitude
  #}
  
  #if(!is.null(map$map$ETa)){
  #  map=map$map$ETa
  #}
  d <- data.frame(x=longitude, y=latitude)
  d2=d
  coordinates(d) <- c("x", "y")
  if(longitude[[1]]<300||latitude[[1]]<300){
  projection(d) <- CRS("+init=epsg:4326") # WGS 84
  CRS.new=CRS( proj4string(map))
  dtrans <- spTransform(d, CRS.new)
  }else{
    dtrans= cbind(longitude, latitude)
    #projection(dtrans) <- proj4string(map)
  }
  
  data=extract(map, dtrans,df=TRUE)
  data=cbind(d2,data)
  data$ID=NULL
  thisname=paste("value",1:(length(data)-2),sep="")
  names(data)=c("longitude","latitude",thisname)
  data
}
#xydata=read.table("C:/Users/GeoKings/Documents/George Owusu/UG PhD/lectures/remote sensing/model/sebkc/inst/extdata/sys/xyvalues.txt",header=TRUE)

#Tsdata2[5]=rev(Tsdata2[5])
#Ts=data$Ts
#Px=read.table("C:/Users/GeoKings/Documents/George Owusu/UG PhD/lectures/remote sensing/model/sebkc/inst/extdata/sys/P.txt",header=T)
#ETOI=read.table("C:/Users/GeoKings/Documents/George Owusu/UG PhD/lectures/remote sensing/model/sebkc/inst/extdata/sys/ETOI.txt",header=T)
# I=ETOI$I
#ETo=ETOI$Eto
#mapmatrix(Px,1,map)
#

#' @title Internal function for doing interpolation
#' 
#' @description Internal function for doing interpolation
#'
#' @param xydata x and y data
#' @param i index
#' @param map geographical map
#'
#' @return map
#' @author George Owusu
#' @export 
#' @examples
#' \dontrun{
#' folder=system.file("extdata","stack",package="sebkc")
#' sebiauto=sebi(folder=folder,welev=317,Tmax=31,Tmin=28)
#' points=sampleRandom(sebiauto$EF,100,sp=TRUE)
#' pt=cbind(points@@coords,points@@data)
#' longitude=pt$x
#' latitude=pt$y
#' value=pt$layer
#' mapmatrix=mapmatrix(xydata,2)
#' map=invdist(longitude,latitude,var=value,ext=sebiauto$EF)
#' }
mapmatrix=function(xydata,i,map){
  
rclmat <- matrix(xydata[i,], ncol=3, byrow=TRUE)
longitude=apply(rclmat, 1,as.numeric)[1,]
latitude=apply(rclmat, 1,as.numeric)[2,]
var=apply(rclmat, 1,as.numeric)[3,]

#rclmat=data.frame(cbind(longitude,latitude,var))
map=invdist(longitude=longitude,latitude=latitude,var=var,ext=map)
map$var
}
#mapmatrix=mapmatrix(xydata,2)
#map=invdist(rclmat,ext=Ts)


