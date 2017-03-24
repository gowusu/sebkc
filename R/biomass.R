
#' @title Inverse distance interpolation
#' 
#' @description This function performs inverse distance interpolation based 
#' on projection and geometry of input spatial data 'ext'. 
#' @inheritParams ETo
#'
#' @param var numeric variable. The length should be same as the length of 
#' latitude and longitude
#' @param ext RasterLayer object. 
#' @param idp Neighbourhood parameter
#'
#' @return var: interpolated values
#' @author George Owusu
#' @export
#'
#' @examples
#'  \dontrun{
#' #Get ext data
#' folder=system.file("extdata","stack",package="sebkc")
#' data=landsat578(data=folder, welev=362)
#' temp=data$Ts
#' latitude=seq(7.544,7.590,0.001)
#' longitude=seq(-1.211,-1.187,0.001)
#' var=seq(20,30,0.2)
#' latitude=latitude[1:length(longitude)]
#' var=var[1:length(longitude)]
#' map=invdist(longitude,latitude,var,ext=temp)
#' }
invdist=function(longitude,latitude,var,ext,idp = .5){
  
 
#if(!is.null(longitude)){
  d <- data.frame(longitude=longitude, latitude=latitude,var=var)
 names(d)=c("longitude","latitude","var")
  coordinates(d) <- c("longitude", "latitude")
  lon1=ifelse(length(longitude[[1]])>1,longitude[[1]][1],longitude[[1]])
  lat1=ifelse(length(latitude[[1]])>1,latitude[[1]][1],latitude[[1]])
  
  if(lon1<300||lat1<300){
  projection(d) <- CRS("+init=epsg:4326") # WGS 84
  CRS.new=CRS( proj4string(ext))
  d <- spTransform(d, CRS.new)
  #data=cbind(dtrans@coords,dtrans@data)
  #print(coordinates(data))
  }
    projection(d) <- proj4string(ext)
    #data=d
    mg <- gstat(id = "var", formula = var~1, data=d, set=list(idp = idp))
    
    var <- interpolate(ext, mg)
  
  ## inverse distance weighted (IDW)
  var
  #list(var=var)
}

#' Reproject WGS 84 longitude and latitude  into a new raster coordinates and vice versa 
#'
#' @details If x is a list of longitudes then new  horizontal  coordinates based on map
#' coordinates. If x is a raster object with  horizontal  coordinates then it will be projected to
#' WGS 84.
#' @param x a longitude or raster or vector object
#' @param y a latitude or raster or vector object
#' @param map Raster* object that will be used to define the new coordinates.
#' @param epsg coordinate reference number. The default is WGS 84: 4326
#' @inheritParams ETo
#' @inheritParams invdist
#' @return New x, y and map
#' @export 
#' @author George Owusu
#' @examples
#' folder=system.file("extdata","stack",package="sebkc")
#' data=landsat578(data=folder, welev=362)
#' temp=data$Ts
#' latitude=seq(7.544,7.590,0.001)
#' longitude=seq(-1.211,-1.187,0.001)
#' var=seq(20,30,0.2)
#' latitude=latitude[1:length(longitude)]
#' var=var[1:length(longitude)]
#' newlonlat=lonlatreproject(longitude,latitude,var=var,ext=temp)
#' longitude=newlonlat$longitude
#' map=newlonlat$map
#' 
#' 
lonlatreproject=function(x,y=NULL,var=NA,map,epsg=4326){
  longitude=x
  latitude=y
  #ext=map
  if(is.null(var)){
    var=NA
  }
  if(is.null(latitude)){
    d=longitude
  }else{
  d <- data.frame(longitude=longitude, latitude=latitude,var=var)
  coordinates(d) <- c("longitude", "latitude")
  }
  if(extent(d)[1]>180){
    data= projectRaster(d, crs=paste("+init=epsg:",epsg,sep=""))
  }else{
  projection(d) <- CRS(paste("+init=epsg:",epsg,sep="")) # WGS 84
  CRS.new=CRS( proj4string(map))
  data <- spTransform(d, CRS.new)
  }
  #data=cbind(dtrans@coords,dtrans@data)
  
  list(x=coordinates(data)[,1],y=coordinates(data)[,2],map=data)
}

#' @title Biomass and yield estimation of crops and plants
#' 
#' @description This function uses	Wageningen procedure 
#' (De Wit, 1968; Slabbers et al., 1979) to compute the 
#' point or spatial possible maximum yield of a standard crop based on a 
#' Rs data. This function can also perform internal interpolation of 
#' data if longitude and latitude are provided. See datails and examples below.
#'
#' @param Tday numeric. Average Day time temperature in degrees Celsius 
#' for the growing period. A Remote sensing surface Temperature can also be used
#' 
#' @param T24 numeric. Average 24-hour (including night hours) temperature 
#' in degrees Celsius for the growing period. 
#' @param Rs A shortwave radiation data [ca^-2 day^-1]. 
#' RasterLayer data can be supplied.
#' @inheritParams ETo
#' @inheritParams sebal
#' @param HI numeric. Harvest Index. It can be list or spatial
#' @param plantdate character. A date of planting in the form "YYYY-mm-dd"
#' @param harvestdate character. A date of harvest in the form "YYYY-mm-dd"
#' @param c30 numeric. A maintenance respiration coefficient (c30) at 30 oC. 
#' it takes either a value of 0.0108 for non-legume crop or 0.0283 for legumes
#' @param adaptability numeric. The  climate adaptability group  number. 
#' It takes values of 1,2,3,or 4. it depends on the temperature profile of 
#' study area. See references (eg. FAO (1981)) for more explanation. 
#' @param crop character. A crop type. For example "maize"
#' @author George Owusu
#' @details 
#'  \describe{
#'  The function can be used in several ways:
#'  \item{\strong{Point Estimation}}{\code{biomass} can compute biomass 
#'  and yield for  one location. Each input data must be given only one value}
#'  \item{\strong{Point Estimation for several locations }}{The function can 
#'  compute biomass and yields for more than one location. Input data must 
#'  be a list or array.}
#'  \item{\strong{Spatial Estimation}}{Input data such as Rs, Tday, T24,
#'  HI or LAI can be spatial.}
#'  \item{\strong{Spatial Estimation with internal interpolation}}{It takes Tday,
#'   for example from satellite image, and do interpolation for 
#'  the rest of the data. Both longitude, latitude, and the variables must be of
#'  the same length. If longitude is not provided no interpolation will be done.}
#'  }
#'  @references 
#' \itemize{
#' \item{}{De Wit, C. T. (1968). Plant production Miscellaneous Papers: 
#' Landbouw Hogeschool, Wageningen.}
#' \item{}{Slabbers, P. J., Herrendorf, V. S., & Stapper, M. (1979). 
#' Evaluation of simplified water-crop yield models. 
#' Agricultural Water Management, 2(2), 95-129. 
#' doi: http://dx.doi.org/10.1016/0378-3774(79)90026-X}
#' }
#' \item{}{FAO (1981). Report on the Agro-ecological Zones Project: 
#' Food and Agriculture Organization of the United Nations.}
#' @return yield and biomass in t/ha
#' @export
#'
#' @examples
#' \dontrun{
#' ################# Point Estimation ###################
#' # Identify the climatic input data on 
#' latitude=11.18 #in decimal degrees
#' Rs= 452 #Average shortwave radiation over the growing period (cal-2day-1)
#' Tday = 26.9 #	Average day time difference (oC) over the growing period
#' T24 =25.3 #Average 24 hour mean temperature over the growing period
#' # Gather crop Information
#' crop="maize"
#' plantdate="2015-05-1" #Start date
#' harvestdate="2015-08-29" #end date 120 days
#' HI=0.4 #Harvest index for maize
#' adaptability=3 #crop climate adaptability group III
#' LAI = 4 #Leaf Area Index
#' pointmod=biomass(Tday,T24,Rs,latitude,longitude=NULL,HI,
#' plantdate,harvestdate,c30=0.0108,adaptability=3,LAI=5,crop=NULL)
#' 
#' ################# Point Estimation for several locations#################
#' latitude=c(10,11.18,12,9)
#' longitude=c(-1.4,-2,-1,-2.6)
#' Rs= c(450,452,400,500)
#' Tday = c(25,26.9,24,20,30)
#' T24=c(24,25.3,23,28)
#' HI=c(0.4,0.4,0.5,0.25)
#' plantdate=c("1-03-2015","1-04-2015","1-05-2015")
#' harvestdate=c("29-08-2015","29-08-2015","29-08-2015")
#' 
#' # computing several points
#' pointSmod=biomass(Tday,T24,Rs,latitude,longitude=longitude,HI,
#' plantdate,harvestdate,c30=0.0108,adaptability=3,LAI=5,crop=NULL)
#' yield=pointSmod$yield #tabular data
#' biomass=pointSmod$output #tabular data
#' #access on biomass
#' pointSmod$output$yield
#' 
#' ############## Spatial Estimation with internal interpolation ###########################
#' #generate spatial data
#' folder=system.file("extdata","stack",package="sebkc")
#' modauto=landsat578(data=folder, welev=362)
#' Tday=modauto$Ts-273.15 #temperature in degree celcius
#' LAI=modauto$LAI
#' latitude=seq(7.544,7.590,0.001) #for interpolation
#' longitude=seq(-1.211,-1.187,0.001) #for interpolation
#' latitude=latitude[1:length(longitude)] #to make sure they are of the same length
#' spatialmodint=biomass(Tday,T24,Rs,latitude,longitude=longitude,HI,
#' plantdate,harvestdate,c30=0.0108,adaptability=3,LAI=LAI,crop=NULL)
#' # Rs, T24, and HI are interpolated
#' plot(spatialmodint$yield)
#' ############## Spatial Estimation  ###########################
#' Rs= 452
#' T24 =25.3
#' HI=0.4 
#' spatialmod=biomass(Tday,T24,Rs,latitude,longitude=longitude,HI,
#' plantdate,harvestdate,c30=0.0108,adaptability=3,
#' LAI=LAI,crop=NULL)
#' plot(spatialmod$yield)
#' }

biomass=function(Tday=NULL,T24,Rs,latitude,longitude=NULL,HI=NULL,
plantdate,harvestdate,c30=0.0108,adaptability=3,LAI=5,crop=NULL){
  if(class(Tday)=="kc"){
    Tday=Tday$DL
    #print(Tday)
  if(!is.null(Tday$Tday)){
    data1=Tday
    T24=data1$T24
    Rs=data1$Rs
    latitude=data1$latitude
    plantdate=data1$plantdate
    harvestdate=data1$harvestdate
    Tday=data1$Tday
    
  }
  }
  inter=NULL
  if(class(Tday)=="RasterLayer"){
    inter=Tday  
  }else{
    if(class(Rs)=="RasterLayer"){
      inter=Rs  
    }else{
      if(class(T24)=="RasterLayer"){
        inter=T24
      }else{
        if(class(LAI)=="RasterLayer") 
          inter=LAI
      }
    }
  }
  
  
  if(is.null(inter)&&(length(plantdate)>1||length(harvestdate)>1||length(Tday)>1
     ||length(T24)>1||length(Rs)>1||length(c30)>1||length(adaptability)>1
     ||length(LAI)>1||length(crop)>1)){
    
outputdata=data.frame(Tday=numeric(),T24=numeric(),Rs=numeric(),
                      latitude=numeric(),longitude=numeric(),HI=numeric(),
  plantdate=numeric(),harvestdate=numeric(),c30=numeric(),
  adaptability=numeric(),LAI=numeric(),yield=numeric(),biomass=numeric())

outputdata=data.frame()
                         
    i=1
    mod=NULL
    while(i<=length(plantdate)){
      
      Tday2=Tday
      if(length(Tday)>1){
        Tday2=Tday[i]
      }
      T242=T24
      if(length(T24)>1){
        T242=T24[i]
      }
      rs2=Rs
      if(length(T24)>1){
        rs2=Rs[i]
      }
      longitude2=longitude
      if(length(longitude)>1){
        longitude2=longitude[i]
      }
      
      latitude2=latitude
      if(length(latitude)>1){
        latitude2=latitude[i]
      }
      HI2=HI
      if(length(HI)>1){
        HI2=HI[i]
      }
      c302=c30
      if(length(c30)>1){
        c302=c30[i]
      }
      adaptability2=adaptability
      if(length(adaptability)>1){
        adaptability2=adaptability[i]
      }
      LAI2=LAI
      if(length(LAI)>1){
        LAI2=LAI[i]
      }
      plantdate2=plantdate
      if(length(plantdate)>1){
        plantdate2=plantdate[i]
      }
      harvestdate2=harvestdate
      if(length(harvestdate)>1){
        harvestdate2=harvestdate[i]
      }
      crop2=crop
      if(length(crop)>1){
        crop2=crop[i]
      }
      pointmod=biomass(Tday=Tday2,T24=T242,Rs=rs2,latitude=latitude2,
        longitude=longitude2,HI=HI2,plantdate=plantdate2,harvestdate=harvestdate2,
        c30=c302,adaptability=adaptability2,
                        LAI=LAI2,crop=crop2)
      if(is.null(longitude)){longitude2=NA}
      outputdata=rbind(outputdata,data.frame(Tday=Tday2,T24=T242,Rs=rs2,
                                             latitude=latitude2,longitude=longitude2,HI=HI2,plantdate=plantdate2,
                                             harvestdate=harvestdate2,
              c30=c302,adaptability=adaptability2,LAI=LAI2,yield=pointmod$yield,
             biomass=pointmod$biomass))
      #print (c(latitude2,longitude2))
      #mod=rbind(mod,pointmod)
      #print (cbind(output,Tday=Tday2,T24=T242,Rs=rs2,latitude=latitude2,longitude=longitude2,HI=HI2,plantdate=plantdate2,harvestdate=harvestdate2,c30=c302,adaptability=adaptability2,LAI=LAI2,crop=crop2))
      
      
      i=i+1
    }
    
    #colnames(mod)=paste(plantdate,"to",harvestdate,sep="")
    
    return(list(output=outputdata,yield=outputdata,biomass=outputdata))
    
  }else{
  
Ndays=as.numeric(as.Date(harvestdate, "%Y/%m/%d")-as.Date(plantdate, "%Y/%m/%d"))
if(is.na(Ndays)){
  Ndays=as.numeric(as.Date(harvestdate, "%d-%m-%Y")-as.Date(plantdate, "%d-%m-%Y"))
  
}
if(is.na(Ndays)){
  Ndays=as.numeric(as.Date(harvestdate, "%d/%m/%Y")-as.Date(plantdate, "%d/%m/%Y"))
  
}
if(is.na(Ndays)){
  Ndays=as.numeric(as.Date(harvestdate, "%m/%d/%Y")-as.Date(plantdate, "%m/%d/%Y"))
}

if(is.na(Ndays)){
  Ndays=as.numeric(as.Date(harvestdate, "%m-%d-%Y")-as.Date(plantdate, "%m-%d-%Y"))
}
if(is.na(Ndays)){
  Ndays=as.numeric(as.Date(harvestdate, "%Y-%m-%d")-as.Date(plantdate, "%Y-%m-%d"))
}
#print(Ndays)
Pm20=20
inter=NULL


if(latitude<0){
  Pmfile=read.table(system.file("extdata","sys","pmwrite2.txt",package="sebkc"),header=TRUE)
  latrad=read.table(system.file("extdata","sys","latradiancesouth.txt",package="sebkc"),header=TRUE)
  boin=read.table(system.file("extdata","sys","bowrite2south.txt",package="sebkc"),header=TRUE)
  bcin=read.table(system.file("extdata","sys","bcwrite2south.txt",package="sebkc"),header=TRUE)  
}else{Pmfile=read.table(system.file("extdata","sys","pmwrite2.txt",package="sebkc"),header=TRUE)
latrad=read.table(system.file("extdata","sys","latradiance.txt",package="sebkc"),header=TRUE)
boin=read.table(system.file("extdata","sys","bowrite2.txt",package="sebkc"),header=TRUE)
bcin=read.table(system.file("extdata","sys","bcwrite2.txt",package="sebkc"),header=TRUE)
}


start_date=as.numeric(format.Date(as.Date(plantdate), "%m"))+1
end_date=as.numeric(format.Date(as.Date(harvestdate), "%m"))+1
latrad=latrad[c(1,start_date:end_date)]
boin=boin[c(1,start_date:end_date)]
bcin=bcin[c(1,start_date:end_date)]

rdata=FALSE

if(class(Tday)=="RasterLayer"){
inter=Tday  
}else{
  if(class(Rs)=="RasterLayer"){
    inter=Rs  
  }else{
    if(class(T24)=="RasterLayer"){
      inter=T24
    }else{
    if(class(LAI)=="RasterLayer") 
      inter=LAI
    }
  }
}
if(class(Tday)=="RasterLayer"){
  rdata=TRUE
  
  if(is.null(longitude)&&length(latitude>1)){
  latitude=mean(latitude)  
  }
}

Pm=NULL

if(length(latitude)>1){
  latrad2 <- NULL
  boin2 <-NULL
  bcin2 <- NULL
  i=1
  while(i<=length(latitude)){
    boin2b <- boin[boin$id==round(latitude[i],1), ]
    boin2=rbind(boin2,boin2b)
    
    latrad2b <- latrad[latrad$id==round(latitude[i],1), ]
    latrad2=rbind(latrad2,latrad2b)
    
    bcin2b <- bcin[bcin$id==round(latitude[i],1), ]
    bcin2=rbind(bcin2,bcin2b)
    
    #latrad2   
 i=i+1
 }
  
  #Ac
  Ac=rowMeans((latrad2[-1]))
  #bo
  bo=rowMeans((boin2[-1]))
  #bc
  bc=rowMeans((bcin2[-1]))
  
  if(!is.null(inter)&&!is.null(longitude)){
    Ac=invdist(longitude,latitude,Ac,inter)
    #Ac=Ac$var
    bc=invdist(longitude,latitude,bc,inter)
    #bc=bc$var
    bo=invdist(longitude,latitude,bo,inter)
    #bo=bo$var
    
    if(length(latitude)==length(T24)){
      T24=invdist(longitude,latitude,T24,inter)
      #T24=T24$var
    }
    if(length(latitude)==length(LAI)){
      LAI=invdist(longitude,latitude,LAI,inter)
      #LAI=LAI$var
      #plot(LAI)
    }
    
    if(length(latitude)==length(HI)){
      HI=invdist(longitude,latitude,HI,inter)
      #HI=HI$var
    }
    
    if(length(latitude)==length(Rs)){
      Rs=invdist(longitude,latitude,Rs,inter)
      #Rs=Rs$var
    }
    
  }
  }else{
  latrad2 <- latrad[latrad$id==round(latitude,1), ]
  boin2 <- boin[boin$id==round(latitude,1), ]
  bcin2 <- bcin[bcin$id==round(latitude,1), ]
  
  #Ac
  Ac=mean(t(latrad2[-1]))
  #bo
  bo=mean(t(boin2[-1]))
  #bc
  
  bc=mean(t(bcin2[-1]))
  
  
  
  }

##########################################
if(length(Tday)>1&&rdata!=TRUE){
  Pm=NULL
  i=1
  while(i<=length(Tday)){

    if(length(adaptability)==length(Tday)){
      adaptability=adaptability[i]
    }
    
      Tday2=Tday[i]
      if(adaptability==1){
        Pm2 <- Pmfile[Pmfile$id==round(Tday2,1), ][[2]]
      }
      
      if(adaptability==2){
        Pm2 <- Pmfile[Pmfile$id==round(Tday2,1), ][[3]]
      }
      
      if(adaptability==3){
        Pm2 <- Pmfile[Pmfile$id==round(Tday2,1), ][[4]]
      }
      
      if(adaptability==4){
        Pm2 <- Pmfile[Pmfile$id==round(Tday2,1), ][[5]]
        
      }
      Pm=rbind(Pm,Pm2)
    #print(Pm)
    #latrad2   
    i=i+1
  }
}




ct=c30*(0.044+(0.0019*T24)+(0.0010*(T24^2)))

#Ac=361
#bo=427
#bc=228
growth_rate=-0.0307*LAI^2 + 0.3549*LAI - 0.0111
fraction=(Ac-(0.5*Rs))/(0.8*Ac)
bgm=(fraction*bo)+((1-fraction)*bc)
#if(LAI!=5){
  bgm=bgm* growth_rate
#}
  
  if(rdata==TRUE){
    #Pm1=Pmfile[1:2]
    #Pm2=Pmfile[c(1,3)]
    #Pm3=Pmfile[c(1,4)]
    #Pm4=Pmfile[c(1,5)]
    if(adaptability==1){
    Pm=subs(round(Tday,1),Pmfile[1:2])
    }
    
    if(adaptability==2){
      Pm=subs(round(Tday,1),Pmfile[c(1,3)])
    }
    if(adaptability==3){
      Pm=subs(round(Tday,1),Pmfile[c(1,4)])
    }
    if(adaptability==4){
      Pm=subs(round(Tday,1),Pmfile[c(1,5)])
    }
    
  }else{

    if(!is.numeric(Pm)){
      
if(adaptability==1){
  
  Pm <- Pmfile[Pmfile$id==round(Tday,1), ][[2]]
  }


if(adaptability==2){
  Pm <- Pmfile[Pmfile$id==round(Tday,1), ][[3]]
}

if(adaptability==3){
  Pm <- Pmfile[Pmfile$id==round(Tday,1), ][[4]]
}


if(adaptability==4){
  Pm <- Pmfile[Pmfile$id==round(Tday,1), ][[5]]
  
}
     # print(Pm)
}
}

  
#Pm=max(Pm,0)

y=((Pm-Pm20)/Pm20)*100

if(rdata==FALSE){
if(bgm>Pm20){
 bgm=bgm+(((y/5)/100)*(fraction*bo)+((y/2)/100)*(1-fraction)*bc) 
}else{
  bgm=bgm-(((y/2)/100)*(fraction*bo)+(y/100)*(1-fraction)*bc) 
}
  
}else{
 bgm= rastercon(bgm>Pm20,bgm+(((y/5)/100)*(fraction*bo)+((y/2)/100)*(1-fraction)*bc),
                bgm-(((y/2)/100)*(fraction*bo)+(y/100)*(1-fraction)*bc))
  }
  
  
Bn=(0.36*bgm)/(1/Ndays+0.25*ct)/907.18474

  By=Bn*HI
  

list(biomass=Bn,yield=By)
}
}

########################



   
   