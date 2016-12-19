#generic function
#default function
#' @title  R computation of Surface Energy Balance Algorithm for Land 
#' and METRIC
#' 
#' @description This R function computes Bastiaanssen (1998) surface energy 
#' balance components of evapotranspiration (ET), Sensible heat(H), 
#' soil heat flux(G) and Net Radiation(Rn).
#'
#' @author George Owusu
#' @inheritParams ETohr
#' @inheritParams ETo
#' @inheritParams sebkcstack
#' @inheritParams weather
#' 
#' @param albedo A RasterLayer data that has albedo values. 
#' You can also provide file path or location on your computer
#' @param Ts A RasterLayer data that indicates radiometric 
#' surface temperature values from remote sensing image, 
#' preferably in Kelvin (K). You can also point to raster file on your 
#' computer.
#' @param NDVI A RasterLayer data that indicates NDVI 
#' (Normalized Difference Vegetation Index) values. You can also point to a 
#' file on your computer.
#' @param SAVI A RasterLayer data that indicates 
#' SAVI(Soil-adjusted Vegetation Index) values. 
#'  You can also point to file on your computer.
#' @param t1 numeric. The length of the calculation period in hour; 1 for hour, 0.5 for 30 minutes 0.25 for 15 minutes
#' @param xycold numeric or "auto". A list of x and y coordinates of a 
#' cold pixel in the form of c(x,y). If it set to "auto", \code{\link{coldTs}} 
#' will be used to compute it.
#' @param xyhot numeric or "auto". A list of x and y coordinates of a 
#' hot pixel in the form c(x,y). If it is set to "auto", \code{\link{hotTs}} 
#' will be used to compute it.
#' @param DOY Numeric or Date [YYYY-mm-dd]. Day of the Year. 
#' If you give data in the form of date [YYYY-mm-dd], it will be converted to DOY
#' @param sunelev Angle of Sun elevation in degrees, 
#' as found in the meta data of the satellite image
#' @param welev Weather station elevation in meters
#' @param zx The height above the weather station where the 
#' wind speed is measured [m]
#' @param u The satellites  overpass time wind speed at the weather station [m/s]
#' @param zomw the roughness length for the weather station surface [m]
#' @param zom the momentum roughness length for each pixel [m]. 
#' You can also point to a raster file on your computer
#' @param LAI Leaf Area Index, dimensionless
#' @param DEM A digital elevation model[m]
#' @param lapse A local lapse rate that is applied to correct DEM [K/m]. 
#' Default value is 0.0065
#' @param Rn24 A 24 hour net radiation [W/m2]. It is needed to estimate daily ET. 
#' It can be computed with \code{\link{ETo}} with even a minimal data of 
#' Tmax and Tmin. 
#'  You can also point to raster file on your computer
#' @param ETr an alfalfa or grass reference ET at the time of satellites  overpass. 
#' It only needed if you set model to "metric" [mm/hr]
#' @param ETr24 daily grass or alfalfa reference Evapotranspiration (mm/day)
#' @param model character. The type of model. It takes either "SEBAL" or "METRIC"
#' @param iter.max maximum iterations of sensible heat calculation.
#' @param folder An original  directory of the images that contains all 
#' landsat 5,7 or 8 bands and metadata.
#' At the moment only Landsat folder is supported 
#' @return 
#' \itemize{
#' \item{LE:} {  latent heat flux [W/m2]}
#' \item{Rn:} { Net Radiation [W/m2]}
#' \item{H: } { Sensible Heat Flux [W/m2]}
#' \item{G:} { Soil heat flux [W/m2]}
#' \item{ETins:} { instantaneous ET if model = "METRIC" [mm/hr]}
#' \item{ET24:} { 24 hour Evapotranspiration only if Rn24 is not NULL [mm/day]}
#' \item{ETrF:} { METRIC ET fraction [dimensionless]}
#' \item{EF:} {  ET fraction [dimensionless]}
#' \item{Ta:} { Air temperature = Ts-dT [K]}
#' \item{iter.coef:} { Parameters from sensible heat iteration}
#' } 
#' @details 
#'  \describe{
#'  \item{\strong{Computing SEBAL or METRIC}}{\code{sebal} can compute either SEBAL(Bastiaanssen et al,1998)
#'  or METRIC (Allen et al 2007). In the case of METRIC, input data such as DEM,
#'  local lapse rate, and ETr must be provided. See Allen et al (2007) for more
#'  details.}
#'  \item{\strong{Using Landsat Folder as input data}}{Though this sebkc 
#'  function is not sensor dependent, one can point the
#'  folder parameter to a Landsat TM(5), ETM(7) or 8 folder that contains 
#'  landsat bands and metadata. One can also point the folder parameter 
#'  to object of  \code{\link{landsat578}} }
#'  \item{\strong{Automatic computation}}{When xycold and xyhot are set to "auto", about 80 per cent of the input
#'  image pixels are used to compute cold and hot pixels. This avoids using
#'  NAs in classification of pixels. }
#'  \item{\strong{Taking long time to Simulate?}}{If you point the folder to landsat files and 
#'  set xycold and xyhot to "auto", it will take about 2 hours for 
#'  the function to run. Without these "autos", it will take about 45
#'  minutes. It is advisable to perform 3 step analyses of using three 
#'  functions as follows: \code{\link{sebkcstack}}, \code{\link{landsat578}},
#'  \code{\link{hotTs}}, \code{\link{coldTs}}. You need not to go through this if you have your 
#'  albedo, Ts, NDVI, SAVI from MODIS, ASTER etc}
#'  }
#'  
#' @references 
#' \itemize{
#' \item{}{Bastiaanssen, W. G. M., Menenti, M., Feddes, R. A., & Holtslag,
#'  A. A. M. 1998. A remote sensing surface energy balance algorithm for  
#'  land (SEBAL), part 1: formulation. Journal of 
#'  Hydrology: 212-213: 198-212.}
#' \item{}{Bastiaanssen, W. G. M., Pelgrum, H., Wang, J., Ma, Y., Moreno,
#'  J. F., Roerink, G. J., & van der Wal, T. 1998. A remote sensing 
#'  surface energy balance algorithm for land (SEBAL).: Part 2: 
#'  Validation. Journal of Hydrology, 212-213(0): 213-229.}
#'  \item{}{Allen, R., Tasumi, M., & Trezza, R. 2007. Satellite-Based 
#' Energy Balance for Mapping Evapotranspiration with Internalized 
#' Calibration (METRIC)-Model. Journal of Irrigation and Drainage 
#' Engineering, 133(4): 380-394.}
#' \item{}{Allen, R., Tasumi, M., Morse, A., Trezza, R., Wright, J., 
#' Bastiaanssen, W., Kramber, W., Lorite, I., & Robison, C. 2007. 
#' Satellite-Based Energy Balance for Mapping Evapotranspiration 
#' with Internalized Calibration (METRIC)-Applications. 
#' Journal of Irrigation and Drainage Engineering, 133(4): 395-406.}
#' }
#' @import raster
#' @examples
#' \dontrun{
#' #use original landsat 8 data by specifying folder path
#' folder=system.file("extdata","stack",package="SEBKc")
#' modauto=sebal(folder = folder,welev = 380,xycold="full",xyhot="full")
#' #plot ET fraction
#' plot(modauto$EF)
#' 
#' #Semi-auto
#' #load individual input files, for example from MODIS, ASTER etc
#' albedo=raster(system.file("extdata","albedo.grd",package="SEBKc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="SEBKc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="SEBKc"))
#' LAI=raster(system.file("extdata","LAI.grd",package="SEBKc"))
#' mod=sebal(albedo=albedo,Ts=Ts,NDVI=NDVI,SAVI=NULL,
#' iter.max=7,xyhot="full",xycold="full",
#' DOY=37,sunelev=50.71154048,welev=317.1,zx=10,
#' u=2,zomw=3,zom=NULL,LAI=LAI,DEM=NULL,
#' lapse=0.0065,Rn24=NULL,ETr=NULL,model="SEBAL")
#' #get latent heat flux
#' mod$LE
#' #evaporation fraction
#' mod$EF
#' 
#' #3 step simulation using landsat 7 with SLC OFF
#' ##step one: prepare the input data 
#' ### point to landsat 7 data with strips 
#' folder=system.file("extdata","slc",package="SEBKc")
#' #####get the data stacked and remove strips
#' rawdata=sebkcstack(folder,remove_cloud = 2,gap_fill = TRUE)
#' ###### check if strips are removed
#' plot(rawdata$data)
#' #get Ts, Albedo, NDVI, SAVI, sunelev, DOY from landsat data 
#' welev=278
#'  data=landsat578(rawdata,welev=welev)
#'  #perform semi-auto simulation mod=sebal(data,welev=welev)
#' # Determine xyhot. Digitize polygon on the Ts map
#' modhot=hotTs(data,welev=welev,extent="digitize",cluster=2)
#' #determine the cold. Digitize polygon on the Ts map
#' modcold=coldTs(data,welev=welev,extent="digitize",cluster=2)
#' xyhot=modhot$xyhot
#' xycold=modcold$xycold
#' modsebal=sebal(folder=data,xyhot=xyhot,xycold=xycold,welev=welev)
#' }
#' @export 
#' @rdname sebal

sebal<-function(albedo,Ts,NDVI,SAVI,welev,xyhot="auto",xycold="auto",DOY=NULL,
sunelev,zx=10,u=2,zomw=2,zom=NULL,LAI=NULL,DEM=NULL,lapse=0.0065,Rn24=NULL,ETr=NULL,
ETr24=NULL,wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",
latitude=NULL,t1=1,time=NULL, Lz=NULL,Lm=NULL,
model="SEBAL",iter.max=7,clip=NULL,folder=NULL)UseMethod ("sebal")

#' @export 
#' @rdname sebal
sebal.default<-function(albedo=NULL,Ts=NULL,NDVI=NULL,SAVI=NULL,welev,xyhot="auto",
xycold="auto",DOY=NULL,sunelev=NULL,zx=10,u=2,zomw=2,zom=NULL,LAI=NULL,DEM=NULL,
lapse=0.0065,Rn24=NULL,ETr=NULL,ETr24=NULL,wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",
latitude=NULL,t1=1,time=NULL, Lz=NULL,Lm=NULL,model="SEBAL",iter.max=7,clip=NULL,folder=NULL)
{
  
  model=toupper(model)
  #if(!is.null(ETr)){
 #   model="METRIC"
 # }
  
  eNB=NULL
  if(is.null(folder)){
    folder="/nothing454782ghf7poi8r.hope"
  }
  
  if(class(folder)=="landsat578"||class(albedo)=="landsat578"){
    file.info2=TRUE
  }else{
    file.info2=file.info(folder)[["isdir"]]
  }
  if(file.info2==TRUE&&!is.na(file.info2)){
    if(class(albedo)!="landsat578"){
    if(is.null(welev)){
        return(print("Please provide the parameter welev: 
                     the elavation of th weather station"))  
    }
    }
    if(class(folder)=="landsat578"||class(albedo)=="landsat578"){
      if(class(albedo)=="landsat578"){
        mod=albedo
        welev=mod$welev
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
      eNB=mod$eNB
      if(is.null(LAI)){
      LAI=mod$LAI
      }
    }
  
  if(is.null(ETr)||is.null(ETr24)||is.null(Rn24)){
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
    Rn24=ETr24$Rn
    ETr24= ETr24$ETo
    
    thistime=sebkc.tryCatch(mod$time)$value
    if(class(thistime)[1]=="simpleError"){
      thistime=time 
    }
    #print(ETr24)
    if(is.null(ETr)){
    ETr=sebkc.tryCatch(ETohr(DOY=thisDOY,airport = airport,wmo=wmo,latitude=latitude,
    z=zx,Krs =Krs,altitude=welev,time= thistime,Lz=Lz,t1=t1,Lm=Lm))$value
    if(class(ETr)[1]=="simpleError"){
      return (print(paste("Invalid DOY [YYYY-mm-dd] or airport or wmo or time or Lz or Lm or t1",ETr))) 
    }
    #return(print(ETr$ETo))
    ETr=ETr$ETo
    }
  }
  }
 
  
  if(!is.numeric(DOY)&&!is.null(DOY)){
    DOY2=DOY
    DOY=strptime(DOY,"%d-%m-%Y")$yday+1
    if(is.na(DOY)){
      DOY=strptime(DOY2,"%Y-%m-%d")$yday+1
    }
  }
  if(class(NDVI)!="RasterLayer"){
    NDVI=raster(NDVI)
  }
  if(class(Ts)!="RasterLayer"){
    Ts=raster(Ts)
  }
  if(class(albedo)!="RasterLayer"){
    albedo=raster(albedo)
  }
  
  if(class(DEM)!="RasterLayer"&&!is.null(DEM)){
    DEM=raster(DEM)
  }
  if(class(SAVI)!="RasterLayer"&&!is.null(SAVI)){
    SAVI=raster(SAVI)
  }
  
  if(class(zom)!="RasterLayer"&&!is.null(zom)){
    zom=raster(zom)
  }
  
  if(!is.numeric(Rn24)){
  if(class(Rn24)!="RasterLayer"&&!is.null(Rn24)){
    Rn24=raster(Rn24)
  }
  }
  
  if(class(xyhot)=="hotTs"){
    xyhot=xyhot$xyhot
  }
  
  if(class(xycold)=="coldTs"){
    xycold=xycold$xycold
  }
  testcold=NULL
  
  if((xycold=="auto"||xycold=="Auto"||xycold=="AUTO"||xycold=="full"||xycold=="Full"||xycold=="FULL")&&!is.numeric(xycold)){
    xcold=tolower(xycold)
    ycold=tolower(xycold)
    testcold=TRUE
  }
  else{
    if(length(xycold)>2){
      testcold=TRUE 
    }else{
  xcold=xycold[1]
  ycold=xycold[2]
    }
  }
  testhot=NULL
  if((xyhot=="auto"||xyhot=="Auto"||xyhot=="AUTO"||xyhot=="FULL"||xyhot=="full"||xyhot=="Full")&&(!is.numeric(xyhot))){
  xhot=tolower(xyhot)
  yhot=tolower(xyhot)
  testhot=TRUE
  }else{
    if(length(xyhot)>2){
      testhot=TRUE 
    }else{
  xhot=xyhot[1]
  yhot=xyhot[2] 
    }
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
  dr=1+(0.033*cos(DOY*((2*pi)/365)))
  cos_teta=cos(((90-sunelev)/360)*2*pi)
  tsw=(0.75+2*0.00001*welev)^2
  ea=0.85*(-log(tsw))^0.09
  stefan=5.67E-08
  print("computing incoming shortwave solar radiation")
  Rs_in=1367*cos_teta*dr*tsw
  sunangle=sunelev
  if(is.null(DEM)){
    pressure=101.325
  }else{
    print("Computing air pressure with DEM")
    pressure=101.325*((293-(0.0065*DEM))/293)^5.26
    Ts=Ts+(lapse*DEM)
  }
  if(is.null(LAI)){
    LAI=(5.434*SAVI)+1.5416
    #LAI=-(log(0.69-SAVI)/0.59)/(0.91)
  }
  
  if(!is.null(testhot)){
    #print(xyhot)
    modhot=hotTs(Ts=Ts,NDVI=NDVI,albedo=albedo,DEM=DEM,
                 extent=xyhot,upper=0.80,lower=0.1)
      xhot=modhot$x
      yhot=modhot$y
    
  }
  if(!is.null(testcold)){
    modcold=coldTs(Ts=Ts,NDVI=NDVI,albedo=albedo,DEM=DEM,
                   extent=xycold,upper=0.95,lower=0.1,sunangle = sunangle)
    
      xcold=modcold$x
      ycold=modcold$y
    
    
  }
  if(is.null(eNB)){
    eNB=rastercon( NDVI<0 & albedo<0.47,0.99,rastercon(LAI>=3,0.98,0.97+(LAI*0.0033)))   
  }
 
  e0=rastercon( NDVI<0 & albedo<0.47,0.985,rastercon(LAI>=3,0.98,0.95+(LAI*0.01)))
  print("Computing Vaporaisation with Surface Temperature")
  vaporisation=(2.501-(0.00236*(Ts-273.15)))*10^6
  print("Computing outgoing longwave radiation")
  RL_out=e0*stefan*(Ts^4)
  Ts_cold=getValues(Ts)[cellFromXY(Ts,c(xcold,ycold))]
  print("Computing incoming longwave radiation")
  RL_in=ea*stefan*(Ts_cold^4)
  #Net radiation
  Rn=((1-albedo)*Rs_in)+RL_in-RL_out-((1-e0)*RL_in)
  #Soil Heat flux
  print("Computing Soil heat flux")
  G_by_Rn=(Ts-273.15)*(0.0038+0.007*albedo)*(1-0.98*NDVI^4)
  G_by_Rn[NDVI<0,]=0.5
  G_by_Rn[Ts-273.15<4&albedo>0.45,]=0.5
  G=G_by_Rn*Rn
  
  print("Computing surface roughness[zom]")
  #sensible Heat
  if(is.null(zom)){
    zom=0.018*LAI #surface roughness
  }
  u200=((uw*log(200))/zomw)/log(zx/zomw)
  u_star=(k*u200)/(log(200/zom))# friction velocity
  rah=log(z2/z1)/(u_star*k)#heat tranpsort
  df = data.frame(a= numeric(), b = numeric(), u_star = numeric(),
                  L = numeric(),pairhot=numeric(),paircold=numeric(),
                  rah_hot=numeric(),
                  rah_cold=numeric(),dThot=numeric(),dTcold=numeric(),
                  Hhot=numeric(),Hcold=numeric(),Ts_hot=numeric(),
                  Ts_cold=numeric(),
                  Rn_hot=numeric(),Rn_cold=numeric(),G_hot=numeric(),
                  G_cold=numeric(),
                  LEcold=numeric(),LEhot=numeric(),ETint_hot=numeric(),
                  ETint_cold=numeric())
  #df
  Ts_hot=getValues(Ts)[cellFromXY(Ts,c(xhot,yhot))]
  Rn_cold=getValues(Rn)[cellFromXY(Rn,c(xcold,ycold))]
  Rn_hot=getValues(Rn)[cellFromXY(Rn,c(xhot,yhot))]
  G_cold=getValues(G)[cellFromXY(G,c(xcold,ycold))]
  G_hot=getValues(G)[cellFromXY(G,c(xhot,yhot))]
  print("Computing sensible heat flux........................")
  
  i=1
  rah_hot_change=-10000
  coldDT=0
  while (i<=iter.max)
  {
    print(paste("Monin-Obukhov length iteration",i,"of",iter.max))
    
    #change in temerature at the hot pixel
    rah_hot=getValues(rah)[cellFromXY(rah,c(xhot,yhot))]
    rah_cold=getValues(rah)[cellFromXY(rah,c(xcold,ycold))]
    #print(paste("rah_cold",rah_cold))
    
    if(i>1){
      rah_hot_change=round((((rah_hot-rah2)/rah2)*100),2)
    }
    #rah_hot change
    hotDT=((Rn_hot-G_hot)*rah_hot)/(p*cp)
    
    #slope of the line
    a=(hotDT-coldDT)/(Ts_hot-Ts_cold)
    
    #b=hotDT-(a*Ts_hot)
    b=(hotDT-a)/Ts_hot
    #b=-a*Ts_cold
    b=hotDT-(a*Ts_hot)
    
    #the actual DT
    dT=(a*Ts)+b
    pair=(1000*pressure)/(1.01*((Ts-dT)*287))
    p=getValues(pair)[cellFromXY(pair,c(xhot,yhot))]
    
    #sensible heat
    H=(pair*cp*dT)/rah
    #LE=Rn-G-H
    H_hot=getValues(H)[cellFromXY(H,c(xhot,yhot))]
    H_cold=getValues(H)[cellFromXY(H,c(xcold,ycold))]
    
    if(model=="METRIC"){
      if(is.null(ETr)){
        return(print("ETr is needed"))
      }
      LE_cold=1.05*ETr*getValues(pair)[cellFromXY(pair,c(xcold,ycold))]
      H_cold=(Rn_cold-LE_cold-G_cold)/3600
      coldDT=(H_cold*rah_cold)/(getValues(pair)[cellFromXY(pair,
                        c(xcold,ycold))]*cp)
    }else{
      coldDT=(H_cold*rah_cold)/(getValues(pair)[cellFromXY(pair,
                                          c(xcold,ycold))]*cp)
    }
    u_star_hot=getValues(u_star)[cellFromXY(u_star,c(xhot,yhot))]
    
    L=-((p*cp*(u_star_hot^3)*Ts_hot)/(k*g*H_hot))
    x_200m=(1-(16*(200/L)))^0.25
    x_2m=(1-(16*(2/L)))^0.25
    x_01m=(1-(16*(0.1/L)))^0.25
    
    w_200m=ifelse(L<0,2*log((1+x_200m)/2)+log((1+x_200m^2)/2)-2*
                    atan(x_200m)+0.5*pi,(ifelse(L>0,-5*(200/L),0)))
    w_2m=ifelse(L<0,2*log((1+x_2m^2)/2),ifelse(L>0,-5*(2/L),0))
    w_01m=ifelse(L<0,2*log((1+x_01m^2)/2),ifelse(L>0,-5*(2/L),0))
    
    u_star=(u200*k)/(log(200/zom)-w_200m)
    rah=(log(z2/z1)-w_2m+w_01m)/(u_star*k)
    rah2=getValues(rah)[cellFromXY(rah,c(xhot,yhot))]
    
    df = rbind(df,data.frame(a= a, b = b, u_star = u_star_hot,L = L,
                             paircold=getValues(pair)[cellFromXY(pair,
                                  c(xcold,ycold))],pairhot=p,
                             rah_hot=rah_hot,rah_cold=rah_cold,
                             dThot=hotDT,dTcold=coldDT,
                             Hhot=H_hot,Hcold=H_cold,Ts_hot=Ts_hot,
                             Ts_cold=Ts_cold,
                             Rn_hot=Rn_hot,Rn_cold=Rn_cold,G_hot=G_hot,
                             G_cold=G_cold,
                             LEcold=Rn_cold-G_cold-H_cold,
                             LEhot=Rn_hot-G_hot-H_hot,ETint_hot=t1*3600*
                          ((Rn_hot-G_hot-H_hot)/(getValues(vaporisation)
                            [cellFromXY(vaporisation,c(xhot,yhot))])),
                             ETint_cold=t1*3600*((Rn_cold-G_cold-H_cold)/
                            (getValues(vaporisation)[cellFromXY
                              (vaporisation,c(xcold,ycold))]))))
    
    i=i+1
  }
  #H=Con(H<0,0,H)
  print ("computing latent heat flux" )
  LE=Rn-G-H
  #LE=Con(LE<0,0,LE)
  ETrF=NULL
  ET24=NULL
  EF=NULL
  
  ETint=(t1*3600*LE)/vaporisation
  if(model=="METRIC"){
    EF=ETint/ETr
    EF[EF>1.1,]=1.1
    EF[EF<0,]=0
    
    if(!is.null(ETr24)){
     ET24=EF*ETr24 
    }
    
  }else{
    EF=1.1*(LE/(Rn-G))
    EF[EF<0,]=0
    EF[EF>1.1,]=1.1
    if(!is.null(Rn24)){   
      ET24=EF*Rn24*0.035
     }
  }

  factor=list(rah=rah, Rn=Rn, H=H, G=G,LE=LE,ETins=ETint,ETr=ETr,
       rah=rah,iter.coef=df,vaporisation=vaporisation,
       pair=pair,dT=dT,Ta=Ts-dT,wm=w_200m,EF=EF,ETrF=EF,ET24=ET24,
       folder=folder,DOY=DOY,G_by_Rn=G_by_Rn,Rs=Rs_in)
  
  factor$call<-match.call()
  
  class(factor)<-model
  factor
}

