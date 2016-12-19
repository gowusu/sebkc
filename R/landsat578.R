#' @title Radiometric correction and processing of Landsat Images
#'
#' @description \code{landsat578} performs radiometric collection of 
#' landsat 5,7,and 8 images. It returns albedo, NDVI, 
#' fractional vegetation cover, surface temperature
#' SAVI,radiance and reflectance of the image
#' @param data An original folder containing Landsat data or
#' Raster \code{\link[raster]{brick}}  or raster 
#' \code{\link[raster]{stack}}. Arrange the bands in ascending order:
#'  c(band1,band2,band3,band4,band5,band6,band7) for Landsat 5 and 7.
#'  And a similar  order for Landsat 8.  On folder see \code{\link{sebkcstack}}
#'  See also details.
#' @param sensor Numeric. Type of landsat satellite; it takes 5,7, or 8
#' @param gain Band-specific sensor gain. It takes either "high" or "low"
#' @inheritParams sebal
#' @inheritParams sebkcstack
#' @param weights albedo bands weights for each band.
#' list as in order c(band1,band2,band3,band4,band5,band6,band7) 
#' for landsat 5 and 7.
#'  The default is "auto" where 
#' sebal weights are used.
#' @param grescale list as in order weights. Band-specific sensor values. 
#' The default for landsat 5 and 7 is from Chander et.al (2009).
#' For Landsat 8 use ML Radiance scale in the metadata;
#' i.e. Band specific multiplicative rescaling factor from the metadata 
#' (MTL file) (RADIANCE_MULT_BAND_x, where x is the band number).
#' @param brescale list. Band-specific sensor values.
#' The default for landsat 5 and 7 is from Chander et.al (2009).
#' for Landsat 8 use AL Radiance scale in the metadata;
#' i.e. Band specific multiplicative rescaling factor from the metadata 
#' (MTL file) (RADIANCE_MULT_BAND_x, where x is the band number). 
#' @param ESUNx list. The mean solar exo-atmospheric irradiance for each band.
#' It is used for reflectance correction. The default is taken from SEBAL manual
#' @param K1 Temperature constants for each Landsat image. 
#' Default is taken from landsat handbook
#' @param K2 Temperature constants for each Landsat image. 
#' Default is taken from landsat handbook
#' @param Rp path radiance in the 10.4-12.5 um band.
#' MODTRAN runs can accurately estimate this. 
#' The default is taken from Allen et al(2007) 
#' to be 0.91,
#' @param Rsky narrow band downward thermal radiation. 
#' MODTRAN runs can accurately estimate this.
#' The default is taken from Allen et al(2007) to be 1.32 
#' @param tNB narrow band transmissivity of air (10.4-12.5) um range.
#' MODTRAN runs can accurately estimate this.
#' The default is taken from Allen et al(2007) to be 0.866 
#' @param date DATE_ACQUIRED of the image, see metadata
#' @param Mp Reflectance scale for Landsat 8. Band specific multiplicative 
#' rescaling factor from the metadata (MTL file) (REFLECTANCE_MULT_BAND_x, 
#' where x is the band number).
#' @param Ap Reflectance scale for Landsat 8. Band specific additive rescaling
#' factor from the metadata (MTL file) (REFLECTANCE_ADD_BAND_x,
#' where x is the band number)
#' @author George Owusu
#' @details If you set the data to landsat folder that contains 
#' landsat band and meta data \code{\link{sebkcstack}} is used to
#' compute the parameters. The user must specify welev and may
#' use default parameters. 
#'  The data must always be a raster brick 
#' or stack. To use the default values please order bind the bands
#'  c(band1,band2,band3,band4,band5,band6,band7). 
#'  The defaults values are for landsat 5 and 7. 
#' For landsat 8, use the meta data. 
#' @return 
#'  \itemize{
#' \item{albedo:} { albedo values}
#' \item{Ts:} { Surface Temperature}
#' \item{NDVI:} { NDVI values}
#' \item{SAVI:} { SAVI values}
#' \item{fc:} { Fractional vegetation cover}
#' \item{radiance:} { radiance values}
#' \item{reflectance:} { TOA reflectance}
#' }
#' @references 
#'  \itemize{
#' \item{}{Chander, G., Markham, B. L., & Helder, D. L. 2009. 
#' Summary of current radiometric calibration coefficients for 
#' Landsat MSS, TM, ETM+, and EO-1 ALI sensors. 
#' Remote Sensing of Environment, 113(5): 893-903.}
#' \item{}{Waters, R. 2002. SEBAL: Surface Energy Balance Algorithms for Land. 
#' Idaho Implementation, Advanced Training and Users Manual. Kimberly, Idaho: 
#' University of Idaho.}
#' \item{}{Survey, U. S. G. 2015. Landsat 8 (L8) data users handbook: 97p.
#'  from https://landsat.usgs.gov/documents/Landsat8DataUsersHandbook.pdf}
#' }
#' @examples
#' \dontrun{
#' #Automatic detection of parameters
#' folder=system.file("extdata","stack",package="SEBKc")
#' modauto=landsat578(data=folder, welev=362)
#' 
#' #Manual estimation with input parameters
#' landsat7_6feb2004=brick(system.file("extdata","landsat7_6feb2004.grd",
#' package="SEBKc"))
#' mod=landsat578(data=landsat7_6feb2004,welev=317.1,sensor=7,gain="high",
#' sunelev=50.71154048,weights="auto",grescale="auto",brescale="auto",
#' ESUN="auto",K1=666.09,K2=1282.71,date="2002-2-6")
#' #acces albedo value
#' albedo=mod$albedo
#' }
#' @export
#' @rdname landsat578

landsat578<-function(data,welev,sensor=NULL,gain="low",
sunelev=NULL,date=NULL,weights="auto",grescale="auto",brescale="auto",
ESUNx="auto",Mp=NULL,Ap=NULL,K1=NULL,K2=NULL,Rp=0.91,Rsky=1.32,tNB=0.886,
DOY=NULL,clip=NULL,folder=NULL) UseMethod ("landsat578")
#' @export
#' @rdname landsat578
landsat578.default<-function(data=NULL,welev,sensor=NULL,gain="low",
                             sunelev=NULL,date=NULL,
weights="auto",grescale="auto",brescale="auto",ESUNx="auto",Mp=NULL,
Ap=NULL,K1=NULL,K2=NULL,Rp=0.91,Rsky=1.32,tNB=0.886,DOY=NULL,clip=NULL,folder=NULL){
  stack=NULL
  
  if(is.null(data)&&!is.null(folder)){
   data=folder 
  }
  if(is.null(data)&&is.null(folder)){
  return(print("Please provide data or folder path"))
      }
  folder=NULL

  if(class(data)[1]=="character"||class(data)=="sebkcstack"){
    folder=data
    if(class(data)=="sebkcstack"){
      file.info2=TRUE
    }else{
      file.info2=file.info(folder)[["isdir"]]
    }
    if(file.info2==TRUE&&!is.na(file.info2)){
      if(class(data)=="sebkcstack"){
        stack=data
      }else{
        stack=sebkcstack(data,clip=clip) 
      }
      data=stack$data
      sensor=stack$sensor
      satellite=sensor
      sunelev=stack$sunelev
      date=stack$date
      #gain="high"
      grescale=stack$grescale
      brescale=stack$brescale
      Ap=stack$Ap
      Ap[is.na(Ap)]=0
      Mp=stack$Mp
      Mp[is.na(Mp)]=0
      if(sensor==8){
      K1=stack$K1
      K2=stack$K2
      K11=stack$K11
      K21=stack$K21
      }
        }
  }
  
satellite=sensor
thisdate=date
if(!is.numeric(DOY)){
DOY=strptime(thisdate,"%Y-%m-%d")$yday+1
if(is.na(DOY)){
  DOY=strptime(thisdate,"%Y-%m-%d")$yday+1
}
}



dr=1+0.033*cos(DOY*((2*pi)/365))
cos_teta=cos(((90-sunelev)/360)*2*pi)
tsw=0.75+2*1e-05*welev
weights5=c(0.293,0.274,0.233,0.157,0.033,0,0.011)
weights7=c(0.293,0.274,0.231,0.156,0.034,0,0.012)
weights8=c(0.293,0.274,0.231,0.156,0.034,0.012,0)
#landsat 8 weights from Silva et al (2016): 
#DOI: http://dx.doi.org/10.1590/1807-1929/agriambi.v20n1p3-8
weights8=c(0.300,0.277,0.233,0.143,0.035,0.012,0)


lowGrescale7=c(1.180709,1.209843,0.942520,0.969291,
0.191220,0.067087,0.066496)
lowBrescale7=c(-7.38,-7.61,-5.94,-6.07,-1.19,-0.07,-0.42)
ESUNx7=c(1997,1812,1533,1039,230.8,1,84.9)

highGrescale7=c(0.77874,0.798819,0.621654,0.639764,
0.12622,0.037205,0.043898)
highBrescale7=c(-6.98,-7.2,-5.62,-5.74,-1.13,3.16,-0.39)

##############################
Grescale5after=c(0.765827,1.448189,1.043976,0.876024,0.120354,
0.055376,0.065551)
Brescale5after=c(-2.29,-4.29,-2.21,-2.39,-0.49,1.18,-0.22)
ESUNx5=c(1983,1796,1536,1031,220,1,83.44)

Grescale5before=c(0.671339,1.322205,1.043976,0.876024,0.120354,
0.055376,0.065551)
Brescale5before=c(-2.19,-4.16,-2.21,-2.39,-0.49,1.18,-0.22)

if(satellite==7||satellite=="ETM+"||satellite=="7"){
ESUN=ESUNx7
weight=weights7
if(is.null(K1)){
K1=666.09
}
if(is.null(K2)){
  K2=1282.71
}
if(gain=="high"){
Grescale=highGrescale7
Brescale=highBrescale7
}
if(gain=="low"){
Grescale=lowGrescale7
Brescale=lowBrescale7
}
}
#thisdate=paste(year,month,day,sep="-")
beforedate = as.POSIXct("1991-12-31",format='%Y-%m-%d')
afterdate = as.POSIXct(thisdate,format='%Y-%m-%d')
timediff=as.numeric(difftime(afterdate,beforedate,units="days"))
if(satellite==5||satellite=="MSS"||satellite=="5"){
ESUN=ESUNx5
weight=weights5
if(is.null(K1)){
  K1=607.76
}
if(is.null(K2)){
  K2=1260.56
}

if(timediff>0){
Grescale=Grescale5after
Brescale=Brescale5after
}
if(timediff<0){
Grescale=Grescale5before
Brescale=Brescale5before
}
}
#Radiance
if(is.numeric(grescale)){
Grescale=grescale

}

if(is.numeric(brescale)){
Brescale=brescale

}
print("Computing Radiance: this may take several minutes........")

radiance=(Grescale*data)+Brescale

#reflection
if(ESUNx!="auto"){
ESUN=ESUNx
}

print("Computing Reflectance")

if(satellite==8||satellite=="8"){
weight=weights8
Mp[7]=0
Ap[7]=0
reflection=((Mp*data)+Ap)/sin(sunelev*(pi/180))
#plot(data)
#print(Ap)

}else{
  reflection=(pi*radiance)/(ESUN*cos_teta*dr)
}
#albedo weights
if(!is.numeric(weights)){
if(weights!="auto"){
weight=weights
}
}
weights=reflection*weight

print("Computing Albedo")
#albedo toa
albedo_toa=sum(weights)

#albedo
albedo=(albedo_toa-0.03)/(tsw^2)

print("Computing NDVI")

NDVI=(reflection[[4]]-reflection[[3]])/(reflection[[4]]+reflection[[3]])
fc=((NDVI-cellStats(NDVI,min))/(cellStats(NDVI,max)-cellStats(NDVI,min)))^2
NDVI=rastercon(NDVI< (-1),-1,rastercon(NDVI>1,1,NDVI))
print("Computing SAVI")
SAVI=(1.5*(reflection[[4]]-reflection[[3]]))/
(0.5+(reflection[[4]]+reflection[[3]]))
print("Computing LAI")
LAI=(5.434*SAVI)+1.5416
print("Computing emissivities")
eNB=rastercon( NDVI<0 & albedo<0.47,0.99,rastercon(LAI>=3,0.98,0.97+
                                                     (LAI*0.0033)))
#Rc=radiance[[6]]

Ts=NULL
#temperature in Kelvin
if(satellite==8)
{
radiance6=radiance[[7]]  
}else{
radiance6=radiance[[6]]  
}
if(!is.null(K1)&&!is.null(K2)){
  print("Computing corrected thermal radiance")
Rc=((radiance6-Rp)/tNB)-(1-eNB)*Rsky
print("Computing Surface temperature")
Ts=K2/(log(((eNB*K1)/Rc)+1))
}

factor<-list(albedo=albedo,Ts=Ts,NDVI=NDVI,SAVI=SAVI,fc=fc,radiance=radiance,
             reflectance=reflection,DOY=DOY,sunelev=sunelev,welev=welev,
             date=date,LAI=LAI,eNB=eNB,time=stack$time,date=stack$date,folder=folder )


factor$call<-match.call()

class(factor)<-"landsat578"
factor
}