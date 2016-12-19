#' @title Computing Two Source Surface Energy Balance (TSEB) in R
#'
#' @description This function estimates evaporation (E), transpiration (T), 
#' and evapotranspiration(ET) of land surface. The computation is based on
#' Norman et al (1995) and Colaizzi et al(2012) models.   
#' @inheritParams sebal
#' @inheritParams ETo 
#' @inheritParams coldTs
#' @inheritParams ETohr
#' @inheritParams sebkcstack
#' @param TA Air temperature[K]. It can take "sebal" or 'Tmax' or 'Tmin' or numeric value.
#' When it is set to sebal, TA = Ts-dT. The 'Tmax' and 'Tmin' correspond to maximum and 
#' minimum temperature.
#' @param viewangle angle of view of the sensor in degrees
#' @param s numeric. Mean leaf size; it is four times the leaf 
#' area divided by the perimeter. The default is 2 (Grace et al., 1981).
#' @param C coefficient. Weighting a coefficient in the equation for leaf 
#' boundary layer. A value of 90 s^0.5/m (Norman et al, 1995).
#' @param xPT coefficient(s). The Priestley-Taylor equation coefficient (Norman et al, 2005).
#' It can take a single coefficient or in the form c(initial, minimum, increment). 
#' Example is xPT=c(1.26,1,0.1). It must be noted this is a decreasing increment 
#' from the initial value .i.e 1.26. The subsequent increment is only meant for 
#' a cell where LEs is negative. 
#' @param clump coefficient. A canopy clumping coefficient used for 
#' calculation of directional  fc (fractional cover of vegetation)
#' @param fg numeric. The fraction of LA1 that is green. If no information is
#'  available on fg, then it is assumed to be unity (1).
#' @param fc fraction (0-1).Directional vegetation cover. 
#' The default uses Norman et al(1995) fc equation 
#' @param hc numeric. it is canopy height[m]
#' @param network character. It takes either "parallel" (Norman et al, 1995)  or 
#' "series" (Colaizzi et al, 2012; Norman et al, 1995).
#' @param series character.The type of series. It takes "xPT" for 
#' Priestley-Taylor equation(Norman, 1995) or 
#' "PM" for Penman-Monteith equation (Colaizzi et al, 2012)
#' @param rc coefficient. Penman-Monteith equation 
#' coefficient (Colaizzi et al,2012). It can take a single coefficient 
#' or in the form c(initial, maximum, increment). 
#' Example is rc=c(0,1000,100). Unlike xPT, this is positively incremented.
#' @author George Owusu
#' @references 
#' \itemize{
#' \item{}{Norman JM, Kustas WP, Humes KS (1995). A two-source approach for estimating 
#' soil and vegetation energy fluxes from observations of directional 
#' radiometric surface temperature. Agric For Meteorol 1995;77:263-93.}
#' \item{}{Colaizzi, Paul; Kustas, William P.; Anderson, Martha C.; 
#' Agam, Nurit; Tolk, Judy A.; Evett, Steven R.; Howell, Terry A.; 
#' Gowda, Prasanna H.; and O'Shaughnessy, Susan A., 
#' "Two-source energy balance model estimates of evapotranspiration using 
#' component and composite surface temperatures" (2012). Publications 
#' from USDA-ARS / UNL Faculty. Paper 1147.
#' http://digitalcommons.unl.edu/usdaarsfacpub/1147}
#' \item{}{Grace, J., 1981. Some effects of wind on plants. 
#' In: J. Grace, E.D. Ford and P.G. Jarvis (Editors), Plants and
#' their Atmospheric Environment. Blackwell Scientific, London, pp. 31-56.}
#' \item{}{Ronglin Tang, Zhao-Liang Li, Yuanyuan Jia, Chuanrong Li, 
#' Kun-Shan Chen, Xiaomin Sun & Jinyong Lou (2012): 
#' Evaluating one- and two-source energy balance models in
#' estimating surface evapotranspiration from Landsat-derived 
#' surface temperature and field
#' measurements, International Journal of Remote Sensing.
#' http://dx.doi.org/10.1080/01431161.2012.716529}
#' \item{}{KUSTAS, W. & ANDERSON, M. 2009. Advances in thermal infrared remote 
#' sensing for land surface modeling. Agricultural and Forest Meteorology, 149,
#'  2071-2081.}
#'}
#' @details 
#' \describe{
#'  \item{\strong{Theory}}{This function separately computes LEs, LEs using
#'  Net radiation (Rn), sensible heat flux (H) and soil heat flux (G).
#'  The Rn was first computed with \code{sebal}. The soil net radiation computed
#'  as  \eqn{Rns=(sebal$Rn*exp(0.9*log(1-fc)))} . The canopy radiation was computed 
#'  as  \eqn{Rnc=sebal$Rn-Rns} . The G was computed from \code{sebal}. 
#'  The canopy sensible heat (Hc) and soil sensible heat (Hs) were computed
#'  based on corresponding canopy Temperature (Tc) and soil temperature (Tsoil).
#'  Two different approaches were employed to sensible heat fluxes. The first
#'   is Norman et al(1995) parallel network. The second is Norman (1995) 
#'   secant series but uses Priestley-Taylor equation (xPT).   }
#'  \item{\strong{Default Parameter}}{These parameters have been 
#'  set based on literature }
#'  }
#' @return 
#' \itemize{
#' \item{LEc:} { canopy latent heat flux [W/m2]}
#' \item{LEs:} { Soil latent heat flux[W/m2]}
#' \item{T24:} { 24 hour Transpiration [mm/day]}
#' \item{E24:} { 24 hour Evaporation [mm/day]}
#' \item{ET24:} { 24 hour Evapotranspiration [mm/day]}
#' \item{EF:} { ET24 (daily Evapotranspiration) Evaporative Fraction. 
#' EF corresponds to the type of scaling method used: see \code{ETohr} for details}
#' \item{EFs:} { E(Evaporation) Evaporative Fraction}
#' \item{EFc: } { T (Transpiration) Evaporative Fraction}
#' } 
#' @examples
#' \dontrun{
#' albedo=raster(system.file("extdata","albedo.grd",package="SEBKc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="SEBKc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="SEBKc"))
#' LAI=raster(system.file("extdata","LAI.grd",package="SEBKc"))
#' Parallel=tseb(Ts=Ts,LAI=LAI,DOY=37,xyhot="full",
#' xycold="full",albedo=albedo,Tmax=31,
#' Tmin=28,RHmax=NULL,RHmin=NULL,zom=NULL,NDVI=NDVI,SAVI=NULL,hc=20,
#' DEM=NULL,viewangle=0,sunelev=50,welev=345,u=2,
#' s=2,C=90,lapse=0.0065,Rn24=NULL,zx=200,zomw=2,
#' xPT=1.3,clump=1,fg=1,fc="auto",network="parallel",
#' latitude=5.6,n=6.5)
#' 
#' series=tseb(Ts=Ts,LAI=LAI,DOY=37,xyhot="full",
#' xycold="full",albedo=albedo,Tmax=31,
#' Tmin=28,RHmax=NULL,RHmin=NULL,zom=NULL,NDVI=NDVI,SAVI=NULL,hc=20,
#' DEM=NULL,viewangle=0,sunelev=50,welev=345,u=2,
#' s=2,C=90,lapse=0.0065,Rn24=NULL,zx=200,zomw=2,
#' xPT=1.3,clump=1,fg=1,fc="auto",network="series",
#' latitude=5.6,n=6.5)
#' 
#' #using landsat folder parameters
#' folder=system.file("extdata","stack",package="SEBKc")
#' tsebautoseries=tseb(folder=folder,welev=317,Tmax=31,Tmin=27,
#' latitude=5.6,n=6.5)
#' tsebautoparallel=tseb(folder=folder,welev=317,Tmax=31,Tmin=27,
#' latitude=5.6,n=6.5,network="parallel")
#' 
#' #3 step simulation using landsat 7 with SLC OFF
#' ##step one: prepare the input data 
#' ### point to landsat 7 data with strips 
#'
#' #####get the data stacked and remove strips
#' rawdata=sebkcstack(folder)
#' ###### check if strips are removed
#' plot(rawdata$data)
#' #get Ts, Albedo, NDVI, SAVI, sunelev, DOY from landsat data 
#' welev=278
#'  data=landsat578(rawdata,welev=welev)
#'  #perform semi-auto simulation 
#' # Determine xyhot. Digitize polygon on the Ts map
#' modhot=hotTs(data,welev=300,extent="digitize",cluster=2)
#' #determine the cold. Digitize polygon on the Ts map
#' modcold=coldTs(data,welev=275,extent="digitize",cluster=2)
#' xyhot=modhot$xyhot
#' xycold=modcold$xycold
#' modTSEB2=tseb(data,Tmax=31,Tmin=28,
#' latitude=5.6,n=6.5)
#' #the same
#' modTSEB2=tseb(data,Tmax=33,Tmin=28,welev=345,xyhot=xyhot,xycold=xycold,
#' latitude=5.6,n=6.5)
#' #the same
#' modTSEB2=tseb(data,Tmax=33,Tmin=28,welev=345,xyhot=modhot,xycold=modcold,
#' latitude=5.6,n=6.5)

#' }
#' 

#' @export 
#' @rdname tseb
tseb=function(albedo,Ts,LAI,DOY,xyhot="auto",xycold="auto",TA='sebal'
,Tmax,Tmin,RHmax=NULL,RHmin=NULL,n=NULL,latitude=NULL,zom=NULL,NDVI,SAVI,hc=2,DEM=NULL,
viewangle=0,sunelev,welev,u=2,s=2,C=90,lapse=0.0065,Rn24=NULL,ETr=NULL,ETr24=NULL,wmo=NULL, 
airport=NULL,Krs=0.16,surface = "grass",time=NULL, Lz=NULL,Lm=NULL,t1=1,zx,zomw,
xPT=1.3,rc=c(70,1000,100),series="xPT",clump=1,fg=1,fc="auto",
network="parallel",clip=NULL,folder=NULL,iter.max = 7) UseMethod ("tseb")
#' @export 
#' @rdname tseb
tseb.default=function(albedo=NULL,Ts=NULL,LAI=NULL,DOY=NULL,
                      xyhot="auto",xycold="auto",TA='sebal',
Tmax,Tmin,RHmax=NULL,RHmin=NULL,n=NULL,latitude=NULL,zom=NULL,
NDVI,SAVI,hc=2,DEM=NULL,viewangle=0,sunelev,welev,u=2,
s=2,C=90,lapse=0.0065,Rn24=NULL,ETr=NULL,ETr24=NULL,wmo=NULL, 
airport=NULL,Krs=0.16,surface = "grass",
time=NULL, Lz=NULL,Lm=NULL,t1=1,zx=200,zomw=2,
xPT=c(1.26,0,0.1),rc=c(50,1000,100),series="PM",clump=1,fg=1,
fc="auto",network="series",clip=NULL,folder=NULL,iter.max = 7){
  PET=NULL
  thisTA=TA
  
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
      folder=mod
    }
    #mod=landsat578(data=folder, welev=welev)
    albedo=mod$albedo
    Ts=mod$Ts
    NDVI=mod$NDVI
    DOY=mod$DOY
    sunelev=mod$sunelev
    SAVI=mod$SAVI
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
      Tmin=ETr24$data$Tmin
      Tmax=ETr24$data$Tmax
      PET=ETr24
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
        if(thisTA=="Tmean"){
        TA=273.15+ETr$data$Tmean
        #print(TA)
        }
        ETr=ETr$ETo
      }
    }
  }
  
#print(ETr24)
zu=zx  
U=u

a1=0.004
b1=0.012
d=0.65*hc
ZM=hc/8
ZM=0.018*LAI
cp=1004.14
iteration=F
mod=sebal(albedo=albedo,Ts=Ts,NDVI=NDVI,SAVI=SAVI,iter.max=iter.max,
xyhot=xyhot,xycold=xycold,DOY=DOY,sunelev=sunelev,welev=welev,
zx=zx,u=u,zomw=zomw,zom=zom,LAI=LAI,DEM=DEM,lapse=lapse,
Rn24=Rn24,ETr=ETr,model="SEBAL",folder=folder)

if(is.null(PET)){
PET=ETo(Tmax=Tmax,Tmin=Tmin,z=zx,uz=u,altitude=welev,
RHmax=RHmax,RHmin=RHmin,n=n,DOY=mod$DOY,latitude=latitude)
}
Rn24=PET$Rn
y=PET$y
vpd=PET$vpd
slope=PET$slope
pair=mod$pair
G=mod$G
vaporisation=mod$vaporisation
#G=0
RA=mod$rah

if(thisTA=='sebal'||thisTA=='SEBAL'){
  TA=mod$Ta
}

if(thisTA=='Tmax'||thisTA=='TMAX'){
  TA=273.15+Tmax
}
if(thisTA=='Tmin'||thisTA=='TMIN'){
  TA=273.15+Tmin
}
#Ts=mod$Ts

pair=mod$pair
wm=mod$wm
ladt=0.5
angle=viewangle
print("Computing fc based on LAI, angle, clump, etc")

if(fc=="auto"||!is.numeric(fc)){
fc=min(1-exp((-ladt*clump*LAI)/(cos(angle*(pi/180)))),1)
}
print("Computing Soil surface net radiation")

Rns=(mod$Rn*exp(0.9*log(1-fc)))

G=mod$G_by_Rn*Rns
print("Computing canopy net radiation")
Rnc=mod$Rn-Rns
print("Computing transport resistance (Rs)and  wind speed above surface (Us)")

Uc=U*(log((hc-d)/ZM)/(log((zu-d)/ZM)-wm))
a=0.28*LAI^(2/3)*hc^(1/3)*s^(-1/3)
Us=Uc*exp(-1*a*(1-(0.05/hc)))
Rs=1/(a1+(b1*Us))
Udzm=Uc*exp(-1*a*(1-((d+ZM)/hc)))
Rx=(C/LAI)*((s/Udzm)^0.5)
steady=NULL
if(network=="series"){
  print("Network series started")
if(series=="PT"||series=="xPT"){
  print("series = xPT or PT:Priestley-Taylor")
  
  thisxPT=xPT
if(is.na(xPT[2])){
  thisxPT=c(xPT,xPT,xPT)
}
  H2=NULL
  i=thisxPT[1]
  while(i>=thisxPT[2]){
    print(paste("Trying PT=",i,"to correct soil surface latent heat"))
    
TClin=((TA/RA)+(Ts/(Rs*(1-fc)))+((Rnc*Rx)/(pair*cp))*
(1-(i*fg*(slope/(slope+y))))*((1/RA)+(1/Rs)+(1/Rx)))/
  ((1/RA)+(1/Rs)+(fc/(Rs*(1-fc))))
Tslin=TClin*(1+(Rs/RA))-((Rnc*Rx)/(pair*cp))*
  (1-(i*fg*(slope/(slope+y))))*(1+(Rs/Rx)+(Rs/RA))-(TA*(Rs/RA))
dTC=((Ts^4)-(fc*TClin^4)-((1-fc)*Tslin^4))/((4*(1-fc)*Tslin^3)*
    (1+(Rs/RA))+(4*fc*TClin^3))
Tc=TClin+dTC
Tsoil=(((Ts^4)-(fc*Tc^4))/(1-fc))^0.25
TAC=((TA/RA)+(Tsoil/Rs)+(Tc/Rx))/((1/RA)+(1/Rs)+(1/Rx))
Hs=pair*cp*((Tsoil-TAC)/Rs)
Hc=pair*cp*((Tc-TAC)/Rx)
H=pair*cp*((TAC-TA)/RA)
if(i==thisxPT[1]){
LEs=Rns-G-Hs
LEc=Rnc-Hc
#H2=H
}else{
  LEc=rastercon(LEs<=0,Rnc-Hc,LEc)
  #H2==rastercon(LEs<=0,H,H2)
  LEs=rastercon(LEs<=0,Rns-G-Hs,LEs)
}
#print(i)
i=i-thisxPT[3]

if(minValue(LEs)>0){
  i=thisxPT[2]-i 
  steady=TRUE
}
}
#LEs2=LEs
#LEs2[LEs2>0,]=NA
#H=H2
##if LEs<0
  if(is.null(steady)){
HsEso=Rns-G

TAClin=((TA/RA)+(Ts/(fc*Rx))-((1-fc)/(fc*Rx))*((HsEso*Rs)/(pair*cp))+
          (HsEso/(pair*cp)))/((1/RA)+(1/Rx)+((1-fc)/(fc*Rx)))
TE=TAClin*(1+(Rx/RA))-((HsEso*Rx)/(pair*cp))-((TA*Rx)/RA)
dTAC=((Ts^4)-(1-fc)*(((HsEso*Rs)/(pair*cp))+(TAClin))^4-(fc*(TE^4)))/
  ((4*fc*(TE^3))*(1+(Rx/RA))+(4*(1-fc))*(((HsEso*Rs)/(pair*cp))+(TAClin))^3)
TAC=TAClin+dTAC
Tsoil2=TAC+((HsEso*Rs)/(pair*cp))
Hs2=pair*cp*((Tsoil2-TAC)/Rs)
LEc=rastercon(LEs<=0,Rnc-Hc,LEc)
LEs=rastercon(LEs<=0,Rns-G-Hs2,LEs)

if(minValue(LEs)<0){
  LEs=rastercon(LEs<=0,0,LEs)
}
}
}
if(series=="PM"||series=="penman-monteith"||series=="pm"){
  thisrc=rc
  if(is.na(rc[2])){
    thisrc=c(rc,rc,rc)
  }
  #H2=NULL
  print("series = PM:Penman-Monteith")
  
  i=thisrc[1]
  while(i<=thisrc[2]){
    print(paste("Trying rc=",i,"to correct soil surface latent heat"))
    
y_star=y*(1+(i/RA))
#print(i)
TClin=((TA/RA)+(Ts/(Rs*(1-fc)))+(((Rx*y_star*Rnc)/(pair*cp*(slope+y_star)))-
((Rx/RA)*(vpd/(slope+y_star))))*((1/RA)+(1/Rs)+(1/Rx)))/((1/RA)+(1/Rs)+(fc/
(Rs*(1-fc))))
Tslin=TClin*(1+(Rs/RA))-(TA*(Rs/RA))-(((Rx*y_star*Rnc)/(pair*cp*(slope+y_star)))-
  ((Rx/RA)*(vpd/(slope+y_star))))*(1+(Rs/Rx)+(Rs/RA))

dTC=((Ts^4)-(fc*TClin^4)-((1-fc)*Tslin^4))/((4*(1-fc)*Tslin^3)*(1+(Rs/RA))+
                                              (4*fc*TClin^3))
Tc=TClin+dTC
Tsoil=(((Ts^4)-(fc*Tc^4))/(1-fc))^0.25
TAC=((TA/RA)+(Tsoil/Rs)+(Tc/Rx))/((1/RA)+(1/Rs)+(1/Rx))

Hs=pair*cp*((Tsoil-TAC)/Rs)
Hc=pair*cp*((Tc-TAC)/Rx)
H=pair*cp*((TAC-TA)/RA)
#LEc=Rnc-Hc

if(i==thisrc[1]){
  LEs=Rns-G-Hs
  LEc=Rnc-Hc
  #H2=H
}else{
  LEc=rastercon(LEs<=0,Rnc-Hc,LEc)
  #H2=rastercon(LEs<=0,H,H2)
  LEs=rastercon(LEs<=0,Rns-G-Hs,LEs)
  #LEc[LEs<=0,]=Rns-G-Hs
  #LEc[LEs<=0,]=Rnc-Hc
}

i=i+thisrc[3]
if(minValue(LEs)>0){
i=i+thisrc[2]  
steady=TRUE
}
  }
  
 
  
#LEs2[LEs2>0,]=NA
#H=H2
##if LEs<0
  if(is.null(steady)){
    
HsEso=Rns-G
          
TAClin=((TA/RA)+(Ts/(fc*Rx))+(HsEso/(pair*cp))-((HsEso*Rs)/(pair*cp))*((1-fc)/(fc*Rx)))/
  ((1/RA)+(1/Rx)+((1-fc)/(fc*Rx)))
TE=TAClin*(1+(Rx/RA))-(TA*(Rx/RA))-((HsEso*Rx)/(pair*cp))
dTAC=((Ts^4)-(1-fc)*(((HsEso*Rs)/(pair*cp))+(TAClin))^4-(fc*(TE^4)))/
  ((4*fc*(TE^3))*(1+(Rx/RA))+(4*(1-fc))*(((HsEso*Rs)/(pair*cp))+(TAClin))^3)
TAC=TAClin+dTAC
Tsoil2=TAC+((HsEso*Rs)/(pair*cp))
Hs2=pair*cp*((Tsoil2-TAC)/Rs)
LEs=rastercon(LEs<=0,Rns-G-Hs2,LEs)
if(minValue(LEs)<0){
LEs=rastercon(LEs<0,0,LEs)
}
}
}
}else{
  thisxPT=xPT
  if(is.na(xPT[2])){
    thisxPT=c(xPT,xPT,xPT)
  }
  #thisxPT=c(xPT[1],xPT[1],xPT[1])
 LEc2=NULL
 LEs2=NULL
 #H2=NULL
 #G2=NULL
 print(paste("network =",netwok))
  i=thisxPT[1]
  while(i>=thisxPT[2]){
    print(paste("Trying PT=",i,"to correct soil surface latent heat"))
    
LEc=i*fg*(slope/(slope+y))*Rnc
Hc=Rnc*(1-(i*fg*(slope/(slope+y))))
Tc=((Hc*RA)/(pair*cp))+TA
Tsoil=(((Ts^4)-(fc*Tc^4))/(1-fc))^0.25
Hs=pair*cp*((Tsoil-TA)/(Rs+RA))
if(i==thisxPT[1]){
LEs=Rns-Hs-G
}else{
LEs=rastercon(LEs<=0,Rns-G-Hs,LEs)
}
Hs=rastercon(LEs<0,Rns-G,Hs)
#Tsoil2=(Hs2/(pair*cp))*(Rs+RA)+TA
Tsoil=rastercon(LEs<0,(Hs/(pair*cp))*(Rs+RA)+TA,Tsoil)
#Tc2=(((Ts^4)-(1-fc)*Tsoil2^4)/fc)^0.25
Tc=rastercon(LEs<0,(((Ts^4)-(1-fc)*Tsoil^4)/fc)^0.25,Tc)
#Hc2=pair*cp*((Tc2-TA)/RA)
Hc=rastercon(LEs<0,pair*cp*((Tc-TA)/RA),Hc)
#LEc2=Rnc-Hc2
LEc=rastercon(LEs<0,Rnc-(pair*cp*((Tc-TA)/RA)),LEc)
Hc=rastercon(Hc>Rnc,Rnc,Hc)
Tc=rastercon(Hc>Rnc,((Hc*RA)/(pair*cp))+TA,Tc)
Tsoil=rastercon(Hc>Rnc,(((Ts^4)-(fc*Tc^4))/(1-fc))^0.25,Tsoil)
Hs=rastercon(Hc>Rnc,pair*cp*((Tsoil-TA)/(Rs+RA)),Hs)
G=rastercon(Hc>Rnc,Rns-Hs,G)
Rstar=(Ts-TA)/(((Tc-TA)/RA)+((Tsoil-TA)/(RA+Rs)))
#H=pair*cp*((Ts-TA)/Rstar)
H=Hs+Hc
if(i==thisxPT[1]){
  LEs2=Rns-G-Hs
  LEc2=Rnc-Hc
  #H2=H
  #G2=G
  #H3=H
}else{
  LEc2=rastercon(LEs2<=0,Rnc-Hc,LEc2)
  LEs2=rastercon(LEs2<=0,Rns-G-Hs,LEs2)
  
  #G2=rastercon(LEs2<=0,G,G2)
  #H2=rastercon(LEs2<=0,H,H2)
}


i=i-thisxPT[3]

  }
  #H=H2
  #G=G2
  LEc=LEc2
  LEs=LEs2
  if(minValue(LEs)<0){
    LEs=rastercon(LEs<=0,0,LEs)
  }
}
#print(LEs)
#print(spplot(stack(LEc,LEs)))
EFc=EFs=EF=NULL
#LE=mod$Rn-G-H
LE=LEc+LEs
#print(spplot(LE-LE2))
print ("computing E, T and ET" )

if(is.null(ETr)){
EFc=1.1*(LEc/Rnc)
EFs=1.1*(LEs/(Rns-G))
#EF=1.1*((LEc+LEs)/((Rns+Rnc)-G))
EF=1.1*(LE/(mod$Rn-G))

#EF=mean(EFc)
EFc[EFc<0,]=0
EFc[EFc>1.1,]=1.1
EFs[EFs<0,]=0
EFs[EFs>1.1,]=1.1

EF[EF<0,]=0
EF[EF>1.1,]=1.1
ET24=NULL
T24=NULL
E24=NULL
if(!is.null(Rn24)){
  T24=max(EFc*Rn24*0.035,0)
  E24=max(EFs*Rn24*0.035,0)
  ET24=max(EF*Rn24*0.035,0)
}
}else{
  T24=(t1*3600*LEc)/vaporisation
  EFc=T24/ETr
  E24=(t1*3600*LEs)/vaporisation
  EFs=E24/ETr
  ET24=(t1*3600*LE)/vaporisation
  EF=ET24/ETr
  if(minValue(EF)<0){
    EF[EF<0,]=0
  }
  if(maxValue(EF)>1.1){
    EF[EF>1.1,]=1.1
  }
  
  if(minValue(EFs)<0){
    EFs[EFs<0,]=0
  }
  if(maxValue(EFs)>1.1){
    EFs[EFs>1.1,]=1.1
  }
  
  if(minValue(EFc)<0){
    EFc[EFc<0,]=0
  }
  if(maxValue(EFc)>1.1){
    EFc[EFc>1.1,]=1.1
  }
  
  if(!is.null(ETr24)){
    E24=EFs*ETr24
    T24=EFc*ETr24
    ET24=EF*ETr24
    }
}
factor<-list(LE=LE,LEc=LEc,LEs=LEs,T24=T24,E24=E24,ET24=ET24,EF=EF,
EFc=EFc,EFs=EFs,ETr=ETr,fc=fc,vaporisation=vaporisation,folder=folder,t1=t1,Rs=mod$Rs)
factor$call<-match.call()

class(factor)<-"tseb"
factor
}