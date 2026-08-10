
#' @title SEBS Estimation of Evapotranspiration Fraction
#'
#' @description \code{sebs} computes Su (2002) 
#' Surface Energy Balance System (SEBS) for estimation of turbulent heat fluxes
#' @inheritParams sebal
#' @inheritParams ETo
#' @inheritParams ETohr
#' @inheritParams weather
#' @inheritParams sebkcstack
#' @details Function \code{\link{ETo}} will be automatically used to 
#' estimate  parameters y (psychrometric constant), 
#' slope (slope of saturation vapour pressure curve) and 
#' vpd (Vapour pressure deficit). For information on data input procedure
#' see 'details' at \code{\link{sebal}}
#' @seealso \code{\link{sebal}} , \code{\link{ETo}}, \code{\link{sebi}}  
#' @references 
#' Su, Z. 2002. The Surface Energy Balance System (SEBS) for estimation 
#' of turbulent heat fluxes. Hydrol. Earth Syst. Sci., 6(1): 85-100.
#' @return EF Evapotranspiration Fraction
#' @author George Owusu
#'
#' @examples
#' \dontrun{
#' #define Input data
#' albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="sebkc"))
#' NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
#' LAI=raster(system.file("extdata","LAI.grd",package="sebkc"))
#' Tmax=31
#' Tmin=28
#' RHmax=84
#' RHmin=63
#' #perform simulation
#' modsebs=sebs(albedo=albedo,Ts=Ts,Tmax=Tmax,Tmin=Tmin,RHmax=RHmax,
#' RHmin=RHmin,NDVI=NDVI,SAVI=NULL,iter.max=7,xyhot="full",
#' xycold="full",DOY=37,sunelev=50.71154048,welev=317.1,zx=10,
#' u=2,zomw=2,zom=NULL,LAI=LAI,Rn24=NULL,model="SEBS")
#'
#' 
#' #use original landsat 7 data by specifying folder path
#' folder=system.file("extdata","stack",package="sebkc")
#' modauto=sebs(folder = folder,welev = 380,Tmax=Tmax,
#' Tmin=Tmin,RHmax=RHmax,RHmin=RHmin)
#' #plot ET fraction
#' plot(modauto$EF)
#' 
#' #get Ts, Albedo, NDVI, SAVI, sunelev, DOY from landsat data 
#' welev=278
#' data=landsat578(rawdata,welev=welev)
#'  #perform semi-auto simulation 
#' #Determine xyhot. Digitize polygon on the Ts map
#' modhot=hotTs(data,welev=300,extent="auto",cluster=2)
#' #determine the cold. Digitize polygon on the Ts map
#' modcold=coldTs(data,welev=275,extent="auto",cluster=2)
#' xyhot=modhot$xyhot
#' xycold=modcold$xycold
#' #use object of \code{\link{coldTs}} or \code{\link{hotTs}} in different ways
#' modsebs=sebs(data,xyhot=modhot,xycold=xycold,welev=welev,
#' Tmax=Tmax,Tmin=Tmin,RHmax=RHmax,RHmin=RHmin)
#'  }
#' @export
#' @rdname sebs 

sebs=function(albedo,Ts,Tmax,Tmin,RHmax=NULL,RHmin=NULL,n=NULL,
      NDVI,SAVI,iter.max=7,xyhot,xycold,DOY,
      sunelev,welev,zx,u,zomw,zom=NULL,LAI=NULL,DEM=NULL,
      lapse=0.0065,Rn24=NULL,wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",
      latitude=NULL,model="SEBS",clip=NULL,folder=NULL) UseMethod ("sebs")
#' @export
#' @rdname sebs 
sebs.default=function(albedo=NULL,Ts,Tmax,Tmin,RHmax=NULL,RHmin=NULL,n=NULL,
NDVI,SAVI,iter.max=7,xyhot="auto",xycold="auto",DOY=NULL,
sunelev,welev,zx=10,u=2,zomw=2,zom=NULL,LAI=NULL,DEM=NULL,lapse=0.0065,Rn24=NULL,
wmo=NULL, airport=NULL,Krs=0.16,surface = "grass",latitude=NULL,
model="SEBS",clip=NULL,folder=NULL){
  
  if(is.null(folder)){
    folder="/nothing454782ghf7poi8r.hope"
  }
  
  if(inherits(folder, "landsat578")||inherits(albedo, "landsat578")){
    file.info2=TRUE
  }else{
    file.info2=file.info(folder)[["isdir"]]
  }
  if(file.info2==TRUE&&!is.na(file.info2)){
    if(is.null(welev)){
      return(message("Please provide the parameter welev: 
                   the elavation of th weather station"))  
    }
    
    if(inherits(folder, "landsat578")||inherits(albedo, "landsat578")){
      if(inherits(albedo, "landsat578")){
        mod=albedo
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
    }
  PET=NULL

    if(!is.null(wmo)||!is.null(airport)){
      mod=sebkc.tryCatch(mod)$value
      thisDOY=sebkc.tryCatch(mod$date)$value
      if(inherits(thisDOY, "simpleError")){
        thisDOY=DOY 
      }
      ETr24=sebkc.tryCatch(ETo(DOY=thisDOY,airport = airport,wmo=wmo,latitude=latitude,
                               z=zx,Krs =Krs,altitude=welev))$value
      if(inherits(ETr24, "simpleError")){
        return (message(paste("Invalid DOY [YYYY-mm-dd] or airport or wmo",ETr24))) 
      }
      
      Tmin=ETr24$data$Tmin
      Tmax=ETr24$data$Tmax
      Rn24=ETr24$Rn
      PET=ETr24

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

if(!is.numeric(DOY)&&!is.null(DOY)){
  DOY2=DOY
  DOY=strptime(DOY,"%d-%m-%Y")$yday+1
  if(is.na(DOY)){
    DOY=strptime(DOY2,"%Y-%m-%d")$yday+1
  }
}

mod=sebal(albedo=albedo,Ts=Ts,NDVI=NDVI,SAVI=SAVI,
iter.max=iter.max,xyhot=xyhot,xycold=xycold,
DOY=DOY,sunelev=sunelev,welev=welev,
zx=zx,u=uw,zomw=zomw,zom=zom,LAI=LAI,DEM=DEM,
lapse=lapse,Rn24=Rn24,model="SEBAL",folder=folder)

if(is.null(PET)){
PET=ETo(Tmax=Tmax,Tmin=Tmin,z=zx,uz=u,altitude=welev,
RHmax=RHmax,RHmin=RHmin,n=n,DOY=mod$DOY,latitude=1)
}
y=PET$y
vpd=PET$vpd
slope=PET$slope
# ---------------------------------------------------------------------------
# Native SEBS sensible heat (Su, 2002). The earlier version took SEBAL's
# anchor-calibrated H and rescaled it, which (a) made SEBS depend on the chosen
# hot/cold pixels - the domain sensitivity SEBS is designed to remove - and
# (b) is algebraically redundant: any relative-evaporation rescaling of a given
# H returns LE/(Rn-G), i.e. exactly the SEBAL result. A genuine SEBS therefore
# has to compute its OWN H, which is done here: per pixel, from the actual
# surface-to-air temperature difference and a roughness length for heat derived
# from Su's kB^-1 excess-resistance parameterisation, with a Monin-Obukhov
# stability iteration. Only Rn, G and air density are reused from sebal(); these
# are independent of the anchor pixels, so the SEBS result does not depend on the
# spatial domain. H is then bounded by the physical dry and wet limits.
# NOTE: correctly formulated but pending eddy-covariance validation before
# quantitative use (see NEWS.md).
# ---------------------------------------------------------------------------
warning("sebs(): the native SEBS sensible-heat implementation (Su 2002 kB^-1) is ",
        "EXPERIMENTAL. It is domain-independent by construction but its kB^-1 excess ",
        "resistance has not yet been calibrated/validated against measured sensible ",
        "heat, so absolute EF/ET values are indicative only. Prefer sebal()/METRIC for ",
        "quantitative work. See NEWS.md.", call. = FALSE)
message("computing SEBS sensible heat (Su 2002)")
tmpl <- mod$Rn
Rnv <- getValues(mod$Rn); Gv <- getValues(mod$G); pairv <- getValues(mod$pair)
Tsv <- getValues(Ts)
if(max(Tsv, na.rm=TRUE) < 200) Tsv <- Tsv + 273.15                #ensure Kelvin
if(is.null(LAI)){ LAI <- (5.434*SAVI)+1.5416 }
LAIv <- getValues(if(inherits(LAI,"RasterLayer")) LAI else raster(LAI))
LAIv[is.na(LAIv) | LAIv < 0.05] <- 0.05

nu <- 1.5e-5; Pr <- 0.71; Cd_leaf <- 0.2; Ct_leaf <- 0.01; zref <- zx
Ta_K   <- max(Tmax) + 273.15                                      #overpass (near-midday) air temperature [K]
zomv   <- pmax(0.018*LAIv, 1e-3)                                  #momentum roughness [m]
fc     <- 1-exp(-0.5*LAIv); fs <- 1-fc                            #canopy / soil fractions
hcan   <- pmax(zomv/0.136, 0.02)                                  #canopy height [m]
d0     <- 0.667*hcan                                              #displacement height [m]
uzref  <- u*log(zref/zomw)/log(zx/zomw)                           #wind at blending height [m/s]
ustar_uh <- k/log(pmax(0.333*hcan/zomv, 1.2))                     #within-canopy wind ratio u*/u(h)
nec    <- Cd_leaf*LAIv/(2*ustar_uh^2)                             #canopy wind extinction
ustar  <- k*uzref/log((zref-d0)/zomv)                             #neutral friction velocity
Re_s   <- pmax(0.009*ustar/nu, 0.1)                               #soil roughness Reynolds number
kBs    <- 2.46*(Re_s^0.25) - log(7.4)                             #Brutsaert soil kB^-1
Ct_s   <- Pr^(-2/3)*Re_s^(-0.5)
kB1    <- (k*Cd_leaf)/(4*Ct_leaf*ustar_uh*(1-exp(-nec/2)))*fc^2 +
          (k*ustar_uh*(zomv/hcan))/Ct_s*fc*fs + kBs*fs^2          #Su (2002) kB^-1
kB1[!is.finite(kB1)] <- 2.3; kB1 <- pmin(pmax(kB1, 0), 12)
z0h    <- pmax(zomv/exp(kB1), 1e-6)                               #roughness length for heat [m]
theta0 <- Tsv                                                     #surface potential temperature ~ Ts
thetaa <- Ta_K + 0.0098*zref                                     #air pot. temperature at blending height
dtheta <- theta0 - thetaa
psi_m <- function(z,L){ x <- suppressWarnings((1-16*z/L)^0.25)
  ifelse(L<0, 2*log((1+x)/2)+log((1+x^2)/2)-2*atan(x)+pi/2, ifelse(L>0, -5*z/L, 0)) }
psi_h <- function(z,L){ x <- suppressWarnings((1-16*z/L)^0.25)
  ifelse(L<0, 2*log((1+x^2)/2), ifelse(L>0, -5*z/L, 0)) }
L  <- rep(1e6, length(Rnv))                                      #neutral first guess
Hv <- rep(0, length(Rnv)); rah <- rep(50, length(Rnv))
for(it in 1:8){                                                  #Monin-Obukhov stability iteration
  ustar <- k*uzref/(log((zref-d0)/zomv) - psi_m(zref-d0,L) + psi_m(zomv,L))
  rah   <- (log((zref-d0)/z0h) - psi_h(zref-d0,L) + psi_h(z0h,L))/(k*ustar)
  rah[!is.finite(rah) | rah < 1] <- 1
  Hv    <- pairv*cp*dtheta/rah
  Hv[!is.finite(Hv)] <- 0
  L     <- -(pairv*cp*(ustar^3)*theta0)/(k*g*ifelse(abs(Hv)<1e-3, 1e-3, Hv))
}
message("computing EF")
avail <- Rnv-Gv
Hwet  <- (avail-(pairv*cp/rah)*(vpd/y))/(1+slope/y)              #wet limit (Penman-Monteith)
Hv    <- pmin(pmax(Hv, Hwet), avail)                            #bound H to [Hwet, Hdry=Rn-G]
EFv   <- 1.1*((avail-Hv)/avail)                                 #EF = 1.1 * LE/(Rn-G)
EFv[avail<=0 | !is.finite(EFv) | EFv<0] <- 0; EFv[EFv>1.1] <- 1.1
EF <- setValues(tmpl, EFv)
H  <- setValues(tmpl, Hv)

ET24=NULL
#print(Rn24)
if(!is.null(Rn24)){
  message("computing ETa" )
  ET24=EF*Rn24*0.035
}
factor<-list(EF=EF,H=H,ET24=ET24,folder=folder)
factor$call<-match.call()
class(factor)<-"sebs"
factor
}




