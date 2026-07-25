## =====================================================================
## sebkc manual - reproduction script
## Runs every workflow in the manual on the BUNDLED data and captures the
## numbers (-> manual/repro_output.txt) and figures (-> manual/figures/*.png),
## so the manual can be written from real output, not guesses.
##
## Run:  source("C:/Users/PC/Documents/GitHub/sebkc/manual/reproduce_manual.R")
## Needs: sebkc (>= 1.0-6, with the modern-R fixes), raster, sp.
## =====================================================================
OUT <- "C:/Users/PC/Documents/GitHub/sebkc/manual"
FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
suppressWarnings(suppressMessages(library(sebkc)))
suppressWarnings(suppressMessages(library(raster)))

figpng <- function(name, expr, w = 900, h = 650) {
  f <- file.path(FIG, paste0(name, ".png")); cat("[FIG]", name, ": ")
  png(f, width = w, height = h, res = 110)
  r <- try(force(expr), silent = TRUE); try(dev.off(), silent = TRUE)
  cat(if (inherits(r, "try-error")) paste("FAIL", conditionMessage(attr(r,"condition"))) else "ok", "\n")
  invisible(r)
}
st <- function(x) if (is.null(x)) "NULL" else sprintf("min %.3f  mean %.3f  max %.3f",
      cellStats(x,"min"), cellStats(x,"mean"), cellStats(x,"max"))
sec <- function(s) cat("\n\n==========", s, "==========\n")

sink(file.path(OUT, "repro_output.txt"), split = TRUE)
cat("sebkc manual reproduction\n")
cat("sebkc", as.character(packageVersion("sebkc")), "|", R.version.string, "\n")

sf <- function(...) system.file("extdata", ..., package = "sebkc")

## ---- 1. Reference ET (FAO-56) ---------------------------------------
sec("1. Reference ET: ETo")
PETmin <- ETo(Tmax = 31, Tmin = 28, latitude = 5.6, DOY = "6-2-2001")
cat("ETo minimal (Accra-like):", round(PETmin$ETo, 3), "mm/day\n")
PET18 <- ETo(Tmax = 21.5, Tmin = 12.3, z = 10, uz = 2.78, altitude = 100,
             RHmax = 84, RHmin = 63, n = 9.25, DOY = 187, latitude = 50.8)
cat("ETo FAO-56 Example 18 (Brussels):", round(PET18$ETo, 3), "mm/day (textbook 3.9)\n")
cat("  Rn:", round(PET18$Rn_mj_m2_day,3), "MJ/m2/day  slope:", round(PET18$slope,4),
    "  psychrometric y:", round(PET18$y,4), "  vpd:", round(PET18$vpd,3), "\n")
sec("1b. Hourly ET: ETohr")
ET1 <- ETohr(Tmean=38,RH=52,DOY=274,Lz=15,t1=1,time=14.5,Lm=16.25,latitude=16.22,uz=3.3,Rs=2.450)
cat("ETohr (FAO-56 Ex19-like):", round(ET1$ETo,3), "mm/hour\n")

## ---- 2. Crop coefficient / FAO-56 water balance ---------------------
sec("2. Crop coefficient (kc): FAO-56 textbook examples")
Ex32 <- kc(ETo=7,P=0,RHmin=20,soil="sandy loam",crop="Broccoli",I=10,kctype="dual",u2=3,kc=c(0.9,0,0),h=1,Zr=0.1)
cat("Ex32 sprinkler: fc =", round(Ex32$output$fc,3), " Kc_adj =", round(Ex32$output$Kc_adj,3), "\n")
Ex33 <- kc(ETo=7,P=0,RHmin=20,soil="sandy loam",crop="Broccoli",I=10,kctype="dual",u2=3,kc=c(0.9,0,0),h=1,Zr=0.1,fw=0.3)
cat("Ex33 furrow (fw=0.3): Kc_adj =", round(Ex33$output$Kc_adj,3), "\n")
Ex34 <- kc(ETo=7,P=0,RHmin=20,soil="sandy loam",crop="Broccoli",I=10,kctype="dual",u2=3,kc=c(0.9,0,0),h=1,Zr=0.1,fw="drip")
cat("Ex34 drip: Kc_adj =", round(Ex34$output$Kc_adj,3), "\n")

sec("2b. Multi-day water balance from bundled irrigation.txt")
d <- read.table(sf("sys","irrigation.txt"), header = TRUE)
cat("input columns:", paste(names(d), collapse=", "), "  days:", nrow(d), "\n")
msingle <- kc(ETo=d$ETo,P=d$P,RHmin=35,soil="sandy loam",crop="Broccoli",I=d$I)
cat("single-Kc output columns:", paste(names(msingle$output), collapse=","), "\n")
print(round(msingle$output[c("Day","ETo","Kc","ETc","TAW","RAW","AW_modelled","Runoff")], 2))
mdual <- kc(ETo=d$ETo,P=d$P,RHmin=35,soil="sandy loam",crop="Broccoli",I=d$I,kctype="dual",
            FC=0.23,p=0.6,WP=0.10,u2=1.6,kc=c(0.3,0,0),Zr=d$Zr)
cat("\ndual-Kc mean ETc:", round(mean(mdual$output$ETc),3), "mm/day\n")
## classic FAO-56 dual-Kc season curve using the bundled Kcb column + stage lengths
mcurve <- kc(ETo=d$ETo,P=d$P,RHmin=35,soil="sandy loam",crop="Broccoli",I=d$I,kctype="dual",
  FC=0.23,p=0.6,WP=0.10,u2=1.6,rc=0,kc=d$Kcb,Zr=d$Zr,lengths=c(4,4,2,2),fw=0.8)
oc <- mcurve$output
cat("Kc_adj range:", paste(round(range(oc$Kc_adj),3),collapse="-"),
    " Kcb range:", paste(round(range(oc$Kcb),3),collapse="-"), "\n")
figpng("fig_kc_curve", { plot(oc$Day,oc$Kc_adj,type="o",pch=19,col="#1f7a8c",ylim=c(0,1.2),
  xlab="Day of season", ylab="crop coefficient", main="Dual crop coefficient over the season (FAO-56)")
  lines(oc$Day,oc$Kcb,type="o",pch=1,lty=2,col="#c1462f")
  legend("topleft", c("Kc (Kcb+Ke, adjusted)","Kcb (basal)"), col=c("#1f7a8c","#c1462f"),
         pch=c(19,1), lty=c(1,2), bty="n") })

## ---- 3. Landsat preprocessing ---------------------------------------
sec("3. Landsat 7 radiometric processing: landsat578")
l7 <- brick(sf("landsat7_6feb2004.grd"))
L <- landsat578(data=l7, welev=317.1, sensor=7, gain="high", sunelev=50.71154048,
                K1=666.09, K2=1282.71, date="2002-2-6")
cat("albedo :", st(L$albedo), "\n")
cat("NDVI   :", st(L$NDVI), "\n")
cat("SAVI   :", st(L$SAVI), "\n")
cat("Ts (K) :", st(L$Ts), "\n")
figpng("fig_albedo", plot(L$albedo, main="Surface albedo (Landsat 7, 6 Feb 2004)"))
figpng("fig_NDVI",   plot(L$NDVI,   main="NDVI (Landsat 7)"))
figpng("fig_Ts",     plot(L$Ts,     main="Surface temperature Ts [K]"))

## ---- 4. Energy-balance inputs + anchor pixels -----------------------
sec("4. Anchor pixels: coldTs / hotTs")
albedo <- raster(sf("albedo.grd")); Ts <- raster(sf("Ts.grd"))
NDVI <- raster(sf("NDVI.grd")); LAI <- raster(sf("LAI.grd"))
cat("bundled Ts [K]:", st(Ts), "\n")
set.seed(1); cold <- coldTs(Ts=Ts,NDVI=NDVI,albedo=albedo,sunangle=50,cluster=7,extent="full",plot=FALSE,seed=1)
cat("cold pixel: x=",round(cold$x),"y=",round(cold$y)," Tscold=",round(cold$Tscold,2),"K\n")
set.seed(1); hot <- hotTs(Ts=Ts,NDVI=NDVI,albedo=albedo,cluster=7,extent="full",plot=FALSE,seed=1)
cat("hot  pixel: x=",round(hot$x),"y=",round(hot$y)," Tshot=",round(hot$Tshot,2),"K\n")

## ---- 5. Surface energy balance models -------------------------------
sec("5. One-source models: SEBAL / SEBS / SEBI / SSEB / S-SEBI")
runmod <- function(tag, m) {
  if (inherits(m, "try-error")) { cat(sprintf("%-8s FAIL %s\n", tag, conditionMessage(attr(m,"condition")))); return(NULL) }
  ef <- if (!is.null(m$EF)) m$EF else m$ETrF
  cat(sprintf("%-8s EF: %s\n", tag, st(ef)))
  if (!is.null(m$ET24)) cat(sprintf("         ET24 [mm/day]: %s\n", st(m$ET24)))
  figpng(paste0("fig_EF_", tag), plot(ef, main=paste(tag, "evaporative fraction")))
  invisible(m)
}
set.seed(1); mSEBAL <- try(sebal(albedo=albedo,Ts=Ts,NDVI=NDVI,SAVI=NULL,iter.max=7,xyhot="full",
  xycold="full",DOY=37,sunelev=50.71154048,welev=317.1,zx=10,u=2,zomw=3,LAI=LAI,model="SEBAL"), silent=TRUE)
runmod("SEBAL", mSEBAL)
if (!inherits(mSEBAL,"try-error")) { cat("  SEBAL Rn:",st(mSEBAL$Rn),"\n  SEBAL G:",st(mSEBAL$G),
  "\n  SEBAL H:",st(mSEBAL$H),"\n  SEBAL LE:",st(mSEBAL$LE),"\n") }
set.seed(1); mSEBS <- try(sebs(albedo=albedo,Ts=Ts,Tmax=31,Tmin=28,RHmax=84,RHmin=63,NDVI=NDVI,SAVI=NULL,
  iter.max=7,xyhot="full",xycold="full",DOY=37,sunelev=50.71,welev=317.1,zx=10,u=2,zomw=2,LAI=LAI,model="SEBS"), silent=TRUE)
runmod("SEBS", mSEBS)
set.seed(1); mSEBI <- try(sebi(albedo=albedo,Ts=Ts,Tmax=31,Tmin=28,RHmax=84,RHmin=63,NDVI=NDVI,SAVI=NULL,
  iter.max=7,xyhot="full",xycold="full",DOY=37,sunelev=50.71,welev=317.1,zx=10,u=2,zomw=2,LAI=LAI,model="sebi"), silent=TRUE)
runmod("SEBI", mSEBI)
set.seed(1); mSSEB <- try(sseb(Ts=Ts,TH="full",TC="full",sunelev=50,NDVI=NDVI,albedo=albedo,ETo="auto",
  x=1.2,Tmax=31,Tmin=28,zx=10,u=2,DOY=37,latitude=5.6), silent=TRUE)
runmod("SSEB", mSSEB)
set.seed(1); mSSEBI <- try(ssebi(Ts=Ts,albedo=albedo,threshold=0,plot=FALSE), silent=TRUE)
runmod("S-SEBI", mSSEBI)

sec("5b. Two-source model: TSEB (parallel and series)")
set.seed(1); mTSp <- try(tseb(Ts=Ts,LAI=LAI,DOY=37,xyhot="full",xycold="full",albedo=albedo,Tmax=31,Tmin=28,
  NDVI=NDVI,SAVI=NULL,hc=20,sunelev=50,welev=345,u=2,zx=200,zomw=2,xPT=1.3,network="parallel",latitude=5.6,n=6.5), silent=TRUE)
if (!inherits(mTSp,"try-error")) { cat("TSEB-parallel  ET24:",st(mTSp$ET24)," E24:",st(mTSp$E24)," T24:",st(mTSp$T24),"\n")
  figpng("fig_TSEB_parallel_ET24", plot(mTSp$ET24, main="TSEB (parallel) ET24 [mm/day]")) } else cat("TSEB parallel FAIL\n")
set.seed(1); mTSs <- try(tseb(Ts=Ts,LAI=LAI,DOY=37,xyhot="full",xycold="full",albedo=albedo,Tmax=31,Tmin=28,
  NDVI=NDVI,SAVI=NULL,hc=20,sunelev=50,welev=345,u=2,zx=200,zomw=2,xPT=1.3,network="series",latitude=5.6,n=6.5), silent=TRUE)
if (!inherits(mTSs,"try-error")) cat("TSEB-series    ET24:",st(mTSs$ET24),"\n") else cat("TSEB series FAIL\n")

## ---- 6. Biomass and Water Deficit Index -----------------------------
sec("6. Biomass (point) and WDI")
bm <- try(biomass(Tday=26.9,T24=25.3,Rs=452,latitude=11.18,HI=0.4,plantdate="2015-05-1",
              harvestdate="2015-08-29",c30=0.0108,adaptability=3,LAI=5,crop="maize"), silent=TRUE)
if (inherits(bm,"try-error")) { cat("biomass FAIL:", conditionMessage(attr(bm,"condition")),"\n")
} else { cat("biomass return names:", paste(names(bm), collapse=","), "\n")
  y <- tryCatch(bm$yield, error=function(e) NULL); if (is.null(y)) y <- tryCatch(bm$output, error=function(e) NULL)
  cat("biomass yield/output:\n"); print(y) }
set.seed(1); mw <- try(wdi(Ts=Ts,Ta=299,NDVI=NDVI,Tsmax=313,Tsmin=304), silent=TRUE)
if (!inherits(mw,"try-error")) { cat("WDI EF:", st(mw$EF), "\n")
  figpng("fig_WDI", plot(mw$EF, main="Water Deficit Index (EF)")) } else cat("WDI FAIL:", conditionMessage(attr(mw,"condition")),"\n")

## ---- 7. Model comparison of mean EF ---------------------------------
sec("7. Mean evaporative fraction across models (same scene)")
efmean <- function(m) if (inherits(m,"try-error")||is.null(m)) NA else {
  ef <- if(!is.null(m$EF)) m$EF else m$ETrF; round(cellStats(ef,"mean"),3) }
comp <- data.frame(model=c("SEBAL","SEBS","SEBI","SSEB","S-SEBI"),
  mean_EF=c(efmean(mSEBAL),efmean(mSEBS),efmean(mSEBI),efmean(mSSEB),efmean(mSSEBI)))
print(comp)

sec("8. sessionInfo")
print(sessionInfo())
cat("\n\nDONE\n"); sink()
cat("\nReproduction complete. Text -> repro_output.txt ; figures ->", FIG, "\n")
