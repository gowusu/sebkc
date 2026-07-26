## reproduce_water_needs.R
## Regenerates every number and figure in the crop-water-needs manual
## ("How Much Water Does My Crop Need?"). Run it and you get exactly the
## tables and plots printed there.
##
##   "C:/Program Files/R/R-4.6.0/bin/Rscript.exe" reproduce_water_needs.R
##
## Produced on R 4.6.0 with sebkc >= 1.0-8.
suppressWarnings(suppressMessages(library(sebkc)))
set.seed(1)
here <- "C:/Users/PC/Documents/GitHub/sebkc/manual"
figdir <- file.path(here, "figures"); dir.create(figdir, showWarnings = FALSE)
outfile <- file.path(here, "repro_water_needs_output.txt")
con <- file(outfile, "w")
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = con) }
sec <- function(t) say("\n========== ", t, " ==========")

## --- Kumasi weather: already updated to the present by weather() -------
lat <- 6.72; lon <- -1.60; alt <- 222
kum <- readRDS(system.file("extdata","kc","kumasi1960_present.rds", package="sebkc"))
kum$DOY <- as.character(kum$DOY); kum$Date <- as.Date(kum$DOY)
RHm <- round(mean(kum$RHmin, na.rm = TRUE), 1)

sec("1. The weather record we are standing on")
say("rows:", nrow(kum), " span:", as.character(min(kum$Date)), "->",
    as.character(max(kum$Date)))
say("mean RHmin used for kc:", RHm, "%")

## --- Climate normals (last 10 full years) for context -----------------
rec <- kum[kum$Date >= as.Date("2015-01-01") & kum$Date <= as.Date("2024-12-31"), ]
rec$mon <- as.integer(format(rec$Date, "%m")); rec$yr <- format(rec$Date, "%Y")
etoDay <- function(w) sapply(seq_len(nrow(w)), function(k)
  suppressMessages(ETo(Tmax=w$Tmax[k],Tmin=w$Tmin[k],RHmax=w$RHmax[k],RHmin=w$RHmin[k],
      n=w$n[k],uz=w$uz[k],DOY=w$DOY[k],latitude=lat,altitude=alt)$ETo))
rec$ETo <- suppressMessages(etoDay(rec))
mon_eto  <- tapply(rec$ETo, rec$mon, mean)                       # mm/day
mon_rain <- tapply(rec$P, list(rec$mon, rec$yr), sum)            # mm/month by yr
mon_rain <- rowMeans(mon_rain, na.rm = TRUE)                     # mean mm/month
mlab <- month.abb
normals <- data.frame(month = mlab,
                      ETo_mm_day = round(as.numeric(mon_eto), 2),
                      rain_mm = round(as.numeric(mon_rain), 0))
sec("2. Kumasi climate normals 2015-2024 (mean)")
capture.output(print(normals, row.names = FALSE)) |> writeLines(con)
print(normals, row.names = FALSE)

## --- season helper ----------------------------------------------------
season_eto <- function(start, ndays) {
  i0 <- which(kum$Date == as.Date(start)); w <- kum[i0:(i0 + ndays - 1), ]
  data.frame(Date = w$Date, ETo = suppressMessages(etoDay(w)), P = w$P)
}

## --- crop parameters (FAO-56 Table 11 & 12) ---------------------------
crops <- list(
  Maize   = list(kc=c(0.30,1.20,0.35), len=c(20,35,40,30), Zr=1.0, p=0.55, h=2.0),
  Tomato  = list(kc=c(0.60,1.15,0.80), len=c(30,40,40,25), Zr=0.9, p=0.40, h=0.6),
  Cabbage = list(kc=c(0.70,1.05,0.95), len=c(40,60,50,15), Zr=0.6, p=0.45, h=0.4),
  Onion   = list(kc=c(0.70,1.05,0.75), len=c(15,25,70,40), Zr=0.45,p=0.30, h=0.4),
  Pepper  = list(kc=c(0.60,1.05,0.90), len=c(30,35,40,20), Zr=0.7, p=0.30, h=0.7))

runcrop <- function(cr, start, label) {
  nd <- sum(cr$len); s <- season_eto(start, nd)
  sink(tempfile())  # swallow kc's per-day chatter
  m <- suppressMessages(suppressWarnings(kc(ETo=s$ETo, P=s$P, RHmin=RHm,
        soil="sandy loam", crop=label, I=0, kctype="single",
        kc=cr$kc, lengths=cr$len, Zr=cr$Zr, p=cr$p, h=cr$h)))
  sink()
  o <- m$output
  list(days=nd, ETo=sum(s$ETo), ETc=sum(o$ETc,na.rm=TRUE),
       rain=sum(s$P), peak=max(o$ETc,na.rm=TRUE),
       s=s, o=o)
}

seasons <- list(dry = "2024-11-15", main = "2025-03-15")
tab <- data.frame()
detail <- list()
for (sn in names(seasons)) for (cn in names(crops)) {
  r <- runcrop(crops[[cn]], seasons[[sn]], cn)
  detail[[paste(cn,sn)]] <- r
  tab <- rbind(tab, data.frame(crop=cn, season=sn, days=r$days,
    ETo=round(r$ETo), ETc=round(r$ETc), perday=round(r$ETc/r$days,1),
    peak=round(r$peak,1), rain=round(r$rain),
    irrigation=max(0, round(r$ETc - r$rain))))
}
sec("3. Crop water needs near Kumasi (planted dry 2024-11-15 / main 2025-03-15)")
capture.output(print(tab, row.names = FALSE)) |> writeLines(con)
print(tab, row.names = FALSE)

## --- featured walk-through: dry-season Tomato -------------------------
ft <- detail[["Tomato dry"]]
sec("4. Featured example: dry-season Tomato, stage by stage")
o <- ft$o; len <- crops$Tomato$len; cum <- cumsum(len)
stage <- cut(o$Day, breaks = c(0, cum), labels = c("initial","development","mid","late"))
byst <- data.frame(
  stage = levels(stage),
  days  = as.integer(table(stage)),
  ETc_mm = round(tapply(o$ETc, stage, sum)),
  rain_mm = round(tapply(o$P, stage, sum)))
byst$irrigation_mm <- pmax(0, byst$ETc_mm - byst$rain_mm)
capture.output(print(byst, row.names = FALSE)) |> writeLines(con)
print(byst, row.names = FALSE)
say("season crop water need ETc:", round(ft$ETc), "mm over", ft$days, "days")
say("peak daily need:", round(ft$peak,1), "mm/day  (design figure for the pump/drip line)")
say("season rain:", round(ft$rain), "mm  ->  irrigation top-up:",
    max(0, round(ft$ETc - ft$rain)), "mm")

## ---------------- figures ---------------------------------------------
## Fig 1: climate normals - ETo vs rainfall
png(file.path(figdir,"cwn_climate.png"), width=1600, height=950, res=200)
op <- par(mar=c(4,4,3,4))
bp <- barplot(normals$rain_mm, names.arg=normals$month, col="#9ecae1",
      border=NA, ylim=c(0,260), ylab="Rainfall (mm/month)",
      main="Kumasi climate normals (2015-2024): rain vs atmospheric demand")
par(new=TRUE)
plot(bp, normals$ETo_mm_day, type="o", pch=19, lwd=2, col="#d95f02",
     axes=FALSE, xlab="", ylab="", ylim=c(0,6))
axis(4, col="#d95f02", col.axis="#d95f02")
mtext("Reference ET, ETo (mm/day)", side=4, line=2.5, col="#d95f02")
legend("topright", c("Rainfall","ETo (demand)"), pch=c(15,19),
       col=c("#9ecae1","#d95f02"), bty="n")
par(op); dev.off()

## Fig 2: featured tomato - daily crop water use vs rainfall
png(file.path(figdir,"cwn_tomato.png"), width=1600, height=950, res=200)
op <- par(mar=c(4,4,3,2))
plot(ft$o$Day, ft$o$ETc, type="h", col="#d95f02", lwd=2,
     xlab="Days after planting", ylab="Water (mm/day)",
     main="Dry-season Tomato near Kumasi: daily crop water use vs rainfall",
     ylim=c(0, max(ft$o$ETc, ft$s$P)))
points(ft$s$Date - min(ft$s$Date) + 1, ft$s$P, type="h", col="#3182bd", lwd=2)
legend("topleft", c("Crop water use (ETc)","Rainfall"),
       col=c("#d95f02","#3182bd"), lwd=3, bty="n")
par(op); dev.off()

## Fig 3: irrigation requirement by crop and season
png(file.path(figdir,"cwn_irrigation.png"), width=1600, height=950, res=200)
op <- par(mar=c(4,4,3,2))
mat <- rbind(dry  = tab$irrigation[tab$season=="dry"],
             main = tab$irrigation[tab$season=="main"])
colnames(mat) <- tab$crop[tab$season=="dry"]
barplot(mat, beside=TRUE, col=c("#e6550d","#31a354"), border=NA,
        ylab="Irrigation top-up (mm/season)",
        main="Irrigation need by crop: dry-season vs main-season planting")
legend("topright", c("Dry-season planting","Main-season planting"),
       fill=c("#e6550d","#31a354"), border=NA, bty="n")
par(op); dev.off()

say("\nfigures written to", figdir)
close(con)
cat("\nDONE. Console log saved to", outfile, "\n")
