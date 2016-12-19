#' @title S-SEBI Evapotranspiration Fraction
#' 
#' @description \code{ssebi} computes Roerink et. al (2000) S-SEBI: A simple
#'  remote sensing algorithm to estimate the surface energy balance.
#' @inheritParams sebal
#' @inheritParams coldTs
#' @param threshold numeric. Albedo inflection point threshold
#' @author George Owusu
#'
#' @return EF, amax,bmax,amin,bmin
#' @references 
#' Roerink, G. J., Su, Z., & Menenti, M. 2000. S-SEBI: 
#' A simple remote sensing algorithm to estimate the surface energy balance.
#' Physics and Chemistry of the Earth, Part B: Hydrology, Oceans and Atmosphere
#' , 25(2): 147-157.
#' @examples
#' \dontrun{
#' #Manual data specification
#' albedo=raster(system.file("extdata","albedo.grd",package="SEBKc"))
#' Ts=raster(system.file("extdata","Ts.grd",package="SEBKc"))
#' mod=ssebi(Ts=Ts,albedo=albedo,cluster=100,threshold=0,plot=TRUE)
#' 
#' #Using landsat folder
#' folder=system.file("extdata","stack",package="SEBKc")
#' ssebiauto=ssebi(folder=folder,welev=278)
#' }
#' @export
#' @rdname ssebi 
ssebi=function(Ts,albedo=NULL,cluster=100,
               threshold=0.2,plot=TRUE,Rn24=NULL,folder=NULL,welev=NULL)
  UseMethod ("ssebi")
#' @export
#' @rdname ssebi 
  
ssebi.default=function(Ts=NULL,albedo=NULL,cluster=100,
            threshold=0.2,plot=TRUE,Rn24=NULL,folder=NULL,welev=NULL){
  
  
  if(is.null(folder)){
    folder="/nothing454782ghf7poi8r.hope"
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
      mod=landsat578(folder=folder, welev=welev)
    }
    #mod=landsat578(data=folder, welev=welev)
    albedo=mod$albedo
    Ts=mod$Ts
    NDVI=mod$NDVI
    DOY=mod$DOY
    sunelev=mod$sunelev
    SAVI=mod$SAVI
    
    }
  
if(class(Ts)!="RasterLayer"){
Ts=raster(Ts)
}
if(class(albedo)!="RasterLayer"&&!is.null(albedo)){
albedo=raster(albedo)
}
max=maxValue(albedo)
min=minValue(albedo)
range=max-min
ave=range/cluster
seq=(seq(min,max,ave))
albedoclass=cut(albedo,breaks=seq)
albedomin=zonal(albedo, albedoclass, fun='min') 
Tsmin=zonal(Ts, albedoclass, fun='min')
Tsmax=zonal(Ts, albedoclass, fun='max') 
albedomax=zonal(albedo, albedoclass, fun='max') 

alTsmindata=merge(albedomin,Tsmin,by="zone")
alTsmaxdata=merge(albedomax,Tsmax,by="zone")
alTsmindata2<- alTsmindata[ which(alTsmindata[2]>0), ]
alTsmaxdata2<- alTsmaxdata[ which(alTsmaxdata[2]>0), ]
if(threshold=="auto"){
inflection=alTsmaxdata2[ which(alTsmaxdata2[3]==max(alTsmaxdata2[[3]])), ]
alTsmaxdata2<- alTsmaxdata2[ which(alTsmaxdata2[1]>=inflection[[1]]), ]
}
if(is.numeric(threshold)){
alTsmaxdata2<- alTsmaxdata2[ which(alTsmaxdata2[2]>threshold), ]
}
modTmax=lm(alTsmaxdata2[[3]]~alTsmaxdata2[[2]])
modTmin=lm(alTsmindata2[[3]]~alTsmindata2[[2]])
amax=coef(modTmax)[[1]]
bmax=coef(modTmax)[[2]]
amin=coef(modTmin)[[1]]
bmin=coef(modTmin)[[2]]
EF=(((bmax*albedo)+amax)-Ts)/(((bmax*albedo)+amax)-((bmin*albedo)+amin))
EF[EF<0]=NA
EF[EF>1]=NA

if(plot==T||plot==TRUE){
#plot(albedo,Ts, maxpixels=5000000,xlab="Albedo",ylab="Surface Tempearture[Ts]")
plot(albedo,Ts,maxpixels=ncell(albedo)/2,xlab="Albedo",ylab="Surface Tempearture[Ts]")

#lines(alTsmindata2[[2]],alTsmindata2[[3]],type="p",col=2)
#lines(alTsmaxdata3[[2]],alTsmaxdata3[[3]],type="p",col=3)
abline(modTmin, col = "forestgreen")
abline(modTmax, col = "red")
dev.new(width=15, height=10)

plot(EF)
}
ET24=NULL
if(!is.null(Rn24)){
  ET24=(EF*Rn24)/2.45
}

factor<-list(EF=EF,amax=amax,bmax=bmax,amin=amin,bmin=bmin,folder=folder)
factor$call<-match.call()
class(factor)<-"ssebi"
factor
}



