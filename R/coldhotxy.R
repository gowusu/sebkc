#' @title  Automatic computation of hot and cold pixels from Radiometric Surface Temperature (Ts)

#' @description This function computes temperature coordinates of a hot and cold pixels based on only
#' Surface Temperature. This function uses ranges of temperature.
#' @param Tmax Numeric. Maximum surface Temperature in map units (eg. Kilvin) that you want computation from. 
#' You can keep changing untill you get your desired values
#' @param Tmin Numeric. Minimum surface Temperature in  map units (eg. Kilvin) that you want computation from. 
#' You can keep changing untill you get your desired values
#' @inheritParams sebal
#' @inheritParams coldTs2
#' @inheritParams coldTs
#' @inheritParams landsat578
#' @details The function first divides the selected area into a 
#' number of clusters. The mean and standard deviations are
#' computed. The cluster means with the lowest standard deviations
#' and qualifying for upper (Ts) and lower (NDVI) limits 
#' (see above) are selected. The provided upper and lower probabilities
#'  are also used to select the final Ts.
#' @seealso \code{\link{coldTs}},  \code{\link{hotTs}}, \code{\link{coldT2}} 
#' @return 
#' \itemize{
#' \item{Ts:} { The selected or input Ts}
#' \item{Tshot:} { The Temperature of hot pixel}
#' \item{Tscold:} { The Temperature of cold pixel}
#' \item{xyhot: } {combined x and y coordinates}
#' \item{candidates:}{ The similar candidates' pixels that can be used}
#' }
#' @examples
#' \dontrun{
#'  folder=system.file("extdata","stack",package="sebkc")
#'   model=coldhot(Tmin=307,Tmax=307,folder=folder,welev=170)
#'   spplot(model$Tshot)
#'   spplot(model$Tscold)
#'   modauto=sebal(folder = folder,welev = 380,xycold=model$xycold,xyhot=model$xyhot)}
#' @author George Owusu
#' @references Allen, R. G., Burnett, B., Kramber, W., Huntington, J., 
#' Kjaersgaard, J., Kilic, A., Kelly, C., & Trezza, R. 2013. Automated 
#' Calibration of the METRIC-Landsat Evapotranspiration Process. 
#' JAWRA Journal of the American Water Resources Association, 49(3): 563-576.
#' @export
#' @rdname coldhot
#'

coldhot<-function(Tmin,Tmax, Ts=NULL,folder=NULL,welev=NULL ) UseMethod ("coldhot")
#' @export
#' @rdname coldhot
coldhot.default=function(Tmin,Tmax,Ts=NULL,folder=NULL,welev=NULL){
  
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
    
    
    Ts=mod$Ts

    }
  
 
  
  #get the hot atemperature values
  temphot=Ts
  #set a threshold for hot temperatures
  temphot[temphot<=Tmax]=NA
  phot=data.frame(rasterToPoints(temphot))
  
  #get the cold atemperature values
  tempcold=Ts
  #set a threshold for cold temperatures
  tempcold[tempcold>Tmin]=NA
  pcold=data.frame(rasterToPoints(tempcold))
  xyhot=c(phot[1,][[1]], phot[1,][[2]])
  xycold=c(pcold[1,][[1]], pcold[1,][[2]])
  
 
  
  factor<-list(xyhot=xyhot,xycold=xycold,Ts=Ts,Tshot=temphot,Tscold=tempcold)
  
  factor$call<-match.call()
  
  class(factor)<-"coldhot"
  factor
  }

#folder=system.file("extdata","stack",package="sebkc")
#model=coldhot(Tmin=307,Tmax=310,folder=folder,welev=170)
#spplot(model$Tshot)
#spplot(model$Tscold)
