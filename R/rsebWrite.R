#' @title Writing sebkc object to a folder
#' 
#' @description This function writes spatial output of the following functions
#' \code{\link{sebal}}, \code{\link{sebi}},\code{\link{ssebi}},
#' \code{\link{wdi}},\code{\link{sebs}},\code{\link{sseb}},
#' \code{\link{tseb}} and \code{\link{landsat578}}  
#'
#' @param object sebkc object 
#' @param folder Folder the files should be written to.If it is set to NULL,
#' it is written to the input folder of sebkc object.
#' @param xy A dataframe of of xy coordinates in in decimal degrees or 
#' meters in the order of c(x,y). If it is not set NULL,
#'  the corresponding values are extracted and written to the folder
#' @param overwrite logical whether the file should be over written
#' @return Writes output to a folder
#' @export
#' @examples
#' \dontrun{
#' folder=system.file("extdata","stack",package="SEBKc")
#' stack=landsat578(data=folder, welev=362)
#' writesebkc(stack)
#' sebaloutput=sebal(folder = folder,welev = 317)
#' writesebkc(sebaloutput)
#' 
#' }
writesebkc<-function(object,folder=NULL,xy=NULL,overwrite=TRUE){
  #longitude=c(-1.62,-1.61,-1.59,-1.59,-1.61,-1.62,-1.58,-1.57)
  #latitude=c(6.7,6.5,6.3,7.0,6.7,6.5,6.3,7.0)
  #xy1 <- data.frame(x=longitude, y=latitude)
  #xy2=NULL
  if(!is.null(xy)){
    xy=data.frame(x=xy[1],y=xy[2])
    if(xy[[1]][1]<200){
      coordinates(xy) <- c("x", "y")
      projection(xy) <- CRS("+init=epsg:4326") # WGS 84
      CRS.new=CRS( proj4string(object[[1]]))
      dtrans <- spTransform(xy, CRS.new)
      xy=dtrans@coords  
    }
    
  }
  
  
  if(is.null(folder)){
    folder=object$folder
  }
  
  if(class(object)=="sebal"){
    writeRaster(object$EF, filename=paste(folder,"EF.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$LE, filename=paste(folder,"LE.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$H, filename=paste(folder,"H.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$Rn, filename=paste(folder,"Rn.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$G, filename=paste(folder,"G.tif",sep="/"), format="GTiff",overwrite=overwrite)
    if(!is.null(object$ET24)){
      writeRaster(object$ET24, filename=paste(folder,"ET24.tif",sep="/"), format="GTiff",overwrite=overwrite)
      if(!is.null(xy)){
      xydata=stack(object$ET24, object$EF,object$LE,object$H,object$Rn,object$G)
      varname=c("ET24","EF","LE","H","Rn","G")
      names(xydata)=varname
      }
    }else{
      if(!is.null(xy)){
      xydata=stack(object$EF,object$LE,object$H,object$Rn,object$G)
      varname=c("EF","LE","H","Rn","G")
      names(xydata)=varname
      }
      
    }
    #xy=modcold$p[1:2]
  }  
  
  if(class(object)=="sebi"||class(object)=="ssebi"||class(object)=="wdi"||
     class(object)=="sebs"||class(object)=="sseb"){
    writeRaster(object$EF, filename=paste(folder,"EF.tif",sep="/"), format="GTiff",overwrite=overwrite)
    if(!is.null(object$ET24)){
      writeRaster(object$ET24, filename=paste(folder,"ET24.tif",sep="/"), format="GTiff",overwrite=overwrite)
    }
    if(!is.null(xy)&&!is.null(object$ET24)){
      xydata=stack(object$ET24, object$EF)
      varname=c("ET24","EF")
      names(xydata)=varname
    }else{
      xydata= object$EF
      varname=c("EF")
      names(xydata)=varname
    }
  }
  
  
  if(class(object)=="tseb"){
    writeRaster(object$EFc, filename=paste(folder,"EFc.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$EFs, filename=paste(folder,"EFs.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$LEc, filename=paste(folder,"LEc.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$LEs, filename=paste(folder,"Les.tif",sep="/"), format="GTiff",overwrite=overwrite)
    if(!is.null(object$ET24)){
      writeRaster(object$T24, filename=paste(folder,"T24.tif",sep="/"), format="GTiff",overwrite=overwrite)
      writeRaster(object$E24, filename=paste(folder,"E24.tif",sep="/"), format="GTiff",overwrite=overwrite)
      writeRaster(object$ET24, filename=paste(folder,"ET24.tif",sep="/"), format="GTiff",overwrite=overwrite)
    }
    
    if(!is.null(xy)&&!is.null(object$ET24)){
      xydata=stack(object$EFc,object$EFs,object$LEc,object$LEs,object$ET24,
                   object$T24, object$E24)
      varname=c("EFcanopy","EFsoil","LEcanopy","LEsoil","ET24","Transpiration24hr", "Evaporation24hr")
      names(xydata)=varname
    }else{
      xydata=stack(object$EFc,object$EFs,object$LEc,object$LEs)
      varname=c("EFcanopy","EFsoil","LEcanopy","LEsoil")
      names(xydata)=varname
    }
  }
  
  if(class(object)=="landsat578"){
    writeRaster(object$albedo, filename=paste(folder,"albedo.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$Ts, filename=paste(folder,"Ts.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$NDVI, filename=paste(folder,"NDVI.tif",sep="/"), format="GTiff",overwrite=overwrite)
    writeRaster(object$SAVI, filename=paste(folder,"SAVI.tif",sep="/"), format="GTiff",overwrite=overwrite)
      writeRaster(object$fc, filename=paste(folder,"fc.tif",sep="/"), format="GTiff",overwrite=overwrite)
      writeRaster(object$radiance, filename=paste(folder,"radiance.tif",sep="/"), format="GTiff",overwrite=overwrite)
      writeRaster(object$reflectance, filename=paste(folder,"reflectance.tif",sep="/"), format="GTiff",overwrite=overwrite)
      writeRaster(object$LAI, filename=paste(folder,"LAI.tif",sep="/"), format="GTiff",overwrite=overwrite)
      if(!is.null(xy)){
        xydata=stack(object$albedo,object$Ts,object$NDVI,object$SAVI,object$LAI,object$fc,
                     object$radiance, object$reflectance)
        varname=c("albedo","Ts","NDVI","SAVI","LAI","fc","radiance", "reflectance")
        names(xydata)=varname
      } 
  }
  
  if(class(object)=="sebkcstack"){
    writeRaster(object$data, filename=paste(folder,"sebkcstack.tif",sep="/"), 
                format="GTiff",bandorder='BIL',overwrite=overwrite)
    if(!is.null(object$band8)){
      writeRaster(object$band8, filename=paste(folder,"band8pan.tif",sep="/"), 
                  format="GTiff",overwrite=overwrite)
    }
    
    if(!is.null(xy)){
      xydata=object$data
      #names(xydata)=varname
    }
  }
  
  if(class(object)=="ETo"){
    writeRaster(object$ETa, filename=paste(folder,"ETa.tif",sep="/"), 
                format="GTiff",overwrite=overwrite)
    
    if(!is.null(xy)){
      if(!is.null(object$xy)){
        xydata=stack(object$ETa,object$EF,object$Rn)
        names(xydata)=c("ETa","EF","Rn")
      }else{
        if(!is.null(object$EF)){
        xydata=stack(object$ETa,object$EF)
        names(xydata)=c("ETa","EF")  
        }else{
          xydata=object$ETa
          names(xydata)=c("ETa") 
        }
      }
     
    }
  }
  
  if(!is.null(xy)){
    xytable=cbind(xy,extract(xydata,xy,df=TRUE))
    filename=paste(class(object),"_xydata.txt",sep="")
    write.table(data.frame(xytable), file=paste(folder,filename,sep="/"),row.names = FALSE)
  }
  
}