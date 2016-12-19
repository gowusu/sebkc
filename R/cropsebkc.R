#' @title Reprojecting and cropping images
#' 
#' 
#' @description This function returns the subset of x object as defined by y extent object. 
#' See \code{\link[raster]{crop}} and \code{\link[raster]{extent}} for more details.
#' The difference between \code{\link[raster]{crop}}  and cropsebkc is the ability for the
#' latter to reproject y object with x object projection.
#'
#' @param x object to be extracted
#' @param y extent object or raster object or polygon from which an Extent 
#'  object can be extracted. A polygon or raster will be reprojected to 
#'  conform x
#' @author George Owusu
#'
#' @return cropped  x and reprojected y
#' @seealso  \code{\link[raster]{crop}} 
#' @export
#'
#' @examples
#' \dontrun{
#' folder=system.file("extdata","stack",package="SEBKc")
#' stack=sebkcstack(folder=folder)
#' data=cropsebkc(stack$data,stack$data)
#' }
#' 
cropsebkc<-function(x,y){
if(is.numeric(y)){
  # extent=y
  raster1y=crop(x, y)
  yproj=y
  #print("yes")
}else{
  if(class(y)=="RasterLayer"){
    xy=xyFromCell(y, 1) 
    if(xy[[1]]<300||xy[[2]]<300){
      projection(y) <- CRS("+init=epsg:4326") # WGS 84
      CRS.new=CRS( proj4string(x))
      yproj=projectRaster(y, crs=CRS.new)
      
    }
    else{
      CRS.new=CRS( proj4string(x))
      yproj=projectRaster(y, crs=CRS.new) 
    }
    raster1y=crop(x, yproj)
  }else{
    #extent
    if(class(y)=="Extent"){
      # print("Yes")
      raster1y=crop(x, y)
      yproj=y
      
      
    }else{
      
      
      #polygons functions here
      xy=coordinates(y)
      if(xy[[1]]<300||xy[[2]]<300){
        projection(y) <- CRS("+init=epsg:4326") # WGS 84
        CRS.new=CRS( proj4string(x))
        yproj <- spTransform(y, CRS.new)
        
      }else{
        projection(y) <- proj4string(x)
        yproj=y 
        
      }
      
      raster1y=crop(x, yproj)
      
    }
  }
  
}
  list(x=raster1y,y=yproj)
}