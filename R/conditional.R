#' Spatial Conditional function on raster object
#' 
#' It computes equivalent of logical if function in R
#'
#' @param condition logical
#' @param trueValue logical
#' @param falseValue logical
#' @author George Owusu
#' @return Return True or False
#' @export
#' @examples
#' \dontrun{
#' NDVI=raster(system.file("extdata","NDVI.grd",package="sebkc"))
#' albedo=raster(system.file("extdata","albedo.grd",package="sebkc"))
#' LAI=raster(system.file("extdata","LAI.grd",package="sebkc"))
#'   eNB=rastercon( NDVI<0 & albedo<0.47,0.99,rastercon(LAI>=3,0.98,0.97+(LAI*0.0033)))
#'   }

rastercon=function(condition, trueValue,falseValue){
    return (condition*trueValue+(!condition)*falseValue)
}

