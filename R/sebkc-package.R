#' sebkc: Surface Energy Balance and Crop Coefficient Evapotranspiration
#' Estimation
#'
#' Computes and integrates surface energy balance components of
#' evapotranspiration into the FAO56 water balance model, providing the
#' SEBAL, METRIC, SEBI, SSEB, SEBS and TSEB models.
#'
#' Generics that the 'raster' package also provides as methods (such as plot,
#' hist, quantile, text, lines and points) are intentionally not imported
#' here, because they are made available through import(raster); importing
#' them again would trigger "replacing previous import" warnings.
#'
#' @keywords internal
#' @importFrom stats na.omit kmeans lm coef cor pf
#' @importFrom graphics abline axis legend par
#' @importFrom grDevices dev.new
#' @importFrom utils read.csv read.table read.delim write.csv write.table data
"_PACKAGE"
