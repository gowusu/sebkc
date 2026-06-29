#' sebkc: Surface Energy Balance and Crop Coefficient Evapotranspiration
#' Estimation
#'
#' Computes and integrates surface energy balance components of
#' evapotranspiration into the FAO56 water balance model, providing the
#' SEBAL, METRIC, SEBI, SSEB, SEBS and TSEB models.
#'
#' @keywords internal
#'
#' Base-package functions used internally. Generics that \pkg{raster} also
#' provides as methods (e.g. plot, hist, quantile, text, lines, points) are
#' intentionally NOT imported here because they arrive via \code{import(raster)};
#' importing them again would trigger "replacing previous import" warnings.
#'
#' @importFrom stats na.omit kmeans lm coef cor
#' @importFrom graphics abline axis legend par
#' @importFrom grDevices dev.new
#' @importFrom utils read.csv read.table read.delim write.csv write.table
"_PACKAGE"
