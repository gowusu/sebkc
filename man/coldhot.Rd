% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coldhotxy.R
\name{coldhot}
\alias{coldhot}
\alias{coldhot.default}
\title{Automatic computation of hot and cold pixels from Radiometric Surface Temperature (Ts)}
\usage{
coldhot(Tmin, Tmax, Ts = NULL, folder = NULL, welev = NULL)

\method{coldhot}{default}(Tmin, Tmax, Ts = NULL, folder = NULL,
  welev = NULL)
}
\arguments{
\item{Tmin}{Numeric. Minimum surface Temperature in  map units (eg. Kilvin) that you want computation from. 
You can keep changing untill you get your desired values}

\item{Tmax}{Numeric. Maximum surface Temperature in map units (eg. Kilvin) that you want computation from. 
You can keep changing untill you get your desired values}

\item{Ts}{A RasterLayer data that indicates radiometric 
surface temperature values from remote sensing image, 
preferably in Kelvin (K). You can also point to raster file on your 
computer.}

\item{folder}{An original  directory of the images that contains all 
landsat 5,7 or 8 bands and metadata.
At the moment only Landsat folder is supported}

\item{welev}{Weather station elevation in meters}
}
\value{
\itemize{
\item{Ts:} { The selected or input Ts}
\item{Tshot:} { The Temperature of hot pixel}
\item{Tscold:} { The Temperature of cold pixel}
\item{xyhot: } {combined x and y coordinates}
\item{candidates:}{ The similar candidates' pixels that can be used}
}
}
\description{
This function computes temperature coordinates of a hot and cold pixels based on only
Surface Temperature. This function uses ranges of temperature.
}
\details{
The function first divides the selected area into a 
number of clusters. The mean and standard deviations are
computed. The cluster means with the lowest standard deviations
and qualifying for upper (Ts) and lower (NDVI) limits 
(see above) are selected. The provided upper and lower probabilities
 are also used to select the final Ts.
}
\examples{
\dontrun{
 folder=system.file("extdata","stack",package="sebkc")
  model=coldhot(Tmin=300,Tmax=310,folder=folder,welev=170)
  spplot(model$Tshot)
  spplot(model$Tscold)
  modauto=sebal(folder = folder,welev = 380,xycold=model$xycold,xyhot=model$xyhot)}
}
\author{
George Owusu
}
\references{
Allen, R. G., Burnett, B., Kramber, W., Huntington, J., 
Kjaersgaard, J., Kilic, A., Kelly, C., & Trezza, R. 2013. Automated 
Calibration of the METRIC-Landsat Evapotranspiration Process. 
JAWRA Journal of the American Water Resources Association, 49(3): 563-576.
}
\seealso{
\code{\link{coldTs}},  \code{\link{hotTs}}, \code{\link{coldT2}}
}

