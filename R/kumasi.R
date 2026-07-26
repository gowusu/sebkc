#' Daily Weather data of Kumasi from 1960 to the present
#'
#' A dataset containing Rain, maximum and minimum temperature, wind speed, sunshine hours,
#' and minimum and maximum Relative Humidity for Kumasi, Ghana. The record from 1960 to
#' 2015 is the historical station series (see \code{\link{kumasi1960_2015}}); from 2016
#' onward it is extended with NASA POWER daily data retrieved through \code{\link{weather}}
#' (latitude 6.72, longitude -1.60). For the appended period the incoming shortwave
#' radiation Rs is converted to sunshine hours \code{n} with the inverse Angstrom-Prescott
#' relation so that the whole series keeps a single, homogeneous set of columns and can be
#' passed to \code{\link{ETo}} and \code{\link{kc}} without special handling. Regenerate or
#' extend it at any time by re-running \code{weather()} for the new dates and re-appending.
#'
#' @name kumasi1960_present
#' @docType data
#' @author George Owusu
#' @format A data frame with 8 variables:
#' \describe{
#'   \item{DOY}{Calendar date [YYYY-mm-dd]}
#'   \item{Tmax}{Maximum Temperature [oC]}
#'   \item{Tmin}{Minimum Temperature [oC]}
#'   \item{RHmin}{Minimum Relative Humidity [per cent]}
#'   \item{RHmax}{Maximum Relative Humidity [per cent]}
#'   \item{n}{sunshine hours}
#'   \item{uz}{wind speed [m/s]}
#'   \item{P}{Daily precipitation [mm]}
#' }
#'
NULL

#' Historical Weather data of Kumasi from 1960 to 2015
#'
#' A dataset containing Rain, maximum and minimum temperature, wind speed, sunshine hours,
#' and minimum and maximum Relative Humidity.
#'
#' @name kumasi1960_2015
#' @docType data
#' @author George Owusu
#' @format A data frame with 20454 rows and 8 variables:
#' \describe{
#'   \item{DOY}{Day of the Year in dates}
#'   \item{Tmax}{Maximum Temperature [oC]}
#'   \item{Tmin}{Minimum Temperature [oC]}
#'   \item{RHmin}{Minimum Relative Humidity [per cent]}
#'   \item{RHmax}{Maximum Relative Humidity [per cent]}
#'   \item{n}{sunshine hours}
#'   \item{uz}{wind speed [m/s]}
#'   \item{P}{Daily precipitation [mm]}
#' }
#' 
NULL

#' A 2.6 RCP downscaled Weather data of Kumasi from 2016 to 2050
#'
#' A dataset containing Rain, maximum and minimum temperature, wind speed, sunshine hours,
#' and minimum and maximum Relative Humidity. The RCPs 2.6 assumes that green house 
#' gas will peak in between 2010-2020 and declines throughout 21st Century. 
#' 
#' @name kumasircp262016_2050
#' @docType data
#' @author George Owusu
#' @format A data frame with 12784 rows and 8 variables:
#' \describe{
#'   \item{DOY}{Day of the Year in dates}
#'   \item{Tmax}{Maximum Temperature [oC]}
#'   \item{Tmin}{Minimum Temperature [oC]}
#'   \item{RHmin}{Minimum Relative Humidity [per cent]}
#'   \item{RHmax}{Maximum Relative Humidity [per cent]}
#'   \item{n}{sunshine hours}
#'   \item{uz}{wind speed [m/s]}
#'   \item{P}{Daily precipitation [mm]}
#' }
#' 
NULL

#' A 4.5 RCP downscaled Weather data of Kumasi from 2016 to 2050
#'
#' A dataset containing Rain, maximum and minimum temperature, wind speed, sunshine hours,
#' and minimum and maximum Relative Humidity. The RCP 4.5 assumes that green house 
#' gas will peak in 2040 and declines throughout 21st Century. 
#' 
#' @name kumasircp452016_2050
#' @docType data
#' @author George Owusu
#' @format A data frame with 12784 rows and 8 variables:
#' \describe{
#'   \item{DOY}{Day of the Year in dates}
#'   \item{Tmax}{Maximum Temperature [oC]}
#'   \item{Tmin}{Minimum Temperature [oC]}
#'   \item{RHmin}{Minimum Relative Humidity [per cent]}
#'   \item{RHmax}{Maximum Relative Humidity [per cent]}
#'   \item{n}{sunshine hours}
#'   \item{uz}{wind speed [m/s]}
#'   \item{P}{Daily precipitation [mm]}
#' }
#' 
NULL

#' A 8.5 RCP downscaled Weather data of Kumasi from 2016 to 2050
#'
#' A dataset containing Rain, maximum and minimum temperature, wind speed, sunshine hours,
#' and minimum and maximum Relative Humidity. The RCP 8.5 assumes that green house 
#' gas will rise throughout 21st Century. 
#' 
#' @name kumasircp852016_2050
#' @docType data
#' @author George Owusu
#' @format A data frame with 12784 rows and 8 variables:
#' \describe{
#'   \item{DOY}{Day of the Year in dates}
#'   \item{Tmax}{Maximum Temperature [oC]}
#'   \item{Tmin}{Minimum Temperature [oC]}
#'   \item{RHmin}{Minimum Relative Humidity [per cent]}
#'   \item{RHmax}{Maximum Relative Humidity [per cent]}
#'   \item{n}{sunshine hours}
#'   \item{uz}{wind speed [m/s]}
#'   \item{P}{Daily precipitation [mm]}
#' }
#' 
NULL