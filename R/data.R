#' Temperature and precipitation in Jyväskylä between 1951-2017.
#'
#' A data.frame containing the daily time series of daily mean temperature and daily precipitation
#' measures in the Jyväskylä airport (25.68 E, 62.40 N) between 1951-2017.
#' Missing values have been replaced with previous day's values.
#'
#' @format A data frame with 24139 rows and 3 variables:
#' \describe{
#'   \item{date}{time stamps for each day}
#'   \item{tas}{daily mean temperature in Celcius}
#'   \item{pr}{daily precipitation in mm/d}
#'   ...
#' }
#' @source \url{https://climexp.knmi.nl/}
"station_Jyvaskyla"


#' ALADIN5.3 EURO-CORDEX 0.11 deg resolution.
#'
#' A dataset containing the time series of 2m daily mean temperature and daily precipitation
#' simulated with CNRM-ALADIN5.3 at 0.11 degree resolution using RCP4.5 emission scenario forcings
#' in a grid box nearest to Jyväskylä airport (25.68 E, 62.40 N).
#'
#' @format A data frame with 54819 rows and 3 variables:
#' \describe{
#'   \item{date}{time stamps for each day}
#'   \item{tas}{daily mean temperature in Celcius}
#'   \item{pr}{daily precipitation in mm/d}
#'   ...
#' }
#' @source \url{https://esgf.llnl.gov/}
"ALADIN_Jyvaskyla"
