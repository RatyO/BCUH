# BCUH provides tools to bias adjust simulated time series with a number of uni- and multivariate methods.
# Copyright (C) 2018 Olle Räty
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#' @include TimeSeries-class.R BiascoTimeSeries-class.R

#' @title A wrapper function for univariate bias correction and delta changes methods.
#'
#' @description \code{biasco} returns bias corrected data as and object of type \code{BiascoTimeSeries}
#'
#' @param obs.in A vector containing the reference data
#' @param ctrl.in A vector of simulated data used for calibration
#' @param scen.in A vector of data to be adjusted
#' @param type String defining the type of adjustment. Currently \code{type} recognizes two values:
#' \itemize{
#'   \item{"abs"}{Absolute scale is used. This approach is typically applied measured in the absolute scale (e.g temperature)}
#'   \item{"ratio"}{The adjustment is made in ratio format. This approach is typically applied to non-continuous variables such as precipitation.}
#' }
#' @param method Method to be used to bias adjust \code{scen.in}. The following options are available:
#'
#' @details The following univariate methods are available:
#'  \itemize{
#'  \item{M1: Simple delta change scaling of time mean}
#'  \item{M2: Delta change scaling on mean and standard deviation. If \code{type} = ratio, then Engen-Skaugen algorithm is used.}{}
#'  \item{M3: Delta change scaling on mean, standard deviation and skewness.  If \code{type} = ratio, M3 corresponds to delta change scaling of mean and standard deviation using.}
#'  \item{M4: Parametric quantile mapping applied in the delta change mode}
#'  \item{M5: Non-parametric quantile mapping applied in the delta change mode}
#'  \item{M6: As M1 but applied in the bias correction mode}
#'  \item{M7: As M2 but applied in the bias correction mode}
#'  \item{M8: As M3 but applied in the bias correction mode}
#'  \item{M9: As M4 but applied in the bias correction mode}
#'  \item{M10: As M5 but applied in the bias correction mode}
#'  }
#'
#' @param ... Additional parameters used to handle the behavior of specific bias correction methods.
#'
#' @return An object of type \code{BiascoTimeSeries},
#'   which contains the adjusted data and parameter values related to a particular method.
#' @seealso \code{\link{biasco2d}} for two-dimensional bias correction of temperature and precipitation
#' @export
#' @examples
#' \dontrun{
#' library(BCUH)
#' library(lubridate)
#' data(station_Jyvaskyla)
#' data("ALADIN_Jyvaskyla")
#' #Extract data only for December
#' obs <- station_Jyvaskyla[!is.na(match(year(station_Jyvaskyla$date),seq(ctrl[1],ctrl[2]))), -1]$tas
#' ctrl <- ALADIN_Jyvaskyla[!is.na(match(year(ALADIN_Jyvaskyla$date),seq(ctrl[1],ctrl[2]))), -1]$tas
#' scen <- ALADIN_Jyvaskyla[!is.na(match(year(ALADIN_Jyvaskyla$date),seq(scen[1],scen[2]))), -1]$tas
#' TS <- biasco(obs.ctrl[,1],ctrl[,1],scen[,1],type="abs",method="M1")
#'
#' }
#'
biasco <- function(obs.in, ctrl.in, scen.in, type = "abs", method = "M1", ...){

  if(is.TS(obs.in)){
    obs <- obs.in
  }else{
    obs <- TS(data = obs.in, type = type)
  }

  if(is.TS(ctrl.in)){
    ctrl <- ctrl.in
  }else{
    ctrl <- TS(data = ctrl.in, type = type)
  }

  if(is.TS(scen.in)){
    scen <- scen.in
  }else{
    scen <- TS(data = scen.in, type = type)
  }

  switch(type,
         abs = {
           BiascoObject <- BC.abs(adj=TS(), obs = obs, ctrl = ctrl, scen = scen, method = method)
           switch(method,
                  M1 = {
                    BiascoObject <- .DcMean(BiascoObject)
                  },
                  M2 = {
                    BiascoObject <- .DcMeanSd(BiascoObject)
                  },
                  M3 = {
                    BiascoObject <- .DcMeanSdSkew(BiascoObject)
                  },
                  M4 = {
                    BiascoObject <- .DcQmEmpir(BiascoObject)
                  },
                  M5 = {
                    BiascoObject <- .DcQmParam(BiascoObject)
                  },
                  M6 = {
                    BiascoObject <- .BcMean(BiascoObject)
                  },
                  M7 = {
                    BiascoObject <- .BcMeanSd(BiascoObject)
                  },
                  M8 = {
                    BiascoObject <- .BcMeanSdSkew(BiascoObject)
                  },
                  M9 = {
                    BiascoObject <- .BcQmEmpir(BiascoObject)
                  },
                  M10 = {
                    BiascoObject <- .BcQmParam(BiascoObject)
                  })
         },

         ratio = {
           BiascoObject <- BC.ratio(adj=TS(), obs = obs, ctrl = ctrl, scen = scen, method = method)
           switch(method,
                  M1 = {
                    BiascoObject <- .DcMean(BiascoObject)
                  },
                  M2 = {
                    BiascoObject <- .DcMeanSd1(BiascoObject)
                  },
                  M3 = {
                    BiascoObject <- .DcMeanSd2(BiascoObject)
                  },
                  M4 = {
                    BiascoObject <- .DcQmEmpir(BiascoObject)
                  },
                  M5 = {
                    BiascoObject <- .DcQmParam(BiascoObject)
                  },
                  M6 = {
                    BiascoObject <- .BcMean(BiascoObject)
                  },
                  M7 = {
                    BiascoObject <- .BcMeanSd1(BiascoObject)
                  },
                  M8 = {
                    BiascoObject <- .BcMeanSd2(BiascoObject)
                  },
                  M9 = {
                    BiascoObject <- .BcQmEmpir(BiascoObject)
                  },
                  M10 = {
                    BiascoObject <- .BcQmParam(BiascoObject)
                  })
         }
  )
  return(as(BiascoObject,"BiascoTimeSeries"))
  # return(as(BiascoObject,"BiascoTimeSeries"))
}

biasco2D <- function(obs.in, ctrl.in, scen.in, method, ...){}
