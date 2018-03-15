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

##' @include TimeSeries-class.R TimeSeriesTP-class.R BiascoTimeSeries-class.R BiascoTimeSeriesTP-class.R

#' @title A wrapper function for univariate bias correction and delta changes methods.
#'
#' @description \code{biasco} returns bias corrected data as an object of type \code{BiascoTimeSeries}. 
#' Several methods of different complexity are available in the package (see details).
#'
#' @param obs.in A vector containing the reference data.
#' @param ctrl.in A vector of simulated data used to calibrate the method
#' @param scen.in A vector of data for the scenario period simulation
#' @param type String defining the type of adjustment. Currently \code{type} recognizes two values:
#' \itemize{
#'   \item{abs: Absolute scale is used. This approach is typically applied measured in the absolute scale (e.g temperature)}
#'   \item{ratio: The adjustment is made in ratio format. This approach is typically applied to non-continuous variables such as precipitation.}
#' }
#' @param method Method to be used to bias adjust \code{scen.in}. The following options are available:
#'
#' @details The following univariate methods have curretnly been implemented in this package:
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
#' @importFrom Rdpack reprompt
#' @export
#' @examples
#' \dontrun{
#'library(lubridate)
#'library(BCUH)
#'
#'data("station_Jyvaskyla")
#'data("ALADIN_Jyvaskyla")
#'
#'ctrl <- c(1971,2000)
#'scen <- c(2071,2100)
#'
#'ind.obs <- which(month(station_Jyvaskyla$date) == 12 &
#'                   year(station_Jyvaskyla$date) %in% seq(ctrl[1],ctrl[2]))
#'ind.ctrl <- which(month(ALADIN_Jyvaskyla$date) == 12 &
#'                    year(ALADIN_Jyvaskyla$date) %in% seq(ctrl[1],ctrl[2]))
#'ind.scen <- which(month(ALADIN_Jyvaskyla$date) == 12 &
#'                    year(ALADIN_Jyvaskyla$date) %in% seq(scen[1],scen[2]))
#'
#'obs.ctrl <- station_Jyvaskyla[ind.obs,]
#'rcm.ctrl <- ALADIN_Jyvaskyla[ind.ctrl,]
#'rcm.scen <- ALADIN_Jyvaskyla[ind.scen,]
#'
#'dc.object <- biasco(obs.ctrl$tas, rcm.ctrl$tas, rcm.scen$tas, type = "abs", method = "M1")
#'bc.object <- biasco(obs.ctrl$tas, rcm.ctrl$tas, rcm.scen$tas, type = "abs", method = "M6")
#'
#'plot(quantile(dat(adj(bc.object)), seq(0,1,0.01)), type = "l",
#'     main="Quantile plot", xlab = "%", ylab = "Celcius", ylim = c(-30,10))
#'lines(quantile(dat(adj(dc.object)), seq(0,1,0.01)), lty = 2)
#'lines(quantile(obs.ctrl$tas, seq(0,1,0.01)), col = "red")
#'lines(quantile(rcm.ctrl$tas, seq(0,1,0.01)), col = "blue")
#'lines(quantile(rcm.scen$tas, seq(0,1,0.01)), col = "blue", lty = 2)
#'legend("topleft", c("M1","M6","Obs","Ctrl","Scen"),
#'       col = c("black","black","red","blue","blue"), lty = c(1,2,1,1,2))
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
                    BiascoObject <- .DcMeanSd(BiascoObject, ...)
                  },
                  M3 = {
                    BiascoObject <- .DcMeanSdSkew(BiascoObject, ...)
                  },
                  M4 = {
                    BiascoObject <- .DcQmEmpir(BiascoObject, ...)
                  },
                  M5 = {
                    BiascoObject <- .DcQmParam(BiascoObject, ...)
                  },
                  M6 = {
                    BiascoObject <- .BcMean(BiascoObject)
                  },
                  M7 = {
                    BiascoObject <- .BcMeanSd(BiascoObject, ...)
                  },
                  M8 = {
                    BiascoObject <- .BcMeanSdSkew(BiascoObject, ...)
                  },
                  M9 = {
                    BiascoObject <- .BcQmEmpir(BiascoObject, ...)
                  },
                  M10 = {
                    BiascoObject <- .BcQmParam(BiascoObject, ...)
                  })
         },

         ratio = {
           BiascoObject <- BC.ratio(adj=TS(), obs = obs, ctrl = ctrl, scen = scen, method = method)
           switch(method,
                  M1 = {
                    BiascoObject <- .DcMean(BiascoObject)
                  },
                  M2 = {
                    BiascoObject <- .DcMeanSd1(BiascoObject, ...)
                  },
                  M3 = {
                    BiascoObject <- .DcMeanSd2(BiascoObject, ...)
                  },
                  M4 = {
                    BiascoObject <- .DcQmEmpir(BiascoObject, ...)
                  },
                  M5 = {
                    BiascoObject <- .DcQmParam(BiascoObject, ...)
                  },
                  M6 = {
                    BiascoObject <- .BcMean(BiascoObject)
                  },
                  M7 = {
                    BiascoObject <- .BcMeanSd1(BiascoObject, ...)
                  },
                  M8 = {
                    BiascoObject <- .BcMeanSd2(BiascoObject, ...)
                  },
                  M9 = {
                    BiascoObject <- .BcQmEmpir(BiascoObject, ...)
                  },
                  M10 = {
                    BiascoObject <- .BcQmParam(BiascoObject, ...)
                  })
         }
  )
  return(as(BiascoObject,"BiascoTimeSeries"))
  # return(as(BiascoObject,"BiascoTimeSeries"))
}

#' @title Joint bias correction of temperature and precipitation.
#'
#' @description A wrapper function for a copula-based bias correction method of temperature and precipitation
#' 
#' @param obs.in Either a 2-column array, matrix or data.frame. of observed temperature and precipitation
#' @param ctrl.in Data used to calibrate the method given in the format as \code{obs.in}
#' @param scen.in Data for the scenario period simulation given in the format as \code{obs.in}
#' @param cond The conditioning order of temperature and precipitation. If \code{cond} = "P" then the precipitation marginal
#' is adjusted first, while the adjustment of temperature is made conditionally on precipitation
#' @param threshold Threshold for wet day precipitation. Everything below \code{threshold} is handled as zero precipitation
#' @param ... No additional parameter available yet
#'
#' @return An object of type BiascoTimeSeriesPT which contains a matrix of bias adjusted temperature and precipitation and also
#' additional information on the used parameter values 
#' @export
#'
#' @examples
#' #' \dontrun{
#' library(lubridate)
#' library(copula)
#' library(BCUH)

#' data("station_Jyvaskyla")
#' data("ALADIN_Jyvaskyla")

#' ctrl <- c(1971,2000)
#' scen <- c(2071,2100)

#' ind.obs <- which(month(station_Jyvaskyla$date) == 12 &
#'                   year(station_Jyvaskyla$date) %in% seq(ctrl[1],ctrl[2]))
#' ind.ctrl <- which(month(ALADIN_Jyvaskyla$date) == 12 &
#'                    year(ALADIN_Jyvaskyla$date) %in% seq(ctrl[1],ctrl[2]))
#' ind.scen <- which(month(ALADIN_Jyvaskyla$date) == 12 &
#'                    year(ALADIN_Jyvaskyla$date) %in% seq(scen[1],scen[2]))

#' obs.ctrl <- station_Jyvaskyla[ind.obs,2:3]
#' rcm.ctrl <- ALADIN_Jyvaskyla[ind.ctrl,2:3]
#' rcm.scen <- ALADIN_Jyvaskyla[ind.scen,2:3]

#' biasco2d.object <- biasco2D(obs.ctrl,rcm.ctrl,rcm.scen, names = c("tas", "pr"))
#' 
#' Visualise the results
#' par(mfrow=c(1,2))

#' plot(obs.ctrl)
#' points(rcm.ctrl,col="red")
#' legend("topleft",c("Obs","Ctrl"),col=c("black","red"),pch=c(1,1))

#' plot(dat(adj(biasco2d.object)))
#' points(rcm.scen,col="red")
#' legend("topleft",c("Adj","Scen"),col=c("black","red"),pch=c(1,1))
#' }
biasco2D <- function(obs.in, ctrl.in, scen.in, names = NULL, cond = "P", threshold = 0.1, ...){

  if(is.null(names)) names <- colnames(obs.in)
  
  if(is.TS2D(obs.in)){
    obs <- obs.in
  }else{
    obs <- TS_2D(data = obs.in, names = names, type = type)
  }
  obs
  if(is.TS2D(ctrl.in)){
    ctrl <- ctrl.in
  }else{
    ctrl <- TS_2D(data = ctrl.in, names = names, type = type)
  }
  
  if(is.TS2D(scen.in)){
    scen <- scen.in
  }else{
    scen <- TS_2D(data = scen.in, names = names, type = type)
  }
  
  biascoObject2D <- BC.joint(obs = obs, ctrl = ctrl, scen = scen)
  biascoObject2D <- .JBC(biascoObject2D, cond = cond, threshold = threshold)
  return(as(biascoObject2D,"BiascoTimeSeriesTP"))
}
