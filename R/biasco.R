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

#' @include TimeSeries-class.R TimeSeriesTP-class.R BiascoTimeSeries-class.R BiascoTimeSeriesTP-class.R

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
#' @details The following univariate methods have currently been implemented in this package:
#'  \itemize{
#'  \item{M1: Simple delta change scaling of time mean}
#'  \item{M2: Delta change scaling on mean and standard deviation. If \code{type} = ratio, then Engen-Skaugen algorithm is used.}{}
#'  \item{M3: Delta change scaling on mean, standard deviation and skewness.  If \code{type} = ratio, M3 corresponds to delta change scaling of mean and standard deviation using.}
#'  \item{M4: Non-parametric quantile mapping applied in the delta change mode}
#'  \item{M5: Parametric quantile mapping applied in the delta change mode}
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
#' @import lubridate
#' @importFrom Rdpack reprompt
#' @export
#' @examples
#' 
#' library(lubridate)
#' library(BCUH)
#' 
#' data("obs_data")
#' data("model_data")
#' 
#' #Create objects which contain the bias corrected data
#' bc6 <- biasco(varsO[,"temp"], varsC[,"temp"], varsC[,"temp"], type = "abs", method = "M6")
#' bc9 <- biasco(varsO[,"temp"], varsC[,"temp"], varsC[,"temp"], type = "abs", method = "M9")
# 
#' #Visualise the results
#' plot(quantile(dat(adj(bc6)), seq(0,1,0.01)), type = "l",
#'      main="Quantile plot", xlab = "%", ylab = "Celcius", ylim = c(10,30))
#' lines(quantile(dat(adj(bc9)), seq(0,1,0.01)), col = "black",lty=2)
#' lines(quantile(varsO[,"temp"], seq(0,1,0.01)), col = "red")
#' lines(quantile(varsC[,"temp"], seq(0,1,0.01)), col = "blue")
#' lines(quantile(varsC[,"temp"], seq(0,1,0.01)), col = "blue", lty = 2)
#' legend("topleft", c("M6","M9","Obs","Ctrl","Scen"),
#'        col = c("black","black","red","blue","blue"), lty = c(1,2,1,1,2))
#' 
#'
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
           BiascoObject <- BC.abs(adj = TS(), obs = obs, ctrl = ctrl, scen = scen, method = method)
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
           BiascoObject <- BC.ratio(adj = TS(), obs = obs, ctrl = ctrl, scen = scen, method = method)
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
#' @param names names of the variables at different columns
#' @param threshold Threshold for wet day precipitation. Everything below \code{threshold} is handled as zero precipitation
#' @param ... Method-specific parameters to be implemented
#'
#' @return An object of type BiascoTimeSeriesPT which contains a matrix of bias adjusted temperature and precipitation and also
#' additional information on the used parameter values 
#' @import lubridate
#' @export
#'
#' @examples
#'
#' library(lubridate)
#' library(BCUH)
#' library(copula)
#' 
#' data("obs_data")
#' data("model_data")
#' biasco2d.object <- biasco2D(varsO,varsC,varsC,cond="P",
#'                             names = c("temp", "prec"), 
#'                             threshold=0.1)
#' 
#' #Visualise the results:
#' 
#' par(mfrow=c(1,2))
#' plot(varsO)
#' points(varsC,col="red")
#' legend("topleft",c("Obs","Ctrl"),col=c("black","red"),pch=c(1,1))
#' 
#' plot(dat(adj(biasco2d.object)))
#' points(varsO,col="red")
#' legend("topleft",c("Adj","Scen"),col=c("black","red"),pch=c(1,1))

biasco2D <- function(obs.in, ctrl.in, scen.in, names = NULL, cond = "P", separate = F, threshold = 0.1, jittering = F, fit.skew= F, ...){

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
  biascoObject2D <- .JBC(biascoObject2D, cond = cond, threshold = threshold, separate = separate, jittering = jittering, fit.skew = fit.skew)
  return(as(biascoObject2D,"BiascoTimeSeriesTP"))
}
