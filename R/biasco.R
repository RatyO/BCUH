#' @title A wrapper function for univariate bias correction and delta changes methods.
#'
#' @description \code{biasco} returns bias corrected data as and object of type \code{BiascoTimeSeries}
#'
#' @param obs.in A vector containing the reference data
#' @param ctrl.in A vector of simulated data used for calibration
#' @param scen.in A vector of data to be adjusted
#' @param type String defining the type of the input. If the variable in on absolute scale, the
#' @param method Method to be used to adjust the
#' @param ... Additional parameters to used in bias correction methods
#'
#' @return An object of type \code{BiascoTimeSeries},
#'   which contains the adjusted data and parameter values related to a particular method
#' @seealso \code{\link{biasco2d}} for two-dimensional bias correction of temperature and precipitation
#' @export
#'
 biasco <- function(obs.in, ctrl.in, scen.in, type="abs", method, ...){

   obs <- TS(obs.in,method = data.type)
   ctrl <- TS(ctrl.in,type = data.type)
   scen <- TS(scen.in,type =data.type)

   bc.object <- .BC()

   BiascoObject <- new("BiascoTimeSeriesAbs")
   switch(data.type,
     abs = {
     BiascoObject <- new("AbsCorrectionObject")
     AbsCorrections(BiascoObject)
     },
     ratio = {
       BiascoObject <- new("RatioCorrectionObject")
       data(adj(BiascoObject)) <- ratioCorrections(obs,ctrl,scen,BiascoObject)
     }
   )

   return(BiascoObject)
 }

biasco2D <- function(obs.in, ctrl.in, scen.in, method, ...){}
