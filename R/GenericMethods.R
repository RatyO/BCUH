#A useful class union to handle missing or NULL input for class methods

#setClassUnion("missingOrNULL",c("missing", "NULL"))

#Getter generics

#' @title get nvar
#'
#' @description \code{nvar} returns the value stored in the slot nvar.
#'
#' @param object An object
#'
#' @export
setGeneric(name = "nvar",
           def = function(object){
             standardGeneric("nvar")
           })

#' @title get type
#'
#' @description \code{type} returns the value stored in the slot type.
#'
#' @param object An object
#'
#' @export
setGeneric(name = "type",
           def = function(object){
             standardGeneric("type")
           })

#' @title get data
#'
#' @description \code{data} returns the values stored in the slot data.
#'
#' @param object An object
#'
#' @export
setGeneric(name = "dat",
           def = function(object){
             standardGeneric("dat")
           })

setGeneric(name = "dat<-",
           def = function(object,i,value){
             standardGeneric("dat<-")
           })

setGeneric(name = "dat[",
           def = function(object,i,...,value){
             standardGeneric("dat[")
           })

#' @title get method
#'
#' @description \code{method} returns the string stored in the slot method.
#'
#' @param object An object
#' @param Name of the method
#' @param ... Additional arguments
#'
#' @export
setGeneric(name = "method",
           def = function(object, value, ...){
             standardGeneric("method")
           })

#' @title get adj
#'
#' @description \code{adj} returns the elements of adj, of type either  in an object
#'
#' @param object An object
#' @param ... Additional arguments
#'
#' @export
setGeneric(name = "adj",
           def = function(object, ...){
             standardGeneric("adj")
           })

#' @title get obs
#'
#' @description \code{obs} returns the content of obs,
#'   which can be either type TimeSeries or TimeSeriesTP and can also be part of other object.
#'
#' @param object An object
#' @param ... Additional arguments
#'
#' @export
setGeneric(name = "obs",
           def = function(object, ...){
             standardGeneric("obs")
           })

#' @title get ctrl
#'
#' @description \code{ctrl} returns the content of ctrl,
#'   which can be either type TimeSeries or TimeSeriesTP and can also be part of other object.
#' @param object An object
#' @param ... Additional arguments
#'
#' @export
setGeneric(name = "ctrl",
           def = function(object, ...){
             standardGeneric("ctrl")
           })

#' @title scen
#'
#' @description \code{scen} returns the content of scen,
#'   which can be either type TimeSeries or TimeSeriesTP and can also be part of other object.
#' @param object An object
#' @param ... Additional arguments
#'
#' @export
setGeneric(name = "scen",
           def = function(object, ...){
             standardGeneric("scen")
           })

#Setter generics
#' Set nvar value
#'
#' @param object An object
#' @param value Value to be set slot nvar
#'
#' @return new value for nvar
#'
#' @keywords internal
#'
setGeneric(name = "nvar<-",
           def = function(object, value){
             standardGeneric("nvar<-")
           })

#' Set a new value for slot type
#'
#' @param object An object
#' @param value Value given for \code{slot}
#'
#' @export
#'
setGeneric(name = "type<-",
           def = function(object, value){
             standardGeneric("type<-")
           })

setGeneric(name = "method<-",
           def = function(object, value){
             standardGeneric("method<-")
           })

setGeneric(name = "adj<-",
           def = function(object, value){
             standardGeneric("adj<-")
           })

setGeneric(name = "obs<-",
           def = function(object, value){
             standardGeneric("obs<-")
           })

setGeneric(name = "ctrl<-",
           def = function(object, value){
             standardGeneric("ctrl<-")
           })

setGeneric(name = "scen<-",
           def = function(object, value){
             standardGeneric("scen<-")
           })

setGeneric(name = "is.TS",
           def = function(object){
             standardGeneric("is.TS")
           })

setGeneric(name = "is.TS2D",
           def = function(object){
             standardGeneric("is.TS2D")
           })

#DC

#' Mean delta change
#'
#' A generic for mean delta change scaling
#'
#' @param .Object An object
#' @param ... Additional parameters required by a specific implementation
#'
#' @return Delta change adjusted time series
#'
#' @keywords internal
setGeneric(".DcMean",
          def = function(.Object, ...){
             standardGeneric(".DcMean")
           })

#' Delta change scaling of mean and standard deviation
#'
#' A generic for delta change scaling of mean and standard deviation.
#'
#' @param .Object An object
#' @param ... Additional parameters required by a specific implementation
#'
#' @return Delta change adjusted time series
#'
#' @keywords internal
setGeneric(name = ".DcMeanSd",
           def = function(.Object, ...){
             standardGeneric(".DcMeanSd")
           })

#' Delta change scaling of mean and standard deviation
#'
#' A generic for delta change scaling of mean and standard deviation.
#'
#' @param .Object An object
#' @param ... Additional parameters required by a specific implementation
#'
#' @return Delta change adjusted time series
#'
#' @keywords internal
setGeneric(name = ".DcMeanSd1",
           def = function(.Object, ...){
             standardGeneric(".DcMeanSd1")
           })

#' Delta change scaling of mean and standard deviation
#'
#' A generic for delta change scaling of mean and standard deviation.
#'
#' @param .Object An object
#' @param ... Additional parameters required by a specific implementation
#'
#' @return Delta change adjusted time series
#'
#' @keywords internal
setGeneric(name = ".DcMeanSd2",
           def = function(.Object, ...){
             standardGeneric(".DcMeanSd2")
           })

#' Delta change scaling of mean, standard deviation and skewness
#'
#' A generic for delta change scaling of mean and standard deviation.
#'
#' @param .Object An object
#' @param ... Additional parameters required by a specific implementation
#'
#' @return Delta change adjusted time series
#'
#' @keywords internal
setGeneric(name = ".DcMeanSdSkew",
           def = function(.Object, ...){
             standardGeneric(".DcMeanSdSkew")
           })

#' Empirical quantile mapping applied in the delta change form
#'
#' A generic for quantile mapping applied in the delta change mode.
#'
#' @param .Object An object
#' @param ... Additional parameters required by a specific implementation
#'
#' @return Delta change adjusted time series
#'
#' @keywords internal
setGeneric(name = ".DcQmEmpir",
           def = function(.Object, ...){
             standardGeneric(".DcQmEmpir")
           })

#' Parametric quantile mapping applied in the delta change form
#'
#' A generic for quantile mapping applied in the delta change mode.
#'
#' @param .Object An object
#' @param ... Additional parameters required by a specific implementation
#'
#' @return Delta change adjusted time series
#'
#' @keywords internal
setGeneric(name = ".DcQmParam",
           def = function(.Object, ...){
             standardGeneric(".DcQmParam")
           })


#' Bias correction of the mean
#'
#' A generic for bias correction of mean
#'
#' @param .Object An object
#' @param ...  Additional parameters required by a specific implementation
#'
#' @return Bias corrected time series
#'
#' @keywords internal
setGeneric(".BcMean",
           def = function(.Object, ...){
             standardGeneric(".BcMean")
           })

#' Bias correction of the mean and standard deviation
#'
#' A generic for bias correction of mean and standard deviation
#'
#' @param .Object An object
#' @param ...  Additional parameters required by a specific implementation
#'
#' @return Bias corrected time series
#'
#' @keywords internal
setGeneric(name = ".BcMeanSd",
           def = function(.Object, ...){
             standardGeneric(".BcMeanSd")
           })

#' Bias correction of the mean and standard deviation
#'
#' A generic for bias correction of mean and standard deviation
#'
#' @param .Object An object
#' @param ...  Additional parameters required by a specific implementation
#'
#' @return Bias corrected time series
#'
#' @keywords internal
setGeneric(name = ".BcMeanSd1",
           def = function(.Object, ...){
             standardGeneric(".BcMeanSd1")
           })

#' Bias correction of the mean and standard deviation
#'
#' A generic for bias correction of mean and standard deviation
#'
#' @param .Object An object
#' @param ...  Additional parameters required by a specific implementation
#'
#' @keywords internal
setGeneric(name = ".BcMeanSd2",
           def = function(.Object, ...){
             standardGeneric(".BcMeanSd2")
           })

#' Bias correction of the mean, standard deviation and skewness
#'
#' A generic for bias correction of mean, stardard deviation and skewness
#'
#' @param .Object An object
#' @param ...  Additional parameters required by a specific implementation
#'
#' @return Bias corrected time series
#'
#' @keywords internal
setGeneric(name = ".BcMeanSdSkew",
           def = function(.Object, ...){
             standardGeneric(".BcMeanSdSkew")
           })


#' Bias correction using empirical quantile mapping
#'
#' A generic for bias correction using empirical quantile mapping
#'
#' @param .Object An object
#' @param ... Additional parameters required by a specific implementation
#'
#' @return Bias corrected time series
#' @keywords internal
setGeneric(name = ".BcQmEmpir",
           def = function(.Object, ...){
             standardGeneric(".BcQmEmpir")
           })

#' Bias correction using parametric quantile mapping
#'
#' A generic for bias correction using parametric quantile mapping
#'
#' @param .Object An object
#' @param ... Additional parameters required by a specific implementation
#'
#' @return Bias corrected time series
#' @keywords internal
setGeneric(name = ".BcQmParam",
           def = function(.Object, ...){
             standardGeneric(".BcQmParam")
           })

#' Joint bias correction of temperature and precipitation
#'
#' A generic for joint bias correction of temperature and precipitation using the algorithm by Li et al. (2014)
#'
#' @param .Object An object
#' @param cond character flag indicating the order of conditioning. If value "T" is given temperature is bias corrected first,
#' while "P" is selected the adjustment of temperature is conditioned on precipitation.
#' @param threshold The value used to divide precipitation time series to wet and dry days.
#' @param ... Additional parameters required by the specific implementation of \code{JBC}.
#'
#' @return Bias corrected time series of both temperature and precipitation.
#' @keywords internal
setGeneric(name = ".JBC",
           def = function(.Object, cond, threshold, ...){
             standardGeneric(".JBC")
           })
