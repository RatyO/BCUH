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
setGeneric(name = "data",
           def = function(object){
             standardGeneric("data")
           })

setGeneric(name = "data<-",
           def = function(object,i,value){
             standardGeneric("data<-")
           })

setGeneric(name = "data[",
           def = function(object,i,...,value){
             standardGeneric("data[")
           })

#' @title get method
#'
#' @description \code{method} returns the string stored in the slot method.
#'
#' @param object An object
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
#' @param object
#' @param value
#'
#' @return
#' @export
#'
#' @examples
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

#DC

setGeneric(".DcMean",
          def = function(.Object, ...){
             standardGeneric(".DcMean")
           })

setGeneric(name = ".DcMeanSd",
           def = function(.Object, ...){
             standardGeneric(".DcMeanSd")
           })

setGeneric(name = ".DcMeanSd1",
           def = function(.Object, ...){
             standardGeneric(".DcMeanSd1")
           })

setGeneric(name = ".DcMeanSd2",
           def = function(.Object, ...){
             standardGeneric(".DcMeanSd2")
           })

setGeneric(name = ".DcMeanSdSkew",
           def = function(.Object, ...){
             standardGeneric(".DcMeanSdSkew")
           })


setGeneric(name = ".DcQmEmpir",
           def = function(.Object, ...){
             standardGeneric(".DcQmEmpir")
           })

setGeneric(name = ".DcQmParam",
           def = function(.Object, ...){
             standardGeneric(".DcQmParam")
           })

#BC
setGeneric(".BcMean",
           def = function(.Object, ...){
             standardGeneric(".BcMean")
           })

setGeneric(name = ".BcMeanSdSkew",
           def = function(.Object, ...){
             standardGeneric(".BcMeanSdSkew")
           })

setGeneric(name = ".BcMeanSd",
           def = function(.Object, ...){
             standardGeneric(".BcMeanSd")
           })

setGeneric(name = ".BcMeanSd1",
           def = function(.Object, ...){
             standardGeneric(".BcMeanSd1")
           })

setGeneric(name = ".BcMeanSd2",
           def = function(.Object, ...){
             standardGeneric(".BcMeanSd2")
           })

setGeneric(name = ".BcQmEmpir",
           def = function(.Object, ...){
             standardGeneric(".BcQmEmpir")
           })

setGeneric(name = ".BcQmParam",
           def = function(.Object, ...){
             standardGeneric(".BcQmParam")
           })

setGeneric(name = ".JBC",
           def = function(.Object, cond, threshold, ...){
             standardGeneric(".JBC")
           })

