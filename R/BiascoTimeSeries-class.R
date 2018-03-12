#' @include GenericMethods.R TimeSeries-class.R

#' export BC
BC <- setClass("BiascoTimeSeries",
         slots = c(adj = "TimeSeries", method = "character", bc.attributes = "list")
)

setMethod(f = "initialize",
          signature = "BiascoTimeSeries",
          function(.Object, adj = TS(), method = NA_character_){
            .Object@adj <- adj
            .Object@method <- method
            .Object@bc.attributes <- list()
            validObject(.Object)
            return(.Object)
          })

#getters

setMethod(f = "adj",
          signature = "BiascoTimeSeries",
          definition = function(object){
            return(object@adj)
          })

setMethod(f = "method",
          signature = "BiascoTimeSeries",
          definition = function(object){
            return(object@method)
          })

#Setters

setMethod(f = "adj<-",
           signature = "BiascoTimeSeries",
           definition = function(object, value){
             object@adj <- value
             validObject(object)
             return(object)
           })

setMethod(f = "method<-",
          signature = "BiascoTimeSeries",
          definition = function(object, value){
            object@method <- value
            validObject(object)
            return(object)
          })

# setMethod(f = "as",
#           signature = "BiascoTimeSeries",
#           definition = function(from,to){
#
#           })

#Show
setMethod(f = "show",
          signature = "BiascoTimeSeries",
          definition = function(object){
            cat("*** Class BiascoTimeSeries, method show *** \n")
            cat("* method = "); print (object@method)
            cat("* bc.attributes = "); print (object@bc.attributes)
            cat("* adj@nvar = "); print (object@adj@nvar)
            cat("* adj@Dim = "); print (object@adj@Dim)
            cat("* adj@data (limited to first 100 values) = \n")
            if(length(object@adj@data)!=0){
              print(formatC(object@adj@data[1:min(length(object@adj@data),100)]),quote=FALSE)
            }else{}
            cat("******* End Show (BiascoTimeSeries) ******* \n")
          })
