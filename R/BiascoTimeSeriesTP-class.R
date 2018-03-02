.BC_2D <- setClass("BiascoTimeSeriesTP",
                slots = c(adj = "TimeSeriesTP", method = "character", bc.attributes = "list")
)

setMethod(f = "initialize",
          signature = "BiascoTimeSeriesTP",
          function(.Object, adj = new("TimeSeriesTP"), method = NA_character_, cond = "T"){
            .Object@adj <- adj
            .Object@method <- method
            .Object@bc.attributes <- list()
            .Object@bc.attributes[["cond"]] <- cond
            validObject(.Object)
            return(.Object)
          })

setMethod(f = "adj",
          signature = "BiascoTimeSeriesTP",
          definition = function(object){
            return(object@adj)
          })

setMethod(f = "method",
          signature = "BiascoTimeSeriesTP",
          definition = function(object){
            return(object@method)
          })

setMethod(f = "adj<-",
          signature = "BiascoTimeSeriesTP",
          definition = function(object, value){
            object@adj <- value
            validObject(object)
            return(object)
          })

setMethod(f = "method<-",
          signature = "BiascoTimeSeriesTP",
          definition = function(object, value){
            object@method <- value
            validObject(object)
            return(object)
          })

setMethod(f = "show",
          signature = "BiascoTimeSeriesTP",
          definition = function(object){
            cat("*** Class BiascoTimeSeriesTP, method show *** \n")
            cat("* method = "); print (object@method)
            cat("* bc.attributes = "); print (object@bc.attributes)
            cat("* adj@nvar = "); print (object@adj@nvar)
            cat("* adj@Dim = "); print (object@adj@Dim)
            cat("* adj@data (limited to first 100 values) = \n")
            if(nrow(object@adj@data)!=0){
              print(format(object@adj@data[1:min(nrow(object@adj@data),100),2]),quote=FALSE)
            }else{}
            cat("******* End Show (BiascoTimeSeriesTP) ******* \n")
          })
