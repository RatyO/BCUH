#' export TS
TS <- setClass("TimeSeries",
         slots = list(data = "numeric"),
         prototype = prototype(data = NA_real_),
         contains = "DataFormat"
)

setMethod(f = "initialize",
          signature = "TimeSeries",
          function(.Object, data="ANY", nvar = NA_integer_, ...){
            .Object@data <- as.numeric(data)

            if(length(.Object@data) == 1 && is.na(.Object@data)){
              d <- NA_integer_
            }else if(is.null(dim(.Object@data))){
              d <- as.integer(length(.Object@data))
            }else{
              d <- as.integer(dim(.Object@data))
            }
            .Object <- callNextMethod(.Object, nvar = as.integer(nvar), Dim = d)
            validObject(.Object)
            return(.Object)
          })

#DataFormat methods and generics

setMethod(f = "data",
          signature = "TimeSeries",
          definition = function(object){
            return(object@data)
          })

setMethod(f = "data<-",
          signature = "TimeSeries",
          definition = function(object,value){
            object@data <- value
            object@Dim <- length(object@data)
            validObject(object)
            return(object)
          })

setMethod(f = "data[",
          signature = "TimeSeries",
          definition = function(object,i,value){
            object@data[i] <- value
            object@Dim <- length(object@data)
            validObject(object)
            return(object)
          })

setMethod(f = "dim",
          signature(x = "TimeSeries"),
          definition = function(x) x@Dim, valueClass = "integer")


# setMethod(f = "dim<-",
#           signature(x = "TimeSeries", value = "integer"),
#           definition = function(x, value=NA_integer_){
#             stop("Dim is calculated internally, the value is not changed.")
#           })

#Show the full class object in a nicer format.
setMethod(f = "show",
          signature = "TimeSeries",
          definition = function(object){
            cat("*** Class TimeSeries, method show *** \n")
            cat("* nvar = "); print (object@nvar)
            cat("* Dim = "); print (object@Dim)
            cat("* data (limited to first 100 values) = \n")
            if(length(object@data)!=0){
              print(formatC(object@data[1:min(length(object@data),100)]),quote=FALSE)
            }else{}
            cat("******* End Show (TimeSeries) ******* \n")
          })
