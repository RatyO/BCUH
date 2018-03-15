#' @include GenericMethods.R DataFormat-class.R

#' export TS_2D
TS_2D <- setClass("TimeSeriesTP",
                slots = list(data = "ANY"),
#                prototype = prototype(data = matrix(c(NA,NA), ncol = 2)),
                contains = "DataFormat"
)

setMethod(f = "initialize",
          signature = "TimeSeriesTP",
          function(.Object, data = matrix(c(NA,NA),ncol=2,dimnames = list(c(1), c("NA", "NA"))),
          nvar = NA_integer_, names = NA_character_, ...){
            
            d <- dim(data)
            if(length(d) != 2) stop("Data not a 2-dimensional table!")
            if(d[2] != 2) stop("Only two columns allowed!")
            
            .Object@data <- as.matrix(data)

            if(is.na(nvar))nvar <- ncol(.Object@data)
            if(!any(is.na(names)) & any(is.null(colnames(data)))){
              colnames(.Object@data) <- name
            }else{
              if(!any(is.na(names))) warning("Overwriting original colnames")
              colnames(.Object@data) <- colnames(data)
            }
            .Object <- callNextMethod(.Object, nvar = nvar, Dim = as.integer(d))
            validObject(.Object)
            return(.Object)
          })

#DataFormat methods and generics
setMethod(f = "dat",
          signature = "TimeSeriesTP",
          definition = function(object){
            return(object@data)
          })

setMethod(f = "dat<-",
          signature = "TimeSeriesTP",
          definition = function(object,value){
            object@data <- value
            object@Dim <- length(object@data)
            validObject(object)
            return(object)
          })

setMethod(f = "dat[",
          signature = "TimeSeriesTP",
          definition = function(object,i,j,value){
            object@data[i,j] <- value
            object@Dim <- length(object@data)
            validObject(object)
            return(object)
          })

setMethod(f = "dim",
          signature(x = "TimeSeriesTP"),
          definition = function(x) x@Dim, valueClass = "integer")


setMethod(f = "is.TS2D",
          signature( "ANY"),
          definition = function(object){
            if(is.object(object) == T){
              if(class(object) == "TimeSeriesTP") return(TRUE)
              return(FALSE)
            }else{
              return(FALSE)
            }
          })

setMethod(f = "show",
          signature = "TimeSeriesTP",
          definition = function(object){
            cat("*** Class TimeSeriesTP, method show *** \n")
            cat("* nvar = "); print (object@nvar)
            cat("* Dim = "); print (object@Dim)
            cat("* data (limited to first 100 values) = \n")
            if(nrow(object@data)!=0){
              print(format(object@data[1:min(nrow(object@data),100),]),quote=FALSE)
            }else{}
            cat("******* End Show (TimeSeriesTP) ******* \n")
          })

