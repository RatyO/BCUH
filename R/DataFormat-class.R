#' @include GenericMethods.R

setClassUnion("missingOrNULL", c("missing", "NULL"))

setClass("DataFormat",
         slots = list(nvar = "integer", Dim = "integer"),
#         prototype = prototype(nvar = NA_integer_, Dim = NA_integer_),
         contains = "VIRTUAL"
)

#initialization

setMethod(f = "initialize",
          signature = "DataFormat",
          definition = function(.Object, nvar = NA_integer_, Dim = NA_integer_){
            .Object@nvar <- nvar
            .Object@Dim <- Dim
            validObject(.Object)
            return(.Object)
          })

#DataFormat methods and generics

setMethod(f = "nvar",
          signature = "DataFormat",
          definition = function(object){
            return(object@nvar)
          })

setMethod(f = "nvar<-",
          signature = "DataFormat",
          definition = function(object, value){
            object@nvar <- value
            validObject(object)
            return(object)
          })
