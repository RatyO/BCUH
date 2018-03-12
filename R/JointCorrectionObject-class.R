#' @include TimeSeriesTP-class.R BiascoTimeSeriesTP-class.R
.BC.joint <- setClass("JointCorrection",
                    contains = "BiascoTimeSeriesTP",
                    slots = c(obs = "TimeSeriesTP", ctrl = "TimeSeriesTP", scen = "TimeSeriesTP"),
                    prototype = prototype(obs = new("TimeSeriesTP"), ctrl = new("TimeSeriesTP"), scen = new("TimeSeriesTP"))
)

#initialize
setMethod(f = "initialize",
          signature = "JointCorrection",
          function(.Object, ..., obs = new("TimeSeriesTP"), ctrl = new("TimeSeriesTP"),
                   scen = new("TimeSeriesTP")){
            .Object <- callNextMethod(.Object, adj = scen)
            .Object@obs <- obs
            .Object@ctrl <- ctrl
            .Object@scen <- scen
            validObject(.Object)
            return(.Object)
          })

#getters

setMethod(f = "obs",
          signature = "JointCorrection",
          definition = function(object){
            return(object@obs)
          })

setMethod(f = "ctrl",
          signature = "JointCorrection",
          definition = function(object){
            return(object@ctrl)
          })

setMethod(f = "scen",
          signature = "JointCorrection",
          definition = function(object){
            return(object@scen)
          })

#Setters

setMethod(f = "obs<-",
          signature = "JointCorrection",
          definition = function(object, value){
            .Object@obs <- value
            validObject(object)
            return(object)
          })

setMethod(f = "ctrl<-",
          signature = "JointCorrection",
          definition = function(object, value){
            object@ctrl <- value
            validObject(object)
            return(object)
          })

setMethod(f = "scen<-",
          signature = "JointCorrection",
          definition = function(object, value){
            object@scen <- value
            validObject(object)
            return(object)
          })

setMethod(f = "show",
          signature = "JointCorrection",
          definition = function(object){
            cat("*** Class JointCorrection object, method show *** \n")
            cat("* method = "); print (object@method)
            cat("* bc.attributes = "); print (object@bc.attributes)
            cat("* adj@nvar = "); print (object@adj@nvar)
            cat("* adj@Dim = "); print (object@adj@Dim)
            cat("* adj@data (limited to first 100 values) = \n")
            if(nrow(object@adj@data)!=0){
              print(format(object@adj@data[1:min(nrow(object@adj@data),100),]),quote=FALSE)
            }else{}
            cat("* obs@nvar = "); print (object@obs@nvar)
            cat("* obs@Dim = "); print (object@obs@Dim)
            cat("* obs@data (limited to first 100 values) = \n")
            if(nrow(object@obs@data)!=0){
              print(format(object@obs@data[1:min(nrow(object@obs@data),100),]),quote=FALSE)
            }else{}
            cat("* ctrl@nvar = "); print (object@ctrl@nvar)
            cat("* ctrl@Dim = "); print (object@ctrl@Dim)
            cat("* ctrl@data (limited to first 100 values) = \n")
            if(nrow(object@ctrl@data)!=0){
              print(format(object@ctrl@data[1:min(nrow(object@ctrl@data),100),]),quote=FALSE)
            }else{}
            cat("* scen@nvar = "); print (object@scen@nvar)
            cat("* scen@Dim = "); print (object@scen@Dim)
            cat("* scen@data (limited to first 100 values) = \n")
            if(nrow(object@scen@data)!=0){
              print(format(object@scen@data[1:min(nrow(object@scen@data),100),]),quote=FALSE)
            }else{}
            cat("******* End Show (JointCorrectionObject) ******* \n")
          })

setMethod(f = ".JBC",
          signature = "JointCorrection",
          definition = function(.Object, cond="P", threshold = 0.1){

            #Create data frames for the data. Update these data frames whenever necessary.
            obs <- list(T=.Object@obs@data[,1],P=.Object@obs@data[,2],
                        margT=NA,margP=NA,pwet=NA,
                        wet=NA,dry=NA,cpar1=NA,cpar2=NA)
            ctrl <- list(T=.Object@ctrl@data[,1],P=.Object@ctrl@data[,2],
                         margT=NA,margP=NA,pwet=NA,
                         wet=NA,dry=NA,cpar1=NA,cpar2=NA)
            scen <- list(T=.Object@scen@data[,1],P=.Object@scen@data[,2],
                         margT=NA,margP=NA,pwet=NA,
                         wet=NA,dry=NA,cpar1=NA,cpar2=NA)
            adj <- list(T=array(NA,dim=c(length(.Object@scen@data[,1]))),
                        P=array(NA,dim=c(length(.Object@scen@data[,2]))),
                        margT=NA,margP=NA,pwet=NA,
                        wet=NA,dry=NA,cpar1=NA,cpar2=NA)
print(threshold)
            tmp <- .ppRain(ref = obs$P, adj = ctrl$P, ref.threshold = threshold)
            ctrl$P <- tmp$adj

            tmp <- .ppRain(ref = obs$P, adj = scen$P, ref.threshold = threshold)
            scen$P <- tmp$adj

            if(tmp$nDry[2]==length(scen$P)){
              warning("All dry days!")
              .Object@adj@data[,1] <- scen$T
              .Object@adj@data[,2] <- scen$P
              return(.Object)
            }

            obs$wet <- which(obs$P>threshold)
            obs$dry <- which(obs$P<=threshold)
            obs$pwet <- length(obs$wet)/length(obs$P)

            ctrl$wet <- which(ctrl$P>0)
            ctrl$dry <- which(ctrl$P==0)

            scen$wet <- which(scen$P>0)
            scen$dry <- which(scen$P==0)

            ctrl$pwet <- length(ctrl$P[which(ctrl$P>0)])/length(ctrl$P)
            scen$pwet <- length(scen$P[which(scen$P>0)])/length(scen$P)

            obs$margP <- .fitMarginal(obs$P[obs$wet],type="gamma")
            obs$margTW <- .fitMarginal(obs$T[obs$wet],type="norm")
            obs$margTD <- .fitMarginal(obs$T[obs$dry],type="norm")

            ctrl$margP <- .fitMarginal(ctrl$P[ctrl$wet],type="gamma")
            ctrl$margTW <- .fitMarginal(ctrl$T[ctrl$wet],type="norm")
            ctrl$margTD <- .fitMarginal(ctrl$T[ctrl$dry],type="norm")

            fit.normO <- copula::fitCopula(normalCopula(dim=2),pobs(matrix(c(obs$T[obs$wet],obs$P[obs$wet]),ncol=2)),
                                   method="itau",start=0,lower=NULL,upper=NULL,
                                   optim.control=list(maxit=100))
            obs$cpar1 <- fit.normO@copula@parameters

            fit.normC <- copula::fitCopula(normalCopula(dim=2),pobs(matrix(c(ctrl$T[ctrl$wet],ctrl$P[ctrl$wet]),ncol=2)),
                                   method="itau",start=0,lower=NULL,upper=NULL,
                                   optim.control=list(maxit=100))
            ctrl$cpar1 <- fit.normC@copula@parameters

            fit.normF <- copula::fitCopula(normalCopula(dim=2),pobs(matrix(c(scen$T[scen$wet],scen$P[scen$wet]),ncol=2)),
                                   method="itau",start=0,lower=NULL,upper=NULL,
                                   optim.control=list(maxit=100))

            scen$cpar1 <- fit.normF@copula@parameters

            if(cond=="P"){
              print(cond)
              adj <- .adjPT(obs,ctrl,scen)
            } else {
              print(cond)
              adj <- .adjTP(obs,ctrl,scen)
            }

            fit.normCor <- copula::fitCopula(normalCopula(dim=2),matrix(c(adj$T[scen$wet],adj$P[scen$wet]),ncol=2),
                                     method="itau",start=NULL,lower=NULL,upper=NULL,
                                     optim.method="BFGS",optim.control=list(maxit=1000))


#            adj$threshold <- sort(adj$P)[floor(length(adj$P)*(1-obs$pwet))]
#            adj$pwet <- length(adj$P[which(adj$P>adj$threshold)])/length(adj$P)
            adj$cpar1 <- fit.normCor@copula@parameters
            print(c(obs$cpar1,ctrl$cpar1,scen$cpar1,adj$cpar1))
            .Object@adj@data[,1] <- adj$T
            .Object@adj@data[,2] <- adj$P
            .Object@bc.attributes[["threshold"]] <- threshold
            .Object@bc.attributes[["cond"]] <- cond
            return(.Object)
          })
