#' @include TimeSeriesTP-class.R BiascoTimeSeriesTP-class.R
BC.joint <- setClass("JointCorrection",
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
          definition = function(.Object, cond="P", threshold = 0.1, jittering = T, jitter.amount = 1e-5, separate = F){

            #Create data frames for the data. Update these data frames whenever necessary.
            obs <- list(T=.Object@obs@data[,1],
                        P=.Object@obs@data[,2],
                        margT=NA,margP=NA,pwet=NA,
                        wet=NA,dry=NA,cpar1=NA,cpar2=NA)
            ctrl <- list(T=.Object@ctrl@data[,1],
                         P=.Object@ctrl@data[,2],
                         margT=NA,margP=NA,pwet=NA,
                         wet=NA,dry=NA,cpar1=NA)
            scen <- list(T=.Object@scen@data[,1],
                         P=.Object@scen@data[,2],
                         margT=NA,margP=NA,pwet=NA,
                         wet=NA,dry=NA,cpar1=NA)
            adj <- list(T=array(NA,dim=c(length(.Object@scen@data[,1]))),
                        P=array(NA,dim=c(length(.Object@scen@data[,2]))),
                        margT=NA,margP=NA,pwet=NA,
                        wet=NA,dry=NA,cpar1=NA)
            
            if(jittering){
              obs$T <- jitter(obs$T, amount = jitter.amount)
              obs$P <- obs$P+runif(length(obs$P))*jitter.amount
              ctrl$T <- jitter(ctrl$T, amount = jitter.amount)
              ctrl$P <- ctrl$P+runif(length(ctrl$P))*jitter.amount
              scen$T <- jitter(scen$T, amount = jitter.amount)
              scen$P <- scen$P+runif(length(scen$P))*jitter.amount
            }
            
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
            ctrl$margP <- .fitMarginal(ctrl$P[ctrl$wet],type="gamma")
            scen$margP <- .fitMarginal(scen$P[scen$wet],type="gamma")
            
            if(separate){
              obs$margTW <- .fitMarginal(obs$T[obs$wet],type="norm")
              obs$margTD <- .fitMarginal(obs$T[obs$dry],type="norm")
              
              ctrl$margTW <- .fitMarginal(ctrl$T[ctrl$wet],type="norm")
              ctrl$margTD <- .fitMarginal(ctrl$T[ctrl$dry],type="norm")
            }else{
              obs$margTW <- obs$margTD <- .fitMarginal(obs$T,type="norm")
              ctrl$margTW <- ctrl$margTD <- .fitMarginal(ctrl$T,type="norm")
            }

            print(obs$margP)
            print(ctrl$margP)
            print(scen$margP)
            
            fit.normO <- copula::fitCopula(copula::normalCopula(dim=2),copula::pobs(matrix(c(obs$T[obs$wet],obs$P[obs$wet]),ncol=2)),
                                   method="itau",start=0,lower=NULL,upper=NULL,
                                   optim.control=list(maxit=100))
            obs$cpar1 <- fit.normO@copula@parameters
            
            fit.normC <- copula::fitCopula(copula::normalCopula(dim=2),copula::pobs(matrix(c(ctrl$T[ctrl$wet],ctrl$P[ctrl$wet]),ncol=2)),
                                   method="itau",start=0,lower=NULL,upper=NULL,
                                   optim.control=list(maxit=100))
            ctrl$cpar1 <- fit.normC@copula@parameters

            # fit.normF <- copula::fitCopula(copula::normalCopula(dim=2),copula::pobs(matrix(c(scen$T[scen$wet],scen$P[scen$wet]),ncol=2)),
            #                        method="itau",start=0,lower=NULL,upper=NULL,
            #                        optim.control=list(maxit=100))
            # 
            # scen$cpar1 <- fit.normF@copula@parameters

            mvd.o <- copula::mvdc(copula = copula::ellipCopula(family = "normal", param = 0),
                          margins = c("norm", "gamma"), 
                          paramMargins = list(list(mean = obs$margTW$estimate[1], sd = obs$margTW$estimate[2]),
                                              list(shape = obs$margP$estimate[1], rate = obs$margP$estimate[2])))
            start <- as.vector(c(obs$margTW$estimate[1],obs$margTW$estimate[2],obs$margP$estimate[1],obs$margP$estimate[2],fit.normO@copula@parameters))
            fit.mvdO <- suppressMessages(copula::fitMvdc(matrix(c(obs$T[obs$wet],obs$P[obs$wet]),ncol=2), mvd.o, method = "Nelder-Mead",
                                                 start=start, optim.control=list(trace = 0, reltol = 1e-4, maxit=100),
                                                 lower = c(-100,0,0,0,-1), upper = c(100,inf,inf,inf,1)))
            
#            obs$cpar1 <- coef(fit.mvdO)[5]
#            obs$margTW$estimate[1] <- coef(fit.mvdO)[1]
#            obs$margTW$estimate[2] <- coef(fit.mvdO)[2]
#            obs$margP$estimate[1] <- coef(fit.mvdO)[3]
#            obs$margP$estimate[2] <- coef(fit.mvdO)[4]

            mvd.c <- copula::mvdc(copula = copula::ellipCopula(family = "normal", param = 0),
                          margins = c("norm", "gamma"), 
                          paramMargins = list(list(mean = ctrl$margTW[1], sd = ctrl$margTW[2]),
                                              list(shape = ctrl$margP[1], rate = ctrl$margP[2])))
            start <- c(ctrl$margTW$estimate[1],ctrl$margTW$estimate[2],ctrl$margP$estimate[1],ctrl$margP$estimate[2],fit.normC@copula@parameters)
            fit.mvdC <- suppressMessages(copula::fitMvdc(matrix(c(ctrl$T[ctrl$wet],ctrl$P[ctrl$wet]),ncol=2), mvd.c, method = "Nelder-Mead",
                                                 start=start, optim.control=list(trace = 0, reltol = 1e-4, maxit=1000)))
            
            # ctrl$cpar1 <- coef(fit.mvdC)[5]
            # ctrl$margTW$estimate[1] <- coef(fit.mvdC)[1]
            # ctrl$margTW$estimate[2] <- coef(fit.mvdC)[2]
            # ctrl$margP$estimate[1] <- coef(fit.mvdC)[3]
            # ctrl$margP$estimate[2] <- coef(fit.mvdC)[4]
            
            if(cond=="P"){
              print(cond)
              adj <- .adjPT(obs,ctrl,scen)
            } else {
              print(cond)
              adj <- .adjTP(obs,ctrl,scen)
            }
            fit.normCor <- copula::fitCopula(copula::normalCopula(dim=2),copula::pobs(matrix(c(adj$T[scen$wet],adj$P[scen$wet]),ncol=2)),
                                     method="itau",start=NULL,lower=NULL,upper=NULL,
                                     optim.method="BFGS",optim.control=list(maxit=1000))
            
            adj$cpar1 <- fit.normCor@copula@parameters
            .Object@adj@data[,1] <- adj$T
            .Object@adj@data[,2] <- adj$P
            .Object@bc.attributes[["threshold"]] <- threshold
            .Object@bc.attributes[["cond"]] <- cond
            .Object@bc.attributes[["copula.param"]] <- adj$cpar1
            .Object@method <- "J1"
            return(.Object)
          })
