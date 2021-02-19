# BCUH provides tools to bias adjust simulated time series with a number of uni- and multivariate methods.
# Copyright (C) 2018 Olle Räty
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundataion, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

#' @include GenericMethods.R Auxfunctions.R TimeSeries-class.R BiascoTimeSeries-class.R

#' An S4 class for an object of type AbsCorrection
#'
#' @slot obs an object of type TimeSeries containing the observed data
#' @slot ctrl an object of type TimeSeries containing the control period data for the simulation to be adjusted
#' @slot scen an object of type TimeSeries containing the scenario period data for the simulation to be adjusted
#' export BC.abs
BC.abs <- setClass("AbsCorrection",
          contains = "BiascoTimeSeries",
          slots = c(obs = "TimeSeries", ctrl = "TimeSeries", scen = "TimeSeries")
#          prototype = prototype(obs = new("TimeSeries"), ctrl = new("TimeSeries"), scen = new("TimeSeries"))
)

#' initialize
#' @rdname initialize
setMethod(f = "initialize",
          signature = "AbsCorrection",
          function(.Object, ..., obs = new("TimeSeries"), ctrl = new("TimeSeries"),
                   scen = new("TimeSeries")){
            .Object <- callNextMethod(.Object, ...)
            .Object@obs <- obs
            .Object@ctrl <- ctrl
            .Object@scen <- scen
            validObject(.Object)
            return(.Object)
          })

#' obs
#' @rdname obs
setMethod(f = "obs",
          signature = "AbsCorrection",
          definition = function(object){
            return(object@obs)
          })

#' ctrl
#' @rdname ctrl
setMethod(f = "ctrl",
          signature = "AbsCorrection",
          definition = function(object){
            return(object@ctrl)
          })

#' scen
#' @rdname scen
setMethod(f = "scen",
          signature = "AbsCorrection",
          definition = function(object){
            return(object@scen)
          })

#' obs<-
#' @rdname obs<-
setMethod(f = "obs<-",
          signature = "AbsCorrection",
          definition = function(object, value){
            .Object@obs <- value
            validObject(object)
            return(object)
          })

#' ctrl<-
#' @rdname ctrl<-
setMethod(f = "ctrl<-",
          signature = "AbsCorrection",
          definition = function(object, value){
            object@ctrl <- value
            validObject(object)
            return(object)
          })

#' scen<-
#' @rdname scen<-
setMethod(f = "scen<-",
          signature = "AbsCorrection",
          definition = function(object, value){
            object@scen <- value
            validObject(object)
            return(object)
          })

#' show
#' @rdname show
setMethod(f = "show",
          signature = "AbsCorrection",
          definition = function(object){
            cat("*** Class AbsCorrection object, method show *** \n")
            cat("* method = "); print (object@method)
            cat("* bc.attributes = "); print (object@bc.attributes)
            cat("* adj@nvar = "); print (object@adj@nvar)
            cat("* adj@Dim = "); print (object@adj@Dim)
            cat("* adj@data (limited to first 100 values) = \n")
            if(length(object@adj@data)!=0){
              print(formatC(object@adj@data[1:min(length(object@adj@data),100)]),quote=FALSE)
            }else{}
            cat("* obs@nvar = "); print (object@obs@nvar)
            cat("* obs@Dim = "); print (object@obs@Dim)
            cat("* obs@data (limited to first 100 values) = \n")
            if(length(object@obs@data)!=0){
              print(formatC(object@obs@data[1:min(length(object@obs@data),100)]),quote=FALSE)
            }else{}
            cat("* ctrl@nvar = "); print (object@ctrl@nvar)
            cat("* ctrl@Dim = "); print (object@ctrl@Dim)
            cat("* ctrl@data (limited to first 100 values) = \n")
            if(length(object@ctrl@data)!=0){
              print(formatC(object@ctrl@data[1:min(length(object@ctrl@data),100)]),quote=FALSE)
            }else{}
            cat("* scen@nvar = "); print (object@scen@nvar)
            cat("* scen@Dim = "); print (object@scen@Dim)
            cat("* scen@data (limited to first 100 values) = \n")
            if(length(object@scen@data)!=0){
              print(formatC(object@scen@data[1:min(length(object@scen@data),100)]),quote=FALSE)
            }else{}
            cat("******* End Show (AbsCorrectionObject) ******* \n")
          })

#' @rdname .DcMean
setMethod(".DcMean","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .DcMean
setMethod(".DcMean",
          signature = "AbsCorrection",
          definition = function(.Object){

            a <- mean(.Object@scen@data,na.rm=T) - mean(.Object@ctrl@data,na.rm=T)

            .Object@adj@data <- .Object@obs@data + a
            .Object@bc.attributes <- list("mean.factor"=a)
            .Object@adj@Dim <- length(.Object@adj@data)
            validObject(.Object)
            return(.Object)
          })

#' @rdname .DcMeanSd
setMethod(".DcMeanSd","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .DcMeanSd
setMethod(".DcMeanSd",
          signature = "AbsCorrection",
          definition = function(.Object, nseq=1){

            m.obs <- mean(.Object@obs@data, na.rm = T)
            m.ctrl <- mean(.Object@ctrl@data, na.rm = T)
            m.scen <- mean(.Object@scen@data, na.rm = T)

            sd.obs <- sd(.Object@obs@data, na.rm = T)
            sd.ctrl <- sd(.Object@ctrl@data, na.rm = T)
            sd.scen <- sd(.Object@scen@data, na.rm = T)

            a <- m.scen - m.ctrl
            b <-  sd.scen/sd.ctrl

            .Object@adj@data <- m.obs + a + (.Object@obs@data-m.obs)*b
            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("mean.factor"=a,
                                          "sd.facor"=b,
                                          "nseq"=nseq)
            .Object@adj@Dim <- length(.Object@adj@data)
            return(.Object)
          })

#' @rdname .DcMeanSdSkew
setMethod(".DcMeanSdSkew","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .DcMeanSdSkew
setMethod(".DcMeanSdSkew",
          signature = "AbsCorrection",
          definition = function(.Object, nseq=1){

            m.obs <- mean(.Object@obs@data, na.rm = T)
            m.ctrl <- mean(.Object@ctrl@data, na.rm = T)
            m.scen <- mean(.Object@scen@data, na.rm = T)

            sd.obs <- sd(.Object@obs@data, na.rm = T)
            sd.ctrl <- sd(.Object@ctrl@data, na.rm = T)
            sd.scen <- sd(.Object@scen@data, na.rm = T)

            skew.obs <- moments::skewness(.Object@obs@data, na.rm = T)
            skew.ctrl <- moments::skewness(.Object@ctrl@data, na.rm = T)
            skew.scen <- moments::skewness(.Object@scen@data, na.rm = T)

            m.target <- m.obs + (m.scen - m.ctrl)
            sd.target <- sd.obs*(sd.scen/sd.ctrl)
            skew.target <- skew.obs + (skew.scen - skew.ctrl)

            .Object@adj@data <- .iterSkewness(.Object@obs@data,skew.target)

            m.adj <- mean(.Object@adj@data)
            sd.adj <- sd(.Object@adj@data)

            .Object@adj@data <- m.target + (.Object@adj@data-m.adj)*(sd.target/sd.adj)
            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("nseq"=nseq)
            .Object@adj@Dim <- length(.Object@adj@data)
            validObject(.Object)
            return(.Object)
          })

#' @rdname .DcQmParam
setMethod(".DcQmParam","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .DcQmParam
setMethod(".DcQmParam",
          signature = "AbsCorrection",
          definition = function(.Object, fit.type = "linear", eps = 1e-5){

#            obs.rand <- sort(.Object@obs@data + runif(.Object@obs@Dim)*eps)
            #To cope with different lengths of ctrl and scen, the linear model
            #is fitted to the empirical quantiles instead of raw values
            ctrl.rand <- quantile(.Object@ctrl@data + runif(.Object@ctrl@Dim)*eps,seq(0,1,1/1000), method = 7)
            scen.rand <- quantile(.Object@scen@data + runif(.Object@scen@Dim)*eps,seq(0,1,1/1000), method = 7)

            data <- data.frame(y=scen.rand,x=ctrl.rand)

            if(fit.type == "linear"){ #Linear regression
              lm.model <- lm(y ~ x, data = data)
              a0 <- lm.model$coefficients[1]
              a1 <- lm.model$coefficients[2]
              tmp <- array(NA,dim=length(.Object@adj@data))
              smaller <- which(.Object@obs@data < min(.Object@ctrl@data))
              larger <- which(.Object@obs@data > max(.Object@ctrl@data))
              other <- seq(1,.Object@obs@Dim)
              
              if(length(smaller) != 0){
                tmp[smaller] <- a0 + (a1-1)*min(.Object@ctrl@data) + .Object@obs@data[smaller]
              }
              
              if(length(larger) != 0){
                tmp[larger] <- a0 + (a1-1)*max(.Object@ctrl@data) + .Object@obs@data[larger]              
              }
              
              if(length(c(smaller,larger)) != 0) other <- other[-c(larger,smaller)]
              tmp[other] <- a0 + a1*.Object@obs@data[other]
              
              .Object@adj@data <- as.numeric(tmp)
              .Object@bc.attributes <- list("fit.type"=fit.type,
                                            fit.parameters=list(fit=lm.model))
            }else if(fit.type == "normal"){ #Parametric fit
              ctrl.param <- MASS::fitdistr(ctrl.rand,"normal")
              scen.param <- MASS::fitdistr(scen.rand,"normal")
              .Object@adj@data <- qnorm(pnorm(.Object@obs@data,ctrl.param$estimate[1],
                                              ctrl.param$estimate[2]),
                                        scen.param$estimate[1],scen.param$estimate[2])
              .Object@bc.attributes <- list("fit.type"=fit.type,
                                            fit.parameters=list(ctrl.param=ctrl.param,
                                                                scen.param=scen.param))
            }else{
              stop("Wrong parametric fit")
            }
            .Object@adj@Dim <- length(.Object@adj@data)
            validObject(.Object)
            return(.Object)
          }
)

#' @rdname .DcQmEmpir
setMethod(".DcQmEmpir","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .DcQmEmpir
setMethod(".DcQmEmpir",
          signature = "AbsCorrection",
          definition = function(.Object, smooth = 0.05, pre.adj = F, post.adj = F, eps = 1e-5){
            #function(TObs,TCtrl,TScen,qsmooth,ladj,verbose=F){

            #1. Initialisations
            # nObs <- .Objeclength(TObs)
            # nCtrl <- length(TCtrl)
            # nScen <- length(TScen)

            TObs2 <- array(NA,dim=c(.Object@obs@Dim))

            TCtrlSmooth <- array(NA,dim=c(.Object@ctrl@Dim))
            TScenSmooth <- array(NA,dim=c(.Object@scen@Dim))

            #2. Insert random numbers to the values
            TObs2 <- .Object@obs@data + eps*runif(.Object@obs@Dim,0,1)

            TScenSort <- sort(.Object@scen@data + eps*runif(.Object@scen@Dim,0,1))
            TCtrlSort <- sort(.Object@ctrl@data + eps*runif(.Object@ctrl@Dim,0,1))

            #3. Smoothing the time series
            ScenSum <- array(NA,dim=c(.Object@scen@Dim+1))
            CtrlSum <- array(NA,dim=c(.Object@ctrl@Dim+1))

            ScenSum[1] <- 0.
            CtrlSum[1] <- 0.

            for (i in 1:.Object@ctrl@Dim){
              CtrlSum[i+1] <- CtrlSum[i]+TCtrlSort[i]
            }
            for (i in 1:.Object@scen@Dim){
              ScenSum[i+1] <- ScenSum[i]+TScenSort[i]
            }

            nsmooth <- floor(.Object@ctrl@Dim*smooth)
            for (i in 2:(.Object@ctrl@Dim+1)){
              j1 <- max(i-nsmooth,2)-1
              j2 <- min(i+nsmooth,.Object@ctrl@Dim+1)
              TCtrlSmooth[i-1] <- (CtrlSum[j2]-CtrlSum[j1])/(j2-j1)
            }
            nsmooth <- floor(.Object@scen@Dim*smooth)
            for (i in 2:(.Object@scen@Dim+1)){
              j1 <- max(i-nsmooth,2)-1
              j2 <- min(i+nsmooth,.Object@scen@Dim+1)
              TScenSmooth[i-1] <- (ScenSum[j2]-ScenSum[j1])/(j2-j1)
            }

            #5. Find scenario period precipitation values
            #   corresponding to each observed precipitation value.
            #
            #   Logic: First find the value corresponding observations
            #   from cmodjs (smoothed and scaled version of the simulated
            #   baseline time series in the baseline period). Scenario
            #   value will then be obtained from the same location in
            #   table pmodjs. Below the lowest and above the highest value
            #   scenario value is extrapolated assuming constant relative change

            .Object@adj@data <- as.numeric(.quantMapT(TObs2,TCtrlSmooth,TScenSmooth))

            #  6. Scale the scenario distribution, so that the change in
            #  the time mean temperature will be 'right', ie. corresponding
            #  to the simulated time mean temperature change between the
            #  baseline and scenario periods.

            if(post.adj){
              MObs <- mean(.Object@obs@data)
              MCtrl <- mean(.Object@ctrl@data)
              MScen <- mean(.Object@scen@data)
              MScenObs <- mean(.Object@adj@data)

              diff <- (MScen-MCtrl)-(MScenObs-MObs)

              .Object@adj@data <- .Object@adj@data+diff
            }

            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("smooth" = smooth, "pre.adj" = pre.adj, "post.adj" = post.adj)
            validObject(.Object)
            return(.Object)
          })

#-----------------------------------------
#BC
#-----------------------------------------

#' @rdname .BcMean
setMethod(".BcMean","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .BcMean
setMethod(".BcMean",
          signature = "AbsCorrection",
          definition = function(.Object){

            a <- mean(.Object@obs@data,na.rm=T) - mean(.Object@ctrl@data,na.rm=T)

            .Object@adj@data <- .Object@scen@data + a
            .Object@bc.attributes <- list()
            .Object@adj@Dim <- length(.Object@adj@data)
            validObject(.Object)
            return(.Object)
          })

#' @rdname .BcMeanSd
setMethod(".BcMeanSd","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .BcMeanSd
setMethod(".BcMeanSd",
          signature = "AbsCorrection",
          definition = function(.Object,nseq=1){

            m.obs <- mean(.Object@obs@data, na.rm = T)
            m.ctrl <- mean(.Object@ctrl@data, na.rm = T)
            m.scen <- mean(.Object@scen@data, na.rm = T)

            sd.obs <- sd(.Object@obs@data, na.rm = T)
            sd.ctrl <- sd(.Object@ctrl@data, na.rm = T)
            sd.scen <- sd(.Object@scen@data, na.rm = T)

            a <- m.obs - m.ctrl
            b <-  sd.obs/sd.ctrl

            .Object@adj@data <- m.scen + a + (.Object@scen@data-m.scen)*b
            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("nseq" = nseq)
            return(.Object)
          })

#' @rdname .BcMeanSdSkew
setMethod(".BcMeanSdSkew","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .BcMeanSdSkew
setMethod(".BcMeanSdSkew",
          signature = "AbsCorrection",
          definition = function(.Object, nseq=1){

            m.obs <- mean(.Object@obs@data, na.rm = T)
            m.ctrl <- mean(.Object@ctrl@data, na.rm = T)
            m.scen <- mean(.Object@scen@data, na.rm = T)

            sd.obs <- sd(.Object@obs@data, na.rm = T)
            sd.ctrl <- sd(.Object@ctrl@data, na.rm = T)
            sd.scen <- sd(.Object@scen@data, na.rm = T)

            skew.obs <- moments::skewness(.Object@obs@data, na.rm = T)
            skew.ctrl <- moments::skewness(.Object@ctrl@data, na.rm = T)
            skew.scen <- moments::skewness(.Object@scen@data, na.rm = T)

            m.target <- m.scen + (m.obs - m.ctrl)
            sd.target <- sd.scen*(sd.obs/sd.ctrl)
            skew.target <- skew.scen + (skew.obs - skew.ctrl)

            .Object@adj@data <- .iterSkewness(.Object@scen@data,skew.target)

            m.adj <- mean(.Object@adj@data)
            sd.adj <- sd(.Object@adj@data)

            .Object@adj@data <- m.target + (.Object@adj@data-m.adj)*(sd.target/sd.adj)
            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("nseq" = nseq)
            validObject(.Object)
            return(.Object)
          })

#' @rdname .BcQmParam
setMethod(".BcQmParam","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .BcQmParam
setMethod(".BcQmParam",
          signature = "AbsCorrection",
          definition = function(.Object, fit.type = "linear", eps = 1e-5){

            #To cope with different lengths of ctrl and scen, the linear regression
            #is fitted to the empirical quantiles instead of raw values
            ctrl.rand <- quantile(.Object@ctrl@data + runif(.Object@ctrl@Dim)*eps,seq(0,1,1/1000),method=7)
            obs.rand <- quantile(.Object@obs@data + runif(.Object@obs@Dim)*eps,seq(0,1,1/1000),method=7)

            data <- data.frame(y=obs.rand,x=ctrl.rand)

            if(fit.type == "linear"){#Linear regression
              lm.model <- lm(y ~ x, data = data)
              a0 <- lm.model$coefficients[1]
              a1 <- lm.model$coefficients[2]
              tmp <- array(NA,dim=length(.Object@adj@data))
              smaller <- which(.Object@scen@data < min(.Object@ctrl@data))
              larger <- which(.Object@scen@data > max(.Object@ctrl@data))
              other <- seq(1,.Object@scen@Dim)
              
              if(length(smaller) != 0){
#                other <- other[-c(smaller)]
                tmp[smaller] <- a0 + (a1-1)*min(.Object@ctrl@data) + .Object@scen@data[smaller]
              }
              
              if(length(larger) != 0){
#                other <- other[-c(larger)]
                tmp[larger] <- a0 + (a1-1)*max(.Object@ctrl@data) + .Object@scen@data[larger]
              }
              
              if(length(c(smaller,larger) != 0)) other <- other[-c(smaller,larger)]
              tmp[other] <- a0 + a1*.Object@scen@data[other]
              .Object@adj@data <- as.numeric(tmp)
              .Object@bc.attributes <- list("fit.type"=fit.type,
                                            fit.parameters=list(fit=lm.model))
            }else if(fit.type == "normal"){#Parametric fit
              ctrl.param <- MASS::fitdistr(ctrl.rand,"normal")
              obs.param <- MASS::fitdistr(obs.rand,"normal")
              .Object@adj@data <- qnorm(pnorm(.Object@scen@data,ctrl.param$estimate[1],ctrl.param$estimate[2]),
                                        obs.param$estimate[1],obs.param$estimate[2])
              .Object@bc.attributes <- list("fit.type"=fit.type,
                                            fit.parameters=list(ctrl.param=ctrl.param,
                                                                scen.param=scen.param))
            }else{
              stop("Wrong parametric fit")
            }
            .Object@adj@Dim <- length(.Object@adj@data)
            validObject(.Object)
            return(.Object)
          }
)

#' @rdname .BcQmEmpir
setMethod(".BcQmEmpir","ANY",function(.Object) print("Missing or wrong input"))

#' @rdname .BcQmEmpir
setMethod(".BcQmEmpir",
          signature = "AbsCorrection",
          definition = function(.Object, smooth = 0.05, pre.adj = F, post.adj = F, eps = 1e-5){

            TObs2 <- array(NA,dim=c(.Object@obs@Dim))
            TCtrlSmooth <- array(NA,dim=c(.Object@ctrl@Dim))
            TObsSmooth <- array(NA,dim=c(.Object@obs@Dim))

            TScen2 <- .Object@scen@data + eps*runif(.Object@scen@Dim,0,1)
            TObsSort <- sort(.Object@obs@data + eps*runif(.Object@obs@Dim,0,1))
            TCtrlSort <- sort(.Object@ctrl@data + eps*runif(.Object@ctrl@Dim,0,1))

            ObsSum <- array(NA,dim=c(.Object@obs@Dim+1))
            CtrlSum <- array(NA,dim=c(.Object@ctrl@Dim+1))

            ObsSum[1] <- 0.
            CtrlSum[1] <- 0.

            for (i in 1:.Object@ctrl@Dim){
              CtrlSum[i+1] <- CtrlSum[i]+TCtrlSort[i]
            }
            for (i in 1:.Object@obs@Dim){
              ObsSum[i+1] <- ObsSum[i]+TObsSort[i]
            }

            nsmooth <- floor(.Object@ctrl@Dim*smooth)
            for (i in 2:(.Object@ctrl@Dim+1)){
              j1 <- max(i-nsmooth,2)-1
              j2 <- min(i+nsmooth,.Object@ctrl@Dim+1)
              TCtrlSmooth[i-1] <- (CtrlSum[j2]-CtrlSum[j1])/(j2-j1)
            }
            nsmooth <- floor(.Object@obs@Dim*smooth)
            for (i in 2:(.Object@obs@Dim+1)){
              j1 <- max(i-nsmooth,2)-1
              j2 <- min(i+nsmooth,.Object@obs@Dim+1)
              TObsSmooth[i-1] <- (ObsSum[j2]-ObsSum[j1])/(j2-j1)
            }

            .Object@adj@data <- as.numeric(.quantMapT(TScen2,TCtrlSmooth,TObsSmooth))

            if(post.adj){
              MObs <- mean(.Object@obs@data)
              MCtrl <- mean(.Object@ctrl@data)
              MScen <- mean(.Object@scen@data)
              MScenObs <- mean(.Object@adj@data)

              diff <- (MScen-MCtrl)-(MScenObs-MObs)

              .Object@adj@data <- .Object@adj@data+diff
            }

            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("smooth" = smooth, 
                                          "pre.adj" = pre.adj, 
                                          "post.adj" = post.adj)
            validObject(.Object)
            return(.Object)
          })
