#' @include GenericMethods.R Auxfunctions.R TimeSeries-class.R BiascoTimeSeries-class.R

#' export BC.ratio
BC.ratio <- setClass("RatioCorrection",
         contains = "BiascoTimeSeries",
         slots = c(obs = "TimeSeries", ctrl = "TimeSeries", scen = "TimeSeries")
#         prototype = prototype(obs = new("TimeSeries"), ctrl = new("TimeSeries"), scen = new("TimeSeries"))
)

setMethod(f = "initialize",
          signature = "RatioCorrection",
          function(.Object, ..., obs = new("TimeSeries"), ctrl = new("TimeSeries"),
                   scen = new("TimeSeries")){
            .Object <- callNextMethod(.Object, ...)
            .Object@obs <- obs
            .Object@ctrl <- ctrl
            .Object@scen <- scen
            validObject(.Object)
            return(.Object)
          })

#getters
setMethod(f = "obs",
          signature = "RatioCorrection",
          definition = function(object){
            return(object@obs)
          })

setMethod(f = "ctrl",
          signature = "RatioCorrection",
          definition = function(object){
            return(object@ctrl)
          })

setMethod(f = "scen",
          signature = "RatioCorrection",
          definition = function(object){
            return(object@scen)
          })

#Setters
setMethod(f = "obs<-",
          signature = "RatioCorrection",
          definition = function(object, value){
            object@obs <- value
            validObject(object)
            return(object)
          })

setMethod(f = "ctrl<-",
          signature = "RatioCorrection",
          definition = function(object, value){
            object@ctrl <- value
            validObject(object)
            return(object)
          })

setMethod(f = "scen<-",
          signature = "RatioCorrection",
          definition = function(object, value){
            object@scen <- value
            validObject(object)
            return(object)
          })

setMethod(f = "show",
          signature = "RatioCorrection",
          definition = function(object){
            cat("*** Class BiascoTimeSeries, method show *** \n")
            cat("* method = "); print (object@method)
            cat("* bc.attributes = "); print (object@bc.attributes)
            cat("* adj@nvar = "); print (object@adj@nvar)
            cat("* adj@Dim = "); print (object@adj@Dim)
            cat("* adj@dat (limited to first 100 values) = \n")
            if(length(object@adj@dat)!=0){
              print(formatC(object@adj@dat[1:min(length(object@adj@dat),100)]),quote=FALSE)
            }else{}
            cat("* obs@nvar = "); print (object@obs@nvar)
            cat("* obs@Dim = "); print (object@obs@Dim)
            cat("* obs@dat (limited to first 100 values) = \n")
            if(length(object@obs@dat)!=0){
              print(formatC(object@obs@dat[1:min(length(object@obs@dat),100)]),quote=FALSE)
            }else{}
            cat("* ctrl@nvar = "); print (object@ctrl@nvar)
            cat("* ctrl@Dim = "); print (object@ctrl@Dim)
            cat("* ctrl@dat (limited to first 100 values) = \n")
            if(length(object@ctrl@dat)!=0){
              print(formatC(object@ctrl@dat[1:min(length(object@ctrl@dat),100)]),quote=FALSE)
            }else{}
            cat("* scen@nvar = "); print (object@scen@nvar)
            cat("* scen@Dim = "); print (object@scen@Dim)
            cat("* scen@dat (limited to first 100 values) = \n")
            if(length(object@scen@dat)!=0){
              print(formatC(object@scen@dat[1:min(length(object@scen@dat),100)]),quote=FALSE)
            }else{}
            cat("******* End Show (BiascoTimeSeries) ******* \n")
          })

setMethod(".DcMean","ANY",function(.Object, ...) print("Missing or wrong input"))
setMethod(".DcMean",
          signature = "RatioCorrection",
          definition = function(.Object, ratio.max = 5){

            mean.obs <- mean(.Object@obs@dat)
            mean.ctrl <- mean(.Object@ctrl@dat)
            mean.scen <- mean(.Object@scen@dat)

            a <- min(mean(.Object@scen@dat,na.rm = T)/mean(.Object@ctrl@dat,na.rm = T), ratio.max)

            .Object@adj@dat <- .Object@obs@dat*a
            .Object@adj@Dim <- length(.Object@adj@dat)
            .Object@bc.attributes <- list("ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)
})

setMethod(".DcMeanSd1","ANY",function(.Object, ...) print("Missing or wrong input"))
setMethod(".DcMeanSd1",
          signature = "RatioCorrection",
          definition = function(.Object, nseq = 1, ratio.max = 5){

            m.obs <- mean(.Object@obs@dat)
            m.ctrl <- mean(.Object@ctrl@dat)
            m.scen <- mean(.Object@scen@dat)

            sd.obs <- sd(.Object@obs@dat)
            sd.ctrl <- sd(.Object@ctrl@dat)
            sd.scen <- sd(.Object@scen@dat)

            m.target <- m.obs*min(m.scen/m.ctrl, ratio.max)
            sd.target <- sd.obs*min(sd.scen/sd.ctrl, ratio.max)

            .Object@adj@dat <- .adjustEs(.Object@obs@dat, m.target, sd.target)
            .Object@adj@Dim <- length(.Object@adj@dat)
            .Object@bc.attributes <- list("nseq" = nseq, "ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)
          })


setMethod(".DcMeanSd2","ANY",function(.Object, ...) print("Missing or wrong input"))
setMethod(".DcMeanSd2",
          signature = "RatioCorrection",
          definition = function(.Object, ratio.max = 5){

          m.obs <- mean(.Object@obs@dat)
          m.ctrl <- mean(.Object@ctrl@dat)
          m.scen <- mean(.Object@scen@dat)

          sd.obs <- sd(.Object@obs@dat)
          sd.ctrl <- sd(.Object@ctrl@dat)
          sd.scen <- sd(.Object@scen@dat)

          m.targ <- m.obs*min(m.scen/m.ctrl, ratio.max)
          sd.targ <- sd.obs*min(sd.scen/sd.ctrl, ratio.max)

          .Object@adj@dat <- .iterAb(.Object@obs@dat, m.targ, sd.targ)
          .Object@adj@Dim <- length(.Object@adj@dat)
          .Object@bc.attributes <- list("ratio.max" = ratio.max)
          validObject(.Object)
          return(.Object)
          })

setMethod(".DcQmEmpir","ANY",function(.Object, ...) print("Missing or wrong input"))
setMethod(".DcQmEmpir",
          signature = "RatioCorrection",
          definition = function(.Object, smooth=0.02, pre.adj = F, post.adj = F, eps = 1e-5){

            #1. Initialisations
            PObs2 <- array(NA,dim=c(.Object@obs@Dim))

            PObsSort <- array(NA,dim=c(.Object@obs@Dim))
            PCtrlSort <- array(NA,dim=c(.Object@ctrl@Dim))
            PScenSort <- array(NA,dim=c(.Object@scen@Dim))

            PCtrlSmooth <- array(NA,dim=c(.Object@ctrl@Dim))
            PScenSmooth <- array(NA,dim=c(.Object@scen@Dim))

            #2. Insert random numbers to the time series
            PScenSort <- .Object@scen@dat + eps*runif(.Object@scen@Dim,0,1)
            PCtrlSort <- .Object@ctrl@dat + eps*runif(.Object@ctrl@Dim,0,1)
            PObs2 <- .Object@obs@dat + eps*runif(.Object@obs@Dim,0,1)

            PObs2 <- .zeros(PObs2,eps)

            #3. Order zero values evenly, sort obs and ctrl time series
            PScenSort <- sort(.zeros(PScenSort,eps))
            PCtrlSort <- sort(.zeros(PCtrlSort,eps))
            PObsSort <- sort(PObs2)

            if(pre.adj){
              mse.adj <- .minimizeMSE(PObsSort,PCtrlSort,PScenSort)
              PCtrlSort <- mse.adj$cal
              PScenSort <- mse.adj$val
            }

            PCtrlSmooth <- .smoothQuantiles(PCtrlSort, smooth)
            PScenSmooth <- .smoothQuantiles(PScenSort, smooth)

            .Object@adj@dat <- as.numeric(.quantMapP(PObs2,PCtrlSmooth,PScenSmooth))
            .Object@adj@dat[which(.Object@adj@dat < eps)] <- 0.0

            if(post.adj){
              m.ctrl <- mean(.Object@ctrl@dat)
              m.scen <- mean(.Object@scen@dat)
              if(m.ctrl > 0){
                Ratio1 <- m.scen/m.ctrl
              }
              else{
                Ratio1 <- 1.
              }
              m.adj <- mean(.Object@adj@dat)
              m.obs <- mean(.Object@obs@dat)
              if(m.obs>0){
                Ratio2 <- m.adj/m.obs
              }
              else{
                Ratio2 <- 1.
              }

              .Object@adj@dat <- .Object@adj@dat*Ratio1/Ratio2
            }
            .Object@adj@Dim <- length(.Object@adj@Dim)
            .Object@bc.attributes <- list("smooth" = smooth, "pred.adj" = pre.adj, "post.adj" = post.adj)
            validObject(.Object)
            return(.Object)
          })

#---------------------------------------
#BC
#---------------------------------------

setMethod(".BcMean","ANY",function(.Object, ...) print("Missing or wrong input"))
setMethod(".BcMean",
          signature = "RatioCorrection",
          definition = function(.Object, ratio.max = 5){

            mean.obs <- mean(.Object@obs@dat)
            mean.ctrl <- mean(.Object@ctrl@dat)
            mean.scen <- mean(.Object@scen@dat)

            a <- min(mean(.Object@obs@dat,na.rm = T)/mean(.Object@ctrl@dat,na.rm = T), ratio.max)

            .Object@adj@dat <- .Object@scen@dat*a
            .Object@adj@Dim <- length(.Object@adj@dat)
            .Object@bc.attributes <- list("ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)
          })


setMethod(".BcMeanSd1","ANY",function(.Object, ...) print("Missing or wrong input"))
setMethod(".BcMeanSd1",
          signature = "RatioCorrection",
          definition = function(.Object, nseq = NULL, ratio.max = 5){

            m.obs <- mean(.Object@obs@dat)
            m.ctrl <- mean(.Object@ctrl@dat)
            m.scen <- mean(.Object@scen@dat)

            sd.obs <- sd(.Object@obs@dat)
            sd.ctrl <- sd(.Object@ctrl@dat)
            sd.scen <- sd(.Object@scen@dat)

            m.target <- m.scen*min(m.obs/m.ctrl, ratio.max)
            sd.target <- sd.scen*min(sd.obs/sd.ctrl, ratio.max)

            .Object@adj@dat <- .adjustEs(.Object@scen@dat, m.target, sd.target)

            .Object@adj@Dim <- length(.Object@adj@dat)
            .Object@bc.attributes <- list("ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)

          })

setMethod(".BcMeanSd2","ANY",function(.Object, ...) print("Missing or wrong input"))
setMethod(".BcMeanSd2",
          signature = "RatioCorrection",
          definition = function(.Object, ratio.max = 5){

            m.obs <- mean(.Object@obs@dat)
            m.ctrl <- mean(.Object@ctrl@dat)
            m.scen <- mean(.Object@scen@dat)

            sd.obs <- sd(.Object@obs@dat)
            sd.ctrl <- sd(.Object@ctrl@dat)
            sd.scen <- sd(.Object@scen@dat)

            m.targ <- m.scen*min(m.obs/m.ctrl, ratio.max)
            sd.targ <- sd.scen*min(sd.obs/sd.ctrl, ratio.max)

            .Object@adj@dat <- .iterAb(.Object@scen@dat, m.targ, sd.targ)
            .Object@adj@Dim <- length(.Object@adj@dat)
            .Object@bc.attributes <- list("ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)
          })

setMethod(".BcQmParam","ANY",function(.Object, ...) print("Missing or wrong input"))
setMethod(".BcQmParam",
          signature = "RatioCorrection",
          definition = function(.Object, q.point = 0.95, threshold = 0.1){

  eps <- 1e-10
  obs <- list(P=.Object@obs@dat,
              margT=NA,margP=NA,pwet=NA,
              wet=NA,dry=NA)
  ctrl <- list(P=.Object@ctrl@dat,
               margT=NA,margP=NA,pwet=NA,
               wet=NA,dry=NA)
  scen <- list(P=.Object@scen@dat,
               margT=NA,margP=NA,pwet=NA,
               wet=NA,dry=NA)
  adj <- list(P=array(NA,dim=c(length(scen$P))),
              margT=NA,margP=NA,pwet=NA,
              wet=NA,dry=NA)

  tmp <- .ppRain(ref = obs$P, adj = ctrl$P, ref.threshold = threshold)
  ctrl$P <- tmp$adj

  tmp <- .ppRain(ref = obs$P, adj = scen$P, ref.threshold = threshold)
  scen$P <- tmp$adj

  #TODO handle rare cases, where no rain is observed or simulated
  # if(tmp$nP>=length(obs$P)){
  #   return(list(P=scen$P,T=scen$T))
  # }

  #Find locations of wet and dry days
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

  p1 <- pgamma(scen$P[scen$wet],shape=ctrl$margP$estimate[1],rate=ctrl$margP$estimate[2])
  p1[which(p1==1)] <- p1[which(p1==1)]-eps
  adj$P[scen$wet] <- qgamma(p1,shape=obs$margP$estimate[1],rate=obs$margP$estimate[2])
  adj$P[scen$dry] <- 0.0

  .Object@adj@dat <- as.numeric(adj$P)
  .Object@bc.attributes <- list("threshold" = threshold)
  return(.Object)

})

setMethod(".BcQmEmpir","ANY",function(.Object, ...) print("Missing or wrong input"))
setMethod(".BcQmEmpir",
          signature = "RatioCorrection",
          definition = function(.Object, smooth=0.02, pre.adj = F, post.adj = F, eps = 1e-5){

  #1. Initialisations
  PScen2 <- array(NA,dim=c(.Object@scen@Dim))

  PObsSort <- array(NA,dim=c(.Object@obs@Dim))
  PCtrlSort <- array(NA,dim=c(.Object@ctrl@Dim))
  PScenSort <- array(NA,dim=c(.Object@scen@Dim))

  PCtrlSmooth <- array(NA,dim=c(.Object@ctrl@Dim))
  PObsSmooth <- array(NA,dim=c(.Object@obs@Dim))

  #2. Insert random numbers to the time series
  PObsSort <- .Object@obs@dat + eps*runif(.Object@obs@Dim,0,1)
  PCtrlSort <- .Object@ctrl@dat + eps*runif(.Object@ctrl@Dim,0,1)
  PScen2 <- .Object@scen@dat + eps*runif(.Object@scen@Dim,0,1)

  PScen2 <- .zeros(PScen2,eps)

  #3. Order zero values evenly, sort obs and ctrl time series
  PObsSort <- sort(.zeros(PObsSort,eps))
  PCtrlSort <- sort(.zeros(PCtrlSort,eps))
  PScenSort <- sort(PScen2)

  if(pre.adj){
    .minimizeMSE(PObsSort,PCtrlSort,PScenSort)
    PCtrlSort <- mse.min$cal
    PScenSort <- mse.min$val
  }

  #6. Smoothing the time series
  PCtrlSmooth <- .smoothQuantiles(PCtrlSort, smooth)
  PObsSmooth <- .smoothQuantiles(PObsSort, smooth)

  .Object@adj@dat <- as.numeric(.quantMapP(PScen2,PCtrlSmooth,PObsSmooth))
  .Object@adj@dat[which(.Object@adj@dat < eps)] <- 0.0

  if(post.adj){
    m.ctrl <- mean(.Object@ctrl@dat)
    m.scen <- mean(.Object@scen@dat)
    if(m.ctrl > 0){
      Ratio1 <- m.scen/m.ctrl
    }
    else{
      Ratio1 <- 1.
    }
    m.adj <- mean(.Object@adj@dat)
    m.obs <- mean(.Object@obs@dat)
    if(m.obs>0){
      Ratio2 <- m.adj/m.obs
    }
    else{
      Ratio2 <- 1.
    }

    .Object@adj@dat <- .Object@adj@dat*Ratio1/Ratio2
  }
  .Object@adj@Dim <- length(.Object@adj@Dim)
  .Object@bc.attributes <- list("smooth" = smooth, "pred.adj" = pre.adj, "post.adj" = post.adj)
  validObject(.Object)
  return(.Object)
})

