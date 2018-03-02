
.BC.ratio <- setClass("RatioCorrection",
         contains = "BiascoTimeSeries",
         slots = c(obs = "TimeSeries", ctrl = "TimeSeries", scen = "TimeSeries"),
         prototype = prototype(obs = new("TimeSeries"), ctrl = new("TimeSeries"), scen = new("TimeSeries"))
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
            cat("******* End Show (BiascoTimeSeries) ******* \n")
          })

setMethod(".DcMean","ANY",function() print("Missing or wrong input"))
setMethod(".DcMean",
          signature = "RatioCorrection",
          definition = function(.Object,ratio.max=5){

            mean.obs <- mean(.Object@obs@data)
            mean.ctrl <- mean(.Object@ctrl@data)
            mean.scen <- mean(.Object@scen@data)

            a <- min(mean(.Object@scen@data,na.rm=T)/mean(.Object@ctrl@data,na.rm=T),ratio.max)

            .Object@adj@data <- .Object@obs@data*a
            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)
})

setMethod(".DcMeanSd1","ANY",function() print("Missing or wrong input"))
setMethod(".DcMeanSd1",
          signature = "RatioCorrection",
          definition = function(.Object, nseq=1, ratio.max=5){

            m.obs <- mean(.Object@obs@data)
            m.ctrl <- mean(.Object@ctrl@data)
            m.scen <- mean(.Object@scen@data)

            sd.obs <- sd(.Object@obs@data)
            sd.ctrl <- sd(.Object@ctrl@data)
            sd.scen <- sd(.Object@scen@data)

            m.target <- m.obs*min(m.scen/m.ctrl, ratio.max)
            sd.target <- sd.obs*min(sd.scen/sd.ctrl, ratio.max)

            .Object@adj@data <- .adjustEs(.Object@obs@data, m.target, sd.target)
            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("nseq" = nseq, "ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)
          })


setMethod(".DcMeanSd2","ANY",function() print("Missing or wrong input"))
setMethod(".DcMeanSd2",
          signature = "RatioCorrection",
          definition = function(.Object, ratio.max = 5){

          m.obs <- mean(.Object@obs@data)
          m.ctrl <- mean(.Object@ctrl@data)
          m.scen <- mean(.Object@scen@data)

          sd.obs <- sd(.Object@obs@data)
          sd.ctrl <- sd(.Object@ctrl@data)
          sd.scen <- sd(.Object@scen@data)

          m.targ <- m.obs*min(m.scen/m.ctrl, ratio.max)
          sd.targ <- sd.obs*min(sd.scen/sd.ctrl, ratio.max)

          .Object@adj@data <- .iterAb(.Object@obs@data, m.targ, sd.targ)
          .Object@adj@Dim <- length(.Object@adj@data)
          .Object@bc.attributes <- list("ratio.max" = ratio.max)
          validObject(.Object)
          return(.Object)
          })

setMethod(".DcQmEmpir","ANY",function() print("Missing or wrong input"))
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
            PScenSort <- .Object@scen@data + eps*runif(.Object@scen@Dim,0,1)
            PCtrlSort <- .Object@ctrl@data + eps*runif(.Object@ctrl@Dim,0,1)
            PObs2 <- .Object@obs@data + eps*runif(.Object@obs@Dim,0,1)

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

            .Object@adj@data <- as.numeric(.quantMapP(PObs2,PCtrlSmooth,PScenSmooth))
            .Object@adj@data[which(.Object@adj@data < eps)] <- 0.0

            if(post.adj){
              m.ctrl <- mean(.Object@ctrl@data)
              m.scen <- mean(.Object@scen@data)
              if(m.ctrl > 0){
                Ratio1 <- m.scen/m.ctrl
              }
              else{
                Ratio1 <- 1.
              }
              m.adj <- mean(.Object@adj@data)
              m.obs <- mean(.Object@obs@data)
              if(m.obs>0){
                Ratio2 <- m.adj/m.obs
              }
              else{
                Ratio2 <- 1.
              }

              .Object@adj@data <- .Object@adj@data*Ratio1/Ratio2
            }
            .Object@adj@Dim <- length(.Object@adj@Dim)
            .Object@bc.attributes <- list("smooth" = smooth, "pred.adj" = pre.adj, "post.adj" = post.adj)
            validObject(.Object)
            return(.Object)
          })

#---------------------------------------
#BC
#---------------------------------------

setMethod(".BcMean","ANY",function() print("Missing or wrong input"))
setMethod(".BcMean",
          signature = "RatioCorrection",
          definition = function(.Object,ratio.max=5){

            mean.obs <- mean(.Object@obs@data)
            mean.ctrl <- mean(.Object@ctrl@data)
            mean.scen <- mean(.Object@scen@data)

            a <- min(mean(.Object@obs@data,na.rm=T)/mean(.Object@ctrl@data,na.rm=T), ratio.max)

            .Object@adj@data <- .Object@scen@data*a
            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)
          })


setMethod(".BcMeanSd1","ANY",function() print("Missing or wrong input"))
setMethod(".BcMeanSd1",
          signature = "RatioCorrection",
          definition = function(.Object,nseq=NULL,ratio.max=5){

            m.obs <- mean(.Object@obs@data)
            m.ctrl <- mean(.Object@ctrl@data)
            m.scen <- mean(.Object@scen@data)

            sd.obs <- sd(.Object@obs@data)
            sd.ctrl <- sd(.Object@ctrl@data)
            sd.scen <- sd(.Object@scen@data)

            m.target <- m.scen*min(m.obs/m.ctrl, ratio.max)
            sd.target <- sd.scen*min(sd.obs/sd.ctrl, ratio.max)

            .Object@adj@data <- .adjustEs(.Object@scen@data, m.target, sd.target)

            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)

          })

setMethod(".BcMeanSd2","ANY",function() print("Missing or wrong input"))
setMethod(".BcMeanSd2",
          signature = "RatioCorrection",
          definition = function(.Object, ratio.max = 5){

            m.obs <- mean(.Object@obs@data)
            m.ctrl <- mean(.Object@ctrl@data)
            m.scen <- mean(.Object@scen@data)

            sd.obs <- sd(.Object@obs@data)
            sd.ctrl <- sd(.Object@ctrl@data)
            sd.scen <- sd(.Object@scen@data)

            m.targ <- m.scen*min(m.obs/m.ctrl, ratio.max)
            sd.targ <- sd.scen*min(sd.obs/sd.ctrl, ratio.max)

            .Object@adj@data <- .iterAb(.Object@scen@data, m.targ, sd.targ)
            .Object@adj@Dim <- length(.Object@adj@data)
            .Object@bc.attributes <- list("ratio.max" = ratio.max)
            validObject(.Object)
            return(.Object)
          })

setMethod(".BcQmParam","ANY",function() print("Missing or wrong input"))
setMethod(".BcQmParam",
          signature = "RatioCorrection",
          definition = function(.Object, q.point = 0.95, threshold = 0.1){

  eps <- 1e-10
  obs <- list(P=.Object@obs@data,
              margT=NA,margP=NA,pwet=NA,
              wet=NA,dry=NA)
  ctrl <- list(P=.Object@ctrl@data,
               margT=NA,margP=NA,pwet=NA,
               wet=NA,dry=NA)
  scen <- list(P=.Object@scen@data,
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

  .Object@adj@data <- as.numeric(adj$P)
  .Object@bc.attributes <- list("threshold" = threshold)
  return(.Object)

})

setMethod(".BcQmEmpir","ANY",function() print("Missing or wrong input"))
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
  PObsSort <- .Object@obs@data + eps*runif(.Object@obs@Dim,0,1)
  PCtrlSort <- .Object@ctrl@data + eps*runif(.Object@ctrl@Dim,0,1)
  PScen2 <- .Object@scen@data + eps*runif(.Object@scen@Dim,0,1)

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

  .Object@adj@data <- as.numeric(.quantMapP(PScen2,PCtrlSmooth,PObsSmooth))
  .Object@adj@data[which(.Object@adj@data < eps)] <- 0.0

  if(post.adj){
    m.ctrl <- mean(.Object@ctrl@data)
    m.scen <- mean(.Object@scen@data)
    if(m.ctrl > 0){
      Ratio1 <- m.scen/m.ctrl
    }
    else{
      Ratio1 <- 1.
    }
    m.adj <- mean(.Object@adj@data)
    m.obs <- mean(.Object@obs@data)
    if(m.obs>0){
      Ratio2 <- m.adj/m.obs
    }
    else{
      Ratio2 <- 1.
    }

    .Object@adj@data <- .Object@adj@data*Ratio1/Ratio2
  }
  .Object@adj@Dim <- length(.Object@adj@Dim)
  .Object@bc.attributes <- list("smooth" = smooth, "pred.adj" = pre.adj, "post.adj" = post.adj)
  validObject(.Object)
  return(.Object)
})

