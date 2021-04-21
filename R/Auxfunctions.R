#' Find correct skewness iteratively (internal)
#' @title .iterSkewness
#'
#' @param data data to be adjusted
#' @param skew.target target skewness
#' @param nseq data set is divided into nseq chunks on which the skewness is adjusted separately
#' @param tol threshold for iteration
#' @param a.lim limits for the exponent used to adjust the skewness
#' @param iter.max max. number of iterations
#'
#' @importFrom moments skewness
#'
#' @return Bias adjusted time series
#' @keywords internal
#'
.iterSkewness <- function(source, skew.target, nseq = 1, tol = 1e-3, 
                          a.lim = c(0.5,9), iter.max = 100){

  if(length(a.lim) < 2)stop("Missing one or both limits for the skewness scaling factor")
#  print(source)
  skew.source <- moments::skewness(source)

  if(skew.source < skew.target){
    b <- min(source)
    source <- source-b
    sign <- 1
  }else{
    skew.target <- -skew.target
    b <- max(source)
    source <- -(source-b)
    sign <- -1
  }

  i <- 1
  source.scaled <- source
  while(abs(skew.source-skew.target)>tol){
    a <- sum(a.lim)/2
    source.scaled <- source^a
    skew.source <- moments::skewness(source.scaled)
#    cat("skew.source: ",skew.source,"skew.target: ",skew.target,"a: ",a,"\n")
    if(skew.source<skew.target){
      a.lim[1] <- a
    }else{
      a.lim[2] <- a
    }
    if(i == iter.max){
      warning(paste("Skewness(",skew.source,") not converging to the target (",skew.target,")"))
      return(source.scaled)
    }
    i <- i + 1
  }
  
  source <- b + source.scaled*sign
  return(source)
}


#' Quantile mapping algorithm (internal)
#'
#'Determine the future temperatures (ScenCor) that
#  correspond to the observed temperatures (Obs).
#  In case that the baseline and scenario period include
#  a different number of days (nObs <> nCtrl):
#    iObs <- round(0.5+nObs*(iCtrl-0.5)/nCtrl)
#    iObs1 <- round(0.5+nObs*(iCtrl-1.5)/nCtrl)
#  ..so the idea is to find more equivalent quantiles.
#  Nevertheless, it is cleanest to use baseline and scenario
#  time series of the same length.
#'
#' @title .quantMapT
#' @param Scen scenario period time series
#' @param Ctrl control period time series
#' @param Obs control period observations
#'
#' @return Bias adjusted time series
#'
#' @keywords internal
#'
.quantMapT <- function(Scen,Ctrl,Obs){

  nScen <- length(Scen)
  nCtrl <- length(Ctrl)
  nObs <- length(Obs)

  ScenCor <- array(NA,dim=c(nScen))

  for (i in 1:nScen){
    if(Scen[i]>Ctrl[nCtrl]){
      iCtrl <- nCtrl+1
    }
    else{
      iCtrl <- min(which(Ctrl>=Scen[i]))
    }

    iObs <- round(0.5+nObs*(iCtrl-0.5)/nCtrl)
    iObs1 <- round(0.5+nObs*(iCtrl-1.5)/nCtrl)

    if(iCtrl==1){
      ScenCor[i] <- Scen[i] + (Obs[iObs]-Ctrl[iCtrl])
    }
    else if(iCtrl==(nCtrl+1)){
      ScenCor[i] <- Scen[i] + (Obs[iObs1]-Ctrl[iCtrl-1])

    }
    else
    {
      a <- (Scen[i]-Ctrl[iCtrl-1])/(Ctrl[iCtrl]-Ctrl[iCtrl-1])
      ScenCor[i] <- (1.-a)*Obs[iObs1]+a*Obs[iObs]
    }
  }
  return(ScenCor)
}


#' Handle zero values (internal)
#'
#' @title .zeros
#' @param Data Time series to be adjusted
#' @param eps small value used to adjust the zero values
#'
#' @return adjusted time series
#'
#' @keywords internal
#'
.zeros <- function(Data,eps){
  zero <- which(Data<eps & Data>=0)
  nzero <- length(zero)
  nsmaller <- array(NA,dim=c(nzero))

  for (i in 1:nzero){
    nsmaller[i] <- length(which(Data[zero[-i]]<Data[zero[i]]))
  }

  for (i in 1:nzero){
    Data[zero[i]] <- eps*(1+2.*nsmaller[i])/(2.*nzero)
  }
  return(Data)
}


#' Quantile mapping algorithm for precipitation (internal)
#'
#'Determine the future temperatures (ScenCor) that
#  correspond to the observed temperatures (Obs).
#  In case that the baseline and scenario period include
#  a different number of days (nObs <> nCtrl):
#    iObs <- round(0.5+nObs*(iCtrl-0.5)/nCtrl)
#    iObs1 <- round(0.5+nObs*(iCtrl-1.5)/nCtrl)
#  ..so the idea is to find more equivalent quantiles.
#  Nevertheless, it is cleanest to use baseline and scenario
#  time series of the same length.
#'
#' @title .quantMapP
#' @param Scen scenario period time series
#' @param Ctrl control period time series
#' @param Obs control period observations
#'
#' @return Adjusted time series for the Scenario period
#'
#' @keywords internal
#'
.quantMapP <- function(Scen,Ctrl,Obs){

  nScen <- length(Scen)
  nCtrl <- length(Ctrl)
  nObs <- length(Obs)

  ScenCor <- array(NA,dim=c(nScen))
  for (i in 1:nScen){
    if(Scen[i]>Ctrl[nCtrl]){
      iCtrl <- nCtrl+1
    }
    else{
      iCtrl <- min(which(Ctrl>=Scen[i]))
    }

    iObs <- round(0.5+nObs*(iCtrl-0.5)/nCtrl)
    iObs1 <- round(0.5+nObs*(iCtrl-1.5)/nCtrl)

    if(iCtrl==1){
      ScenCor[i] <- Scen[i]/Ctrl[iCtrl]*Obs[iObs]
    }
    else if(iCtrl==nCtrl+1){
      ScenCor[i] <- (Scen[i]/Ctrl[iCtrl-1])*Obs[iObs1]
    }
    else
    {
      a <- (Scen[i]-Ctrl[iCtrl-1])/(Ctrl[iCtrl]-Ctrl[iCtrl-1])
      ScenCor[i] <- (1.-a)*Obs[iObs1]+a*Obs[iObs]
    }

    if(Scen[i]>=0)ScenCor[i] <- max(ScenCor[i],0.)

  }

  return(ScenCor)

}

#' Engen-Skaugen algorithm (internal)
#'
#' Adjust the cv interatively using the Engen-Skaugen algorithm.
#'
#' @title .adjustEs
#' @param data Data to be adjusted
#' @param mean.target Traget value for the mean
#' @param sd.target Target value for standard deviation
#' @param tol Tolerance for the convergence
#' @param max.iter Max. number of iterations, after which the execution is stopped
#'
#' @return Adjusted time series for the scenario period.
#'
#' @keywords internal
#'
.adjustEs <- function(data, mean.target, sd.target, tol = 1e-3, max.iter = 100){

  data.unadj <- data
  mean.data <- mean(data)
  sd.data <- sd(data)

  i <- 1
  while(abs(mean.target/mean.data - 1) > tol | abs(sd.target/sd.data - 1) > tol ){
    data <- data*(mean.target/mean.data)
    sd.data <- (mean.target/mean.data)*sd.data
    data <- mean.target + (data - mean.target)*(sd.target/sd.data)
    data[which(data < 0)] <- 0.
    mean.data <- mean(data)
    sd.data <- sd(data)
    if(i == max.iter){
      warning("coefficient of variability not converging to the target value!")
      data <- data*mean.target/mean.data
      return(data)
    }
    i <- i + 1
  }
  return(data)
}

#' Power transformation algorithm (internal)
#'
#' Find iteratively the correct exponent in the power transformation of a ratio variable.
#'
#' @title .iterAb
#' @param data Input data
#' @param mean.target Target values for the mean
#' @param sd.target Target value for standard deviations
#' @param b.lim Limits for the exponent value
#' @param tol Tolerance for the convergence
#' @param max.iter  Max. number of iterations. Used to stop the execution,
#'   if the algorithm does not converge to the target values.
#'
#' @return Adjusted time series for the scenario period.
#'
#' @keywords internal
#'
.iterAb <- function(data, mean.target, sd.target, b.lim = c(0.001, 5), tol = 1e-3, max.iter = 100){
  b <- 1
  cv.target <- sd.target/max(mean.target,0.001)
  mean.data <- mean(data)
  sd.data <- sd(data)
  cv <- sd.data/mean.data

  data.unadj <- data

  i <- 1
  while(abs(cv/cv.target-1) > tol){

    if(cv<cv.target){
      b.lim[1] <- b
      b <- mean(b.lim)
    }else{
      b.lim[2] <- b
      b <- mean(b.lim)
    }
      data <- data.unadj^b

      mean.data <- mean(data)
      sd.data <- sd(data)
      cv <- sd.data/mean.data

      if(i == max.iter){
        warning("coefficient of variability not converging to the target value!")
        data <- data*mean.target/mean.data
        return(data)
      }
      i <- i + 1
    }

  data <- data*mean.target/mean.data
  return(data)
}

#' Pre-process precipitation data for parametric quantile mapping.
#'
#' @title .ppRain
#' @param ref Reference data
#' @param adj Data to be adjusted
#' @param ref.threshold Threshold for wet-day values.
#'
#' @return A list containing the number of dry days both in the reference and calibration data sets.
#'
#' @keywords internal
#'

.ppRain <- function(ref, adj, ref.threshold = 0.1){
  er <- ecdf(ref)
  ea <- ecdf(adj)
  pr <- er(0.1)
  pa <- ea(0.1)
  rdry <- ceiling(length(ref)*pr)
  adry <- ceiling(length(adj)*pr)
  adj.sorted <- sort(adj)
  adj.sorted.index <- sort(adj, index.return = TRUE)$ix
  ath <- adj.sorted[adry]
  if(pa>pr & !all(ref <= ref.threshold)){
    ref.sorted <- sort(ref)
    start <- ceiling(length(adj)*pr)
    end <- ceiling(length(adj)*pa)
    adj.sorted[start:end] <- runif(n=(end-start+1),
                                   min=ref.threshold,
                                   max=ref.sorted[end])
  }
  if(pa<pr) adj.sorted[1:ceiling(length(adj)*pr)] <- 0
  adj[adj.sorted.index] <- adj.sorted
  
  return(list("nDry" = c(rdry, adry), 
              "adj.threshold" = ath, "adj" = adj))
}


#' Minimize MSE between the reference and simulated time series (internal)
#'
#' @title .minimizeMSE
#' @param ref Reference data
#' @param cal Calibration data
#' @param val Scenario data
#'
#' @return Adjusted calibration and scenario data.
#'
#' @keywords internal
#'
.minimizeMSE <- function(ref,cal,val){
  omsum <- 0
  m2sum <- 0

  for (i in 1:length(cal)){
    j <- round(0.5+length(ref)*(i-0.5)/length(cal))
    omsum <- omsum + cal[i]*ref[j]
    m2sum <- m2sum + cal[i]^2
  }

  cal <- cal*omsum/m2sum
  val <- val*omsum/m2sum
  return(list(cal=cal,val=val))
}

.smoothQuantiles <- function(data, smooth) {
  smoothed <- array(NA,dim=c(length(data)))
  cum.data <- array(NA,dim=c(length(data)+1))
  cum.data[1] <- 0.
  for (i in 1:length(data)){
    cum.data[i+1] <- cum.data[i] + data[i]
  }

  nsmooth <- floor(length(data)*smooth)
  for (i in 2:(length(data)+1)){
    j1 <- max(i-nsmooth,2)-1
    j2 <- min(i+nsmooth,length(data)+1)
    smoothed[i-1] <- (cum.data[j2]-cum.data[j1])/(j2-j1)
  }
  return(smoothed)
}

#' Fit parametric marginal distribution to the data (internal)
#'
#'#Two options are possible:
#' -gaussian marginal for temperature
#' -gamma marginal for precipitation
#'
#' @title .fitMarginal
#' @param data Data used to fit the selected distribution
#' @param type Type of distribution to be fitted
#'
#' @return Parameters of the fitted distribution
#'
#' @keywords internal
#'
.fitMarginal <- function(data,type){
#  require("fitdistrplus")
  if(type=="gamma"){
    params <- fitdistrplus::fitdist(data, type, method="mle", 
                                    optim.method="Nelder-Mead")
  } else {
    params <- fitdistrplus::fitdist(data, type, method="mle", 
                                    optim.method="Nelder-Mead")
  }
  return(params)
}

#' Conditional adjustment of temperature with respect to precipitation (internal).
#'
#' @title .adjPT
#' @param obs Observed time series.
#' @param ctrl Control period time series.
#' @param scen Scenario period time series.
#' @param eps Small value used to handle zero and unity probabilities.
#'
#' @return Adjusted time series of temperature and precipitation.
#'
#' @keywords internal
#'
# .adjPT <- function(obs,ctrl,scen,eps=1e-5){
# 
#   adj <- list(T=array(NA,dim=c(length(scen$T))),
#               P=array(NA,dim=c(length(scen$P))))
# 
#   p2 <- p1 <- pgamma(scen$P[scen$wet],shape=ctrl$margP$estimate[1],rate=ctrl$margP$estimate[2])
#   p2[which(p2==1)] <- p2[which(p2==1)]-eps
#   adj$P[scen$wet] <- qgamma(p2,shape=obs$margP$estimate[1],rate=obs$margP$estimate[2])
#   adj$P[scen$dry] <- 0.0
# 
#   # Calculate conditional probability of temperature given precipitation
#   # for future simulations with parameters from historical simulation
# 
#   Fp <- array(NA,dim=c(length(scen$P)))
#   Ft <- array(NA,dim=c(length(scen$T)))
# 
#   p3 <- Ft
# 
#   Fp[scen$wet] <- pgamma(scen$P[scen$wet],shape=ctrl$margP$estimate[1],rate=ctrl$margP$estimate[2])
#   Fp[scen$dry] <- scen$P[scen$dry]
# 
#   Ft[scen$wet] <- pnorm(scen$stdTW,mean=ctrl$margTW$estimate[1],sd=ctrl$margTW$estimate[2])
#   Ft[scen$dry] <- pnorm(scen$stdTD,mean=ctrl$margTD$estimate[1],sd=ctrl$margTD$estimate[2])
# 
#   h <- pmax(pmin(pnorm((qnorm(Ft[scen$wet])-ctrl$cpar1*qnorm(Fp[scen$wet]))/sqrt(1-ctrl$cpar1^2)),1-eps),eps)
#   p3[scen$wet] <- unlist(h)
#   p3[scen$dry] <- Ft[scen$dry]
# 
#   #Finally, calculate the bias corrected temperature, conditioned on precipitation
#   #with the observed dependence structure
#   u1 <- pgamma(adj$P[scen$wet],shape=obs$margP$estimate[1],rate=obs$margP$estimate[2])
#   u1[which(u1==1)] <- 1-eps
#   u2 <- p3[scen$wet]
# 
#   tmp <- pnorm(qnorm(u2)*sqrt(1-obs$cpar1^2)+obs$cpar1*qnorm(u1))
#   #  tmp[which(tmp==1)] <- tmp[which(tmp==1)]-1e-5
#   adj$T[scen$wet] <- qnorm(tmp,mean=obs$margTW$estimate[1],sd=obs$margTW$estimate[2])
#   adj$T[scen$dry] <- qnorm(p3[scen$dry],mean=obs$margTD$estimate[1],sd=obs$margTD$estimate[2])
# 
#   return(adj)
# }

.adjPT <- function(obs,ctrl,scen,fit.skew,eps=1e-15){
  
  adj <- list(T = array(NA, dim = c(length(scen$T))),
              P = array(NA, dim = c(length(scen$P))))
  
  p2 <- pgamma(scen$P[scen$wet], shape=ctrl$margP$estimate[1], 
               rate=ctrl$margP$estimate[2])
  p2[which(p2==1)] <- p2[which(p2==1)]-eps
  
  adj$P[scen$wet] <- qgamma(p2, shape=obs$margP$estimate[1], 
                            rate=obs$margP$estimate[2])
  adj$P[scen$dry] <- 0.0

  # Calculate conditional probability of temperature given precipitation
  # for future simulations with parameters from historical simulation

  Fp <- array(NA,dim=c(length(scen$P)))
  Ft <- array(NA,dim=c(length(scen$T)))

  p3 <- Ft

  Fp[scen$wet] <- pgamma(scen$P[scen$wet], shape=ctrl$margP$estimate[1], 
                         rate=ctrl$margP$estimate[2])
  Fp[scen$dry] <- NA
  if(fit.skew){
    Ft[scen$wet] <- psn(scen$T[scen$wet], xi=ctrl$margTW$estimate[1], 
                        omega=ctrl$margTW$estimate[2], 
                        alpha=ctrl$margTW$estimate[3])
    Ft[scen$dry] <- psn(scen$T[scen$dry], xi=ctrl$margTD$estimate[1], 
                        omega=ctrl$margTD$estimate[2], 
                        alpha=ctrl$margTD$estimate[3])
  }else{
    Ft[scen$wet] <- pnorm(scen$T[scen$wet], mean=ctrl$margTW$estimate[1],
                          sd=ctrl$margTW$estimate[2])
    Ft[scen$dry] <- pnorm(scen$T[scen$dry], mean=ctrl$margTD$estimate[1],
                          sd=ctrl$margTD$estimate[2])
  }

  h <- pmax(pmin(pnorm((qnorm(Ft[scen$wet])-ctrl$cpar1*qnorm(Fp[scen$wet]))/sqrt(1-ctrl$cpar1^2)),(1-eps)),eps)
  plot(qnorm(Ft[scen$wet])-ctrl$cpar1*qnorm(Fp[scen$wet]))
#  p3[scen$wet] <- unlist(h)
#  p3[scen$dry] <- Ft[scen$dry]

  #Finally, calculate the bias corrected temperature, conditioned on precipitation
  #with the observed dependence structure
  u1 <- pgamma(adj$P[scen$wet], shape=obs$margP$estimate[1], 
               rate=obs$margP$estimate[2])
  u1[which(u1==1)] <- 1-eps
  u2 <- unlist(h) #p3[scen$wet]
  
  tmp <- pnorm(qnorm(u2)*sqrt(1-obs$cpar1^2)+obs$cpar1*qnorm(u1))
  plot(tmp)
  if(fit.skew){
    adj$T[scen$wet] <- qsn(tmp, xi=obs$margTW$estimate[1], omega=obs$margTW$estimate[2], 
                           alpha=obs$margTW$estimate[3])
    adj$T[scen$dry] <- qsn(Ft[scen$dry], xi=obs$margTD$estimate[1], omega=obs$margTD$estimate[2], 
                           alpha=obs$margTD$estimate[3])
  }else{
    adj$T[scen$wet] <- qnorm(tmp, mean=obs$margTW$estimate[1], sd=obs$margTW$estimate[2])
    adj$T[scen$dry] <- qnorm(Ft[scen$dry], mean=obs$margTD$estimate[1], sd=obs$margTD$estimate[2])
  }
  #  tmp[which(tmp==1)] <- tmp[which(tmp==1)]-1e-5

  return(adj)
}

#' Conditional adjustment of precipitation with respect to temperature (internal).
#'
#' @title .adjTP
#' @param obs A list containing the observed time series of temperature and precipitation.
#' @param ctrl A list containing the simulated time series of temperature and precipitation in the control period.
#' @param scen A list containing the simulated time series of temperature and precipitation in the scenario period.
#' @param eps Small value used to handle zero and unity probabilities.
#'
#' @return Adjusted time series of temperature and precipitation.
#'
#' @keywords internal
#'
# .adjTP <- function(obs,ctrl,scen,eps=1e-5){
# 
#   adj <- list(T=array(NA,dim=c(length(scen$T))),
#               P=array(NA,dim=c(length(scen$P))),
#               pwet=NA,
#               wet=NA,
#               dry=NA,
#               cpar1=NA)
# 
#   pw <- pmax(pmin(pnorm(scen$stdTW,mean=ctrl$margTW$estimate[1],sd=ctrl$margTW$estimate[2]),1-eps),eps)
#   pd <- pmax(pmin(pnorm(scen$stdTD,mean=ctrl$margTD$estimate[1],sd=ctrl$margTD$estimate[2]),1-eps),eps)
# 
#   adj$T[scen$wet] <- qnorm(pw,mean=obs$margTW$estimate[1],sd=obs$margTW$estimate[2])
#   adj$T[scen$dry] <- qnorm(pd,mean=obs$margTD$estimate[1],sd=obs$margTD$estimate[2])
# 
#   Fp <- pmin(pgamma(scen$P[scen$wet],shape=ctrl$margP$estimate[1],rate=ctrl$margP$estimate[2]),1-eps)
#   Ft <- pmax(pmin(pnorm(scen$T[scen$wet],mean=ctrl$margTW$estimate[1],sd=ctrl$margTW$estimate[2]),1-eps),eps)
#   p3 <- u2 <- h <- pnorm((qnorm(Fp)-ctrl$cpar1*qnorm(Ft))/sqrt(1-ctrl$cpar1^2))
# 
#   adj$wet <- scen$wet
#   adj$dry <- scen$dry
# 
#   #Finally, calculate the bias corrected precipitation, conditioned on temperature
#   #with the observed depence structure.
#   #  u2 <- p3[adj$wet]
#   u1 <- pnorm(adj$T[adj$wet],mean=obs$margTW$estimate[1],sd=obs$margTW$estimate[2])
# 
#   tmp <- pnorm(qnorm(u2)*sqrt(1-obs$cpar1^2)+obs$cpar1*qnorm(u1))
#   adj$P[adj$wet] <- qgamma(tmp,shape=obs$margP$estimate[1],rate=obs$margP$estimate[2])
#   adj$P[adj$dry] <- 0.0
# 
#   return(adj)
# }


.adjTP <- function(obs,ctrl,scen,fit.skew,eps=1e-15){

  adj <- list(T = array(NA, dim = c(length(scen$T))),
              P = array(NA, dim = c(length(scen$P))))
  
  if(fit.skew){
    pw <- pmax(pmin(psn(scen$T[scen$wet], xi = ctrl$margTW$estimate[1], 
                        omega = ctrl$margTW$estimate[2], 
                        alpha = ctrl$margTW$estimate[3]), (1-eps)), eps)
    pd <- pmax(pmin(psn(scen$T[scen$dry], xi = ctrl$margTD$estimate[1], 
                        omega = ctrl$margTD$estimate[2], 
                        alpha = ctrl$margTD$estimate[3]), (1-eps)), eps)
    
    adj$T[scen$wet] <- qsn(pw, xi=obs$margTW$estimate[1], 
                           omega=obs$margTW$estimate[2], 
                           alpha=obs$margTW$estimate[3])
    adj$T[scen$dry] <- qsn(pd, xi=obs$margTD$estimate[1], 
                           omega=obs$margTD$estimate[2], 
                           alpha=obs$margTD$estimate[3])
  }else{
    pw <- pnorm(scen$T[scen$wet], mean=ctrl$margTW$estimate[1], 
                          sd=ctrl$margTW$estimate[2])
    pd <- pnorm(scen$T[scen$dry], mean=ctrl$margTD$estimate[1], 
                          sd=ctrl$margTD$estimate[2])
    
    adj$T[scen$wet] <- qnorm(pw, mean=obs$margTW$estimate[1], 
                             sd=obs$margTW$estimate[2])
    adj$T[scen$dry] <- qnorm(pd, mean=obs$margTD$estimate[1], 
                             sd=obs$margTD$estimate[2])
  }
  
  Fp <- array(NA,dim=c(length(scen$P)))
  Ft <- array(NA,dim=c(length(scen$T)))
  
  Fp[scen$wet] <- pmin(pgamma(scen$P[scen$wet], shape=ctrl$margP$estimate[1], 
                    rate=ctrl$margP$estimate[2]), (1-eps))
  
#  p3 <- Fp
  
  if(fit.skew){
    Ft <- pmax(pmin(psn(scen$T[scen$wet], xi=ctrl$margTW$estimate[1], 
                        omega=ctrl$margTW$estimate[2], 
                        alpha=ctrl$margTW$estimate[3]), (1-eps)), eps)
  }else{
    Ft[scen$wet] <- pnorm(scen$T[scen$wet], mean=ctrl$margTW$estimate[1], 
                          sd=ctrl$margTW$estimate[2])
  }
  
  h <- pnorm((qnorm(Fp[scen$wet]) - ctrl$cpar1*qnorm(Ft[scen$wet]))/sqrt(1 - ctrl$cpar1^2))
  
#  p3[scen$wet] <- unlist(h)
#  p3[scen$dry] <- Fp[scen$dry]
  
  #Finally, calculate the bias corrected precipitation, conditioned on temperature
  #with the observed dependence structure.
  u2 <- unlist(h)#p3[adj$wet]
  if(fit.skew){
    u1 <- psn(adj$T[scen$wet], xi=obs$margTW$estimate[1], 
              omega=obs$margTW$estimate[2], 
              alpha=obs$margTW$estimate[3])
  }else{
    u1 <- pnorm(adj$T[scen$wet], mean=obs$margTW$estimate[1], 
                sd=obs$margTW$estimate[2])
  }
  tmp <- pnorm(qnorm(u2)*sqrt(1-obs$cpar1^2)+obs$cpar1*qnorm(u1))
  adj$P[scen$wet] <- qgamma(tmp, shape=obs$margP$estimate[1], 
                           rate=obs$margP$estimate[2])
  adj$P[scen$dry] <- 0.0
  
  return(adj)
}
