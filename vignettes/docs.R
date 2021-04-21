## ----global_options, include=FALSE--------------------------------------------
# R output pre blocks are styled by default to indicate output
knitr::opts_chunk$set(comment=NA)

## -----------------------------------------------------------------------------
library(BCUH)

data("obs")
data("model")

## -----------------------------------------------------------------------------
bc6 <- biasco(obs=varsO[,'temp'], ctrl=varsC[,'temp'], scen=varsC[,'temp'], type="abs", method="M6")
bc9 <- biasco(obs=varsO[,'temp'], ctrl=varsC[,'temp'], scen=varsC[,'temp'], type="abs", method="M9")


## -----------------------------------------------------------------------------
plot(quantile(varsO[,'temp'], seq(0,1,0.01)), type="l",
    main="Quantile plot", xlab="%", ylab="Celcius", ylim=c(-30,10))
lines(quantile(dat(adj(bc6)), seq(0,1,0.01)), lty=1, col="orange")
lines(quantile(dat(adj(bc9)), seq(0,1,0.01)), lty=2, col="orange")
lines(quantile(varsC[,'temp'], seq(0,1,0.01)), col="blue")
legend("topleft", c("Obs","M6","M9","Ctrl"),
      col=c("black","orange","orange","blue"), lty=c(1,1,2,1))

## -----------------------------------------------------------------------------
bc6 <- biasco(obs=varsO[,'prec'], ctrl=varsC[,'prec'], scen=varsC[,'prec'], type="ratio", method="M6")
bc9 <- biasco(obs=varsO[,'prec'], ctrl=varsC[,'prec'], scen=varsC[,'prec'], type="ratio", method="M9")

## -----------------------------------------------------------------------------
plot(quantile(varsO[,'temp'], seq(0,1,0.01)), type="l",
    main="Quantile plot", xlab="%", ylab="mm/d", ylim=c(-30,10))
lines(quantile(dat(adj(bc6)), seq(0,1,0.01)), lty=1, col="orange")
lines(quantile(dat(adj(bc9)), seq(0,1,0.01)), lty=2, col="orange")
lines(quantile(varsC[,'temp'], seq(0,1,0.01)), col="blue")
legend("topleft", c("Obs","M6","M9","Ctrl"),
      col=c("black","orange","orange","blue"), lty=c(1,1,2,1))

## -----------------------------------------------------------------------------
#test <- biasco2D(obs.in=varsO, ctrl.in=varsC, scen.in=varsC, cond = "T")
#plot(dat(adj(test)))
#points(varsO,col="red")

