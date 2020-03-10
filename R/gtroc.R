#' Nonparametric Estimation of Distributions and Diagnostic Accuracy Based on Group-Tested Results with Differential Misclassification
#' @import isotone pracma
#' @description This package concerns the problem of estimating a continuous distribution in a diseased or
#' nondiseased population when only group-based test results on the disease status are available. The problem
#'  is challenging in that individual disease statuses are not observed and testing results are often subject
#'  to misclassification, with further complication that the misclassification may be differential as the group
#'  size and the number of the diseased individuals in the group vary. We propose a method to construct
#'  nonparametric estimation of the distribution and obtain its asymptotic properties. The performance of the
#'  distribution estimator is evaluated under various design considerations concerning group sizes and
#'  classiffication errors. The method is exemplified with data from the National Health and Nutrition
#'  Examination Survey study to estimate the distribution and diagnostic accuracy of C-reactive protein in
#'  blood samples in predicting chlamydia incidence.
#' @return None
#'
#' @examples
#' plot_crayons()
#'
#' @export
#' #' @details
#' This is a ...
#' @author Karl W Broman, \email{broman@@wisc.edu}
#' @references \url{https://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors}
#' @seealso \code{\link{brocolors}}
#' @keywords hplot
#' ...
#' @importFrom

###############################################################
#### source the main functions
#source("./MainFunctions/sen_fun.R")
#source("./MainFunctions/loglik_prevalence.R")
#source("./MainFunctions/cdf_estimate.R")
#source("./MainFunctions/ROC_estimate.R")

###############################################################
gtroc<-function( pi0 = 0.9,
                 pi1 = 0.9,
                 lambda = 0.02,
                 tempDat,
                 nodeNum = 100,
                 ...
){
  K <- ncol(tempDat)-1
  n <- nrow(tempDat)
  groupM <- tempDat[,1]
  groupPos <- which(groupM==1)
  groupNeg <- which(groupM==0)
  n0 <- sum(groupM==0)
  n1 <- sum(groupM==1)
  Xmat <- tempDat[,-1]
  Xgroup <- split(Xmat, 1:n)

  ############################################################################################
  ######the ROC curve estimate  based on the estimated isotonic CDF
  #' @description the ROC curve estimate  based on the estimated isotonic CDF
  ROCEst <- function(u, F0.val, F1.val)
  {
    invLoc <- min(which(F0.val>=1-u))
    roc <- 1-F1.val[invLoc]
    return(roc)
  }
  ROCEst <- Vectorize(ROCEst, "u")


  ###################### Estimate the prevalence and the distribution functions
  temp.estRes <- CDFEst(Xgroup=Xgroup, Mvec=groupM, pi0=pi0, pi1=pi1, K=K, lambda=lambda)
  p.est <- temp.estRes$prevalEst
## use the derived equation to calculate the variance
  pi1.vec <- apply(matrix(1:K, ncol=1), 1, SenFun, pi1=pi1, K=K, lambda=lambda)
  temp.exp1 <- rep(NA, K)
  for(dd in 1:K)
  {
    temp.exp2 <- dd*p.est^(dd-1)*(1-p.est)^(K-dd)-(K-dd)*p.est^dd*(1-p.est)^(K-1-dd)
    temp.exp1[dd] <- pi1.vec[dd]*choose(K,dd)*temp.exp2
  }
  logLik_der1 <- K*(1-pi0)*(1-p.est)^(K-1)- sum(temp.exp1)
  p.sd <- sqrt(p.est*(1-p.est)/(n*logLik_der1^2))
## estimate the distribution function
  F0.fun <- temp.estRes$F0.est
  F1.fun <- temp.estRes$F1.est
  F0.val <- temp.estRes$F0.pava
  F1.val <- temp.estRes$F1.pava

####### the estimated CDF values at some given points (quantile)
  quan <- quantile(as.vector(as.matrix(Xmat)), seq(from=0.1, to=0.9, by=0.1))
  F0.AllEst <- apply(matrix(quan, ncol=1), 1, F0.fun)
  F1.AllEst <- apply(matrix(quan, ncol=1), 1, F1.fun)

#################### Estimate the ROC curve
  aa <- seq(from=1/nodeNum, to=1-1/nodeNum, by=1/nodeNum)
  ROC.value  <- apply(matrix(aa, ncol=1), MARGIN=1, ROCEst, F0.val=F0.val, F1.val=F1.val)

################## Estimate the AUC
  AUC.value <- integrate(ROCEst, 0, 1, subdivisions=5000, F0.val=F0.val, F1.val=F1.val)$value
  ROCAOC<- list(ROC.value, AUC.value)
  return(ROCAOC)
}


###########################################################################################
########## Calculate the sensitivity with the dilution effect
#' Calculate the sensitivity with the dilution effect
#' @description Calculate the sensitivity with the dilution effect
SenFun <- function(pi1, K, dnum, lambda)
{
  res <- pi1*dnum/(dnum+lambda*(K-dnum))
  return(res)
}

###################################################################################
###### estimate the distribution function in the diseased (case) and non-disease (control) population
#' estimate the distribution function in the diseased (case) and non-disease (control) population
#' @description estimate the distribution function in the diseased (case) and non-disease (control) population
CDFEst <- function(Xgroup, Mvec, pi0, pi1, K, lambda)
{
  n <- length(Xgroup)
  n0 <- sum(Mvec==0)
  n1 <- n-n0
  Xvec <- unlist(Xgroup)
  Xvec.neg <- unlist(Xgroup[Mvec==0])

  ### Step 1: estimate the prevalence based on the screening test observations
  #p.est <- optimize(f=Loglik_GT, interval=c(1e-6,1-1e-6), maximum=T, n0=n0, n1=n1, K=K, pi0=pi0, pi1=pi1, lambda=lambda, tol=.Machine$double.eps)$maximum
  p.est <- uniroot(f=Loglik_der1, interval=c(1e-6,1-1e-6), n0=n0, n1=n1, K=K, pi0=pi0, pi1=pi1, lambda=lambda, tol=.Machine$double.eps)$root
  ### Step 2: estimate the mixture CDF and the conditional CDF given D=0
  Fmix.est <- function(z)  return(mean(Xvec<=z))
  Fneg.est <- function(z)  return(sum(Xvec.neg<=z)/(n*K))

  ###
  d1.vec <- 1:K
  temp1_pi1 <- apply(matrix(1:K, ncol=1), 1, SenFun, pi1=pi1, K=K, lambda=lambda)
  tempA <- sum((1-temp1_pi1)*choose(K-1,d1.vec-1)*(p.est^d1.vec)*((1-p.est)^(K-d1.vec)))
  if(K>1)
  {
    temp2_pi1 <- apply(matrix(1:(K-1), ncol=1), 1, SenFun, pi1=pi1, K=K, lambda=lambda)
    d2.vec <- 1:(K-1)
    tempB <- pi0*(1-p.est)^K+sum((1-temp2_pi1)*choose(K-1,d2.vec)*(p.est^d2.vec)*((1-p.est)^(K-d2.vec)) )
  }
  else
  {
    tempB <- pi0*(1-p.est)^K
  }
  tempC <- p.est*tempB-tempA*(1-p.est)

  ### Step 3: obtain the initial nonparametric estimators of F0 and F1
  tempF0.est <- function(z)
  {
    temp3 <- p.est*Fneg.est(z)-tempA*Fmix.est(z)
    res.est <- temp3/tempC
    return(res.est)
  }
  tempF1.est <- function(z)
  {
    temp5 <- p.est*tempB*Fmix.est(z)-p.est*(1-p.est)*Fneg.est(z)
    temp6 <- p.est*tempC
    res.est <- temp5/temp6
    return(res.est)
  }

  ### Step 4: use the PAV algorithm to adjust the initial CDF estimators
  X.order <- sort(Xvec, decreasing=F)
  F0.vec <- apply(matrix(X.order, ncol=1), MARGIN=1, tempF0.est)
  F1.vec <- apply(matrix(X.order, ncol=1), MARGIN=1, tempF1.est)
  F0.pava <- gpava(z=X.order, y=F0.vec, ties="primary")$x
  F0.pava <- apply(matrix(F0.pava,ncol=1), MARGIN=1, function(z) min(max(z,0),1))
  F1.pava <- gpava(z=X.order, y=F1.vec, ties="primary")$x
  F1.pava <- apply(matrix(F1.pava,ncol=1), MARGIN=1, function(z) min(max(z,0),1))
  ### generate the step function based on the pava for F0
  F0.est <- stepfun(x=X.order, y=c(0,F0.pava), right=FALSE)
  F1.est <- stepfun(x=X.order, y=c(0,F1.pava), right=FALSE)

  return(list(prevalEst=p.est, F0.est=F0.est, F1.est=F1.est, F0.pava=F0.pava, F1.pava=F1.pava))
}

###########################################################################################
### the loglikelihood function for the prevalence
#' the loglikelihood function for the prevalence
#' @description the loglikelihood function for the prevalence
Loglik_GT <- function(n0, n1, p, K, pi0, pi1, lambda)
{
  d.vec <- 1:K
  pi1.vec <- apply(matrix(d.vec, ncol=1), 1, SenFun, pi1=pi1, K=K, lambda=lambda)
  temp1 <- pi1.vec*choose(K,d.vec)*(p^d.vec)*((1-p)^(K-d.vec))
  prob0 <- 1-(1-pi0)*(1-p)^K-sum(temp1)
  prob1 <- 1-prob0
  res <- n0*log(prob0)+n1*log(prob1)
  return(res)
}

### the first derivative of the loglikelihood function for the prevalence
#' the first derivative of the loglikelihood function for the prevalence
#' @description the first derivative of the loglikelihood function for the prevalence
Loglik_der1 <- function(n0, n1, p, K, pi0, pi1, lambda)
{
  d.vec <- 1:K
  pi1.vec <- apply(matrix(d.vec, ncol=1), 1, SenFun, pi1=pi1, K=K, lambda=lambda)
  temp1 <- pi1.vec*choose(K,d.vec)*(p^d.vec)*((1-p)^(K-d.vec))
  prob0 <- 1-(1-pi0)*(1-p)^K-sum(temp1)
  res <- n0/prob0-n1/(1-prob0)
  return(res)
}






