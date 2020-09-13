#' @title The function to fit time varying generalized autoregressive conditional heteroscedastic (p,q) model for count data
#' @description Takes the mean-zero time series as input and returns MCMC samples
#' @references Karmakar and Roy (2020)
#'     "Bayesian modelling of time-varying conditional heteroscedasticity"
#'
#' @param data is the time-series of count-valued data
#' @param order1 is the order of the time varying AR part, parameter p
#' @param order2 is the order of the time varying conditional variance, parameter q
#' @param knot is the number of equidistant knots for B-spline
#' @param norder is the order of B-splines
#' @param P is the order in tvAR() for initialization
#' @param Total_sample is the total number of iterations of MCMC
#' @param burn is the number of burn-in MCMC samples

#' @return fit.tvGARCHMCMCcomboN returns a list of the posterior samples of lambda0, mean, AR and CH coefficient
#' functions (sig2lp, Mfn, Afn and Bfn)\cr

#Assume order1>0 if order2>0
fit.tvGARCHMCMCcomboN <- function(data, order1 = 5, order2 = 0, knot = 6, norder = 4, P=10, Total_itr = 10000, burn = 5000){
  library(fda)
  library(pracma)
  set.seed(1)
  HMC = function (U, grad_U, epsilon, L = 30, current_q, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  HMCA = function (U, grad_U, epsilon, L = 30, current_q, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      q = q * (q >= 0)
      q = (q > 1) + q * (q <= 1)
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  HMCS = function (U, grad_U, epsilon, L = 30, current_q, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      q = q * (q > 0)
      q[q<=0] = 1e-10
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  order <- max(order1, order2)
  time <- (1:length(data)) / length(data)
  
  J       <- knot + norder - 1 
  timesp  <- bsplineS(time,  breaks=seq(0,1,1/knot), norder = norder)
  timespI <- timesp
  
  if(order>0){timespI <- timesp[-(1:order), ]}
  #timesp  <- matrix(rep(timespI, order), nrow = nrow(timespI))
  
  timespIder <- bsplineS(time,  breaks=seq(0,1,1/knot), norder = norder, nderiv = 1)
  timespIder <- timespIder[-(1:order), ]
  
  Umu <- function(x){
    mut     <- array(timespI %*% array(exp(x)))
    comp2 = array(At * X)
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * X))}
    
    vart <- array(mut) + array(comp2)
    
    if(order2>0){
      temp <- sigma2lat
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order+i-1):i])
        temp    <- c(temp, vart[i])
      } 
    }
    
    compo   <- Y/vart + log(vart)
    return(sum(compo)/2 + sum(x^2) / (2*100))
  }
  
  grad_Umu <- function(x){
    mut     <- array(timespI %*% array(exp(x)))
    comp2 = array(At * X)
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * X))}
    
    vart <- array(mut) + array(comp2)
    
    if(order2>0){
      temp <- sigma2lat
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order+i-1):i])
        temp    <- c(temp, vart[i])
      } 
    }
    
    compo   <- - Y/vart^2 + 1/vart
    
    compo   <- t(timespI) %*% array(compo) 
    return(array(compo)*exp(x)/2 + x / 100)
  }
  
  UA <- function(x){
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- timespI %*% x[(j-1)*J+1:J]*M[j+1]#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
    }
    comp2 = At * X
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * X))}
    
    vart <- array(mut) + array(comp2)
    
    if(order2>0){
      temp <- sigma2lat
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order+i-1):i])
        temp    <- c(temp, vart[i])
      } 
    }
    
    compo   <- Y/vart + log(vart)
    return(sum(compo)/2 + sum(x^2) / (2*100))
  }
  
  grad_UA <- function(x, deltaBc, Mc, temp){
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- timespI %*% x[(j-1)*J+1:J]*M[j+1]#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
    }
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- timespI %*% deltaBc[(j-1)*J+1:J]*Mc[order1+j+1]#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
    }
    comp2 = At * X
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * X))}
    
    vart <- array(mut) + array(comp2)
    
    if(order2>0){
      #temp <- sigma2lat
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order+i-1):i])
        temp    <- c(temp, vart[i])
      } 
    }
    
    if(order1==1){compo   <- array(t(timespI) %*% array((- Y/vart^2 + 1/vart)*array(X)))*M[order+1]}
    if(order1 > 1){
      compo <- NULL
      for(j in 1:order1){
        compo   <- c(compo, array(t(timespI) %*% array((- Y/vart^2 + 1/vart)*X[, j]))*M[j+1])
      }
    }
    #expxdev <- diag(exp(x)/sum(exp(x))) - tcrossprod(exp(x), exp(x)) /  (sum(exp(x)))^2
    #return(array(expxdev %*% compo) + (x) / (100))
    return(array(compo)/2)
  }
  
  UB <- function(x){
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- timespI %*% x[(j-1)*J+1:J]*M[order1+j+1]#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
    }
    comp2 = At * X
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * X))}
    
    vart <- array(mut) + array(comp2)
    
    if(order2>0){
      temp <- sigma2lat
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order+i-1):i])
        temp    <- c(temp, vart[i])
      } 
    }
    
    compo   <- Y/vart + log(vart)
    return(sum(compo)/2)# + sum(x^2) / (2*100))
  }
  
  grad_UB <- function(x, deltaAa, Mm, tempsig){
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- timespI %*% deltaAa[(j-1)*J+1:J]*Mm[j+1] #/  sum(exp(x))
    }
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- timespI %*% x[(j-1)*J+1:J]*Mm[order1+j+1]#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
    }
    comp2 = At * X
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * X))}
    
    vart <- array(mut) + array(comp2)
    
    sigma2m <- Bt
    
    if(order2>0){
      temp <- tempsig
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order2+i-1):i])
        sigma2m[i, ] <- temp[(order2+i-1):i]
        temp    <- c(temp, vart[i])
      } 
    }
    
    if(order2==1){compo   <- array(t(timespI) %*% array((- Y/vart^2 + 1/vart)*sigma2m))*M[order1+order2+1]}
    if(order2 > 1){
      compo <- NULL
      for(j in 1:order2){
        compo   <- c(compo, array(t(timespI) %*% array((- Y/vart^2 + 1/vart)*sigma2m[, j]))*M[order1+j+1])
      }
    }
    #expxdev <- diag(exp(x)/sum(exp(x))) - tcrossprod(exp(x), exp(x)) /  (sum(exp(x)))^2
    #return(array(expxdev %*% compo) + (x) / (100))
    return(array(compo)/2)
  }
  
  US <- function(x, deltaAa, deltaBb, Mm){
    ret = Inf
    if(!is.nan(x)){
    if(x>=0){
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- timespI %*% deltaAa[(j-1)*J+1:J]*Mm[j+1] #/  sum(exp(x))
    }
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- timespI %*% deltaBb[(j-1)*J+1:J]*Mm[order1+j+1]#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
    }
    vartcan <- length(Y)
    temp <- x
    if(order2>0){
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vartcan[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order2+i-1):i])
        temp    <- c(temp, vart[i])
        #sigma2m[i, ] <- temp[(order2+i-1):i]
      } 
    }
    compo   <- sum(Y/vartcan + log(vartcan)) + sum(log(x) + (data[1:order2])/x)
    ret <- sum(compo)/2 - sum(dgamma(x, 0.1, 0.1, log = T))
    }
    }
    return(ret)# + sum(x^2) / (2*100))
  }
  
  grad_US <- function(x, deltaAa, deltaBb, Mm){
    ret <- jacobian(US, x, deltaAa = deltaAa, deltaBb = deltaBb, Mm=Mm)
    return(array(ret))
  }
  
  UM <- function(x){
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- timespI %*% deltaA[(j-1)*J+1:J]*exp(x[j+1]) /  sum(exp(x))
    }
    
    comp2 = At * X
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * X))}
    
    vart <- array(mut) + array(comp2)
    
    if(order2>0){
      temp <- sigma2lat
      Bt <- matrix(0, length(data) - order, order2)
      for(j in 1:order2){
        Bt[, j] <- timespI %*% deltaB[(j-1)*J+1:J]*exp(x[order1+j+1]) /  sum(exp(x))
      }
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order+i-1):i])
        temp    <- c(temp, vart[i])
      } 
    }
    
    compo   <- Y/vart + log(vart)
    return(sum(compo)/2 + sum(x^2) / (2*100))
  }
  
  grad_UM <- function(x, deltaAa, deltaBb, tempsig){
    At <- matrix(0, length(data) - order, order1)
    Atder <- matrix(0, length(data) - order, order1)
    
    for(j in 1:order1){
      At[, j] <- timespI %*% deltaAa[(j-1)*J+1:J]*exp(x[j+1]) /  sum(exp(x))#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
      Atder[, j] <- timespI %*% deltaAa[(j-1)*J+1:J]
    }
    
    comp2 = array(At * X)
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = rowSums(At * X)}
    vart <- array(mut) + array(comp2)
    
    if(order2>0){
      sigma2m <- Bt
      temp <- tempsig#sigma2lat
      Bt <- matrix(0, length(data) - order, order2)
      Btder <- matrix(0, length(data) - order, order2)
      for(j in 1:order2){
        Bt[, j] <- timespI %*% deltaBb[(j-1)*J+1:J]*exp(x[order1+j+1]) /  sum(exp(x))#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
        Btder[, j] <- timespI %*% deltaBb[(j-1)*J+1:J]
      }
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order+i-1):i])
        temp    <- c(temp, vart[i])
        sigma2m[i, ] <- temp[(order+i-1):i]
      } 
    }
    part1 <- - Y/vart^2 + 1/vart
    
    comp2der <- array(crossprod(part1, Atder*X))
    if(order2>0){comp2der <- array(crossprod(part1, cbind(Atder*X, Btder*sigma2m)))}
    
    #comp2der <- c(0, comp2der)
    expxdev <- diag(exp(x)/sum(exp(x))) - tcrossprod(exp(x), exp(x)) /  (sum(exp(x)))^2
    expxdev <- matrix(expxdev[, -1], nrow = nrow(expxdev))
    return(array(expxdev %*% comp2der)/2 + (x) / (100))
  }
  
  Ucombo <- function(x){
    temp    <- x[order2]
    ret = Inf
    if(!is.nan(temp)){
    if(temp>=0){
    deltaAc <- x[1:length(deltaA) + order2]
    deltaBc <- x[1:length(deltaB) + order2 + length(deltaA)]
    deltaMc <- x[1:length(deltaM) + order2 + length(deltaB) + length(deltaA)]
    
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- timespI %*% deltaAc[(j-1)*J+1:J]*exp(deltaMc[j+1]) /  sum(exp(deltaMc))
    }
    
    comp2 = At * X
    if(order1==0){comp2 = rep(0, length(mut))}
    if(order1>1){comp2 = array(rowSums(At * X))}
    
    vart <- array(mut) + array(comp2)
    
    if(order2>0){
      #temp <- sigma2lat
      Bt <- matrix(0, length(data) - order, order2)
      for(j in 1:order2){
        Bt[, j] <- timespI %*% deltaBc[(j-1)*J+1:J]*exp(deltaMc[order1+j+1]) /  sum(exp(deltaMc))
      }
      vart <- rep(0, length(data)-order) 
      for(i in 1:(length(data)-order)){
        vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order+i-1):i])
        temp    <- c(temp, vart[i])
      } 
    }
    
    compo   <- Y/vart + log(vart)
    ret <- sum(compo)/2 + sum(deltaMc^2) / (2*100)
     } 
    }
    return(ret)# + sum(x^2) / (2*100))
  }
  
  grad_Ucombo <- function(x){
    temp    <- x[order2]
    deltaAc <- x[1:length(deltaA) + order2]
    deltaBc <- x[1:length(deltaB) + order2 + length(deltaA)]
    deltaMc <- x[1:length(deltaM) + order2 + length(deltaB) + length(deltaA)]
    
    Mc <- exp(deltaMc)/sum(exp(deltaMc))
    
    ret1 <- grad_US(temp, deltaAc, deltaBc, Mc)
    ret2 <- grad_UA(deltaAc, deltaBc, Mc, temp)
    ret3 <- grad_UB(deltaBc, deltaAc, Mc, temp)
    ret4 <- grad_UM(deltaMc, deltaAc, deltaBc, temp)
    ret <- c(ret1, ret2, ret3, ret4)
    return(array(ret))
  }
  
  HMC_combo = function (U, grad_U, epsilon, L = 30, current_q, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      
      q0 <- q[order2]
      q1 <- q[1:length(deltaA) + order2]
      q2 <- q[1:length(deltaB) + order2 + length(deltaA)]
      q3 <- q[1:length(deltaM) + order2 + length(deltaB) + length(deltaA)]
      
      #q = q * (q > 0)
      q0 = q0 * (q0 > 0) + (q0<=0)*1e-10
      q1 = q1 * (q1 >= 0)
      q1 = (q1 > 1) + q1 * (q1 <= 1)
      
      q2 = q2 * (q2 >= 0)
      q2 = (q2 > 1) + q2 * (q2 <= 1)
      # Make a full step for the momentum, except at end of trajectory
      
      q <- c(q0, q1, q2, q3)
      if (j!=L) p = p - epsilon * grad_U(q)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if(is.na(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  Y <- data^2
  if(order>0){Y <- (data[-(1:order)])^2}
  X <- NULL
  if(order1 > 0){
    for(j in 1:order1){
      ind <- (order-j+1):(length(data) - j)
      X <- cbind(X, (array(data[ind]))^2)
    } 
  }
  At <- 0
  #est=summary(garch(data, order = c(1,1), coef = NULL, itmax = 1000,eps = NULL, grad = c("analytical","numerical"),series = NULL, trace = F))$coef
  
  #P = 10
  fit <- tvAR(data^2, p = P) 
  #design <- ginv(t(timesp[-1:15]) %*% timesp[-1:15]) %*% t( timesp[-1:15])
  
  p1 <- fit$coefficients[, 1]
  p2 <- fit$coefficients[, 2]
  p3 <- rowSums(fit$coefficients[, 3:P])
  
  p23 <- p2 / p3 + 1
  
  bcoef  <- 1/p23
  mucoef <- p1*(1-bcoef)
  acoef  <- p2
  
  fit1    <- lm(mucoef~timesp[-c(1:P),]-1)
  temp    <- fit1$coefficients
  temp[which(is.na(temp))] <- 0.1
  temp[(temp < 0)] <- 0.1
  deltamu <- log(array(temp))
  
  fit2    <- lm(acoef~timesp[-c(1:P),]-1)
  temp    <- fit2$coefficients
  temp[which(is.na(temp))] <- 0
  temp[(temp < 0)] <- 0
  deltaA  <- array(temp)
  
  fit3    <- lm(bcoef~timesp[-c(1:P),]-1)
  temp    <- fit3$coefficients
  temp[which(is.na(temp))] <- 0
  temp[(temp < 0)] <- 0
  deltaB  <- array(temp)
  
  if((order1+order2)>0){
    deltaM  <- rep(0, order1+order2+1)#rnorm(order1+order2+1)
    deltaM[1] <- -30
    M       <- exp(deltaM) / sum(exp(deltaM))
  }
  mut     <- timespI %*% exp(deltamu)
  
  deltaA <- deltaA / M[2]
  deltaB <- deltaB / M[3]
  
  deltaA[(deltaA>1)] <- 1
  deltaB[(deltaB>1)] <- 1
  
  if(order1>0){
    At <- matrix(0, length(data) - order, order1)
    for(j in 1:order1){
      At[, j] <- timespI %*% deltaA[(j-1)*J+1:J]*M[j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
    }
  }
  
  if(order2>0){
    Bt <- matrix(0, length(data) - order, order2)
    for(j in 1:order2){
      Bt[, j] <- timespI %*% deltaB[(j-1)*J+1:J]*M[order1+j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
    }
  }
  sigma2lat <- NULL
  if(order2>0)
  {
    sigma2lat <- 1#(deltaA[1]*M[2])/(1-(deltaB[1]*M[3]))
  }
  delta <- c(sigma2lat, deltaA, deltaB, deltaM)
  temp <- sigma2lat
  if(order2>0){
    vart <- rep(0, length(data)-order) 
    for(i in 1:(length(data)-order)){
      vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order+i-1):i])
      temp    <- c(temp, vart[i])
    } 
  }
  
  itr <- 0
  arA <- 0
  arB <- 0
  arM <- 0
  armu <- 0
  arS  <- 0
  arcom <- 0
  sdA <- 1e-6
  sdB <- 1e-10
  sdM <- 1e-6
  sdmu <- 1e-4
  sdS  <- 1e-4
  sdcom <- 1e-6
  Als <- list()
  Bls <- list()
  Mls <- list()
  Mlsder <- list()
  Alsder <- list()
  Blsder <- list()
  siglatls <- list()
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  Pred <- rep(0, Total_itr)
  while(itr < Total_itr){
    itr <- itr + 1
    temp    <- HMC(Umu, grad_Umu, sdmu, L = 30, deltamu, armu)
    deltamu <- temp$up
    armu     <- temp$arc
    
    mut <- timespI %*% exp(deltamu)
    
    
    if(order>0){
      
      # if(order1>0){
      #   temp   <- HMCA(UA, grad_UA, sdA, L = 30, deltaA, arA)
      #   #print(sum(temp$up-deltaA)^2)
      #   deltaA <- temp$up
      #   arA    <- temp$arc
      #   
      #   At <- matrix(0, length(data) - order, order1)
      #   for(j in 1:order1){
      #     At[, j] <- timespI %*% deltaA[(j-1)*J+1:J]*M[j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
      #   }
      # }
      
      if(order2==0){
        temp   <- HMC(UM, grad_UM, sdM, L = 30, deltaM, arM)
        #print(sum(temp$up-deltaM)^2)
        deltaM <- temp$up
        arM    <- temp$arc
        
        M       <- exp(deltaM) / sum(exp(deltaM))
        
        }
      if(order2>0){
        temp   <- HMC_combo(Ucombo, grad_Ucombo, sdcom, L = 30, delta, arcom)
        #print(sum(temp$up-deltaA)^2)
        delta <- temp$up
        arcom <- temp$arc
        sigma2lat <- delta[order2]
        deltaA <- delta[1:length(deltaA) + order2]
        deltaB <- delta[1:length(deltaB) + order2 + length(deltaA)]
        deltaM <- delta[1:length(deltaM) + order2 + length(deltaB) + length(deltaA)] 
        
        M       <- exp(deltaM) / sum(exp(deltaM))
        
        At <- matrix(0, length(data) - order, order1)
        for(j in 1:order1){
          At[, j] <- timespI %*% deltaA[(j-1)*J+1:J]*M[j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
        }
        
        Bt <- matrix(0, length(data) - order, order2)
        for(j in 1:order2){
          Bt[, j] <- timespI %*% deltaB[(j-1)*J+1:J]*M[order1+j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
        }
        comp2 = At * X
        if(order1==0){comp2 = rep(0, length(mut))}
        if(order1>1){comp2 = array(rowSums(At * X))}
        
        vart <- array(mut) + array(comp2)
        
        if(order2>0){
          temp <- sigma2lat
          vart <- rep(0, length(data)-order) 
          for(i in 1:(length(data)-order)){
            vart[i] <- mut[i] + sum(At[i, ]*X[i, ])+ sum(Bt[i, ]*temp[(order2+i-1):i])
            temp    <- c(temp, vart[i])
          } 
        }
      }
    }
    
    print(Pred[itr] <- mean((data^2-temp)^2))
    siglatls[[itr]] <- sigma2lat
    
    if(itr %% 100 == 0){
      if(order > 0){
        if(order2==0){
          
          ar <- arM/ itr
          cat(ar, "acceptance rate for M")
          if(ar<.60){sdM <- sdM * (.1)}
          if(ar>.90){sdM <- sdM * (10)} 
        }
        
        ar <- arA/ itr
        cat(ar, "acceptance rate for A")
        if(ar<.60){sdA <- sdA * (.1)}
        if(ar>.90){sdA <- sdA * (10)}
        
        if(order2>0){
          ar <- arcom/ itr
          cat(ar, "acceptance rate for combo")
          if(ar<.60){sdcom <- sdcom * (.1)}
          if(ar>.90){sdcom <- sdcom * (10)}
          
          # ar <- arS/ itr
          # cat(ar, "acceptance rate for Latentsig")
          # if(ar<.60){sdS <- sdS * (.1)}
          # if(ar>.90){sdS <- sdS * (10)}  
        }
      }
      ar <- armu/ itr
      cat(ar, "acceptance rate for mu")
      if(ar<.60){sdmu <- sdmu * (.1)}
      if(ar>.90){sdmu <- sdmu * (10)}
    }
    
    mutder <- timespIder %*% exp(deltamu)
    if(order1>0){
      Atder <- matrix(0, length(data) - order, order1)
      for(j in 1:order1){
        Atder[, j] <- timespIder %*% deltaA[(j-1)*J+1:J]*M[j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
      }
      Als[[itr]] <- At 
      Alsder[[itr]] <- Atder
    }
    
    if(order2>0){
      Btder <- matrix(0, length(data) - order, order2)
      for(j in 1:order2){
        Btder[, j] <- timespIder %*% deltaB[(j-1)*J+1:J]*M[order1+j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
      }
      Bls[[itr]] <- Bt 
      Blsder[[itr]] <- Btder
    }
    
    Mls[[itr]] <- mut
    Mlsder[[itr]] <- mutder
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
    
    # plot(Atder)
    #plot(Bt[, order2])
    #points(At[, 1], col = 2)
    #points(At[, 2], col=3)
    #print(sigma2lat)
  }
  close(pb)
  if(order1>0){out <- list(sig0 = siglatls[(burn+1):Total_itr], pred = Pred[(burn+1):Total_itr], Afn = Als[(burn+1):Total_itr], Mfn = Mls[(burn+1):Total_itr], Afnder = Alsder[(burn+1):Total_itr], Mfnder = Mlsder[(burn+1):Total_itr])}
  if(order2>0){out <- list(sig0 = siglatls[(burn+1):Total_itr], pred = Pred[(burn+1):Total_itr], Afn = Als[(burn+1):Total_itr], Bfn = Bls[(burn+1):Total_itr], Mfn = Mls[(burn+1):Total_itr], Afnder = Alsder[(burn+1):Total_itr], Bfnder = Blsder[(burn+1):Total_itr], Mfnder = Mlsder[(burn+1):Total_itr])}
  if(order==0){out <- list(Mfn = Mls[(burn+1):Total_itr], Mfnder = Mlsder[(burn+1):Total_itr])}
  return(out)
}
