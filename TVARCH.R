#' @title The function to fit time varying autoregressive conditional heteroscedastic (p,q) model for count data
#' @description Takes the mean-zero time series as input and returns MCMC samples
#' @references Karmakar and Roy (2020)
#'     "Bayesian modelling of time-varying conditional heteroscedasticity"
#'
#' @param data is the time-series of count-valued data
#' @param order is the order of the time varying AR part, parameter p
#' @param knot is the number of equidistant knots for B-spline
#' @param norder is the order of B-splines
#' @param Total_sample is the total number of iterations of MCMC
#' @param burn is the number of burn-in MCMC samples

#' @return fit.tvARMCMCunbd returns a list of the posterior samples of mean and AR
#' functions (Mfn and Afn)\cr

fit.tvARMCMCunbd <- function(data, order = 5, knot = 4, norder = 4, Total_itr = 20000, burn = 10000){
  library(fda)
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
  
  HMCM = function (U, grad_U, epsilon, L = 30, current_q, arc)
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
      #q = q * (q >= 0)
      
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
  
  time <- (1:length(data)) / length(data)
  
  J       <- knot + 3 
  timesp  <- bsplineS(time,  breaks=seq(0,1,1/knot))
  timespI <- timesp
  
  if(order>0){timespI <- timesp[-(1:order), ]}
  #timesp  <- matrix(rep(timespI, order), nrow = nrow(timespI))
  
  timespIder <- bsplineS(time,  breaks=seq(0,1,1/knot), nderiv = 1)
  timespIder <- timespIder[-(1:order), ]
  
  Umu <- function(x){
    mut     <- array(timespI %*% array(exp(x)))
    comp2 = array(At * X)
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * X))}
    lambdat <- array(mut) + array(comp2)
    compo   <- log(lambdat) + Y / lambdat
    return(sum(compo)/2 + sum(x^2) / (2*100))
  }
  
  grad_Umu <- function(x){
    mut     <- array(timespI %*% array(exp(x)))
    comp2 = array(At * X)
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * X))}
    lambdat <- array(mut) + array(comp2)
    compo   <- t(timespI) %*% array(1/lambdat-Y/lambdat^2) 
    return(array(compo)*exp(x) + x / 100)
  }
  
  UA <- function(x){
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- timespI %*% x[(j-1)*J+1:J]*M[j+1]#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
    }
    comp2 = At * X
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = rowSums(At * X)}
    lambdat <- array(mut) + array(comp2)
    compo   <- log(lambdat) + Y / lambdat
    return(sum(compo)/2)# + sum(x^2) / (2*100))
  }
  
  grad_UA <- function(x){
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- timespI %*% x[(j-1)*J+1:J]*M[j+1]#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
    }
    comp2 = At * X
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = rowSums(At * X)}
    lambdat <- array(mut) + array(comp2)
    if(order==1){compo   <- array(t(timespI) %*% array(array(X)*(1/lambdat-Y/lambdat^2)))*M[order+1]}
    if(order > 1){
      compo <- NULL
      for(j in 1:order){
        compo   <- c(compo, array(t(timespI) %*% array(X[, j]*(1/lambdat-Y/lambdat^2)))*M[j+1])
      }
    }
    #expxdev <- diag(exp(x)/sum(exp(x))) - tcrossprod(exp(x), exp(x)) /  (sum(exp(x)))^2
    #return(array(expxdev %*% compo) + (x) / (100))
    return(array(compo))
  }
  
  UM <- function(x){
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- timespI %*% deltaA[(j-1)*J+1:J]*exp(x[j+1]) /  sum(exp(x))
    }
    comp2 = array(At * X)
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = rowSums(At * X)}
    lambdat <- array(mut) + array(comp2)
    compo   <- log(lambdat) + Y / lambdat
    return(sum(compo)/2 + sum(x^2) / (2*100))
  }
  
  grad_UM <- function(x){
    At <- matrix(0, length(data) - order, order)
    Atder <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      At[, j] <- timespI %*% deltaA[(j-1)*J+1:J]*exp(x[j+1]) /  sum(exp(x))#exp(x[(j-1)*J+1:J]) /  sum(exp(x))
      Atder[, j] <- timespI %*% deltaA[(j-1)*J+1:J]
    }
    comp2 = array(At * X)
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = rowSums(At * X)}
    lambdat <- array(mut) + array(comp2)
    
    part1 <- (1/lambdat-Y/lambdat^2)
    
    comp2der <- array(crossprod(part1, Atder*X))
    
    #comp2der <- c(0, comp2der)
    expxdev <- diag(exp(x)/sum(exp(x))) - tcrossprod(exp(x), exp(x)) /  (sum(exp(x)))^2
    expxdev <- matrix(expxdev[, -1], nrow = nrow(expxdev))
    return(array(expxdev %*% comp2der) + (x) / (100))
  }
  Y <- data^2
  if(order>0){Y <- (data[-(1:order)])^2}
  X <- NULL
  if(order > 0){
    for(j in 1:order){
      ind <- (order-j+1):(length(data) - j)
      X <- cbind(X, (array(data[ind]))^2)
    } 
  }
  At <- 0
  
  deltaA <- rep(0, ncol(timespI))
  fit <- tvAR(data^2, p = order)
  #design <- ginv(t(timesp[-1:15]) %*% timesp[-1:15]) %*% t( timesp[-1:15]
  P <- order
  fit1    <- lm(fit$coefficients[, 1]~timesp[-c(1:P),]-1)
  temp    <- fit1$coefficients
  temp[which(is.na(temp))] <- 0.1
  temp[(temp < 0)] <- 0.1
  deltamu <- log(array(temp))

  for(j in 1:order){
    fit2    <- lm(fit$coefficients[, 1+j]~timesp[-c(1:P),]-1)
    temp    <- fit2$coefficients
    temp[which(is.na(temp))] <- 0
    temp[(temp < 0)] <- 0
    deltaA[(j-1)*J+1:J]  <- array(temp)
  }
  
  #deltamu <- rnorm(ncol(timesp))
  #deltaA  <- runif(ncol(timespI)*order)
  
  if(order>0){
    deltaM  <- rnorm(order+1)
    deltaM[1] <- -30
    M       <- exp(deltaM) / sum(exp(deltaM))
    mut     <- timespI %*% exp(deltamu)
    At <- matrix(0, length(data) - order, order)
    for(j in 1:order){
      temp <- deltaA[(j-1)*J+1:J] / M[j+1]
      temp[temp>1] <- 1
      deltaA[(j-1)*J+1:J] <- temp
      At[, j] <- timespI %*% deltaA[(j-1)*J+1:J]*M[j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
    }
  }
  
  itr <- 0
  arA <- 0
  arM <- 0
  armu <- 0
  sdA <- 1e-6
  sdM <- 1e-6
  sdmu <- 1e-4
  Als <- list()
  Mls <- list()
  Mlsder <- list()
  Alsder <- list()
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  Pred <- rep(0, Total_itr)
  while(itr < Total_itr){
    itr <- itr + 1
    temp    <- HMC(Umu, grad_Umu, sdmu, L = 30, deltamu, armu)
    deltamu <- temp$up
    armu     <- temp$arc
    
    mut <- timespI %*% exp(deltamu)
    
    if(order>0){
      temp   <- HMCA(UA, grad_UA, sdA, L = 30, deltaA, arA)
      #print(sum(temp$up-deltaA)^2)
      deltaA <- temp$up
      arA    <- temp$arc
      
      At <- matrix(0, length(data) - order, order)
      for(j in 1:order){
        At[, j] <- timespI %*% deltaA[(j-1)*J+1:J]*M[j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
      }
      
      temp   <- HMCM(UM, grad_UM, sdM, L = 30, deltaM, arM)
      #print(sum(temp$up-deltaM)^2)
      deltaM <- temp$up
      arM    <- temp$arc
      
      M       <- exp(deltaM) / sum(exp(deltaM))
      At <- matrix(0, length(data) - order, order)
      for(j in 1:order){
        At[, j] <- timespI %*% deltaA[(j-1)*J+1:J]*M[j+1]# exp(deltaA[(j-1)*J+1:J]) /  sum(exp(deltaA))
      }
    }
    
    comp2 = array(At * X)
    if(order==0){comp2 = rep(0, length(mut))}
    if(order>1){comp2 = array(rowSums(At * X))}
    lambdat <- array(mut) + array(comp2)
    
    #print(Pred[itr] <- mean((Y-lambdat)^2))
    
    if(itr %% 100 == 0){
      if(order > 0){
        ar <- arA/ itr
        cat(ar, "acceptance rate for A")
        if(ar<.60){sdA <- sdA * (.1)}
        if(ar>.90){sdA <- sdA * (10)}
        
        ar <- arM/ itr
        cat(ar, "acceptance rate for M")
        if(ar<.60){sdM <- sdM * (.1)}
        if(ar>.90){sdM <- sdM * (10)} 
      }
      
      ar <- armu/ itr
      cat(ar, "acceptance rate for mu")
      if(ar<.60){sdmu <- sdmu * (.1)}
      if(ar>.90){sdmu <- sdmu * (10)}
    }
    Atder <- timespIder %*% deltaA[1:J]*M[2]
    mutder <- timespIder %*% exp(deltamu)
    if(order>0){
      Als[[itr]] <- At 
      Alsder[[itr]] <- Atder
    }
    Mls[[itr]] <- mut
    Mlsder[[itr]] <- mutder
    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, itr)
    
    #plot(At[, order])
    # plot(Atder)
    #plot(mut)
    #points(At[, 1], col = 2)
    #points(At[, 2], col=3)
  }
  close(pb)
  out <- list(pred = Pred[(burn+1):Total_itr], Afn = Als[(burn+1):Total_itr], Mfn = Mls[(burn+1):Total_itr], Afnder = Alsder[(burn+1):Total_itr], mufnder = Mlsder[(burn+1):Total_itr])
  if(order==0){out <- list(pred = Pred[(burn+1):Total_itr], Mfn = Mls[(burn+1):Total_itr], mufnder = Mlsder[(burn+1):Total_itr])}
  return(out)
}