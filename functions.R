######################## Sub-routines ###############
#
#####################################################

### Calculate mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

### projection function
proj <- function(x){
  x <- as.matrix(x)
  return(x%*%solve(t(x)%*%x)%*%t(x))
}

rortho = function(n) {
  return(cbind(qr.Q(qr(array(runif(n), dim = c(n, n))))))
}

theta <- function(a,b) {
  acos( sum(a*b)/( sqrt(sum(a*a))*sqrt(sum(b*b)) ) )
}

## generate AR(1) correlation matrix
ar1cor <- function(p,rho){
  if(rho==0) sigma = diag(p)
  else sigma = rho^abs(outer(1:p,1:p,'-'))
  return(sigma)
}

## find the power of a matrix
matpower <- function(A,n){
  A = round((A+t(A))/2,7)
  tmp = eigen(A)
  return(tmp$vectors%*%diag((tmp$values)^n)%*%t(tmp$vectors))
}


## discretize the continuous response
discretize = function(y,h){
  n = length(y)
  m = floor(n/h)
  y = y + 0.00001*mean(y)*rnorm(n)
  yord = y[order(y)]
  divpt = numeric();
  for(i in 1:(h-1)) divpt = c(divpt,yord[i*m+1])
  y1 = rep(0,n)
  y1[y < divpt[1]] = 1
  y1[y >= divpt[h-1]] = h
  for(i in 2:(h-1)) y1[(y >= divpt[i-1]) & (y < divpt[i])] = i
  return(y1)
}


## Gaussian RBF kernel of Distance matrix
kern <- function(D){
  n <- dim(D)[1]
  Dupper <- D[upper.tri(D,diag = FALSE)]
  sigma2 <- sum(Dupper^2)/choose(n,2)
  gamma <- 1/(2*sigma2)
  return(exp(-gamma*D))
}


## kernel matrix for SIR
sirmat = function(Z,y,h){
  ydis =  discretize(y,h)
  ylab = unique(ydis)
  prob = table(ydis)/length(ydis)
  exy = numeric()
  for(i in 1:length(prob)) exy = rbind(exy,apply(Z[ydis==ylab[i],,drop=FALSE],2,mean))
  sirmat = t(exy)%*%diag(prob)%*%exy
  return(sirmat)
}


#### Computing pairwise distance matrices ############
#
#######################################################

#### Pairwise Euclidean Distances
edist <- function(Y){
  D <- as.matrix(dist(Y,method="euclidean",diag=TRUE,upper=TRUE))
  return(D)
}

#### Pairwise Wasserstein Distances
wdist <- function(Y){  ## Y is a list
  n <- length(Y)
  
  ##upper triangular index
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  
  dist <- function(i,j){
    return(wasserstein1d(Y[[i]],Y[[j]]))
  }
  
  kupper <- mapply(dist, i=ind[,1], j=ind[,2])
  k <- matrix(0,n,n)
  k[upper.tri(k,diag = FALSE)]=kupper^2
  k <- k+t(k)
  return(k)
}


### Pairwise Wasserstein distance on smoothed histograms
den_dist <- function(Y){
  n <- length(Y)
  
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  
  dist <- function(i,j){
    return(dist4den(Y[[i]], Y[[j]]))
  }
  
  kupper <- mapply(dist, i=ind[,1], j=ind[,2])
  k <- matrix(0,n,n)
  k[upper.tri(k,diag = FALSE)]=kupper^2
  k <- k+t(k)
  return(k)
}


## Pairwise Fourier distance
FourierDistance <- function(x, y, n=(floor(length(x)/2)+1)) {
  fft1 <- fft(x)
  fft2 <- fft(y)
  d <- sqrt(sum(Mod(fft1[1:n] - fft2[1:n]) ^ 2))
  
  return(d)
}


tsdist <- function(Y){
  n <- length(Y)
  
  ##upper triangular index
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  
  dist <- function(i,j){
    return(FourierDistance(Y[[i]], Y[[j]], n=20) )
  }
  
  kupper <- mapply(dist, i=ind[,1], j=ind[,2])
  k <- matrix(0,n,n)
  k[upper.tri(k,diag = FALSE)]=kupper^2
  k <- k+t(k)
  return(k)
}

###### Metric Cook's Distance  #####
mcd <- function(X,Dy){
  
  ## find MMDS scores
  fit <- cmdscale(Dy, eig=TRUE, k=2) 
  Sy <- fit$points[,1] 
  
  ## fit linear model Sy
  fit.lm <- lm(Sy ~ X)
  mcd <- cooks.distance(fit.lm)
  mcd_df <- data.frame(i=1:length(mcd),mcd)
  
  return(list(MCD = mcd, MCD.index=mcd_df, Sy = Sy))
}




######### Surrogate-assisted Methods ############
##
#################################################

## sa-OLS estimator
sa_ols <- function(X,D,N=1000){
  
  n=dim(X)[1]
  p = dim(X)[2]
  
  ## Standardize X
  mu = apply(X,2,mean)
  signrt = matpower(var(X),-1/2)
  Z = (X-mu)%*%signrt
  
  ## generate unit vectors
  set.seed(1)
  G = rmvnorm(N, mean=rep(0,n),sigma=diag(n))
  U = G / apply(G,1,function(x) Norm(x))
  
  ## surrogate responses
  Ys <- D%*%t(U)
  
  ## beta estimate for each response
  beta <- apply(Ys, MARGIN=2, cov, x=Z)
  M <- beta%*%t(beta)/N
  
  evals <- eigen(M)$values
  bbt <- eigen(M)$vectors
  beta.final <- signrt%*%bbt
  
  return(list(dir = beta.final, M=M, evalues=evals))
}


## sa-SIR estimator
sa_sir <- function(X,D,h=5,N=1000){
  
  n=dim(X)[1]
  p = dim(X)[2]
  
  ## Standardize X
  mu = apply(X,2,mean)
  signrt = matpower(var(X),-1/2)
  Z = (X-mu)%*%signrt
  
  ## generate unit vectors 
  set.seed(1)
  G = rmvnorm(N, mean=rep(0,n),sigma=diag(n))
  U = G / apply(G,1,function(x) Norm(x))
  
  ## surrogate response
  Ys <- D%*%t(U)
  
  mhat <- lapply(1:N, function(k) sirmat(Z,Ys[,k],h))
  Mhat <- do.call(cbind, mhat)
  M <- (Mhat%*%t(Mhat))/N
  
  evals <- eigen(M)$values
  Bhat <- signrt%*%eigen(M)$vectors
  return(list(dir=Bhat ,M=M, evalues=evals))
}



############### Frechet Methods  ####
#
#####################################

## Frechet-OLS estimator
fols <- function(X,D){
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  ## Standardize X
  mu = apply(X,2,mean)
  signrt <- matpower(var(X),-1/2)
  Z <- (X-mu)%*%signrt
  
  ## universal kernel for D
  y <- kern(D)
  
  beta <- apply(y, MARGIN=2, cov, x=Z)
  M <- beta%*%t(beta)/n
  
  evals <- eigen(M)$values
  bbt <- eigen(M)$vectors
  beta.final <- signrt%*%bbt
  return(list(dir=beta.final, M=M, evalues = evals))
}


## Frechet-SIR estimator

fsir <- function(X,D,h=5) {
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  ## Standardize X
  mu=apply(X,2,mean)
  signrt=matpower(var(X),-1/2)
  Z=(X-mu)%*%signrt
  
  ## Universal kernel for D
  y <- kern(D)
  
  ## find sir candidate matrices
  smeans <- numeric()
  for(j in 1:n)  smeans = cbind(smeans, sirmat(Z,y[j,],h) )
  
  Mhat <- smeans%*%t(smeans)/n
  evals <- eigen(Mhat)$values
  bbt <- eigen(Mhat)$vectors
  beta.final <- signrt%*%bbt
  
  return(list(dir = beta.final,M=Mhat, evalues = evals))
  
}



