library(mvtnorm)
library(fdadensity)
library(transport)
library(pdist)
library(pracma)
library(dr)
library(dplyr)
library(frechet)
require(tictoc)
require(MASS)
require(magrittr)
require(dplyr)
library(parallel)
library(FCPS)
library(boot)
library(cdcsis)
library(ggplot2)
library(tidyr)
library(stringr)
library(statip)
library(plotly)
library(latex2exp)
library(gridExtra)
library(cowplot)
library(tictoc)
library(astsa)
library(psd)
library(nlme)
library(fda)
library(fda.usc)
library(celltrackR)
library(rstudioapi)
library(rrpack)
library(EnvStats)

## Set working directory and call functions
rm(list=ls())
setwd(dirname(getActiveDocumentContext()$path))
source("functions.R")
ncores <- detectCores()


##----------------- Example 1 ----------------#####

####### Simulations ###########
n <- 100  
p <- 10    
rho <- 0.5
sig <- ar1cor(p,rho)

b1 <- c(1,1,rep(0,p-2))
beta <- b1 
d <- 1

nsim <- 500 # no. of simulation runs
fnorm <- mclapply(1:nsim, function(q){
  
  set.seed(q)
  X <- rmvnorm(n, rep(0,p), sig)
  
  eps <- rnorm(n,0,sd=0.5)
  eps[c(10,15)] <- rnorm(2,5,sd=2)
  Y <- 0.5 + X%*%b1 + eps
  
  Bols <- coef(lm(Y ~ X))[-1]
  Bsir <- dr(Y~X, method="sir", nslices=5)$evectors[,1]
  
  Dy <- as.matrix(dist(Y, method="euclidean", diag=TRUE, upper=TRUE)) 
 
  MCD <- mcd(X,Dy)
  idx <- which(MCD$MCD >= 4/(n-(p+1)))
  
  if(length(idx)==0){
    Xi <- X 
    Yi <- Y
  } else{
    Xi <- X[-idx,]
    Yi <- Y[-idx]
  }
  
  bols <- coef(lm(Yi ~ Xi))[-1]
  bsir <- dr(Yi~Xi, method="sir", nslices=5)$evectors[,1]
  
  f1 <- theta(Bols,beta)
  f2 <- theta(bols,beta)
  f3 <- theta(Bsir,beta)
  f4 <- theta(bsir,beta)
  
  return(c(f1,f2,f3,f4))
}, mc.cores = 1)


ff <- matrix(unlist(fnorm),nsim,byrow=TRUE)
results <- rbind(apply(ff, 2, mean),apply(ff,2,sd)/sqrt(nsim) )
colnames(results) <- c("OLS","OLS*","SIR","SIR*")
rownames(results) <- c("Mean","Sd")
round(results,4)



##### --------------- Example 2 -------------###
n <- 100  
p <- 10   
rho <- 0.5
sig <- ar1cor(p,rho)

b1 <- c(1,-1,rep(0,p-2))
q=2
r=1

nsim <- 500 # no. of simulation runs

fnorm <- mclapply(1:nsim, function(m){
  
  set.seed(m)
  B = cbind(b1)%*%rbind(rortho(q)[1:r,])
  X <- rmvnorm(n, rep(0,p), sig)
  
  sig_e <- matrix(c(1,0.7,0.7,1),nrow=2,byrow=T)
  eps <- rmvnorm(n, rep(0,q), sig_e)
  eps[10,1] <- rnorm(1,5,2)
  eps[15,2] <- rnorm(1,-5,2)
  
  Y <- sin(X%*%B) + eps
  
  Bols <- rrr(Y,X,maxrank = 2)$coef[,1] 
  Bsir <- dr(Y~X, method="sir", nslices=5)$evectors[,1]
  
  Dy <- as.matrix(dist(Y, method="euclidean", diag=TRUE, upper=TRUE)) 
  
  MCD <- mcd(X,Dy)
  idx <- which(MCD$MCD >= 4/(n-(p+1)))
  
  if(length(idx)==0){
    Xi <- X 
    Yi <- Y
  } else{
    Xi <- X[-idx,]
    Yi <- Y[-idx,]
  }
  
  bols <- rrr(Yi,Xi,maxrank=2)$coef[,1] 
  bsir <- dr(Yi~Xi, method="sir", nslices=5)$evectors[,1]
  
  f1 <- theta(Bols,b1)
  f2 <- theta(bols,b1)
  f3 <- theta(Bsir,b1)
  f4 <- theta(bsir,b1)
  
  return(c(f1,f2,f3,f4))
},mc.cores = 1)


ff <- matrix(unlist(fnorm),nsim,byrow=TRUE)
results <- rbind(apply(ff, 2, mean),apply(ff,2,sd)/sqrt(nsim) )
colnames(results) <- c("OLS","OLS*","SIR","SIR*")
rownames(results) <- c("Mean","Sd")
round(results,4)





##-------------- Example 3 -------------##
n <- 100 
p <- 10 
m <- 50

b1 <- c(1,-1,rep(0,p-2))/sqrt(2)

nsim <- 500 # no. of simulation runs
fnorm <- mclapply(1:nsim, function(q){
  
  set.seed(q)
  X <- matrix(NA,n,p)
  for(h in 1:n) X[h,] <- runif(p,0,1)
  
  Y1 <- list()
  for(v in 1:n ) Y1[[v]] <- runif(m,-5*abs(X[v,]%*%b1),5*abs(X[v,]%*%b1))
  Y1[[10]] <-  rexp(m,abs(X[10,]%*%b1)) 
  Y1[[15]] <-  rnormMix(m,mean1=0,sd1=abs(X[15,]%*%b1),
                        mean2=(X[15,]%*%b1),sd2=5,p.mix = 0.5)
  
  Dy <- wdist(Y1)
  
  FOLS <- fols(X,Dy)$dir[,1]
  SAOLS <- sa_ols(X,Dy,N=1000)$dir[,1]
  FSIR <- fsir(X,Dy,h=5)$dir[,1]
  SASIR <- sa_sir(X,Dy,h=5,N=1000)$dir[,1]
  
  ## Find outliers 
  MCD <- mcd(X,Dy)
  idx <- which(MCD$MCD >= 4/(n-(p+1)))
  
  if(length(idx)==0){
    Xi <- X 
    Yi <- Y1
  } else{
    Xi <- X[-idx,]
    Yi <- Y1[-idx]
  }
  
  Dyi <- Dy[-idx,-idx]
  
  fols1 <- fols(Xi,Dyi)$dir[,1]
  saols1 <- sa_ols(Xi,Dyi,N=1000)$dir[,1]
  fsir1 <- fsir(Xi,Dyi,h=5)$dir[,1]
  sasir1 <- sa_sir(Xi,Dyi,N=1000)$dir[,1]
  
  f1 <- theta(FOLS,b1)
  f2 <- theta(fols1,b1)
  f3 <- theta(SAOLS,b1)
  f4 <- theta(saols1,b1)
  f5 <- theta(FSIR,b1)
  f6 <- theta(fsir1,b1)
  f7 <- theta(SASIR,b1)
  f8 <- theta(sasir1,b1)
  
  return(c(f1,f2,f3,f4,f5,f6,f7,f8))
},mc.cores = 1)

ff <- matrix(unlist(fnorm),nsim,byrow=TRUE)
results <- rbind(apply(ff, 2, mean),apply(ff,2,sd)/sqrt(nsim) )
colnames(results) <- c("fOLS","fOLS*","saOLS","saOLS*","fSIR","fSIR*","saSIR","saSIR*")
rownames(results) <- c("Mean","Sd")
round(results,4)



##---------------- Example 4  ----------############
n <- 100 
p <- 10 
m <- 30
rho <- 0
sig <- ar1cor(p,rho)

b1 <- c(1,-1,rep(0,p-2))/sqrt(2)

nsim <- 500  # no. of simulation runs
fnorm <- mclapply(1:nsim, function(q){
  
  set.seed(q)
  X <- rmvnorm(n, rep(0,p), sig)
  
  TT <- sort(round(runif(m,0,10), 4)) # Functional data support
  mu <- 10*cos(pi + pi*TT/5)
  mu_t <- matrix(mu,n,m, byrow=TRUE)
  
  Y <- list()
  for (i in seq_len(n)) {
    Y[[i]] <- mu_t[i,] + c(5*sin(X[i,]%*%b1)) + rnorm(1)
  }
  Y[[10]][seq(5,25,5)] <- Y[[10]][seq(5,25,5)] + c(5*exp(2+X[seq(5,25,5),]%*%b1))
  Y[[15]] <- mu_t[15,] + c(5*sin(X[15,]%*%b1)) - 2*rnorm(m,0,10)
  
  Dy <- tsdist(Y)
  
  SAOLS <- sa_ols(X,Dy,N=1000)$dir[,1]
  SASIR <- sa_sir(X,Dy,h=5,N=1000)$dir[,1]
  
  MCD <- mcd(X,Dy)
  idx <- which(MCD$MCD >= 4/(n-(p+1)))
  
  if(length(idx)==0){
    Xi <- X 
    Yi <- Y
  } else{
    Xi <- X[-idx,]
    Yi <- Y[-idx]
  }
  
  Dyi <- Dy[-idx,-idx]
  saols1 <- sa_ols(Xi,Dyi,N=1000)$dir[,1]
  sasir1 <- sa_sir(Xi,Dyi,h=5,N=1000)$dir[,1]
  
  f1 <- theta(SAOLS,b1)
  f2 <- theta(saols1,b1)
  f3 <- theta(SASIR,b1)
  f4 <- theta(sasir1,b1)
  
  return(c(f1,f2,f3,f4))
}, mc.cores = 1)

ff <- matrix(unlist(fnorm),nsim,byrow=TRUE)
results <- rbind(apply(ff, 2, mean),apply(ff,2,sd)/sqrt(nsim) )
colnames(results) <- c("saOLS","saOLS*","saSIR","saSIR*")
rownames(results) <- c("Mean","Sd")
round(results,4)

