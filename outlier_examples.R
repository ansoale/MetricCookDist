library(mvtnorm)
library(fdadensity)
library(transport)
library(pracma)
library(dr)
library(dplyr)
library(frechet)
require(tictoc)
require(MASS)
require(magrittr)
require(dplyr)
library(parallel)
library(ggplot2)
library(tidyr)
library(stringr)
library(statip)
library(plotly)
library(latex2exp)
library(gridExtra)
library(tictoc)
library(nlme)
library(fda)
library(rstudioapi)
library(rrpack)
library(EnvStats)
library(devtools)
library(roxygen2)

## Set working directory and call functions
rm(list=ls())
setwd(dirname(getActiveDocumentContext()$path))
source("functions.R")
ncores <- detectCores()


##--------------- Example 1 ------------------#####
## set parameters
set.seed(1)
n <- 100 
p <- 10    
rho <- 0.5
sig <- ar1cor(p,rho)
b1 <- c(1,1,rep(0,p-2))
beta <- cbind(b1)
d <- 1

## generate X
X <- rmvnorm(n, rep(0,p), sig)

## generate errors
set.seed(1)
eps <- rnorm(n,0,sd=0.5)
eps[c(10,15)] <- rnorm(2,5,sd=2)

## generate Y
Y <- 0.5 + X%*%b1 + eps

## calculate pairwise distances
Dy <- as.matrix(dist(Y,method="euclidean",diag=TRUE,upper=TRUE))

## Calculate surrogate response and metric Cook's distances
MCD <- mcd(X, Dy)

## OLS estimate based on the surrogate response
bhat <- coef(lm(MCD$Sy~X))[-1]


## Plot Cook's distances
cols <- rep(1,n)
cols[c(10,15)] <- 2

par(mfrow=c(2,2))
plot(Y~X%*%b1, xlab=expression(X~beta),col=cols, main="(A)")
plot(MCD$Sy~X%*%bhat, xlab=expression(X~hat(beta)), ylab=expression(S^Y),
     main="(B)", col=cols, pch=18)
plot(lm(Y ~ X), 4, main="(C)")
plot(mcd~i, data=MCD$MCD.index, type="h", ylab="Cook's distance",
     xlab="Obs. number", main="(D) Metric Cook's distance", ylim=c(0,1))
text(x=MCD$MCD.index$i[c(10,15)],y=MCD$MCD.index$mcd[c(10,15)],
     labels=c(10,15),pos=3)




##### -------------- Example 2 -----------###
set.seed(1)
n <- 100 
p <- 10   
rho <- 0.5
sig <- ar1cor(p,rho)
X <- rmvnorm(n, rep(0,p), sig)

b1 <- c(1,-1,rep(0,p-2))
q=2
r=1
B = cbind(b1)%*%rbind(rortho(q)[1:r,])

sig_e <- matrix(c(1,0.7,0.7,1),nrow=2,byrow=T)
eps <- rmvnorm(n, rep(0,q), sig_e)
eps[10,1] <- rnorm(1,5,2)
eps[15,2] <- rnorm(1,-5,2)
 
Y <- sin(X%*%B) + eps


## 3D scatter plot of responses vs sufficient predictor
Ydf <- data.frame(y=Y,x=X%*%b1)
colnames(Ydf) <- c("y1","y2","x")
Ydf$col <- "n"
Ydf$col[c(10,15)] <- "y"

fig <- plot_ly(Ydf, x=~x, y=~y1, z=~y2, color=~col, colors=c('black','#BF382A'), 
                type='scatter3d', mode="markers") %>%
  layout(scene = list(
    xaxis=list(title="\u03B2<sup>T</sup>x",zerolinewidth=0),
    yaxis=list(title="y<sub>1</sub>",zerolinewidth=0),
    zaxis=list(title="y<sub>2</sub>",zerolinewidth=0))) 

fig <- hide_colorbar(fig) %>% layout( showlegend = TRUE, 
                                        legend=list(title=list(text='<b> Outlier </b>')))
fig

## Calculate pairwise distances 
Dy <- as.matrix(dist(Y, method="euclidean", diag=TRUE, upper=TRUE))

## Calculate surrogate response and metric Cook's distances
MCD <- mcd(X,Dy)

## Plot Cook's distances
par(mfrow=c(2,2))
plot(lm(Y[,1] ~ X), 4, main="(A)")
plot(lm(Y[,2] ~ X), 4, main="(B)")
plot(mcd~i, data=MCD$MCD.index, type="h", ylab="Cook's distance",
     xlab="Obs. number", main="(c) Metric Cook's distance", ylim=c(0, 0.3))
text(x=MCD$MCD.index$i[c(10,15)],y=MCD$MCD.index$mcd[c(10,15)],
     labels=c(10,15),pos=3)



##---------------- Example 3 ----------------##
n <- 100  
p <- 10  
m <- 50

b1 <- c(1,-1,rep(0,p-2))/sqrt(2)
d <- 1

set.seed(123)
X <- matrix(NA,n,p)
for(h in 1:n) X[h,] <- runif(p,0,1)

Y <- list()
for(v in 1:n ) Y[[v]] <- runif(m,-5*abs(X[v,]%*%b1),5*abs(X[v,]%*%b1))
Y[[10]] <-  rexp(m,abs(X[10,]%*%b1)) 
Y[[15]] <-  rnormMix(m,mean1=0,sd1=abs(X[15,]%*%b1),
              mean2=(X[15,]%*%b1),sd2=5,p.mix = 0.5)

### plot densities for each response
dens0 <- lapply(Y, density)
dens <- data.frame(
  x = unlist(lapply(dens0, "[[", "x")),
  y = unlist(lapply(dens0, "[[", "y"))
)

xb <- round(c(X%*%b1),4)
cols <- rep("n",100)
cols[c(10,15)] <- "y"

Pred <-  xb[rep(seq_len(length(xb)), each=length(dens0[[1]]$x))]
Cols <-  cols[rep(seq_len(length(cols)), each=length(dens0[[1]]$x))]
dens_df <- cbind(dens,Pred,Cols)

axx <- list(title = "y")
axy <- list(title = "density")
axz <- list(title = list(text='\u03B2<sup>T</sup>x'))

fig <- plot_ly(dens_df, x=~x, y=~y, z=~Pred, split=~Pred, color=~Cols, 
               colors=c('black','#BF382A'),type='scatter3d', mode='lines') %>% 
              layout(scene=list(xaxis=axx,yaxis=axy,zaxis=axz))
fig <- hide_colorbar(fig) %>% layout(showlegend = FALSE)
fig


## Compute pairwise distances 
Dy <- wdist(Y)

## Calculate surrogate response and metric Cook's distances
MCD <- mcd(X,Dy)

## Plot Cook's distances
par(mfrow=c(1,1))
plot(mcd~i, data=MCD$MCD.index, type="h", ylab="Cook's distance",
     xlab="Obs. number", main="(D) Metric Cook's distance", ylim=c(0,1))
text(x=MCD$MCD.index$i[c(10,15)],y=MCD$MCD.index$mcd[c(10,15)],
     labels=c(10,15),pos=3)



##--------------- Example 4 -------------------############
set.seed(123)
n <- 100
p <- 10
rho <- 0
sig <- ar1cor(p,rho)
X <- rmvnorm(n, rep(0,p), sig)

b1 <- c(1,-1,rep(0,p-2))/sqrt(2)
m <- 30
TT <- sort(round(runif(m,0,10), 4)) # Functional data support
mu <- 10*cos(pi + pi*TT/5)
mu_t <- matrix(mu,n,m, byrow=TRUE)

Y <- list()
for (i in seq_len(n)) {
  Y[[i]] <- mu_t[i,] + c(5*sin(X[i,]%*%b1)) + rnorm(1)
}
Y[[10]][seq(5,25,5)] <- Y[[10]][seq(5,25,5)] + c(5*exp(2+X[seq(5,25,5),]%*%b1))
Y[[15]] <- mu_t[15,] + c(5*sin(X[15,]%*%b1)) - 2*rnorm(m,0,10)


## Plot response functions
matplot(do.call(cbind,Y),type='l',xlab="times (t)", ylab=expression(y(t)))

## Compute pairwise distances
Dy <- tsdist(Y)

## Calculate surrogate response and metric Cook's distances
MCD <- mcd(X,Dy)

## Plot Cook's distances
par(mfrow=c(1,2))
matplot(do.call(cbind,Y),type='l',xlab="times (t)", 
        main="Functional responses", ylab=expression(y(t)))
plot(mcd~i, data=MCD$MCD.index, type="h", ylab="Cook's distance",
     xlab="Obs. number", main="(D) Metric Cook's distance", ylim=c(0,1))
text(x=MCD$MCD.index$i[c(10,15)],y=MCD$MCD.index$mcd[c(10,15)],
     labels=c(10,15),pos=3)
