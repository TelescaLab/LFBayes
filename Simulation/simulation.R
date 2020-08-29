### Call this file in terminal by Rscript simulation.R SEED
args = commandArgs(trailingOnly = TRUE)
args <- as.numeric(args[1])
setwd("/Users/johnshamshoian/Documents/R_projects/LFBayes/Simulation")
if(!dir.exists("output")) dir.create("output")
library(splines)
library(MASS)
library(LFBayes)
library(pracma)
library(dplyr)
if(args[1] <= 500){
  q1s <- 3
  q2s <- 3
  splinenum <- 10
}
if(args[1] > 500 & args[1] <= 1000){
  q1s <- 3
  q2s <- 3
  splinenum <- 12
}
if(args[1] > 1000 & args[1] <= 1500){
  q1s <- 3
  q2s <- 3
  splinenum <- 14
}
if(args[1] > 1500 & args[1] <= 2000){
  q1s <- 3
  q2s <- 3
  splinenum <- 16
}
if(args[1] > 2000 & args[1] <= 2500){
  q1s <- 6
  q2s <- 6
  splinenum <- 10
}
if(args[1] > 2500 & args[1] <= 3000){
  q1s <- 6
  q2s <- 6
  splinenum <- 12
}
if(args[1] > 3000 & args[1] <= 3500){
  q1s <- 6
  q2s <- 6
  splinenum <- 14
}
if(args[1] > 3500 & args[1] <= 4000){
  q1s <- 6
  q2s <- 6
  splinenum <- 16
}

set.seed(args %% 500)
cat("Seed is", args %% 500, "splinenum is", splinenum, "q is", q1s)
alpha <- 0.1
errorvar <- .025
SS <- 20
TT <- 20
t <- seq(from = 0, to = 1, length.out = TT)
s <- seq(from = 0, to = 1, length.out = SS)
n <- 60
tt <- list()
tt[[1]] <- 1:(TT*SS)
tt <- rep(tt, n)
p1 <- 12
p2 <- 12
q1 <- 3
q2 <- 3
Bt <- bs(t, df = p1, intercept = TRUE)
Bs <- bs(s, df = p2, intercept = TRUE)

Loading.Brown.Bridge <- function(t, p, k){
  B <- bs(t, df = p, intercept = TRUE)
  Loading <- matrix(nrow = p, ncol = k)
  for(i in 1:k){
    eigval <- 1/(i^2*pi^2)
    #if(i %% 2 == 0){
    #  psi <- sqrt(2) * -sqrt(2)*cos((i+1)*pi*t)
    #}else{
    psi <- sqrt(2)*sin(i*pi*t)
    #}
    
    Loading[,i] <- sqrt(eigval)*solve(t(B)%*%B)%*%t(B)%*%psi
  }
  return(Loading)
}

Matern.Cov <- function(s){
  rho <- 0.5
  sigmasq <- 1
  Matern.Cov <- matrix(nrow = length(s), ncol = length(s))
  for(i in 1:length(s)){
    for(j in 1:length(s)){
      d <- abs(s[i] - s[j])
      Matern.Cov[i, j] <- sigmasq * (1 + sqrt(3) * d / rho) * exp(- sqrt(3) * d / rho)
      #Matern.Cov[i, j] <- (1 + sqrt(5) * d / rho + 5 * d^2 / (3*rho^2)) * exp(-sqrt(5) * d / rho)
    }
  }
  Matern.Cov
}

Loading.Matern <- function(s, p, k, B){
  rho <- 0.5
  sigmasq <- 1
  Matern.Cov <- matrix(nrow = length(s), ncol = length(s))
  for(i in 1:length(s)){
    for(j in 1:length(s)){
      d <- abs(s[i] - s[j])
      Matern.Cov[i, j] <- sigmasq * (1 + sqrt(3) * d / rho) * exp(- sqrt(3) * d / rho)
    }
  }
  Loading <- matrix(nrow = p, ncol = k)
  evec <- eigen(Matern.Cov)$vectors
  eval <- eigen(Matern.Cov)$values
  for(i in 1:k){
    Loading[,i] <- sqrt(eval[i]) * solve(t(B)%*%B)%*%t(B)%*%evec[,i]
  }
  Loading
}


GenerateH <- function(q1,q2){
  H <- matrix(0,q1, q2)
  for(i in 1:q1){
    for(j in 1:q2){
      H[i,j] <- exp(-(sqrt(.01*i) + sqrt(.1*j)))
    }
  }
  H <- diag(c(H))
  H
}

GenerateMu1 <- function(s,t){
  mu <- matrix(0, nrow= length(t),ncol=length(s))
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      mu[i,j] <- sqrt(1/(5*sqrt(s[j]+1)))*sin(5*t[i])
    }
  }
  c(mu)
}
H <- GenerateH(q1, q2)

mu1 <- GenerateMu1(s,t)
Lambda <- Loading.Matern(t, p1, q1, Bt)
Gamma <- Loading.Brown.Bridge(s, p2, q2)
Cov <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar * diag(SS * TT)
Bt1 <- bs(t, df = splinenum, intercept = TRUE)
Bs1 <- bs(s, df = splinenum, intercept = TRUE)

iter <- 5000 # Number of iterations
burnin <- 1000 # Burn-in iterations
thin <- 1 # Thinning for each chain
nchain <- 1 # Number of chains
neig <- 3 # Number of eigenfunctions for inference

x <- mvrnorm(n, mu  = as.vector(mu1), Sigma = Cov)
sx <- sd(x)
mx <- mean(x)
x <- (x-mx)/sx
Smooth_scaled_cov <- (Cov - errorvar * diag(SS * TT)) / sx^2
mu <- (mu1 - mx)/sx
Marg.Long <- getMarginalLong(Smooth_scaled_cov,SS,TT)
Marg.Func <- getMarginalFunc(Smooth_scaled_cov,SS,TT)
m1 <- eigen(Marg.Long)$vectors[,1:3] # Marginal longitudinal eigenvectors
m2 <- eigen(Marg.Func)$vectors[,1:3] # Marginal functional eigenvectors
y <- lapply(1:n, function(i) x[i,])
missing <- list()
for(ii in 1:n){
  missing[[ii]] <- numeric(0)
}
X <- rep(1,n)
dim(X) <- c(n,1)
cat("\n")
MCMC <- mcmcWeakChains(y, missing, X, Bs1, Bt1,
                       q1s, q2s, iter, thin, burnin, nchain)
MCMC_eigen <- eigenLFChains(Bs1, Bt1, MCMC, neig,
                            iter, burnin, nchain, s, t, alpha)

signLong <- rep(1,3)
signLong[1] <- if(sum((MCMC_eigen$eigvecLongmean[,1] + m1[,1])^2) < sum((MCMC_eigen$eigvecLongmean[,1] - m1[,1])^2)) -1 else 1
signLong[2] <- if(sum((MCMC_eigen$eigvecLongmean[,2] + m1[,2])^2) < sum((MCMC_eigen$eigvecLongmean[,2] - m1[,2])^2)) -1 else 1
signLong[3] <- if(sum((MCMC_eigen$eigvecLongmean[,3] + m1[,3])^2) < sum((MCMC_eigen$eigvecLongmean[,3] - m1[,3])^2)) -1 else 1

signFunc <- rep(1,3)
signFunc[1] <- if(sum((MCMC_eigen$eigvecFuncmean[,1] + m2[,1])^2) < sum((MCMC_eigen$eigvecFuncmean[,1] - m2[,1])^2)) -1 else 1
signFunc[2] <- if(sum((MCMC_eigen$eigvecFuncmean[,2] + m2[,2])^2) < sum((MCMC_eigen$eigvecFuncmean[,2] - m2[,2])^2)) -1 else 1
signFunc[3] <- if(sum((MCMC_eigen$eigvecFuncmean[,3] + m2[,3])^2) < sum((MCMC_eigen$eigvecFuncmean[,3] - m2[,3])^2)) -1 else 1

# Evaluate
results <- numeric(15)
results[1] <- sum((mu - as.numeric(MCMC_eigen$postmean))^2) / sum(mu^2)
results[2] <- all(mu < MCMC_eigen$upper & mu > MCMC_eigen$lower)

# Eignfn in longitudinal direction
results[3] <- sum((m1[,1] - signLong[1] * MCMC_eigen$eigvecLongmean[,1])^2)
results[4] <- sum((m1[,2] - signLong[2] * MCMC_eigen$eigvecLongmean[,2])^2)
results[5] <- sum((m1[,3] - signLong[3] * MCMC_eigen$eigvecLongmean[,3])^2)

if(signLong[1] == -1){
  results[6] <- all(m1[,1] > -MCMC_eigen$eigvecLongupper[,1] & m1[,1] < -MCMC_eigen$eigvecLonglower[,1])
} else{
  results[6] <- all(m1[,1] < MCMC_eigen$eigvecLongupper[,1] & m1[,1] > MCMC_eigen$eigvecLonglower[,1])
}
if(signLong[2] == -1){
  results[7] <- all(m1[,2] > -MCMC_eigen$eigvecLongupper[,2] & m1[,2] < -MCMC_eigen$eigvecLonglower[,2])
} else{
  results[7] <- all(m1[,2] < MCMC_eigen$eigvecLongupper[,2] & m1[,2] > MCMC_eigen$eigvecLonglower[,2])
}
if(signLong[3] == -1){
  results[8] <- all(m1[,3] > -MCMC_eigen$eigvecLongupper[,3] & m1[,3] < -MCMC_eigen$eigvecLonglower[,3])
} else{
  results[8] <- all(m1[,3] < MCMC_eigen$eigvecLongupper[,3] & m1[,3] > MCMC_eigen$eigvecLonglower[,3])
}

# Eignfn in functional direction
results[9] <- sum((m2[,1] - signFunc[1] * MCMC_eigen$eigvecFuncmean[,1])^2)
results[10] <- sum((m2[,2] - signFunc[2] * MCMC_eigen$eigvecFuncmean[,2])^2)
results[11] <- sum((m2[,3] - signFunc[3] * MCMC_eigen$eigvecFuncmean[,3])^2)

if(signFunc[1] == -1){
  results[12] <- all(m2[,1] > -MCMC_eigen$eigvecFuncupper[,1] & m2[,1] < -MCMC_eigen$eigvecFunclower[,1])
} else{
  results[12] <- all(m2[,1] < MCMC_eigen$eigvecFuncupper[,1] & m2[,1] > MCMC_eigen$eigvecFunclower[,1])
}
if(signFunc[2] == -1){
  results[13] <- all(m2[,2] > -MCMC_eigen$eigvecFuncupper[,2] & m2[,2] < -MCMC_eigen$eigvecFunclower[,2])
} else{
  results[13] <- all(m2[,2] < MCMC_eigen$eigvecFuncupper[,2] & m2[,2] > MCMC_eigen$eigvecFunclower[,2])
}
if(signFunc[3] == -1){
  results[14] <- all(m2[,3] > -MCMC_eigen$eigvecFuncupper[,3] & m2[,3] < -MCMC_eigen$eigvecFunclower[,3])
} else{
  results[14] <- all(m2[,3] < MCMC_eigen$eigvecFuncupper[,3] & m2[,3] > MCMC_eigen$eigvecFunclower[,3])
}

# result 1: relative squared error of mean
# result 2: was mean covered?
# result 3: squared error of first longitudinal eigenfunction
# result 4: squared error of second longitudinal eigenfunction
# result 5: squared error of third longitudinal eigenfunction
# result 6: was first longitudinal eigenfunction covered?
# result 7: was second longitudinal eigenfunction covered?
# result 8: was third longitudinal eigenfunction covered?
# result 9: squared error of first functional eigenfunction
# result 10: squared error of second functional eigenfunction
# result 11: squared error of third functional eigenfunction
# result 12: was first functional eigenfunction covered?
# result 13: was second functional eigenfunction covered?
# result 14: was third functional eigenfunction covered?
save(results, file = paste0("output/results", args %% 500, ".RData"))

