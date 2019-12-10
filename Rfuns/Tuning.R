library(splines)
library(MASS)
library(fdapace)
library(doParallel)
library(doSNOW)
library(coda)
library(MCMCglmm)
library(LFBayes)
#setwd("/Users/John/Downloads/LongFunc Code/ChenCode")
#setwd("/Users/John/Documents/Johnstuff/LFBayes/Rfuns")
setwd("E:/Rcpp stuff/LFBayes/LFBayes/Rfuns")
source("MarginalFPCA.R")
source("ProductFPCA.R")

errorvar <- .05
t <- seq(from = 0, to = 1, length.out = 20)
s <- seq(from = 0, to = 1, length.out = 10)
n <- 30
tt <- list()
tt[[1]] <- 1:(length(t)*length(s))
tt <- rep(tt, n)
p1 <- 10
p2 <- 5
q1 <- 5
q2 <- 5
Bt <- bs(t, df = p1, intercept = TRUE)
Bs <- bs(s, df = p2, intercept = TRUE)
Bt1 <- bs(t, df = p1+5, intercept = TRUE)
Bs1 <- bs(s, df = p2+3, intercept = TRUE)
productcov <- function(eig, scores){
  n <- dim(scores)[1]
  neig <- dim(eig)[2]
  v <- numeric(dim(eig)[1])
  mycov <- matrix(0, nrow = length(v), ncol = length(v))
  for(l in 1:neig){
    v = eig[,l]
    mycov <- mycov + outer(v, v) * 1/n * as.numeric(scores[,l]%*%scores[,l])
  }
  mycov
}
marginalcov <- function(eig, scores){
  n <- dim(scores)[1]
  neig <- dim(eig)[2]
  v <- numeric(dim(eig)[1])
  mycov <- matrix(0, nrow = length(v), ncol = length(v))
  for(j in 1:neig){
    v <- eig[,j]
    mycov <- mycov + outer(v,v) * 1/n * as.numeric(scores[,j]%*%scores[,j])
  }
  
  
  mycov
}
Brown.Bridge.Coef <- function(B, t, l){
  eigval <- 1/(l^2*pi^2)
  psi <- sqrt(2)*sin(l*pi*t)
  return(sqrt(eigval)*solve(t(B)%*%B)%*%t(B)%*%psi)
}
Brown.Motion.Coef <- function(B, s, l){
  eigval <- 1/((l - 1/2)^2*pi^2)
  psi <- sqrt(2)*sin((l-1/2)*pi*s)
  return(sqrt(eigval)*solve(t(B)%*%B)%*%t(B)%*%psi)
}
Loading.Brown.Bridge <- function(t, p, k){
  B <- bs(t, df = p, intercept = TRUE)
  Loading <- matrix(nrow = p, ncol = k)
  for(i in 1:k){
    Loading[,i] <- Brown.Bridge.Coef(B, t, i)
  }
  return(Loading)
}
Loading.Brown.Motion <- function(t, p, k){
  B <- bs(t, df = p, intercept = TRUE)
  Loading <- matrix(nrow = p, ncol = k)
  for(i in 1:k){
    Loading[,i] <- Brown.Motion.Coef(B, t, i)
  }
  return(Loading)
}
Matern.Func <- function(t1, t2){
  d <- abs(t1 - t2)
  rho = 0.5
  return((1 + sqrt(5) * d / rho + 5 * d^2 / (3*rho^2)) * exp(-sqrt(5) * d / rho))
}
Matern.Cov <- function(s){
  Matern.Cov <- matrix(nrow = length(s), ncol = length(s))
  for(i in 1:length(s)){
    for(j in 1:length(s)){
      Matern.Cov[i, j] <- Matern.Func(s[i], s[j])
    }
  }
  Matern.Cov
}

Loading.Matern <- function(s, p, k, B){
  Loading <- matrix(nrow = p, ncol = k)
  Matern <- Matern.Cov(s)
  evec <- eigen(Matern)$vectors
  eval <- eigen(Matern)$values
  for(i in 1:k){
    Loading[,i] <- sqrt(eval[i]) * solve(t(B)%*%B)%*%t(B)%*%evec[,i]
  }
  Loading
}
Lambda <- Loading.Matern(t, p1, q1, Bt)
#Lambda <- Loading.Brown.Bridge(t, p1, q1)
Gamma <- Loading.Brown.Motion(s, p2, q2)
Brown.Motion.Cov <- Bs%*%Gamma%*%t(Gamma)%*%t(Bs)
Matern.Cov <- Bt%*%Lambda%*%t(Lambda)%*%t(Bt)
Cov.Strong <- kronecker(Brown.Motion.Cov, Matern.Cov)
#H <- diag(rgamma(q1*q2, shape = 1, rate = 1))
H <- diag(5*5)
setwd("E:/Rcpp stuff/LFBayes/LFBayes/Rfuns")
load("H.RData")
Cov.Weak <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar*diag(length(s) * length(t))
pc.j = NULL
pc.k = NULL
fpca.op1 = list(dataType = "Sparse", maxK = pc.j, FVEthreshold = .9999, nRegGrid = 20)
fpca.op2 = list(dataType = "Sparse", maxK = pc.k, FVEthreshold = .9999, nRegGrid = 10)
mu1 <- kronecker(sqrt(1/(5*sqrt(s)+1)), sin(5*t))
image(Cov.Weak[,200:1], col = heat.colors(100))

iterations <- 50
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(4)
registerDoSNOW(cl)
system.time(n30_s.05_WSep_HSep_Phi<-foreach(index=1:iterations,.combine=rbind, .packages = c("MASS", "LFBayes", "fdapace"), .options.snow = opts)%dopar%{
  print(index)
  #x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
  x <- mvrnorm(n, mu  = as.vector(mu1), Sigma = Cov.Weak)
  sx <- sd(x)
  mx <- mean(x)
  x <- (x-mx)/sx
  Smooth_scaled_cov <- (Cov.Weak - errorvar * diag(length(s) * length(t))) / sx^2
  mu <- (mu1 - mx)/sx
  y <- lapply(1:n, function(i) x[i,])
  missing <- list()
  for(i in 1:n){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,n)
  dim(X) <- c(n,1)
  
  q1 <- 8
  q2 <- 8
  mcmc <- mcmcWeak(y, missing, X, Bs1, Bt1, q1, q2, 25000, 1, 5000)
  
  resMarginal <- MarginalFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
  resProduct <- ProductFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resProductCov <- productcov(resProduct$eig, resProduct$scores)
  resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
  yt <- t(x)
  EmpCov <- 1 / (n-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
  EmpMean <- colMeans(x)
  results <- numeric(19)
  
  results[1] <- max(abs(eigen(Smooth_scaled_cov - resProductCov)$values)) / (length(s) * length(t))
  results[2] <- max(abs(eigen(Smooth_scaled_cov - resMarginalCov)$values)) / (length(s) * length(t))
  results[3] <- max(abs(eigen(Smooth_scaled_cov - resPACE$smoothedCov)$values)) / (length(s) * length(t))
  results[4] <- max(abs(eigen(Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(length(s) * length(t))))$values)) / (length(s) * length(t))
  results[5] <- max(abs(eigen(Smooth_scaled_cov - mcmc$postcov)$values)) / (length(s) * length(t))

  m1 <- eigen(Brown.Motion.Cov)$vectors[,1:3]
  results[6] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
  results[7] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
  results[8] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
  results[9] <-min(sum((mcmc$eigvecLongmean[,3] - m1[,1])^2), sum((mcmc$eigvecLongmean[,3] + m1[,1])^2))
  results[10] <-min(sum((mcmc$eigvecLongmean[,2] - m1[,2])^2), sum((mcmc$eigvecLongmean[,2] + m1[,2])^2))
  results[11] <-min(sum((mcmc$eigvecLongmean[,1] - m1[,3])^2), sum((mcmc$eigvecLongmean[,1] + m1[,3])^2))
  
  m2 <- eigen(Matern.Cov)$vectors[,1:3]
  
  results[12] <-min(sum((resProduct$psi[,1] - m2[,1])^2), sum((resProduct$psi[,1] + m2[,1])^2))
  results[13] <-min(sum((resProduct$psi[,2] - m2[,2])^2), sum((resProduct$psi[,2] + m2[,2])^2))
  results[14] <-min(sum((resProduct$psi[,3] - m2[,3])^2), sum((resProduct$psi[,3] + m2[,3])^2))
  results[15] <-min(sum((mcmc$eigvecFuncmean[,3] - m2[,1])^2), sum((mcmc$eigvecFuncmean[,3] + m2[,1])^2))
  results[16] <-min(sum((mcmc$eigvecFuncmean[,2] - m2[,2])^2), sum((mcmc$eigvecFuncmean[,2] + m2[,2])^2))
  results[17] <-min(sum((mcmc$eigvecFuncmean[,1] - m2[,3])^2), sum((mcmc$eigvecFuncmean[,1] + m2[,3])^2))
  
  results[18] <- sum((mu - resPACE$mu)^2)/(length(s) * length(t))
  results[19] <- sum((mu - as.numeric(mcmc$postmean))^2)/(length(s) * length(t))
  gc()
  #save(results, file = paste("n30_s.05_WSep_H1_NoPen" ,index,".RData", sep = ""))
  results
  
})[3]
stopCluster(cl)
save(n30_s.05_WSep_H.1_, file = "n30_s.05_WSep_H1_NoPen.RData")