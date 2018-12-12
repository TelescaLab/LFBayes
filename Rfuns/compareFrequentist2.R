library(splines)
library(MASS)
library(fdapace)
library(doParallel)
library(doSNOW)
library(coda)
library(MCMCglmm)
library(LFBayes)
#setwd("/Users/John/Downloads/LongFunc Code/ChenCode")
setwd("/Users/John/Documents/Johnstuff/LFBayes/Rfuns")

source("MarginalFPCA.R")
source("ProductFPCA.R")

errorvar <- .05
t <- seq(from = 0, to = 1, length.out = 20)
s <- seq(from = 0, to = 1, length.out = 20)
n <- 90
tt <- list()
tt[[1]] <- 1:(length(t)*length(s))
tt <- rep(tt, n)
p1 <- 10
p2 <- 6
q1 <- 5
q2 <- 3
Bt <- bs(t, df = p1, intercept = TRUE)
Bs <- bs(s, df = p2, intercept = TRUE)
Bt1 <- bs(t, df = p1+2, intercept = TRUE)
Bs1 <- bs(s, df = p2+2, intercept = TRUE)
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

CosCov <- function(s){
  alpha <- 2
  CosCov <- matrix(0,nrow=length(s),ncol=length(s))
  for(k in 1:50){
    CosCov <- k^(-2*alpha)*outer(cos(k*pi*s),cos(k*pi*s)) + CosCov
  }
  CosCov
}

Loading.CosCov <- function(s, p, k, B){
  Loading <- matrix(nrow = p, ncol = k)
  CosCov <- CosCov(s)
  evec <- eigen(CosCov)$vectors
  eval <- eigen(CosCov)$values
  for(i in 1:k){
    Loading[,i] <- sqrt(eval[i]) * solve(t(B)%*%B)%*%t(B)%*%evec[,i]
  }
  Loading
}
Lambda <- Loading.Matern(t, p1, q1, Bt)
#Lambda <- Loading.Brown.Bridge(t, p1, q1)
#Gamma <- Loading.Brown.Bridge(s, p2, q2)
Gamma <- Loading.CosCov(s,p2,q2,Bs)
#Brown.Bridge.Cov <- Bs%*%Gamma%*%t(Gamma)%*%t(Bs)
CosCov.Recon <- Bs%*%Gamma%*%t(Gamma)%*%t(Bs)
Matern.Cov <- Bt%*%Lambda%*%t(Lambda)%*%t(Bt)
Cov.Strong <- kronecker(CosCov.Recon, Matern.Cov)
#H <- diag(rgamma(q1*q2, shape = 1, rate = 1))


#H <- diag(q1*q2)
#setwd("/Users/John/Documents/Johnstuff/LFBayes/Rfuns")
#load("H2.RData")
#H <- diag(q1*q2)
H <- diag(rgamma(q1*q2,.1,.1))
#H <- diag(sin(1:(q1*q2)))
Cov.Weak <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar*diag(length(s) * length(t))
pc.j = NULL
pc.k = NULL
fpca.op1 = list(dataType = "Sparse", maxK = pc.j, FVEthreshold = .9999, nRegGrid = length(t))
fpca.op2 = list(dataType = "Sparse", maxK = pc.k, FVEthreshold = .9999, nRegGrid = length(s))
mu1 <- kronecker(sqrt(1/(5*sqrt(s)+1)), sin(5*t))
image(Cov.Weak[,400:1], col = heat.colors(100))
#image(Cov.Strong[,200:1],col=heat.colors(100))  

iterations <- 100
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(5)
registerDoSNOW(cl)
system.time(n30_50<-foreach(index=1:iterations,.combine=rbind, .packages = c("MASS", "LFBayes", "fdapace"), .options.snow = opts)%dopar%{
  print(index)
  #x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
  x <- mvrnorm(n, mu  = as.vector(mu1), Sigma = Cov.Weak)
  sx <- sd(x)
  mx <- mean(x)
  x <- (x-mx)/sx
  Smooth_scaled_cov <- (Cov.Weak - errorvar * diag(length(s) * length(t))) / sx^2
  mu <- (mu1 - mx)/sx
  #scaled_cov <- Cov/sx^2
  y <- lapply(1:n, function(i) x[i,])
  missing <- list()
  for(i in 1:n){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,n)
  dim(X) <- c(n,1)

  q1 <- 8
  q2 <- 6
  mcmc <- mcmcWeak(y, missing, X, Bs1, Bt1, q1, q2, 5000, 5, 1000)

  resMarginal <- MarginalFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
  resProduct <- ProductFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resProductCov <- productcov(resProduct$eig, resProduct$scores)
  resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
  yt <- t(x)
  EmpCov <- 1 / (n-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
  EmpMean <- colMeans(x)
  results <- numeric(20)
  
  results[1] <- max(abs(eigen(Smooth_scaled_cov - resProductCov)$values)) / (length(s) * length(t))
  results[2] <- max(abs(eigen(Smooth_scaled_cov - resMarginalCov)$values)) / (length(s) * length(t))
  results[3] <- max(abs(eigen(Smooth_scaled_cov - resPACE$smoothedCov)$values)) / (length(s) * length(t))
  results[4] <- max(abs(eigen(Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(length(s) * length(t))))$values)) / (length(s) * length(t))
  results[5] <- max(abs(eigen(Smooth_scaled_cov - mcmc$postcov)$values)) / (length(s) * length(t))
  #results[6] <- max(abs(eigen(Smooth_scaled_cov - mcmc$postcov_SE)$values)) / (length(s) * length(t))
  m1 <- eigen(getMarginalLong(Smooth_scaled_cov,20,20))$vectors[,1:3]
  results[7] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
  results[8] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
  results[9] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
  results[10] <-min(sum((mcmc$eigvecLongmean[,3] - m1[,1])^2), sum((mcmc$eigvecLongmean[,3] + m1[,1])^2))
  results[11] <-min(sum((mcmc$eigvecLongmean[,2] - m1[,2])^2), sum((mcmc$eigvecLongmean[,2] + m1[,2])^2))
  results[12] <-min(sum((mcmc$eigvecLongmean[,1] - m1[,3])^2), sum((mcmc$eigvecLongmean[,1] + m1[,3])^2))

  m2 <- eigen(getMarginalFunc(Smooth_scaled_cov,20,20))$vectors[,1:3]

  results[13] <-min(Re(sum((resProduct$psi[,1] - m2[,1])^2)), Re(sum((resProduct$psi[,1] + m2[,1])^2)))
  results[14] <-min(Re(sum((resProduct$psi[,2] - m2[,2])^2)), Re(sum((resProduct$psi[,2] + m2[,2])^2)))
  results[15] <-min(Re(sum((resProduct$psi[,3] - m2[,3])^2)), Re(sum((resProduct$psi[,3] + m2[,3])^2)))
  results[16] <-min(Re(sum((mcmc$eigvecFuncmean[,3] - m2[,1])^2)), Re(sum((mcmc$eigvecFuncmean[,3] + m2[,1])^2)))
  results[17] <-min(Re(sum((mcmc$eigvecFuncmean[,2] - m2[,2])^2)), Re(sum((mcmc$eigvecFuncmean[,2] + m2[,2])^2)))
  results[18] <-min(Re(sum((mcmc$eigvecFuncmean[,1] - m2[,3])^2)), Re(sum((mcmc$eigvecFuncmean[,1] + m2[,3])^2)))

  results[19] <- sum((mu - resPACE$mu)^2)/(length(s) * length(t))
  results[20] <- sum((mu - as.numeric(mcmc$postmean))^2)/(length(s) * length(t))
  save(results, file= paste0("n",n,"_20long20func",index,".RData"))
  gc()
  results
})[3]
stopCluster(cl)

round(mean(n30_50[,20]),3)
short_results <- matrix(0,nrow=10,ncol=20)
for(i in c(1:10)){
  load(paste0("n30_ccc",i,".RData"))
  short_results[i,] <- results
}
short_results<-short_results[c(1:10),]
n30_results <- matrix(0,nrow=250,ncol = 20)
for(i in 1:250){
  load(paste0("n30_",i,".RData"))
  n30_results[i,] <- results
}
n30_results <- matrix(as.numeric(n30_results),nrow=250,ncol=20)

round(mean(n30_results[,5]),3)


n60_results <- matrix(0,nrow=250,ncol = 20)
for(i in c(1:19,21,23:30,32:42,44:209,211:250)){
  load(paste0("n60_",i,".RData"))
  n60_results[i,] <- results
}
n60_results <- matrix(as.numeric(n60_results),nrow=250,ncol=20)
rows <- apply(n60_results,1,function(i) all(i==0))
n60_results<-n60_results[!rows,]
round(mean(n60_results[,5]),3)

n90_results <- matrix(0,nrow=250,ncol = 20)
for(i in c(1:48,50:67,69:87,89:91,93:99,101:118,120:136,138,140:165,
           167:171,173,175:180,182:185,187:199,201,203:207,209:213,215:220,222:250)){
  load(paste0("n90_",i,".RData"))
  n90_results[i,] <- results
}
n90_results <- matrix(as.numeric(n90_results),nrow=250,ncol=20)
rows <- apply(n90_results,1,function(i) all(i==0))
n90_results<-n90_results[!rows,]

round(mean(n90_results[,5]),3)


mlongb <- getMarginalLong(mcmc$postcov,20,20)
mfuncb <- getMarginalFunc(mcmc$postcov,20,20)
mlongt <- getMarginalLong(Smooth_scaled_cov,20,20)
mfunct <- getMarginalFunc(Smooth_scaled_cov,20,20)
mlongf <- getMarginalLong(resProductCov,20,20)
mfuncf <- getMarginalFunc(resProductCov,20,20)
mlonge <- getMarginalLong(EmpCov,20,20)
mfunce <- getMarginalFunc(EmpCov,20,20)
image(mlongt[,20:1],col = heat.colors(100))
image(mlongb[,20:1],col = heat.colors(100))
image(mlongf[,20:1],col = heat.colors(100))
image(mlonge[,20:1],col = heat.colors(100))
plot(eigen(mlongf)$vectors[,1],type="l")
lines(-resProduct$phi[,1],col="red")
plot(eigen(mlongf)$vectors[,2],type="l")
lines(-resProduct$phi[,2],col="red")
plot(eigen(mlongf)$vectors[,3],type="l")
lines(resProduct$phi[,3],col="red")
plot(eigen(mfuncf)$vectors[,1],type="l")
lines(-resProduct$psi[,1],col="red")
plot(eigen(mfuncf)$vectors[,2],type="l")
lines(-resProduct$psi[,2],col="red")
plot(eigen(mfuncf)$vectors[,3],type="l")
lines(resProduct$psi[,3],col="red")

j <- 14
plot(mlongt[j,],type="l")
lines(mlongb[j,],col="blue")
lines(mlongf[j,],col="red")
lines(mlonge[j,],col="green")

j <- 199
plot(Smooth_scaled_cov[j,],type="l")
lines(mcmc$postcov[j,],col="blue")
lines(resProductCov[j,],col="red")
lines(EmpCov[j,],col="green")

j <- 2
plot(mfunct[j,],type="l")
lines(mfuncb[j,],col="blue")
lines(mfuncf[j,],col="red")
lines(mfunce[j,],col="green")

mylist <- list()
mylist[[1]] <- mfunct
mylist[[2]] <- mfuncb
mylist[[3]] <- mfuncf
mylist[[4]] <- mfunce
image(mfunct, col= heat.colors(100),zlim=c(min(unlist(mylist)),max(unlist(mylist))))
image(mfuncb, col= heat.colors(100),zlim=c(min(unlist(mylist)),max(unlist(mylist))))
image(mfuncf, col= heat.colors(100),zlim=c(min(unlist(mylist)),max(unlist(mylist))))
image(mfunce, col= heat.colors(100),zlim=c(min(unlist(mylist)),max(unlist(mylist))))

image(mlongt, col= heat.colors(100))
image(mlongb, col= heat.colors(100))
image(mlongf, col= heat.colors(100))
image(mlonge, col= heat.colors(100))




errorvar <- .01
Cov.Weak <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar*diag(length(s) * length(t))
iterations <- 500
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(6)
registerDoSNOW(cl)
system.time(n90_s.01_Sep<-foreach(index=1:iterations,.combine=rbind, .packages = c("MASS", "LFBayes", "fdapace"), .options.snow = opts)%dopar%{
  print(index)
  #x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
  x <- mvrnorm(30, mu  = as.vector(mu1), Sigma = Cov.Weak)
  sx <- sd(x)
  mx <- mean(x)
  x <- (x-mx)/sx
  Smooth_scaled_cov <- (Cov.Weak - errorvar * diag(length(s) * length(t))) / sx^2
  mu <- (mu1 - mu)/sx
  y <- lapply(1:n, function(i) x[i,])
  missing <- list()
  for(i in 1:n){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,n)
  dim(X) <- c(n,1)
  
  q <- 8
  mcmc <- mcmcWeak(y, missing, X, Bs1, Bt1, q, q, 25500, 1, 5000)
  
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
  results
})[3]
stopCluster(cl)
save(n90_s.01_Sep, file = "n90_s.01_Sep.RData")


errorvar <- .001
Cov.Weak <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar*diag(length(s) * length(t))
iterations <- 500
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(6)
registerDoSNOW(cl)
system.time(n90_s.001_Sep<-foreach(index=1:iterations,.combine=rbind, .packages = c("MASS", "LFBayes", "fdapace"), .options.snow = opts)%dopar%{
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
  
  q <- 8
  mcmc <- mcmcWeak(y, missing, X, Bs1, Bt1, q, q, 25000, 1, 5000)
  
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
  results
})[3]
stopCluster(cl)
save(n90_s.001_Sep, file = "n90_s.001_Sep.RData")

halft <- function(sig){
  nu <- 1
  A <- 10
  (1+1/nu*(sig/A)^2)^(-(nu+1)/2)
}

sample.x = runif(100000,0,20)
accept = c()
for(i in 1:length(sample.x)){
  U = runif(1, 0, 1)
  if(dunif(sample.x[i], 0, 20)*50*U <= halft(sample.x[i])) {
    
    accept[i] = 'Yes'
  }
  else if(dunif(sample.x[i],0,20)*50*U > halft(sample.x[i])) {
    accept[i] = 'No'
  }
}
T = data.frame(sample.x, accept = factor(accept, levels= c('Yes','No')))

We can plot the results along with the true distribution with the following code.

hist(T[,1][T$accept=='Yes'], breaks = seq(0,1,0.01), freq = FALSE, main = 'Histogram of X', xlab = 'X')
lines(sample.x, dbeta(x,6,3))
