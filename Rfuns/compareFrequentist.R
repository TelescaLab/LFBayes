library(splines)
library(MASS)
library(fdapace)
library(doParallel)
library(doSNOW)
library(coda)
library(MCMCglmm)
setwd("E:/Rcpp stuff/ChenCode")
source("MarginalFPCA.R")
source("ProductFPCA.R")

errorvar <- 1
t <- seq(from = 0, to = 1, length.out = 20)
s <- seq(from = 0, to = 1, length.out = 10)
n <- 50
tt <- list()
tt[[1]] <- 1:(length(t)*length(s))
tt <- rep(tt, n)
tt <- list()
tt[[1]] <- 1:(length(t)*length(s))
tt <- rep(tt, n)
p1 <- 10
p2 <- 5
q1 <- 3
q2 <- 3
Bt <- bs(t, df = p1, intercept = TRUE)
Bs <- bs(s, df = p2, intercept = TRUE)
Bt1 <- bs(t, df = 15, intercept = TRUE)
Bs1 <- bs(s, df = 10, intercept = TRUE)
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
mymode <- function(x){
  posterior.mode(mcmc(x))
}
covMode <- function(D){
  apply(D, c(1,2), mymode)
}
covMedian <- function(D){
  apply(D, c(1,2), median)
}
covMean <- function(D){
  apply(D, c(1,2), mean)
}
Lambda <- Loading.Matern(t, p1, q1, Bt)
#Lambda <- Loading.Brown.Bridge(t, p1, q1)
Gamma <- Loading.Brown.Motion(s, p2, q2)
Brown.Motion.Cov <- Bs%*%Gamma%*%t(Gamma)%*%t(Bs)
Brown.Bridge.Cov <- Bt%*%Lambda%*%t(Lambda)%*%t(Bt)
Cov.Strong <- kronecker(Brown.Motion.Cov, Brown.Bridge.Cov)
H <- diag(rgamma(q1*q2, shape = 10, rate = 10))
Cov.Weak <- 50*kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar*diag(length(s) * length(t))
SmoothedCov <- Cov.Weak - errorvar * diag(length(s) * length(t))
pc.j = NULL
pc.k = NULL
fpca.op1 = list(dataType = "Sparse", maxK = pc.j, FVEthreshold = .9999, nRegGrid = 20)
fpca.op2 = list(dataType = "Sparse", maxK = pc.k, FVEthreshold = .9999, nRegGrid = 10)



iterations <- 30
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(10)
registerDoSNOW(cl)
system.time(matrix3<-foreach(index=1:iterations,.combine=cbind, .packages = c("MASS", "LongFunc", "fdapace"), .options.snow = opts)%dopar%{
  print(index)
  x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
  y <- lapply(1:n, function(i) x[i,])
  missing <- list()
  for(i in 1:n){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,n)
  dim(X) <- c(n,1)

  q <- 6
  mcmc <- mcmcWeak(y, missing, X, Bs1, Bt1, q, q, 10000, 10)

  D <- getStatistics(kronecker(Bs1, Bt1), mcmc, 1, rep(0, length(t) * length(s)), 2000, 1)
  Vmean <- covMean(D$covC)
  Vmedian <- covMedian(D$covC)
  Vmode <- covMode(D$covC)
  resMarginal <- MarginalFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
  resProduct <- ProductFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resProductCov <- productcov(resProduct$eig, resProduct$scores)
  resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
  yt <- t(x)
  EmpCov <- 1 / (n-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
  EmpMean <- colMeans(x)
  results <- numeric(25)

  results[1] <- max(abs(eigen(Cov.Weak - resProductCov - resPACE$sigma2*diag(length(s) * length(t)))$values)) / (length(s) * length(t))
  results[2] <- max(abs(eigen(Cov.Weak - resMarginalCov - resPACE$sigma2*diag(length(s) * length(t)))$values)) / (length(s) * length(t))
  results[3] <- max(abs(eigen(Cov.Weak - resPACE$fittedCov)$values)) / (length(s) * length(t))
  results[4] <- max(abs(eigen(Cov.Weak - EmpCov)$values)) / (length(s) * length(t))
  results[5] <- max(abs(eigen(Cov.Weak - D$postcov)$values)) / (length(s) * length(t))
  results[6] <- max(abs(eigen(Cov.Weak - Vmedian)$values)) / (length(s) * length(t))
  results[7] <- max(abs(eigen(Cov.Weak - Vmode)$values)) / (length(s) * length(t))

  #results[6] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(resProductCov + resPACE$sigma2*diag(length(s) * length(t))))$values)) / (length(s) * length(t))
  #results[7] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(resMarginalCov + resPACE$sigma2*diag(length(s) * length(t))))$values)) / (length(s) * length(t))
  #results[8] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(resPACE$smoothedCov))$values)) / (length(s) * length(t))
  #results[9] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(EmpCov))$values)) / (length(s) * length(t))
  #results[10] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(D$postcov))$values)) / (length(s) * length(t))
  #results[11] <- sum((Cov.Weak - resProductCov)^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[12] <- sum((Cov.Weak - resMarginalCov)^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[13] <- sum((Cov.Weak - resPACE$smoothedCov)^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[14] <- sum((Cov.Weak - EmpCov)^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[15] <- sum((Cov.Weak - D$postcov)^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[16] <- sum((Cov.Weak - EmpCovBayes)^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[16] <- sum((cov2cor(Cov.Weak) - cov2cor(resProductCov))^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[17] <- sum((cov2cor(Cov.Weak) - cov2cor(resMarginalCov))^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[18] <- sum((cov2cor(Cov.Weak) - cov2cor(resPACE$smoothedCov))^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[19] <- sum((cov2cor(Cov.Weak) - cov2cor(EmpCov))^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[20] <- sum((cov2cor(Cov.Weak) - cov2cor(D$postcov))^2) / ((length(s) * length(t)) * (length(s) * length(t) + 1)/2)
  #results[21] <- sum(resProduct$mu^2)/(length(s) * length(t))
  #results[22] <- sum(resMarginal$mu^2)/(length(s) * length(t))
  results[23] <- sum(resPACE$mu^2)/(length(s) * length(t))
  results[24] <- sum(EmpMean^2)/(length(s) * length(t))
  results[25] <- sum(D$postmean^2)/(length(s) * length(t))
  results
})[3]
stopCluster(cl)

iter <- 9500
thetamat <- matrix(nrow = n, ncol = 15*10)
for(i in 1:15){
  for(j in 1:10){
    thetamat[,15 * (j - 1)+ i] <- mcmc$theta[[iter]][i,j,]
  }
}
q1 <- 20
q2 <- 10
etamat <- matrix(nrow = n, ncol =q1*q2)
for(i in 1:q1){
  for(j in 1:q2){
    etamat[,q1 * (j - 1) + i]<- mcmc$eta[[iter]][i,j,]
  }
}
A <- (kronecker(mcmc$Gamma[,,iter], mcmc$Lambda[,,iter])%*%solve(mcmc$HC[,,iter])%*%t(kronecker(mcmc$Gamma[,,iter],mcmc$Lambda[,,iter])) + diag(kronecker(1/mcmc$sigma2[,iter],1/mcmc$sigma1[,iter])))
B <- kronecker(mcmc$Gamma[,,iter], mcmc$Lambda[,,iter])%*%cov(etamat)%*%t(kronecker(mcmc$Gamma[,,iter], mcmc$Lambda[,,iter])) + diag(kronecker(1/mcmc$sigma2[,iter],1/mcmc$sigma1[,iter]))
C <- kronecker(Bs, Bt)%*%A %*%t(kronecker(Bs, Bt)) + 1/mcmc$varphi[iter]*diag(length(s)*length(t))
mylist2 <- list()
mylist2[[1]] <- cov(thetamat)
mylist2[[2]] <- A
mylist2[[3]] <- B
image(A[,(15*10):1], zlim = c(min(unlist(mylist2)),max(unlist(mylist2))), col = heat.colors(100))
image(cov(thetamat)[,(15*10):1], zlim = c(min(unlist(mylist2)),max(unlist(mylist2))), col = heat.colors(100))
image(B[,(15*10):1], zlim = c(min(unlist(mylist2)),max(unlist(mylist2))), col = heat.colors(100))
getCovWeak <- function(i){
  kronecker(Bs%*%mcmc$Gamma[,,i], Bt%*%mcmc$Lambda[,,i])%*%cov(etamat) %*% t(kronecker(Bs%*%mcmc$Gamma[,,i], Bt%*%mcmc$Lambda[,,i])) +
    kronecker(Bs%*%diag(1/mcmc$sigma2[,i])%*%t(Bs),Bt%*%diag(1/mcmc$sigma1[,i])%*%t(Bt)) + diag(1/mcmc$varphi[i], length(t)*length(s))
}
matrix<-foreach(index=1:2,.combine=cbind, .packages = c("MASS", "LongFunc", "fdapace"))%dopar%{
  x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
}
stopCluster(cl)
x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
y <- lapply(1:n, function(i) x[i,])
missing <- list()
for(i in 1:n){
  missing[[i]] <- numeric(0)
}
X <- rep(1,n)
dim(X) <- c(n,1)
q <- sample(3:6, 1)
mcmc <- mcmcWeak(y, missing, X, Bs, Bt, q, q, 10000, 5)
D <- getStatistics(kronecker(Bs, Bt), mcmc, 1, rep(0, length(t) * length(s)), 2000, 1)
resMarginal <- MarginalFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
resProduct <- ProductFPCA(x, n, length(s), length(t), fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
resProductCov <- productcov(resProduct$eig, resProduct$scores)
resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .99))
yt <- t(x)
EmpCov <- 1 / (n-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
EmpMean <- colMeans(x)
results <- numeric(15)
results[1] <- max(abs(eigen(Cov.Weak - resProductCov - resPACE$sigma2*diag(length(s) * length(t)))$values)) / (length(s) * length(t))
results[2] <- max(abs(eigen(Cov.Weak - resMarginalCov - resPACE$sigma2*diag(length(s) * length(t)))$values)) / (length(s) * length(t))
results[3] <- max(abs(eigen(Cov.Weak - resPACE$fittedCov)$values)) / (length(s) * length(t))
results[4] <- max(abs(eigen(Cov.Weak - EmpCov)$values)) / (length(s) * length(t))
results[5] <- max(abs(eigen(Cov.Weak - D$postcov)$values)) / (length(s) * length(t))
results[6] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(resProductCov + resPACE$sigma2*diag(length(s) * length(t))))$values)) / (length(s) * length(t))
results[7] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(resMarginalCov + resPACE$sigma2*diag(length(s) * length(t))))$values)) / (length(s) * length(t))
results[8] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(resPACE$fittedCov))$values)) / (length(s) * length(t))
results[9] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(EmpCov))$values)) / (length(s) * length(t))
results[10] <- max(abs(eigen(cov2cor(Cov.Weak) - cov2cor(D$postcov))$values)) / (length(s) * length(t))
results[11] <- sum(resProduct$mu^2)
results[12] <- sum(resMarginal$mu^2)
results[13] <- sum(resPACE$mu^2)
results[14] <- sum(EmpMean^2)
results[15] <- sum(D$postmean^2)

lines(resPACE$smoothedCov[1,], type = "l")
resPACE$smoothedCov[1:5,1:5]
resProductCov[1:5,1:5]
resMarginalCov[1:5,1:5]
Cov.Weak[1:5,1:5]

library(doSNOW)
  library(foreach)
  cl<-makeCluster(detectCores()-1) #change the 2 to your number of CPU cores
  registerDoSNOW(cl)
  matrix<-foreach(i=1:10,.combine=cbind) %dopar% {
    a
  }



  set.seed(1)
  pts <- seq(0, 1, by=0.05)
  sampWiener <- Wiener(n, pts)

  Ly <- lapply(1:20, function(i) sampWiener[i,])
  Lt <- lapply(1:20, function(i) pts)
  res <- FPCA(Ly, Lt,
              list(dataType='Dense', error=TRUE, kernel='epan', verbose=TRUE, useBinnedData = "OFF"))
  plot(res) # The design plot covers [0, 1] * [0, 1] well.
  CreateCovPlot(res, 'Fitted')
plot(res$mu, type = "l")
plot(colMeans(sampWiener), type = "l")
lines(res$mu, col = "red")
plot(res$mu, type = "l")
length(res$mu)
dim(res$fittedCov)
res$fittedCoresults[1:5,1:5]
cov(sampWiener)[1:5,1:5]



x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000
ptime <- system.time({
r <- foreach(i = 1:1000, .combine=cbind, .packages = "MASS") %dopar% {
        mvrnorm(1, mu = c(0,0), Sigma = diag(2))
         ind <- sample(100, 100, replace=TRUE)
         result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
         coefficients(result1)
       }
   })[3]
ptime

library(doParallel)
registerDoParallel(cores=11)
foreach(i=1:3) %dopar% sqrt(i)




library(doRNG)
set.seed(123)
a <- foreach(i=1:10,.combine=cbind, .packages = "MASS") %dopar% {
  x <- mvrnorm(2, mu  = rep(0, 2), Sigma = diag(2))
  x

}
a
set.seed(123)
b <- foreach(i=1:2,.combine=cbind) %dopar% {rnorm(5)}
identical(a,b)
