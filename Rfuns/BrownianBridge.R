library(splines)
library(MASS)
n <- 30
t <- seq(from = 0, to = 1, length.out = 20)
s <- seq(from = 0, to = 1, length.out = 10)
p1 <- 10
p2 <- 5
q1 <- 3
q2 <- 3
Bt <- bs(t, df = p1, intercept = TRUE)
Bs <- bs(s, df = p2, intercept = TRUE)
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
Matern.L<-Bs%*%Lambda%*%t(Lambda)%*%t(Bs)
MaternCov <- Matern.Cov(s)
Lambda <- Loading.Brown.Bridge(t, p1, q1)
Gamma <- Loading.Brown.Motion(s, p2, q2)
Brown.Motion.Cov <- Bs%*%Gamma%*%t(Gamma)%*%t(Bs)
Brown.Bridge.Cov <- Bt%*%Lambda%*%t(Lambda)%*%t(Bt)
Cov.Strong <- kronecker(Brown.Motion.Cov, Brown.Bridge.Cov)
IMSE <- matrix(0, nrow = 50, ncol = 3)
Cov.Norm <- matrix(0, nrow = 50, ncol = 3)
EmpCov.Norm <- matrix(0, nrow = 50, ncol = 3)
E <- matrix(0, nrow = 50, ncol = 3)
H <- diag(rgamma(q1*q2, shape = 10, rate = 10))
Cov.Weak <- 50*kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + .5*diag(length(s) * length(t))
factors <- c(3,6,9)
kseq <- c(35,50)
for(j in 2:2){
  for(k in 15:50){
    print(k)
    #H <- diag(rgamma(q1*q2, shape = 10, rate = 10))
    x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
    #EmpMean <- colMeans(x)
    #lines(EmpMean, col = "blue")
    #plot(EmpMean, type = "l")
    y <- lapply(1:n, function(i) x[i,])
    missing <- list()
    for(i in 1:n){
      missing[[i]] <- numeric(0)
    }
    X <- rep(1,n)
    dim(X) <- c(n,1)
    mcmc <- mcmcWeak(y, missing, X, Bs, Bt, 5, 5, 10000, 3)
    D <- getStatistics(kronecker(Bs, Bt), mcmc, 1, rep(0, length(t) * length(s)), 2000, 1)
    yt <- t(x)
    EmpCov <- 1 / n  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
    #IMSE[k,j] <- D$IMSE
    n90q369[[1]][k,j] <- D$IMSE
    n90q369[[2]][k,j] <- max(abs(eigen(D$postcov-Cov.Weak)$values)) / (length(t) * length(s))
    n90q369[[3]][k,j] <- max(abs(eigen(EmpCov-Cov.Weak)$values)) / (length(t) * length(s))
    n90q369[[4]][k,j] <- max(abs(eigen(EmpCov-D$postcov)$values)) / (length(t) * length(s))
    sum((D$postcov-Cov.Weak)^2)/400
    sum((EmpCov-Cov.Weak)^2)/400
    #Cov.Norm[k,j] <- max(abs(eigen(D$postcov-Cov.Weak)$values)) / (length(t) * length(s))
    #EmpCov.Norm[k,j] <- max(abs(eigen(EmpCov-Cov.Weak)$values)) / (length(t) * length(s))
    #E[k,j] <- max(abs(eigen(EmpCov-D$postcov)$values)) / (length(t) * length(s))
    #save(IMSE, file = "IMSE.RData")
    #save(Cov.Norm, file = "Cov.Norm.RData")
    #save(EmpCov.Norm, file = "EmpCov.Norm.RData")
    #save(E, file = "E.RData")
  }
}
mcmc <- mcmcWeak(y, missing, X, Bs, Bt, 6, 6, 3000)

temp <- list(IMSE = IMSE, Cov.Norm = Cov.Norm, EmpCov.Norm = EmpCov.Norm, E = E)
save(temp, file = "temp.RData")
H <- diag(rgamma(q1*q2, shape = 10, rate = 10))
Cov.Weak <- 50*kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + .25*diag(length(s) * length(t))
x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
#EmpMean <- colMeans(x)
#lines(EmpMean, col = "blue")
#plot(EmpMean, type = "l")
y <- lapply(1:n, function(i) x[i,])
plot(y[[1]][90:150], type = "l")
missing <- list()
for(i in 1:n){
  missing[[i]] <- numeric(0)
}
X <- rep(1,n)
dim(X) <- c(n,1)
mcmc <- mcmcWeak(y, missing, X, Bs, Bt, 7, 7, 50000)

getCovWeak <- function(i){
  kronecker(Bs1%*%mcmc$Gamma[,,i], Bt1%*%mcmc$Lambda[,,i])%*%solve(mcmc$HC[,,i]) %*% t(kronecker(Bs1%*%mcmc$Gamma[,,i], Bt1%*%mcmc$Lambda[,,i])) +
    kronecker(Bs1%*%diag(1/mcmc$sigma2[,i])%*%t(Bs1),Bt1%*%diag(1/mcmc$sigma1[,i])%*%t(Bt1)) + diag(1/mcmc$varphi[i], length(t)*length(s))
}

fromt <- 5000
iter <- 25000
j <- 0
#v <- numeric(length(seq(from = 1, to = 10000, by = 100)))
V <- array(dim = c(length(s) * length(t),length(s) * length(t),length(seq(from = fromt, to = iter, by = 10))))
S <- matrix(0, nrow = length(t) * length(s), ncol = length(t) * length(s))
for(i in seq(from = fromt, to = iter, by = 10)){
  j <- j + 1
  #v[j] <- getCovWeak(i)[50,50]
  V[,,j] <- getCovWeak(i)
  #print(j)
  S <- V[,,j] + S
}
S <- S / (length(seq(from = fromt, to = iter, by = 10)))
#myfun <- function(x){
#  posterior.mode(mcmc(x))
#}
#system.time(Vmode <- vapply(D$covC, myfun, c(1,2)))
#system.time(Vmode <- apply(D$covC, c(1,2),myfun))

#Vmedian <- apply(D$covC, c(1,2), median)
#Vmean <- apply(V, c(1,2), mean)
A<-sapply(1:2000, function(i){max(abs(eigen(V[,,i] - Cov.Weak)$values))})/200
#lines(A1, col ="red")
plot(A, type = "l")
lines(A2, col = "blue")
sapply(1:100, function(i){max(abs(eigen(V2[,,i] - Cov.Weak)$values))})
A3 <- sapply(1:100, function(i){sum((V[,,i]-Cov.Weak)^2)})
plot(V[52,52,], type = "l")
S <- S / length(seq(from = 5000, to = 7500, by = 10))
yt <- t(x)
EmpCov <- 1 / n  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
max(abs(eigen(Vmean - Cov.Weak)$values)) / 200
max(abs(eigen(Vmedian - Cov.Weak)$values)) / 200
max(abs(eigen(Vmode - Cov.Weak)$values)) / 200
max(abs(eigen(Cov.Weak - EmpCov)$values)) /200
max(abs(eigen(S - EmpCov)$values))
mcmcEigen <- eigenLF(Bs1, Bt1, mcmc, 4, 200)
k <- 120
quantile(V[k,k,], c(.025,.975))
EmpCov[k,k]
Cov.Weak[k,k]
x.axis <- t
y.axis <- s
z1 <- round(quantile(mcmcEigen$eigvalFunc[20,], c(0.025,.975)),2)
z2 <- round(quantile(mcmcEigen$eigvalFunc[19,], c(0.025,.975)),2)
z3 <- round(quantile(mcmcEigen$eigvalFunc[18,], c(0.025,.975)),2)
z4 <- round(quantile(mcmcEigen$eigvalLong[20,], c(0.025,.975)),2)
z5 <- round(quantile(mcmcEigen$eigvalLong[19,], c(0.025,.975)),2)
z6 <- round(quantile(mcmcEigen$eigvalLong[18,], c(0.025,.975)),2)
library(ggplot2)
library(gridExtra)
plot1 <- qplot(x.axis, mcmcEigen$eigvecFuncmean[,4], geom = "line", xlab = "", ylab = "Eigenfunction 1") + annotate("text", x = Inf, y = Inf, label = paste0("(",sprintf("%.2f", z1[1]),", ", sprintf("%.2f", z1[2]), ")"), vjust = 1, hjust = 1) + geom_ribbon(aes(ymin=mcmcEigen$eigvecFunclower[,4], ymax=mcmcEigen$eigvecFuncupper[,4]), alpha=0.2)
plot2 <- qplot(x.axis, mcmcEigen$eigvecFuncmean[,3], geom = "line", xlab = "", ylab = "Eigenfunction 2") + annotate("text", x = Inf, y = Inf, label = paste0("(",sprintf("%.2f", z2[1]),", ", sprintf("%.2f", z2[2]), ")"), vjust = 1, hjust = 1) + geom_ribbon(aes(ymin=mcmcEigen$eigvecFunclower[,3], ymax=mcmcEigen$eigvecFuncupper[,3]), alpha=0.2)
plot3 <- qplot(x.axis, -1*mcmcEigen$eigvecFuncmean[,2], geom = "line", xlab = "", ylab = "Eigenfunction 3") + annotate("text", x = Inf, y = Inf, label = paste0("(",sprintf("%.2f", z3[1]),", ", sprintf("%.2f", z3[2]), ")"), vjust = 1, hjust = 1) + geom_ribbon(aes(ymin=-1*mcmcEigen$eigvecFunclower[,2], ymax=-1*mcmcEigen$eigvecFuncupper[,2]), alpha=0.2)
grid.arrange(plot1, plot2, plot3, nrow = 1, ncol = 3)

plot4 <- qplot(y.axis, mcmcEigen$eigvecLongmean[,4], geom = "line", xlab = "", ylab = "Eigenfunction 1") + annotate("text", x = Inf, y = Inf, label = paste0("(",sprintf("%.2f", z4[1]),", ", sprintf("%.2f", z4[2]), ")"), vjust = 1, hjust = 1) + geom_ribbon(aes(ymin=mcmcEigen$eigvecLonglower[,4], ymax=mcmcEigen$eigvecLongupper[,4]), alpha=0.2)
plot5 <- qplot(y.axis, -1*mcmcEigen$eigvecLongmean[,3], geom = "line", xlab = "", ylab = "Eigenfunction 2") + annotate("text", x = Inf, y = Inf, label = paste0("(",sprintf("%.2f", z5[1]),", ", sprintf("%.2f", z5[2]), ")"), vjust = 1, hjust = 1) + geom_ribbon(aes(ymin=-1*mcmcEigen$eigvecLonglower[,3], ymax=-1*mcmcEigen$eigvecLongupper[,3]), alpha=0.2)
plot6 <- qplot(y.axis, -1*mcmcEigen$eigvecLongmean[,2], geom = "line", xlab = "", ylab = "Eigenfunction 3") + annotate("text", x = Inf, y = Inf, label = paste0("(",sprintf("%.2f", z6[1]),", ", sprintf("%.2f", z6[2]), ")"), vjust = 1, hjust = 1) + geom_ribbon(aes(ymin=-1*mcmcEigen$eigvecLonglower[,2], ymax=-1*mcmcEigen$eigvecLongupper[,2]), alpha=0.2)
grid.arrange(plot4, plot5, plot6, nrow = 1, ncol = 3)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, nrow = 2, ncol = 3)
grid.arrange(arrangeGrob(plot1, plot2, plot3, top = "Brownian Bridge", nrow = 1), arrangeGrob(plot4, plot5, plot6, top = "Brownian Motion", nrow = 1), nrow = 2)
EmpMean <- colMeans(x)
plot(EmpMean, type = "l")
lines(mcmcEigen$postmean, col = "blue")
plot(eigen(M)$vectors[,3], type = "l", col = "blue")
for(j in seq(from = 1, to = 3900, length.out = 500)){
  lines(mcmcEigen$eigvecLong[,2,j])
}
plot(S[400,], type = "l")
lines(Cov.Strong[400,], col = "blue")
mcmc <- mcmcWeak()
matplot(t(y), type = "l")
Loading.Brown.Motion(s, 12, 6)
plot(sqrt(eigval) * psi, type = "l", ylim = c(-1,1))
lines(B%*%getCol(B, t, 1), col = "red")
lines(B%*%getCol(B, t, 2), col = "blue")
lines(B%*%getCol(B, t, 3), col = "green")
lines(B%*%getCol(B, t, 4), col = "purple")
lines(B%*%getCol(B, t, 5), col = "yellow")
Lambda <- Loading.Brown.Bridge(t, 10, 5)
B%*%Lambda
dim(Lambda)
dim(B)


brown.bridge.cov <- B%*%Lambda%*%t(Lambda)%*%t(B)
dim(brown.bridge.cov)
psi1 <- sqrt(2)*sin(1*pi*t)
psi2 <- sqrt(2)*sin(2*pi*t)
psi3 <- sqrt(2)*sin(3*pi*t)
psi4 <- sqrt(2)*sin(4*pi*t)
psi5 <- sqrt(2)*sin(5*pi*t)
eigval1 <- 1/(1^2*pi^2)
eigval2 <- 1/(2^2*pi^2)
eigval3 <- 1/(3^2*pi^2)
eigval4 <- 1/(4^2*pi^2)
eigval5 <- 1/(5^2*pi^2)
truth <- eigval1 * psi1 %*%t(psi1) + eigval2 * psi2 %*% t(psi2) + eigval3 * psi3 %*% t(psi3) + eigval4 * psi4 %*% t(psi4) + eigval5 * psi5%*%t(psi5)
truth[1:5,1:5]
brown.bridge.cov[1:5,1:5]
plot(brown.bridge.cov[4,], type = "l")
lines(truth[4,], col = "blue")
max(abs(eigen(brown.bridge.cov-truth)$values))
t(psi1)%*%psi1

y <- mvrnorm(1, mu = rep(0,40), Sigma = brown.bridge.cov)
plot(y, type = "l")



## sapply(*, "array") -- artificial example
(v <- structure(10*(5:8), names = LETTERS[1:4]))
f2 <- function(x, y) outer(rep(x, length.out = 3), y)
(a2 <- sapply(v, f2, y = 2*(1:5), simplify = "array"))
a.2 <- vapply(v, f2, outer(1:3, 1:5), y = 2*(1:5))

V2small <- V2[1:5,1:5,1:50]
Vmode <- apply(V2, c(1,2), myfun)
posterior.mode(mcmc(V2[1,1,]))
myfun <- function(x){
  posterior.mode(mcmc(x))
}


