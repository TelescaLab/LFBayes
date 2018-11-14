library(fdapace)
t + 1
s[length(s)] - s[1]
c(sapply(seq(from = s[1], to = length(s), by = s[length(s)]-s[1])-s[1], function(i){t+i}))
tt <- seq(from = t[1], to = length(s), by = t[2] - t[1])
tt <- list()
tt[[1]] <- seq(from = 1, to = 400, length.out = 400)
tt <- rep(tt, n)
length(t)*length(s)


FPCADense <- FPCA(y, tt)
Cov.Weak[1:5,1:5]
Cov.Weak <- Cov.Weak + diag(1, length(t)*length(s))
FPCADense$smoothedCov[1:5,1:5]
x <- mvrnorm(n, mu  = rep(0, length(t)*length(s)), Sigma = Cov.Weak)
#EmpMean <- colMeans(x)
#lines(EmpMean, col = "blue")
#plot(EmpMean, type = "l")
y <- lapply(1:n, function(i) x[i,])
FPCADense <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .95))
FPCADense$fittedCov[1:5,1:5]
length(FPCADense$mu)
max(abs(eigen(FPCADense$fittedCov-Cov.Weak)$values)) / (length(t) * length(s))

max(abs(eigen(aa-Cov.Weak)$values)) / (length(t) * length(s))
max(abs(eigen(bb-Cov.Weak)$values)) / (length(t) * length(s))

max(abs(eigen(FPCADense$smoothedCov-Cov.Weak)$values)) / (length(t) * length(s))
max(abs(eigen(EmpCov-Cov.Weak)$values)) / (length(t) * length(s))

v = kronecker(res$psi[,1], res$phi[,1])
outer(v,v) * res$scores
length(res$phi[,1])

prodcov <- function(phi, psi, scores){
  n <- dim(scores)[1]
  j <- dim(phi)[2]
  k <- dim(psi)[2]
  v <- numeric(dim(psi)[1]*dim(phi)[1])

  mycov <- matrix(0, nrow = length(v), ncol = length(v))
  for(l in 1:k){
    for(i in 1:j){
      v = kronecker(phi[,i], psi[,l])
      mycov <- mycov + outer(v, v) * 1/n * as.numeric(scores[,i + j * (l - 1)]%*%scores[,i + j * (l - 1)])
    }
  }
  mycov
}

marginalcov <- function(phi, psi, scores){
  n <- dim(scores)[1]
  neig <- dim(res$phi)[2]
  v <- numeric(dim(psi)[1] * dim(phi)[1])
  mycov <- matrix(0, nrow = length(v), ncol = length(v))
  for(j in 1:neig){
    i <- ceiling(j / dim(res$psi)[2])
    v <- kronecker(phi[,j], psi[,i])
    mycov <- mycov + outer(v,v) * 1/n * as.numeric(scores[,j]%*%scores[,j])
  }
  mycov
}
res1 <- ProductFPCA(x, n, 20, 20, fpca.op1, fpca.op2, pc.j, pc.k)
res2 <- MarginalFPCA(x, n, 20, 20, fpca.op1, fpca.op2, pc.j, pc.k)
aa <- marginalcov(res2$phi, res2$psi, res2$scores)
bb <- prodcov(res1$phi, res1$psi, res1$scores)
aa[1:5,1:5]


Ly <- lapply(1:n, function(i) X.age.year[i,])
Lt <- lapply(1:n, function(i) 1:2464)
fertilitypace <- FPCA(Ly, Lt, list(useBinnedData = "OFF", FVEthreshold = .95))
