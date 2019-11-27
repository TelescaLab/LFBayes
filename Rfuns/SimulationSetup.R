### Source me

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
Loading.Brown.Motion <- function(t, p, k){
  B <- bs(t, df = p, intercept = TRUE)
  Loading <- matrix(nrow = p, ncol = k)
  for(i in 1:k){
    eigval <- 1/((k - 1/2)^2*pi^2)
    psi <- sqrt(2)*sin((k-1/2)*pi*t)
    Loading[,i] <- sqrt(eigval)*solve(t(B)%*%B)%*%t(B)%*%psi
  }
  return(Loading)
}
Matern.Cov <- function(s){
  rho = 0.5
  Matern.Cov <- matrix(nrow = length(s), ncol = length(s))
  for(i in 1:length(s)){
    for(j in 1:length(s)){
      d <- abs(s[i] - s[j])
      Matern.Cov[i, j] <- (1 + sqrt(5) * d / rho + 5 * d^2 / (3*rho^2)) * exp(-sqrt(5) * d / rho)
    }
  }
  Matern.Cov
}

Loading.Matern <- function(s, p, k, B){
  rho = 0.5
  Matern.Cov <- matrix(nrow = length(s), ncol = length(s))
  for(i in 1:length(s)){
    for(j in 1:length(s)){
      d <- abs(s[i] - s[j])
      Matern.Cov[i, j] <- (1 + sqrt(5) * d / rho + 5 * d^2 / (3*rho^2)) * exp(-sqrt(5) * d / rho)
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


Loading.CosCov <- function(s, p, q, B){
  alpha <- 2
  CosCov <- matrix(0,nrow=length(s),ncol=length(s))
  for(k in 1:50){
    CosCov <- k^(-2*alpha)*outer(cos(k*pi*s),cos(k*pi*s)) + CosCov
  }
  Loading <- matrix(nrow = p, ncol = q)
  evec <- eigen(CosCov)$vectors
  eval <- eigen(CosCov)$values
  for(i in 1:q){
    print(i)
    Loading[,i] <- sqrt(eval[i]) * solve(t(B)%*%B)%*%t(B)%*%evec[,i]
  }
  Loading
}
GenerateNonsep <- function(s,t){
  Cov <- matrix(0, nrow = length(t)*length(s), ncol = length(t)*length(s))
  S <- length(s)
  T <- length(t)
  for(i in 1:S){
    for(ii in 1:S){
      for(j in 1:T){
        for(jj in 1:T){
          Cov[(i-1)*T + j, (ii - 1)*T + jj] <- 1/((t[j]-t[jj])^2+1) * exp(-((s[i]-s[ii])^2)/((t[j]-t[jj])^2+1))
        }
      }
    }
  }
  Cov
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
GenerateMu2 <- function(s,t){
  mu <- matrix(0, nrow = length(t),ncol=length(s))
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      mu[i,j] <- 5*sqrt(1-(s[j]-.5)^2-(t[i]-.5)^2)
    }
  }
  c(mu)
}
GenerateMu3 <- function(s,t){
  mu <- matrix(0, nrow = length(t),ncol=length(s))
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      mu[i,j] <- sqrt(1+sin(pi*s[j]) + cos(pi*t[i]))
    }
  }
  c(mu)
}

### Run for automated array job using multicore
Run_Arrayjob <- function(numsimulation, Cov, Mu, Errorvar, N, SS, TT, Bs1, Bt1, q1s, q2s, iter,
                          thin, info_thin, burnin, nchain, neig, fpca.op1, fpca.op2, pc.j, pc.k,
                          str_filename, str_directory){
  system.time(simulation <-foreach(index=1:numsimulation,.combine=rbind,
                                   .packages = c("MASS", "LFBayes", "fdapace"),
                                   .export = c("MarginalFPCA", "ProductFPCA", "marginalcov", "productcov", "mFPCA.lm"),
                                   .options.snow = opts)%dopar%{
    print(index)
    #x <- mvrnorm(n, mu  = rep(0, TT*SS), Sigma = Cov.Weak)
    x <- mvrnorm(N, mu  = as.vector(Mu), Sigma = Cov)
    sx <- sd(x)
    mx <- mean(x)
    x <- (x-mx)/sx
    Smooth_scaled_cov <- (Cov - Errorvar * diag(SS * TT)) / sx^2
    mu <- (Mu - mx)/sx
    Marg.Long <- getMarginalLong(Smooth_scaled_cov,SS,TT)
    Marg.Func <- getMarginalFunc(Smooth_scaled_cov,SS,TT)
    #scaled_cov <- Cov/sx^2
    y <- lapply(1:N, function(i) x[i,])
    missing <- list()
    for(i in 1:N){
      missing[[i]] <- numeric(0)
    }
    X <- rep(1,N)
    dim(X) <- c(N,1)
    tt <- list()
    tt[[1]] <- 1:(TT*SS)
    tt <- rep(tt, N)
    MCMC <- mcmcWeakChains(y, missing, X, Bs1, Bt1, q1s, q2s, iter, thin, burnin, nchain)
    MCMC_eigen <- eigenLFChains(Bs1, Bt1, MCMC, neig, iter, burnin, nchain)
    infoD <- calculate_DIC(t(x), X, MCMC, Bs1, Bt1, burnin, info_thin)
    infoB <- calculate_BIC(t(x), X, MCMC, Bs1, Bt1, burnin, info_thin)
    
    resMarginal <- MarginalFPCA(x, N, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
    resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
    resProduct <- ProductFPCA(x, N, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
    resProductCov <- productcov(resProduct$eig, resProduct$scores)
    resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
    yt <- t(x)
    EmpCov <- 1 / (N-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
    EmpMean <- colMeans(x)
    results <- numeric(26)
    results[1] <- infoD$DIC
    results[2] <- infoB$BIC1
    results[3] <- infoB$BIC2
    results[4] <- sum((Smooth_scaled_cov - resProductCov)^2) / sum(Smooth_scaled_cov^2)
    results[5] <- sum((Smooth_scaled_cov - resMarginalCov)^2) / sum(Smooth_scaled_cov^2)
    results[6] <- sum((Smooth_scaled_cov - resPACE$smoothedCov)^2) / sum(Smooth_scaled_cov^2)
    results[7] <- sum((Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(SS * TT)))^2) / sum(Smooth_scaled_cov^2)
    results[8] <- sum((Smooth_scaled_cov - MCMC_eigen$postcov)^2) / sum(Smooth_scaled_cov^2)
    
    m1 <- eigen(Marg.Long)$vectors[,1:3] # Marginal longitudinal eigenvectors
    
    results[9] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
    results[10] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
    results[11] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
    results[12] <-min(sum((MCMC_eigen$eigvecLongmean[,3] - m1[,1])^2), sum((MCMC_eigen$eigvecLongmean[,3] + m1[,1])^2))
    results[13] <-min(sum((MCMC_eigen$eigvecLongmean[,2] - m1[,2])^2), sum((MCMC_eigen$eigvecLongmean[,2] + m1[,2])^2))
    results[14] <-min(sum((MCMC_eigen$eigvecLongmean[,1] - m1[,3])^2), sum((MCMC_eigen$eigvecLongmean[,1] + m1[,3])^2))
    
    m2 <- eigen(Marg.Func)$vectors[,1:3] # Marginal functional eigenvectors
    
    results[15] <-min(Re(sum((resProduct$psi[,1] - m2[,1])^2)), Re(sum((resProduct$psi[,1] + m2[,1])^2)))
    results[16] <-min(Re(sum((resProduct$psi[,2] - m2[,2])^2)), Re(sum((resProduct$psi[,2] + m2[,2])^2)))
    results[17] <-min(Re(sum((resProduct$psi[,3] - m2[,3])^2)), Re(sum((resProduct$psi[,3] + m2[,3])^2)))
    results[18] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,3] - m2[,1])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,3] + m2[,1])^2)))
    results[19] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,2] - m2[,2])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,2] + m2[,2])^2)))
    results[20] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,1] - m2[,3])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,1] + m2[,3])^2)))
    
    # Compare marginal covariances
    
    results[21] <- sum((getMarginalLong(MCMC_eigen$postcov,SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
    results[22] <- sum((getMarginalLong(resProductCov, SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
    results[23] <- sum((getMarginalFunc(MCMC_eigen$postcov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
    results[24] <- sum((getMarginalFunc(resProductCov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
    
    results[25] <- sum((mu - EmpMean)^2)/ sum(mu^2)
    results[26] <- sum((mu - as.numeric(MCMC_eigen$postmean))^2) / sum(mu^2)
    #save(results, file= paste0("splinenum",splinenum,"_20long20splines",index,".RData"))
    gc()
    results
  })[3]
  stopCluster(cl)
  setwd(str_directory)
  save(simulation, file = str_filename)
}

Single_Job <- function(Cov, Mu, Errovar, N, SS, TT, Bs1, Bt1, q1s, q2s, iter, thin, info_thin, burnin, nchain, neig,fpca.op1, fpca.op2, pc.j, pc.k){
  print(dim(Cov))
  x <- mvrnorm(N, mu  = as.vector(Mu), Sigma = Cov)
  sx <- sd(x)
  mx <- mean(x)
  x <- (x-mx)/sx
  Smooth_scaled_cov <- (Cov - Errorvar * diag(SS * TT)) / sx^2
  mu <- (Mu - mx)/sx
  Marg.Long <- getMarginalLong(Smooth_scaled_cov,SS,TT)
  Marg.Func <- getMarginalFunc(Smooth_scaled_cov,SS,TT)
  y <- lapply(1:N, function(i) x[i,])
  missing <- list()
  for(i in 1:N){
    missing[[i]] <- numeric(0)
  }
  X <- rep(1,N)
  dim(X) <- c(N,1)
  tt <- list()
  tt[[1]] <- 1:(TT*SS)
  tt <- rep(tt, N)
  MCMC <- mcmcWeakChains(y, missing, X, Bs1, Bt1, q1s, q2s, iter, thin, burnin, nchain)
  MCMC_eigen <- eigenLFChains(Bs1, Bt1, MCMC, neig, iter, burnin, nchain)
  infoD <- calculate_DIC(t(x), X, MCMC, Bs1, Bt1, burnin, info_thin)
  infoB <- calculate_BIC(t(x), X, MCMC, Bs1, Bt1, burnin, info_thin)
  resMarginal <- MarginalFPCA(x, N, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
  resProduct <- ProductFPCA(x, N, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
  resProductCov <- productcov(resProduct$eig, resProduct$scores)
  resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
  yt <- t(x)
  EmpCov <- 1 / (N-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
  EmpMean <- colMeans(x)
  results <- numeric(26)
  results <- numeric(26)
  results[1] <- infoD$DIC
  results[2] <- infoB$BIC1
  results[3] <- infoB$BIC2
  results[4] <- sum((Smooth_scaled_cov - resProductCov)^2) / sum(Smooth_scaled_cov^2)
  results[5] <- sum((Smooth_scaled_cov - resMarginalCov)^2) / sum(Smooth_scaled_cov^2)
  results[6] <- sum((Smooth_scaled_cov - resPACE$smoothedCov)^2) / sum(Smooth_scaled_cov^2)
  results[7] <- sum((Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(SS * TT)))^2) / sum(Smooth_scaled_cov^2)
  results[8] <- sum((Smooth_scaled_cov - MCMC_eigen$postcov)^2) / sum(Smooth_scaled_cov^2)
  
  m1 <- eigen(Marg.Long)$vectors[,1:3] # Marginal longitudinal eigenvectors
  
  results[9] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
  results[10] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
  results[11] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
  results[12] <-min(sum((MCMC_eigen$eigvecLongmean[,3] - m1[,1])^2), sum((MCMC_eigen$eigvecLongmean[,3] + m1[,1])^2))
  results[13] <-min(sum((MCMC_eigen$eigvecLongmean[,2] - m1[,2])^2), sum((MCMC_eigen$eigvecLongmean[,2] + m1[,2])^2))
  results[14] <-min(sum((MCMC_eigen$eigvecLongmean[,1] - m1[,3])^2), sum((MCMC_eigen$eigvecLongmean[,1] + m1[,3])^2))
  
  m2 <- eigen(Marg.Func)$vectors[,1:3] # Marginal functional eigenvectors
  
  results[15] <-min(Re(sum((resProduct$psi[,1] - m2[,1])^2)), Re(sum((resProduct$psi[,1] + m2[,1])^2)))
  results[16] <-min(Re(sum((resProduct$psi[,2] - m2[,2])^2)), Re(sum((resProduct$psi[,2] + m2[,2])^2)))
  results[17] <-min(Re(sum((resProduct$psi[,3] - m2[,3])^2)), Re(sum((resProduct$psi[,3] + m2[,3])^2)))
  results[18] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,3] - m2[,1])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,3] + m2[,1])^2)))
  results[19] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,2] - m2[,2])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,2] + m2[,2])^2)))
  results[20] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,1] - m2[,3])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,1] + m2[,3])^2)))
  
  # Compare marginal covariances
  
  results[21] <- sum((getMarginalLong(MCMC_eigen$postcov,SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
  results[22] <- sum((getMarginalLong(resProductCov, SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
  results[23] <- sum((getMarginalFunc(MCMC_eigen$postcov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
  results[24] <- sum((getMarginalFunc(resProductCov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
  
  results[25] <- sum((mu - EmpMean)^2)/ sum(mu^2)
  results[26] <- sum((mu - as.numeric(MCMC_eigen$postmean))^2) / sum(mu^2)
  results
}
### Run for more simple output
Array_Job <- function(Cov, Mu, Errorvar, N, SS, TT, Bs1, Bt1, q1s, q2s, iter,
                          thin, info_thin, burnin, nchain, neig, fpca.op1, fpca.op2, pc.j, pc.k){
                                     x <- mvrnorm(N, mu  = as.vector(Mu), Sigma = Cov)
                                     sx <- sd(x)
                                     mx <- mean(x)
                                     x <- (x-mx)/sx
                                     Smooth_scaled_cov <- (Cov - Errorvar * diag(SS * TT)) / sx^2
                                     print(dim(Smooth_scaled_cov))
                                     print(Smooth_scaled_cov[1:5,1:5])
                                     mu <- (Mu - mx)/sx
                                     Marg.Long <- getMarginalLong(Smooth_scaled_cov,SS,TT)
                                     #Marg.Func <- getMarginalFunc(Smooth_scaled_cov,SS,TT)
                                     #scaled_cov <- Cov/sx^2
                                     y <- lapply(1:N, function(i) x[i,])
                                     missing <- list()
                                     for(i in 1:N){
                                       missing[[i]] <- numeric(0)
                                     }
                                     X <- rep(1,N)
                                     dim(X) <- c(N,1)
                                     tt <- list()
                                     tt[[1]] <- 1:(TT*SS)
                                     tt <- rep(tt, N)
                                     MCMC <- mcmcWeakChains(y, missing, X, Bs1, Bt1, q1s, q2s, iter, thin, burnin, nchain)
                                     MCMC_eigen <- eigenLFChains(Bs1, Bt1, MCMC, neig, iter, burnin, nchain)
                                     infoD <- calculate_DIC(t(x), X, MCMC, Bs1, Bt1, burnin, info_thin)
                                     infoB <- calculate_BIC(t(x), X, MCMC, Bs1, Bt1, burnin, info_thin)
                                     
                                     resMarginal <- MarginalFPCA(x, N, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
                                     resMarginalCov <- marginalcov(resMarginal$eig, resMarginal$scores)
                                     resProduct <- ProductFPCA(x, N, SS, TT, fpca.op1, fpca.op2, pc.j, rep(pc.k, pc.j))
                                     resProductCov <- productcov(resProduct$eig, resProduct$scores)
                                     resPACE <- FPCA(y, tt, list(useBinnedData = "OFF", FVEthreshold = .9999, dataType = "Dense"))
                                     yt <- t(x)
                                     EmpCov <- 1 / (N-1)  * (yt - rowMeans(yt))%*%t(yt - rowMeans(yt))
                                     EmpMean <- colMeans(x)
                                     results <- numeric(26)
                                     results[1] <- infoD$DIC
                                     results[2] <- infoB$BIC1
                                     results[3] <- infoB$BIC2
                                     results[4] <- sum((Smooth_scaled_cov - resProductCov)^2) / sum(Smooth_scaled_cov^2)
                                     results[5] <- sum((Smooth_scaled_cov - resMarginalCov)^2) / sum(Smooth_scaled_cov^2)
                                     results[6] <- sum((Smooth_scaled_cov - resPACE$smoothedCov)^2) / sum(Smooth_scaled_cov^2)
                                     results[7] <- sum((Smooth_scaled_cov - (EmpCov - resPACE$sigma2 * diag(SS * TT)))^2) / sum(Smooth_scaled_cov^2)
                                     results[8] <- sum((Smooth_scaled_cov - MCMC_eigen$postcov)^2) / sum(Smooth_scaled_cov^2)
                                     
                                     m1 <- eigen(Marg.Long)$vectors[,1:3] # Marginal longitudinal eigenvectors
                                     
                                     results[9] <- min(sum((resProduct$phi[,1] - m1[,1])^2), sum((resProduct$phi[,1] + m1[,1])^2))
                                     results[10] <-min(sum((resProduct$phi[,2] - m1[,2])^2), sum((resProduct$phi[,2] + m1[,2])^2))
                                     results[11] <-min(sum((resProduct$phi[,3] - m1[,3])^2), sum((resProduct$phi[,3] + m1[,3])^2))
                                     results[12] <-min(sum((MCMC_eigen$eigvecLongmean[,3] - m1[,1])^2), sum((MCMC_eigen$eigvecLongmean[,3] + m1[,1])^2))
                                     results[13] <-min(sum((MCMC_eigen$eigvecLongmean[,2] - m1[,2])^2), sum((MCMC_eigen$eigvecLongmean[,2] + m1[,2])^2))
                                     results[14] <-min(sum((MCMC_eigen$eigvecLongmean[,1] - m1[,3])^2), sum((MCMC_eigen$eigvecLongmean[,1] + m1[,3])^2))
                                     
                                     m2 <- eigen(Marg.Func)$vectors[,1:3] # Marginal functional eigenvectors
                                     
                                     results[15] <-min(Re(sum((resProduct$psi[,1] - m2[,1])^2)), Re(sum((resProduct$psi[,1] + m2[,1])^2)))
                                     results[16] <-min(Re(sum((resProduct$psi[,2] - m2[,2])^2)), Re(sum((resProduct$psi[,2] + m2[,2])^2)))
                                     results[17] <-min(Re(sum((resProduct$psi[,3] - m2[,3])^2)), Re(sum((resProduct$psi[,3] + m2[,3])^2)))
                                     results[18] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,3] - m2[,1])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,3] + m2[,1])^2)))
                                     results[19] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,2] - m2[,2])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,2] + m2[,2])^2)))
                                     results[20] <-min(Re(sum((MCMC_eigen$eigvecFuncmean[,1] - m2[,3])^2)), Re(sum((MCMC_eigen$eigvecFuncmean[,1] + m2[,3])^2)))
                                     
                                     # Compare marginal covariances
                                     
                                     results[21] <- sum((getMarginalLong(MCMC_eigen$postcov,SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
                                     results[22] <- sum((getMarginalLong(resProductCov, SS, TT) - Marg.Long)^2) / sum(Marg.Long^2)
                                     results[23] <- sum((getMarginalFunc(MCMC_eigen$postcov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
                                     results[24] <- sum((getMarginalFunc(resProductCov,SS,TT) - Marg.Func)^2) / sum(Marg.Func^2)
                                     
                                     results[25] <- sum((mu - EmpMean)^2)/ sum(mu^2)
                                     results[26] <- sum((mu - as.numeric(MCMC_eigen$postmean))^2) / sum(mu^2)
                                     #save(results, file= paste0("splinenum",splinenum,"_20long20splines",index,".RData"))
                                     gc()
                                     results
}

