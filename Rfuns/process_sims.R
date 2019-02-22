getwd()
results_30 <- matrix(nrow = 100,ncol = 20)
for(i in 1:100){
  load(paste0("n30_20long20func",i,".RData"))
  results_30[i,]<-results
}

results_60 <- matrix(nrow = 100,ncol = 20)
for(i in 1:100){
  load(paste0("n60_20long20func",i,".RData"))
  results_60[i,]<-results
}
results_90 <- matrix(nrow = 100,ncol = 20)
for(i in 1:100){
  load(paste0("n90_20long20func",i,".RData"))
  results_90[i,]<-results
}
results_120 <- matrix(nrow = 100,ncol = 20)
for(i in 1:100){
  load(paste0("n120_20long20func",i,".RData"))
  results_120[i,]<-results
}

round(quantile(sims60[,35],c(.5,.1,.9)),3)
round(quantile(sims60[,36],c(.5,.1,.9)),3)
round(quantile(sims60[,5],c(.5,.1,.9)),3)



round(mean(results_30[,1]),3)
round(mean(results_60[,1]),3)
round(mean(results_90[,1]),3)
round(mean(results_120[,1]),3)

round(mean(results_30[,3]),3)
round(mean(results_60[,3]),3)
round(mean(results_90[,3]),3)
round(mean(results_120[,3]),3)

round(mean(results_30[,5]),3)
round(mean(results_60[,5]),3)
round(mean(results_90[,5]),3)
round(mean(results_120[,5]),3)

round(mean(results_30[,7]),3)
round(mean(results_60[,7]),3)
round(mean(results_90[,7]),3)
round(mean(results_120[,7]),3)

round(mean(results_30[,10]),3)
round(mean(results_60[,10]),3)
round(mean(results_90[,10]),3)
round(mean(results_120[,10]),3)

round(mean(results_30[,8]),3)
round(mean(results_60[,8]),3)
round(mean(results_90[,8]),3)
round(mean(results_120[,8]),3)

round(mean(results_30[,11]),3)
round(mean(results_60[,11]),3)
round(mean(results_90[,11]),3)
round(mean(results_120[,11]),3)

round(mean(results_30[,13]),3)
round(mean(results_60[,13]),3)
round(mean(results_90[,13]),3)
round(mean(results_120[,13]),3)

round(mean(results_30[,16]),3)
round(mean(results_60[,16]),3)
round(mean(results_90[,16]),3)
round(mean(results_120[,16]),3)

round(mean(results_30[,14]),3)
round(mean(results_60[,14]),3)
round(mean(results_90[,14]),3)
round(mean(results_120[,14]),3)

round(mean(results_30[,17]),3)
round(mean(results_60[,17]),3)
round(mean(results_90[,17]),3)
round(mean(results_120[,17]),3)

round(mean(results_30[,19]),3)
round(mean(results_60[,19]),3)
round(mean(results_90[,19]),3)
round(mean(results_120[,19]),3)

round(mean(results_30[,20]),3)
round(mean(results_60[,20]),3)
round(mean(results_90[,20]),3)
round(mean(results_120[,20]),3)



tempsum <- t(mcmc2$Lambda%*%mcmc2$eta[[500]][,,1])
S1 <- rep(1, dim(mcmc2$Lambda)[1])
S2 <- rep(1, dim(mcmc2$Gamma)[1])
D1 <- rep(1, dim(mcmc2$Lambda)[2])
D2 <- rep(1, dim(mcmc2$Gamma)[2])
Sig <- matrix(1, nrow = dim(mcmc2$Lambda)[1],ncol=dim(mcmc2$Gamma)[1])
phi1 <- matrix(1, nrow = dim(mcmc$Lambda)[1],ncol = dim(mcmc$Lambda)[2])
C <- updateLambdaSmoothD(mcmc2$eta[[500]], mcmc2$Gamma[,,500], S1, S2, D1, mcmc2$theta[[500]])
D <- updateLambdaSmoothDSig(mcmc2$eta[[500]], mcmc2$Gamma[,,500], Sig, D1, mcmc2$theta[[500]])
C[1:5,1:5]
D[1:5,1:5]
set.seed(1)
A <- updateGammaSmoothD(mcmc2$eta[[500]], mcmc2$Lambda[,,500], S1, S2, D1, mcmc2$theta[[500]])
B <- updateGammaSmoothDSig(mcmc2$eta[[500]], mcmc2$Lambda[,,500], Sig, D1, mcmc2$theta[[500]])
A[1:5,1:5]
B[1:5,1:5]
dim(A)
dim(B)

mcmcASD <- list()
mcmcASD$Lambda[[1]] <- mcmc1ASD$Lambda[[1]]
mcmcASD$Gamma[[1]] <- mcmc1ASD$Gamma[[1]]
mcmcASD$H[[1]] <- mcmc1ASD$H[[1]]
mcmcASD$Sigma[[1]] <- mcmc1ASD$Sigma[[1]]
mcmcASD$Beta[[1]] <- mcmc1ASD$Beta[[1]]
mcmcASD$Lambda[[2]] <- mcmc2ASD$Lambda[[1]]
mcmcASD$Gamma[[2]] <- mcmc2ASD$Gamma[[1]]
mcmcASD$H[[2]] <- mcmc2ASD$H[[1]]
mcmcASD$Sigma[[2]] <- mcmc2ASD$Sigma[[1]]
mcmcASD$Beta[[2]] <- mcmc2ASD$Beta[[1]]
mcmcASD$Lambda[[3]] <- mcmc3ASD$Lambda[[1]]
mcmcASD$Gamma[[3]] <- mcmc3ASD$Gamma[[1]]
mcmcASD$H[[3]] <- mcmc3ASD$H[[1]]
mcmcASD$Sigma[[3]] <- mcmc3ASD$Sigma[[1]]
mcmcASD$Beta[[3]] <- mcmc3ASD$Beta[[1]]

mcmcASD$Theta[[1,1]] <- mcmc1ASD$Theta[[1,1]]
