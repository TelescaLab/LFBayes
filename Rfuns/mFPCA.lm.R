# mFPCA.lm Linear model fitting from the output of marginal FPCA
#
# Input:
#       big.X Functional data matrix, dim = (n.t*n.s) x n, n=no.of funct and stional data, n.t, n.s number of values of args. 
#   eig.funct Matrix of eigen functions, dim = (n.t*n.s) x r, 
#       model Subset of 1:r, indicating what eig.funct are used in the lm model
#
mFPCA.lm <- function(big.X,eig.funct,model=1:dim(eig.funct)[2],label=1:dim(big.X)[2],printing=F,ploting=printing){
  n <- dim(big.X)[2]
  r <- dim(eig.funct)[2]
  R2 <- numeric(n)
  SSR <- numeric(n)
  SST <- numeric(n)
  coefs <- matrix(0,nrow=n,ncol=length(model))
  
  mu <- apply(big.X,1,mean) 
  eig.f <- as.matrix(eig.funct[,model])
  model.char <- names(as.data.frame(eig.funct))[model]
  
  for (i in 1:n){
    Xc.i <- big.X[,i]-mu
    lm.i <- lm(Xc.i ~ eig.f - 1)
    summ.lm.i <- summary(lm.i)
    R2[i] <- summ.lm.i$r.squared
    SSR[i] <- sum(summ.lm.i$residuals^2)
    SST[i] <- SSR[i]/(1-R2[i])
    coefs[i,] <- lm.i$coefficients
    if(printing){
      print(label[i])
      print((SST[i]-sum(Xc.i^2))/sum(Xc.i^2))
      print(summ.lm.i)
    } 
  }
  perct.model <- (1-sum(SSR)/sum(SST))*100
  if (ploting){
    op<-par(mfrow=c(1,1))
    boxplot(R2)
    bp.out <- which( R2 %in% boxplot.stats(R2, do.out = TRUE)$out )
    if (length(bp.out>0)) {
      text(1,R2[bp.out],label[bp.out],pos=4)
    }
  }
  
  #cat("model: ",model.char, "\n")
  #cat("Pctg.var.expalined by the model: ", round(perct.model,2), " %","\n")
  
  return(list(model=model,model.char=model.char,perct.model=perct.model,R2=R2,SSR=SSR,SST=SST,coefs=coefs))
}
