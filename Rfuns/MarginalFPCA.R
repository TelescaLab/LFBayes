
# The main function for Marginal Functional Principal Component Analysis
# Reference: Modeling Function-Valued Stochastic Processes, With Applications to Fertility Dynamics
#            by Kehui Chen, Pedro Delicado and Hans-Georg M\"{u}ller
# The code works for dense functional data, missing values are allowed.
# Written by kehui chen 10/09/2015, based on the original code by Pedro Delicado

#### Usage  MarginalFPCA(X.age.year, n, num.years, num.ages, fpca.op1, fpca.op2, pc.j, pc.k)
#### Input:
     # X.age.year: your data, a matrix with dimension n by (num.years*num.ages)
		 # n: sample size
		 # num.years: dimension 1
		 # num.ages: dimension 2
		 # fpca.op1: optns for fpca in age direction, check the help function of FPCA in tPACE for details
		 # fpca.op2: optns for fpca in year direction
		 # pc.j: number of components for first step, chosen by FVE if empty
		 # pc.k: a vector of 1*pc.j, number of components for the second step, chosen by FVE if empty

#### output a list of the following values
     # Xest: The estimated X.age.year from marginal FPCA
		 # mu: The estimated mean function
		 # scores: loadings \chi_{jk}
		 # res.psi: the FPCA output for psi
		 # res.phi: a list contains the FPCA output for second steps
		 # eig: the estimated product functions
		 # pc.j: number of components for first step
		 # psi: a matrix of num.ages * pc.j, contains eigenfunctions from step 1
		 # pc.k: a vector of 1*pc.j, number of components for the second step,
		 # phi: a matrix of num.years*(pc.j*pc.k), contains all eigenfunctions phi from step 2,
		 #       combining eigenfunctions for each $\xi_j$ processes.
		 # FVE: the variance explained by each product function
		 # VarOrdered: Varance explaned by each term. The terms are ordered by
		 #         var(\chi_{jk}). One can select the best model by truncating at a desired level of FVE, and         use names(VarOrdered) to see the corresponding model terms.

#### example code
 #      source('MarginalFPCA.R')
 #      load('X_age_year.RData')
 # 			n = dim(X.age.year)[1]
 # 			num.years = 56
 # 			num.ages = 44
 # 			pc.j = 3
 # 			pc.k = 3
 # 			res <- MarginalFPCA(X.age.year, n, num.years, num.ages, fpca.op1, fpca.op2, pc.j, rep(pc.k,pc.j))
 #


source("mFPCA.lm.R")  # projection
library(fdapace)  # PACE functions

MarginalFPCA <- function(X.age.year, n, num.years, num.ages, fpca.op1, fpca.op2, pc.j = NULL, pc.k = NULL){
# step 1: centering
      mean.X = apply(X.age.year, 2, mean, na.rm = TRUE)
      mean.Xmat = matrix(rep(mean.X,n), nrow = n, ncol = num.years*num.ages, byrow = TRUE)
      X.age.year.c = X.age.year-mean.Xmat
      X.age.c = matrix(t(X.age.year.c), nrow = n*num.years, ncol = num.ages, byrow = TRUE)

# step 2: computing marginal eigenfunctions psi

			# use PACE functions so that it can handle missing values

			if (sum(is.na(X.age.c)) == 0){
				tList <- lapply(1:(n*num.years), function(x) 1:num.ages)
			  yList <- lapply(1:(n*num.years), function(x) X.age.c[x,])
		  } else {
				fpca.op1$dataType = 'DenseWithMV'
				tList <- vector('list', length  = n*num.years)
				yList <- vector('list', length = n*num.years)
				tvec <- 1:num.ages
		  	for (i in 1:(n*num.years)){
		  		ind <- which(!is.na(X.age.c[i,]))
					yList[[i]] <- X.age.c[i,ind]
					tList[[i]] <- tvec[ind]
		  	}
		  }
			res.psi <- FPCA(yList, tList, optns = fpca.op1)
			psi <- res.psi$phi
			if (is.null(pc.j)){
				pc.j <- dim(psi)[2]
	    } else {
	    	psi <- psi[,1:pc.j]
	    }

			xiEst <- res.psi$xiEst


 # step 3
         res.phi <- vector('list', length = pc.j)
				 pc.k.2 <- numeric(pc.j)
      for (j in 1:pc.j){
         fpc <- matrix(xiEst[,j],nrow=n,ncol=num.years,byrow=TRUE)
	 			 tList <- lapply(1:n, function(x) 1:num.years)
	 			 yList <- lapply(1:n, function(x) fpc[x,])
	 			 res <- FPCA(yList, tList, optns = fpca.op2)
				 res.phi[[j]] <- res
	 			 eigenfunctions <- res$phi
				 if (is.null(pc.k)){
				 	 pc.k.2[j] <- dim(eigenfunctions)[2]
				 } else {
				 	 pc.k.2[j] <- pc.k[j]
					 eigenfunctions <- eigenfunctions[,1:pc.k[j]]
				 }

         if (j == 1){
             phi <- eigenfunctions
         } else {
             phi <- cbind(phi,eigenfunctions)
         }
      }

			pc.k <- pc.k.2
      eig.Age.Year <- matrix(0, nrow= num.years*num.ages, ncol= sum(pc.k))

      i = 0
			id <- 1:pc.j
      for (j in 1:pc.j){
         for (k in 1:pc.k[j]){
            i = i +1
						#eig.Age.Year[,i] <- as.numeric(outer(psi[,j],phi[,(j-1)*pc.k[j]+k],"*"))
            eig.Age.Year[,i] <- as.numeric(outer(psi[,j],phi[,sum(pc.k[id<j])+k],"*"))

         }
      }

			lmfull.simu = mFPCA.lm(t(X.age.year.c),eig.Age.Year, ploting=F)
			scores <- lmfull.simu$coefs
      Xest = mean.Xmat+scores%*%t(eig.Age.Year)

			df.eig.Age.Year <- as.data.frame(eig.Age.Year)
     	r <- dim(df.eig.Age.Year)[2]
      mFPCA.Age.Year.pctvar <- numeric(r)
      for (i in 1:r){
        mFPCA.Age.Year.pctvar[i] <-
				     round(mFPCA.lm(t(X.age.year.c),df.eig.Age.Year,model=i,ploting=F)$perct.model,2)
      }

      best.model <- order(-mFPCA.Age.Year.pctvar)
      VarOrdered <- mFPCA.Age.Year.pctvar[best.model]
			namesV = vector('list', length = dim(eig.Age.Year)[2])
			i = 0
			for (j in 1:pc.j){
				for (k in 1:pc.k[j]){
					i = i+1
					namesV[i] = as.character(paste('Psi',j, 'Phi',j,k, sep = ""))
				}
			}
		  names(VarOrdered) = namesV[best.model]
			return(list(Xest = Xest, mu = mean.X, scores = scores, res.psi = res.psi, res.phi = res.phi, eig = eig.Age.Year,
				          pc.j = pc.j, psi = psi, pc.k = pc.k, phi = phi, FVE = mFPCA.Age.Year.pctvar, VarOrdered = VarOrdered))
}



