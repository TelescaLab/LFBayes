
# The main function for Product Functional Principal Component Analysis
# Reference: Modeling Function-Valued Stochastic Processes, With Applications to Fertility Dynamics
#            by Kehui Chen, Pedro Delicado and Hans-Georg M\"{u}ller
# The code works for dense functional data, missing values are allowed.
# Written by kehui chen 10/09/2015, based on the original code by Pedro Delicado

#### Usage  ProductFPCA(X.age.year, n, num.years, num.ages, fpca.op1, fpca.op2, pc.j, pc.k)
#### Input:
     # X.age.year: your data, a matrix with dimension n by (num.years*num.ages)
		 # n: sample size
		 # num.years: dimension 1
		 # num.ages: dimension 2
		 # fpca.op1: optns for fpca in age direction, check the help function of FPCA in tPACE for details
		 # fpca.op2: optns for fpca in year direction
		 # pc.j: number of components for age direction, chosen by FVE if empty
		 # pc.k: number of components for year direction, chosen by FVE if empty


#### output a list of the following values
     # Xest: The estimated X.age.year from product FPCA
		 # mu: The estimated mean function
		 # scores: loadings \chi_{jk}
		 # res.psi: the FPCA output for psi
		 # res.phi: the FPCA output for phi
		 # eig: the estimated product functions
		 # pc.j: number of components for first step
		 # psi: a matrix of num.ages * pc.j, contains eigenfunctions psi for ages
		 # pc.k: number of components for the second step
		 # phi: a matrix of num.years*pc.k, contains all eigenfunctions phi for years
		 # FVE: the variance explained by each product function
		 # VarOrdered: Varance explaned by each term. The terms are ordered by
		 #         var(\chi_{jk}). One can select the best model by truncating at a desired level of FVE, and         use names(VarOrdered) to see the corresponding model terms.

#### example code
 #      source('ProductFPCA.R')
 #      load('X_age_year.RData')
 # 			n = dim(X.age.year)[1]
 # 			num.years = 56
 # 			num.ages = 44
 # 			pc.j = 3
 # 			pc.k = 3
 # 			fpca.op1 = list(dataType = "Dense")
 #		  fpca.op2 = list(dataType = "Dense")
 #     res <- ProductFPCA(X.age.year, n, num.years, num.ages, fpca.op1, fpca.op2, pc.j, pc.k)
 #


source("mFPCA.lm.R")  # projection
library(fdapace)  # PACE functions

ProductFPCA <- function(X.age.year, n, num.years, num.ages, fpca.op1, fpca.op2, pc.j = NULL, pc.k = NULL){
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
				psi <- psi[, 1:pc.j]
			}


 # step 3: computing marginal eigenfunctions phi
      X.year.c = matrix(X.age.year.c,  nrow = n*num.ages, ncol = num.years)
			if (sum(is.na(X.year.c)) == 0){
				tList <- lapply(1:(n*num.ages), function(x) 1:num.years)
			  yList <- lapply(1:(n*num.ages), function(x) X.year.c[x,])
		  } else {
				fpca.op2$dataType = 'DenseWithMV'
				tList <- vector('list', length  = n*num.ages)
				yList <- vector('list', length = n*num.ages)
				tvec <- 1:num.years
		  	for (i in 1:(n*num.ages)){
		  		ind <- which(!is.na(X.year.c[i,]))
					yList[[i]] <- X.year.c[i,ind]
					tList[[i]] <- tvec[ind]
		  	}
		  }

			res.phi <- FPCA(yList, tList, optns = fpca.op2)
			phi <- res.phi$phi
			if (is.null(pc.k)){
				pc.k <- dim(phi)[2]
			} else {
				phi <- phi[, 1:pc.k]
			}



      eig.Age.Year <- matrix(0, nrow= num.years*num.ages, ncol= pc.j*pc.k)

      i = 0
      for (j in 1:pc.j){
         for (k in 1:pc.k){
            i = i +1
            eig.Age.Year[,i] <- as.numeric(outer(psi[,j],phi[, k],"*"))
         }
      }

			lmfull.simu = mFPCA.lm(t(X.age.year.c),eig.Age.Year, ploting=F)
			scores <- lmfull.simu$coefs
      Xest = mean.Xmat + scores%*%t(eig.Age.Year)

			df.eig.Age.Year <- as.data.frame(eig.Age.Year)
     	r <- dim(df.eig.Age.Year)[2]
      mFPCA.Age.Year.pctvar <- numeric(r)
      for (i in 1:r){
        mFPCA.Age.Year.pctvar[i] <-
				     round(mFPCA.lm(t(X.age.year.c),df.eig.Age.Year,model=i,ploting=F)$perct.model,2)
      }

      best.model <- order(-mFPCA.Age.Year.pctvar)
			VarOrdered <- mFPCA.Age.Year.pctvar[best.model]
			namesV = as.character(t(outer(paste("Psi",1:pc.j,sep=""),paste("Phi",1:pc.k,sep=""),paste,sep=".")))
		  names(VarOrdered) = namesV[best.model]
			return(list(Xest = Xest,mu = mean.X, scores = scores, res.psi = res.psi, res.phi = res.phi,
				         eig = eig.Age.Year, pc.j = pc.j, psi = psi, pc.k = pc.k, phi = phi, FVE = mFPCA.Age.Year.pctvar, VarOrdered = VarOrdered))
}



