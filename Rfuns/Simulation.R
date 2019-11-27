##############################################################################################
# Load relevant libraries

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
source("SimulationSetup.R")
##############################################################################################
# True data generating parameters

Errorvar <- .025
SS <- 10
TT <- 20
t <- seq(from = 0, to = 1, length.out = TT)
s <- seq(from = 0, to = 1, length.out = SS)
p1 <- 6
p2 <- 6
q1 <- 4
q2 <- 3
Bt <- bs(t, df = p1, intercept = TRUE)
Bs <- bs(s, df = p2, intercept = TRUE)
H <- GenerateH(q1, q2)

##############################################################################################

Mu <- GenerateMu1(s,t)
Lambda <- Loading.Matern(t, p1, q1, Bt)
Gamma <- Loading.Brown.Bridge(s, p2, q2)
Cov <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + Errorvar * diag(SS * TT)

##############################################################################################
# Uncomment to run case 2

#mu2 <- GenerateMu2(s,t)
#Lambda <- Loading.Brown.Bridge(t, p1, q1)
#Gamma <- Loading.CosCov(s,p2,q2,Bs)
#Cov <- kronecker(Bs%*%Gamma, Bt%*%Lambda)%*%H%*%t(kronecker(Bs%*%Gamma, Bt%*%Lambda)) + errorvar * diag(SS * TT)

##############################################################################################
# Uncomment to run case 3

#mu3 <- GenerateMu3(s,t)
#Cov <- GenerateNonsep(s,t)

##############################################################################################
# Options for PACE

pc.j = NULL
pc.k = NULL
fpca.op1 = list(dataType = "Sparse", maxK = pc.j, FVEthreshold = .9999, nRegGrid = TT)
fpca.op2 = list(dataType = "Sparse", maxK = pc.k, FVEthreshold = .9999, nRegGrid = SS)

##############################################################################################
# Options for MCMC

iter <- 5000 # Number of iterations
burnin <- 1000 # Burnin iterations
thin <- 1 # Thinning for each chain
nchain <- 4 # Number of chains
neig <- 3 # Number of eigenfunctions for inference
info_thin <- 10 # Thinning for information criteria (shortens computation time)

##############################################################################################
# Design parameters used in estimation

N <- 30 # Number of subjects
p1s <- 9
p2s <- 8 # Number of splines used in estimation
q1s <- 6 # Number of latent factors for functional dimension
q2s <- 6 # Number of latent factors for longitudinal dimension
Bt1 <- bs(t, df = p1s, intercept = TRUE)
Bs1 <- bs(s, df = p2s, intercept = TRUE)
set.seed(1)

##############################################################################################
# Setting up multi-core use
# Ignore if running a single job

setwd("/Users/John/Documents/Johnstuff/LFBayes/Rfuns")
numsimulation <- 50
pb <- txtProgressBar(max = numsimulation, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl <- makeCluster(5)
registerDoSNOW(cl)

##############################################################################################
### Array job
# Uses multiple cores to run multiple jobs in parallel
# Saves output to a user-specified directory

Array_Job(numsimulation, Cov, Mu, Errorvar, N, SS, TT, Bs1, Bt1, q1s, q2s, iter, thin, info_thin, burnin, nchain, neig, fpca.op1, fpca.op2, pc.j, pc.k, "test.RData", getwd())

##############################################################################################
### Single job
# Uses all cores to run a single job
# Saves output in R object

output <- Single_Job(Cov, Mu, Errorvar, N, SS, TT, Bs1, Bt1, q1s, q2s, iter, thin, info_thin, burnin, nchain, neig, fpca.op1, fpca.op2, pc.j, pc.k)
