W <- calculate_WAIC(t(x), X, mcmc, Bs1, Bt1, 5000)

for(i in 1:17){
  xt <- matrix(x[i,],nrow=44,ncol=56)
  xt <- xt[-c(seq(from = 2, to = 44, by = 2)), -c(seq(from = 2, to = 56, by = 2))]
  xthin[i,]<-c(xt)
}


p1 <- 44
p2 <- 50
q1 <- 20
q2 <- 20
bs1 <- bs(seq(from=0,to=1,length.out=56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from=0,to=1,length.out=44), df = p1, intercept = TRUE)
X <- rep(1,17)
dim(X) <- c(17,1)
MCMC1_5000 <- mcmcWeakChains(y,missing,X,bs1,bt1,q1,q2,5000,1,1000, 1)
bt1thin <- predict(bt1, seq(from = 0, to = 1, length.out = 44)[seq(from=1,to=43,by=2)])
bs1thin <- predict(bs1, seq(from = 0, to = 1, length.out = 56)[seq(from = 1, to = 55, by = 2)])
W1 <- calculate_WAIC(t(xthin), X, MCMC1_5000, bs1thin, bt1thin, 1000)
save(MCMC1_5000, file = "MCMC1_5000.RData")
save(W1, file = "W1.RData")

p1 <- 44
p2 <- 50
q1 <- 6
q2 <- 6
bs1 <- bs(seq(from=0,to=1,length.out=56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from=0,to=1,length.out=44), df = p1, intercept = TRUE)
X <- rep(1,17)
dim(X) <- c(17,1)
MCMC2_5000 <- mcmcWeakChains(y,missing,X,bs1,bt1,q1,q2,5000,1,5000, 1)
bt1thin <- predict(bt1, seq(from = 0, to = 1, length.out = 44)[seq(from=1,to=43,by=2)])
bs1thin <- predict(bs1, seq(from = 0, to = 1, length.out = 56)[seq(from = 1, to = 55, by = 2)])
W2test <- calculate_WAIC(t(xthin), X, MCMC2_5000, bs1thin, bt1thin, 1000)
save(MCMC2_5000, file = "MCMC2_5000.RData")
save(W2, file = "W2.RData")

p1 <- 16
p2 <- 20
q1 <- 12
q2 <- 12
bs1 <- bs(seq(from=0,to=1,length.out=56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from=0,to=1,length.out=44), df = p1, intercept = TRUE)
X <- rep(1,17)
dim(X) <- c(17,1)
MCMC3_5000 <- mcmcWeakChains(y,missing,X,bs1,bt1,q1,q2,5000,1,5000, 1)
bt1thin <- predict(bt1, seq(from = 0, to = 1, length.out = 44)[seq(from=1,to=43,by=2)])
bs1thin <- predict(bs1, seq(from = 0, to = 1, length.out = 56)[seq(from = 1, to = 55, by = 2)])
W3test <- calculate_WAIC(t(xthin), X, MCMC3_5000, bs1thin, bt1thin, 1000)
save(MCMC3_5000, file = "MCMC3_5000.RData")
save(W3, file = "W3.RData")

p1 <- 16
p2 <- 20
q1 <- 6
q2 <- 6
bs1 <- bs(seq(from=0,to=1,length.out=56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from=0,to=1,length.out=44), df = p1, intercept = TRUE)
X <- rep(1,17)
dim(X) <- c(17,1)
MCMC4_5000 <- mcmcWeakChains(y,missing,X,bs1,bt1,q1,q2,5000,1,5000, 1)
bt1thin <- predict(bt1, seq(from = 0, to = 1, length.out = 44)[seq(from=1,to=43,by=2)])
bs1thin <- predict(bs1, seq(from = 0, to = 1, length.out = 56)[seq(from = 1, to = 55, by = 2)])
W4 <- calculate_WAIC(t(xthin), X, MCMC4_5000, bs1thin, bt1thin, 1000)
save(MCMC4_5000, file = "MCMC4_5000.RData")
save(W4, file = "W4.RData")

