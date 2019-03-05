setwd("/Users/John/Documents/Johnstuff/ERP")
load("ERPdat.RData")
head(ERPdat)
library(LFBayes)

library(splines)
# TD group first
uniqueid <- unique(ERPdat$id[which(ERPdat$group == "TD")])
datTD <- ERPdat[which(ERPdat$group == "TD"
                      & !(ERPdat$id %in% c(217,333))
                      & ERPdat$channel %in% c(4, 118, 123, 124)
                      & !(ERPdat$id == 162 & ERPdat$channel == 118)
                      & !(ERPdat$id == 184 & ERPdat$channel == 124)
                      & ERPdat$trial_window >= 5
                      & ERPdat$trial_window <= 60),]

datASD <- ERPdat[which(ERPdat$group == "ASD"
                       & !(ERPdat$id %in% c(102, 118, 185))
                       & ERPdat$channel %in% c(4, 118, 123, 124)
                       & ERPdat$trial_window >= 5
                       & ERPdat$trial_window <= 60),]


# Mean center data for each trial window within each patient and channel
subjects <- unique(datASD$id)
#for(k in subjects){
#  chann <- chann <- unique(datASD[which(datASD$id == k),]$channel)
#  for(i in chann){
#    tri <- unique(datASD[which(datASD$id == k
#                               & datASD$channel == i),]$trial_window)
#    for(j in tri){
#      indicies <- which(datASD$id == k
#                        & datASD$channel == i
#                        & datASD$trial_window == j)
#      datASD$diff[indicies] <- datASD$diff[indicies] - mean(datASD$diff[indicies])
#    }
#  }
#  print(k)
#}

plot(datASD$diff[which(datASD$id==subjects[1] & datASD$channel == 4 & datASD$trial_window == 5)], type = "l")
# Figure out how to average across channels

mywin <- 6
plot(datASD$diff[which(datASD$id == subjects[1] & datASD$channel == 118 & datASD$trial_window == mywin)], type = "l", ylim = c(-20,0))
lines(datASD$diff[which(datASD$id == subjects[1] & datASD$channel == 124 & datASD$trial_window == mywin)])
lines(datASD$diff[which(datASD$id == subjects[1] & datASD$channel == 123 & datASD$trial_window == mywin)])
lines(datASD$diff[which(datASD$id == subjects[1] & datASD$channel == 4 & datASD$trial_window == mywin)])

mymean <- (datASD$diff[which(datASD$id == subjects[1] & datASD$channel == 4 & datASD$trial_window == mywin)] +
             datASD$diff[which(datASD$id == subjects[1] & datASD$channel == 123 & datASD$trial_window == mywin)] +
             datASD$diff[which(datASD$id == subjects[1] & datASD$channel == 124 & datASD$trial_window == mywin)] +
             datASD$diff[which(datASD$id == subjects[1] & datASD$channel == 118 & datASD$trial_window == mywin)]) / 4

lines(mymean, col = "red")

# Average across channels
id <- numeric(0)
win <- numeric(0)
diff <- numeric(0)
A <- numeric(0)
B <- numeric(0)
C <- numeric(0)
D <- numeric(0)
y <- list()
missing <- list()
i <- 1
for(k in subjects){
  getwin <- Reduce(union, list(datASD$trial_window[which(datASD$id == k & datASD$channel == 4)],
                               datASD$trial_window[which(datASD$id == k & datASD$channel == 118)],
                               datASD$trial_window[which(datASD$id == k & datASD$channel == 123)],
                               datASD$trial_window[which(datASD$id == k & datASD$channel == 124)]))
  missing[[i]] <- Reduce(setdiff, sort(getwin), 5:60) - 4
  for(g in sort(getwin)){
    A <- datASD$diff[which(datASD$id== k & datASD$trial_window == g & datASD$channel== 4)]
    B <- datASD$diff[which(datASD$id== k & datASD$trial_window == g & datASD$channel== 118)]
    C <- datASD$diff[which(datASD$id== k & datASD$trial_window == g & datASD$channel== 123)]
    D <- datASD$diff[which(datASD$id== k & datASD$trial_window == g & datASD$channel== 124)]
    #print(length(A))
    #print(length(B))
    #print(length(C))
    #print(length(D))
    
    if(length(A) != 0 & length(A) != 37){
      print(paste0("PROBLEM CHILD IS ", k,"A"))
    }
    if(length(B) != 0 & length(B) != 37){
      print(paste0("PROBLEM CHILD IS ", k,"B"))
    }
    if(length(C) != 0 & length(C) != 37){
      print(paste0("PROBLEM CHILD IS ", k,"C"))
    }
    if(length(D) != 0 & length(D) != 37){
      print(paste0("PROBLEM CHILD IS ", k,"D",g))
    }
    if(length(A) == 0){
      A <- rep(0, 37)
    }
    if(length(B) == 0){
      B <- rep(0, 37)
    }
    if(length(C) == 0){
      C <- rep(0, 37)
    }
    if(length(D) == 0){
      D <- rep(0, 37)
    }
    tempdiff <- (A + B + C + D) / 4
    id <- c(id, rep(k, 37))
    win <- c(win, rep(g, 37))
    diff <- c(diff, tempdiff)
  }
  y[[i]] <- tempdiff
  print(k)
  i <- i + 1
}
mydf <- data.frame(id = id, trial_window = win, diff = diff)
y <- lapply(1:length(subjects), function(i) mydf$diff[which(mydf$id == subjects[i])])
meany <- mean(unlist(y))
sdy <- sd(unlist(y))
y <- lapply(1:length(subjects), function(i) (y[[i]] - meany) / sdy)
X <- rep(1,length(subjects))
dim(X) <- c(length(subjects),1)
observed <- list()
for(i in 1:length(subjects)){
  observed[[i]] <- numeric(0)
  for(m in setdiff(1:56, missing[[i]])){
    observed[[i]] <- c(observed[[i]], ((m-1)*37):(m*37 - 1))
  }
}

iter <- 5000
nchain <- 4
burnin <- 1000

p1 <- ceiling(37/4)
p2 <- ceiling(56/4)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
load("ASDERPFebruaryMCMC1.RData")
#ASDERPFebruaryMCMC1 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
#save(ASDERPFebruaryMCMC1, file = "ASDERPFebruaryMCMC1.RData")
#ASDERPFebruaryDIC1 <- calculate_DIC_Missing(y, X, observed, ASDERPFebruaryMCMC1, bs1, bt1, burnin, 10)
#save(ASDERPFebruaryDIC1, file = "ASDERPFebruaryDIC1.RData")
ASDERPFebruaryBIC1 <- calculate_BIC_Missing(y, X, observed, ASDERPFebruaryMCMC1, bs1, bt1, burnin, 10)
save(ASDERPFebruaryBIC1, file = "ASDERPFebruaryBIC1.RData")
ASDERPFebruaryeigenMCMC1 <- eigenLFChains(bs1, bt1, ASDERPFebruaryMCMC1, 3, iter, burnin, nchain)
save(ASDERPFebruaryeigenMCMC1, file = "ASDERPFebruaryeigenMCMC1.RData")
rm(ASDERPFebruaryMCMC1)
rm(ASDERPFebruaryeigenMCMC1)

p1 <- ceiling(37/3)
p2 <- ceiling(56/3)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
ASDERPFebruaryMCMC2 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
save(ASDERPFebruaryMCMC2, file = "ASDERPFebruaryMCMC2.RData")
ASDERPFebruaryDIC2 <- calculate_DIC_Missing(y, X, observed, ASDERPFebruaryMCMC2, bs1, bt1, burnin, 10)
save(ASDERPFebruaryDIC2, file = "ASDERPFebruaryDIC2.RData")
ASDERPFebruaryBIC2 <- calculate_BIC_Missing(y, X, observed, ASDERPFebruaryMCMC2, bs1, bt1, burnin, 10)
save(ASDERPFebruaryBIC2, file = "ASDERPFebruaryBIC2.RData")
ASDERPFebruaryeigenMCMC2 <- eigenLFChains(bs1, bt1, ASDERPFebruaryMCMC2, 3, iter, burnin, nchain)
save(ASDERPFebruaryeigenMCMC2, file = "ASDERPFebruaryeigenMCMC2.RData")
rm(ASDERPFebruaryMCMC2)
rm(ASDERPFebruaryeigenMCMC2)

p1 <- ceiling(37/2)
p2 <- ceiling(56/2)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
ASDERPFebruaryMCMC3 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
save(ASDERPFebruaryMCMC3, file = "ASDERPFebruaryMCMC3.RData")
ASDERPFebruaryDIC3 <- calculate_DIC_Missing(y, X, observed, ASDERPFebruaryMCMC3, bs1, bt1, burnin, 10)
save(ASDERPFebruaryDIC3, file = "ASDERPFebruaryDIC3.RData")
ASDERPFebruaryBIC3 <- calculate_BIC_Missing(y, X, observed, ASDERPFebruaryMCMC3, bs1, bt1, burnin, 10)
save(ASDERPFebruaryBIC3, file = "ASDERPFebruaryBIC3.RData")
ASDERPFebruaryeigenMCMC3 <- eigenLFChains(bs1, bt1, ASDERPFebruaryMCMC3, 3, iter, burnin, nchain)
save(ASDERPFebruaryeigenMCMC3, file = "ASDERPFebruaryeigenMCMC3.RData")
rm(ASDERPFebruaryMCMC3)
rm(ASDERPFebruaryeigenMCMC3)

p1 <- ceiling(37*2/3)
p2 <- ceiling(56*2/3)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
ASDERPFebruaryMCMC4 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
save(ASDERPFebruaryMCMC4, file = "ASDERPFebruaryMCMC4.RData")
ASDERPFebruaryDIC4 <- calculate_DIC_Missing(y, X, observed, ASDERPFebruaryMCMC4, bs1, bt1, burnin, 1)
save(ASDERPFebruaryDIC4, file = "ASDERPFebruaryDIC4.RData")
ASDERPFebruaryBIC4 <- calculate_BIC_Missing(y, X, observed, ASDERPFebruaryMCMC4, bs1, bt1, burnin, 10)
save(ASDERPFebruaryBIC4, file = "ASDERPFebruaryBIC4.RData")
ASDERPFebruaryeigenMCMC4 <- eigenLFChains(bs1, bt1, ASDERPFebruaryMCMC4, 3, iter, burnin, nchain)
save(ASDERPFebruaryeigenMCMC4, file = "ASDERPFebruaryeigenMCMC4.RData")
rm(ASDERPFebruaryMCMC4)
rm(ASDERPFebruaryeigenMCMC4)

p1 <- ceiling(37*3/4)
p2 <- ceiling(56*3/4)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
ASDERPFebruaryMCMC5 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
save(ASDERPFebruaryMCMC5, file = "ASDERPFebruaryMCMC5.RData")
ASDERPFebruaryDIC5 <- calculate_DIC_Missing(y, X, observed, ASDERPFebruaryMCMC5, bs1, bt1, burnin, 10)
save(ASDERPFebruaryDIC5, file = "ASDERPFebruaryDIC5.RData")
ASDERPFebruaryBIC5 <- calculate_BIC_Missing(y, X, observed, ASDERPFebruaryMCMC5, bs1, bt1, burnin, 10)
save(ASDERPFebruaryBIC5, file = "ASDERPFebruaryBIC5.RData")
ASDERPFebruaryeigenMCMC5 <- eigenLFChains(bs1, bt1, ASDERPFebruaryMCMC5, 3, iter, burnin, nchain)
save(ASDERPFebruaryeigenMCMC5, file = "ASDERPFebruaryeigenMCMC5.RData")
rm(ASDERPFebruaryMCMC5)
rm(ASDERPFebruaryeigenMCMC5)

p1 <- ceiling(37*4/5)
p2 <- ceiling(56*4/5)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
ASDERPFebruaryMCMC6 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
save(ASDERPFebruaryMCMC6, file = "ASDERPFebruaryMCMC6.RData")
ASDERPFebruaryDIC6 <- calculate_DIC_Missing(y, X, observed, ASDERPFebruaryMCMC6, bs1, bt1, burnin, 10)
save(ASDERPFebruaryDIC6, file = "ASDERPFebruaryDIC6.RData")
ASDERPFebruaryBIC6 <- calculate_BIC_Missing(y, X, observed, ASDERPFebruaryMCMC6, bs1, bt1, burnin, 10)
save(ASDERPFebruaryBIC6, file = "ASDERPFebruaryBIC6.RData")
ASDERPFebruaryeigenMCMC6 <- eigenLFChains(bs1, bt1, ASDERPFebruaryMCMC6, 3, iter, burnin, nchain)
save(ASDERPFebruaryeigenMCMC6, file = "ASDERPFebruaryeigenMCMC6.RData")
rm(ASDERPFebruaryMCMC6)
rm(ASDERPFebruaryeigenMCMC6)
