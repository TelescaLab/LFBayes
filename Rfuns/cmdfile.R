setwd("/Users/John/Documents/Johnstuff/ERP")
library(LFBayes)
load("ERPdat.RData")
head(ERPdat)


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
subjects <- unique(datTD$id)
#for(k in subjects){
#  chann <- chann <- unique(datTD[which(datTD$id == k),]$channel)
#  for(i in chann){
#    tri <- unique(datTD[which(datTD$id == k
#                               & datTD$channel == i),]$trial_window)
#    for(j in tri){
#      indicies <- which(datTD$id == k
#                        & datTD$channel == i
#                        & datTD$trial_window == j)
#      datTD$diff[indicies] <- datTD$diff[indicies] - mean(datTD$diff[indicies])
#    }
#  }
#  print(k)
#}

plot(datTD$diff[which(datTD$id==135 & datTD$channel == 4 & datTD$trial_window == 5)], type = "l")
# Figure out how to average across channels

mywin <- 6
plot(datTD$diff[which(datTD$id == 135 & datTD$channel == 118 & datTD$trial_window == mywin)], type = "l", ylim = c(-20,0))
lines(datTD$diff[which(datTD$id == 135 & datTD$channel == 124 & datTD$trial_window == mywin)])
lines(datTD$diff[which(datTD$id == 135 & datTD$channel == 123 & datTD$trial_window == mywin)])
lines(datTD$diff[which(datTD$id == 135 & datTD$channel == 4 & datTD$trial_window == mywin)])

mymean <- (datTD$diff[which(datTD$id == 135 & datTD$channel == 4 & datTD$trial_window == mywin)] +
             datTD$diff[which(datTD$id == 135 & datTD$channel == 123 & datTD$trial_window == mywin)] +
             datTD$diff[which(datTD$id == 135 & datTD$channel == 124 & datTD$trial_window == mywin)] +
             datTD$diff[which(datTD$id == 135 & datTD$channel == 118 & datTD$trial_window == mywin)]) / 4

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
  getwin <- Reduce(union, list(datTD$trial_window[which(datTD$id == k & datTD$channel == 4)],
                               datTD$trial_window[which(datTD$id == k & datTD$channel == 118)],
                               datTD$trial_window[which(datTD$id == k & datTD$channel == 123)],
                               datTD$trial_window[which(datTD$id == k & datTD$channel == 124)]))
  missing[[i]] <- Reduce(setdiff, sort(getwin), 5:60) - 4
  for(g in sort(getwin)){
    A <- datTD$diff[which(datTD$id== k & datTD$trial_window == g & datTD$channel== 4)]
    B <- datTD$diff[which(datTD$id== k & datTD$trial_window == g & datTD$channel== 118)]
    C <- datTD$diff[which(datTD$id== k & datTD$trial_window == g & datTD$channel== 123)]
    D <- datTD$diff[which(datTD$id== k & datTD$trial_window == g & datTD$channel== 124)]
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
p1 <- ceiling(37*2/3)
p2 <- ceiling(56*2/3)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
load("ERPFebruaryMCMC4.RData")
#ERPFebruaryMCMC4 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
#save(ERPFebruaryMCMC4, file = "ERPFebruaryMCMC4.RData")
#ERPFebruaryDIC4 <- calculate_DIC_Missing(y, X, observed, ERPFebruaryMCMC4, bs1, bt1, burnin, 10)
#save(ERPFebruaryDIC4, file = "ERPFebruaryDIC4.RData")
ERPFebruaryBIC4 <- calculate_BIC_Missing(y, X, observed, ERPFebruaryMCMC4, bs1, bt1, burnin, 10)
save(ERPFebruaryBIC4, file = "ERPFebruaryBIC4.RData")
ERPFebruaryeigenMCMC4 <- eigenLFChains(bs1, bt1, ERPFebruaryMCMC4, 3, iter, burnin, nchain)
save(ERPFebruaryeigenMCMC4, file = "ERPFebruaryeigenMCMC4.RData")
rm(ERPFebruaryMCMC4)
rm(ERPFebruaryeigenMCMC4)

p1 <- ceiling(37*3/4)
p2 <- ceiling(56*3/4)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
ERPFebruaryMCMC5 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
save(ERPFebruaryMCMC5, file = "ERPFebruaryMCMC5.RData")
ERPFebruaryDIC5 <- calculate_DIC_Missing(y, X, observed, ERPFebruaryMCMC5, bs1, bt1, burnin, 10)
save(ERPFebruaryDIC5, file = "ERPFebruaryDIC5.RData")
ERPFebruaryBIC5 <- calculate_BIC_Missing(y, X, observed, ERPFebruaryMCMC5, bs1, bt1, burnin, 10)
save(ERPFebruaryBIC5, file = "ERPFebruaryBIC5.RData")
ERPFebruaryeigenMCMC5 <- eigenLFChains(bs1, bt1, ERPFebruaryMCMC5, 3, iter, burnin, nchain)
save(ERPFebruaryeigenMCMC5, file = "ERPFebruaryeigenMCMC5.RData")
rm(ERPFebruaryMCMC5)
rm(ERPFebruaryeigenMCMC5)

p1 <- ceiling(37*4/5)
p2 <- ceiling(56*4/5)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
ERPFebruaryMCMC6 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
save(ERPFebruaryMCMC6, file = "ERPFebruaryMCMC6.RData")
ERPFebruaryDIC6 <- calculate_DIC_Missing(y, X, observed, ERPFebruaryMCMC6, bs1, bt1, burnin, 10)
save(ERPFebruaryDIC6, file = "ERPFebruaryDIC6.RData")
ERPFebruaryBIC6 <- calculate_BIC_Missing(y, X, observed, ERPFebruaryMCMC6, bs1, bt1, burnin, 10)
save(ERPFebruaryBIC6, file = "ERPFebruaryBIC6.RData")
ERPFebruaryeigenMCMC6 <- eigenLFChains(bs1, bt1, ERPFebruaryMCMC6, 3, iter, burnin, nchain)
save(ERPFebruaryeigenMCMC6, file = "ERPFebruaryeigenMCMC6.RData")
rm(ERPFebruaryMCMC6)
rm(ERPFebruaryeigenMCMC6)