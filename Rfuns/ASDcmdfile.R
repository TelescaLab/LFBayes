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

#counter <- 1
#longmeans <- matrix(0, ncol = length(subjects), nrow = 56)
#for(s in unique(mydf$id)){
#  print(s)
#  for(t in 5:60){
#    longmeans[t-4,counter] <- mean(mydf$diff[mydf$id==s & mydf$trial_window==t])
#  }
#  counter <- counter + 1
#}

#counter <- 0
#myid <- 7
#bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = 15, intercept = TRUE)
#bs1 <- predict(bs1, seq(from = 0, to = 1, length.out= 56)[-which(is.na(longmeans[,myid]))])
#reg1 <- lm(na.omit(longmeans[,myid]) ~ bs1-1)
#plot(na.omit(longmeans[,myid]),type="l")
#lines(reg1$fitted.values,col="red")

iter <- 5000
nchain <- 4
burnin <- 1000

p1 <- ceiling(37)
p2 <- ceiling(56)
q1 <- floor(p1/2)
q2 <- floor(p2/2)
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
ASDERPFebruaryMCMC7 <- mcmcWeakChains(y, missing, X, bs1, bt1, q1, q2, iter, 1, burnin, nchain)
save(ASDERPFebruaryMCMC7, file = "ASDERPFebruaryMCMC7.RData")
ASDERPFebruaryDIC7 <- calculate_DIC_Missing(y, X, observed, ASDERPFebruaryMCMC7, bs1, bt1, burnin, 10)
save(ASDERPFebruaryDIC7, file = "ASDERPFebruaryDIC7.RData")
ASDERPFebruaryBIC7 <- calculate_BIC_Missing(y, X, observed, ASDERPFebruaryMCMC7, bs1, bt1, burnin, 10)
save(ASDERPFebruaryBIC7, file = "ASDERPFebruaryBIC7.RData")
ASDERPFebruaryeigenMCMC7 <- eigenLFChains(bs1, bt1, ASDERPFebruaryMCMC7, 3, iter, burnin, nchain)
save(ASDERPFebruaryeigenMCMC7, file = "ASDERPFebruaryeigenMCMC7.RData")
rm(ASDERPFebruaryMCMC7)
rm(ASDERPFebruaryeigenMCMC7)
