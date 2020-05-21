library(tidyverse)
library(LFBayes)
library(splines)
library(loo)
load("E:/Rcpp stuff/ERP/ERPdat.RData")
ASD <- ERPdat %>%
  filter(group == 'ASD' & trial_window >= 5) %>%
  mutate(trial_window = trial_window - 4) %>%
  group_by(id, trial_window, erp_time) %>% 
  summarize(y = mean(diff))

Missing_data_ASD <- ASD %>% group_by(id) %>% summarize(unique_elements = n_distinct(trial_window))
sum(Missing_data_ASD[,2]) / (56 * 37)
p1 <- ASD %>%
  filter(id == 113 & trial_window == 18) %>%
  ggplot(aes(x = erp_time)) + 
  geom_line(aes(y = y)) +
  xlab("ERP time") +
  ylab("Condition differentiation") +
  ggtitle("Patient ID = 113, Trial = 18") +
  theme_bw()
p2 <- ASD %>%
  filter(id == 115 & trial_window == 5) %>%
  ggplot(aes(x = erp_time)) + 
  geom_line(aes(y = y)) +
  xlab("ERP time") +
  ylab("Condition differentiation") +
  ggtitle("Patient ID = 115, Trial = 5") +
  theme_bw()
p3 <- ASD %>%
  filter(id == 204 & trial_window == 35) %>%
  ggplot(aes(x = erp_time)) + 
  geom_line(aes(y = y)) +
  xlab("ERP time") +
  ylab("Condition differentiation") +
  ggtitle("Patient ID = 204, Trial = 35") +
  theme_bw()

ASD_subjects <- pull(ASD %>%
  ungroup(id, trial_window, erp_time) %>%
  distinct(id))
TD <- ERPdat %>%
  filter(group == 'TD' & trial_window >= 5) %>%
  mutate(trial_window = trial_window - 4) %>%
  group_by(id, trial_window, erp_time) %>% 
  summarize(y = mean(diff))

Missing_data_TD <- TD %>% group_by(id) %>% summarize(unique_elements = n_distinct(trial_window))
sum(Missing_data_TD[,2]) / (56 * 34)
ASD_full <- tibble(id = rep(ASD_subjects, each = 56),
                   trial_window = rep(1:56, length(ASD_subjects)))
ASD_missing <- ASD %>%
  group_by(id) %>%
  distinct(trial_window) %>%
  setdiff(x = ASD_full)

missing_list <- lapply(1:length(ASD_subjects), function(i)
  pull(ASD_missing %>%
  filter(id == ASD_subjects[i]) %>%
  select(trial_window)))

yy <- lapply(1:length(ASD_subjects), function(i)
  pull(ASD %>%
       filter(id == ASD_subjects[i]) %>%
       select(y)))
mx <- mean(unlist(yy))
sx <- sd(unlist(yy))
yy <- lapply(1:length(ASD_subjects), function(i) (yy[[i]]-mx)/sx)
iter <- 5000
nchain <- 1
burnin <- 1000
p1 <- 14
p2 <- 50
q1 <- 10
q2 <- 10
bs1 <- bs(seq(from = 0, to = 1, length.out = 56), df = p2, intercept = TRUE)
bt1 <- bs(seq(from = 0, to = 1, length.out = 37), df = p1, intercept = TRUE)
X <- rep(1,length(ASD_subjects))
dim(X) <- c(length(ASD_subjects),1)
mcmc <- mcmcWeakChains(yy, missing_list, X, bs1, bt1, q1, q2, iter, 1, burnin, 1)
L <- loglik(unlist(yy), X, bs1, bt1, missing_list, mcmc$Theta, mcmc$Varphi, iter, burnin)
rel_n_eff <- relative_eff(exp(L), chain_id = rep(1, each = 4000))
loo(L, r_eff = rel_n_eff, cores = 4)
