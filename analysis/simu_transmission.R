#### SIMULATION OF TRANSMISSION PAIRS

### Simulate SNP cutoff
library(scales)
library(truncnorm)
library(reshape2)
library(ggplot2)

### Home
home <- "~/Documents/snp_cutoff/"
plots <- "~/Documents/plots"
theme_set(theme_bw(base_size = 24))

###*** Functions ***############################################################################################################################################################
setwd(home)
source("analysis/simu_transmission_fn.R")

###*** Times ***############################################################################################################################################################
times <- read.csv("../data/sample_time_dist2.csv")[,-1] ### t distribution from FC
times_c <- 180 # 6 months to match FC

###*** mutation rate ***###################################################################################################################################################
mu = 5 / 365 # 5 SNPS per year
# 5.21*10^(-9)*2842618*365 = 5.4 from Campbell 2018
mu_e = 5.41 / 365
mu_l = 3.79 / 365
mu_h = 6.37 / 365
# since 99.7% of normal values fall within 3 SD
mn = mu_l + (mu_h - mu_l)/2
sdmn = mn/3

###*** RUNS ***###################################################################################################################################################
setwd("plots")
npat = 1000
nruns = 10000

###*** CONSTANT 6mo***###################################################################################################################################################
s_constant <- simu_runs(180, mu, npat, nruns,0,0)
p <- ggplot(s_constant$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,position = "identity") + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("constant_6mo_SNP_dist_limits.pdf")
p1 <- ggplot(s_constant$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,alpha = 0.4,position = "identity") + 
  scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("grp_constant_6mo_SNP_dist_limits.pdf")

###*** VARYING TIMES ***###################################################################################################################################################
s_varying <- simu_runs(times, mu, npat, nruns,0,1)
p <- ggplot(s_varying$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,position = "identity") + 
   scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability") + facet_wrap(~variable) 
ggsave("varying_SNP_dist_limits.pdf")
p1 <- ggplot(s_varying$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,alpha = 0.4,position = "identity") + 
  scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("grp_varying_SNP_dist_limits.pdf")

###*** VARYING TIMES & MU ***###################################################################################################################################################
s_varyingtmu <- simu_runs(times, c(mn,sdmn), npat, nruns,0,0)
p <- ggplot(s_varyingtmu$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,position = "identity") + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("varyingtmu_SNP_dist_limits.pdf")
p1 <- ggplot(s_varyingtmu$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,alpha = 0.4,position = "identity") + 
  scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("grp_varyingtmu_SNP_dist_limits.pdf")

###*** CONSTANT & SAME MODEL FOR SOURCE/RECIPIENT ***###################################################################################################################################################
s_constant_same <- simu_runs(times_c, mu, npat, nruns,1,0)
p <- ggplot(s_constant_same$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,position = "identity") + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("constant_same_SNP_dist_limits.pdf")
p1 <- ggplot(s_constant_same$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,alpha = 0.4,position = "identity") + 
  scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("grp_constant_same_SNP_dist_limits.pdf")

###*** VARYING & SAME MODEL FOR SOURCE/RECIPIENT ***###################################################################################################################################################
s_vary_same <- simu_runs(times, c(mn,sdmn), npat, nruns,1,0)
p <- ggplot(s_vary_same$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,position = "identity") + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("vary_same_SNP_dist_limits.pdf")
p1 <- ggplot(s_vary_same$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,alpha = 0.4,position = "identity") + 
  scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("grp_vary_same_SNP_dist_limits.pdf")

###*** CONSTANT 3mo***###################################################################################################################################################
s_constant3 <- simu_runs(90, mu, npat, nruns,0,0)
p <- ggplot(s_constant3$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,position = "identity") + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("constant_3mo_SNP_dist_limits.pdf")
p1 <- ggplot(s_constant3$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,alpha = 0.4,position = "identity") + 
  scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")
ggsave("grp_constant_3mo_SNP_dist_limits.pdf")

###*** CONSTANT 12mo***###################################################################################################################################################
s_constant12 <- simu_runs(365, mu, npat, nruns,0,0)
p <- ggplot(s_constant12$store,aes(x=value,fill = variable))  + geom_histogram(aes(y = ..density..),binwidth = 1,position = "identity") + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance",limits = c(0,45)) + scale_y_continuous("Probability")
ggsave("constant_12mo_SNP_dist_limits.pdf")
p1 <- ggplot(s_constant12$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,alpha = 0.4,position = "identity") + 
  scale_x_continuous("SNP distance",limits = c(0,100)) + scale_y_continuous("Probability")
ggsave("grp_constant_12mo_SNP_dist_limits.pdf")

###*** CONSTANT 3years***###################################################################################################################################################
s_constant3y <- simu_runs(3*365, mu, npat, nruns,0,0)
p <- ggplot(s_constant3y$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,position = "identity") + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance",limits = c(0,100)) + scale_y_continuous("Probability")
ggsave("constant_3yr_SNP_dist_limits.pdf")
p1 <- ggplot(s_constant3y$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,alpha = 0.4,position = "identity") + 
  scale_x_continuous("SNP distance",limits = c(0,100)) + scale_y_continuous("Probability")
ggsave("grp_constant_3yr_SNP_dist_limits.pdf")

### To get underlying data use:
# underlying data
gg <- ggplot_build(p)
gg_data <- gg$data[[1]]
gg_data <- gg_data[,c("x","density","PANEL")]
gg_data <- dcast(gg_data, x ~ PANEL, value.var = "density" )
# thresholds: last non-zero (subtracted 1 as first row zero)
c(tail(which(gg_data[,c(2)]!=0),1),tail(which(gg_data[,c(3)]!=0),1),tail(which(gg_data[,c(4)]!=0),1),tail(which(gg_data[,c(5)]!=0),1)) - 1

### group together
constant_all <- as.data.frame(rbind(s_constant$store,s_constant3$store, s_constant12$store))
constant_all$time <- c(rep(6,dim(s_constant$store)[1]),rep(3,dim(s_constant$store)[1]),rep(12,dim(s_constant$store)[1]))
ggplot(constant_all, aes(x=value, fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,alpha = 0.7,position = "identity") + 
  scale_x_continuous("SNP distance",limits = c(0,40)) + scale_y_continuous("Probability") + facet_wrap(~time)
ggsave("grp_constant_all_SNP_dist_limits.pdf")              
                       
