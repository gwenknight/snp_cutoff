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
nruns = 1000

###*** CONSTANT ***###################################################################################################################################################
s_constant <- simu_runs(times_c, mu, npat, nruns,0)
p <- ggplot(s_constant$store,aes(x=value)) + geom_histogram(aes(y = ..density..),binwidth = 1) + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance") + scale_y_continuous("Probability")
ggsave("constant_SNP_dist_limits.pdf")

###*** VARYING TIMES ***###################################################################################################################################################
s_varying <- simu_runs(times, mu, npat, nruns,0)
p <- ggplot(s_varying$store,aes(x=value)) + geom_histogram(aes(y = ..density..),binwidth = 1) + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance") + scale_y_continuous("Probability")
ggsave("varying_SNP_dist_limits.pdf")

###*** VARYING TIMES & MU ***###################################################################################################################################################
s_varyingtmu <- simu_runs(times, c(mn,sdmn), npat, nruns,0)
p <- ggplot(s_varyingtmu$store,aes(x=value)) + geom_histogram(aes(y = ..density..),binwidth = 1) + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance") + scale_y_continuous("Probability")
ggsave("varyingtmu_SNP_dist_limits.pdf")

###*** CONSTANT & SAME MODEL FOR SOURCE/RECIPIENT ***###################################################################################################################################################
s_constant_same <- simu_runs(times_c, mu, npat, nruns,1)
p <- ggplot(s_constant_same$store,aes(x=value)) + geom_histogram(aes(y = ..density..),binwidth = 1) + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance") + scale_y_continuous("Probability")
ggsave("constant_same_SNP_dist_limits.pdf")

