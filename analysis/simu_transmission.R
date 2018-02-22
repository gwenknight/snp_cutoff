#### SIMULATION OF TRANSMISSION PAIRS

### Simulate SNP cutoff
library(scales)
library(truncnorm)
library(reshape2)

### Home
home <- "~/Documents/snp_cutoff/"
plots <- "~/Documents/plots"
theme_set(theme_bw(base_size = 24))

###*** Functions ***############################################################################################################
# Sample from own distribution
rMydist <- function(n,x) {
  sample(seq(1,length(x),1), size = n, prob = x, replace=T)
}
# Sample from distribution 
rMydist_plus <- function(n,x,n_mut) {
  sample(seq(1,length(x)+n_mut,1), size = n, prob = c(matrix(0,1,n_mut),x), replace=T) # zero probability of less than n_muts
}
## General model for how distance from transmitted changes over time
gen_mod_day <- function(day,l_x){
  x <- seq(0,l_x,1)
  aa_g <-  0.99338 #mod_a_mod[1]
  bb_g <- -0.02207 #mod_a_mod[2]
  out <- exp(exp(aa_g + bb_g*day) + p_aa * x)
  return(out)
}

###*** Times ***############################################################################################################################################################
times <- read.csv("../data/sample_time_dist2.csv")[,-1]

###*** mutation rate ***###################################################################################################################################################
mu = 5 / 365 # 5 SNPS per year
# 5.21*10^(-9)*2842618*365 = 5.4 from Campbell 2018
mu_e = 5.41 / 365
mu_l = 3.79 / 365
mu_h = 6.37 / 365
# since 99.7% of normal values fall within 3 SD
mn = mu_l + (mu_h - mu_l)/2
sdmn = mn/3

#### RUNS
store_all <- c()
store_limits <- c()

for(j in 1:1000){
  
  ### STORE / PREP
  n = 1000 # number of pairs of samples
  
  store <- matrix(0,n,2)
  
  # SAMPLE mutation rates to use
  mu_v <- rtruncnorm(n,a = 0, b = Inf, mn,sdmn) # lower truncated to be positive
  
  # SAMPLE times to use
  ht <- hist(times,breaks = seq(0,500,1),plot = FALSE)
  times_v <- rMydist(n,ht$counts) 
  
  for(i in 1:n){
    sample_gap <- times_v[i]
    store[i,1] <- sample_gap # store time
    
    # SOURCE
    dis_source <- gen_mod_day(sample_gap,150) 
    dis_source_from_transm <- rMydist(1,dis_source)
    
    # RECEIVER
    mut_number <- rpois(1,mu_v[i] * sample_gap) # number of mutations likely to have occured since sample
    dis_receiver_from_transm <- rexp(1,1/mut_number) # 1/mean (rate in this definition of exp)
    
    # DISTANCE SOURCE <- TRANS -> RECEIVER
    mut_dist <- dis_source_from_transm + dis_receiver_from_transm
    store[i,2] <- mut_dist # store mut distance
  }
  
  #plot(store[,1], store[,2], xlab = "Time Distance (in days)", ylab = "SNP distance")
  h<-hist(store[,2],breaks = seq(1,400,1),plot=FALSE)
  h_cum <- cumsum(h$density)
  #plot(seq(1,39,1),h_cum)
  l1<-length(which(h_cum < 0.7)) # 70% 1 SNP differences
  l2<-length(which(h_cum < 0.9)) # 90% 4 SNPS different
  l3<-length(which(h_cum < 0.95)) # 95% 6 SNPS different 
  l4<-length(which(h_cum < 0.99)) # 99% 14 SNPS different
  
  ### Store
  store_all <- rbind(store_all,cbind(j,store))
  store_limits <- rbind(store_limits,c(j,l1,l2,l3,l4))
  
}

store_limits <- as.data.frame(store_limits)
colnames(store_limits) <- c("run","l70","l90","l95","l99")
m_store_limits <- melt(store_limits,id.vars = "run")
ggplot(m_store_limits,aes(x=value)) + geom_histogram(stat="count") + facet_wrap(~variable)
