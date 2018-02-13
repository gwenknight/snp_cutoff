### Simulate SNP cutoff
library(scales)
setwd("~/Documents/snp_cutoff/")
theme_set(theme_bw(base_size = 24))

# needed function
rMydist_plus <- function(n,x,n_mut) {
  sample(seq(1,length(x)+n_mut,1), size = n, prob = c(matrix(0,1,n_mut),x), replace=T) # zero probability of less than n_muts
}
rMydist <- function(n,x) {
  sample(seq(1,length(x),1), size = n, prob = x, replace=T)
}

n = 1000 # number of pairs of samples

# mutation rate
mu = 5 / 365 # 5 SNPS per year

# snp dist within a host
dist0 <- read.csv("data/same_day_snp_dist.csv")[,-1]
# acts like a probability for the distance between sampled and transmitted strain

store <- matrix(0,n,2)

for(i in 1:n){
  sample_gap <- round(runif(1,min = 1, max = 365),0) # rounded to days
  store[i,1] <- sample_gap # store time
  mut_number <- rpois(1,mu * sample_gap) # number of mutations likely to have occured since sample
  
  mut_dist <- rMydist(1,dist0,mut_number) # samples mutation distribution from dist0, assuming minimum of mut_number
  store[i,2] <- mut_dist # store mut distance
}

plot(store[,1], store[,2], xlab = "Time Distance (in days)", ylab = "SNP distance")

#####*** Sample actual patient times MATCHED PAIRS ***###################################################################################################
times <- read.csv("data/sample_time_dist.csv")[,-1]
n <- length(times)
store <- matrix(0,n,2)

for(i in 1:n){
  sample_gap <- times[i] # use data times 
  
  store[i,1] <- sample_gap # store time
  mut_number <- rpois(1,mu * sample_gap) # number of mutations likely to have occured since sample
  
  mut_dist <- rMydist_plus(1,dist0,mut_number) # samples mutation distribution from dist0, assuming minimum of mut_number
  store[i,2] <- mut_dist # store mut distance
}

setwd("../plots")
pdf("pairs_snp_dist.pdf",width=10,height = 4)
plot(store[,1], store[,2], xlab = "Time Distance (in days)", ylab = "SNP distance")
dev.off()

h<-hist(store[,2],breaks = seq(1,max(store[,2],1)))
h_cum <- cumsum(h$density)
plot(seq(1,175,1),h_cum[1:175])
length(which(h_cum < 0.9)) # SNP cutoff of 7
length(which(h_cum < 0.95)) # SNP cutoff of 15

store <- as.data.frame(store)
colnames(store) <- c("time","liSNPdis")

ggplot(store, aes(x = liSNPdis)) + geom_histogram(aes(y=..density..), binwidth = 1) + 
  scale_x_continuous("SNP distance")
setwd("../plots")
ggsave("sample_time_store_output.pdf")

ggplot(dd_wo, aes(x = liSNPdis)) + geom_histogram(binwidth = 1, alpha = 0.4, aes(y=..density..)) + 
  geom_histogram(data = store,aes(x=liSNPdis,y=..density..), binwidth = 1, fill = "red", alpha = 0.4) + 
  scale_x_continuous(limits = c(0,60))
# are those dd_wo not overlapping with red the non-transmission ones? 
  
#####*** NON-MATCHED pairs***###################################################################################################
distbeta <- read.csv("data/22_snp_dist.csv")[,-1]
storebeta <- matrix(0,n,2)

storebeta[,1] <- times
storebeta[,2] <- rMydist(n,distbeta) # sample n times from the distbeta distribution


setwd("../plots")
pdf("pairs_snp_dist.pdf",width=10,height = 4)
plot(storebeta[,1], storebeta[,2], xlab = "Time Distance (in days)", ylab = "SNP distance")
dev.off()

h<-hist(storebeta[,2],breaks = seq(1,max(storebeta[,2],1)))
h_cum <- cumsum(h$density)
plot(seq(1,175,1),h_cum[1:175])
length(which(h_cum < 0.9)) # SNP cutoff of 191
length(which(h_cum < 0.95)) # SNP cutoff of 449

storebeta <- as.data.frame(storebeta)
colnames(storebeta) <- c("time","liSNPdis")

ggplot(storebeta, aes(x = liSNPdis)) + geom_histogram(aes(y=..density..), binwidth = 1) + 
  scale_x_continuous("SNP distance")
setwd("../plots")
ggsave("beta_sample_time_store_output.pdf")

ggplot(dd_wo, aes(x = liSNPdis)) + geom_histogram(binwidth = 1, alpha = 0.4, aes(y=..density..)) + 
  geom_histogram(data = storebeta,aes(x=liSNPdis,y=..density..), binwidth = 1, fill = "red", alpha = 0.4) + 
  scale_x_continuous(limits = c(0,60))
# are those dd_wo not overlapping with red the non-transmission ones? 

