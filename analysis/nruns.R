# how many simulations? 

library(reshape2)
setwd("~/Documents/snp_cutoff/")
source("analysis/simu_transmission_fn.R")

# Choose which mapping 
for(kk in 4){
  print(c("kk = ", kk))
  if(kk == 1){map <- "CC22"} # ST22 strain HO 5096 0412 mapping data 
  if(kk == 2){map <- "CC30"}
  if(kk == 3){map <- "CC22_core"}
  if(kk == 4){map <- "CC30_core"}
  
  ## mu = substitution_rate
  if(map == "CC22"){mu = 4 / 365}
  if(map == "CC30"){mu = 5 / 365}
  if(map == "CC22_core"){mu = 3 / 365}
  if(map == "CC30_core"){mu = 3 / 365}
  
  ## Data for time zero transferred variability
  ## Which sheet? 
  if(map == "CC22"){sheet_num = 1}
  if(map == "CC30"){sheet_num = 2}
  if(map == "CC22_core"){sheet_num = 3}
  if(map == "CC30_core"){sheet_num = 4}
  dataS1 = read.xls(dataS1_file, sheet = sheet_num, header = T) 
  # Removing outliers
  dataS1 = dataS1[-which(grepl("outlier",dataS1$Note)==TRUE),]
  # Same day 
  dd0 = dataS1[which(dd$TimeGap == 0),] # 104 pairs
  h <- hist(dd0$SNPs)
  t0prob_dist <- h$counts/sum(h$counts)
  
  # Parameters from model fit
  param_general_fit <- read.csv(paste0("~/Documents/snp_cutoff/output/param_general_fit_",map,".csv"))[,2]
  
  # Simulation population 
  #vnpats = c(100,500,1000,5000,10000)
  #vnruns = c(100000,100000,100000,100000,100000,
  #           200000,200000,200000,200000,200000)
  
  param_general_fit <- read.csv(paste0("~/Documents/snp_cutoff/output/param_general_fit_",map,".csv"))[,2]
  #param_general_fit <- read.csv(paste0("~/Dropbox/MRC SD Fellowship/Research/SNP cutoffs/output/param_general_fit_",map,".csv"))[,2]
  
  # Simulation population 
  #vnpats = c(100,500,1000,5000,10000)
  nrep = 10
  vnruns <- c(rep(1e4, nrep), rep(5e4, nrep), rep(1e5, nrep), rep(2e5, nrep), rep(5e5, nrep))
  
  #store99 <- c()
  store95 <- c()
  
  #for(i in 1:length(vnpats)){
  #  i = 2
  
  for(j in 1:length(vnruns)){
    print(c(j))
    npat = 459 #vnpats[i]
    nruns = vnruns[j]
    
    ndays = 180 # time between samples
    
    # Run the simulation 
    ss <- simu_runs(ndays,mu,npat,nruns, param_general_fit, t0prob_dist)
    
    # Maximum number of SNPs needed to capture 99% of the transmission events
    #store99 <- c(store99, max(ss$store_limits[which(ss$store_limits$variable == "99%"),"value"]))
    store95 <- c(store95, max(ss$store_limits[which(ss$store_limits$variable == "95%"),"value"]))
  }
  
  
  #store <- cbind(store99,store95)
  write.csv(store95, paste0("~/Documents/snp_cutoff/output/store_run_check_pat459_",map,".csv"))
  
  
  # vnpatss <- rep(vnpats, 5)
  # vnrunss <- rep(vnruns, each = 5)
  
  # plot(vnpatss, store[,1])
  # points(vnrunss, store[,1], col = "red")
  # 
  # plot(vnpatss,vnrunss,col = rainbow(length(store))[rank(store)])
  
  # S <- cbind(vnpatss, vnrunss, store)
  # 
  # m <- matrix(store, ncol = 5)
  
  # rows m = n runs. cols = n pats 
  ### for cc22_core:  500 runs and 10,000 patients gives consistent results
  # [,1] [,2] [,3] [,4] [,5]
  # [1,]   27   26   19   17   17
  # [2,]   29   24   21   18   18
  # [3,]   29   23   22   18   18
  # [4,]   36   25   23   19   18
  # [5,]   38   27   24   19   18
  
  
  ### for cc22:  500 runs and 10,000 patients gives consistent results
  # [,1] [,2] [,3] [,4] [,5]
  # [1,]   25   21   20   18   18
  # [2,]   29   25   21   19   19
  # [3,]   38   26   22   19   18
  # [4,]   39   26   24   19   19
  # [5,]   37   27   25   19   19  
  
  
}



### Analyse output
store22 <- read.csv(paste0("~/Documents/snp_cutoff/output/store_run_check_pat459_CC22.csv"))[,-1]
store30 <- read.csv(paste0("~/Documents/snp_cutoff/output/store_run_check_pat459_CC30.csv"))[,-1]
store22_c <- read.csv(paste0("~/Documents/snp_cutoff/output/store_run_check_pat459_CC22_core.csv"))[,-1]
store30_c <- read.csv(paste0("~/Documents/snp_cutoff/output/store_run_check_pat459_CC30_core.csv"))[,-1]

vnruns

out <- as.data.frame(cbind(c(store22, store30, store22_c, store30_c), c(rep(22,length(vnruns)),rep(30,length(vnruns)),rep("22 core",length(vnruns)),rep("30 core",length(vnruns)))))
colnames(out) <- c("SNPthreshold", "Map")
out$Nruns <- rep(vnruns,4)

ggplot(out, aes(x=Nruns, y = SNPthreshold)) + geom_count() + facet_wrap(~ Map) + scale_x_continuous("Number of simulations") + scale_y_discrete("SNP threshold") 
ggsave("output/Nruns_threshold.pdf")

# store22 <- read.csv(paste0("~/Dropbox/MRC SD Fellowship/Research/SNP cutoffs/output/store_run_check_CC22.csv"))[,-1]
# store30 <- read.csv(paste0("~/Dropbox/MRC SD Fellowship/Research/SNP cutoffs/output/store_run_check_CC30.csv"))[,-1]
# store22_c <- read.csv(paste0("~/Dropbox/MRC SD Fellowship/Research/SNP cutoffs/output/store_run_check_CC22_core.csv"))[,-1]
# store30_c <- read.csv(paste0("~/Dropbox/MRC SD Fellowship/Research/SNP cutoffs/output/store_run_check_CC30_core.csv"))[,-1]

# vnpats = c(100,500,1000,5000,10000)
# vnruns = c(100,500,1000,5000,10000)
# vnpatss <- rep(vnpats, 5)
# vnrunss <- rep(vnruns, each = 5)
# 
# plot(vnpatss, store22[,1], xlab = "npats",ylab = "store99")
# points(vnrunss, store22[,1], col = "red")
# plot(vnpatss, store30[,1], xlab = "npats",ylab = "store99")
# points(vnrunss, store30[,1], col = "red")
# plot(vnpatss, store22_c[,1], xlab = "npats",ylab = "store99")
# points(vnrunss, store22_c[,1], col = "red")
# plot(vnpatss, store30_c[,1], xlab = "npats",ylab = "store99")
# points(vnrunss, store30_c[,1], col = "red")
# 
# plot(vnpatss, store22[,2], xlab = "npats",ylab = "store99")
# points(vnrunss, store22[,2], col = "red")
# plot(vnpatss, store30[,2], xlab = "npats",ylab = "store99")
# points(vnrunss, store30[,2], col = "red")
# plot(vnpatss, store22_c[,2], xlab = "npats",ylab = "store99")
# points(vnrunss, store22_c[,2], col = "red")
# plot(vnpatss, store30_c[,2], xlab = "npats",ylab = "store99")
# points(vnrunss, store30_c[,2], col = "red")



plot(vnpatss,vnrunss,col = rainbow(length(store))[rank(store)])

S <- cbind(vnpatss, vnrunss, store22)

m22<- matrix(store22[,2], ncol = 5)

# rows m = n runs. cols = n pats 
### for cc22_core:  500 runs and 10,000 patients gives consistent results
# [,1] [,2] [,3] [,4] [,5]
# [1,]   27   26   19   17   17
# [2,]   29   24   21   18   18
# [3,]   29   23   22   18   18
# [4,]   36   25   23   19   18
# [5,]   38   27   24   19   18


### for cc22:  500 runs and 10,000 patients gives consistent results
# [,1] [,2] [,3] [,4] [,5]
# [1,]   25   21   20   18   18
# [2,]   29   25   21   19   19
# [3,]   38   26   22   19   18
# [4,]   39   26   24   19   19
# [5,]   37   27   25   19   19  


# ### Analyse output - fixed patient 
# store22 <- read.csv(paste0("~/Documents/snp_cutoff/output/store_run_check_fix_patCC22_core.csv)[,-1]
# store30 <- read.csv(paste0("~/Dropbox/MRC SD Fellowship/Research/SNP cutoffs/output/store_run_check_CC30.csv"))[,-1]
# store22_c <- read.csv(paste0("~/Dropbox/MRC SD Fellowship/Research/SNP cutoffs/output/store_run_check_CC22_core.csv"))[,-1]
# store30_c <- read.csv(paste0("~/Dropbox/MRC SD Fellowship/Research/SNP cutoffs/output/store_run_check_CC30_core.csv"))[,-1]
# 
# vnpats = c(100,500,1000,5000,10000)
# vnruns = c(100,500,1000,5000,10000)
# vnpatss <- rep(vnpats, 5)
# vnrunss <- rep(vnruns, each = 5)
# 
# plot(vnpatss, store22[,1], xlab = "npats",ylab = "store99")
# points(vnrunss, store22[,1], col = "red")
# plot(vnpatss, store30[,1], xlab = "npats",ylab = "store99")
# points(vnrunss, store30[,1], col = "red")
# plot(vnpatss, store22_c[,1], xlab = "npats",ylab = "store99")
# points(vnrunss, store22_c[,1], col = "red")
# plot(vnpatss, store30_c[,1], xlab = "npats",ylab = "store99")
# points(vnrunss, store30_c[,1], col = "red")
# 
# plot(vnpatss, store22[,2], xlab = "npats",ylab = "store99")
# points(vnrunss, store22[,2], col = "red")
# plot(vnpatss, store30[,2], xlab = "npats",ylab = "store99")
# points(vnrunss, store30[,2], col = "red")
# plot(vnpatss, store22_c[,2], xlab = "npats",ylab = "store99")
# points(vnrunss, store22_c[,2], col = "red")
# plot(vnpatss, store30_c[,2], xlab = "npats",ylab = "store99")
# points(vnrunss, store30_c[,2], col = "red")
# 
# 
# 
# plot(vnpatss,vnrunss,col = rainbow(length(store))[rank(store)])
# 
# S <- cbind(vnpatss, vnrunss, store22)
# 
# m22<- matrix(store22[,2], ncol = 5)
# 
# # rows m = n runs. cols = n pats 
# ### for cc22_core:  500 runs and 10,000 patients gives consistent results
# # [,1] [,2] [,3] [,4] [,5]
# # [1,]   27   26   19   17   17
# # [2,]   29   24   21   18   18
# # [3,]   29   23   22   18   18
# # [4,]   36   25   23   19   18
# # [5,]   38   27   24   19   18
# 
# 
# ### for cc22:  500 runs and 10,000 patients gives consistent results
# # [,1] [,2] [,3] [,4] [,5]
# # [1,]   25   21   20   18   18
# # [2,]   29   25   21   19   19
# # [3,]   38   26   22   19   18
# # [4,]   39   26   24   19   19
# # [5,]   37   27   25   19   19  
# 
# 
# }
# 
# ### Analyse output
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#   