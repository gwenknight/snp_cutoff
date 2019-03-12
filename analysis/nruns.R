# how many simulations? 

library(reshape2)
setwd("~/Documents/snp_cutoff/")
source("analysis/simu_transmission_fn.R")

# Choose which mapping 
for(kk in 2:4){
  print(c("kk = ", kk))
  if(kk == 1){map <- "CC22"} # ST22 strain HO 5096 0412 mapping data 
  if(kk == 2){map <- "CC30"}
  if(kk == 3){map <- "CC22_core"}
  if(kk == 4){map <- "CC30_core"}
  
  # Parameters from model fit
  param_general_fit <- read.csv(paste0("~/Documents/snp_cutoff/output/param_general_fit_",map,".csv"))[,2]
  
  # Simulation population 
  #vnpats = c(100,500,1000,5000,10000)
  vnruns = c(100000,100000,100000,100000,100000,
             200000,200000,200000,200000,200000)
  
  store99 <- c()
  store95 <- c()
  
  #for(i in 1:length(vnpats)){
    for(j in 1:length(vnruns)){
      print(c(i,j))
      npat = 500 #vnpats[i]
      nruns = vnruns[j]
      
      ndays = 180 # time between samples
      
      mu = 4.874424 / 365 #per day
      
      # Run the simulation 
      ss <- simu_runs(ndays,mu,npat,nruns, param_general_fit)
      
      # Maximum number of SNPs needed to capture 99% of the transmission events
      store99 <- c(store99, max(ss$store_limits[which(ss$store_limits$variable == "99%"),"value"]))
      store95 <- c(store95, max(ss$store_limits[which(ss$store_limits$variable == "95%"),"value"]))
    }
  
  
  store <- cbind(store99,store95)
  write.csv(store, paste0("~/Documents/snp_cutoff/output/store_run_check_pat500_",map,".csv"))
  
  vnpatss <- rep(vnpats, 5)
  vnrunss <- rep(vnruns, each = 5)
  
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
  