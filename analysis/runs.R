setwd("~/Documents/snp_cutoff/")

# Simulation population 
npat = 459 # number of patient samples in cohort 1 - matches above analysis
ndays = 180 # time between samples



# Choose which mapping 
#map <- "CC22" # ST22 strain HO 5096 0412 mapping data 
#map <- "CC30"
#map <- "CC22_core"
map <- "CC30_core"

#### CC22  ### NEEDS TO CHANGE FOR EACH MAP: FRANCESC HOW ARE YOU DOING IN THE ABOVE? 
#mu = substitution_rate_cc22/365 # mutation rate for CC22: 4.874424 / 365 per day
if(map == "CC22"){mu = 4 / 365}
if(map == "CC30"){mu = 5 / 365}
if(map == "CC22_core"){mu = 3 / 365}
if(map == "CC30_core"){mu = 3 / 365}


# Parameters from model fit
param_general_fit <- read.csv(paste0("output/param_general_fit_",map,".csv"))[,2]

nruns_v <- c(rep(1e4, 10), rep(5e4, 10), rep(1e5, 10), rep(2e5, 10), rep(5e5, 10))
#nruns_v <- c(rep(5e4, 100))

m <- c()

for(i in 1:length(nruns_v)){
  print(i)
  nruns <- nruns_v[i]
  # Run the simulation 
  ss <- simu_runs(ndays,mu,npat,nruns, param_general_fit)
  
  # Maximum number of SNPs needed to capture 95% or 99% of the transmission events
  m <- rbind(m,c(nruns, max(ss$store_limits[which(ss$store_limits$variable == "95%"),"value"]), # = below this, capture 95% of all transmissions
  max(ss$store_limits[which(ss$store_limits$variable == "99%"),"value"]))) # = 
  
}
write.csv(m, paste0("1903",map,".csv"))

# Choose which mapping 
map <- "CC22" # ST22 strain HO 5096 0412 mapping data 
#map <- "CC30"
#map <- "CC22_core"
#map <- "CC30_core"

#### CC22  ### NEEDS TO CHANGE FOR EACH MAP: FRANCESC HOW ARE YOU DOING IN THE ABOVE? 
#mu = substitution_rate_cc22/365 # mutation rate for CC22: 4.874424 / 365 per day
if(map == "CC22"){mu = 4 / 365}
if(map == "CC30"){mu = 5 / 365}
if(map == "CC22_core"){mu = 3 / 365}
if(map == "CC30_core"){mu = 3 / 365}


# Parameters from model fit
param_general_fit <- read.csv(paste0("output/param_general_fit_",map,".csv"))[,2]


m <- c()

for(i in 1:length(nruns_v)){
  print(i)
  nruns <- nruns_v[i]
  # Run the simulation 
  ss <- simu_runs(ndays,mu,npat,nruns, param_general_fit)
  
  # Maximum number of SNPs needed to capture 95% or 99% of the transmission events
  m <- rbind(m,c(nruns, max(ss$store_limits[which(ss$store_limits$variable == "95%"),"value"]), # = below this, capture 95% of all transmissions
                 max(ss$store_limits[which(ss$store_limits$variable == "99%"),"value"]))) # = 
  
}
write.csv(m, paste0("1903",map,".csv"))

# Choose which mapping 
#map <- "CC22" # ST22 strain HO 5096 0412 mapping data 
map <- "CC30"
#map <- "CC22_core"
#map <- "CC30_core"

#### CC22  ### NEEDS TO CHANGE FOR EACH MAP: FRANCESC HOW ARE YOU DOING IN THE ABOVE? 
#mu = substitution_rate_cc22/365 # mutation rate for CC22: 4.874424 / 365 per day
if(map == "CC22"){mu = 4 / 365}
if(map == "CC30"){mu = 5 / 365}
if(map == "CC22_core"){mu = 3 / 365}
if(map == "CC30_core"){mu = 3 / 365}


# Parameters from model fit
param_general_fit <- read.csv(paste0("output/param_general_fit_",map,".csv"))[,2]


m <- c()

for(i in 1:length(nruns_v)){
  print(i)
  nruns <- nruns_v[i]
  # Run the simulation 
  ss <- simu_runs(ndays,mu,npat,nruns, param_general_fit)
  
  # Maximum number of SNPs needed to capture 95% or 99% of the transmission events
  m <- rbind(m,c(nruns, max(ss$store_limits[which(ss$store_limits$variable == "95%"),"value"]), # = below this, capture 95% of all transmissions
                 max(ss$store_limits[which(ss$store_limits$variable == "99%"),"value"]))) # = 
  
}
write.csv(m, paste0("1903",map,".csv"))

# Choose which mapping 
#map <- "CC22" # ST22 strain HO 5096 0412 mapping data 
#map <- "CC30"
map <- "CC22_core"
#map <- "CC30_core"

#### CC22  ### NEEDS TO CHANGE FOR EACH MAP: FRANCESC HOW ARE YOU DOING IN THE ABOVE? 
#mu = substitution_rate_cc22/365 # mutation rate for CC22: 4.874424 / 365 per day
if(map == "CC22"){mu = 4 / 365}
if(map == "CC30"){mu = 5 / 365}
if(map == "CC22_core"){mu = 3 / 365}
if(map == "CC30_core"){mu = 3 / 365}


# Parameters from model fit
param_general_fit <- read.csv(paste0("output/param_general_fit_",map,".csv"))[,2]


m <- c()

for(i in 1:length(nruns_v)){
  print(i)
  nruns <- nruns_v[i]
  # Run the simulation 
  ss <- simu_runs(ndays,mu,npat,nruns, param_general_fit)
  
  # Maximum number of SNPs needed to capture 95% or 99% of the transmission events
  m <- rbind(m,c(nruns, max(ss$store_limits[which(ss$store_limits$variable == "95%"),"value"]), # = below this, capture 95% of all transmissions
                 max(ss$store_limits[which(ss$store_limits$variable == "99%"),"value"]))) # = 
  
}
write.csv(m, paste0("1903",map,".csv"))


#### Explore
setwd("~/Documents/snp_cutoff/")
mcc22 <- read.csv(paste0("1903","CC22",".csv"))[,-1]
mcc30 <- read.csv(paste0("1903","CC30",".csv"))[,-1]
mcc22c <- read.csv(paste0("1903","CC22_core",".csv"))[,-1]
mcc30c <- read.csv(paste0("1903","CC30_core",".csv"))[,-1]


mm <- as.data.frame(rbind(mcc22,mcc30,mcc22c,mcc30c))
mm$names <- c(rep(22,dim(mcc22)[1]),rep(30,dim(mcc30)[1]),rep(220,dim(mcc22c)[1]),rep(300,dim(mcc30c)[1]))

colnames(mm) <- c("nrun","95","99","names")
mmm <- melt(mm[,c("nrun","95","names")], id.vars = c("nrun","names"))

ggplot(mmm, aes(x=nrun, y = value)) + geom_point(aes(col = variable)) + 
  facet_wrap(~names, scales = 'free') + scale_x_continuous("Number of simulations") + scale_y_continuous("Predicted SNP threshold")
ggsave("output/1903_nrun.pdf")

#### Implies 95% more stable - 99% subject to stochastic fluctuation in final transmissions
# Need 50,000 to get stable top limit
max(mcc22[,2])
max(mcc30[,2])
max(mcc22c[,2])
max(mcc30c[,2])

