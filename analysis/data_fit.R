### Code to fit exponential curve to number of SNPs distribution by time

##########################################################################################################
###                                         1. INPUT FILES                                            ####
##########################################################################################################

#working_dir = "/Users/francesccoll/fellowship/2.projects/19.mrsa_hicf_project/within_host/"
working_dir = "~/Documents/snp_cutoff/"

setwd(working_dir)

#dataS1_file = "/Users/francesccoll/fellowship/2.projects/19.mrsa_hicf_project/manuscript/supplementary_data/SupplementaryData1.xlsx"
dataS1_file = "data/SupplementaryData1.xlsx"

require(gdata)
require(ggplot2)
require(svglite)

##########################################################################################################
###                                   2. ST22 strain HO 5096 0412 mapping data                        ####
##########################################################################################################

dataS1 = read.xls(dataS1_file, sheet = 1, header = T)

dim(dataS1)
# [1] 1557   14

# Removing outliers
dataS1 = dataS1[-which(grepl("outlier",dataS1$Note)==TRUE),]
dim(dataS1)

# Only CC22 data
dd = dataS1[which(dataS1$CC1 == 22),];

##########################################################################################################
###                                         2a. MODEL FIT                                              ####
##########################################################################################################

u<-sort(unique(dd$TimeGap))
coef_store<-c(); data_store<-c()

# Run through every time gap. If more than 10 data points then fit exponential to the data and store coefficients
for(i in 1:length(u)){
  # For this time gap grab the specific data
  w<-which(dd$TimeGap == u[i])
  if(length(w)>10){ # need > 10? for fit to be ok
    print(i)
    # Histogram of distribution 
    h <- hist(dd[w,]$SNPs,breaks=seq(0,max(dd[w,"SNPs"]),1),plot = T)
    data <- data.frame( x = h$breaks[-length(h$breaks)], y = h$counts)
    # Fit exponential to distribution of SNPs at each time point
    mod <- nls(y ~ exp(a + b * x), data = data, start = list(a = log(1), b = -1),
               control = list(maxiter = 500))
    # Can plot as go or all together below
    #plot(data$x,data$y)
    #lines(data$x, predict(mod, list(x = data$x)))
    
    # store coefficients
    data_store <- rbind(data_store,cbind(u[i],data,predict(mod, list(x = data$x))))
    coef_store <- rbind(coef_store, c(u[i],coef(mod)))
    }
}

colnames(data_store) <- c("day","x","y","m")

coef_store <- as.data.frame(coef_store)
colnames(coef_store) <- c("day","a","b")

ggplot(data_store,aes(x=x, y=y)) + geom_point() + geom_line(aes(x=x,y=m),col="red") + 
  facet_wrap(~day,scales="free")
ggsave("output/fit_exponential_per_day.pdf",width = 12, height = 12)

### What is the trend in the coefficients? 
plot(coef_store$day,coef_store$a,ylim = c(0,4))
# Assuming expoential link between day and coefficient 
mod_a <- nls(a ~ exp(aa + bb * day), data = coef_store, start = list(aa = 1.4, bb = -1),
           control = list(maxiter = 500))
lines(coef_store$day, predict(mod_a, list(x = coef_store$day)),col="red")
mod_a_mod <- coef(mod_a)

pdf("output/a_exponential.pdf")
plot(coef_store$day,coef_store$a,ylim = c(0,4),xlab = "Time between samples",ylab="Intercept") 
lines(coef_store$day, predict(mod_a, list(x = coef_store$day)),col="red")
dev.off()

pdf("output/b_flat.pdf")
plot(coef_store[,1],coef_store[,3],xlab = "Time between samples",ylab = "Slope")
# Assume all v similar - take mean without the v small outliers greater than -8 
mean(coef_store[,3])
p_aa <- mean(coef_store[which(coef_store[,3] > (-4)),3]) # remove outliers
abline(h = p_aa,col="red",lty="dashed")
dev.off()

## General model for the fit
data_store$mod_general <- exp(exp(mod_a_mod[1] + mod_a_mod[2]*data_store[,"day"]) + p_aa * data_store[,"x"])
  
# plot the general model for the fit (blue) against individual (red)
ggplot(data_store,aes(x=x, y=y)) + geom_point() + geom_line(aes(x=x,y=m),col="red") + 
  facet_wrap(~day,scales="free") + geom_line(aes(x=x,y=mod_general),col="blue") + 
   scale_x_continuous("SNP distance") + scale_y_continuous("Count")
ggsave("output/fit_exponential_gen&ind_per_day.pdf",width = 12, height = 12)

# Parameters
mod_a_mod[1]
mod_a_mod[2]
p_aa <- mean(coef_store[which(coef_store[,3] > (-4)),3]) # remove outliers


##########################################################################################################
###                                         2b. MODEL TO GENERATE FIT                                  ####
##########################################################################################################
# Also in simu_transmission_fn.R

## General model 
gen_mod_day <- function(day,l_x){
  ### day = number of days from first sample
  ### l_x = maximum number of SNPs
  x <- seq(0,l_x,1)
  aa_g <-  1.016402 #mod_a_mod[1]
  bb_g <- -0.022452 #mod_a_mod[2]
  p_aa <- -0.976289 # mean 
  out <- exp(exp(aa_g + bb_g*day) + p_aa * x)
  return(out)
}



##########################################################################################################
###                                   3. ST30 strain MRSA252 mapping data                               ####
##########################################################################################################

dataS1 = read.xls(dataS1_file, sheet = 2, header = T)

dim(dataS1)
# [1] 1557   12

# Removing outliers
dataS1 = dataS1[-which(grepl("outlier",dataS1$Note)==TRUE),]
dim(dataS1) # 1510 12

# Only CC22 data
dd = dataS1[which(dataS1$CC1 == 22),];

##########################################################################################################
###                                         3a. MODEL FIT                                              ####
##########################################################################################################

u<-sort(unique(dd$TimeGap))
coef_store<-c(); data_store<-c()

# Run through every time gap. If more than 10 data points then fit exponential to the data and store coefficients
for(i in 1:length(u)){
  # For this time gap grab the specific data
  w<-which(dd$TimeGap == u[i])
  if(length(w)>10){ # need > 10? for fit to be ok
    print(i)
    # Histogram of distribution 
    h <- hist(dd[w,]$SNPs,breaks=seq(0,max(dd[w,"SNPs"]),1),plot = T)
    data <- data.frame( x = h$breaks[-length(h$breaks)], y = h$counts)
    # Fit exponential to distribution of SNPs at each time point
    mod <- nls(y ~ exp(a + b * x), data = data, start = list(a = 1.5, b = -0.2), # Needed different starting values
               control = list(maxiter = 500))
    # Can plot as go or all together below
    #plot(data$x,data$y)
    #lines(data$x, predict(mod, list(x = data$x)))
    
    # store coefficients
    data_store <- rbind(data_store,cbind(u[i],data,predict(mod, list(x = data$x))))
    coef_store <- rbind(coef_store, c(u[i],coef(mod)))
  }
}

colnames(data_store) <- c("day","x","y","m")

coef_store <- as.data.frame(coef_store)
colnames(coef_store) <- c("day","a","b")

ggplot(data_store,aes(x=x, y=y)) + geom_point() + geom_line(aes(x=x,y=m),col="red") + 
  facet_wrap(~day,scales="free")
ggsave("output/fit_exponential_per_day_st30.pdf",width = 12, height = 12)

### What is the trend in the coefficients? 
plot(coef_store$day,coef_store$a,ylim = c(0,4))
# Assuming expoential link between day and coefficient 
mod_a <- nls(a ~ exp(aa + bb * day), data = coef_store, start = list(aa = 1.4, bb = -1),
             control = list(maxiter = 500))
lines(coef_store$day, predict(mod_a, list(x = coef_store$day)),col="red")
mod_a_mod <- coef(mod_a)

pdf("output/a_exponential_st30.pdf")
plot(coef_store$day,coef_store$a,ylim = c(0,4),xlab = "Time between samples",ylab="Intercept") 
lines(coef_store$day, predict(mod_a, list(x = coef_store$day)),col="red")
dev.off()

# Flat line? 
pdf("output/b_flat_st30.pdf")
plot(coef_store[,1],coef_store[,3],xlab = "Time between samples",ylab = "Slope")
# Assume all v similar - take mean without the v small outliers greater than -8 
mean(coef_store[,3])
p_aa <- mean(coef_store[,3]) 
abline(h = p_aa,col="red",lty="dashed")
dev.off()

# No longer flat line for b: time increases the distance: BUT
ll <- lm(coef_store[,3] ~ coef_store[,1])

pdf("output/b_line_st30.pdf")
plot(coef_store[,1],coef_store[,3],xlab = "Time between samples",ylab = "Slope")
# Assume all v similar - take mean without the v small outliers greater than -8 
abline(lm(coef_store[,3] ~ coef_store[,1]), col = "red", lty = "dashed")
dev.off()

## General model for the fit
data_store$mod_general <- exp(exp(mod_a_mod[1] + mod_a_mod[2]*data_store[,"day"]) + 
                          (ll$coefficients[1] + ll$coefficients[2]*data_store[,"day"]) * data_store[,"x"])

# plot the general model for the fit (blue) against individual (red)
ggplot(data_store,aes(x=x, y=y)) + geom_point() + geom_line(aes(x=x,y=m),col="red") + 
  facet_wrap(~day,scales="free") + geom_line(aes(x=x,y=mod_general),col="blue") + 
  scale_x_continuous("SNP distance") + scale_y_continuous("Count")
ggsave("output/fit_exponential_gen&ind_per_day_st30.pdf",width = 12, height = 12)

# Parameters
mod_a_mod[1]
mod_a_mod[2]
ll$coefficients[1]
ll$coefficients[2]

##########################################################################################################
###                                         3b. MODEL TO GENERATE FIT                                  ####
##########################################################################################################
# Also in simu_transmission_fn.R

## General model 
gen_mod_day_st30 <- function(day,l_x){
  ### day = number of days from first sample
  ### l_x = maximum number of SNPs
  x <- seq(0,l_x,1)
  aa_g <-  1.016402 #mod_a_mod[1]
  bb_g <- -0.022452 #mod_a_mod[2]
  p_aa <- -0.127701 #ll$coefficients[1]
  p_bb <-  0.002448 #ll$coefficients[1]
  out <- exp(exp(aa_g + bb_g*day) + (p_aa + p_bb*day) * x)
  return(out)
}

