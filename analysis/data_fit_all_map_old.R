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
#### CHOOSE HERE - which mapping

#map <- "CC22" # ST22 strain HO 5096 0412 mapping data 
#map <- "CC30"
#map <- "CC22_core"
map <- "CC30_core"

##########################################################################################################

## Which sheet? 
if(map == "CC22"){sheet_num = 1}
if(map == "CC30"){sheet_num = 2}
if(map == "CC22_core"){sheet_num = 3}
if(map == "CC30_core"){sheet_num = 4}

dataS1 = read.xls(dataS1_file, sheet = sheet_num, header = T) 

dim(dataS1)
# [1] 1557   14

# Removing outliers
dataS1 = dataS1[-which(grepl("outlier",dataS1$Note)==TRUE),]
dim(dataS1)

# Only CC22 data
dd = dataS1[which(dataS1$CC1 == 22),];

##########################################################################################################
###                                         1a. MODEL FIT                                              ####
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
    mod <- nls(y ~ exp(a + b * x), data = data, start = list(a = 1.5, b = -0.2), # works for CC22/CC30
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
ggsave(paste0("output/fit_exponential_per_day_",map,".pdf"),width = 12, height = 12)

### What is the trend in the coefficients? 
## for a? 
plot(coef_store$day,coef_store$a,ylim = c(0,4))
# Assuming expoential link between day and coefficient 
mod_a <- nls(a ~ exp(aa + bb * day), data = coef_store, start = list(aa = 1.4, bb = -0.2),
           control = list(maxiter = 500))
lines(coef_store$day, predict(mod_a, list(x = coef_store$day)),col="red")
mod_a_mod <- coef(mod_a)

pdf(paste0("output/a_exponential_",map,".pdf"))
plot(coef_store$day,coef_store$a,ylim = c(0,4),xlab = "Time between samples",ylab="Intercept") 
lines(coef_store$day, predict(mod_a, list(x = coef_store$day)),col="red")
dev.off()


## What is the trend for b? 
plot(coef_store[,1],coef_store[,3],xlab = "Time between samples",ylab = "Slope")
# Assume linear in the non-outliers
w<-which(coef_store[,3] > (-4))
lb <- lm(coef_store[w,3] ~ coef_store[w,1])
pred_b <- coefficients(lb)[1] + coefficients(lb)[2]*coef_store[,1]

pdf(paste0("output/b_flat_",map,".pdf"))
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
ggsave(paste0("output/fit_exponential_gen&ind_per_day_",map,".pdf"),width = 12, height = 12)

# Parameters - output
mod_a_mod[1]
mod_a_mod[2]
p_aa 

param_general_fit = c(mod_a_mod[1],mod_a_mod[2],p_aa)

write.csv(param_general_fit,paste0("output/param_general_fit_",map,".csv"))

##########################################################################################################
###                                         1b. MODEL TO GENERATE FIT                                 ####
##########################################################################################################
# Also in simu_transmission_fn.R

## General model 
gen_mod_day <- function(day,l_x, param_general_fit){
  ### day = number of days from first sample
  ### l_x = maximum number of SNPs
  ### param_general_fit = parameters from fit to data
  
  x <- seq(0,l_x,1)
  aa_g <- param_general_fit[1]
  bb_g <- param_general_fit[2]
  p_aa <- param_general_fit[3]
  out <- exp(exp(aa_g + bb_g*day) + p_aa * x)
  return(out)
}

### Prediction for 6months = 180 days
day = 180
x = seq(0,10,1)

gen_pred_180 <- exp(exp(mod_a_mod[1] + mod_a_mod[2]*180) + p_aa * x)
plot(x, gen_pred_180)

ggplot(data_store,aes(x=x, y=y, group = day)) + geom_point(aes(col = day)) + 
  geom_line(aes(x=x,y=mod_general, group = day, col = day)) + 
  scale_x_continuous("SNP distance") + scale_y_continuous("Count")


