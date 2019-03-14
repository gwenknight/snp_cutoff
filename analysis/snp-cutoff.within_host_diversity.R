

# Intructions on using this R script:
#   1. Change working_dir variable to include the full path of your working directory
#   2. Change full path to SupplementaryData1.xlsx
#   3. The R package gdata needs to be installed, which is used to load .xlsx files
#   4. The R package ggplot2 needs to be installed, which is used for plotting
#   5. The R package lme4 needs to be installed, which is used to run linear mixed models
#   6. The R package svglite needs to be installed, which is used to save plots as SVG files


# The script is structured in X sections where:
#   1. The input files are defined and loaded in INPUT FILES section


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

dataS1 = read.xls(dataS1_file, sheet = 1, header = T)

dim(dataS1)
# [1] 1557   14



##########################################################################################################
###                                     2. PERCENTAGE OF MIXED STRAINS                                ####
##########################################################################################################

# Number of individuals with more than one isolate

length(unique(dataS1$AnonPtNo))
# [1] 459

# Number of individuals with more than one strain (as defined by having isolates from different CCs or
# isolates from the same CC but different clades, labelled as outliers)

length(unique(dataS1$AnonPtNo[which(grepl("outlier",dataS1$Note)==TRUE)]))
# [1] 23

# NOTE: 23/459 (5%) individuals with consecutive isolates had different MRSA strains 

# Removing outliers

dataS1 = dataS1[-which(grepl("outlier",dataS1$Note)==TRUE),]
dim(dataS1)
# [1] 1510   14

##########################################################################################################
###                           3. SNP DISTANCES AMONG ISOLATES COLLECTED ON THE SAME DAY               ####
##########################################################################################################

# Isolates from the same patient collected on the same day will be used to calculate the cloud of diversity

dataS1sd = dataS1[which(dataS1$TimeGap==0),];

# An uneven number of isolates per patient collected on the same day could bias such calculation. Thus, for 
# individuals with more than two isolates, the one with the maximum SNP distance is kept, to keep only one
# isolate pair per patient

keepInd = vector()

individuals = unique(as.vector(dataS1sd$AnonPtNo))

for(i in 1:length(individuals))
{
  tmp = which(dataS1sd$AnonPtNo==individuals[i])
  if(length(tmp)==1)
  {
    keepInd = c(keepInd,tmp)
  } else
  {
    tmp2 = which(dataS1sd$SNPs[tmp]==max(dataS1sd$SNPs[tmp]))
    keepInd = c(keepInd,tmp[tmp2[1]])
  }
}

# Keeping one isolate pair per patient (the one with the maximum SNP distance)

dataS1sd_max = dataS1sd[keepInd,]
dim(dataS1sd_max)
# [1] 82 14


##########################################################################################################
###                               4. EMPIRICAL DISTRIBUTION OF CLOUD OF DIVERSITY                     ####
##########################################################################################################

# The "cloud of diversity" seems to follow an exponential distribution  

quantile(dataS1sd$SNPs)
# 0%  25%  50%  75% 100% 
# 0    1    3    5   82 

quantile(dataS1sd_max$SNPs)
# 0%  25%  50%  75% 100% 
# 0    1    3    5   82 

quantile(dataS1sd_max$SNPs, probs = 0.95)
# 95% 
# 14.95 

quantile(dataS1sd$SNPs, probs = 0.95)
# 95% 
# 14.85 

## Plots

# Empirical distribution of the cloud of diversity across all CCs

plot_width = 6; plot_height = 5;

plot_cloud_of_diversity = function(data, text_x_offset, plot_title)
{
  size_dot = 1; size_axis_lines = 0.3; text_y_offset = 4; font = "Times"; dot_color = "dimgray";
  axis_text_size = 15; axis_title_size = 20; ann_text_size = 5;
  
  co_y = round(as.numeric(quantile(data$SNPs, probs = 0.95)));
  
  co_x = which(data$SNPs <= co_y); co_x = co_x[length(co_x)];
  
  g1 <- ggplot(data, aes(x=seq(1,nrow(data),1), y=SNPs)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
    geom_point(shape = 21, colour = dot_color, fill = dot_color, size = size_dot) + 
    ylim(0, 100) + 
    geom_segment(aes(x= 0, y = co_y, xend = co_x, yend = co_y), linetype="dashed", size=size_axis_lines) + 
    geom_segment(aes(x= co_x, y = 0, xend = co_x, yend = co_y), linetype="dashed", size=size_axis_lines) + 
    annotate("text", x = 0 + text_x_offset, y = co_y + text_y_offset, label = paste("95 percentile =",co_y,"SNPs",sep=" "), family=font, size=ann_text_size) + 
    ylab("Number of SNPs") + 
    xlab("Patients") +
    ggtitle(plot_title) +
    theme(text = element_text(family = font)) +
    theme(axis.text = element_text(size=axis_text_size, color="black"), axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))
  
  return(g1)
}

text_x_offset = 15;

dataS1sd = dataS1sd_max[order(as.numeric(dataS1sd_max$SNPs)),];

plot_title = "Empirical Cloud of Diversity"

g1 = plot_cloud_of_diversity(dataS1sd,text_x_offset,plot_title)

plot_file = "empirical_clould_of_diversity.allCCs.pdf";

ggsave(plot_file, plot = g1, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")



# Empirical distribution of the cloud of diversity for CC22 isolates only

dataS1sdCC22 = dataS1sd_max[which(dataS1sd_max$CC1 == 22),];

dataS1sdCC22 = dataS1sdCC22[order(as.numeric(dataS1sdCC22$SNPs)),];

text_x_offset = 12; 

plot_title = "Empirical Cloud of Diversity (CC22)";

g2 = plot_cloud_of_diversity(dataS1sdCC22,text_x_offset,plot_title)

plot_file = "empirical_clould_of_diversity.CC22.pdf";

ggsave(plot_file, plot = g2, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")




##########################################################################################################
###                                   6. APPLYING LINEAR MIXED MODELS                                 ####
##########################################################################################################

# Linear mixed models are applied to calculate the SNP accumulation rate and to model the "cloud of diversity"

# The number of SNPs between MRSA isolates (SNPs) is modelled as a function of the time gap (TimeGap) between isolates

# The intercept (that is, number of SNPs at time 0) is interpreted as the "cloud of diversity" at time 0 and assumbed to 
# vary by patient (AnonPtNo, random variable)

# Since most of MRSA isolates belong to CC22, the linear mixed model was applied to CC22 isolates only
# as well as to all MRSA CCs too

library(lme4)

dataS1cc22 = dataS1[which(dataS1$CC1==22),];

lmer1_cc22  =  lmer(SNPs  ~ TimeGap   + (1|AnonPtNo), data=dataS1cc22)

summary(lmer1_cc22)

lmer1_all  =  lmer(SNPs  ~ TimeGap   + (1|AnonPtNo), data=dataS1)

summary(lmer1_all)


##########################################################################################################
###                                    5. PLOT SNP DISTANCES OVER TIME                                ####
##########################################################################################################

# Binning data point by TimeGap in months

bins_from = seq(0,330,30); bins_to = seq(30,360,30); bin_month = seq(1,12,1);

dataS1$bin = NA;

for(b in 1:length(bin_month))
{
  tmp = which(dataS1$TimeGap>=bins_from[b] & dataS1$TimeGap<bins_to[b])
  if(length(tmp)>0){ dataS1$bin[tmp] = as.character(bin_month[b]); }
}

# Converting month bin label to factor

dataS1$bin = factor(dataS1$bin,seq(1,12,1))

# Creating X labels

xlab <- paste(levels(dataS1$bin),"\n(N=",table(dataS1$bin),")",sep="")

size_axis_lines = 0.3; axis_x_text_size = 8; axis_y_text_size = 15; axis_title_size = 20;

plot_title = "Number of SNPs over time"

# Boxplot of binned SNP distances per month

boxplot = ggplot(dataS1,aes(x = bin, y = SNPs)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",
                                                                     size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
  geom_boxplot(varwidth = TRUE, notch=FALSE, outlier.size = 0.5) + 
  scale_y_continuous(limits=c(0,25)) +
  scale_x_discrete(labels=xlab) + 
  ylab("Number of SNPs") + 
  xlab("Time distance in months\n(Pairwise Comparisons)") +
  ggtitle(plot_title) +
  theme(text = element_text(family = font)) +
  theme(axis.text.x = element_text(size=axis_x_text_size, color="black"), axis.text.y = element_text(size=axis_y_text_size, color="black"), 
        axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))


boxplot_file = "snp_distances_over_time.allCCs.svg";

ggsave(boxplot_file, plot = boxplot, device = "svg", width = 6, height = 5, dpi = 300, units = "in")


# Boxplot of binned SNP distances per month (CC22 only)

dataS1cc22 = dataS1[which(dataS1$CC1==22),];

xlab <- paste(levels(dataS1cc22$bin),"\n(N=",table(dataS1cc22$bin),")",sep="")

boxplot = ggplot(dataS1cc22,aes(x = bin, y = SNPs)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black",
                                                                     size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines)) +
  geom_boxplot(varwidth = TRUE, notch=FALSE, outlier.size = 0.5) + 
  scale_y_continuous(limits=c(0,25)) +
  scale_x_discrete(labels=xlab) + 
  ylab("Number of SNPs") + 
  xlab("Time distance in months\n(Pairwise Comparisons)") +
  ggtitle(plot_title) +
  theme(text = element_text(family = font)) +
  theme(axis.text.x = element_text(size=axis_x_text_size, color="black"), axis.text.y = element_text(size=axis_y_text_size, color="black"), 
        axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))


boxplot_file = "snp_distances_over_time.CC22.svg";

ggsave(boxplot_file, plot = boxplot, device = "svg", width = 6, height = 5, dpi = 300, units = "in")





##########################################################################################################
###                                 7. CALCULATION OF THE MRSA SUBSTITUTION RATE                      ####
##########################################################################################################

# Using CC22 MRSA isolates only

timegap_coefficient_cc22 = coef(summary(lmer1_cc22))[2,1]
# [1] 0.01335459

# As the independent variable (TimeGap) unit is days, the units of the coefficient are SNPs per genome per day

substitution_rate_cc22 = timegap_coefficient_cc22*365
substitution_rate_cc22
# [1] 4.874424

# The units of substitution rate are now SNPs per genome per year


# Using MRSA isolates from all CCs

timegap_coefficient_all = coef(summary(lmer1_all))[2,1]

substitution_rate_all = timegap_coefficient_all*365

substitution_rate_all
# [1] 4.091438



##########################################################################################################
###                                 8. MODELED DISTRIBUTION OF CLOUD OF DIVERSITY                     ####
##########################################################################################################

# Extracting varying intercepts for all patients

beta0_cc22 = as.vector(unlist(coef(lmer1_cc22)$AnonPtNo["(Intercept)"]));

quantile(beta0_cc22)
# 0%        25%        50%        75%       100% 
# -0.9801063  1.4102880  3.1013106  6.3307336 51.2781168 

quantile(beta0_cc22,prob=0.95)
# 95% 
# 18.09833 

beta0_all = as.vector(unlist(coef(lmer1_all)$AnonPtNo["(Intercept)"]));

quantile(beta0_all)
# 0%        25%        50%        75%       100% 
# -0.5326214  2.2241942  4.2607700  7.8025847 56.8797944 

quantile(beta0_all,prob=0.95)
# 95% 
# 19.43245 

## Plots

# Modelled distribution of the cloud of diversity across all CCs

beta0_all_df = as.data.frame(cbind(seq(1,length(beta0_all),1),sort(beta0_all)))
colnames(beta0_all_df) = c("X","SNPs")

text_x_offset = 60;

plot_title = "Modelled Cloud of Diversity";

g3 = plot_cloud_of_diversity(beta0_all_df,text_x_offset,plot_title)

plot_file = "modelled_clould_of_diversity.allCCs.pdf";

ggsave(plot_file, plot = g3, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")

# Modelled distribution of the cloud of diversity for CC22 isolates only

beta0_cc22_df = as.data.frame(cbind(seq(1,length(beta0_cc22),1),sort(beta0_cc22)))
colnames(beta0_cc22_df) = c("X","SNPs")

text_x_offset = 55;

plot_title = "Modelled Cloud of Diversity (CC22)";

g4 = plot_cloud_of_diversity(beta0_cc22_df,text_x_offset,plot_title)

plot_file = "modelled_clould_of_diversity.CC22.pdf";

ggsave(plot_file, plot = g4, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")


##########################################################################################################
####****                          9. SIMULATION MODEL OF TRANSMISSION                             ****####
##########################################################################################################
# Load in the required functions
# These are the distribution sampler / the curve for the number of SNPs over time / the simulation model
source("analysis/simu_transmission_fn.R")
require(reshape2)

## Simulation population 
npat = 459 # number of patient samples in cohort 1
nruns = 200000 # number of transmission samples - can be increased to reduce variation in number of SNPs

ndays = 180 # time between samples

## Choose which mapping: All of CC22/CC30 or core genome only
#map <- "CC22" # ST22 strain HO 5096 0412 mapping data 
#map <- "CC30"
#map <- "CC22_core"
map <- "CC30_core"

## mu = substitution_rate
if(map == "CC22"){mu = 4 / 365}
if(map == "CC30"){mu = 5 / 365}
if(map == "CC22_core"){mu = 3 / 365}
if(map == "CC30_core"){mu = 3 / 365}

## Parameters from model fit
param_general_fit <- read.csv(paste0("output/param_general_fit_",map,".csv"))[,2]

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

## Run the simulation 10 times to give a range on the maximum, some random variation expected, dependening on the number of runs used. 
m <- rep(0,10)
for(i in 1:10){
  ss <- simu_runs(ndays,mu,npat,nruns, param_general_fit, t0prob_dist)
  
  ## Maximum number of SNPs needed to capture 95% or 99% of the transmission events
  m[i] <- max(ss$store_limits[which(ss$store_limits$variable == "95%"),"value"]) # = below this, capture 95% of all transmissions
  
}


## Plot the output from the last run
# g5 = ggplot(ss$store_limits, aes(x=value, fill = variable)) + geom_histogram(aes(y=..density..), binwidth = 1, position = "identity") +   
#   facet_wrap(~variable) + guides(fill=FALSE) + scale_y_continuous(paste0("Density across ", nruns, " simulations")) + scale_x_continuous("Number of SNPs")
# 
# plot_file = paste0("simulation_model_distribution_of_SNPs_",map,".pdf");
# 
# ggsave(plot_file, plot = g5, device = "pdf", width = plot_width, height = plot_height, dpi = 300, units = "in")


## Results: 
max(m)
range(m)

##            95% of transmission events (max [range])
## map = CC22      = 15 (14 - 15)
## map = CC30      = 51 (48 - 51)
## map = CC22_core = 13 (12 - 13)
## map = CC30_core = 13 (12 - 13)

