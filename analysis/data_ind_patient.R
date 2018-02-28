### Data analysis of CC22 MRSA 

library(ggplot2)
library(fitdistrplus)

# Grab data
dd <- read.csv("~/Documents/snp_cutoff/data/mrsa_pros.csv")
ggplot(dd, aes(x = liDayGap,y=liSNPdis)) + geom_point(alpha=0.05) 
#ggplot(dd, aes(x = liDayGap,y=liSNPdis)) + geom_point(alpha=0.05) + scale_x_continuous(limits = c(0,30)) + scale_y_continuous(limits = c(0,50))

ggplot(dd, aes(x = liDayGap)) + geom_histogram(binwidth=1) 
table(d)

###**** Same day samples ***########################################################################################################################################################################################################
w<-which(dd$liDayGap == 0)
length(w) # 85 pairs taken on same day
length(unique(dd[w,"lsPt"])) # from 69 patients (matches Francesc ppt)
dd_w <- dd[w,]

# if choose "main" strain then also no SNP difference so
dls <- dd_w$liSNPdis
dls <- c(dls, matrix(0,1,length(dls)))

h <- hist(dls,breaks=seq(0,175,1))
setwd("data")
write.csv(h$density, "same_day_snp_dist.csv")
h_cum <- cumsum(h$density)
plot(seq(1,175,1),h_cum)
length(which(h_cum < 0.7)) # 70% 0 SNP differences
length(which(h_cum < 0.9)) # 90% 4 SNPS different
length(which(h_cum < 0.95)) # 95% 12 SNPS different 

ggplot(dd_w, aes(x = liSNPdis)) + geom_histogram(breaks=seq(0,90,1)) + scale_x_continuous("SNP distance") + 
  scale_y_continuous("# of isolates")
ggsave("Same_day_SNP_distn.pdf")

#ggplot(dd_w, aes(x = liSNPdis)) + geom_histogram(binwidth = 1) + facet_wrap(~lsPt)


## Do we want to take into account that most people had v few differences? 
# no. think we want to treat each pair as an independent sample. 
#for(i in 1:length(unique(dd_w$lsPt]))){
#  ww <- which(dd_w$lsPt == i)
#  ll <- length(ww)
#}

###**** Other day samples ***########################################################################################################################################################################################################
# Try and capture this distribution using mutation rate? and initial within host distribution? 
w<-which(dd$liDayGap > 0)
w<-which(dd$liDayGap == 4)
length(w) # 1096 pairs not on same day
length(unique(dd[w,"lsPt"])) # from 322 patients (less than Francesc ppt?)
dd_wo <- dd[w,]

ho <- hist(dd_wo$liSNPdis,breaks=seq(0,max(dd_wo$liSNPdis),1))
h_cumo <- cumsum(ho$density)
plot(seq(1,175,1),h_cumo[1:175])
length(which(h_cumo < 0.6)) # < 30% chance get 0 SNPs
length(which(h_cumo < 0.95)) # 12 SNPS for 95%

setwd("data")
write.csv(dd_wo$liDayGap, "sample_time_dist.csv")

ggplot(dd_wo, aes(x = liSNPdis)) + geom_histogram(binwidth = 1) 
ggplot(dd_wo, aes(x = liDayGap,y=liSNPdis)) + geom_point(alpha=0.05) 
ggplot(dd_wo[1:100,], aes(x = liSNPdis)) + geom_histogram(binwidth = 1) + facet_wrap(~lsPt)

### COMPARE
ggplot(dd_wo, aes(x = liSNPdis)) + geom_histogram(binwidth = 1) + 
  geom_histogram(data = dd_w, aes(x = liSNPdis, fill = "red"), binwidth = 1) +
  scale_x_continuous(limits = c(0,50))

###**** Fit distributions to the data ***########################################################################################################################################################################################################
f1<-fitdist(h$density,"exp",method="mle")
f2<-fitdist(h$density,"gamma",method="mme")

summary(f1)
summary(f2)

b_best<-bootdist(f1)
print(f1)
plot(f1)
summary(f1)

###**** Between patient samples ***###################################################################################################################
ddt <- read.csv("~/Documents/snp_cutoff/data/mrsa_pros.first_isolate.snp_dis.csv")

ddt22 <- ddt[intersect(which(ddt$CC1 == 22),which(ddt$CC2 == 22)),]
w<-which(ddt22$AnonPtNo1 == ddt22$AnonPtNo2) # by row
length(w) # no same patient comparisons

h <- hist(ddt22$SNPs,breaks=seq(0,725,1))
setwd("../data")
write.csv(h$density, "22_snp_dist.csv")

h_cum <- cumsum(h$density)
plot(seq(1,725,1),h_cum)
length(which(h_cum < 0.9)) # 193 SNPS for 90%
length(which(h_cum < 0.95)) # 444 SNPS for 95%

