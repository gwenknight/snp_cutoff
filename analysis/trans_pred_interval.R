# De Silva Transmission Prediction Interval 
library(ggplot2)
theme_set(theme_bw())

####**** FOR CC22 MRSA DATA ****************************************************************************************************************************************************************#######
# s(t) (number of SNP differences between two samples taken a time t apart).
# s ~ Pois(mu*(t+2u))

# Mu can be estimated from genomic data / literature [De Silva ref gets it from time scaled phylogenies]
# 5 SNPs per genome per year [Coll from this data]

mu <- 5

# u can be estimated as the time to MRCA of isolates taken from within the
# same host at the same time (approximately). 
# s0 = SNP diversity in a single host = 2*mu*u
# [s0 is known in De Silva ref from two datasets, as is mu]. 
# So u can be estimated (u = s0 / (2mu)).

s0 = 4.25 # calculated from mrsa_pros data in data_ind_patient.R

u = s0 / (2*mu) # time to MRCA

# s(t) can then be sampled from Pois(mu*(t+2u)) to generate 99% prediction intervals at each time t.

### Calculate exactly from distribution 
# z score at 99% = 2.576

vec_t <- seq(0,11,0.1)
s <- as.data.frame(matrix(0,length(vec_t),2))
sm <- as.data.frame(matrix(0,length(vec_t),1))
s$t <- vec_t
sm$t <- vec_t

for(i in 1:length(vec_t)){
  t = vec_t[i]
  lambda = mu * (t + 2*u)
  
  s[i,1:2] <- c(lambda + 2.576 * sqrt(lambda), lambda - 2.576 * sqrt(lambda))
  sm[i,1] <- lambda
  
  ## Could calculate by sampling
  #   rand_pois_sample <- rpois(100000,lambda)  
  #   q <- quantile(rand_pois_sample,probs = c(1,99)/100)
  #   qm <- mean(rand_pois_sample)
  #   s[i,1:2] <- q
  #   sm[i,1] <- qm
  
}

colnames(s)<- c("lo","hi","t")

# to generate polygon
sp1 <- as.data.frame(s[,c("lo","t")])
colnames(sp1)<- c("y","x")
sp2<-as.data.frame(cbind(rev(s$hi),rev(s$t)))
colnames(sp2)<- c("y","x")
sp <- rbind(sp1,sp2)
sp[which(sp$y < 0),"y"] <- 0 # 0 lowest SNP difference

# round to give actual SNP difference
sp$y <- floor(sp$y)

sm <- as.data.frame(sm)
colnames(sm) <- c("mean","t")

setwd("../output")

ggplot(sp, aes(x = x, y = y)) + geom_polygon( alpha = 0.5) + geom_line(data = sm,aes(x = t,y = mean), linetype = 2) 
  scale_x_continuous("Time between samples (years)") + scale_y_continuous("SNPs between samples") + 
  annotate("Text",x=c(7,2,9), y=c(40,60,5), label=c("Transmission supported","Transmission unlikely","Consider sample error")) + 
  ggtitle("MRSA CC22 data")
ggsave("transmission_pred_interval_mrsa.pdf")

w<-which(sp$x < 1.1)
ggplot(sp[w,], aes(x = x, y = y)) + geom_polygon( alpha = 0.5) + geom_line(data = sm,aes(x = t,y = mean), linetype = 2) + 
  scale_x_continuous("Time between samples (years)", lim = c(0,1)) + scale_y_continuous("SNPs between samples", lim = c(0,20)) + 
  annotate("Text",x=c(0.5,0.25,0.85), y=c(8,16,0.5), label=c("Transmission supported","Transmission unlikely","Consider sample error")) + 
  ggtitle("MRSA CC22 data")
ggsave("transmission_pred_interval_mrsa_zoom.pdf")


####**** FOR N Gonorrhoeae data...? ****************************************************************************************************************************************************************#######

mu <-3.55 # (95% credibility interval 3.27-3.83) SNPs/genome/year.

# u can be estimated as the time to MRCA of isolates taken from within the
# same host at the same time (approximately). 
# s0 = SNP diversity in a single host = 2*mu*u
# [s0 is known in De Silva ref from two datasets, as is mu]. 
# So u can be estimated (u = s0 / (2mu)).

s0 = 2 #???

u = s0 / (2*mu) # time to MRCA

# s(t) can then be sampled from Pois(mu*(t+2u)) to generate 99% prediction intervals at each time t.

### Calculate exactly from distribution 
# z score at 99% = 2.576

vec_t <- seq(0,11,0.2)
s <- as.data.frame(matrix(0,length(vec_t),2))
sm <- as.data.frame(matrix(0,length(vec_t),1))
s$t <- vec_t
sm$t <- vec_t

for(i in 1:length(vec_t)){
  t = vec_t[i]
  lambda = mu * (t + 2*u)
  
  s[i,1:2] <- c(lambda + 2.576 * sqrt(lambda), lambda - 2.576 * sqrt(lambda))
  sm[i,1] <- lambda
  
  ## Could calculate by sampling
  #   rand_pois_sample <- rpois(100000,lambda)  
  #   q <- quantile(rand_pois_sample,probs = c(1,99)/100)
  #   qm <- mean(rand_pois_sample)
  #   s[i,1:2] <- q
  #   sm[i,1] <- qm
  
}

colnames(s)<- c("lo","hi","t")

# to generate polygon
sp1 <- as.data.frame(s[,c("lo","t")])
colnames(sp1)<- c("y","x")
sp2<-as.data.frame(cbind(rev(s$hi),rev(s$t)))
colnames(sp2)<- c("y","x")
sp <- rbind(sp1,sp2)
sp[which(sp$y < 0),"y"] <- 0 # 0 lowest SNP difference

sm <- as.data.frame(sm)
colnames(sm) <- c("mean","t")

ggplot(sp, aes(x = x, y = y)) + geom_polygon( alpha = 0.5) + geom_line(data = sm,aes(x = t,y = mean), linetype = 2) + 
  scale_x_continuous("Time between samples (years)") + scale_y_continuous("SNPs between samples") + 
  annotate("Text",x=c(7.7,2,9), y=c(27,46,3), label=c("Transmission supported","Transmission unlikely","Consider sample error")) + 
  ggtitle("N gonorrhoeae data")
ggsave("transmission_pred_interval_ngon.pdf")
