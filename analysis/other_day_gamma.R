### Other day samples - fit gamma

dd <- read.csv("~/Documents/snp_cutoff/data/mrsa_pros.csv")
u<-sort(unique(dd$liDayGap))
coef_store<-c(); data_store<-c()
library(fitdistrplus)

# trial
i = 14
w<-which(dd$liDayGap == u[i])
h <- hist(dd[w,]$liSNPdis,breaks=seq(0,max(dd[w,"liSNPdis"]),1),plot = FALSE)
data <- data.frame( x = h$breaks[-length(h$breaks)], y = h$counts)
plot(data$x,data$y/sum(data$y))
k = 1
s = 1.2
y <- 1/((s^k)*gamma(k))*data$x^(k-1)*exp(-(data$x/s))
lines(data$x,y, col="red")

fit.gamma <- fitdist(data$y/sum(data$y), distr = "gamma", 
                     method = "mle",start = list(shape = 1, rate = 1.2),lower = c(0,0))
fit.exp <- fitdist(data$y/sum(data$y), distr = "exp", method = "mle")
fit.exp <- fitdistr(data$y/sum(data$y), densfun = "exponential")
fit.poisson <- fitdistr(data$y/sum(data$y), distr = "poisson",densfun = "poisson", method = "mle")
gofstat(list(fit.exp,fit.poisson))

coef_store<-c(); data_store<-c()
for(i in 1:length(u)){
  w<-which(dd$liDayGap == u[i])
  if(length(w)>10){ # need > 10? for fit to be ok
    print(i)
    h <- hist(dd[w,]$liSNPdis,breaks=seq(0,max(dd[w,"liSNPdis"]),1),plot = FALSE)
    #data <- data.frame( x = h$breaks[-length(h$breaks)], y = h$counts/sum(h$counts))
    data <- data.frame( x = h$breaks[-length(h$breaks)], y = h$counts)
    #fit.gamma <- fitdist(data$y/sum(data$y), distr = "gamma", 
    #                     method = "mle",start = list(shape = 3, rate = 0.5),lower = c(0,0))
    fit.gamma <- fitdist(data$y, distr = "gamma", 
                         method = "mle",start = list(shape = 1, rate = 0.5),lower = c(0,0))
    
    plot(data$x,data$y)
    lines(data$x, dgamma(data$x, coef(fit.gamma)[1], coef(fit.gamma)[2]))
    data_store <- rbind(data_store,cbind(u[i],data,dgamma(data$x, coef(fit.gamma)[1], coef(fit.gamma)[2])))
    coef_store <- rbind(coef_store, c(u[i],coef(fit.gamma)))
    }
}

colnames(data_store) <- c("day","x","y","m")

coef_store <- as.data.frame(coef_store)
colnames(coef_store) <- c("day","shape","rate")

ggplot(data_store,aes(x=x, y=y)) + geom_point() + geom_line(aes(x=x,y=m),col="red") + 
  facet_wrap(~day,scales="free")

plot(coef_store$day,coef_store$shape) #,ylim = c(0,4))
# Assuming expoential link between day and coefficient 
mod_a <- nls(a ~ exp(aa + bb * day), data = coef_store, start = list(aa = 1.4, bb = -1),
           control = list(maxiter = 500))
lines(coef_store$day, predict(mod_a, list(x = coef_store$day)),col="red")
mod_a_mod <- coef(mod_a)

plot(coef_store$day,coef_store$rate)
# Assume all v similar - take mean
mean(coef_store[,3])
p_aa <- mean(coef_store[which(coef_store[,3] > (-8)),3]) # remove outliers


data_store$mod_general <- exp(exp(mod_a_mod[1] + mod_a_mod[2]*data_store[,"day"]) + p_aa * data_store[,"x"])
  
ggplot(data_store,aes(x=x, y=y)) + geom_point() + geom_line(aes(x=x,y=m),col="red") + 
  facet_wrap(~day,scales="free") + geom_line(aes(x=x,y=mod_general),col="blue")
ggsave("exponential_other_day.pdf",width = 12, height = 12)

## General model 
gen_mod_day <- function(day,l_x){
  x <- seq(0,l_x,1)
  aa_g <-  0.99338 #mod_a_mod[1]
  bb_g <- -0.02207 #mod_a_mod[2]
  out <- exp(exp(aa_g + bb_g*day) + p_aa * x)
  return(out)
}

mod_gen_plot <- c()
for(i in 1:35){
  days <- seq(0,350,10)
  mod_gen_plot <- rbind(mod_gen_plot,cbind(x,days[i],gen_mod_day(i,50)))
}
mod_gen_plot <- as.data.frame(mod_gen_plot)
colnames(mod_gen_plot) <- c("snp","day","count")

ggplot(mod_gen_plot,aes(x=snp,y=count,colour=factor(day))) + geom_line()
# don't change much! 

