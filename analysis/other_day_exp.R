### Other day samples - fit exponential

dd <- read.csv("~/Documents/snp_cutoff/data/mrsa_pros.csv")
u<-sort(unique(dd$liDayGap))
coef_store<-c(); data_store<-c()

# trial
i = 23
w<-which(dd$liDayGap == u[i])
h <- hist(dd[w,]$liSNPdis,breaks=seq(0,max(dd[w,"liSNPdis"]),1),plot = TRUE)
ggplot(dd[w,],aes(x=liSNPdis)) + geom_histogram(stat="count") + scale_x_continuous("SNP distance") + scale_y_continuous("Count")
ggsave("day_23_ggplot.pdf")
data <- data.frame( x = h$breaks[-length(h$breaks)], y = h$counts)
pdf("day_23.pdf")
plot(data$x,data$y)
lines(data$x,exp(1.665 + -0.976 * data$x), col="red")
dev.off()

for(i in 1:length(u)){
  w<-which(dd$liDayGap == u[i])
  if(length(w)>10){ # need > 10? for fit to be ok
    print(i)
    h <- hist(dd[w,]$liSNPdis,breaks=seq(0,max(dd[w,"liSNPdis"]),1),plot = FALSE)
    data <- data.frame( x = h$breaks[-length(h$breaks)], y = h$counts)
    mod <- nls(y ~ exp(a + b * x), data = data, start = list(a = log(h$counts[1]), b = -1),
               control = list(maxiter = 500))
    
    plot(data$x,data$y)
    lines(data$x, predict(mod, list(x = data$x)))
    data_store <- rbind(data_store,cbind(u[i],data,predict(mod, list(x = data$x))))
    coef_store <- rbind(coef_store, c(u[i],coef(mod)))
    }
}

colnames(data_store) <- c("day","x","y","m")

coef_store <- as.data.frame(coef_store)
colnames(coef_store) <- c("day","a","b")

ggplot(data_store,aes(x=x, y=y)) + geom_point() + geom_line(aes(x=x,y=m),col="red") + 
  facet_wrap(~day,scales="free")

plot(coef_store$day,coef_store$a,ylim = c(0,4))
# Assuming expoential link between day and coefficient 
mod_a <- nls(a ~ exp(aa + bb * day), data = coef_store, start = list(aa = 1.4, bb = -1),
           control = list(maxiter = 500))
lines(coef_store$day, predict(mod_a, list(x = coef_store$day)),col="red")
mod_a_mod <- coef(mod_a)

pdf("a_exponential.pdf")
plot(coef_store$day,coef_store$a,ylim = c(0,4),xlab = "Time between samples",ylab="Intercept") 
lines(coef_store$day, predict(mod_a, list(x = coef_store$day)),col="red")
dev.off()

pdf("b_flat.pdf")
plot(coef_store[,1],coef_store[,3],xlab = "Time between samples",ylab = "Slope")
# Assume all v similar - take mean
mean(coef_store[,3])
p_aa <- mean(coef_store[which(coef_store[,3] > (-8)),3]) # remove outliers
abline(h = p_aa,col="red",lty="dashed")
dev.off()


data_store$mod_general <- exp(exp(mod_a_mod[1] + mod_a_mod[2]*data_store[,"day"]) + p_aa * data_store[,"x"])
  
ggplot(data_store,aes(x=x, y=y)) + geom_point() + geom_line(aes(x=x,y=m),col="red") + 
  facet_wrap(~day,scales="free") + geom_line(aes(x=x,y=mod_general),col="blue") + 
   scale_x_continuous("SNP distance") + scale_y_continuous("Count")
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

