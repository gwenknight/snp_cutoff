# Mutation rate
# From Campbell et al 2018
# per site per year
m1 = 1.22 * 10^(-6) #(0.60, 1.86)
m2 = 1.30 * 10^(-6) #(1.20, 1.40)
m3 = 1.40 * 10^(-6) #(1.04, 1.80)
m4 = 3.30 * 10^(-6) #(2.50, 4.00)
m5 = 1.21 * 10^(-6)
m6 = 2.99 * 10^(-6)

mean_e = sum(m1,m2,m3,m4,m5,m6) / 6 
low_e =  sum(0.6,1.2,1.04,2.5,0,0)* 10^(-6) / 4 
high_e = sum(1.76,1.4,1.8,4,0,0)* 10^(-6) / 4 

d_mu <- as.data.frame(cbind(seq(1,6,1),
                            c(m1,m2,m3,m4,m5,m6),
                            c(0.6,1.2,1.04,2.5,0,0)* 10^(-6),
                            c(1.76,1.4,1.8,4,0,0)* 10^(-6)))
colnames(d_mu) <- c("exp","mean","low","high")

ggplot(d_mu, aes(x=exp,y=mean)) + geom_errorbar(aes(ymin = low,ymax = high)) + 
  scale_x_continuous("Estimate") + scale_y_continuous("Mutation rate \n(per site per year)") + 
  geom_hline(yintercept = mean_e,col="red") + geom_hline(yintercept = low_e,col="red",lty=2) + 
  geom_hline(yintercept = high_e,col="red",lty=2)


### Convert to per year and whole genome
mean_e = mean_e * 2842618 
low_e = low_e * 2842618 
high_e = high_e * 2842618 
mean_e
low_e
high_e
