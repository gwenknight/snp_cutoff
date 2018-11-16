### CHECK OLD NEW

mu = 5/365
npat = 1000
nruns = 10000


# OLD
s_constant <- simu_runs_orig(180, mu, npat, nruns,0,0)
ggplot(s_constant$store,aes(x=value,fill = variable)) + geom_histogram(aes(y = ..density..),binwidth = 1,position = "identity") + 
  facet_wrap(~variable) + scale_x_continuous("SNP distance",limits = c(0,30)) + scale_y_continuous("Probability")

w<-which(s_constant$store$variable == "l99")
range(s_constant$store[w,"value"])
quantile(s_constant$store[w,"value"],probs=c(.025,.975))


# NEW
ss <- simu_runs(180,mu,npat,nruns)
ggplot(ss$store_limits, aes(x=value, fill = variable)) + geom_histogram(aes(y=..density..), binwidth = 1, position = "identity") +   
  facet_wrap(~variable)

w<-which(ss$store_limits$variable == "l99")
range(ss$store_limits[w,"value"])
quantile(ss$store_limits[w,"value"],probs=c(.025,.975))
