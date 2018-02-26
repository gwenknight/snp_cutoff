# Functions for simulating transmission

###*** Functions ***############################################################################################################
# Sample from own distribution
rMydist <- function(n,x) {
  sample(seq(1,length(x),1), size = n, prob = x, replace=T)
}
# Sample from distribution 
rMydist_plus <- function(n,x,n_mut) {
  sample(seq(1,length(x)+n_mut,1), size = n, prob = c(matrix(0,1,n_mut),x), replace=T) # zero probability of less than n_muts
}
## General model for how distance from transmitted changes over time
gen_mod_day <- function(day,l_x){
  x <- seq(0,l_x,1)
  aa_g <-  0.99338 #mod_a_mod[1]
  bb_g <- -0.02207 #mod_a_mod[2]
  out <- exp(exp(aa_g + bb_g*day) + p_aa * x)
  return(out)
}


#### RUNS
simu_runs <- function(timesv, mu, npat, nruns, same){
  # Store
  store_all <- c()
  store_limits <- c()
  
  for(j in 1:nruns){
    
    ### STORE / PREP
    n = npat # number of pairs of samples
    
    store <- matrix(0,n,3)
    
    # SAMPLE mutation rates to use
    if(length(mu)>1){
      mu_v <- rtruncnorm(n,a = 0, b = Inf, mu[1],mu[2]) # lower truncated to be positive
    } else{
      # CONSTANT mutation rate
      mu_v <- matrix(mu,1,n)}
    
    # SAMPLE times to use
    if(length(timesv > 1)){
      when_transm <- ceiling(runif(n,0,180)) # time to second sample
      ht <- hist(times,breaks = seq(0,500,1),plot = FALSE)
    } else {
      # CONSTANT times
      when_transm <- matrix(times_c,1,n)
      first_sample <- matrix(times_c,1,n)}
    
    for(i in 1:n){
      store[i,1:2] <- c(when_transm[i], first_sample[i]) # store time
      
      # SOURCE
      if(length(timesv > 1)){
        first_sample_i <- rMydist(1,ht$counts[which(ht$breaks < when_transm[i])]) # sample time depends on when_transm
      } else { first_sample_i <- first_sample[i]}
      
      if(same == 1){  
        mut_number <- rpois(1,mu_v[i] * first_sample_i) # number of mutations likely to have occured since sample
        dis_source_from_transm <- rexp(1,1/mut_number)
      } else {
        dis_source <- gen_mod_day(first_sample_i,150) 
        dis_source_from_transm <- rMydist(1,dis_source)}
      
      # RECIPIENT
      mut_number <- rpois(1,mu_v[i] * when_transm[i]) # number of mutations likely to have occured since sample
      dis_receiver_from_transm <- rexp(1,1/mut_number) # 1/mean (rate in this definition of exp)
      
      # DISTANCE SOURCE <- TRANS -> RECIPIENT
      mut_dist <- dis_source_from_transm + dis_receiver_from_transm
      store[i,3] <- mut_dist # store mut distance
    }
    
    #plot(store[,1], store[,2], xlab = "Time Distance (in days)", ylab = "SNP distance")
    h<-hist(store[,3],breaks = seq(0,300,1),plot=FALSE)
    h_cum <- cumsum(h$density)
    #plot(seq(1,39,1),h_cum)
    l1<-length(which(h_cum < 0.7)) # 70% 1 SNP differences
    l2<-length(which(h_cum < 0.9)) # 90% 4 SNPS different
    l3<-length(which(h_cum < 0.95)) # 95% 6 SNPS different 
    l4<-length(which(h_cum < 0.99)) # 99% 14 SNPS different
    
    ### Store
    store_all <- rbind(store_all,cbind(j,store))
    store_limits <- rbind(store_limits,c(j,l1,l2,l3,l4))
    
  }
  
  # output
  store_limits <- as.data.frame(store_limits)
  colnames(store_limits) <- c("run","l70","l90","l95","l99")
  m_store_limits <- melt(store_limits,id.vars = "run")
  
  # underlying data
  gg <- ggplot_build(p)
  gg_data <- gg$data[[1]]
  gg_data <- gg_data[,c("x","density","PANEL")]
  gg_data <- dcast(gg_data, x ~ PANEL, value.var = "density" )
  
  # RETURN
  return(list(store = m_store_limits, gg_data = gg_data))
}
