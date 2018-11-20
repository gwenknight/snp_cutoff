###### Functions for simulating transmission ####

###*** Functions ***############################################################################################################

## Sample, with replacement, from your own distribution
rMydist <- function(n,x) {
  # n = number of samples
  # x = list of probabilities
  sample(seq(1,length(x),1), size = n, prob = x, replace=T)
}

## General model for how distance from transmitted changes over time
## parameters taken from data_fit.R 
## THIS IS FOR ST22 strain HO 5096 0412 mapping data   
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

## General model 
# THIS IS FOR ST30 strain MRSA252 mapping data         
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

#### Transmission simulation
simu_runs <- function(timesv, mu, npat, nruns, gen_mod = gen_mod_day){
  ### timesv = time when sample in days from transmission
  ### mu = mutation rate
  ### npat = number of patients
  ### nruns = number of runs
  ### same = if source same as recipient (1)
  ### randt = if random timings of sampling
  
  ### STORE / PREP
  store <- matrix(0,nruns,npat)
  store_limits <- matrix(0,nruns, 5)
  
  for(j in 1:nruns){
    ### SOURCE PATIENT
    # Distribution of SNPs timesv days after transmission
    dis_source <- gen_mod(timesv,150) 
    
    # Sample from this distribution to give number of SNPs different at this time point in SOURCE patients
    dis_source_from_transm <- rMydist(npat,dis_source)
    
    ### RECIPIENT PATIENT (of transmitted strain, timesv previously)
    # number of mutations likely to have occured since sample
    mut_number <- rpois(npat, mu * timesv)
    
    # how many mutations away if assume above number of mutations is the mean of an exponential distribution (1/mean = rate in this definition of exp)
    dis_receiver_from_transm <- c()
    for(i in 1:npat){
    dis_receiver_from_transm <- c(dis_receiver_from_transm,rexp(1,1/mut_number[i]))}

    
    # DISTANCE SOURCE <- TRANS -> RECIPIENT
    mut_dist <- dis_source_from_transm + dis_receiver_from_transm
    
    # STORE
    store[j,] <- mut_dist # store mutation distance
    
    # Calculate cut offs
    h<-hist(mut_dist,breaks = seq(0,300,1),plot=FALSE) 
    h_cum <- cumsum(h$density)
    l1<-length(which(h_cum < 0.7)) # 70% of transmission events captured
    l2<-length(which(h_cum < 0.9)) # 90% of transmission events captured
    l3<-length(which(h_cum < 0.95)) # 95% of transmission events captured
    l4<-length(which(h_cum < 0.99)) # 99% of transmission events captured
    # STORE
    store_limits[j,] <- c(j,l1,l2,l3,l4)
  }
  
  # OUTPUT
  store_limits <- as.data.frame(store_limits)
  colnames(store_limits) <- c("run","70%","90%","95%","99%")
  m_store_limits <- melt(store_limits,id.vars = "run")
  
  # RETURN
  # limits cutoffs / distribution at these limits / example distance
  return(list(store_limits = m_store_limits, store = store))
}