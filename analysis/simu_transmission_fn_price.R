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
gen_mod_day_price <- function(day, param_general_fit, l_x = 150){
  ### day = number of days from first sample
  ### l_x = maximum number of SNPs
  ### param_general_fit = parameters from fit to data
  
  x <- seq(0,l_x,1)
  aa_g <- param_general_fit[1]
  bb_g <- param_general_fit[2]
  p_aa <- param_general_fit[4]
  out <- exp(exp(aa_g + bb_g*day) + p_aa * x)
  return(out)
}

#### Transmission simulation
# for price - only difference is that it runs a different gen_mod_day
simu_runs_price <- function(timesv, mu, npat, nruns, param_general_fit, t0dist, max_t0 = 0, gen_mod = gen_mod_day_price){
  ### timesv = time when sample in days from transmission
  ### mu = mutation rate
  ### npat = number of patients
  ### nruns = number of runs
  ### param_general_fit = parameters of generalised function 
  ### t0dist = distribution within source patient at time zero
  ### max_t0 = indicator as to whether to sample from t0dist or to take the maximum number of SNPs. 0 = baseline = sample. If change to 1 then use maximum number
  ### gen_mod = function type for generalised function
  
  ### STORE / PREP
  store <- matrix(0,nruns,npat)
  
  store_limits <- matrix(0,nruns,5)
  
  ### SOURCE PATIENT
  # Distribution of SNPs timesv days after transmission
  dis_source <- gen_mod(timesv,param_general_fit) # Fixed

  for(j in 1:nruns){
    # Sample from SOURCE distribution to give number of SNPs different at this time point in SOURCE patients
    dis_source_from_transm <- rMydist(npat,dis_source)
    
    ### RECIPIENT PATIENT (of transmitted strain, timesv previously)
    # number of mutations likely to have occured since sample
    mut_number <- rpois(npat, mu * timesv)
    r <- 1/mut_number
    
    # RECIPIENT: how many mutations away if assume above number of mutations is 
    # the mean of an exponential distribution (1/mean = rate in this definition of exp)
    # dis_receiver_from_transm <- c()
    # for(i in 1:npat){
    # dis_receiver_from_transm <- c(dis_receiver_from_transm,rexp(1,1/mut_number[i]))}
     
    if(max_t0 == 0){
    dis_receiver_from_transm <- rexp(n=1*length(r), rate=r) + rMydist(npat,t0dist) # sample from time zero distribution 
    }else{dis_receiver_from_transm <- rexp(n=1*length(r), rate=r) + rep(length(t0prob_dist), npat)} # take maximum possible SNP distance
    
    # DISTANCE SOURCE <- TRANS -> RECIPIENT
    mut_dist <- dis_source_from_transm + dis_receiver_from_transm
    
    # STORE
    store[j,] <- c(mut_dist) # store mutation distance
    
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

##########################################################################################################
#### SUB-SAMPLE function 

sub_sample_price <- function(data_all){
  
  host_ids = unique(as.vector(data_all$AnonymisedPatientId))
  length(host_ids)
  
  # Vector to store pairwise isolate comparisons to keep in each iteration
  host_isolates_kept = vector()
  
  # Across all hosts, sub-sample randomly to de-duplicate dataset
  for(h in 1:length(host_ids))
  {
    hhh = which(data_all$AnonymisedPatientId == host_ids[h])
    # Sample only one CC
    host_clonal_complexes = as.character(unique(c(as.vector(data_all$ClonalComplex1[hhh]), as.vector(data_all$ClonalComplex2[hhh]))))
    host_clonal_complex = sample(host_clonal_complexes,1)
    # Select comparison of sampled CC
    hhh = intersect(which(data_all$AnonymisedPatientId == host_ids[h]),which(data_all$ClonalComplex1 == host_clonal_complex))
    # Get all available collection dates
    host_collection_dates = unique(c(as.vector(data_all$CollectionDate1[hhh]),as.vector(data_all$CollectionDate2[hhh])))
    # If only one available collection data/sample > keep two random isolates
    if(length(host_collection_dates)==1)
    {
      host_isolates = vector()
      hhhd1 = which(data_all$AnonymisedPatientId == host_ids[h] & data_all$CollectionDate1 == host_collection_dates[1] & data_all$ClonalComplex1 == host_clonal_complex)
      if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(data_all$SequencingTag1[hhhd1])); }
      hhhd1 = which(data_all$AnonymisedPatientId == host_ids[h] & data_all$CollectionDate2 == host_collection_dates[1] & data_all$ClonalComplex2 == host_clonal_complex)
      if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(data_all$SequencingTag2[hhhd1])); }
      host_isolates = unique(host_isolates)
      host_isolates = sample(host_isolates, 2);
      host_isolates_kept = c(host_isolates_kept, host_isolates)
    } else
    {
      # Else, for each collection date, randomly sample one isolate
      for(d in 1:length(host_collection_dates))
      {
        host_isolates = vector()
        hhhd1 = which(data_all$AnonymisedPatientId == host_ids[h] & data_all$CollectionDate1 == host_collection_dates[d] & data_all$ClonalComplex1 == host_clonal_complex)
        if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(data_all$SequencingTag1[hhhd1])); }
        hhhd1 = which(data_all$AnonymisedPatientId == host_ids[h] & data_all$CollectionDate2 == host_collection_dates[d] & data_all$ClonalComplex2 == host_clonal_complex)
        if(length(hhhd1)>0){ host_isolates = c(host_isolates, as.vector(data_all$SequencingTag2[hhhd1])); }
        host_isolates = unique(host_isolates)
        host_isolates = sample(host_isolates, 1);
        host_isolates_kept = c(host_isolates_kept, host_isolates)
      }
    }
  }
  
  # print(paste("Number of isolates sub-sampled: ",length(host_isolates_kept), sep=""))
  
  # Keeping comparisons including isolates sub-sampled
  iii1 = which(!is.na(match(data_all$SequencingTag1, host_isolates_kept)))
  iii2 = which(!is.na(match(data_all$SequencingTag2, host_isolates_kept)))
  iii = iii1[which(!is.na(match(iii1,iii2)))]
  # print(paste("Number of pairwise comparisons sub-sampled: ",length(iii), sep=""))
  # print(paste("Number of patients sub-sampled: ",length(unique(data_all$AnonymisedPatientId[iii])), sep=""))
  # tmp = match(host_ids, data_all$AnonymisedPatientId[iii])
  # print(paste("Missing sub-sampled patients: ",paste(host_ids[which(is.na(tmp))], collapse = ";"), sep=""))
  
  ### Running linear mixed model
  data_all_sub = data_all[iii,]
  data_all_sub <- droplevels(data_all_sub)
  data_all_sub$AnonymisedPatientId=as.factor(data_all_sub$AnonymisedPatientId)
  
  return(data_all_sub)
  
}