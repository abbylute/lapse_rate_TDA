
# Topoclimatic Dissimilarity Approach (TDA)
#--------------------------------------------------
# A.C. Lute (Oct, 2019)
# for more details regarding the methods, application, and evaluation of TDA, see 
# Lute and Abatzoglou (2020) .......

# requires: tidyr, dplyr, parallel, data.table


# INPUTS:
# covdf: data.frame of covariate values for each station, organized like this:
#        each row represents one station
#        first column is some station identifier, named 'ID'
#        second column is the station elevations, named 'elev'
#        third column is temperature values, named 'temp'
#        subsequent columns are covariates

# samplesize: number of stations in each sample

# maxsamples: maximum number of samples to evaluate for any sample size.
#        the default is the total number of combinations of stations of size 'samplesize' drawn from covdf

# cores: number of cores to run at once for parallelization, must be an integer, default is 1

# OUTPUTS:
# for each sample evaluated,
# TD: the Topoclimatic Dissimilarity metric (smaller values indicate greater sample similarity)
# lapserate: the lapse rate (°C/km) computed for each sample using simple linear regression




calc_TD = function(covdf, samplesize, maxsamples=ncol(combn(covdf$ID,samplesize)), cores=1){
  
  # Load required packages:
  if(!all(c('dplyr','tidyr','parallel','data.table') %in% rownames(installed.packages())))
    stop("Error: calc_TD requires the following packages:

         dplyr
         tidyr
         parallel
         data.table

         One or more of these packages are currently not installed or not available.
         Run install.packages(packagename) to install them.")
  
  library(dplyr)
  library(tidyr)
  library(parallel)
  library(data.table)
  

  
  covcols <- which(!names(covdf) %in% c('ID','temp'))
  covs <- names(covdf) [covcols] # names of covariates
  
  # Check the format of covdf
  if(!all(c('ID','elev','temp') %in% names(covdf)))
    stop("Error: covdf must contain columns 'ID','elev', and 'temp', for station identifiers, elevations, and temperature values, respectively")

  if(!all(apply(covdf[covcols], 2, is.numeric)))
    stop("Error: covariates must only have values of class numeric")
  
  
  # Check samplesize
  nstations <- nrow(covdf) # total number of stations
  if(!is.numeric(samplesize))
    stop("Error: 'samplesize' must be a numeric whole number")
  if(samplesize > nstations)
    stop("Error: 'samplesize' must be equal to or less than the number of stations (rows) in covdf")
  
  
  # Create weights
  # weights are the sqrt of absolute value of the correlation between temperature and each covariate
  wts <- tibble::enframe(apply(covdf[,covcols], 2, function(x) sqrt(abs(cor.test(covdf$temp,x)$estimate))), name='covnm', value='covcor')
  
  
  # Scale the covariates
  covdf_sc <- covdf %>% dplyr::select(covs) %>%  mutate_each(scale)
  covdf_sc <- cbind(covdf[,-covcols], covdf_sc)
  
  
  # Compute all possible combinations of stations of size 'samplesize'
  samples <- t(combn(as.character(covdf_sc$ID), samplesize)) 
  
  
  # Limit the number of samples if necessary:
  if (nrow(samples) > maxsamples){
    samples <- samples[sample(1:nrow(samples), maxsamples, replace=F),]
  }
  
  
  # Calculate the TD metric for each sample:
  cl <- makeCluster(cores)
  clusterExport(cl, c("TDA","covdf","covdf_sc","wts","covs"), envir=environment())
  
  print(paste0('running TDA for sample size ', samplesize))
  tda_df <- parApply(cl, samples, 1, function(i) TDA(samp=i, covdf.=covdf, covdf_sc.=covdf_sc, wts.=wts, covs.=covs))
  stopCluster(cl)
  
  tda_df <- rbindlist(tda_df)
  
} # end runTDA


# The TDA function is used in line 97 of calc_TD().
# Since it is called within parApply(), TDA takes one sample at a time ('samp') from 'samples' and computes
# the TD metric
TDA = function(samp, covdf., covdf_sc., wts., covs.){
  library(dplyr)
  library(tidyr)
  
  nn <- length(samp)
  tab <- covdf_sc %>% filter(ID %in% samp)
  tabmeta <- tab %>% dplyr::select(ID, elev)
  
  # Get covariate ranges:
  cov_rngs <- tab %>% 
    pivot_longer(covs, names_to='covnm', values_to='value') %>% 
    left_join(tabmeta, by=c('ID')) %>% 
    group_by(covnm, add=F) %>% 
    dplyr::summarise(cov_range=diff(range(value))) %>% 
    left_join(wts, by='covnm')
  
  # Rescale elevation ranges so that larger ranges get smaller numbers 
  evr <- as.vector(dist(covdf_sc$elev, method = "manhattan")) # all possible elevation ranges
  cov_rngs$cov_range[which(cov_rngs$covnm=='elev')] <- max(evr)-cov_rngs$cov_range[which(cov_rngs$covnm=='elev')]
  rm(evr); gc()
  
  # Compute similarity 
  samp_td <- cov_rngs %>% 
    mutate(cov_range_wt = cov_range*abs(covcor)) %>% 
    summarise('td'=sum(cov_range_wt))
  
  # Compute the lapse rate
  samp_lr <- covdf %>% 
    filter(ID %in% samp) %>%
    summarise(lapse=summary(lm(temp~elev))$coefficients[2]*1000) 
  
  out <- data.frame('TD'=samp_td[[1]], 'lapserate'=samp_lr[[1]])
  

}






