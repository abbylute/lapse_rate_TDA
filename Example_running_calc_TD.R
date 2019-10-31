# Example of running calc_TD.R

# Setup
  mydir <- 'FinalCode/ForPublication/'
  
  library(ggplot2)
  library(parallel)
  source(paste0(mydir,'calc_TD.R'))

  # load example_data.txt, which contains station covariates. 
  # mydf should contain values representing one time period 
  # and one type of temperature variable (e.g. tmax or tmin)
  mydf <- read.table(paste0(mydir,'example_data.txt'), header=T)

# Set options
  # set the number of stations in each sample, ie the sample size
  # (see Figure 10 in Lute and Abatzoglou (2020) for guidance 
  # on choice of sample size)
  sampsz <- 2
  
  # set the maximum number of samples to evaluate
  # the default number of samples to evaluate for this sample size is
  ncol(combn(mydf$ID, sampsz))
  # you can set 'maxsamps' to limit the number of samples evaluated
  maxsamps <- 15000
  
  # set the number of cores to use for computing TD. Use detectCores() 
  # to determine how many cores are available.
  ncores <- 4

# Run calc_TD to calculate the Topoclimatic Dissimilarity metric
  td <- calc_TD(covdf=mydf, samplesize=sampsz, maxsamples=maxsamps, cores=ncores)

# calc_TD returns a data.frame with nrow=number of samples evaluated and two columns.
  # the first column, 'TD', is the topoclimatic dissimilarity metric value (lower values are more similar samples), 
  # and the second column, 'lapserate', is the lapse rate (°C/km) calculated for that sample using simple linear regression 
  

# Visualize how lapse rate uncertainty (spread) changes with TD  
  ggplot(td) +
    geom_boxplot(aes(x=TD, y=lapserate, group=cut_width(TD,2)))
  
  ggplot(td) +
    geom_boxplot(aes(x=TD, y=lapserate, group=cut_number(TD,8)))
  