#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Post-processing on msAbund models
# Using spAbundance (Jeff's package)
# Wendy Leuenberger
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load packages ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(tidyverse)
library(magrittr)
library(spAbundance)
library(pracma)
library(spOccupancy)
library(coda)
library(abind)
library(WendyFunctions)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Define Functions

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# extractorRData is now in WendyFunctions package
# From https://stackoverflow.com/questions/65964064/programmatically-extract-an-object-from-collection-of-rdata-files
#' extractorRData <- function(file, object) {
#'   #' Function for extracting an object from a .RData file created by R's save() command
#'   #' Inputs: RData file, object name
#'   E <- new.env()
#'   load(file=file, envir=E)
#'   return(get(object, envir=E, inherits=F))
#' }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set TRUE/FALSE and pull file names ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# CHANGE HERE AS NEEDED
# I'm using a unique identifier to link all files from particular runs 
# or chains from particular runs with identical data
# Set those ID's here and we'll use them throughout the post-processing
# The goal is to reduce the number of places where these need to be set to 
# reduce human error and mixing up results with species/group names
# ID <- '-ResidentSpecialistHesPap24Spp-'; Run <- paste0('run', 1.1:5.1)
# ID <- '-ResidentSpecialistLycNymPieRio23Spp-'; Run <- paste0('run', 7.1:10.1)
# ID <- '-ResidentSpecialistLycNymPieRio21Spp-'; Run <- paste0('run', 1:4)
# ID <- '-ResidentGeneralistHesPap15Spp-'; Run <- paste0('run', 4.1:8.1)
# ID <- '-ResidentGeneralistLycNymPieRio19Spp-'; Run <- paste0('run', 4.1:7.1)
ID <- '-ResidentUnivoltineHesPap11Spp-'; Run <- paste0('run', 6.1:9.1)
# ID <- '-ResidentUnivoltineLycNymPieRio23Spp-'; Run <- paste0('run', c(2, 3, 5))
# ID <- '-ResidentUnivoltineLycNymPieRio23Spp-'; Run <- paste0('run', 11.1:15.1)
# ID <- '-Migratory20Spp-'; Run <- paste0('run', c(1:2))
# ID <- '-Migratory21Spp-'; Run <- paste0('run', 4.1:8.1)

# Run altogether
{
# Pull the file path as well so it can be used throughout
Path <- file.path(getwd(), 'Output/1_runs') %>% paste0('/')
# Path <- file.path('Y:/ButterflyChapterOne/Output/1_runs') %>% paste0('/')
# Set which run numbers we're interested in.
# Set next to ID to lessen confusion
# Hard to automate because might have multiple runs, so leave as manual
# Run <- c('run1', 'run2', 'run3', 'run4', 'run5')
IDRun <- paste0(ID, Run)
# Combine chains
CombineChains <- TRUE
# Print convergence report
ConvergenceCheck <- FALSE
# Save medians as future initial values
SaveMedians <- TRUE
# Subsample posterior to 1000 iterations
Subsample <- TRUE # OK to always leave this as true
# Print convergence report on the subsample
ConvergenceCheckSubsample <- FALSE
# Save out predictions
OutPred <- TRUE
# Calculate Abundance area under the curve
CalcAbundance <- TRUE
# Calculate Richness
CalcRichness <- TRUE

# Calculate Index
CalcIndexAny <-
  any(c(CalcAbundance,
        CalcRichness) ==
        TRUE) # DON'T CHANGE! Calculation, not a toggle

# Set seed
set.seed(seed = 761)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Get files and load out objects ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Grab the possible files that fit this unique ID. 
Files <- data.frame(Path = Path,
                    IDRun = IDRun) %>%
  mutate(PossibleFiles = paste0(Path, IDRun, '.RData'))
# Pull out the ID's from those files and put in a column. These should be the
# same as the ID we set above, but here's another way to get it and link it 
# at all steps. 
# Could maybe be used to pull out a broader set of files and keep the ID's 
# linked...idk though
Files %<>% 
  mutate(ID = str_extract(PossibleFiles, '\\-[:alnum:]+\\-'))

# Right now I only have 1 file, but I will have more in the future. 
# Here's a note thaMt I have to determine what to do about that. 
# But it's hard to do without results!
out.list <- vector('list', length = nrow(Files))
# ff <- 1 # For working on local computer
# Note - pull one run from HPCC to local computer to get data.list
for(ff in 1:nrow(Files)){
  print(paste('File', ff, 'of', nrow(Files), 'files'))
  out.list[[ff]] <- extractorRData(file = Files$PossibleFiles[ff],
                                   object = IDRun[ff])
  data.list <- extractorRData(file = Files$PossibleFiles[ff],
                              object = 'data.list')
}
# if(dim(Files)[1] == 1){
#   print('One file')
#   # Load if there's only a single file
#   out.1 <- extractorRData(file = Files$PossibleFiles,
#                           object = IDRun)
#   data.list <- extractorRData(file = Files$PossibleFiles,
#                               object = 'data.list')
# } else {
#   print(paste(nrow(Files), 'files'))
#   # Load multiple files. Maybe something like this?
  # out.1 <- extractorRData(file = Files$PossibleFiles[1],
  #                         object = IDRun[1])
  # out.2 <- extractorRData(file = Files$PossibleFiles[2],
  #                         object = IDRun[2])
  # out.3 <- extractorRData(file = Files$PossibleFiles[3],
  #                         object = IDRun[3])
#   # out.4 <- extractorRData(file = Files$PossibleFiles[4],
#   #                         object = IDRun[4])
#   # out.5 <- extractorRData(file = Files$PossibleFiles[5],
#   #                         object = IDRun[5])
#   data.list <- extractorRData(file = Files$PossibleFiles[1],
#                               object = 'data.list')
#   # out.3 <- extractorRData(file = Files$PossibleFiles[2],
#   # object = IDRun[3])
#   # Some might be old chains not interested in - how to filter out those? \n
#   # Condense back to one row? File instead of Files?')
# }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Combine chains ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(CombineChains == TRUE & nrow(Files) > 1) {
  print('Combining chains')
  out.full <- list()
  # Combine all MCMC samples together 
  out.full$n.chains <- length(out.list)
  # beta.comm.samples = community-level effects
  out.full$beta.comm.samples <- mcmc(
    do.call(
      rbind, lapply(out.list, function(one.chain) (one.chain$beta.comm.samples))
    )
  )
  # Community-level variances
  out.full$tau.sq.beta.samples <-mcmc(
    do.call(
      rbind, lapply(out.list, function(one.chain) (one.chain$tau.sq.beta.samples))
    )
  )
  # Species-level regression coefficients
  out.full$beta.samples <- mcmc(
    do.call(
      rbind, lapply(out.list, function(one.chain) (one.chain$beta.samples))
    )
  )
  # Species-specific NB overdispersion parameter (lower = more overdispersion)
  out.full$kappa.samples <- mcmc(
    do.call(
      rbind, lapply(out.list, function(one.chain) (one.chain$kappa.samples))
    )
  )
  # Overall random effect variances
  out.full$sigma.sq.mu.samples <- mcmc(
    do.call(
      rbind, lapply(out.list, function(one.chain) (one.chain$sigma.sq.mu.samples))
    )
  )
  # Species-specific random effect values
  out.full$beta.star.samples <- mcmc(
    do.call(
      rbind, lapply(out.list, function(one.chain) (one.chain$beta.star.samples))
    )
  )
  

    # Get RHat values for the main parameters 
    out.full$rhat <- out.list[[1]]$rhat
    # beta.comm
    tmplist <- mcmc.list(
      lapply(
        out.list, function(one.chain) (one.chain$beta.comm.samples)
        )
      )
    out.full$rhat$beta.comm <- 
      as.vector(gelman.diag(tmplist, autoburnin = FALSE)$psrf[, 2])
    
    # tau.sq.beta
    tmplist <- mcmc.list(
      lapply(
        out.list, function(one.chain) (one.chain$tau.sq.beta.samples)
      )
    )
    out.full$rhat$tau.sq.beta <- 
      as.vector(gelman.diag(tmplist, autoburnin = FALSE)$psrf[, 2])
    
    # beta
    tmplist <- mcmc.list(
      lapply(
        out.list, function(one.chain) (one.chain$beta.samples)
      )
    )
    out.full$rhat$beta <- 
      as.vector(gelman.diag(tmplist, autoburnin = FALSE)$psrf[, 2])
    
    # kappa
    tmplist <- mcmc.list(
      lapply(
        out.list, function(one.chain) (one.chain$kappa.samples)
      )
    )
    out.full$rhat$kappa <- 
      as.vector(gelman.diag(tmplist, autoburnin = FALSE)$psrf[, 2])
    
    # sigma.sq.mu
    tmplist <- mcmc.list(
      lapply(
        out.list, function(one.chain) (one.chain$sigma.sq.mu.samples)
      )
    )
    out.full$rhat$sigma.sq.mu <- 
      as.vector(gelman.diag(tmplist, autoburnin = FALSE)$psrf[, 2])

    # Run time
    # Setting the "overall" run time to be the longest run time across the three chains
    out.full$run.time <- do.call(
      rbind, 
      lapply(out.list, function(one.chain) (one.chain$run.time))
      ) %>%
      as_tibble %>%
      filter(elapsed == max(elapsed))
    
  # Get ESS values for the main parameters 
  out.full$ESS <- list()
  out.full$ESS$beta.comm <- effectiveSize(out.full$beta.comm.samples)
  out.full$ESS$tau.sq.beta <- effectiveSize(out.full$tau.sq.beta.samples)
  out.full$ESS$beta <- effectiveSize(out.full$beta.samples)
  out.full$ESS$kappa <- effectiveSize(out.full$kappa.samples)
  out.full$ESS$sigma.sq.mu <- effectiveSize(out.full$sigma.sq.mu.samples)
  
  # Other stuff 
  # This stuff is used under the hood for use with summaries, plotting, 
  # prediction, WAIC,  PPCs. 
  # Names of random effects, which are needed to keep track of stuff if 
  # doing any prediction
  out.full$re.level.names <- out.list[[1]]$re.level.names
  # Design matrix for fixed effects
  out.full$X <- out.list[[1]]$X
  # Design matrix for random effects
  out.full$X.re <- out.list[[1]]$X.re
  # The count data
  out.full$y <- out.list[[1]]$y
  # Information on the call to the function
  out.full$call <- out.list[[1]]$call
  # Overall number of samples run per chain
  out.full$n.samples <- out.list[[1]]$n.samples
  # Names of columns in C
  out.full$x.names <- out.list[[1]]$x.names
  # Species names
  out.full$sp.names <- out.list[[1]]$sp.names
  # Number of posterior samples saved for each chain
  out.full$n.post <- out.list[[1]]$n.post
  # Thinning rate
  out.full$n.thin <- out.list[[1]]$n.thin
  # Amount of burn-in
  out.full$n.burn <- out.list[[1]]$n.burn
  # Number of chains
  out.full$n.chains <- nrow(Files)
  # Distribution used for abundance (NB or Poisson)
  out.full$dist <- out.list[[1]]$dist
  # Names of columns with random effects
  out.full$re.cols <- out.list[[1]]$re.cols
  # Logical indicating if there were any random effects
  out.full$muRE <- out.list[[1]]$muRE
  
  # Make sure the new object has class msAbund 
  class(out.full) <- 'msAbund'
  
  save(out.full, data.list, file = paste0('Output/2_combined-chains/', ID, 
                               'combined-chains', 
                               '.RData'))
} else {
  if(nrow(Files) == 1){
    print("Only one chain, can't combine")
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nrow(Files) > 1) {
  print('Grabbing combined chains')
  if (exists('out.full')) {
    out <- out.full
  } else {
    out <- extractorRData(
      file = paste0(
        'Output/2_combined-chains/',
        ID,
        'combined-chains',
        '.RData'
      ),
      object = 'out.full'
    )
    data.list <- extractorRData(
      file = paste0(
        'Output/2_combined-chains/',
        ID,
        'combined-chains',
        '.RData'
      ),
      object = 'data.list'
    )
  }
} else {
  print('Only one chain')
  out <- out.list[[1]]
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Run convergence report ####

# Call a .Rmd process to check convergence and make traceplots

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(ConvergenceCheck == TRUE) {
  rmarkdown::render(
    'Output/Convergence/Convergence.Rmd',
    output_file = paste0(ID, 'convergence-check-pt1')
    )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Save medians ####

# Save median values of parameters from previous runs to 
# use as starting values for future runs

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(SaveMedians == TRUE) {
  # Get the medians for all parameter estimates. You would calculate these from your
  # first model run.
  beta.comm.meds <- apply(out$beta.comm.samples, 2, median)
  tau.sq.beta.meds <- apply(out$tau.sq.beta.samples, 2, median)
  beta.meds <- apply(out$beta.samples, 2, median)
  kappa.meds <- apply(out$kappa.samples, 2, median)
  beta.star.meds <- apply(out$beta.star.samples, 2, median)
  sigma.sq.mu.meds <- apply(out$sigma.sq.mu.samples, 2, median)
  
  save(beta.comm.meds,
      tau.sq.beta.meds,
      beta.meds,
      kappa.meds,
      beta.star.meds,
      sigma.sq.mu.meds,
    file = paste0('Output/3_medians/', ID,
                  'medians',
                  '.RData')
  )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Subsample posterior to 1000 ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Subsample == TRUE){
  # obtain 1000 evenly spaced samples from the out object
  # Go from 1 to the length of the samples for 
  indx <- seq(from = 1, to = dim(out$beta.comm.samples)[1], 
              length.out = 1000)
  # Make sure they're whole numbers
  # Doesn't matter if there's rounding error, just generally spaced
  indx %<>% round()
  
  # Subsample out object
  out$beta.comm.samples <- out$beta.comm.samples[indx, , drop = FALSE]
  out$tau.sq.beta.samples <- out$tau.sq.beta.samples[indx, , drop = FALSE]
  out$beta.samples <- out$beta.samples[indx, , drop = FALSE]
  out$kappa.samples <- out$kappa.samples[indx, , drop = FALSE]
  # out$y.rep.samples <- out$y.rep.samples[indx, , , , drop = FALSE]
  # out$mu.samples <- out$mu.samples[indx, , , , drop = FALSE]
  out$sigma.sq.mu.samples <- out$sigma.sq.mu.samples[indx, , drop = FALSE]
  out$beta.star.samples <- out$beta.star.samples[indx, , drop = FALSE]
  out$n.post <- length(indx) / out$n.chains
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Run convergence report ####

# Call a .Rmd process to check convergence and make traceplots

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(ConvergenceCheckSubsample == TRUE) {
  rmarkdown::render(
    'Output/Convergence/Convergence.Rmd',
    output_file = paste0(ID, 'convergence-check-subsampled')
  )
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculate indices ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(CalcIndexAny == TRUE | OutPred == TRUE){
  print('Preparing for predictions or index calculations')
  # Predict
  # Don't need species; the msAbund function already knows about species
  # Collect data to feed into the predict function
  # Predict at each site during each week (in period of modeling) during each year.
  week.0.linear <- data.list$covs$week.cov %>% scale %>% unique # All of my weeks unique values
  n.weeks <- length(week.0.linear) # Number of unique week values
  week.0.quadratic <- week.0.linear^2
  year.0 <- data.list$covs$year.cov %>% unique %>% sort
  n.years <- length(year.0)
  # Predict at counties and sites
  CountySite <- tibble(county.cov = data.list$covs$county.cov,
                       # Use NABA for predictions (all others set at 0)
                       # NABA is spread out throughout Midwest, so best for calculations
                       ia.ind = 0, 
                       il.ind = 0, 
                       mi.ind = 0,
                       oh.ind = 0,
                       openco.cov = data.list$covs$openco.cov,
                       opensi.cov = 0,
                       effort.cov = 1) %>%
    distinct() 
  # How many counties/sites and weeks/years
  J.0 <- nrow(CountySite)
  n.rep.0 <- length(week.0.linear) * length(year.0)
  # Create a J.0 x n.rep.0 matrix x p.abund matrix
  # Intercept
  int.0.full <- matrix(1, J.0, n.rep.0)
  # Predict relative abundance for a single week
  week.rep.0.full <- matrix(1, J.0, n.rep.0)
  # Replicates are ordered by year, then week within year.
  # Order of columns: week 1 year 1, week 2 year 1, ....
  week.0.full <- matrix(rep(rep(week.0.linear, each = J.0), times = n.years), J.0, n.rep.0)
  week.2.0.full <- week.0.full^2
  year.indx.full <- matrix(rep(year.0, each = n.weeks * J.0), J.0, n.rep.0)
  # Grab names of random and fixed effects
  REnames <- out$X.re %>% attr(which = 'dimnames') %>% extract2(3) %>% unique
  FEnames <- out$X %>% attr(which = 'dimnames') %>% extract2(3)
  Nparams <- length(REnames) + length(FEnames)
  # array with 3 dimensions for X (X is design matrix), 
  # J = site/county; 2nd dimension number weeks * number years; 3rd = number of covariates
  X.0 <- array(NA, dim = c(J.0, n.rep.0, Nparams))
  X.0[, , 1] <- int.0.full # All 1's
  X.0[, , 2] <- CountySite$county.cov #(use same order as in model)
  X.0[, , 3] <- CountySite$ia.ind
  X.0[, , 4] <- CountySite$il.ind
  X.0[, , 5] <- CountySite$mi.ind
  X.0[, , 6] <- CountySite$oh.ind
  X.0[, , 7] <- CountySite$openco.cov
  X.0[, , 8] <- CountySite$opensi.cov
  X.0[, , 9] <- week.0.full
  X.0[, , 10] <- week.2.0.full
  X.0[, , 11] <- year.indx.full
  X.0[, , 12] <- log(CountySite$effort.cov)
  # X.0[, , 13] <- CountySite$nsurveys.cov
  # Need scale() in there - exactly as in the formula
  dimnames(X.0)[[3]] <- c('(Intercept)', 'county.cov', 
                          'ia.ind', 'il.ind', 'mi.ind', 'oh.ind',
                          'openco.cov', 'opensi.cov',
                          'scale(week.cov)', 'I(scale(week.cov)^2)', 'year.cov',
                          'log(effort.cov)')
  
  # Predict in pieces to avoid a memory error 
  n.samples <- out$n.post * out$n.chains
  n.spp <- out$sp.names %>% n_distinct
  # Play around with the 10 - bigger might work
  vals <- split(1:J.0, ceiling(seq_along(1:J.0) / 10))
  # Number of pieces to predict over
  n.pieces <- length(vals)
  out.pred <- 0; abund.indx <- 0
  # Save data needed for combining files later - dimensions, etc.
  save(n.samples, n.spp, J.0, vals, n.years, n.pieces, n.weeks, data.list,
       file = 'Output/DataAndDimensions.RData')
  # For troubleshooting purposes:
  # ll <- 1
  # Outpred ####
  if(OutPred == TRUE){
    print('Making output predictions on model object')
    for (ll in 1:n.pieces) {
      rm(out.pred, abund.indx)
      gc()  # Returns memory to the system after large object removal
      J.curr <- length(vals[[ll]])
      # abund.indx <- array(0, dim = c(n.samples, n.spp, J.curr, n.years))
      print(paste("Currently on prediction piece ", ll, " out of ", n.pieces, sep = '')) 
      # log(mu[species, site, week, year]) = intercept + sitelevel covariates + 
      # Use covariates from site (give predict one row corresponding to original data)
      # msAbund knows that it's by species, so don't need to specify
      
      out.pred <- predict(out, X.0[vals[[ll]], , ])
      # Add names to mu.0.samples
      dimnames(out.pred[[1]])[[2]] <- out$sp.names
      # Add names to y.0.samples
      dimnames(out.pred[[2]])[[2]] <- out$sp.names
      save(out.pred, file = paste('Output/4_out-pred/', ID,  
                                  'outpred-piece-', 
                                  ll, '.RData', sep = '')) 
    } # ll
  } # OutPred
  # abund.indx ####
  if(CalcAbundance == TRUE){
    out.pred <- 0; abund.indx <- 0
    print('Calculating abundance indices')
    # Put abundance.indx into one array
    abund.indx <-  array(0, dim = c(n.samples, n.spp, J.0, n.years))
    # For troubleshooting purposes:
    # ll <- 1; jj <- 1; kk <- 1; ii <- 1; tt <- 1
    # n.pieces <- 1; n.samples <- 1; J.curr <- 1; n.years <- 1; n.spp <- 1
    for (ll in 1:n.pieces) {
      # Housekeeping - removing and creating objects
      rm(out.pred, abund.indx.tmp)
      gc()  # Returns memory to the system after large object removal
      J.curr <- length(vals[[ll]])
      abund.indx.tmp <- array(0, dim = c(n.samples, n.spp, J.curr, n.years))
      print(paste("Currently on abundance piece ", ll, " out of ", n.pieces, sep = '')) 
      
      # Load out.pred for this piece
      load(file = paste('Output/4_out-pred/', ID,  
                        'outpred-piece-', 
                        ll, '.RData', sep = '')) 
      
      # Calculate abundance index
      for (kk in 1:n.samples) {
        # Don't print out every sample because there's a bunch
        # Only print out it k is evenly divisible by 50
        if(kk == 1 | ( mod(kk, 50) == 0 ) ) {
          print(paste("Currently on sample ", kk, " of ", n.samples, sep = ''))
        }
        for (jj in 1:J.curr) {
          for (tt in 1:n.years) {
            curr.indx <- ((tt - 1) * n.weeks + 1):(tt * n.weeks)
            for (ii in 1:n.spp) {
              # Mean abundance - continuous value. Best for area under curve,
              # but use y.0 to add variability to richness/diversity calculations
              abund.indx.tmp[kk, ii, jj, tt] <- trapz(1:n.weeks, out.pred$mu.0.samples[kk, ii, jj, curr.indx])
              # Would use abund.indx for map of relative abundance (not output from post hoc linear reg)
            } # ii
          } # tt
        } # jj
      } # kk
      # Put into a summary indx
      abund.indx[, , vals[[ll]],] <- abund.indx.tmp
      # Save current results
      # save(abund.indx, file = paste('Output/5_abundance-indx/', ID,  
      #                               'abundance-piece-', 
      #                               ll, '.RData', sep = '')) 
    } # ll
    dimnames(abund.indx)[[2]] <- dimnames(out.pred$mu.0.samples)[[2]]
    save(abund.indx, file = paste('Output/5_abund-indx/', 
                                ID,
                                'abundance-indx.RData', 
                                sep = ''))
    # Pull out each species and save separately. Keeps species name and more 
    # manageable size
    for ( ii in 1:dim(abund.indx)[2] ) {
      Species.indx <- abund.indx[ , ii, , ]
      Species <- dimnames(out.pred$mu.0.samples)[[2]][ii]
      save(Species.indx, file = paste('Output/5_abund-indx/-', 
                                      Species,
                                      '-abundance-indx.RData', 
                                      sep = ''))
    } # ii species
      
  } # CalcAbundance
  
  # Richness, total abundance, spp. abundance y0 ####
  if(CalcRichness == TRUE){
    print('Calculating richness indices')
    # For troubleshooting purposes:
    # ll <- 1; jj <- 1; kk <- 1; ii <- 1; tt <- 1
    # n.pieces <- 1; n.samples <- 1; J.curr <- 1; n.years <- 1; n.spp <- 1
    for (ll in 1:n.pieces) {
      # Housekeeping - removing and creating objects
      rm(out.pred, richness.indx)
      gc()  # Returns memory to the system after large object removal
      J.curr <- length(vals[[ll]])
      y.0.AUC <- array(0, dim = c(n.samples, n.spp, J.curr, n.years))
      richness <- array(0, dim = c(n.samples, J.curr, n.years))
      TotalAbundance <- array(0, dim = c(n.samples, J.curr, n.years))
      print(paste("Currently on richness piece ", ll, " out of ", n.pieces, sep = '')) 
      
      # Load out.pred for this piece
      load(file = paste('Output/4_out-pred/', ID,  
                        'outpred-piece-', 
                        ll, '.RData', sep = '')) 
      
      # Calculate richness and total abundance
      for (kk in 1:n.samples) {
        # Don't print out every sample because there's a bunch
        # Only print out it k is evenly divisible by 50
        if(kk == 1 | ( mod(kk, 50) == 0 ) ) {
          print(paste("Currently on sample ", kk, " of ", n.samples, sep = ''))
        }
        for (jj in 1:J.curr) {
          for (tt in 1:n.years) {
            curr.indx <- ((tt - 1) * n.weeks + 1):(tt * n.weeks)
            # In out.pred, y.0.samples - NB predicted value, predicted integer (not mean, so discrete)
            # Richness: if y.0.samples > 0 for county, yr etc., then it's there, sum up species
            # y's - might have prediction of giant number due to negbin
            # More variation, if low probability, have a chance at getting observations >= 1 even if mean expectation is <1
            # 
            # y.0 has weekly estimates
            # Calculate how many individuals are present per species with AUC 
            # Adding up all abundances for all weeks together - Area under curve
            for (ii in 1:n.spp) {
              y.0.AUC[kk, ii, jj, tt] <- trapz(
                1:n.weeks, out.pred$y.0.samples[kk, ii, jj, curr.indx]
                )
            } # ii
            
            # Add up the number of species that were present at any point that year
            # Species that are present at any point in the year will have a value >0
            richness[kk, jj, tt] <-  sum(y.0.AUC[kk, , jj, tt] > 0)
            
            # Add all species together for total abundance overall (kk, jj, tt)
            TotalAbundance[kk, jj, tt] <- sum(y.0.AUC[kk, , jj, tt])
            
          } # tt
        } # jj
      } # kk
      dimnames(y.0.AUC)[[2]] <- dimnames(out.pred$y.0.samples)[[2]]
      
      # Save current results
      save(richness, file = paste('Output/6_richness/', ID,  
                                  'richness-piece-', 
                                  ll, '.RData', sep = '')) 
      # Save max abundance - combined abundance of all species in the group
      save(TotalAbundance, file = paste('Output/7_TotalAbund/', ID,
                                        'total-abund-piece', 
                                        ll, '.RData', sep = ''))
      # Save abundance per species to use for diversity calculation later
      save(y.0.AUC, file = paste('Output/8_y0AUC/', ID,
                                 'y0-AUC-piece', ll, '.RData', sep = ''))
    } # ll
  } # CalcRichness, MaxAbundance
} # CalcIndexAny and OutPred
