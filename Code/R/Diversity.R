#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Combined metrics - all species all groups together 
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

# library(doParallel) # For running each group on its own node

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set TRUE/FALSE and pull file names ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# CHANGE HERE AS NEEDED
# I'm using a unique identifier to link all files from particular runs 
# or chains from particular runs with identical data
# Set those ID's here and we'll use them throughout the post-processing
# The goal is to reduce the number of places where these need to be set to 
# reduce human error and mixing up results with species/group names
IDs <- c(
  '-ResidentSpecialistHesPap24Spp-',
  '-ResidentSpecialistLycNymPieRio23Spp-',
  '-ResidentGeneralistHesPap15Spp-',
  '-ResidentGeneralistLycNymPieRio19Spp-',
  '-ResidentUnivoltineHesPap11Spp-',
  '-ResidentUnivoltineLycNymPieRio23Spp-',
  '-Migratory21Spp-'
)
IDShort <- data.frame(ID = IDs,
                      Short = c(
                        'RSHP',
                        'RSLNPR',
                        'RGHP',
                        'RGLNPR',
                        'RUHP',
                        'RULNPR',
                        'Migrant'
                      ))

# File path
# RemoteHPCC <- FALSE
RemoteHPCC <- TRUE
if(RemoteHPCC == TRUE){
  path = 'Y:/ButterflyChapterOne/'
} else {
  path = ''
}


# Load object with dimensions 
load(file = 'Output/DataAndDimensions.RData')
# Load object with county information
load(file = paste0(path, 'Data/CountyKey.RData'))
# Source code for species traits
source('Code/R/SpeciesFunctionalTraits.R')
# Load SpeciesLists for info on which species are in which group
# Only changes if abundances change - if which species are rare/common would be 
# affected
load(file = paste0(path, 'Output/12_AbundanceTrends/SpeciesLists.RData'))

n.samples <- 1000
Pieces <- 1:n.pieces
week.0.linear <- data.list$covs$week.cov %>% scale %>% unique # All of my weeks unique values
n.weeks <- length(week.0.linear) # Number of unique week values



# Toggles if needed ####
CalcCombinedAbundance <- FALSE
CalcSeparateProportions <- FALSE
PostHocHprime <- FALSE
WithFunctions <- FALSE
CreateTraitDiversity <- TRUE
GroupTraitAbundance <- FALSE
GroupTraitHPrime <- FALSE
AnalyzeTraitDiversity <- TRUE
# Set seed
set.seed(seed = 761)

# Registering a parallel backend, using the 'multicore' functionality
# registerDoParallel(cores = as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE')[1]))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculate total abundance for diversity ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(CalcCombinedAbundance == TRUE){
  print('Combining total abundance among groups')
  # Pull the file path for the total abundance folder
  Path <- paste0(path, 'Output/7_TotalAbund/') 
  # All of the possible files in one dataframe
  Files <- expand.grid(ID = IDs, 
                       Piece = Pieces) %>% 
    mutate(PossibleFiles = paste0(Path, ID, 'total-abund-piece', Piece, '.RData')) %>% 
    left_join(IDShort)
  # Loop through each piece.
  # TotalAbundance objects are small because all species are combined (1 less dimension)
  # But we'll need to run the next calculations on y.0.AUC's, which are larger
  # So let's keep it split into pieces. Smaller pieces might queue faster on HPCC
  for (ll in 1:n.pieces) {
    print(paste("Currently on combining abundance piece ", ll, " out of ", n.pieces, sep = '')) 
    # Housekeeping - removing and creating objects
    rm(AbundanceCombined, AbundanceSeparated)
    gc()  # Returns memory to the system after large object removal
    # Just grab the files for a given piece
    FilePieces <- Files %>% 
      filter(Piece == ll)
    # Which set of counties are we looking at (numbers/index)
    J.curr <- length(vals[[ll]])
    # Create the empty array to store everything in
    AbundanceCombined <- array(0, dim = c(n.samples, J.curr, n.years))
    # Pull out each of the functional groups
    n.groups <- n_distinct(FilePieces$Short)
    groups <- unique(FilePieces$Short)
    AbundanceSeparated <- array(0, dim = c(n.samples, J.curr, n.years, n.groups))
    for (gg in 1:n.groups) {
      AbundanceSeparated[, , , gg] <-
        extractorRData(file = FilePieces$PossibleFiles[FilePieces$Short == groups[gg]],
                       object = 'TotalAbundance')
    }
    
    # Add them all together for each iteration, county, and year
    # kk <- 1; jj <- 1; tt <- 1 # For troubleshooting
    for (kk in 1:n.samples) {
      # Don't print out every sample because there's a bunch
      # Only print out it k is evenly divisible by 50
      if(kk == 1 | ( mod(kk, 50) == 0 ) ) {
        print(paste("Currently on sample ", kk, " of ", n.samples, sep = ''))
      }
      for (jj in 1:J.curr) {
        for (tt in 1:n.years) {
          curr.indx <- ((tt - 1) * n.weeks + 1):(tt * n.weeks)
          # Add up across the groups 
          # Leave the last index blank in AbundanceSeparated
          AbundanceCombined[kk, jj, tt] <- sum(AbundanceSeparated[kk, jj, tt, ])
        } # tt
      } # jj
    } # kk
    # Save current results
    # Save combined abundance to use for diversity calculation later
    save(AbundanceCombined, 
         file = paste('Output/9_CombinedAbundance/CombinedAbundance-piece',
                      ll, '.RData', sep = ''))
    
  } # ll
} # CalcCombinedAbundance


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculate Proportions of each group ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(CalcSeparateProportions == TRUE) {
  print('Calculating proportions of spp within groups using combined abundance')
  # Grab the species-specific predicted abundance file paths
  Path <- paste0(path, 'Output/8_y0AUC/')
  # Path for each file
  Files <- expand.grid(ID = IDs,
                       Piece = Pieces) %>%
    mutate(PossibleFiles = paste0(Path, ID, 'y0-AUC-piece', Piece, '.RData')) %>%
    left_join(IDShort)
  HprimeAll <- array(0, dim = c(n.samples, J.0, n.years))
  # In pieces for memory
  for (ll in 1:n.pieces) {
    print(paste("Currently on diversity piece ", ll, " out of ", n.pieces, sep = '')) 
    rm(SppAbundance, Proportions, pi_x_ln_pi, sum_pixlnpi, Hprime)
    gc()
    # The overall abundance of all groups combined (denominator)
    FilesCombinedAbundance <- file.path(
      paste0(path,
        'Output/9_CombinedAbundance/CombinedAbundance-piece',
        ll,
        '.RData'
      )
    )
    AbundanceCombined <-
      extractorRData(file = FilesCombinedAbundance,
                     'AbundanceCombined')
    # Just the files from this piece
    FilePieces <- Files %>%
      filter(Piece == ll)
    # Which set of counties are we looking at (numbers/index)
    J.curr <- length(vals[[ll]])
    # Number of groups and their names
    n.groups <- n_distinct(FilePieces$Short)
    groups <- unique(FilePieces$Short)
    SppAbundance <- vector(mode = 'list', length = n.groups)
    Proportions <- vector(mode = 'list', length = n.groups)
    pi_x_ln_pi <- vector(mode = 'list', length = n.groups)
    sum_pixlnpi <-
      array(0, dim = c(n.samples, J.curr, n.years, n.groups))
    Hprime <- array(0, dim = c(n.samples, J.curr, n.years))
    # SppAbundance <- array(0, dim = c(n.samples, J.curr, n.years, n.groups))
    for (group.gg in 1:n.groups) {
      SppAbundance[[group.gg]] <-
        extractorRData(file = FilePieces$PossibleFiles[FilePieces$Short == groups[group.gg]],
                       object = 'y.0.AUC')
      Proportions[[group.gg]] <-
        array(0, dim = c(n.samples, dim(SppAbundance[[group.gg]])[2], J.curr, n.years))
      pi_x_ln_pi[[group.gg]] <-
        array(0, dim = c(n.samples, dim(SppAbundance[[group.gg]])[2], J.curr, n.years))
    } # group.gg
    # Add them all together for each iteration, county, and year
    # kk <- 1; ii <- 1; jj <- 1; tt <- 1 # For troubleshooting
    for (kk in 1:n.samples) {
      # Don't print out every sample because there's a bunch
      # Only print out it k is evenly divisible by 50
      if (kk == 1 | (mod(kk, 50) == 0)) {
        print(paste("Currently on sample ", kk, " of ", n.samples, sep = ''))
      }
      for (jj in 1:J.curr) {
        for (tt in 1:n.years) {
          curr.indx <- ((tt - 1) * n.weeks + 1):(tt * n.weeks)
          # Add up across the groups
          # Leave the last index blank in AbundanceSeparated
          for (gg in 1:n.groups) {
            for (ii in 1:dim(SppAbundance[[group.gg]])[2]) {
              Proportions[[group.gg]][kk, ii, jj, tt] <-
                SppAbundance[[group.gg]][kk, ii, jj, tt] /
                AbundanceCombined[kk, jj, tt]
              # Proportion * ln proportion
              # Sum these all up to get H'
              pi_x_ln_pi[[group.gg]][kk, ii, jj, tt] <-
                Proportions[[group.gg]][kk, ii, jj, tt] *
                log(Proportions[[group.gg]][kk, ii, jj, tt])
            } # ii species
            sum_pixlnpi[kk, jj, tt, group.gg] <-
              sum(pi_x_ln_pi[[group.gg]][kk, , jj, tt],
                  na.rm = TRUE)
          } # gg group
          Hprime[kk, jj, tt] <-
            -1 * sum(sum_pixlnpi[kk, jj, tt, ])
        } # tt
      } # jj
    } # kk
    # Combine diversity values into one array instead of in pieces
    HprimeAll[ , vals[[ll]], ] <- Hprime
    # } # gg group
  } # ll
  # Save combined richness for analyses later
  save(HprimeAll, 
       file = 'Output/11_CombinedHprime/CombinedHprime.RData')
  
} # CalcProportions

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# PostHoc Regression on Diversity data ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(PostHocHprime == TRUE) {
  print('Running Hprime post hoc LMs')
  # Load file
  # For all species in one index. Deprecated because file is too big
  if(exists('HprimeAll') == FALSE) {
    HprimeAll <- extractorRData(file = paste0(path, 'Output/11_CombinedHprime/CombinedHprime.RData'),
                                  object = 'HprimeAll')
  }
  
  DateCreated <- file.info(paste0(path, 'Output/11_CombinedHprime/CombinedHprime.RData'))
  # When I fixed the year order problem (Thanks Jeff!!!)
  # file.info('Code/R/PostProcess_MsAbund.R')$mtime
  if(DateCreated$mtime < "2024-05-30 15:00:54 EDT") {
    year.0.wrong <- c(3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
                      21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 1, 2, 4, 6, 32)
    HprimeAll <- HprimeAll[ , , order(year.0.wrong)]
  }
  # Create empty data frame to store post hoc results
  post.hoc.df.Hprime <- data.frame(
    trend.est = NA,
    int.est = NA,
    trend.prob.pos = NA,
    trend.low = NA,
    trend.high = NA
  )
  J <- ncol(HprimeAll)
  n.samples <- nrow(HprimeAll)
  unique.years <- 1:dim(HprimeAll)[3]
  plot.years <- unique.years
  years <- rep(1:length(unique.years), each = J) 
  n.years <- length(unique(years))
  sites <- rep(1:J, times = n.years)
  # Format data for postHocLM
  y <- matrix(HprimeAll, nrow = n.samples, ncol = J * n.years)
  # Use log transformation because the data are small and near zero
  # Intercept and random effect variance jump around a lot 
  new.data.list <- list(y = log(y), 
                        covs = data.frame(years, sites)) 
  out.posthoc <- postHocLM(formula = ~ years + (1 | sites), 
                           data = new.data.list, 
                           n.chains = 1, # keep at one
                           verbose = TRUE)
  out.posthoc %>% summary
  # Diagnostics
  # new.data.list$y %>% hist
  # plot(out.posthoc$beta.samples, density = FALSE)
  # plot(out.posthoc$beta.star.samples, density = FALSE)
  # plot(out.posthoc$beta.samples[,1] + out.posthoc$beta.star.samples, density = FALSE)
  # str(out.posthoc)
  # out.posthoc$y.hat.samples # predicted values
  # plot(c(out.posthoc$y.hat.samples - new.data.list$y), c(out.posthoc$y.hat.samples))
  # MedianResid <- apply(out.posthoc$y.hat.samples - new.data.list$y, 2, median)
  # HprimeAll[1, 1:10, 1]
  # y[1, 1:10]
  
  # Save results and generate summary figure 
  prob.pos <- apply(out.posthoc$beta.samples, 2, function(a) mean(a > 0))
  beta.star.mean <- out.posthoc$beta.star.samples %>% colMeans()
  site.ranef <- mean(beta.star.mean)
  post.hoc.df.Hprime$trend.prob.pos <- prob.pos[2]
  # y.log.med.full <- matrix(apply(log(y), 2, median), J, n.years)
  # abund.indx.log.site.avg <- apply(log(abund.indx), c(1, 3), mean)
  # abund.indx.log.quants <- apply(abund.indx.log.site.avg, 2, quantile, probs = c(0.025, 0.5, 0.975))
  # If not logged, use these two instead
  # Hprime.log.quants <- apply(HprimeAll[ , , ], 3, quantile, probs = c(0.025, 0.5, 0.975))
  # beta.quants <- apply(out.posthoc$beta.samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
  # Logged values and back transformations
  Hprime.log.quants <- apply(log(HprimeAll[ , , ]), 3, quantile, probs = c(0.025, 0.5, 0.975))
  beta.quants <- apply(exp(out.posthoc$beta.samples), 2, quantile, probs = c(0.025, 0.5, 0.975))
  beta.log.quants <- apply(out.posthoc$beta.samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
  post.hoc.df.Hprime$trend.est <- beta.log.quants[2, 2]
  post.hoc.df.Hprime$int.est <- beta.log.quants[2, 1]
  post.hoc.df.Hprime$trend.low <- beta.log.quants[1, 2]
  post.hoc.df.Hprime$trend.high <- beta.log.quants[3, 2]
  unique.years <- unique(out.posthoc$X[, 'years'])
  plot.min <- min(Hprime.log.quants)
  plot.max <- max(Hprime.log.quants)
  plot.vals <- seq(from = min(unique.years), to = max(unique.years), length.out = 100)
  
  # Generate min and max y values
  for (ss in 1:1000) {
    tmp <- out.posthoc$beta.samples[ss, 1] + out.posthoc$beta.samples[ss, 2] * plot.vals
    if (min(tmp) < plot.min) {
      plot.min <- min(tmp)
    }
    if (max(tmp) > plot.max) {
      plot.max <- max(tmp)
    }
  }
  Trend <- post.hoc.df.Hprime$int.est + 
    post.hoc.df.Hprime$trend.est * plot.vals 
  
  # TODO: Fix and make line the mean
  Lines <- data.frame(year = plot.vals,
                      med = Trend)
  # sp.name.plot <- simpleCap(str_replace(curr.sp.name, '-', ' '))
  plot.df <- data.frame(med = Hprime.log.quants[2, ], 
                        low = Hprime.log.quants[1, ],
                        high = Hprime.log.quants[3, ], 
                        year = plot.years)
  # ggplot structure
  my.plot <- ggplot(data = plot.df, aes(x = year, y = med))
  # 1 line per sample (takes awhile)
  for (ss in sample(1:1000, 200)) {
    my.plot <- my.plot + geom_abline(slope = out.posthoc$beta.samples[ss, 2],
                                     intercept = out.posthoc$beta.samples[ss, 1],
                                     color = 'darkorchid4', alpha = 0.1)
  }
  # Add points, error bars, axis labels
  my.plot <- my.plot +
    geom_point(size = 3) +
    geom_segment(aes(x = year, y = low, xend = year, yend = high)) +
    theme_classic(base_size = 18) +
    # labs(x = 'Year', y = 'Hprime')
    labs(x = 'Year', y = 'Shannon Diversity (log scale)')
  # Add stats and species name
  # my.plot <- my.plot +
  #   ggtitle(label = paste("Hprime: P(trend < 0) = ", round(1 - prob.pos[2], 2), 
  #                         "\n% Change/Year: ", 
  #                         round((beta.quants[2, 2] - 1) * 100, 2), " (", 
  #                         round((beta.quants[1, 2] - 1) * 100, 2), ", ", 
  #                         round((beta.quants[3, 2] - 1) * 100, 2), ")", sep = ''))
  my.plot <- my.plot +
    ggtitle(label = paste("Diversity: P(trend < 0) = ", round(1 - prob.pos[2], 2), 
                          "\n% Change/Year: ", 
                          round((beta.quants[2, 2]), 3), " (", 
                          round((beta.quants[1, 2]), 3), ", ", 
                          round((beta.quants[3, 2]), 3), ")", sep = ''))
  # Add trendline
  my.plot <- my.plot + 
    geom_line(data = Lines, linewidth = 1)
  # lines(plot.vals, tmp, col = 'black', lty = 1, lwd = 5)
  
  pdf(paste('Output/trend-figures/diversity-trend.pdf', sep = ''),
      width = 7, height = 7)
  print(my.plot)
  dev.off()
  
  save(out.posthoc, new.data.list, 
       file = paste('Output/trend-figures/Hprime-posthoc-for-figure.rda',
                    sep = ''))
  
}  


# With functions ####
if (WithFunctions == TRUE) {
  DiversityLM <- ButterflyPostHoc(Data = HprimeAll[, , ], Logged = TRUE)
  DiversityTrendPlot <- ButterflyTrendPlot(
    Data = HprimeAll[, , ],
    Model = DiversityLM,
    Logged = FALSE,
    Name = 'Hprime',
    ylabel = 'Shannon Diversity (log scale)',
    SaveFig = FALSE
  )
}


# Group diversity ####

# Calculate group diversity ####
if ( CreateTraitDiversity == TRUE ) {
  # Pull out the relevant files
  Path <- paste0(path, 'Output/8_y0AUC/')
  SpeciesFileList <- list.files(path = Path, pattern = "[[:upper:]]{5,}-*")
  ColSp <- list.files(path = Path, pattern = 'COL-SP')
  SppFiles <- tibble(Path, SpeciesFileList) %>%
    mutate(SpeciesNames = str_extract(string = SpeciesFileList,
                                      pattern = '[[:alnum:]]{6,}')) %>%
    add_row(Path = Path,
            SpeciesFileList = ColSp,
            SpeciesNames = 'COL-SP') %>%
    mutate(FullPath = paste0(Path, '/', SpeciesFileList))
  
  # Calculate total abundance for each group
  if ( GroupTraitAbundance == TRUE ) {
    
    # Make an empty list to hold the group abundance arrays
    AbundanceList <- vector('list', length = length(SpeciesLists))
    names(AbundanceList) <- names(SpeciesLists)
    
    for (group.hh in 1:length(AbundanceList)) {
      print(paste(
        'Adding species for group',
        group.hh,
        'of',
        length(AbundanceList)
      ))
      # Group files
      ThisGroup <- SpeciesLists[[group.hh]]
      GroupFiles <- SppFiles %>%
        filter(SpeciesNames %in% ThisGroup)
      
      GroupAbundance <- array(0, dim = c(n.samples, J.0, n.years))
      # Add each of the species together
      for (species.ss in 1:nrow(GroupFiles)) {
        ThisSpecies <- extractorRData(file = GroupFiles$FullPath[species.ss], object = 'SppArray')
        
        # Add to the Group so far
        GroupAbundance <- GroupAbundance + ThisSpecies
        
        # Remove this species to make sure the next one gets loaded
        rm(ThisSpecies)
      } # species.ss
      
      # Save GroupAbundance to the RichnessList for later use
      AbundanceList[[group.hh]] <- GroupAbundance
      # Save individual group
      save(
        GroupAbundance,
        file = paste0(
          path,
          'Output/9_CombinedAbundance/-',
          names(SpeciesLists)[group.hh],
          '-Abundance.RData'
        )
      )
      # Save incrementally in case process is interrupted
      save(AbundanceList,
           file = paste0(path, 'Output/9_CombinedAbundance/AbundanceList.RData'))
    } # group.hh
    # load(file = paste0(path, 'Output/9_CombinedAbundance/AbundanceList.RData'))
    
  } # GroupTraitAbundance
  
  ## Group proportions ####
  if ( GroupTraitHPrime == TRUE ) {
    # Calculate diversity proportions for each species
    # Need the group abundance first.
    # Do it for each group
    if (exists('AbundanceList') == FALSE) {
      load(file = paste0(path, 'Output/9_CombinedAbundance/AbundanceList.RData'))
    }
    AbundanceList %>% glimpse
    
    # Make an empty list to hold the group diversity arrays
    DiversityList <- vector('list', length = length(SpeciesLists))
    names(DiversityList) <- names(SpeciesLists)
    
    
    for (group.hh in 1:length(AbundanceList)) {
      # Group files
      ThisGroup <- SpeciesLists[[group.hh]]
      GroupFiles <- SppFiles %>%
        filter(SpeciesNames %in% ThisGroup)
      
      
      print(paste(
        'Calculating proportions for group',
        group.hh,
        'of',
        length(AbundanceList)
      ))
      # Total abundance
      TotalAbundance = AbundanceList[group.hh]
      
      # Do it for each species
      for (species.ss in 1:nrow(GroupFiles)) {
        ThisSpecies <- extractorRData(file = GroupFiles$FullPath[species.ss], object = 'SppArray')
        
        SppProportion <- array(dim = dim(ThisSpecies))
        # Calculate the proportion
        for (number.nn in 1:length(ThisSpecies)) {
          SppProportion[number.nn] = ThisSpecies[number.nn] / TotalAbundance[[1]][number.nn]
          try(if (SppProportion[number.nn] > 1)
            stop(paste(
              'Proportion greater than one detected at location',
              number.nn
            )))
        }
        
        # Proportion * ln proportion
        # Sum these all up to get H'
        pi_x_ln_pi <-
          SppProportion *
          log(SppProportion)
        
        if (species.ss == 1) {
          sum_pi_x_ln_pi <- pi_x_ln_pi
        } else {
          for (number.nn in 1:length(pi_x_ln_pi)) {
            if (is.na(pi_x_ln_pi[number.nn]) == TRUE) {
              sum_pi_x_ln_pi[number.nn] <- sum_pi_x_ln_pi[number.nn]
            } else {
              if (is.na(sum_pi_x_ln_pi[number.nn]) == TRUE) {
                sum_pi_x_ln_pi[number.nn] <- pi_x_ln_pi[number.nn]
              } else {
                sum_pi_x_ln_pi[number.nn] <-
                  sum_pi_x_ln_pi[number.nn] +
                  pi_x_ln_pi[number.nn]
              }
            }
            
          }
        }
        
        
        # Remove this species to make sure the next one gets loaded
        rm(ThisSpecies)
      } # species.ss
      # Calculate group diversity
      DiversityList[[group.hh]] <- -1 * (sum_pi_x_ln_pi)
      # Save incrementally in case it stops
      save(DiversityList,
           file = paste0(path, 'Output/11_CombinedHprime/DiversityList.RData'))
      load(paste0(path, 'Output/11_CombinedHprime/DiversityList.RData'))
    } # group.hh
    
  } # GroupTraitHPrime
  
  # Replace NA's with 0
  # See NA's
  # DiversityList %>% 
  #   map(~ sum(is.na(.)))
  # The lowest possible value is 0
  # DiversityList %>% 
  #   lapply(FUN = min, na.rm = TRUE) 
  # Analyze Trait Diversity ####
  if ( AnalyzeTraitDiversity == TRUE ) {
    # Load diversity
    if (exists('DiversityList') == FALSE) {
      load(paste0(path, 'Output/11_CombinedHprime/DiversityList.RData'))
    }
    
    # Replace NA's with 0's in diversity metrics
    # Note: Nothing will happen if there aren't NA's
    # Note that there are NA's in the beginning
    DiversityList %>%
      purrr::map(~ sum(is.na(.)))
    # Replace the NA's with 0, the lowest possible value
    # . indicates the elements of DiversityList that are passed to map
    DiversityList %<>% 
      purrr::map(~ ifelse(is.na(.), 0, .))
    # Check that there are no longer NA's
    DiversityList %>%
      purrr::map(~ sum(is.na(.)))
    
    # Trait comparisons
    ComparisonList <- list(c('Migratory', 'Residents') %>% sort,
                           c('Generalists', 'Specialists') %>% sort,
                           c('Univoltine', 'Multivoltine') %>% sort,
                           c('Rare', 'Common') %>% sort,
                           c("OverwinterEgg", "OverwinterLarva",
                             "OverwinterAdult", "OverwinterPupa",
                             "OverwinterPL") %>% sort)
    
    # Format output
    OutputList <- vector('list', length = length(ComparisonList))
    post.hoc.df.div <- data.frame(
      trend.est = numeric(),
      trend.prob.pos = numeric(),
      trend.025 = numeric(),
      trend.25 = numeric(),
      trend.75 = numeric(),
      trend.975 = numeric(),
      Name = character()
    )
    
    set.seed(1249)
    # Run analyses
    for ( ll.compare in 1:length(ComparisonList) ) {
      OutputList[[ll.compare]] <- ButterflyDiversity( 
        WhichLists = ComparisonList[[ll.compare]]
      )
      post.hoc.df.div %<>% 
        add_row(OutputList[[ll.compare]]$SummaryStats)
    }
    pdf('Output/11_CombinedHprime/DiversityEvalPlots.pdf', width = 8, height = 6)
    for( nn.list in 1:length(OutputList) ) {
      y.mat <- matrix(OutputList[[nn.list]]$Model$y,
                      nrow = 1000, byrow = TRUE)
      Thousand <- sample(1:16768000, size = 10000)
      hist(y.mat[Thousand], main = ComparisonList[[nn.list]])
      plot(y.mat[Thousand] - OutputList[[nn.list]]$Model$y.hat.samples[Thousand] ~
             OutputList[[nn.list]]$Model$y.hat.samples[Thousand])
      plot(OutputList[[nn.list]]$Model$beta.samples, density = FALSE)
      plot(OutputList[[nn.list]]$Model$beta.star.samples[,1:6], density = FALSE)
    }
    dev.off()
    # OL <- list(OutputList = OutputList[[1]],
    #      ComparisonList = ComparisonList)
    # save(OL, file = 'Output/11_CombinedHprime/MigResForJeff.RData')
    # load(file = 'Output/11_CombinedHprime/MigResForJeff.RData')
    
    post.hoc.df.div %<>% 
      mutate(MeanCompare = ifelse(Name %>% str_detect('-'), 
                                  'Compare', 'Mean'),
             trend.prob.neg = 1 - trend.prob.pos,
             Metric = 'Diversity')
    
    # Save to file
    save(post.hoc.df.div,
         file = paste0(path, 'Output/11_CombinedHprime/CommunityDiversity.RData'))
    # load(file = paste0(path, 'Output/11_CombinedHprime/CommunityDiversity.RData'))
    
    NList <- SpeciesLists %>% lapply(n_distinct) %>% as_tibble %>% t
    Ns <- tibble(Name = rownames(NList),
                 n = NList[,1])
    
    GroupTable <- post.hoc.df.div %>% 
      filter(MeanCompare == 'Mean') %>% 
      left_join(Ns) %>% 
      mutate(Name = factor(Name, levels = c("Migratory", "Residents", 
                                            "Generalists", "Specialists",
                                            "Multivoltine", "Univoltine",
                                            "Common", "Rare",
                                            "OverwinterEgg", "OverwinterLarva", 
                                            "OverwinterPL", "OverwinterPupa",
                                            "OverwinterAdult"))) %>% 
      # mutate(percent.trend = trend.est * 100,
      #        percent.025 = trend.025 * 100,
      #        percent.975 = trend.975 * 100) %>% 
      mutate(CI95 = paste0(
        '(', trend.025 %>% round(2), ' to ', trend.975 %>% round(2), ')')
      ) %>% 
      select(Name, n, trend.est, CI95, trend.prob.neg) %>% 
      arrange(Name)
    
    # write.csv(GroupTable, file = 'Output/11_CombinedHprime/GroupMeanDiv.csv',
                # row.names = FALSE)
    
    
    post.hoc.df.div %>% 
      filter(MeanCompare == 'Compare') %>% 
      select(Name, trend.est, trend.025, trend.975, trend.prob.neg)
    # Eggs, PL, Pupa, Adult, Larva
    # Eggs and PL declining others increase
    # Pupa and Adult increasing less than Larva
    # Eggs declining more
    # Larvae declining more
    # PL not different than adult, likely not changing but not sure
    # Pupa not different than adult, maybe not changing
    
    #
  }
  
} # CreateTraitDiversity


# Check for NA's
DL <- DiversityList %>% 
  lapply(FUN = is.na) %>% lapply(sum)
DiversityList %>% 
  map(is.na) %>% map(sum)
DiversityList %>% 
  map(~ sum(is.na(.)))
DiversityList %>% 
  lapply(FUN = max, na.rm = TRUE) 
