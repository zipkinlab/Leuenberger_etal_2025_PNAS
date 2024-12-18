#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Species richness - all species all groups together 
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
n.weeks <- 23
# week.0.linear <- data.list$covs$week.cov %>% scale %>% unique # All of my weeks unique values
# n.weeks <- length(week.0.linear) # Number of unique week values


# Toggles if needed
CalcCombinedRichness <- FALSE

PostHocRichness <- FALSE

RichnessByTrait <- FALSE
CalcSppy0 <- FALSE
CreateTraitRichness <- FALSE
AnalyzeTraitRichness <- TRUE

# Set seed
set.seed(seed = 761)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculate Richness combined among groups ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if( CalcCombinedRichness == TRUE ) { 
  print('Combining richness among groups')
  # Pull the file path for the total abundance folder
  Path <- paste0(path, 'Output/6_richness/')
  # All of the possible files in one dataframe
  Files <- expand.grid(ID = IDs, 
                       Piece = Pieces) %>% 
    mutate(PossibleFiles = paste0(Path, ID, 'richness-piece-', Piece, '.RData')) %>% 
    left_join(IDShort)
  RichnessAll <- array(0, dim = c(n.samples, J.0, n.years))
  # Loop through each piece.
  # TotalAbundance objects are small because all species are combined (1 less dimension)
  # But we'll need to run the next calculations on y.0.AUC's, which are larger
  # So let's keep it split into pieces. Smaller pieces might queue faster on HPCC
  RichnessCombined <- 0; RichnessSeparated <- 0  # Placeholders to remove in loop
  for (ll in 1:n.pieces) {
    print(paste("Currently on combining richness piece ", ll, " out of ", n.pieces, sep = '')) 
    # Housekeeping - removing and creating objects
    rm(RichnessCombined, RichnessSeparated)
    gc()  # Returns memory to the system after large object removal
    # Just grab the files for a given piece
    FilePieces <- Files %>% 
      filter(Piece == ll)
    # Which set of counties are we looking at (numbers/index)
    J.curr <- length(vals[[ll]])
    # Create the empty array to store everything in
    RichnessCombined <- array(0, dim = c(n.samples, J.curr, n.years))
    # Pull out each of the functional groups
    n.groups <- n_distinct(FilePieces$Short)
    groups <- unique(FilePieces$Short)
    RichnessSeparated <- array(0, dim = c(n.samples, J.curr, n.years, n.groups))
    for (gg in 1:n.groups) {
      RichnessSeparated[, , , gg] <-
        extractorRData(file = FilePieces$PossibleFiles[FilePieces$Short == groups[gg]],
                       object = 'richness')
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
          RichnessCombined[kk, jj, tt] <- sum(RichnessSeparated[kk, jj, tt, ])
        } # tt
      } # jj
    } # kk
    # Combine richness values into one array instead of in pieces
    RichnessAll[ , vals[[ll]], ] <- RichnessCombined
  } # ll
    # Save combined richness for analyses later
    save(RichnessAll, 
         file = 'Output/10_CombinedRichness/CombinedRichness.RData')
} # Combined Richness

# TODO: Note that Migrant has NA in iteration 1000 due to iterations not being a multiple of thin

## Running post hoc regression ####
if(PostHocRichness == TRUE) {
  print('Running Richness post hoc LMs')
  # Load file
  if(exists('RichnessAll') == FALSE) {
    RichnessAll <- extractorRData(file = paste0(path, 'Output/10_CombinedRichness/CombinedRichness.RData'),
                             object = 'RichnessAll')
  }
  
  DateCreated <- file.info(paste0(path, 'Output/10_CombinedRichness/CombinedRichness.RData'))
  # When I fixed the year order problem (Thanks Jeff!!!)
  # file.info('Code/R/PostProcess_MsAbund.R')$mtime
  if(DateCreated$mtime < "2024-05-30 15:00:54 EDT") {
    year.0.wrong <- c(3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
                      21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 1, 2, 4, 6, 32)
    RichnessAll <- RichnessAll[ , , order(year.0.wrong)]
  }
  
  # Transition paragraph (Results -> Discussion)
  # First five year median compared to last five year median
  RichnessAll[ , , 1:5] %>% median
  RichnessAll[ , , 28:32] %>% median
  
  # Create empty data frame to store post hoc results
  post.hoc.df.richness <- data.frame(
    trend.est = NA,
    int.est = NA,
    trend.prob.pos = NA,
    trend.low = NA,
    trend.high = NA
  )
  J <- ncol(RichnessAll)
  n.samples <- nrow(RichnessAll)
  unique.years <- 1:dim(RichnessAll)[3]
  plot.years <- unique.years
  years <- rep(1:length(unique.years), each = J) 
  n.years <- length(unique(years))
  sites <- rep(1:J, times = n.years)
  # Format data for postHocLM
  y <- matrix(RichnessAll, nrow = n.samples, ncol = J * n.years)
  new.data.list <- list(y = log(y), 
                        covs = data.frame(years, sites)) 
  out.posthoc <- postHocLM(formula = ~ years + (1 | sites), 
                           data = new.data.list, 
                           n.chains = 1, 
                           verbose = TRUE)
  out.posthoc %>% summary
  # plot(out.posthoc$beta.samples, density = FALSE)
  # plot(out.posthoc$beta.star.samples, density = FALSE)
  
  # RichnessAll[1, 1:10, 1]
  # y[1, 1:10]
  # Save results and generate summary figure 
  prob.pos <- apply(out.posthoc$beta.samples, 2, function(a) mean(a > 0))
  post.hoc.df.richness$trend.prob.pos <- prob.pos[2]
  # y.log.med.full <- matrix(apply(log(y), 2, median), J, n.years)
  # abund.indx.log.site.avg <- apply(log(abund.indx), c(1, 3), mean)
  # abund.indx.log.quants <- apply(abund.indx.log.site.avg, 2, quantile, probs = c(0.025, 0.5, 0.975))
  # If not log transformed, use these next lines
  # richness.log.quants <- apply(RichnessAll[, , ], 3, quantile, probs = c(0.025, 0.5, 0.975))
  # beta.quants <- apply(out.posthoc$beta.samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
  # Log transformed code
  richness.log.quants <- apply(log(RichnessAll[, , ]), 3, quantile, probs = c(0.025, 0.5, 0.975))
  beta.quants <- apply(exp(out.posthoc$beta.samples), 2, quantile, probs = c(0.025, 0.5, 0.975))
  beta.log.quants <- apply(out.posthoc$beta.samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
  post.hoc.df.richness$trend.est <- beta.log.quants[2, 2]
  post.hoc.df.richness$int.est <- beta.log.quants[2, 1]
  post.hoc.df.richness$trend.low <- beta.log.quants[1, 2]
  post.hoc.df.richness$trend.high <- beta.log.quants[3, 2]
  unique.years <- unique(out.posthoc$X[, 'years'])
  plot.min <- min(richness.log.quants)
  plot.max <- max(richness.log.quants)
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
  Trend <- post.hoc.df.richness$int.est + post.hoc.df.richness$trend.est * plot.vals
  
  # TODO: Fix and make line the mean
  Lines <- data.frame(year = plot.vals,
                      med = Trend)
  # sp.name.plot <- simpleCap(str_replace(curr.sp.name, '-', ' '))
  plot.df <- data.frame(med = richness.log.quants[2, ], 
                        low = richness.log.quants[1, ],
                        high = richness.log.quants[3, ], 
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
    # labs(x = 'Year', y = 'Richness')
    labs(x = 'Year', y = 'Richness (log scale)')
  # Add stats and species name
  # my.plot <- my.plot +
  #   ggtitle(label = paste("Richness: P(trend < 0) = ", round(1 - prob.pos[2], 2), 
  #                         "\n% Change/Year: ", 
  #                         round((beta.quants[2, 2] - 1) * 100, 2), " (", 
  #                         round((beta.quants[1, 2] - 1) * 100, 2), ", ", 
  #                         round((beta.quants[3, 2] - 1) * 100, 2), ")", sep = ''))
  my.plot <- my.plot +
    ggtitle(label = paste("Richness: P(trend < 0) = ", round(1 - prob.pos[2], 2), 
                          "\n% Change/Year: ", 
                          round((beta.quants[2, 2]), 2), " (", 
                          round((beta.quants[1, 2]), 2), ", ", 
                          round((beta.quants[3, 2]), 2), ")", sep = ''))
  # Add trendline
  my.plot <- my.plot + 
    geom_line(data = Lines, linewidth = 1)
  # lines(plot.vals, tmp, col = 'black', lty = 1, lwd = 5)
  
  # pdf(paste('Output/trend-figures/richness-trend.pdf', sep = ''),
  # width = 7, height = 7)
  print(my.plot)
  # dev.off()
  
  save(out.posthoc, 
       file = paste('Output/trend-figures/richness-posthoc-for-figure.rda',
                    sep = ''))
  
}  

# load('Output/trend-figures/richness-posthoc-for-figure.rda')
# load('Y:/ButterflyChapterOne/Output/trend-figures/richness-posthoc-for-figure.rda')
  
RichnessLM <- ButterflyPostHoc(Data = RichnessAll[,,], 
                                Logged = TRUE)
RichnessTrendPlot <- ButterflyTrendPlot(
  Data = RichnessAll[,,],
  Model = RichnessLM,
  Logged = TRUE,
  Name = 'Richness',
  ylabel = 'Richness (log scale)'#,
  # SaveFig = FALSE
)

if ( RichnessByTrait == TRUE ) {
  ## Reformat post-processing arrays ####
  if ( CalcSppy0 == TRUE ) {
    print('Formatting y.0 arrays to be species-specific')
    # Grab the species-specific predicted abundance file paths
    Path <- paste0(path, 'Output/8_y0AUC/')
    # Path for each file
    Files <- expand.grid(ID = IDs,
                         Piece = Pieces) %>%
      mutate(PossibleFiles = paste0(Path, ID, 'y0-AUC-piece', Piece, '.RData'),
             NSpp = str_extract(ID, '\\d{2}')) %>%
      left_join(IDShort)
    # Group information
    n.groups <- n_distinct(Files$Short)
    groups <- unique(Files$Short)
    for ( group.gg in 1:n.groups ) {
      FilePieces <- Files %>% 
        filter(Short == groups[group.gg])
      SppAbundance <-  array(0, dim = c(n.samples, FilePieces$NSpp[1], J.0, n.years))
      print(paste('Starting pulling in pieces for group', group.gg, 'of', n.groups))
      for ( piece.ll in 1:n.pieces ) {
        # Just the files from this piece
        File <- FilePieces %>%
          filter(Piece == piece.ll)
        # Which set of counties are we looking at (numbers/index)
        J.curr <- vals[[piece.ll]]
        SppAbundance[ , , J.curr, ]<-
          extractorRData(file = File$PossibleFiles,
                         object = 'y.0.AUC')
      } # piece.ll
      print(paste('Saving individual species for group', group.gg,
                  'of', n.groups))      
      GetNames <- 
        extractorRData(file = File$PossibleFiles[1],
                       object = 'y.0.AUC') %>% dimnames
      dimnames(SppAbundance)[[2]] <- GetNames[[2]]
      for ( species.ss in 1:dim(SppAbundance)[2] ) {
        SppName <- dimnames(SppAbundance)[[2]][species.ss]
        SppArray <- SppAbundance[ , SppName, , ] 
        save(SppArray, 
             file = paste0(path, 'Output/8_y0AUC/', 
                           SppName, '-y0-AUC.RData'))
      } # species.ss
    } # group.gg
  } # CalcSppy0
  
  # Calculate group richness ####
  if ( CreateTraitRichness == TRUE ) {
    # Make an empty list to hold the richness arrays
    RichnessList <- vector('list', length = length(SpeciesLists))
    names(RichnessList) <- names(SpeciesLists)
    
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
    
    # Calculate richness for each group
    for ( group.rr in 1:length(RichnessList) ) {
      print(paste('Adding species for group', group.rr, 'of', length(RichnessList)))
      # Group files
      ThisGroup <- SpeciesLists[[group.rr]]
      GroupFiles <- SppFiles %>% 
        filter(SpeciesNames %in% ThisGroup)
      
      GroupRichness <- array(0, dim = c(n.samples, J.0, n.years))
      # Add each of the species together
      for ( species.ss in 1:nrow(GroupFiles) ) {
        ThisSpecies <- extractorRData(file = GroupFiles$FullPath[species.ss],
                                      object = 'SppArray')
        ThisSpecies <- ifelse(ThisSpecies > 0, 1, 0)
        
        # Add to the Group so far
        GroupRichness <- GroupRichness + ThisSpecies
        
        # Remove this species to make sure the next one gets loaded
        rm(ThisSpecies)
      } # species.ss
      
      # Save GroupRichness to the RichnessList for later use
      RichnessList[[group.rr]] <- GroupRichness
      # Save incrementally in case process is interrupted
      # save(RichnessList, file = paste0(path, 'Output/10_CombinedRichness/RichnessList.RData'))
    } # group.rr
    # load(file = paste0(path, 'Output/10_CombinedRichness/RichnessList.RData'))
    
  } # CreateTraitRichness
  
  # Analyze trait richness ####
  if ( AnalyzeTraitRichness == TRUE ) {
    ComparisonList <- list(c('Migratory', 'Residents') %>% sort,
                           c('Generalists', 'Specialists') %>% sort,
                           c('Univoltine', 'Multivoltine') %>% sort,
                           c('Rare', 'Common') %>% sort,
                           c("OverwinterEgg", "OverwinterLarva",
                             "OverwinterAdult", "OverwinterPupa",
                             "OverwinterPL") %>% sort)
    # To check for homoscedasticity in single group models
    # These groups have heteroscedastic residuals due to unequal variances 
    # ComparisonList <- list('Migratory', 'Residents',
    #                        'Univoltine', 'Multivoltine')
    OutputList <- vector('list', length = length(ComparisonList))
    post.hoc.df.rich <- data.frame(
      trend.est = numeric(),
      trend.prob.pos = numeric(),
      trend.025 = numeric(),
      trend.25 = numeric(),
      trend.75 = numeric(),
      trend.975 = numeric(),
      Name = character()
    )
    
    if( exists('RichnessList') == FALSE ) {
      load(file = paste0(path, 'Output/10_CombinedRichness/RichnessList.RData'))
    }
    
    set.seed(1249)
    for ( ll.compare in 1:length(ComparisonList) ) {
      OutputList[[ll.compare]] <- ButterflyRichness(
        WhichLists = ComparisonList[[ll.compare]]
      )
      post.hoc.df.rich %<>% 
        add_row(OutputList[[ll.compare]]$SummaryStats)
    }
    pdf('Output/10_CombinedRichness/RichnessEvalPlots.pdf', width = 8, height = 6)
    for( nn.list in 1:length(OutputList) ) {
      plot(OutputList[[nn.list]]$Model$beta.samples, density = FALSE)
      plot(OutputList[[nn.list]]$Model$beta.star.samples[,1:6], density = FALSE)
    }
    dev.off()
    
    Thousand <- 1:100
    plot(OutputList[[3]]$Model$y.hat.samples[Thousand] ~
           OutputList[[3]]$Model$y[Thousand])
    NN <- sample(1:1000, size = 20)
    y.hat <- unlist(OutputList[[1]]$Model$y.hat.samples)
    
    pdf('Output/10_CombinedRichness/RichnessEvalPlots.pdf', width = 8, height = 6)
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
    # save(OL, file = 'Output/10_CombinedRichness/MigResForJeff.RData')
    # load(file = 'Output/10_CombinedRichness/MigResForJeff.RData')
    y.mat <- matrix(OL$OutputList$Model$y,
                    nrow = 1000, byrow = TRUE)
    Thousand <- sample(1:16768000, size = 10000)
    hist(y.mat[Thousand], main = ComparisonList[[1]])
    plot(y.mat[Thousand] - OL$OutputList$Model$y.hat.samples[Thousand] ~
           OL$OutputList$Model$y.hat.samples[Thousand])
    plot(OL$OutputList$Model$beta.samples, density = FALSE)
    plot(OL$OutputList$Model$beta.star.samples[,1:6], density = FALSE)
    
    
    cbind(OutputList[[1]]$Model$y.hat.samples[NN, 1],
          OutputList[[1]]$Model$y.hat.samples[1, NN],
          OutputList[[1]]$Model$y.hat.samples[NN],
          y.mat[NN],
          OutputList[[1]]$Model$y[NN])
    
    post.hoc.df.rich %<>% 
      mutate(MeanCompare = ifelse(Name %>% str_detect('-'), 
                                  'Compare', 'Mean'),
             trend.prob.neg = 1 - trend.prob.pos,
             Metric = 'Richness')
    
    # Save to file
    # save(post.hoc.df.rich,
    #      file = paste0(path, 'Output/10_CombinedRichness/CommunityRichness.RData'))
    # load(file = paste0(path, 'Output/10_CombinedRichness/CommunityRichness.RData'))
    
    NList <- SpeciesLists %>% lapply(n_distinct) %>% as_tibble %>% t
    Ns <- tibble(Name = rownames(NList),
                 n = NList[,1])
    
    GroupTable <- post.hoc.df.rich %>% 
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
    
    # DON'T USE THIS ONE.
    # NOW CODED IN MULTIMETRICFIGURES.R
    # write.csv(GroupTable, file = 'Output/10_CombinedRichness/GroupMeanRich.csv',
    #             row.names = FALSE)
    
    
    
    # Multivoltine and Univoltine declining, not different
    # Generalists and Specialists declining, not different
    # Migrants constant, Residents declining (different)
    # Common declining, Rare declining even more
    # Adults declining
      # Eggs declining more
      # Larvae declining more
      # PL not different than adult, likely not changing but not sure
      # Pupa not different than adult, maybe not changing
    
    
  } # AnalyzeTraitRichness
  
} # RichnessByTrait


# Convert spp/year into percents ####
RichnessStartEnd <- tibble(Name = names(RichnessList),
                           Start = NA,
                           End = NA)
for ( nn.list in 1:length(RichnessList) ) {
  RichnessStartEnd[nn.list, 'Start'] <- 
    RichnessList[[nn.list]][ , , 1:5] %>% median
  RichnessStartEnd[nn.list, 'End'] <- 
    RichnessList[[nn.list]][ , , 28:32] %>% median
}
post.hoc.df.rich %>% 
  filter(MeanCompare == 'Mean') %>% 
  mutate(TotalDecline = trend.est * 32) %>% 
  left_join(RichnessStartEnd) %>% 
  select(Name, trend.est, TotalDecline, Start) %>% 
  mutate(#PropMedian = End / Start,
         PropNum = (Start + TotalDecline) / Start,
         PercentLost = (1-PropNum) * 100,
         AnnualPercLost = trend.est/Start * 100) %>% 
  arrange(PropNum)


RichnessList %>% lapply(summary)

Zeros <- data.frame(Name = names(RichnessList),
                    NZeros = 0,
                    PropZeros = 0)

for ( list.ll in 1:length(RichnessList) ) {
  print(names(RichnessList)[list.ll])
  Total <- RichnessList[[list.ll]] %>% length
  N = sum(RichnessList[[list.ll]] == 0)
  Zeros$NZeros[list.ll] = N
  Zeros$PropZeros[list.ll] = N / Total
}
Zeros
# Zeros %>% write.csv('Output/10_CombinedRichness/Zeros.csv', 
#                     row.names = FALSE)

# Map of change in richness ####
# Not useful because no space * time component in model
RichnessAll %>% dim
# Median of first 5 years and last 5 years?
FirstLastYears <- 5
# Grab the total number of years
LastYr <- dim(RichnessAll)[3]
# Calculate the median for each county 
#   across all the 1000 iterations in the first and last years
FirstYears <- RichnessAll[, , 1:FirstLastYears] %>% 
  apply(MARGIN = 2, FUN = median)
LastYears <-  RichnessAll[, , (LastYr - FirstLastYears + 1):LastYr] %>%
  apply(MARGIN = 2, FUN = median)
# Calculate the change in richness
Diff <- LastYears - FirstYears
# Median richness by county
CountyMedians <- tibble(FirstYears = FirstYears,
                        LastYears = LastYears,
                        Difference = Diff)

# ggplot map code
# Load metadata
load(file = 'Data/CountyKey.RData')
load(file = 'Output/DataAndDimensions.RData')

# County and shapefile type data
# Which states are included
States <- c('iowa', 'wisconsin', 'michigan', 'indiana', 
            'illinois', 'ohio', 'missouri', 'minnesota')
# Which counties are included
Counties <- CountyKey$state.county.name
# Grab the polygon info for the above states
mw.states <- map_data('state') %>%
  filter(region %in% States)
# Grab every county in the above states
mw.counties <- map_data('county') %>% 
  filter(region %in% States) %>% 
  # Some counties in different states have the same name
  # Make a new column that's unique counties. 
  # Also matches the format of state.county.name
  mutate(state.county.name = paste(region, subregion, sep = ',')) 
mw.counties %<>% 
  mutate(state.county.name = 
           case_when(
             state.county.name == 'illinois,de kalb' ~ 'illinois,dekalb',
             state.county.name == 'illinois,du page' ~ 'illinois,dupage',
             state.county.name == 'illinois,la salle' ~ 'illinois,lasalle',
             state.county.name == 'illinois,st clair' ~ 'illinois,st. clair',
             state.county.name == 'indiana,la porte' ~ 'indiana,laporte',
             state.county.name == 'indiana,de kalb' ~ 'indiana,dekalb',
             state.county.name == 'iowa,obrien' ~ "iowa,o'brien",
             state.county.name == 'michigan,st clair' ~ 'michigan,st. clair',
             state.county.name == 'minnesota,st louis' ~ 'minnesota,st. louis',
             state.county.name == 'missouri,st charles' ~ 'missouri,st. charles',
             state.county.name == 'missouri,st clair' ~ 'missouri,st. clair',
             # state.county.name == 'west virginia,cabell' ~ 'west virginia,cabell',
             .default = state.county.name
           ))
# TODO: why is there a West Virginia county in here?
Counties[!Counties %in% mw.counties$state.county.name]
mw.counties %>%
  select(region, subregion, state.county.name) %>% 
  distinct %>% 
  filter(region == 'west virginia') 

# sf objects
# ggplot settings
gpStates <- geom_polygon(
  data = mw.states,
  aes(x = long, y = lat, group = region),
  col = 'black',
  fill = 'white')
gpAllCounties <- geom_polygon(
  data = mw.counties,
  aes(x = long, y = lat, group = state.county.name),
  col = 'black', fill = 'white') 
tbw <- theme_bw(base_size = 18)
tv <- theme_void(base_size = 18)
tlt <- theme(legend.text = element_text(size = 12),
             legend.title = element_text(size = 14))

CountyKeyMedians <- CountyKey %>% arrange(county.ind) %>% 
  bind_cols(CountyMedians)
mw.counties %<>% left_join(CountyKeyMedians) %>% 
  filter(!is.na(Difference))
ggplot() +
  gpStates + 
  gpAllCounties +
  geom_polygon(data = mw.counties,
               aes(x = long, y = lat,
                   group = state.county.name,
                   # fill = FirstYears),
                   # fill = LastYears),
                   fill = Difference),
               col = 'black') +
  scale_fill_gradient2(low = '#E76254', 
                       high = '#1E466E',
                       midpoint = 0) +
  tv +
  labs(fill = "Change in\nrichness")
  # MetBrewer::scale_fill_met_c('Hiroshige', 
                              # midpoint = 0) +
  # MetBrewer::scale_fill_met_c('OKeeffe1') +
  # tv
  
ncol <- 100

## Make a vector with n colors
cols <- RColorBrewer:::brewer.pal(11,"PuOr")  # OR c("purple","white","orange")  
rampcols <- colorRampPalette(colors = cols, space="Lab")(ncol)
rampcols[(n/2) + 1] <- rgb(t(col2rgb("green")), maxColorValue=256) 

## Make a vector with n+1 breaks
rampbreaks <- seq(0, 100, length.out = ncol+1)

## Try it out
heatmap(data_matrix, Rowv = NA, Colv = NA, scale="none",
        col = rampcols, breaks = rampbreaks)