#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Abundance calculations and trends
# Use this code after Code/R/PostProcess_MsAbund.R
# Wendy Leuenberger
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load packages ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Brackets run the entire section with one entry. 
# Added for convenience
{ 
  # Load packages
  library(tidyverse)
  library(magrittr)
  library(spAbundance)
  library(pracma)
  library(spOccupancy)
  library(coda)
  library(abind)
  library(WendyFunctions)
  library(MetBrewer)
  library(ggpp)  # Text annotations on plot
  library(ggpubr)  # Multipanel figures
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ggplot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

{
  # ggplot settings
  gp <- geom_point()
  tbw <-  theme_bw(base_size = 18)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# FUNCTIONS ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Updated WendyFunctions

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set-up ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# CHANGE HERE AS NEEDED
# I'm using a unique identifier to link all files from particular runs 
# or chains from particular runs with identical data
# Set those ID's here and we'll use them throughout the post-processing
# The goal is to reduce the number of places where these need to be set to 
# reduce human error and mixing up results with species/group names
{
  # Run all of the starting stuff (not toggles)
  IDs <- c(
    '-ResidentSpecialistHesPap24Spp-',
    '-ResidentSpecialistLycNymPieRio23Spp-',
    '-ResidentGeneralistHesPap15Spp-',
    '-ResidentGeneralistLycNymPieRio19Spp-',
    '-ResidentUnivoltineHesPap11Spp-',
    '-ResidentUnivoltineLycNymPieRio23Spp-',
    '-Migratory21Spp-'
  )
  IDShort <- data.frame(
    ID = IDs,
    Short = c('RSHP', 'RSLNPR', 'RGHP', 'RGLNPR', 'RUHP', 'RULNPR', 'Migrant'),
    n.spp = c(25, 23, 15, 19, 11, 23, 20)
  )
  IDShort %<>%
    mutate(GroupFamily = str_extract(IDs, '[:alpha:]+'))
  
  # File path
  # Note: as of 9/9, the trends in local computer aren't updated
  # Use HPCC version
  RemoteHPCC <- FALSE
  RemoteHPCC <- TRUE
  if (RemoteHPCC == TRUE) {
    path = 'Y:/ButterflyChapterOne/'
  } else {
    path = ''
  }
  
  # Load object with dimensions
  load(file = paste0(path, 'Output/DataAndDimensions.RData'))
  # Load object with county information
  load(file = paste0(path, 'Data/CountyKey.RData'))
  
  # Source code for species traits
  source('Code/R/SpeciesFunctionalTraits.R')
}

# Toggles ####
{
  # Species toggles
  # Run the regressions/first plots with all the individual species
  ModelIndividualSpecies <- FALSE
  # Run post processing on species trends
  # Compare older trends with newer trends
  # Likely never needed
  CompareTrends <- FALSE
  # Make PDF of each species
  MakeSpeciesPDF <- FALSE
  # Make species-specific plots
  SpeciesPlots <- FALSE
  # Summary stats and processing of trends
  SpeciesTrendSummary <- TRUE
  # Line plot (in presentation and paper)
  SpeciesLinePlot <- TRUE
  # Needed for FannedTrends
  SpeciesFanPlot <- TRUE
}

# Species-specific parts needed for median group models
ModelGroups <- TRUE
if ( ModelGroups == TRUE ) {
  SpeciesTrendSummary <- TRUE
  SpeciesFanPlot <- TRUE
}
  
# TODO: Remove references to dataNew species names once all groups are rerun

# Add spp names to indices here because I hadn't saved them with the objects
# load('Data/FormatData.RData')
# Remove random seed that R automatically associates with a workspace
# Causes results to not be random if just loading
# rm(.Random.seed)

# CodesInput <- Groups %>%
#   filter(Code %in% unique(dataNew$species)) %>%
#   select(Code, GroupFamily) %>% 
#   arrange(GroupFamily, Code) %>% 
#   mutate(
#     GroupFamily = str_replace_all(
#       string = GroupFamily, pattern = ' ', replace = '')
#   )


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Number of iterations per group ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Chains = data.frame(GroupFamily = IDShort$GroupFamily,
                    # N chains manually taken from PostProcess_MsAbund.R
                    Nchains = c(5,  # RSHP 
                                4,  # RSLNPR
                                5,  # RGHP
                                4,  # RGLNPR
                                4,  # RUHP
                                5,  # RULNPR
                                5), # Migratory
                    Niter = c(300,  # RSHP 
                              200,  # RSLNPR
                              460,  # RGHP
                              420,  # RGLNPR
                              1300, # RUHP
                              300,  # RULNPR
                              250), # Migratory
                    Nruns = 2,
                    batch.length = 100
                    ) %>% 
  mutate(TotalIter = Nchains * Niter * batch.length * Nruns)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load all ####
# Requires too much RAM for local computer (11.7 GB for array?)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OneFile <- FALSE
if(OneFile == TRUE){
  print('Putting all species into one array')
  # Pull the file path for the total abundance folder
  Path <- file.path(getwd(), 'Output/5_abund-indx') %>% paste0('/')
  # All of the possible files in one dataframe
  Files <- tibble(ID = IDs) %>% 
    mutate(PossibleFiles = paste0(Path, ID, 'abundance-indx.RData')) %>% 
    left_join(IDShort)
  # Counters for where each species begins and ends in the resulting array
  # Start at 1, then 1 + the end of the last
  # End at the cumulative sum
  Files %<>% 
    mutate(SppStart = lag(cumsum(n.spp), default = 0) + 1,
           SppEnd = cumsum(n.spp))
  all.spp <- sum(Files$n.spp)
  AllSpecies <- array(0, dim = c(n.samples, all.spp, J.0, n.years))
  n.groups <- n_distinct(Files$Short)
  groups <- unique(Files$Short)
  
  for (gg in 1:n.groups) {
    Start = Files$SppStart[Files$Short == groups[gg]]
    End = Files$SppEnd[Files$Short == groups[gg]]
    abund.indx <-
      extractorRData(file = Files$PossibleFiles[Files$Short == groups[gg]],
                     object = 'abund.indx')
    AllSpecies[, Start:End , , ] <- abund.indx
    dimnames(AllSpecies)[[2]][Start:End] <- dimnames(abund.indx)[[2]]
    
  }
  save(AllSpecies, file = paste('Output/5_abund-indx/all-species-abundance-indx.RData',
                                sep = ''))
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Extract all files ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{
  # Extract species files
  Path <- paste0(path, 'Output/5_abund-indx')
  
  FileList <- list.files(path = Path, pattern = "-[[:upper:]]{6,}-*")
  ColSp <- list.files(path = Path, pattern = 'COL-SP')
  Satca1 <- list.files(path = Path, pattern = 'SATCA1')
  Satca3 <- list.files(path = Path, pattern = 'SATCA3')
  Files <- tibble(Path, FileList) %>%
    mutate(SpeciesNames =  str_extract(string = FileList, pattern = '[[:upper:]]{6,}')) %>%
    add_row(
      Path = Path,
      FileList = c(ColSp, Satca1, Satca3),
      SpeciesNames = c('COL-SP', 'SATCA1', 'SATCA3')
    ) %>%
    mutate(FullPath = paste0(Path, '/', FileList))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Species ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n.years <- 32
if (ModelIndividualSpecies == TRUE) {
  set.seed(seed = 48933)
  EachSpp <- nrow(Files)
  LMList <- vector(mode = 'list', length = EachSpp)
  PlotData <- data.frame(Median = NA,
                         Mean = NA,
                         year = rep(1:n.years, times = EachSpp),
                         Species = rep(Files$SpeciesNames, each = n.years))
  Trends <- data.frame(trend.est = NA, 
                       int.est = NA,
                       int.est.median = NA,
                       trend.prob.pos = NA, 
                       trend.025 = NA,
                       trend.25 = NA,
                       trend.75 = NA,
                       trend.975 = NA,
                       Species = Files$SpeciesNames)
  BetaYears <- vector('list', length = EachSpp)
  ii <- 1
  for (ii in 1:EachSpp) {
    print(paste('Analyzing species', ii, 'of', EachSpp, 'species'))
    # Remove previous Species.indx and free up space
    if(ii != 1) {
      rm(Species.indx)
      gc()
    }
    # Pull up file
    Species.indx <- extractorRData(Files$FullPath[ii],
                                   object = 'Species.indx')
    # Run the linear regression
    SpeciesLM <- ButterflyPostHoc(Data = Species.indx, 
                                  Logged = TRUE)
    # Save the beta year values for each iteration for the group trends
    BetaYears[[ii]] <- SpeciesLM$BetaYears
    names(BetaYears)[[ii]] <- Files$SpeciesNames[ii]
    # I didn't get names saved initially, here's how to load them. 
    # for(ii in 1:136){
    #   names(BetaYears)[ii] <- Files$SpeciesNames[ii]
    # }
    
    # Get trends only - not plotting anymore
    # ButterflyTrendPlot can also look at a rough plot.
    #   but it's kinda buggy
    SpeciesTrendPlot <- ButterflyTrendPlot(
      Data = Species.indx,
      Model = SpeciesLM$ModelOutput,
      Logged = TRUE,
      Name = Files$SpeciesNames[ii],
      RunFigCode = FALSE
    )
    
    # TryPlot <- SpeciesTrendPlot[[2]] %>%
    #   select(med, means, year) %>% 
    #   rename(Median = med,
    #          Mean = means) %>% 
    #   pivot_longer(-year, names_to = 'Type', values_to = 'Value')
    # SpeciesTrendPlot[[1]]
    # PlotLines <- tibble(Type = c('Median', 'Mean'),
    #                     Intercept = c(SpeciesTrendPlot[[1]]$int.est.median,
    #                                   SpeciesTrendPlot[[1]]$int.est),
    #                     Slope = SpeciesTrendPlot[[1]]$trend.est)
    # 
    # ggplot(TryPlot, aes(x = year, y = Value)) + 
    #   geom_point() +
    #   facet_wrap(~Type) +
    #   geom_abline(data = PlotLines, aes(slope = Slope, intercept = Intercept))
    
    # Save the model output
    # Too big, don't save
    # LMList[[ii]] <- SpeciesLM
    # Save the trends data
    # These are model outputs - based on the quantiles of beta.samples
    Trends[ii, 1:9] <- SpeciesTrendPlot[[1]]
    # Save the data for the plot
    # These are data (intermediate output) - based on log(Data) instead 
    # These have annual estimates rather than one value as for beta.samples
    Start <- (ii-1) * n.years + 1
    PlotData[Start : (Start + n.years - 1), ] <- 
      SpeciesTrendPlot[[2]][,c('med', 'means', 'year')] %>% 
      cbind(Species = Files$SpeciesNames[ii])
    
    
  } # ii species
  
  
  # save Trends, PlotData, SpeciesLM
  # save(LMList, file = paste0(path, 'Output/12_AbundanceTrends/AllIndividualSpp-LMs.RData'))
  save(BetaYears, file = paste0(path, 'Output/12_AbundanceTrends/AllIndividualSpp-BetaYears.RData'))
  save(PlotData, file = paste0(path, 'Output/12_AbundanceTrends/AllIndividualSpp-PlotData.RData'))
  save(Trends, file = paste0(path, 'Output/12_AbundanceTrends/AllIndividualSpp-Trends.RData'))
} # All individual species

  ### Compare trends ####
  if ( CompareTrends == TRUE ) {
    # Compare original one with runs that were extended by one week
    rm(Trends)
    load(
      paste0(
        path,
        'Output/12_AbundanceTrends/AllIndividualSpp-Trends_20240530.RData'
      )
    )
    OldTrends <- Trends
    rm(Trends)
    load(paste0(
      path,
      'Output/12_AbundanceTrends/AllIndividualSpp-Trends.RData'
    ))
    NewTrends <- Trends
    TrendCompare <- tibble(
      Species = OldTrends$Species,
      OldTrend = OldTrends$trend.est,
      NewTrend = NewTrends$trend.est
    ) %>%
      mutate(Difference = NewTrend - OldTrend)
    TrendCompare$Difference %>% hist
    # Similar but not identical - seems good!
    # Look at the 10 species with the biggest change (5 per direction)
    TrendCompare %>% slice_min(Difference, n = 5)
    TrendCompare %>% slice_max(Difference, n = 5)
  } # Compare trends
  
  if (SpeciesPlots == TRUE) {
    
    ### Species-specific plots ####
    # Make updated figures of mean, median, and trend for each species
    # trend.est in Trends is the logged value, so needs to be (exp(trend) - 1) * 100
    # for percent change
    load(paste0(
      path,
      'Output/12_AbundanceTrends/AllIndividualSpp-PlotData.RData'
    ))
    # PlotData$Species <- NULL
    TryPlot <- PlotData %>%
      rename(Year = year) %>%
      pivot_longer(!Year:Species, names_to = 'Type', values_to = 'Value')
    
    PlotLines <- Trends %>%
      select(trend.est, int.est, int.est.median, Species) %>%
      rename(Mean = int.est, Median = int.est.median) %>%
      pivot_longer(!c(trend.est, Species),
                   names_to = 'Type',
                   values_to = 'Intercept')
    
    # Calculate species in decreasing order for trends
    DecreasingTrend <- Trends %>%
      arrange(desc(trend.est)) %>%
      use_series(Species)
    
    PlotList <- vector(mode = 'list', length = length(Files$SpeciesNames))
    for (nn.species in 1:length(Files$SpeciesNames)) {
      TrendsSpp <- Trends %>%
        # order alphabetically
        # filter(Species == Files$SpeciesNames[nn.species])
        # order by decreasing trend
        filter(Species == DecreasingTrend[nn.species])
      PlotList[[nn.species]] <-  ggplot(TryPlot %>%
                                          filter(Species == TrendsSpp$Species),
                                        aes(x = Year, y = Value)) +
        gp +
        facet_wrap( ~ Type) +
        geom_abline(
          data = PlotLines %>%
            filter(Species == TrendsSpp$Species),
          aes(slope = trend.est, intercept = Intercept)
        ) +
        tbw +
        labs(y = 'Relative abundance estimate (log scale)') +
        ggtitle(
          label = paste(
            TrendsSpp$Species,
            ": P(trend < 0) = ",
            round(1 - TrendsSpp$trend.prob.pos, 2),
            "\n% Change/Year: ",
            round((exp(
              TrendsSpp$trend.est
            ) - 1) * 100, 2),
            " (",
            round((exp(
              TrendsSpp$trend.025
            ) - 1) * 100, 2),
            ", ",
            round((exp(
              TrendsSpp$trend.975
            ) - 1) * 100, 2),
            ")",
            sep = ''
          )
        )
    }
    
    if (MakeSpeciesPDF == TRUE) {
      pdf('Output/12_AbundanceTrends/SpeciesTrendPlots.pdf')
      for (nn.species in 1:length(Files$SpeciesNames)) {
        print(PlotList[[nn.species]])
      }
      dev.off()
    }
  } # Species-specific plots
    
    ## Summarize Trends and back transform ####
    if ( SpeciesTrendSummary == TRUE ) {
    
      # Remove Trends or else columns get added twice
      rm(Trends)
      if (exists('Trends') == FALSE) {
        load(paste0(
          path,
          'Output/12_AbundanceTrends/AllIndividualSpp-Trends.RData'
        ))
      }
      
    
    
      # Add the functional group data to the trends
      Trends %<>% left_join(Groups, join_by(Species == Code))
      
      # if(exists('PlotData') == FALSE){
      #   load(paste0(path, 'Output/12_AbundanceTrends/AllIndividualSpp-PlotData.RData'))
      # }
      Trends %<>%
        mutate(
          DiffZero = ifelse(
            (Trends$trend.025 < 0 & Trends$trend.975 < 0) |
              (Trends$trend.025 > 0 & Trends$trend.975 > 0),
            'Different',
            'Not different'
          ),
          ProbDecreasing = 1 - trend.prob.pos
        )
      
      Transform <- tibble(LogTrend = seq(-0.05, 0.05, by = 0.01))
      Transform %>%
        mutate(ExpTrend = exp(LogTrend),
               Percent = (ExpTrend - 1) * 100)
      # Number Increasing, Stable, Decreasing
      NumberTrends <- Trends %>%
        group_by(Group) %>%
        summarize(
          Increasing = sum(trend.prob.pos > 0.975),
          Decreasing = sum(trend.prob.pos <= 0.025),
          Stable = sum(DiffZero == 'Not different')
        )
      # All spp combined
      
      
      # Look at one species with exact CI's
      Trends %>%
        filter(trend.prob.pos == 0.025 |
                 trend.prob.pos == 0.975)
      # ATRARO should be considered decreasing (97.5 is -0.0002)
      Trends %>%
        summarize(
          Increasing = sum(trend.prob.pos > 0.975),
          Decreasing = sum(trend.prob.pos <= 0.025),
          Stable = sum(DiffZero == 'Not different')
        )
      # Uncomment for dput (to reorder list of species by abundance)
      # Trends %>%
      #   arrange(trend.est) %>%
      #   use_series(Species) %>%
      #   dput()
      Trends$Species %<>%
        factor(levels = c(
          "SPEYAPH", "CERPEG", "CHLOGOR", "NYMVAU", "NYMMIL", "PHYSEL", 
          "HESLEO", "EUPHPHA", "LYCAPHL", "SATACA", "EUPYCON", "BOLOSEL", 
          "LYCAHYL", "THYLIN", "POAHOB", "SATYEUR", "SATEDW", "BOLOBEL", 
          "POAMAS", "POLYPRO", "SATCA1", "NATIOL", "ERYICE", "POLITHE", 
          "HEMISO", "LYCAHEL", "SATTIT", "SATCA3", "AMBVIA", "FENTAR", 
          "SATLIP", "ERYMAR", "ENOANT", "SATYAPP", "POLIORI", "WALEGE", 
          "ACHLYC", "ATRARO", "VANVIR", "CYLGEM", "CALYCEC", "ERYHOR", 
          "ANATLOG", "SPEYCYB", "ERYBRI", "CHLONYC", "VANCAR", "SPEYATL", 
          "THOPYL", "PONPRO", "PIENAP", "COLCES", "NYMANT", "EUPYDIO", 
          "NASLHE", "EUPYVES", "COL-SP", "COLINT", "POLIMYS", "STAPHAY", 
          "AMBHEG", "POLYFAU", "PIERAP", "MEGICYM", "STRMEL", "LIMARC", 
          "DANPLE", "COETUL", "PAPCAN", "CALLAUG", "LYCAEPI", "ASTCLY", 
          "POLYCOM", "EVECOM", "THOBAT", "PHOLCAT", "ERYJUV", "HESSAS", 
          "PIEVIR", "EURLIS", "LYCMEL", "ASTCEL", "HERSOS", "PARRMAL", 
          "PROBYS", "EUPYBIM", "CALLGRY", "CHLOHAR", "BATPHI", "POAZAB", 
          "LYCADIO", "PHYBAT", "EURNIC", "PANOCO", "ATRYHIA", "VANATA", 
          "CALLHEN", "EPACLA", "CELLAD", "LYCIDA", "POAVIA", "POLIPEC", 
          "LIMART", "EURYMAR", "EUPTCLA", "PAPPOL", "ANCNUM", "PHYTHA", 
          "PAPGLA", "HESOTT", "ERYPER", "JUNCOE", "SPEYIDA", "HYLPHY", 
          "LIBCAR", "CALMUT", "LYCADOR", "CALBOR", "ERYBAP", "ANTHMID", 
          "CALLNIP", "POMVER", "PAPCRE", "GLALYG", "EUPYDUK", "CARTPAL", 
          "OARPOW", "ATACAM", "EUCOLY", "PAPTRO", "CALLPOL", "POLYINT", 
          "PHOESEN", "CALLIRU", "PYRCOM", "HESMET"
        ) %>% rev()
        )
      
      # Put on real scale and times 100 for figures
      Trends %<>%
        mutate(
          trend.025.real = (exp(trend.025) - 1) * 100,
          trend.est.real = (exp(trend.est) - 1) * 100,
          trend.975.real = (exp(trend.975) - 1) * 100
        )
    } # SpeciesTrendSummary

# Scratch for results section
Trends %>% 
  filter(DiffZero == 'Different') %>% 
  select(Species, Count, trend.est.real, trend.025.real, trend.975.real, ProbDecreasing) %>%
  # select(Species, trend.est.real, Family, MigrationVerdict, Voltinism, HostPlant) %>%
  arrange(trend.est.real)
Trends %>% 
  filter(DiffZero == 'Not different',
         trend.est.real < 0,
         # ProbDecreasing > 0.9) %>%
         ProbDecreasing >= 0.8, ProbDecreasing < 0.9) %>%
         # ProbDecreasing < 0.8) %>%
  select(Species, Count, trend.est.real, trend.025.real, trend.975.real, ProbDecreasing) %>%
  # select(Species, trend.est.real, Family, MigrationVerdict, Voltinism, HostPlant) %>%
  arrange(trend.est.real) 

# PercentTrendTable (Table S4) is in Results.Rmd
  

  ### Species Line Plots of trends ####
  if ( SpeciesLinePlot == TRUE ) {
    # Get y axis ranges
    PlotObject <- # Use this to get range
      ggplot(Trends, aes(x = Species, y = trend.est.real, color = Group)) +
      geom_point() +
      geom_linerange(aes(min = trend.025.real, max = trend.975.real))
    # If numbers change, will need a new y range for the y axis if removing stable ones first
    yrange <- layer_scales(PlotObject)$y$range$range
    
    ggplot(Trends, # ggplot(AllSpp,
           aes(x = Species, y = trend.est.real, color = Group)) +
      geom_point() +
      geom_linerange(aes(min = trend.025.real, max = trend.975.real)) +
      # For a blank plot, uncomment these lines and comment the two above
      # geom_point(alpha = 0) +
      # geom_linerange(aes(min = trend.025, max = trend.975), alpha = 0) +
      theme_bw(base_size = 24) +
      theme(
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none'
      ) +
      facet_grid( ~ Group, scales = 'free_x', space = 'free') +
      labs(y = 'Annual Percent Change (95% CI)') +
      geom_hline(yintercept = 0) +
      scale_color_manual(values = MetBrewer::met.brewer("Egypt", 4, direction = 1)) +
      # ggpp::geom_text_npc(data = NumberTrends,
      #                     aes(npcx = 'left', npcy = 'top',
      #                         label = paste0('Increasing: ', Increasing)),
      #                     size = 5) +
      # ggpp::geom_text_npc(data = NumberTrends,
      #                     aes(npcx = 'left', npcy = 'top',
      #                         label = paste0('Increasing: ', Increasing, '\n',
      #                                       'Decreasing: ', Decreasing)),
      #                     size = 5) +
      # ggpp::geom_text_npc(data = NumberTrends,
      #                     aes(npcx = 'left', npcy = 'top',
      #                         label = paste0('Increasing: ', Increasing, '\n',
      #                                       'Decreasing: ', Decreasing, '\n',
      #                                       'Stable: ', Stable)),
      #                     size = 5) +
      ylim(yrange[1], yrange[2])
    # ggsave('Output/trend-figures/Trends.pdf', height = 7.52, width = 13.33, units = 'in')
    # ggsave('Output/trend-figures/TrendsBlank.pdf', height = 7.52, width = 13.33, units = 'in')
    # ggsave('Output/trend-figures/TrendsIncreasing.pdf', height = 7.52, width = 13.33, units = 'in')
    
    ## Without the groups ####
    NumberTrendAll <- Trends %>%
      summarize(
        Increasing = sum(trend.prob.pos > 0.975),
        Decreasing = sum(trend.prob.pos <= 0.025),
        Stable = sum(DiffZero == 'Not different')
      )
    colors <- met.brewer('Benedictus', n = 8)
    # colors <- RColorBrewer::brewer.pal('RdYlBu', n = 8)
    # my.colors <- colors[c(8, 6, 3, 1)]
    # colScale <- scale_color_continuous(name = 'Trend', values = my.colors)
    
    # Changed to vertical based on Mollie's advice
    Lines <- ggplot(Trends, # ggplot(AllSpp,
           # aes(y = Species, x = trend.est.real
           aes(x = Species, y = trend.est.real
               #, color = Group
               , color = 1-trend.prob.pos
               )) +
      geom_point() +
      # geom_linerange(aes(xmin = trend.025.real, xmax = trend.975.real)) +
      geom_linerange(aes(ymin = trend.025.real, ymax = trend.975.real)) +
      # For a blank plot, uncomment these lines and comment the two above
      # geom_point(alpha = 0) +
      # geom_linerange(aes(min = trend.025, max = trend.975), alpha = 0) +
      theme_bw(base_size = 12) +  # Paper
      # theme_bw(base_size = 24) +  # Presentation
      theme(
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()#,
        # legend.position = 'bottom'
      ) +
      # labs(x = 'Annual percent change (95% CI)',
      labs(y = 'Annual percent change\n(95% CI)',
           color = 'Probability of decline') +
      geom_hline(yintercept = 0) +
      # geom_vline(xintercept = 0) +
      # ggpp::geom_text_npc(data = NumberTrends,
      #                     aes(npcx = 'left', npcy = 'top',
      #                         label = paste0('Increasing: ', Increasing)),
      #                     size = 5) +
      # ggpp::geom_text_npc(data = NumberTrends,
      #                     aes(npcx = 'left', npcy = 'top',
      #                         label = paste0('Increasing: ', Increasing, '\n',
      #                                        'Decreasing: ', Decreasing)),
      #                     size = 5) +
      # ggpp::geom_text_npc(data = NumberTrendAll,
      #                     aes(npcx = 'left', npcy = 'top',
      #                         label = paste0('Increasing: ', Increasing, '\n',
      #                                        'Decreasing: ', Decreasing, '\n',
      #                                        'No clear trend: ', Stable)),
      #                     size = 3) +
      # All species, red to blue
      scale_color_gradient2(midpoint = 0.5,
                            low = colors[8], mid = 'grey', high = colors[1],
                            limits = c(0,1)) +
      # To remove species overlapping 0, uncomment the lines below
      # scale_color_gradientn(breaks = c(0, 0.975, 1), 
      #                       colors = c('white', 'white', colors[1], colors[1]),
      #                       values = c(0, 0.9749, 0.975, 1)
      #                       ) +
      # xlim(yrange[1], yrange[2])
      ylim(yrange[1], yrange[2])
    # ggsave(filename = 'Output/trend-figures/TrendsNoGroups.jpg',
           # width = 180, height = 75, units = 'mm')
    # ggsave(filename = 'Output/trend-figures/TrendsNoGroupsColor.jpg',
    #        width = 180, height = 90, units = 'mm')
    # ggsave(filename = 'Output/trend-figures/TrendsNoGroupsColor2Blank.jpg',
    #        width = 13.3, height = 7.5, units = 'in')
    # ggsave(filename = 'Output/trend-figures/TrendsNoGroupsColor2Decreasing.jpg',
    #        width = 13.3, height = 7.5, units = 'in')
    # ggsave(filename = 'Output/trend-figures/TrendsNoGroupsColor2All.jpg',
    #        width = 13.3, height = 7.5, units = 'in')
    # ggsave(filename = 'Output/trend-figures/TrendsNoGroupsColor2.jpg',
    #        width = 180, height = 90, units = 'mm')
    
    
    Trends %>%
      dplyr::select(Species,
                    trend.est,
                    trend.est.real,
                    trend.prob.pos,
                    ProbDecreasing)
    
    ggplot(Trends, aes(x = ProbDecreasing)) +
      geom_histogram() +
      labs(y = 'Number of species', x = 'Probability of declining trend') +
      tbw
    # ggsave('Output/trend-figures/ProbDecline.pdf', width = 13.33, height = 7.52)
  } # SpeciesLinePlot
  
  ### Fanned figures of change over time ####
  if ( SpeciesFanPlot == TRUE ) {
    Trends %>%
      arrange(trend.est) %>%
      dplyr::select(Species, trend.est, ProbDecreasing)
    
    FannedTrends <- Trends %>%
      select(Species, ends_with('.real'), DiffZero, ProbDecreasing)
    
    FannedTrendsZero <- FannedTrends %>%
      mutate(StandardizedAbund = 1, Year = 1992)
    FannedTrends32 <- FannedTrends %>%
      mutate(StandardizedAbund = (1 + trend.est.real / 100) ^ 32,
             Year = 2023)
    FannedTrends <- FannedTrendsZero %>%
      bind_rows(FannedTrends32)
    FannedTrends %<>%
      mutate(
        Trend = case_when(
          ProbDecreasing <= 0.025 ~ 'Increasing',
          ProbDecreasing > 0.025 & ProbDecreasing < 0.975 &
            trend.est.real > 0 ~ 'No change; above zero',
          ProbDecreasing > 0.025 & ProbDecreasing < 0.975 &
            trend.est.real < 0 ~ 'No change; below zero',
          ProbDecreasing >= 0.975 ~ 'Decreasing'
        )
      )
    FannedTrends$Trend %<>%
      factor(levels = c(
        'Increasing',
        'No change; above zero',
        'No change; below zero',
        'Decreasing'
      ))
    
    colors <- met.brewer('Benedictus', n = 8)
    my.colors <- colors[c(8, 6, 3, 1)]
    names(my.colors) = c('Increasing',
                         'No change; above zero',
                         'No change; below zero',
                         'Decreasing')
    colScale <- scale_color_manual(name = 'Trend', values = my.colors)
    Fanned <- ggplot(FannedTrends,
           aes(
             x = Year,
             y = StandardizedAbund,
             # color = Trend,
             color = ProbDecreasing,
             group = Species
           )) +
      geom_line(alpha = 0.7) + 
      theme_bw(base_size = 12) +  # for small (90 mm) figure
      # tbw +
      scale_color_gradient2(midpoint = 0.5,
                            low = colors[8], mid = 'grey', high = colors[1],
                            limits = c(0,1)) +
      # colScale +
      # scale_color_grey() +  # If color is DiffZero
      theme(legend.position = 'bottom', # panel.grid = element_blank(),
            # panel.grid.major.y = element_blank(),
            panel.grid = element_blank()) +
      ylim(0, 1.6) +
      xlim(1992, 2023) +
      # guides(color = guide_legend(nrow = 2)) +
      labs(y = 'Standardized abundance',
           color = 'Probability of decline')
    # ggsave(filename = 'Output/trend-figures/StandardizedAbundance.jpg',
    #                     height = 90, width = 90, units = 'mm')
    # ggsave(filename = 'Output/trend-figures/StandardizedAbundance.pdf',
    #                     height = 6, width = 7, units = 'in')
  } # SpeciesFanPlot
  
ggarrange(Lines, Fanned, 
          labels = c('a', 'b'),
          # font.label = list(size = 12, color = 'black', face = 'bold', family = NULL),
          # vjust = 0.5,
          widths = c(2,1),
          common.legend = TRUE,
          legend = 'bottom') #+
  # theme(plot.margin = margin(0.2, 0, 0, 0, "cm"))
ggsave(filename = 'Output/trend-figures/SpeciesTrendCombinedVertical.jpg',
       width = 180, height = 90, units = 'mm')
ggsave(filename = 'Output/trend-figures/SpeciesTrendCombinedHorizontal2.jpg',
       width = 180, height = 90, units = 'mm')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Community Prep ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Total abundance per species ####
# Use to calculate the ten rarest and most common species
# Brackets to run the section
{
  # Calculate the area under the curve of the median values
  # Index of abundance to see how rare/common species are
  if (exists('PlotData') == FALSE) {
    load(paste0(
      path,
      'Output/12_AbundanceTrends/AllIndividualSpp-PlotData.RData'
    ))
  }
  
  SppAbundances <- PlotData %>%
    mutate(RealScaleMedian = exp(Median)) %>%
    group_by(Species) %>%
    summarize(AUC_Abundance = trapz(x = year, y = RealScaleMedian))
  
  SppAbundances %<>%
    arrange(AUC_Abundance) %>%
    mutate(AUC_Sum = cumsum(AUC_Abundance),
           SmallToLarge = 1:nrow(SppAbundances))
  
  # save(SppAbundances, file = 'Output/12_AbundanceTrends/AbundanceIndexCommonToRare.RData')
}

## Lists of species in each functional group ####
# Maybe could be done using the data frame inside the loop, but it's nice to
# to explicitly have the species in a list.
{
  # Bracket to run all of it more quickly
  SpeciesLists <- vector(mode = 'list', length = 13)
  
  SpeciesLists[[1]] <- Groups %>%
    filter(Group == 'Migratory') %>%
    use_series(Code)
  names(SpeciesLists)[1] <-  'Migratory'
  
  SpeciesLists[[2]] <- Groups %>%
    filter(MigrationVerdict == 'Resident') %>%
    use_series(Code)
  names(SpeciesLists)[2] <-  'Residents'
  
  SpeciesLists[[3]] <- Groups %>%
    filter(HostPlant == '>1') %>%
    use_series(Code)
  names(SpeciesLists)[3] <-  'Generalists'
  
  SpeciesLists[[4]] <- Groups %>%
    filter(HostPlant == '1') %>%
    use_series(Code)
  names(SpeciesLists)[4] <-  'Specialists'
  
  SpeciesLists[[5]] <- Groups %>%
    filter(Voltinism == 'U') %>%
    use_series(Code)
  names(SpeciesLists)[5] <-  'Univoltine'
  
  SpeciesLists[[6]] <- Groups %>%
    filter(Voltinism %in% c('B', 'M')) %>%
    use_series(Code)
  names(SpeciesLists)[6] <-  'Multivoltine'
  
  # Rare and common species based on above code chunks
  if (exists('Trends') == FALSE) {
    load(paste0(
      path,
      'Output/12_AbundanceTrends/AllIndividualSpp-Trends.RData'
    ))
    # Add the functional group data to the trends
    Trends %<>% left_join(Groups, join_by(Species == Code))
  }
  
  if (exists('SppAbundances') != TRUE) {
    load(file = 'Output/12_AbundanceTrends/AbundanceIndexCommonToRare.Rdata')
  }
  # Take the rarest species
  SpeciesLists[[7]] <- SppAbundances %>%
    slice_min(AUC_Abundance, n = 34) %>% use_series(Species)
  names(SpeciesLists)[7] <- 'Rare'
  
  # And the most common species
  SpeciesLists[[8]] <- SppAbundances %>%
    slice_max(AUC_Abundance, n = 34) %>% use_series(Species)
  names(SpeciesLists)[8] <- 'Common'
  
  # All species
  SpeciesLists[[9]] <- SppAbundances$Species
  names(SpeciesLists)[9] <- 'AllSpp'
  
  # Overwintering strategy
  SpeciesLists[[10]] <- Trends %>%
    filter(Overwinter == 'E') %>%
    use_series(Species)
  names(SpeciesLists)[10] <- 'OverwinterEgg'
  SpeciesLists[[11]] <- Trends %>%
    filter(Overwinter == 'L') %>%
    use_series(Species)
  names(SpeciesLists)[11] <- 'OverwinterLarva'
  SpeciesLists[[12]] <- Trends %>%
    filter(Overwinter == 'A') %>%
    use_series(Species)
  names(SpeciesLists)[12] <- 'OverwinterAdult'
  SpeciesLists[[13]] <- Trends %>%
    filter(Overwinter == 'P') %>%
    use_series(Species)
  names(SpeciesLists)[13] <- 'OverwinterPupa'
  SpeciesLists[[14]] <- Trends %>%
    filter(Overwinter == 'PL') %>%
    use_series(Species)
  names(SpeciesLists)[14] <- 'OverwinterPL'
}

# Save SpeciesLists for Richness component
# save(SpeciesLists, 
#      file = paste0(path, 'Output/12_AbundanceTrends/SpeciesLists.RData'))

### Trends by Abundance figure ####
FannedTrends %>% 
  left_join(SppAbundances) %>% 
  mutate(Quarter = case_when(
    SmallToLarge <= 34 ~ '0-25%',
    SmallToLarge > 34 & SmallToLarge <= 68 ~ '25-50%',
    SmallToLarge > 68 & SmallToLarge <= 102 ~ '50-75%',
    SmallToLarge > 102 ~ '75-100%'
  )) %>% 
  ggplot(., aes(x = Year, y = StandardizedAbund,
                           color = ProbDecreasing, group = Species)) +
  geom_line(alpha = 0.7) + 
  theme_bw(base_size = 12) +
  # colScale +
  scale_color_gradient2(midpoint = 0.5,
                        low = colors[8], mid = 'grey', high = colors[1],
                        limits = c(0,1)) +
  
  # scale_color_grey() +  # If color is DiffZero
  theme(
    legend.position = 'bottom',
    # panel.grid = element_blank(),
    # panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()) +
  ylim(0, 1.6) +
  xlim(1992, 2023) +
  # guides(color = guide_legend(nrow = 2)) +
  labs(y = 'Standardized abundance',
       color = 'Probability of decline') + 
  facet_grid(~ Quarter) +
  ggtitle('Grouped by abundance')

ggsave(filename = 'Output/trend-figures/FannedByAbundance.jpg',
       width = 180, height = 90, units = 'mm')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Community models ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set up
# Working with the saved BetaYear values
load(paste0(path, 'Output/12_AbundanceTrends/AllIndividualSpp-BetaYears.RData'))
BetaYears %>% length
BetaYears[[1]] %>% median
GroupResults <- vector('list', length = length(SpeciesLists))

## Each community ####
set.seed(511)
for(n.group in 1:length(SpeciesLists)) {
  ThisList <- BetaYears[
    (names(BetaYears) %in% SpeciesLists[[n.group]])
    ]
  y <- matrix(unlist(ThisList),
              nrow = 1000,
              ncol = length(ThisList))
  # Quick check
  y %>% as_tibble %>% colMeans()

  # Looks good
  covs <- data.frame(spp = 1:length(ThisList))
  new.data.list <- list(y = y, covs = covs)
  GroupResults[[n.group]] <- postHocLM(formula = ~ 1,
                     data = new.data.list, n.chains = 1, 
                     verbose = TRUE)
  names(GroupResults)[n.group] <- names(SpeciesLists)[n.group]
}

# Look at results
for ( n.group in 1:length(GroupResults)) {
  cat('\n************************************** \n\n', 
      names(GroupResults)[n.group], 
      '\n\n************************************** \n')
  GroupResults[[n.group]] %>% summary %>% print
}

# Save the all species model
AllSpp <- GroupResults[['AllSpp']]
Beta <- AllSpp$beta.samples[,1]
post.hoc.df.all.abund <- data.frame(
  trend.est = quantile(Beta, 0.5, names = FALSE),
  trend.prob.pos = mean(Beta > 0),
  trend.025 = quantile(Beta, 0.025, names = FALSE),
  trend.25 = quantile(Beta, 0.25, names = FALSE),
  trend.75 = quantile(Beta, 0.75, names = FALSE),
  trend.975 = quantile(Beta, 0.975, names = FALSE),
  Name = 'All',
  Metric = 'Abundance trend'
)

# save(post.hoc.df.all.abund,
#      file = paste0(path, 'Output/12_AbundanceTrends/AllAbundance.RData'))


## Comparison of communities ####
if ( exists('ButterflyComparePostHoc') == FALSE ) {
  stop('Grab ButterflyComparePostHoc from WendyFunctions')
}

ComparisonList <- list(c('Migratory', 'Residents') %>% sort,
                       c('Generalists', 'Specialists') %>% sort,
                       c('Univoltine', 'Multivoltine') %>% sort,
                       c('Rare', 'Common') %>% sort,
                       c("OverwinterEgg", "OverwinterLarva",
                         "OverwinterAdult", "OverwinterPupa",
                         "OverwinterPL") %>% sort)
OutputList <- vector('list', length = length(ComparisonList))
post.hoc.df <- data.frame(
  trend.est = numeric(),
  trend.prob.pos = numeric(),
  trend.025 = numeric(),
  trend.25 = numeric(),
  trend.75 = numeric(),
  trend.975 = numeric(),
  Name = character()
)
set.seed(288)
for ( ll.compare in 1:length(ComparisonList) ) {
  OutputList[[ll.compare]] <- ButterflyComparePostHoc(
    WhichLists = ComparisonList[[ll.compare]]
  )
  post.hoc.df %<>% 
    add_row(OutputList[[ll.compare]]$SummaryStats)
}

post.hoc.df %<>% 
  mutate(MeanCompare = ifelse(Name %>% str_detect('-'), 
                              'Compare', 'Mean'),
         trend.prob.neg = 1 - trend.prob.pos,
         Metric = 'Abundance')

NList <- SpeciesLists %>% lapply(n_distinct) %>% as_tibble %>% t
Ns <- tibble(Name = rownames(NList),
             n = NList[,1])

GroupTable <- post.hoc.df %>% 
  filter(MeanCompare == 'Mean') %>% 
  left_join(Ns) %>% 
  mutate(Name = factor(Name, levels = c("Migratory", "Residents", 
                                        "Generalists", "Specialists",
                                        "Multivoltine", "Univoltine",
                                        "Common", "Rare",
                                        "OverwinterEgg", "OverwinterLarva", 
                                        "OverwinterPL", "OverwinterPupa",
                                        "OverwinterAdult"))) %>% 
  mutate(percent.trend = trend.est * 100,
         percent.025 = trend.025 * 100,
         percent.975 = trend.975 * 100) %>% 
  select(Name, n, starts_with('percent'), trend.prob.neg) %>% 
  arrange(Name)

# write.csv(GroupTable, file = 'Output/12_AbundanceTrends/GroupMean.csv',
#             row.names = FALSE)

# Save file to disk
# save(post.hoc.df,
#      file = paste0(path, 'Output/12_AbundanceTrends/Communities.RData'))
     load(file = paste0(path, 'Output/12_AbundanceTrends/Communities.RData'))


post.hoc.df %>% 
  # filter(MeanCompare == 'Compare') %>% 
  select(Name, trend.est, trend.prob.pos, trend.prob.neg)
post.hoc.df %>% 
  filter(MeanCompare == 'Mean') %>% 
  select(Name, trend.est)
post.hoc.df %>% 
  filter(trend.prob.pos > 0.8 | 
           trend.prob.pos < 0.2) %>% 
  arrange(MeanCompare) %>% 
  select(Name, trend.est, trend.prob.pos)
post.hoc.df %>% 
  filter(trend.prob.pos > 0.975 | 
           trend.prob.pos < 0.025) %>% 
  arrange(MeanCompare) %>% 
  select(Name, trend.est, trend.prob.pos)
  
post.hoc.df %>% 
  filter(MeanCompare == 'Mean') %>% 
  select(Name, trend.est, trend.prob.neg) %>% 
  mutate(ProportionLeft = (1 + trend.est) ^ 32 %>% round(2))

# Scratch ####
Trends %>% use_series(trend.est.real) %>% summary
Trends %>% 
  filter(trend.est.real < -2) %>% 
  arrange(trend.est.real) %>% 
  select(Species, trend.est.real, trend.prob.pos)

dataNew %>% 
  filter(program == 'NFJ') %>% 
  # filter(yr < 1998) %>%
  # filter(yr > 2019) %>%
  group_by(GUEventID) %>% 
  filter(count != 0) %>% 
  summarize(NSpp = n_distinct(species),
            NCount = sum(count)) %>% 
  ungroup %>% 
  summarize(MeanNSpp = median(NSpp),
            MeanCount = median(NCount))
