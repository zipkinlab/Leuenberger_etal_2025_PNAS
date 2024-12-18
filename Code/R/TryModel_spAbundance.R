#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run integrated community model for chapter 1
# Using spAbundance (Jeff's package)
# Wendy Leuenberger
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#------------------------------------------------------------#
# Set up: Load data and packages
#------------------------------------------------------------#

# Load packages
library(tidyverse)
library(magrittr)
library(spAbundance)

# Source data formatting code
# Lines 11 - 123 are for setting changes. Remainder can run based 
# on those settings
# Programs, species, years, months, ScientificName vs. Code, one species
# Preps all the way to standata
# Determine age of file
info <- file.info('Data/FormatData.RData')
code_info <- file.info('Code/R/FormatDataCh1.R')
# Calculate if file is <7 days old or not
LessThanOneWeek <- Sys.time() - info$mtime < lubridate::ddays(7)
# Calculate if code is older than data object 
CheckCode <- (code_info$mtime - info$mtime) < 0

# # Force it to be FALSE
# LessThanOneWeek <- FALSE
# CheckCode <- FALSE

# # Force it to be TRUE
# LessThanOneWeek <- TRUE
# CheckCode <- TRUE

# Load the R workspace if the RData file is <7 days old
# If the file is >7days old, run the source file and save the workspace
# If file is older than the code used to create it, rerun
if (LessThanOneWeek == TRUE & CheckCode == TRUE) {
  print('Loading data')
  load('Data/FormatData.RData')
  # Remove random seed that R automatically associates with a workspace
  # Causes results to not be random if just loading
  rm(.Random.seed)
} else {
  print('Sourcing data')
  source('Code/R/FormatDataCh1.R')
  save.image('Data/FormatData.RData')
}

# Use median values as initial values?
# UseMedians <- TRUE
UseMedians <- TRUE

# Process as if for Stan
# At some point, change into spAbundance-specific processing

# Set data subset
# Save original data for super quick retrieval
dataOriginal <- dataNew
dataNew <- dataOriginal

# NA's in one data point
dataNew %>% 
  filter(is.na(count))
# Check for NA's
CheckNA <- any(is.na(dataNew$count))
if(CheckNA == TRUE){
  NumNA <- is.na(dataNew$count) %>% sum
  dataNew$count[is.na(dataNew$count)] <- 0
  if(NumNA > 1){
    stop('More than the 1 known NA')
  } 
}

# Format data for spAbundance ---------------------------------------------
# StanList <- c(SurveyList,
#                SiteList,
#                OtherList)
# Check for NA's in Family
if(any(is.na(Groups$Family))){
  stop("NA's in Family")
}

# Subset
GroupType <- 'Family'
# GroupType <- 'Funct.Group'

# Note: n.batch * batch.size = total iter
# batch.size = 100
# Set here for easy access and ability to change code in only one place when 
# switching groups
# ; basically runs two separate lines of code.

# ThisGroup <- 'Migratory';                         n.batch = 250
# ThisGroup <- 'Resident Specialist HesPap';        n.batch = 300
# ThisGroup <- 'Resident Specialist LycNymPieRio';  n.batch = 200
# ThisGroup <- 'Resident Generalist';               n.batch = 100
# ThisGroup <- 'Resident Generalist HesPap';        n.batch = 460
# ThisGroup <- 'Resident Generalist LycNymPieRio';  n.batch = 420
# ThisGroup <- 'Resident Univoltine';               n.batch = 100
ThisGroup <- 'Resident Univoltine HesPap';        n.batch = 1300
# ThisGroup <- 'Resident Univoltine LycNymPieRio';  n.batch = 300
if(GroupType == 'Family'){
  print('Functional Group split by family')
  dataNew %<>% 
    filter(species %in% Groups$Code[Groups$GroupFamily == ThisGroup])
} else {
  print('Functional Group')
  dataNew %<>%
    filter(species %in% Groups$Code[Groups$Group == ThisGroup])
}

#------------------------------------------------------------#
# Check that all species have sufficient observations 
#------------------------------------------------------------#  

AtLeast50 <- dataNew %>% 
  group_by(species) %>% 
  summarize(Observations = sum(count, na.rm = TRUE)) %>%
  filter(Observations == min(Observations)) %>% 
  use_series(Observations)

print(AtLeast50)

if(AtLeast50 < 50){
  stop('At least one species has <50 observations')
}

# TOP TEN SPECIES ONLY
TopTen <- FALSE
if(TopTen == TRUE){
  TenMostAbundant <- dataNew %>% 
    group_by(species) %>% 
    summarize(Observations = sum(count)) %>% 
    # Uncomment for top ten
    # slice_max(order_by = Observations, n = 10)
    # Uncomment for random selection of species
    slice_sample(n = 14)
  dataNew <- dataNew %>% filter(species %in% TenMostAbundant$species)
}

#------------------------------------------------------------#
# Data processing that requires the final species list
#------------------------------------------------------------#  
# Processing that requires the final species list

#Add indicators for state monitoring programs (using NABA as a reference level)
dataNew$ia.ind <- ifelse(dataNew$program=='Iowa Butterfly Survey Network',1,0)
dataNew$il.ind <- ifelse(dataNew$program=='Illinois Butterfly Monitoring Network',1,0)
dataNew$mi.ind <- ifelse(dataNew$program=='Michigan Butterfly Network',1,0)
dataNew$oh.ind <- ifelse(dataNew$program=='Ohio Leps',1,0)

#Survey-level covariates in matrix
X_survey <- as.matrix(dataNew[,c('ia.ind','il.ind','mi.ind','oh.ind')])

SpeciesYear <- dataNew %>% 
  select(species, yr) %>% 
  distinct %>% 
  arrange(desc(species), yr) %>% 
  mutate(species_year = paste(species, yr, sep = '_'),
         species_year_cov = species_year %>% as_factor %>% as.numeric)
SpeciesCounty <- dataNew %>% 
  select(species, county.ind) %>% 
  distinct %>% 
  arrange(desc(species), county.ind) %>% 
  mutate(species_county = paste(species, county.ind, sep = '_'),
         species_county_cov = species_county %>% as_factor %>% as.numeric)

SpeciesSite <- dataNew %>% 
  select(species, site.ind) %>% 
  distinct %>% 
  arrange(desc(species), site.ind) %>% 
  mutate(species_site = paste(species, site.ind, sep = '_'),
         species_site_cov = species_site %>% as_factor %>% as.numeric)

dataNew %<>% 
  left_join(SpeciesYear) %>% 
  left_join(SpeciesCounty) %>% 
  left_join(SpeciesSite)

#------------------------------------------------------------#
# Combine ecological covariates into long form (nrows = n_counties * n_years * n_weeks + duplicate surveys in a week)
#------------------------------------------------------------# 
uyears <- sort(unique(dataNew$yr))
ucounties <- sort(unique(dataNew$county.ind))
n_counties <- length(ucounties)
uweeks <- sort(unique(dataNew$wk))
uspecies <- sort(unique(dataNew$species))

cyw <- expand.grid(wk=uweeks,county.ind=ucounties,yr=uyears,species=uspecies)
cyw$yr.ind <- cyw$yr - min(cyw$yr) + 1
#Standarize week
wk.m <- mean(cyw$wk)
wk.sd <- sd(cyw$wk)
cyw$wk.st <- (cyw$wk-wk.m)/wk.sd
# Create an index for unique combinations of county, year, and week
cyw$cyw.ind <- 1:nrow(cyw)

# Join to the observation dataframe
dataNew %<>% left_join(cyw) 

# Change county indicator to start at one
CountyInd <- tibble(county.ind = cyw$county.ind %>% unique,
                    co.ind = 1:length(unique(cyw$county.ind)))
cyw %<>% left_join(CountyInd)

# Species indicator
SpeciesInd = tibble(species = uspecies,
                    species_cov = 1:length(uspecies))
dataNew %<>% left_join(SpeciesInd)

# CountyKey <- dataNew %>% 
#   dplyr::select(county.ind, state.county.name) %>% 
#   distinct
# save(CountyKey, file = 'Data/CountyKey.RData')

# Format data into a list for stan
StanList <- list(year_cov = dataNew$yr.ind,
                 n_years = n_distinct(dataNew$yr.ind),
                 week_st = dataNew$wk.st,
                 y_count = dataNew$count,
                 n_surveys = n_distinct(dataNew$GUEventID),
                 surveys_s = paste0(dataNew$GUEventID, "_", dataNew$species),
                 n_surv_spp = nrow(dataNew),
                 n_cov_alpha = ncol(X_survey),
                 X_program = X_survey,
                 n_species = n_distinct(dataNew$species_cov),
                 species_cov = dataNew$species_cov,
                 n_species_years = n_distinct(dataNew$species_year_cov),
                 species_year_cov = dataNew$species_year_cov,
                 n_species_counties = n_distinct(dataNew$species_county_cov),
                 species_county_cov = dataNew$species_county_cov,
                 n_species_sites = n_distinct(dataNew$species_site_cov),
                 species_site_cov = dataNew$species_site_cov,
                 county_cov = dataNew$county.ind,
                 n_counties = n_distinct(dataNew$county.ind),
                 site_cov = dataNew$site_ind,
                 n_sites = n_distinct(dataNew$site_ind),
                 # crop_county_cov = dataNew$crop_county,
                 open_county_cov = dataNew$open_county,
                 open_site_cov = dataNew$open_site,
                 effort = dataNew$effort)


#------------------------------------------------------------#
# spAbundance-specific processing
#------------------------------------------------------------# 

# Species names
sp.names <- dataNew$species %>% unique %>% sort()
# Number of years
n.years <- StanList$n_years
# Number counties
n.counties <- StanList$n_counties
# Number of sites
n.sites <- StanList$n_sites
# Number of species
n.species <- StanList$n_species
# Number of weeks
n.weeks <- n_distinct(StanList$week_st)
# Number of additional surveys (more than one per week per site)
ExtraSurveys <- dataNew %>%
  group_by(yr, site_ind, species, wk) %>% 
  filter(duplicated(site_ind)) %>% 
  ungroup %>% 
  group_by(yr, site_ind, wk) 
ExtraSurveys %<>%
  group_by(yr, site_ind, wk, species) %>% 
  mutate(wkPlus = wk + (seq_along(site_ind) / 20) ) 
extra.surveys = ExtraSurveys$wkPlus %>% unique
n.extra.surveys = ExtraSurveys$wkPlus %>% n_distinct

dataNew %<>% 
  left_join(ExtraSurveys) %>% 
  mutate(wkPlus = ifelse(is.na(wkPlus), wk, wkPlus))

n.weeks.plus <- dataNew$wkPlus %>% n_distinct

# Total size of potential data points
n.total <- n.years * n.counties * n.species * (n.weeks + n.extra.surveys)
# Get linking indexes
sp.indx.obs <- StanList$species_cov
county.indx.obs <- StanList$county_cov
year.indx.obs <- StanList$year_cov
site.indx.obs <- StanList$site_cov

# TODO: Decide on approach to duplicate surveys in week
# Convert week to a factor then back to numeric to get it as a linking index
week.indx.obs <- as.numeric(factor(StanList$week_st))
week.plus.indx.obs <- as.numeric(factor(dataNew$wkPlus))
week.key <- cbind(week.plus.indx.obs, wkPlus = dataNew$wkPlus) %>% 
  unique

# Total number of observations
n.obs.total <- length(StanList$y_count)
# Get all observational level data in a data frame
obs.data.df <- data.frame(years = year.indx.obs,
                          sites = site.indx.obs,
                          surveys = StanList$surveys_s,
                          counties = StanList$county_cov,
                          species = sp.indx.obs,
                          wk = week.indx.obs,
                          wkplus = week.plus.indx.obs,
                          counts = StanList$y_count,
                          StanList$X_program,
                          lat = dataNew$lat,
                          long = dataNew$long,
                          program = dataNew$program,
                          # crop = dataNew$crop_county,
                          openco = dataNew$open_county,
                          opensi = dataNew$open_site,
                          effort = dataNew$effort)
# Total number of butterflies observed at a given site, in a given year, during
# a given week. Some of the state-wide surveys have multiple observations per week,
# which is accounted for by using the number of surveys in a given week as a variable
# in the model.
count.df <- obs.data.df %>%
  # group_by(years, sites, counties, species, wk) %>%
  group_by(years, sites, counties, species, wk, wkplus) %>%
  summarize(sum.counts = sum(counts),
            n.surveys = n_distinct(surveys),
            ia.ind = unique(ia.ind),
            il.ind = unique(il.ind),
            mi.ind = unique(mi.ind),
            oh.ind = unique(oh.ind),
            # crop.cov = unique(crop),
            openco.cov = unique(openco),
            opensi.cov = unique(opensi),
            effort.cov = unique(effort)) %>%
  ungroup()

# Create arrays of data for use in spAbundance
y <- array(NA, dim = c(n.species, n.sites, n.years, n.weeks.plus))
p.abund <- StanList$n_cov_alpha
program.covs <- array(NA, dim = c(n.sites, n.years, n.weeks.plus, p.abund))
year.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
week.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
week.plus.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
site.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
county.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
survey.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
# crop.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
openco.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
opensi.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
effort.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
nsurveys.cov <- array(NA, dim = c(n.sites, n.years, n.weeks.plus))
for (i in 1:nrow(count.df)) {
  # Put the number of butterflies of a given species for a given site/year/week
  # into its spot in the y array
  # y is species x site x yr x wkplus in shape
  y[count.df$species[i], count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$sum.counts[i]
  if (count.df$species[i] == 1) {
    # Put the values for the site and survey level covariates into their respective arrays
    # These are site x yr x wkplus in shape
    # Only needs to be done once per site/yr/week, so only do it for species 1
    program.covs[count.df$sites[i], count.df$years[i], count.df$wkplus[i], 1] <- count.df$ia.ind[i]
    program.covs[count.df$sites[i], count.df$years[i], count.df$wkplus[i], 2] <- count.df$il.ind[i]
    program.covs[count.df$sites[i], count.df$years[i], count.df$wkplus[i], 3] <- count.df$mi.ind[i]
    program.covs[count.df$sites[i], count.df$years[i], count.df$wkplus[i], 4] <- count.df$oh.ind[i]
    year.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$years[i]
    week.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$wk[i]
    week.plus.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$wkplus[i]
    site.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$sites[i]
    county.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$counties[i]
    survey.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$n.surveys[i]
    # crop.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$crop.cov[i]
    openco.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$openco.cov[i]
    opensi.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$opensi.cov[i]
    effort.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$effort.cov[i]
    nsurveys.cov[count.df$sites[i], count.df$years[i], count.df$wkplus[i]] <- count.df$n.surveys[i]
  }
}

# Investigate
# # All species at site 1, year 1, week 1
# y[ ,1,1,1]
# # Species 1 at all sites, year 1, week 1
# y[1, ,1,1]
# # Species 1 at site 1, all years, week 1
# y[1,1, ,1]
# # Species 1 at site 1, year 1, all weeks
# y[1,1,1, ]
# Lots of NA's - surveys that didn't occur in a given site/year/wk
# 0 represents surveys that occurred where a species was not observed
# Spot check with some actual surveys. Looks good
# # Check without week.plus
# y[1, 21:27, 10:13, 7]
# dataNew %>%
#   filter(species == 'ATACAM',
#          wk == 7,
#          yr.ind %in% 10:13,
#          site_ind %in% 21:27) %>%
#   select(usiteID, yr, wk, count, site_ind)
# Check with week.plus
CheckWkPlus <- dataNew %>% 
  filter(usiteID == 'IA-27001',
         species_cov == 1,
         wk == 15) %>% 
  select(usiteID, site_ind, yr, yr.ind, wk, wkPlus, species, species_cov, count) %>% 
  left_join(
    (week.key %>% as_tibble %>% filter(wkPlus >= 15 & wkPlus < 16))
  )
CheckWkPlus
# Pull out based on indicators
y[1, 425, 19, 22:26]
# Looks good!

# Collapse the data into a species x site x "replicate" three-dimensional array.
# Order of the replicates is week, then year within week.
y <- array(y, dim = c(n.species, n.sites, n.years * n.weeks.plus))
ia.ind <- array(program.covs[, , , 1], dim = c(n.sites, n.years * n.weeks.plus))
il.ind <- array(program.covs[, , , 2], dim = c(n.sites, n.years * n.weeks.plus))
mi.ind <- array(program.covs[, , , 3], dim = c(n.sites, n.years * n.weeks.plus))
oh.ind <- array(program.covs[, , , 4], dim = c(n.sites, n.years * n.weeks.plus))
year.cov <- array(year.cov, dim = c(n.sites, n.years * n.weeks.plus))
week.cov <- array(week.cov, dim = c(n.sites, n.years * n.weeks.plus))
week.plus.cov <- array(week.plus.cov, dim = c(n.sites, n.years * n.weeks.plus))
site.cov <- array(site.cov, dim = c(n.sites, n.years * n.weeks.plus))
county.cov <- array(county.cov, dim = c(n.sites, n.years * n.weeks.plus))
survey.cov <- array(survey.cov, dim = c(n.sites, n.years * n.weeks.plus))
# crop.cov <- array(crop.cov, dim = c(n.sites, n.years * n.weeks.plus))
openco.cov <- array(openco.cov, dim = c(n.sites, n.years * n.weeks.plus))
opensi.cov <- array(opensi.cov, dim = c(n.sites, n.years * n.weeks.plus))
effort.cov <- array(effort.cov, dim = c(n.sites, n.years * n.weeks.plus))
nsurveys.cov <- array(nsurveys.cov, dim = c(n.sites, n.years * n.weeks.plus))

# for(ii in 1:max(count.df$sites)){
#   Check <- y[,ii,] %>% all(is.na)
#   if(Check == TRUE){
#     stop('All NA')
#   }
# }

# Remove missing values to make size of stored objects smaller. In spAbundance,
# this works for non-spatial models to reduce the size of model objects.
y <- matrix(y, n.species, n.sites * n.years * n.weeks.plus)
# Each species has a slice, but only need to look at species 1 because all will
# have the same NA's
# To test this:
identical(is.na(y[1,]), is.na(y[2,]))
y.na <- is.na(y[1, ])
y <- y[, !y.na]
# Only one "slice" aka a dataframe for program data. Can remove the exact same
# positions because they're surveys that didn't happen
ia.ind <- c(ia.ind)[!y.na]
il.ind <- c(il.ind)[!y.na]
mi.ind <- c(mi.ind)[!y.na]
oh.ind <- c(oh.ind)[!y.na]
year.cov <- c(year.cov)[!y.na]
week.cov <- c(week.cov)[!y.na]
week.plus.cov <- c(week.plus.cov)[!y.na]
site.cov <- c(site.cov)[!y.na]
county.cov <- c(county.cov)[!y.na]
# crop.cov <- c(crop.cov)[!y.na]
openco.cov <- c(openco.cov)[!y.na]
opensi.cov <- c(opensi.cov)[!y.na]
effort.cov <- c(effort.cov)[!y.na]
nsurveys.cov <- c(nsurveys.cov)[!y.na]
# Subtract one from the survey variable such that the intercept is equal to
# abundance for a single survey.
# Currently starts at 1
c(survey.cov) %>% table
survey.cov <- c(survey.cov)[!y.na] - 1

# Get data in format for msAbund
abund.covs <- list(ia.ind = ia.ind,
                   il.ind = il.ind,
                   mi.ind = mi.ind,
                   oh.ind = oh.ind,
                   year.cov = year.cov,
                   week.cov = week.cov,
                   week.plus.cov = week.plus.cov,
                   site.cov = site.cov,
                   county.cov = county.cov,
                   survey.cov = survey.cov,
                   # crop.cov = crop.cov,
                   openco.cov = openco.cov,
                   opensi.cov = opensi.cov,
                   effort.cov = effort.cov,
                   nsurveys.cov = nsurveys.cov)
# Add species names
rownames(y) <- sp.names
data.list <- list(y = y,
                  covs = abund.covs)

# Checks
data.list$covs$week.cov %>% table
data.list$covs$week.plus.cov %>% table



# save(data.list, file = 'Data/spAbundance-data.rda')

# Make ID to save objects and load initial values -------------------------

ID <- paste0('-', ThisGroup %>% str_replace_all(" ", ""),
             StanList$n_species, "Spp-")

# Prepare initial values based on medians from previous run ---------------
if(UseMedians == TRUE){
  load(
    paste0('Output/3_medians/', ID,
                'medians',
                '.RData')
  )
  # Rerun the model with new initial values. Note I add a very small amount of noise
  beta.comm.inits <- beta.comm.meds + runif(length(beta.comm.meds), -0.1, 0.1)
  
  tau.sq.beta.inits <- tau.sq.beta.meds + runif(length(tau.sq.beta.meds), 0.1, 0.1)
  # Ensure all tau.sq.beta values are positive
  tau.sq.beta.inits <- ifelse(tau.sq.beta.inits < 0, -tau.sq.beta.inits, tau.sq.beta.inits)
  
  beta.inits <- beta.meds + runif(length(beta.meds), -0.1, 0.1)
  
  kappa.inits <- kappa.meds + runif(length(kappa.meds), -0.1, 0.1)
  # Ensure all kappa values are positive
  kappa.inits <- ifelse(kappa.inits < 0, -kappa.inits, kappa.inits)
  
  beta.star.inits <- beta.star.meds + runif(length(beta.star.meds), -0.1, 0.1)
  
  sigma.sq.mu.inits <- sigma.sq.mu.meds + runif(length(sigma.sq.mu.meds), 0.1, 0.1)
  # Ensure all sigma.sq.mu values are positive
  sigma.sq.mu.inits <- ifelse(sigma.sq.mu.inits < 0, -sigma.sq.mu.inits, sigma.sq.mu.inits)
  
  inits <- list(beta.comm = beta.comm.inits,
                tau.sq.beta = tau.sq.beta.inits,
                # Note beta needs to be a matrix
                beta = matrix(beta.inits, nrow = StanList$n_species), 
                kappa = kappa.inits,
                sigma.sq.mu = sigma.sq.mu.inits,
                beta.star = beta.star.inits)
}


# Run the model in spAbundance --------------------------------------------
# Settings for 47000 should run in 156 hours, about the max for a week

# n.batch moved! 
# Moved to top where the functional group is set because it changes with group
n.batch # take a look at the value
# n.batch = 1880; batch.length = 25; n.burn = 23000 
# n.batch <- 1300  # 1200 1500All running on same core

batch.length <- 100  # 25
n.burn <- n.batch * batch.length * 0.5  # 15000 20000
n.thin <- 25
n.chains <- 1  # Leave here and submit three times to HPCC to get 3 parallel chains
# Jeff will send code to posthoc combine them

# n.batch <- 50; batch.length <- 2; n.burn = 1; n.thin = 1; n.chains = 1

paste('\n ************************************** \n \n Real data', '\n',
      'Functional group =', ThisGroup, '\n',
      'Species =', StanList$n_species, '\n',
      'Data length =', length(StanList$year_cov), 'data points \n',
      'Data size =', object.size(StanList) %>% format(units = 'MB', standard = 'SI'), '\n',
      'Total iter =', n.batch * batch.length, '\n', 
      'Burn-in =', n.burn, '\n',
      'Model = msAbund', '\n',
      'Initial values =', UseMedians, '\n \n',
      '**************************************
      ') %>% cat

tuning <- list(beta = 0.1, kappa = 0.1, beta.star = 0.1)

if(UseMedians == TRUE){
  print('Using medians as initial values')
  out <- msAbund(formula = ~ (1 | county.cov) + # Will estimate one per spp
                   # Site.cov takes a long time. Uncomment if needed though
                   # (1 | site.cov) + # Will estimate one per spp
                   # (1 | species.cov) + # Intercept is ranef of species
                   (1 | year.cov) + # Will estimate one per spp
                   # program (each as individual covs)
                   ia.ind + il.ind + mi.ind + oh.ind + 
                   # crop.cov + 
                   openco.cov + opensi.cov + 
                   # number of surveys per week
                   # Done: change indexing to remove this term
                   # scale(nsurveys.cov) +
                   scale(week.cov) + 
                   # (x|group) random slope of x within group with correlated intercept
                   (scale(week.cov) | year.cov) + 
                   I(scale(week.cov)^2) + 
                   (I(scale(week.cov)^2) | year.cov) +
                   log(effort.cov),
                 data = data.list,
                 # Initial values based on previous medians
                 inits = inits,
                 tuning = tuning,
                 # Overdispersion will be species specific but not linked to a community distribution
                 family = 'NB',
                 # Dramatically decrease size of final object
                 save.fitted = FALSE,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.report = 5,
                 n.chains = n.chains)
} else {
  print('No initial values')
  out <- msAbund(formula = ~ (1 | county.cov) + # Will estimate one per spp
                   # Site.cov takes a long time. Uncomment if needed though
                   # (1 | site.cov) + # Will estimate one per spp
                   # (1 | species.cov) + # Intercept is ranef of species
                   (1 | year.cov) + # Will estimate one per spp
                   # program (each as individual covs)
                   ia.ind + il.ind + mi.ind + oh.ind + 
                   # crop.cov + 
                   openco.cov + opensi.cov + 
                   # number of surveys per week
                   # Done: change indexing to remove this term
                   # scale(nsurveys.cov) +
                   scale(week.cov) + 
                   # (x|group) random slope of x within group with correlated intercept
                   (scale(week.cov) | year.cov) + 
                   I(scale(week.cov)^2) + 
                   (I(scale(week.cov)^2) | year.cov) +
                   log(effort.cov),
                 data = data.list,
                 tuning = tuning,
                 # Overdispersion will be species specific but not linked to a community distribution
                 family = 'NB',
                 # Dramatically decrease size of final object
                 save.fitted = FALSE,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.report = 5,
                 n.chains = n.chains)
}


ID2 <- paste0(ID, 'run',
              length(list.files(path = file.path(getwd(), 'Output/1_runs'),
                                pattern = ID,
                                full.names = FALSE)) + 1)

# Reprint snapshot at the end (prevent scrolling through everything)
paste('\n ************************************** \n \n Real data', '\n',
      'Functional group =', ThisGroup, '\n',
      'Species =', StanList$n_species, '\n',
      'Data length =', length(StanList$year_cov), 'data points \n',
      'Data size =', object.size(StanList) %>% format(units = 'MB', standard = 'SI'), '\n',
      'Batch =', n.batch, '\n', 
      'Burn-in =', n.burn, '\n',
      'Model = msAbund', '\n',
      'Initial values =', UseMedians, '\n \n',
      '**************************************
      ') %>% cat

# Print output file name in slurm.out file
paste('\n ************************************** \n \n',
      'Output File Name:', ID2, '\n \n',
      '**************************************
      ') %>% cat

assign(ID2, out)
assign(ID, ID)

# Add data to the object as well
data.object = "data.list"
assign(data.object, data.list)
# dimnames(data.list$y)[[1]]

save(list = c(ID2, data.object, ID), file = file.path(getwd(), 'Output/1_runs', paste0(ID2, '.RData')))
