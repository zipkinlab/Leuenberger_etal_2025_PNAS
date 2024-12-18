#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Update integrated community model for chapter 1
# Using spAbundance (Jeff's package)
# spOccupancy for updateMCMC
# Wendy Leuenberger
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# rm(list = ls())
# Make sure you pull the latest versions from CRAN for both of these guys
# if you haven't done so recently.
library(spAbundance)
library(spOccupancy)
library(tidyverse)
library(magrittr)
library(WendyFunctions)

# Load old model 
# ID <- '-Migratory21Spp-';                     n.batch = 250
# ID <- '-ResidentUnivoltineLycNymPieRio23Spp-';    n.batch = 300
# ID <- '-ResidentGeneralistLycNymPieRio19Spp-';  n.batch = 420
# ID <- '-ResidentSpecialistLycNymPieRio23Spp-';  n.batch = 200
ID <- '-ResidentUnivoltineHesPap11Spp-';  n.batch = 1300
# ID <- '-ResidentGeneralistHesPap15Spp-';  n.batch = 460
# ID <- '-ResidentSpecialistHesPap24Spp-';  n.batch = 300
Run <- 9.1 
NewRun <- Run + 0.1
IDRun <- paste0(ID, 'run', Run)
Path <- file.path('Output/1_runs') %>% paste0('/')
FilePath <- paste0(Path, IDRun, '.RData')

out <- extractorRData(file = FilePath,
                      object = IDRun)
data.list <- extractorRData(file = FilePath,
                            object = 'data.list')

# ThisGroup <- 'Migratory';                         n.batch = 250
# ThisGroup <- 'Resident Specialist HesPap';        n.batch = 300
# ThisGroup <- 'Resident Specialist LycNymPieRio';  n.batch = 200
# ThisGroup <- 'Resident Generalist';               n.batch = 100
# ThisGroup <- 'Resident Generalist HesPap';        n.batch = 460
# ThisGroup <- 'Resident Generalist LycNymPieRio';  n.batch = 420
# ThisGroup <- 'Resident Univoltine';               n.batch = 100
# ThisGroup <- 'Resident Univoltine HesPap';        n.batch = 1300
# ThisGroup <- 'Resident Univoltine LycNymPieRio';  n.batch = 300

paste('\n ************************************** \n \n Real data', '\n',
      'Functional group =', ID, '\n',
      'Model = msAbund', '\n',
      'Run =', Run, '\n \n',
      '**************************************
      ') %>% cat


# Update the previous model run -------------------------------------------
# NOTE: there is a stupid bug I just found in updateMCMC that happens when 
#       running models with random slopes/intercepts and only 1 chain 
#       which I'm pretty sure is what you're doing. So you will need the 
#       following line in order to get the  model to get updateMCMC to work properly.
p.occ <- dim(out$X)[length(dim(out$X))]

# Then this will work...
out.new <- updateMCMC(object = out, n.batch = n.batch, n.burn = 0, n.thin = 25,
                      keep.orig = FALSE, verbose = TRUE, n.report = 50, 
                      save.fitted = FALSE)

# plot(out.new, 'beta', density = FALSE)
# Notice the keep.orig argument. You can either set it to TRUE or FALSE depending
# on if you want to keep the original samples

ID2 <- paste0(ID, 'run', NewRun)

# Reprint snapshot at the end (prevent scrolling through everything)
paste('\n ************************************** \n \n Real data', '\n',
      'Functional group =', ID, '\n',
      # 'Species =', StanList$n_species, '\n',
      # 'Data length =', length(StanList$year_cov), 'data points \n',
      # 'Data size =', object.size(StanList) %>% format(units = 'MB', standard = 'SI'), '\n',
      'Batch =', n.batch, '\n', 
      'Run =', NewRun, '\n',
      'Model = msAbund', '\n',
      # 'Initial values =', UseMedians, '\n \n',
      '**************************************
      ') %>% cat

# Print output file name in slurm.out file
paste('\n ************************************** \n \n',
      'Output File Name:', ID2, '\n \n',
      '**************************************
      ') %>% cat

assign(ID2, out.new)
assign(ID, ID)

# Add data to the object as well
data.object = "data.list"
assign(data.object, data.list)
# dimnames(data.list$y)[[1]]

save(list = c(ID2, data.object, ID), file = file.path(getwd(), 'Output/1_runs', paste0(ID2, '.RData')))
