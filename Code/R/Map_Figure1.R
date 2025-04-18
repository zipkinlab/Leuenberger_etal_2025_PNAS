# Code to make map of number of surveys at each survey location

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

# Load the R workspace if the RData file is <7 days old
# If the file is >7days old, run the source file and save the workspace
# If file is older than the code used to create it, rerun
if (LessThanOneWeek == TRUE & CheckCode == TRUE) {
  load('Data/FormatData.RData')
  # Load the packages loaded in FormatDataMethod2.R
  library(magrittr)
  library(tidyverse)
  library(ggplot2)
  library(tidybayes)
  library(sf)
  library(sp)
  library(maps)
  # library(maptools)  # Will be deprecated... need map2SpatialPolygons or alternative
  # Different color palettes
  library(RColorBrewer)
  # library(wesanderson)
  library(ggsci)
} else {
  source('Code/R/FormatDataCh1.R')
  save.image('Data/FormatData.RData')
}

# Load information about county names
load(file = 'Data/CountyKey.RData')

# Create map
# Find the number of surveys and a place and format data
counts.by.lat.long <- dataNew %>%
  # filter(Code %in% CaseStudySpp$Code) %>%
  group_by(program, usiteID, lat, long) %>%
  # summarize(Count = sum(count)) %>%
  summarize(Surveys = n_distinct(GUEventID)) %>% 
  ungroup() %>%
  mutate(program = ifelse(program == 'NFJ',
                          'NABA',
                          ifelse(program == 'Ohio Leps',
                                 'OH',
                                 ifelse(program == 'Illinois Butterfly Monitoring Network',
                                        'IL', ifelse(program == 'Iowa Butterfly Survey Network',
                                                     'IA', ifelse(program == 'Michigan Butterfly Network', 'MI', NA)))))) %>%
  mutate(program = factor(program, levels = c('IL', 'IA', 'MI', 'OH', 'NABA'), ordered = TRUE))
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
coords.sf <- st_as_sf(counts.by.lat.long,
                      coords = c('long', 'lat'),
                      crs = 4326)
coords.sf <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Arrange coords.sf with NABA last in the dataset
coords.sf %<>% 
  arrange(program) 
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Full data
mw.states <- usa %>%
  filter(ID %in% c('iowa', 'wisconsin', 'michigan', 'indiana',
                   'illinois', 'ohio', 'missouri', 'minnesota'))
my.cols <- c('#440154FF', pal_tron()(5))
my.cols <- my.cols[1:5]
names(my.cols) <- c('NABA', 'IA', 'IL', 'MI', 'OH')

# 2-column width plot   
ggplot() +
  geom_sf(data = mw.states, col = 'black', fill = 'white', alpha = 1) +
  # geom_sf(data = coords.sf, aes(col = program, size = Count)) +
  geom_sf(data = coords.sf, aes(col = program, size = Surveys)) +
  theme_bw(base_size = 18) +
  scale_color_manual(values = my.cols) +
  guides(color = guide_legend(override.aes = list(size = 2),
                              order = 1)) +
  # scale_size(breaks = c(1, 1000, 5000, 10000), labels = c('1', '1,000', '5,000', '10,000')) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(color = 'Program')

# Save plot to hard drive.
ggsave(filename = 'Output/trend-figures/PNAS_Figure1.eps',
       device = 'eps',
       width = 178, height = 135, units = 'mm')