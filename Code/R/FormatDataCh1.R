#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Butterfly data prep
# Chapter one
# Updated data format (everything on one sheet)
# Integrated with updated NABA data
# Wendy Leuenberger
# Date created: 5/12/2023
# Date updated: 01/09/2024
# Updated by: Michael Belitz
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load packages ####
library(magrittr)
library(tidyverse)
library(waldo)
library(readxl)
# For county/state classification
library(sf)
library(maps)
library(tigris)
# Different color palettes
library(RColorBrewer)
# library(wesanderson)
library(ggsci)
# Processing dates/times
library(lubridate)

# Functions ####

#' Determine counties for each coordinate ####

#' MWB's attempt to determine county for each coordinate 
#' my attempt isn't as efficient in computational speed as the old function, 
#' but it does a better job matching points to counties
#' input parameter is a dataframe with the columns GULongitude and GULatitude

latlon2counties <- function(pointsDF){
  
  state.fips <- maps::state.fips %>%  # data frame with state fips, abreviations, and names
    dplyr::select(fips, abb, polyname)
  
  co <- tigris::counties(resolution = "5m") %>%  # simple feature of counties in the United States
    dplyr::select(STATEFP, COUNTYFP, NAME, geometry) %>%
    mutate(STATEFP = as.integer(STATEFP)) %>% 
    left_join(state.fips, by = c("STATEFP" = "fips"), relationship = "many-to-many") %>% 
    mutate(ids = paste0(str_to_lower(word(polyname, start = 1, end = 1, sep = ":")), 
                        ",", str_to_lower(NAME))) %>% 
    # filter(abb %in% States) %>% 
    st_as_sf() %>% 
    st_transform(crs = "+proj=longlat +datum=WGS84")
  
  pdf <- pointsDF %>% # make point dataframe into simple feature
    st_as_sf(coords = c("GULongitude", "GULatitude"),
             crs = "+proj=longlat +datum=WGS84")
  
  pdf <- pdf %>% # create unique identifier for each record, since one-to-many joins are possible with complex shapefiles
    mutate(uid = 1:nrow(pdf))
  
  coords_df <- data.frame( # make coorinate dataframe since simple features will remove coordinate columns
    uid = pdf$uid,
    GULongitude = st_coordinates(pdf)[,1],
    GULatitude = st_coordinates(pdf)[,2]
  )
  
  point_counties <- pdf %>%  # spatial join point data frame with county dataframe
    st_join(co, join = st_within) %>% 
    st_drop_geometry() %>% 
    distinct(uid, .keep_all = TRUE)
  
  point_counties <- left_join(point_counties, coords_df) # add coordinates
  
  return(point_counties)
  
}

# Fill in the state if NA ####
# Some States are NA. Use this function to fill in the state if
# the county is one of the ones recorded in the state
# Won't work if the survey with state = NA is the only one in that
# county or if the survey is in a different state
FillState <- function(Data, StateCode){
  StateCounties <- Data %>% 
    # Grab just rows recorded as in the particular state
    filter(State == StateCode) %>%
    # Grab all the counties in that state with recorded surveys
    use_series(state.county) %>% 
    # Take unique values so that it is simpler
    unique
  # Find the unique values and assign them to that state if the 
  # county has other recorded surveys in that county/state
  Data$State[is.na(Data$State) &
               Data$state.county %in% StateCounties] <- StateCode
  return(Data)
}

# Standard error function
se <- function(x, na.rm = FALSE){ 
  sqrt(var(x, na.rm = na.rm) / length(x))
}



# Load data ####

# butterfly names
Names <- read_csv('Data/BflyNames.csv') %>% 
  mutate(Taxon = paste(Genus, Species)) %>% 
  distinct(Taxon, .keep_all = TRUE) %>% 
  # Adding family to split larger functional groups
  select(Taxon, UMD_Code, 'Common Name', Family)


# Start with NABA Data
naba <- read_csv('Data/NABAAllData_1977toOct2022.csv') %>% 
  mutate(Data = "nabaUpdate")

naba <- naba %>% 
  select(ID, GU_Loc_NUM, CountNameGU, ObsYear, ObsMonth, ObsDay, State, Country,
         Latit_GU, Long_GU, UMD_CODE, `Scientific  Name`, Num_Butts_GU,
         Party_Hours_GU) %>% 
  mutate(Program = "NFJ",
         UniqueSurveyID = NA,
         NABACode = NA,
         Temp = NA,
         Wind = NA,
         GUEventID = paste(Program, GU_Loc_NUM, ObsMonth, ObsDay, ObsYear, sep = "-"),
         GULocationID = paste(Program, GU_Loc_NUM, sep = "-"),
         EventType = "Count") %>% 
  rename(GULatitude = Latit_GU,
         GULongitude = Long_GU,
         Day = ObsDay,
         Month = ObsMonth,
         Year = ObsYear,
         ScientificName = `Scientific  Name`,
         Code = UMD_CODE,
         Count = Num_Butts_GU,
         Duration = Party_Hours_GU) %>% 
  select(ID, UniqueSurveyID, Program, GUEventID, GULocationID, EventType,
         GULatitude, GULongitude, Country, State, Day, Month, Year, 
         ScientificName, Code, NABACode, Count, Duration, Temp, Wind)

# now read in BMN data
#read in excel sheets
IA_IL_MI_sheets <- excel_sheets("Data/IA IL MI Pollardbase data through 10.04.2023.xlsx")
OH_sheets <- excel_sheets("Data/GU Ohio 2022 Flat File All Observations.xlsx")
coords <- excel_sheets("Data/Pb-route-latlon-IA-IL-MI.xlsx")

# pull out survey info
IA_IL_MI_sheets

IA_surveys <- read_excel(path = "Data/IA IL MI Pollardbase data through 10.04.2023.xlsx",
                         sheet = IA_IL_MI_sheets[2]) %>% 
  mutate(Program = "Iowa Butterfly Survey Network",
         State = "IA")
IL_surveys <- read_excel(path = "Data/IA IL MI Pollardbase data through 10.04.2023.xlsx",
                         sheet = IA_IL_MI_sheets[4]) %>% 
  mutate(Program = "Illinois Butterfly Monitoring Network",
         State = "IL")
MI_surveys <- read_excel(path = "Data/IA IL MI Pollardbase data through 10.04.2023.xlsx",
                         sheet = IA_IL_MI_sheets[6]) %>% 
  mutate(Program = "Michigan Butterfly Network",
         State = "MI")

IA_IL_MI_surveys <- bind_rows(IA_surveys, IL_surveys, MI_surveys) %>% 
  mutate(Break_minutes = if_else(Break_minutes == 9999, true = 0, false = Break_minutes)) %>% 
  mutate(Duration_hours = (Duration_minutes - Break_minutes) / 60,
         eventDate = as.Date(Date_Start_Time),
         year = year(Date_Start_Time),
         Temp = `Temp _start`,
         Wind = Wind_start) %>% 
  rename(Survey_ID = `Survey ID`, Site_Name = `Site Name`, 
         Number_Taxa = `Number taxa`, Total_Individuals = `Total individuals`) %>% 
  select(Program, State, Survey_ID, eventDate, year, Date_Start_Time, Duration_hours, 
         Site_Name, Route, Number_Taxa,Total_Individuals, Temp, Wind) %>% 
  mutate(Route = if_else(Site_Name == "Richfield Park" & !is.na(Site_Name), "Richfield County Park Route 1", Route),# manually change names of 2 sites/routes to match sites sheet 
         Site_Name = if_else(Route == "Parfet Prairie Route 1", "Parfet Prairie", Site_Name),# manually change names of 2 sites/routes to match sites sheet 
         Duration_hours = if_else(Duration_hours < 0, true = 0, false = Duration_hours)) # manually change duration for durations that don't make sense


# read in site information
IA_IL_MI_Sites <- read_excel(path = "Data/Pb-route-latlon-IA-IL-MI.xlsx",
                             sheet = coords[2]) %>% 
  select(ID, `Route ID`, Site, `Route Name`, County, Latitude, Longitude, `Location source`) %>% 
  rename(Route_Name = `Route Name`,
         Route_ID = `Route ID`) %>% 
  distinct(Site, Route_Name, .keep_all = TRUE) # remove three site by route combinations with double coordinates

# join IA_IL_MI survey data with site information
IA_IL_MI_surveys <- left_join(IA_IL_MI_surveys, IA_IL_MI_Sites, by = c("Site_Name" = "Site",
                                                                       "Route" = "Route_Name")) 

# OH_surveys # Ohio data has survey with observation data
OH_obs <- read_excel(path = "Data/GU Ohio 2022 Flat File All Observations.xlsx",
                     sheet = OH_sheets[2])

IA_obs<- read_excel(path = "Data/IA IL MI Pollardbase data through 10.04.2023.xlsx",
                    sheet = IA_IL_MI_sheets[3]) %>% 
  mutate(Program = "Iowa Butterfly Survey Network") 
IA_obs$`Survey ID` <- as.character(IA_obs$`Survey ID`)

IL_obs <- read_excel(path = "Data/IA IL MI Pollardbase data through 10.04.2023.xlsx",
                     sheet = IA_IL_MI_sheets[5])%>% 
  mutate(Program = "Illinois Butterfly Monitoring Network")

MI_obs <- read_excel(path = "Data/IA IL MI Pollardbase data through 10.04.2023.xlsx",
                     sheet = IA_IL_MI_sheets[7])  %>% 
  mutate(Program = "Michigan Butterfly Network")

IA_IL_MI_obs <- bind_rows(IA_obs, IL_obs, MI_obs) %>% 
  rename(Survey_ID = `Survey ID`, Site_Name = `Site Name`) %>% 
  mutate(Route = if_else(Site_Name == "Richfield Park" & !is.na(Site_Name), "Richfield County Park Route 1", Route),# manually change names of 2 sites/routes to match sites sheet 
         Site_Name = if_else(Route == "Parfet Prairie Route 1", "Parfet Prairie", Site_Name)) # manually change names of 2 sites/routes to match sites sheet

# combine IA, IL, and MI with obs
IA_IL_MI_data <- left_join(IA_IL_MI_obs, IA_IL_MI_surveys,
                           by = c("Survey_ID", "Site_Name", "Program", "Route"))
IA_IL_MI_zeros <- filter(IA_IL_MI_surveys, Total_Individuals == 0)

OH_data <- OH_obs %>% 
  mutate(Program = "Ohio Leps") %>% 
  filter(UMDCode != "NONE")
OH_zeros <- filter(OH_obs, UMDCode == "NONE") %>% 
  mutate(Program = "Ohio Leps")

rm(IA_IL_MI_obs, IA_IL_MI_surveys, IA_obs, MI_obs, IL_obs, OH_obs,
   IA_surveys, MI_surveys, IL_surveys, IA_IL_MI_sheets, OH_sheets)
gc()


# what data doesn't have coordinates 
noCoords <- filter(IA_IL_MI_data, is.na(Latitude) | is.na(Longitude)) %>% 
  distinct(Site_Name, State, County) # everything has coordinates :)

# Add an error in case we add data and don't have the coordinates
if(dim(noCoords)[1] != 0){
  stop("Some data doesn't have coordinates")
}

rm(noCoords)
gc()

# format data for cleaning script
IA_IL_MI_data <-  IA_IL_MI_data %>% 
  left_join(Names, by = "Taxon") %>% 
  mutate(GULocationID = paste(State, ID.y, sep = "-")) %>% 
  mutate(GUEventID = paste(GULocationID, Survey_ID, sep = "-"),
         EventType = "Pollard",
         Country = "USA",
         Day = lubridate::mday(eventDate),
         Month = lubridate::month(eventDate),
         NABACode = NA) %>% 
  rename(UniqueSurveyID = Survey_ID,
         ScientificName = Taxon,
         Code = UMD_Code,
         Year = year,
         Duration = Duration_hours,
         ID = ID.y,
         GULatitude = Latitude,
         GULongitude = Longitude) %>% 
  select(ID, UniqueSurveyID, Program, GUEventID, GULocationID, EventType, 
         GULatitude, GULongitude, Country, State, Day, Month, Year, 
         ScientificName, Code, NABACode, Count, Duration, Temp, Wind)
IA_IL_MI_data$GULatitude <- as.double(IA_IL_MI_data$GULatitude)
IA_IL_MI_data$GULongitude <- as.double(IA_IL_MI_data$GULongitude)
IA_IL_MI_data$GULocationID <- as.character(IA_IL_MI_data$GULocationID)


IA_IL_MI_zeros <-  IA_IL_MI_zeros %>% 
  mutate(GULocationID = paste(State, ID, sep = "-")) %>% 
  mutate(GUEventID = paste(GULocationID, Survey_ID, sep = "-"),
         EventType = "Pollard",
         Country = "USA",
         Day = lubridate::mday(eventDate),
         Month = lubridate::month(eventDate),
         NABACode = NA) %>% 
  rename(UniqueSurveyID = Survey_ID,
         Year = year,
         Duration = Duration_hours,
         GULatitude = Latitude,
         GULongitude = Longitude) %>% 
  select(UniqueSurveyID, Program, GUEventID, GULocationID, EventType, 
         GULatitude, GULongitude, Country, State, Day, Month, Year, 
         NABACode, Duration, Temp, Wind)
IA_IL_MI_zeros$GULatitude <- as.double(IA_IL_MI_zeros$GULatitude)
IA_IL_MI_zeros$GULongitude <- as.double(IA_IL_MI_zeros$GULongitude)
IA_IL_MI_zeros$GULocationID <- as.character(IA_IL_MI_zeros$GULocationID)

OH_data <- OH_data %>% 
  mutate(State = "OH",
         GULocationID = paste(State, SiteID, sep = "-"),
         GUEventID = paste(GULocationID, SeqID, sep = "-"),
         EventType = "Pollard",
         Country = "USA",
         Day = lubridate::mday(SiteDate),
         Month = lubridate::month(SiteDate),
         Year = lubridate::year(SiteDate),
         Duration = `GU Duration(min)`/60,
         ScientificName = paste(Genus, Species),
         NABACode = NA,
         ID = NA) %>% 
  rename(UniqueSurveyID = SeqID,
         Code = UMDCode,
         GULatitude = `Latitude GU`,
         GULongitude = `Longitude GU`,
         Wind = StartWindMPH,
         Temp = StartTemp,
         Count = Total) %>% 
  select(ID, UniqueSurveyID, Program, GUEventID, GULocationID, EventType, 
         GULatitude, GULongitude, Country, State, Day, Month, Year, 
         ScientificName, Code, NABACode, Count, Duration, Temp, Wind)
OH_data$UniqueSurveyID <- as.character(OH_data$UniqueSurveyID)
OH_data$Wind <- as.character(OH_data$Wind)

OH_zeros <- OH_zeros %>% 
  mutate(State = "OH",
         GULocationID = paste(State, SiteID, sep = "-"),
         GUEventID = paste(GULocationID, SeqID, sep = "-"),
         EventType = "Pollard",
         Country = "USA",
         Day = lubridate::mday(SiteDate),
         Month = lubridate::month(SiteDate),
         Year = lubridate::year(SiteDate),
         Duration = `GU Duration(min)`/60,
         ScientificName = paste(Genus, Species),
         NABACode = NA,
         ID = NA) %>% 
  rename(UniqueSurveyID = SeqID,
         Code = UMDCode,
         GULatitude = `Latitude GU`,
         GULongitude = `Longitude GU`,
         Wind = StartWindMPH,
         Temp = StartTemp,
         Count = Total) %>% 
  select(UniqueSurveyID, Program, GUEventID, GULocationID, EventType, 
         GULatitude, GULongitude, Country, State, Day, Month, Year, 
         NABACode, Duration, Temp, Wind)
OH_zeros$UniqueSurveyID <- as.character(OH_zeros$UniqueSurveyID)
OH_zeros$Wind <- as.character(OH_zeros$Wind)
# combine pollard
pollard_data <- bind_rows(OH_data, IA_IL_MI_data) 

# there's a few zero count days in the observation sheet for some reason. transfer them to zeros
extra_zeros <- filter(pollard_data, Count == 0)

pollard_data <- pollard_data%>% 
  filter(Count > 0)

pollard_zeros <- bind_rows(OH_zeros, IA_IL_MI_zeros, extra_zeros) %>% 
  mutate(Count = 0)

# join with zeros data
pollard_data <- pollard_data %>% 
  bind_rows(pollard_zeros)

rm(extra_zeros, IA_IL_MI_data, IA_IL_MI_Sites, IA_IL_MI_zeros,
   coords, OH_zeros, OH_data, pollard_zeros)
gc()

# are the colnames of the dataframes the same?
compare(names(naba), names(pollard_data))
# no differences combine!
Bfly <- naba %>% 
  bind_rows(pollard_data)

rm(naba, pollard_data)
gc()

# Combine Colias eurytheme (COLEU2) and Colias philodice (COLPHI)
# Clouded and orange sulphurs cannot be reliably identified in the field,
# even by experts. So we lump them into COL-SP
COLSProws1 <- Bfly %>% 
  filter(Code %in% c('COLEU2', 'COLPHI', 'COL-SP')) %>% 
  nrow

Bfly %<>% 
  mutate(Code = case_when(
    Code %in% c('COLEU2', 'COLPHI') ~ 'COL-SP',
    .default = Code
  )) 

COLSProws2 <- Bfly %>% 
  filter(Code == 'COL-SP') %>% 
  nrow

if(COLSProws1 != COLSProws2){
  stop('Orange and Clouded sulphur combination isnt working right')
}

# Ordinal day and week relative to March 1
DayWeek <- read_csv('Data/DayWeek.csv')
# Functional groups
Groups <- read_csv('Data/ButterflyTraitGroups.csv')

# Add new groups to split by family
# Migratory doesn't have many species, so keep it together as one group
Groups %<>% 
  mutate(FamilyGroup = case_when(
    Family %in% c('Hesperiidae', 'Papilionidae') ~ 'HesPap',
    Family %in% c('Lycaenidae', 'Nymphalidae', 
                  'Pieridae', 'Riodinidae') ~ 'LycNymPieRio'
  )) %>% 
  mutate(GroupFamily = paste(Group, FamilyGroup),
         GroupFamily = case_when(
           Group == 'Migratory' ~ 'Migratory',
           .default = GroupFamily
         ))

# Filter data ####
Programs <- c("NFJ", "Ohio Leps", "Iowa Butterfly Survey Network",  
              "Illinois Butterfly Monitoring Network",
              "Michigan Butterfly Network" )

Years <- 1992:2023

Months <- 5:9

States <- c("IA", "IL", "IN", "MI", "MO", "MN", "OH", "WI") 
# keeping MO for now, need to decide to keep later (WL ok with it in 1/10/24)
# extent will be related to states not latitude
# WL 1/10/24: No points in states without a structured survey AND south of the 
# southern tip of IL (~37 degrees)
#LatitudeMin <- 37
#LatitudeMax <- 48

# Set States ####

# Some programs have surveys in other states
# TRUE = Only states with the programs
# FALSE = Include surveys in other states from approved programs
# There's one survey in MN in 2022, but I'm ignoring that one
# We defined the summer breeding grounds to include 545 counties 
# in eight US states (Illinois, Indiana, Iowa, Michigan, Minnesota, 
# Missouri, Ohio and Wisconsin)
# StatesStrict <- FALSE
# 
# if(StatesStrict == TRUE){
#   States <- c('OH', 'IL', 'IA', 'MI', NA)
# } else {
#   States <- c('OH', 'IL', 'IN', 'IA', 'MI', 'WI', 'MN', 'MO',
#               # Ontario, Canada
#               'ON',
#               NA) 
# }

# Summarize data by Subspecies or by Species? ####
SummaryLevel <- 'Species'
# SummaryLevel <- 'SubSpecies'

# NOTE: Currently Canada is removed before processing. 
# Would need county-level troubleshooting to correct.

Canada <- FALSE

#------------------------------------------------------------#
# Can run as source based on above settings ------------------------------------
#------------------------------------------------------------#

# Take a look ####
head(Bfly)
Bfly %>% glimpse
Names %>% glimpse
DayWeek %>% head
Bfly$Program %>% table
Groups %>% head

# Data Processing ####
# Remove unnecessary columns ####
Bfly[,c('Completed', '...1')] <- NULL
# Observer includes people's names. Let's remove this column
Bfly$Observer <- NULL

# Remove NABACode if they are all NAs
if(all(is.na(Bfly$NABACode))){
  Bfly$NABACode <- NULL
}


# Remove duplicated surveys ####
# MWB Note -- I'm not sure what the next 5 lines of code do? 
Bfly %<>% 
  mutate(Check = str_extract(GUEventID, '[:digit:]+$'),
         Same = ifelse(Check == UniqueSurveyID, TRUE, FALSE))
Bfly %<>% 
  filter(Same == TRUE | is.na(Same))

# Add Week values relative to March 1 ####
DayWeekMerge <- DayWeek %>% 
  select(`Day of month`, Month, OrdinalDay, Mar2FebWeeks) %>% 
  rename(Day = `Day of month`)
Bfly %<>% 
  left_join(DayWeekMerge) 

# Standardize OH/Ohio for state ####

# State names in the case study data
# Bfly %>%
#   filter(Program %in% Programs) %>%
#   select(State) %>% table

# Remove Canada 
if(Canada == FALSE){
  Bfly %<>% 
    filter(!Country %in% c('CND', 'Canada', 'CAN'))
}

# Add County names ####
# are there observations without coordinates?
missingCoords <- Bfly %>% 
  filter(is.na(GULongitude) | is.na(GULatitude))
unique(missingCoords$State)
# there are a few observations without coordinates but they all appear to be outside the study region
# removing obs with missing coordinates

Bfly <- Bfly %>% 
  filter(!is.na(GULongitude),
         !is.na(GULatitude)) %>% 
  filter(Country != "MEX") %>% # Mexico is not in study region and won't have state counties, so remove
  mutate(GULongitude = if_else(GULongitude > 0, 
                               true = GULongitude * -1,
                               false = GULongitude))
Bfly <- latlon2counties(Bfly) %>% 
  dplyr::select(-uid, -COUNTYFP,-STATEFP, -NAME, -abb, -polyname) %>% 
  dplyr::rename(state.county = ids)

# There are some NA's. I used google maps to "ground-truth" these points
Bfly %>% 
  filter(is.na(state.county)) %>% 
  select(GULocationID, Program, state.county, GULatitude, GULongitude) %>% 
  unique %>% dim

# there are state counties with NAs! why!?
stateCountyNas <- Bfly %>% 
  filter(is.na(state.county)) %>% 
  distinct(GULocationID, GULongitude, GULatitude)
# plot on map and find out where
us <- rnaturalearth::ne_countries(country = "United States of America", returnclass = "sf") 

ggplot() +
  geom_sf(us, mapping = aes(), fill = NA) +
  geom_point(stateCountyNas, mapping = aes(x = GULongitude, y = GULatitude))
# looks like the occur outside of US border, These could be manually fixed,
# but since they are outside study area, skipping for now and just removing
# obs with state.county == NA
Bfly <- Bfly %>% 
  filter(!is.na(state.county))

# Make indicator variables to match Erin's data
# May or may not be necessary
# Right now, state.county.ind and county.ind are the same. 
# That's ok, but may change if I use Erin's indicators
Bfly %<>% 
  mutate(state.county.ind = state.county %>% as.factor %>% as.numeric,
         county.ind = state.county %>% as.factor %>% as.numeric,
         site.ind = GULocationID %>% as.factor %>% as.numeric) 

# Check states for is.na(State) ####
Bfly %>% 
  filter(is.na(State)) %>% 
  select(Program, state.county) %>% 
  unique 

# Fill in the values assuming there are surveys in that county and
# the survey is within the expected state
#Bfly %<>% 
#  FillState(., StateCode = 'IL') %>% 
#  FillState(., StateCode = 'MI') %>% 
#  FillState(., StateCode = 'IA') %>% 
#  FillState(., StateCode = 'OH') 

# Calculate if there are still NAs in state
StateNA <- is.na(Bfly$State) %>% sum == 0 

# Throw a stop if there are still NAs in states
if(StateNA != TRUE){
  stop('Some states are NAs')
}

# Filter for Analyses ####
MidwestData <- Bfly %>% 
  filter(Program %in% Programs,
         Year %in% Years,
         Month %in% Months,
         State %in% States) 

# Proportion of surveys before 1992 for Methods
DateProportion <- Bfly %>% 
  filter(Program %in% Programs,
         Month %in% Months,
         State %in% States)
Pre <- DateProportion %>% 
  filter(Year < 1992) %>%
  # If using 1993, there is 1.3% of records before that date
  # filter(Year < Years[1]) %>% 
  select(GUEventID) %>% 
  n_distinct
All <- DateProportion %>% 
  select(GUEventID) %>% 
  n_distinct
PercentYear <- DateProportion %>% 
  select(GUEventID, Year) %>% 
  group_by(Year) %>% 
  summarize(nSurveysYear = n_distinct(GUEventID),
            nPropYear = nSurveysYear / All) %>% 
  mutate(CumulativePercent = cumsum(nPropYear))
ProportionPreDate <- Pre / All

# Check for NA's in States# Check for NA's in state.county
StateCountyNA <- MidwestData %>% 
  filter(is.na(state.county)) %>% 
  select(GULocationID, Program, state.county, GULatitude, GULongitude) %>% 
  unique

# Redo the state.county.ind and county.ind 
# These are used when creating AlignedTo Erins
MidwestData %<>% 
  mutate(state.county.ind = state.county %>% as.factor %>% as.numeric,
         county.ind = state.county %>% as.factor %>% as.numeric)


# Put an error if there are any NAs in state.county.ind
SCNA <- MidwestData %>% 
  filter(Country == 'USA', 
         is.na(state.county.ind)) %>% 
  sum == 0

if(SCNA != TRUE){
  MidwestData %>% 
    filter(is.na(state.county)) %>% 
    select(GULocationID, Program, state.county.ind, GULatitude, GULongitude) %>% 
    unique
  stop('Some state.county values for US sites are NA')
}

# FIGURE OUT CANADA STATE.COUNTY ENTRIES
# Not planning to include Canada at any point (1/10/24)
# MidwestData$state.county %>% table
# MidwestData %>% filter(Country == 'CND')

# Filter out survey and species case study data
MidwestSurveys <- MidwestData %>% 
  select(Program, GUEventID, GULocationID, EventType, 
         GULatitude, GULongitude, Country, State,
         Day, Month, Mar2FebWeeks, Year, Duration, Temp, Wind, OrdinalDay,
         state.county, state.county.ind, county.ind, site.ind) %>% 
  unique

# Code to look for the species with >=50 obs and see if there's any new ones
# to add to our analyses (1/26/2024)
# MidwestData %>% 
#   group_by(Code) %>% 
#   summarize(Sum = sum(Count)) %>% 
#   filter(Sum >= 50) %>% 
#   left_join(Groups) %>% 
#   write_csv(file = 'Output/NewCounts.csv')

# Cut down to only species we'll work with
MidwestDataSpp <- MidwestData %>% 
  filter(Code %in% Groups$Code)
# All species in Groups$Code are in MidwestData$Code
# Except for the other orange/clouded sulphur, which 
# have already been collapsed into COL-SP
# Groups$Code[!Groups$Code %in% SppMW]

# Look at state breakdown by program
MidwestSurveys %>% 
  group_by(Program, State) %>% 
  summarize(N = n())


# Check for same event or location IDs within dataset ####
# Location IDs
#LocIDUnique <- MidwestSurveys %>% 
#  select(Program, GULocationID) %>% 
#  unique %>%
#  ungroup %>% 
#  group_by(GULocationID) %>% 
#  filter(n() > 1) 
# Event IDs
EventIDUnique <- MidwestSurveys %>% 
  select(Program, GUEventID) %>% 
  unique %>%
  ungroup %>% 
  group_by(GUEventID) %>% 
  filter(n()>1) 
if(#nrow(LocIDUnique) > 0 |
  nrow(EventIDUnique) > 0){
  stop('Duplicated location or event IDs')
}

# Some surveys have multiple recordings of a species ####
Spp <- MidwestDataSpp$Code %>% unique
DF <- tibble(Program = NA, GUEventID = NA, ScientificName = NA,
             Code = NA, Count = NA)

for(ss in 1:length(Spp)){
  Events <- MidwestDataSpp %>% 
    filter(Code == Spp[ss]) %>% 
    group_by(GUEventID) %>% 
    filter(n()>1) %>% 
    select(Program, GUEventID, ScientificName, Code, Count)
  DF %<>% bind_rows(Events)
}
# Remove first row
DF %<>% filter(!is.na(Program))

# Number of affected surveys/counts
Affected <- DF %>% 
  group_by(Code) %>% 
  summarize(NumSurveys = n_distinct(GUEventID),
            NumRows = n(),
            NumCount = sum(Count),
            MeanCount = mean(Count)) %>% 
  mutate(RemoveRows = NumRows - NumSurveys)
Affected

# How many rows should we end up with after summing these
ExpectedRowsByCode <- nrow(MidwestDataSpp) - 
  sum(Affected$RemoveRows)

# Sum up these multiple entries by EventID 
NRowByCode <- MidwestDataSpp %>% 
  # Currently keeps subspecies separate (two CELLAD)
  group_by(Program, GUEventID, Code) %>% 
  summarize(Count = sum(Count)) %>% 
  nrow

# Expected and actual rows should line up (TRUE)
OK <- ExpectedRowsByCode == NRowByCode
if(OK != TRUE){
  stop("Expected rows don't line up with actual rows")
}

# Sum up multiple entries by SummaryLevel (set earlier)
if(SummaryLevel == 'Species'){
  MidwestDataSpp %<>%
    group_by(Program, GUEventID, GULocationID, 
             EventType, GULatitude, GULongitude,
             state.county, state.county.ind, county.ind, site.ind,
             Country, State, Day, Month, Year, 
             Code,
             Duration, Temp, Wind, OrdinalDay, Mar2FebWeeks) %>% 
    summarize(Count = sum(Count)) %>% 
    ungroup
} else {
  MidwestDataSpp %<>% 
    mutate(OriginalCode = Code,
           Code = ScientificName) %>% 
    group_by(Program, GUEventID, GULocationID, 
             EventType, GULatitude, GULongitude,
             state.county, state.county.ind, county.ind, site.ind,
             Country, State, Day, Month, Year, 
             Code, OriginalCode,
             Duration, Temp, Wind, OrdinalDay, Mar2FebWeeks) %>% 
    summarize(Count = sum(Count)) %>% 
    ungroup
}

# Fill in Zeros for species that weren't observed ####
MidwestDataSpp %>% 
  select(GUEventID, Code) %>%
  table() %>% 
  as.data.frame %>% 
  filter(Freq != 1) %>% 
  select(Freq) %>% 
  table
# Lots of zeros, none larger than that
# Throw an error if there are more than one row of species in an event
Problem <- MidwestDataSpp %>% 
  select(GUEventID, Code) %>%
  table() %>% 
  as.data.frame %>% 
  filter(Freq > 1) %>% 
  dim
if(Problem[1] > 0){
  stop('Multiple rows of a species during a survey')
}

# # Grab the unique surveys
# # Note that a survey without any of the selected species won't be here
# # Haven't determined yet if there are other surveys that didn't
# # have those species compared to surveys that I don't have observations
# # for
# Uniques <- MidwestSurveys$GUEventID %>% unique
# # Grab the unique species (not CELLAD twice)
# EachSpp <- Univoltine$Code %>% unique
# # Create a dataframe to hold the survey/spp combos that aren't recorded
# AddTo <- data.frame(Code = NA,
#                     GUEventID = NA,
#                     Count = NA)
# start_time <- Sys.time()
# # Loop through the unique surveys
# for(uu in 1:length(Uniques)){
#   # Pull out just one survey
#   US <- Uniques[uu]
#   # Pull out the observations for that survey
#   CSUS <- MidwestDataSpp %>% filter(GUEventID == US)
#   # Loop through the species
#   for(ss in 1:length(EachSpp)){
#     # Check if the species is already present
#     if(!(EachSpp[ss] %in% CSUS$Code)){
#       # If the species isn't present, add the spp name, survey, and count of 0 
#       # to the holder data frame
#       AddTo %<>% add_row(Code = EachSpp[ss],
#                          GUEventID = Uniques[uu],
#                          Count = 0)
#     }
#   }
#   print(uu)
# }
# end_time <- Sys.time()
# end_time - start_time

# # Filter out the first row (NA placeholder) and rename column
# # Original column name didn't seem to work in the loop, but it may
# # have been due to something else. This was an easy solution anyways though
# AddTo %<>% 
#   filter(!is.na(GUEventID))

# # Add the AddTo data to the original data frame and check if
# # all species/surveys have data. All should be 1
# MidwestDataSpp %>% 
#   bind_rows(AddTo) %>% 
#   select(GUEventID, Code) %>% table %>% table
# 
# Expand.grid approach
# Uses MidwestSurveys, so don't need to worry about a subset of species
# not being observed in a survey
AllZeros <- expand_grid(GUEventID = MidwestSurveys$GUEventID %>% unique,
                        Code = MidwestDataSpp$Code %>% unique,
                        Count = 0)
# # Add the AllZeros data to the original data frame and check if
# # all species/surveys have data. All should be 1 or 2
MidwestDataSpp %>%
  bind_rows(AllZeros) %>%
  select(GUEventID, Code) %>% table %>% table

# For each species, how many surveys was it observed on? 
# And how many were observed?
MidwestDataSpp %>% 
  group_by(Code) %>% 
  summarize(Events = n_distinct(GUEventID),
            Count = sum(Count)) %>% 
  print(n=length(MidwestDataSpp$Code %>% unique))

# Add rows and save
# Add survey information 
AllZeros %<>% left_join(MidwestSurveys)

# Join observations with zeros ####
MidwestDataSpp %<>% bind_rows(AllZeros)

# Add the zeros and observations together
MidwestDataSpp %<>%
  group_by(Program, GUEventID, GULocationID, 
           EventType, GULatitude, GULongitude,
           state.county, state.county.ind, county.ind, site.ind,
           Country, State, Day, Month, Year, 
           Code,
           Duration, Temp, Wind, 
           OrdinalDay,
           Mar2FebWeeks) %>% 
  summarize(Count = sum(Count)) %>% 
  ungroup

# Should all have ones
Solution <- MidwestDataSpp %>% 
  select(GUEventID, Code) %>% table %>% table %>% 
  as.data.frame()
if(nrow(Solution) > 1 |
   Solution$Freq[1] != nrow(MidwestSurveys) * n_distinct(MidwestDataSpp$Code)){
  stop("Adding zero's didn't work right")
}

# Add Family here
MidwestDataSpp %<>% 
  left_join(Names, join_by(Code == UMD_Code))

# Only write data if there have been changes to this document
this_info <- file.info('Code/R/FormatDataCh1.R')
data_info <- file.info('Data/CleanedData/DataNewNoPreds.csv')
# Calculate if code is older than data object 
# If positive, code is more recent (later time) and data needs to be rewritten 
CheckData <- (this_info$mtime - data_info$mtime)
if(CheckData > 0){
  data.table::fwrite(x = MidwestDataSpp,
                     file = "Data/CleanedData/DataNewNoPreds.csv")
}

# Format open and crop cover covariates ####
# percent open and crop covariates are static - do not change among years.
# Duration does though

# Join county and site level covariates
# Covariates
# Updated to exclude urban and water from percent open
# CountyCov <- read_csv("Data/Covariates_County_expanded.csv")
# SiteCov <- read_csv('Data/totalSites_propOpen.csv')
CountyCov <- read_csv("Data/CountyCovariates_percOpenNoWaterNoUrb.csv") %>% 
  select(ID, propLC)
SiteCov <- read_csv('Data/totalSites_propOpen_noWaternoUrb.csv')

SiteCov %<>% 
  mutate(perc.open.site = propOpen*100) %>% 
  select(GULocationID, perc.open.site)
# SiteCov %<>% 
#   mutate(perc.open.site = propOpen*100) %>% 
#   select(GULocationID, perc.open.site)
# CountyCov %<>% 
#   rename(perc.open.county = percOpen,
#          perc.crop.county = percCrop,
#          state.county = ID)
CountyCov %<>% 
  mutate(perc.open.county = propLC * 100,
         state.county = ID) %>% 
  select(state.county, perc.open.county)

# do we have all county covariates?
uniqueCounties <- distinct(MidwestDataSpp, state.county)
uniqueCounties <- left_join(uniqueCounties, CountyCov)
# yep all county covariates are gathered

# now looking at sites-specific variates
head(SiteCov)

uniqueSites <- distinct(MidwestDataSpp, GULocationID, .keep_all = TRUE) %>% 
  select(GULocationID, state.county, GULatitude, GULongitude)
uniqueSites <- left_join(uniqueSites, SiteCov)
# all sites have covariates 

# join midwest data with covariates
MidwestDataSpp <- left_join(MidwestDataSpp, CountyCov, by = "state.county")
MidwestDataSpp <- left_join(MidwestDataSpp, SiteCov, by = "GULocationID")


CovNA <- MidwestDataSpp %>% 
  filter(is.na(perc.open.county))

if(nrow(CovNA) != 0){
  stop("There are NA's in county covariates")
}

# Effort
# There's some points with an effort of 0. Replace these values with the median 
# value of the non-zero points
# How many have 0's
# Save as an object so that it can be referenced in Methods.Rmd
DurationZero <- MidwestDataSpp %>% 
  select(GUEventID, Duration) %>% 
  distinct %>% 
  filter(Duration == 0) %>% 
  nrow %>% 
  divide_by(n_distinct(MidwestDataSpp$GUEventID))
# nrow(DurationZero) / n_distinct(MidwestDataSpp$GUEventID)  

# Calculate Median to substitute for zeros
DurationMedian <- MidwestDataSpp$Duration[MidwestDataSpp$Duration > 0] %>% 
  median(na.rm = TRUE) # 1.16
# Some values are huge, most are 1.5 or less
MidwestDataSpp$Duration %>% summary
MidwestDataSpp %>% 
  filter(Duration > 180) %>%
  select(GUEventID, Duration) %>% 
  distinct
# Effort by program (uncomment to look at mean/max before capping effort)
# MidwestDataSpp %>% 
#   group_by(Program) %>% 
#   summarize(meanDuration = mean(Duration),
#             maxDuration = max(Duration))

# replace large NABA (NFJ) values with max value of 15 party-hours 
#   NABA has min of 6 party-hours and 4 people in recent years
# replace large Pollard values with max value of 5 hours, 
# replace 0s with median 1.16
MidwestDataSpp <- MidwestDataSpp %>% 
  mutate(Duration = case_when(
    Program == 'NFJ' & Duration > 15 ~ 15,
    Program != 'NFJ' & Duration > 5 ~ 5,
    Duration == 0 ~ DurationMedian,
    .default = Duration
  ))

# Uncomment to see new values of effort, with caps
# MidwestDataSpp %>% 
#   group_by(Program) %>% 
#   summarize(meanDuration = mean(Duration),
#             maxDuration = max(Duration))

# Make comparable to Erin's ####

# Erin's columns: usiteID, program, state.county, lat, long, yr, wk, monarch, duration, site.ind, county.ind, perc.open
ColumnCompare <- tibble(
  Erins = c('usiteID', 'program', 'state.county', 
            'lat', 'long', 'yr', 'wk', 
            'monarch', 'duration',
            'site.ind', 'county.ind', 'perc.open'),
  Mine = c('GULocationID', 'Program', 'state.county.ind',
           'GULatitude', 'GULongitude', 'Year', 'Mar2FebWeeks',
           'Count', 'Not included yet', 
           'site.ind', 'county.ind', 'Not included yet'))
ColumnCompare

dataNew <- MidwestDataSpp %>% 
  filter(Country != 'CND') %>% 
  select(GULocationID, Program, Country, state.county.ind,
         state.county,
         GULatitude, GULongitude,
         Year, Mar2FebWeeks, Count,
         perc.open.county, #perc.crop.county, 
         perc.open.site,
         Duration,
         site.ind, county.ind, Family, Code, GUEventID) %>% 
  rename(usiteID = GULocationID, 
         program = Program,
         state.county.name = state.county,
         state.county = state.county.ind,
         lat = GULatitude,
         long = GULongitude,
         yr = Year,
         wk = Mar2FebWeeks,
         count = Count,
         species = Code,
         country = Country,
         # Unstandardized covariates
         # crop_county_uns = perc.crop.county,
         open_county_uns = perc.open.county,
         open_site_uns = perc.open.site,
         effort_uns = Duration)

# Look for species with less data
NotEnough <- dataNew %>% 
  group_by(species) %>% 
  summarize(Observations = sum(count)) %>% 
  arrange(Observations) %>% 
  filter(Observations < 50)
# Northern metalmark (CALBOR) only has 2 observations. 
# Henry's elfin (CALLHEN) has 48 observations.
# Remove both for now, but CALLHEN might be enough
# Names %>% 
#   filter(UMD_Code %in% c('CALBOR', 'CALLHEN'))

# Remove the species with 50 or fewer observations
dataNew %<>% 
  filter(!species %in% NotEnough$species)

# Note: as of 1/10/24, all species have enough data

# Generic data prep for fitting the models ----------------------------------
# Any processing that depends on the number of species in the model is done
# in the TryModelICM.R code

# Change site indicator to start at 1
SiteInd <- tibble(site.ind = dataNew$site.ind %>% unique,
                  site_ind = 1:length(unique(dataNew$site.ind)))
dataNew %<>% left_join(SiteInd) 

# Standardize crop_county, open_county, open_site
dataNew %<>% 
  mutate(#crop_county = (crop_county_uns - mean(crop_county_uns)) /
           # sd(crop_county_uns),
         open_county = (open_county_uns - mean(open_county_uns)) / 
           sd(open_county_uns),
         open_site = (open_site_uns - mean(open_site_uns, na.rm = TRUE)) / 
           sd(open_site_uns, na.rm = TRUE),
         open_site = ifelse(is.na(open_site), 0, open_site)
  )

# Standardize effort
# There are some points with a value of 0 for effort. All of them have at least
# one butterfly recorded. Replace the zeros with the mean value that isn't 0.
# All 5 programs have some NA's
# dataNew %>% 
#   group_by(program) %>% 
#   summarize(NotNA = (!is.na(effort_uns)) %>% sum,
#             IsNA = (is.na(effort_uns)) %>% sum)
# 

dataNew %<>% 
  mutate(effort = (effort_uns / mean(effort_uns, na.rm = TRUE)),
         # Replace NA's with 1
         effort = ifelse(is.na(effort),
                         1,
                         effort)) 

# save CleanedData (earlier in code, with check of whether code has changed 
# since CleanedData was last saved)
#fwrite(x = dataNew, file = "Data/CleanedData/DataNew.csv")