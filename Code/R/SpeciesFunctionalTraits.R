## Functional groups ####
Groups <- read_csv(paste0(path, 'Data/ButterflyTraitGroups.csv'))
# Overwintering group
# Data are from LepTraits, though accessed through the SBUS paper
# Overwinter <- read_csv(paste0(path, 'Data/traits_w_disturb.nov2023.csv'))
Overwinter <- read_csv('Data/traits_w_disturb.nov2023.csv') %>% 
  rename(Code = code)
Groups %<>% 
  left_join(Overwinter[,c('Code', 'diapause')])
# Species without overwintering codes
NoOverwinter <- Groups %>% 
  filter(is.na(diapause))
# Some have codes in Mike's other file
# Overwinter2 <- read_csv(paste0(path, 'Data/TraitsButterfly.csv'))
Overwinter2 <- read_csv('Data/TraitsButterfly.csv') %>% 
  rename(Diapause = `overwintering stage`,
         ScientificName = 'Species')
Overwinter2$Diapause %<>% str_replace('P/L', 'PL')
Groups %<>% 
  left_join(Overwinter2[,c('ScientificName', "Diapause")])

# Note: diapause recordings aren't the same
Groups %>% 
  select(diapause, Diapause) #%>% View

Groups %<>%
  mutate(Overwinter = ifelse(
    # If the diapause column from Mike's SBUS paper is NA
    is.na(diapause), 
    # Then grab the value from Mike's other dataset
    Diapause, 
    # Otherwise (not NA), use the one from Mike's SBUS paper
    diapause))
# Add in a couple of NA's
Groups %<>% 
  mutate(Overwinter = case_when(
    Code == 'PHYBAT' ~ 'L', # https://www.butterfliesandmoths.org/species/Phyciodes-batesii
    Code == 'NASLHE' ~ 'P', # https://bugguide.net/node/view/27931 
    Code == 'HEMISO' ~ 'A', # Discussion with Mike
    .default = Overwinter
  ))

# Make csv to discuss with Mike
OverwinterDecision <- Groups %>% 
  select(Code, ScientificName, diapause, Diapause, Overwinter) %>% 
  rename(diapause_SBUS = diapause,
         diapause_2 = Diapause,
         Decision = Overwinter) 
# write_csv(OverwinterDecision, 'Data/OverwinterDecision.csv')
# Take out the intermediate columns
Groups %<>% 
  select(-diapause, -Diapause)

# Remove the species that were collapsed into COL-SP 
# They can't be identified in the field
Groups %<>% 
  filter(!Code %in% c('COLEU2', 'COLPHI'))

# Supplement Table S7 - functional groupings
# Groups %>% 
#   mutate(Prevalence = case_when(Code %in% SpeciesLists[['Rare']] ~ 'Rare',
#                                 Code %in% SpeciesLists[['Common']] ~ 'Common',
#                                 .default = '')) %>% 
#   arrange(ScientificName) %>% 
#   select(ScientificName, Prevalence) %>% 
#   write.csv('Output/12_AbundanceTrends/Prevalence.csv',
#     row.names = FALSE)
