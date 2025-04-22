## Data

### Butterfly observation data
Proprietary; see below for data availability statement
### Butterfly observation data
Proprietary; see below for data availability statement
1. [NABAAllData_1977toOct2022.csv](https://github.com/zipkinlab/Archived-data/blob/master/Leuenberger_etal_2025_PNAS/NABAAllData_1977toOct2022.csv): Data from NABA surveys. Includes survey information (location, date, GPS coordinates, party hours [time spent surveying x number of people]) as well as butterflies observed (species name, species code, number observed). 
2. [IA IL MI Pollardbase data through 10.04.2023.xlsx](https://github.com/zipkinlab/Archived-data/blob/master/Leuenberger_etal_2025_PNAS/IA%20IL%20MI%20Pollardbase%20data%20through%2010.04.2023.xlsx): Butterfly surveys (location, date, GPS coordinates, duration) and observations (species name, species code, number observed) from the Illinois Butterfly Monitoring Network, the Iowa Butterfly Survey Network, and the Michigan Butterfly Network. 
3. [Pb-route-latlon-IA-IL-MI.xlsx](https://github.com/zipkinlab/Archived-data/blob/master/Leuenberger_etal_2025_PNAS/Pb-route-latlon-IA-IL-MI.xlsx): Survey-level information (route, state, county, GPS location) from the Illinois Butterfly Monitoring Network, the Iowa Butterfly Survey Network, and the Michigan Butterfly Network
4. [GU Ohio 2022 Flat File All Observations.xlsx](https://github.com/zipkinlab/Archived-data/blob/master/Leuenberger_etal_2025_PNAS/GU%20Ohio%202022%20Flat%20File%20All%20Observations.xlsx): Butterfly surveys (location, date, GPS coordinates, duration) and observations (species name, species code, number observed) from the Ohio Lepidopterists. 

### Covariate data
1. [CountyCovariates_percOpenNoWaterNoUrb.csv](Data/CountyCovariates_percOpenNoWaterNoUrb.csv): County-level data on proportion of open land, excluding water and urban land covers.
2. [totalSites_propOpen_noWaternoUrb.csv](Data/totalSites_propOpen_noWaternoUrb.csv): Site-level data on proportion of open land, excluding water and urban land covers. Includes location ID and program.

### Miscellaneous data
1. [BflyNames.csv](Data/BflyNames.csv): Contains family, NABA Code (`UMD_Code`), scientific name, common name, and additional information about each species in the NABA dataset. 
2. [ButterflyTraitGroups.csv](Data/ButterflyTraitGroups.csv): Initial functional trait groupings for butterfly species, based on LepTraits and expert knowledge. Split species into four groups: migratory, resident univoltine, resident generalist, and resident specialist. Includes a count of observations, but those numbers are not the final count. 
3. [sbus.traits_w_disturb.nov2023.csv](Data/sbus.traits_w_disturb.nov2023.csv): Trait data from [Edwards et al. 2025](https://www.science.org/doi/10.1126/science.adp4671), used for overwintering (diapause) data. 
4. [TraitsButterflyCrossley.csv](https://github.com/mcrossley3/NorthAmericanButterflies/tree/master): Trait data from [Crossley et al. 2025](https://onlinelibrary.wiley.com/doi/10.1111/gcb.15582), used to check a couple of overwintering stages (diapause).
5. [CountyKey.RData](Data/CountyKey.RData): County indicator and county name, used to keep track of which county is which during post-processing. Originally created in [TryModel_spAbundance.R](Code/R/TryModel_spAbundance.R), line 214 (commented out)