# Three decades of declines restructure butterfly communities in the Midwestern United States

### Wendy Leuenberger, [Jeffrey W. Doser](https://www.jeffdoser.com/), Michael W. Belitz, Leslie Ries, Nick M. Haddad, Wayne E. Thogmartin, and [Elise F. Zipkin](https://zipkinlab.org/)

### In revision

### Please contact the first author for questions about the code or data: Wendy Leuenberger (leuenbe9@msu.edu)

-------------------------------------------------------------------------------

## Abstract:

Insects are declining worldwide. These declines have been documented across taxonomic groups and are worrisome given ecosystem services provided by insects. Long-term data have illuminated butterfly declines across geographic regions. However, critical questions remain on how butterfly declines are distributed across species and functional groups, limiting effective conservation. Here we show unprecedented changes in butterfly biodiversity resulting from 32 years of species-levels declines throughout the Midwestern United States. No species increased over the three-decade study period and abundance declined across every functional group (e.g., rare, common, migratory, resident; annual mean trend: -0.9 to -2.3% per year). Species richness declined across all but one functional group, with concomitant increases in evenness (e.g., abundance among species) in several groups resulting from steeper losses in abundance for common species (abundance: -1.9% per year; richness: -0.04% per year) as compared to rare species (abundance: -0.9% per year; richness: -1.33% per year). Our results paint a bleaker picture than other butterfly studies (Breed et al. 2013, Halsch et al. 2021, Edwards et al. 2025), likely due to our long time series of data and ability to include rare species. Such widespread declines undoubtedly affect other trophic levels and ecosystem services. Our results indicate that focusing risk assessments and management interventions only on rare species is insufficient given broad declines across species, which have fundamentally restructured butterfly communities in the region. As such, conservation efforts should shift focus to species assemblages and entire communities when possible.


-------------------------------------------------------------------------------

## Repository Directory

1. [Code](Code): Contains all code files used for the project
2. [Data](Data): Contains publicly available data files use for the project

-------------------------------------------------------------------------------

## R Code and Workflow

1. Data preparation: 
    1. [FormatDataCh1.R](Code/R/FormatDataCh1.R): R code used to take raw data and create cleaned objects ready for analyses. Output is used in TryModel_spAbundance.R 
    2. [SpeciesFunctionalTraits.R](Code/R/SpeciesFunctionalTraits.R): R code used to link species with their traits. Used in Abundances.R, Diversity.R, Richness.R, and MultiMetricFigures.R
2. Model running: 
    1. [TryModel_spAbundance.R](Code/R/TryModel_spAbundance.R): R code used to select a specific group of butterflies for analyses, perform some final data preparation, and begin the model fitting process. Models take a long time to finish, so these models are continued in UpdateModel_spAbundance.R
    2. [UpdateModel_spAbundance.R](Code/R/UpdateModel_spAbundance.R): R code used to continue the model runs after TryModel_spAbundance.R.
3. Post processing:
    1. [PostProcess_MsAbund.R](Code/R/PostProcess_MsAbund.R): Code that takes several model runs, puts them together, checks convergence, and does some initial abundance and richness processing. This code is run separately for each group that is used for modeling. These outputs aren't ready for interpretation; needs to be combined across groups in the next step.
4. Results: 
    1. [Abundances.R](Code/R/Abundances.R): Calculates species-specific, group, and overall patterns in abundance. Includes code to make figure 3 (abundance trends). 
    2. [Richness.R](Code/R/Richness.R): Calculates species-specific, group, and overall patterns in richness. 
    3. [Diversity.R](Code/R/Diversity.R): Calculates species-specific, group, and overall patterns in evenness (Shannon diversity). 
5. Synthesis:
    1. [MultiMetricFigures.R](Code/R/MultiMetricFigures.R): Takes output from Abundances.R, Richness.R, and Diversity.R to make composite figures (Figure 2) and tables. 
    2. [Map_Figure1.R](Code/R/Map_Figure1.R): Takes survey data and creates a map of the number of surveys per location (Figure 1). 
6. Batch files for super computer:
    Note: likely specific to MSU's ICER HPCC at the time of running. It might be a helpful starting place but wouldn't be ready to use as-is. 
    1. [run_BCSmodel_msAbund.sb](run_BCSmodel_msAbund.sb): Runs the code in [TryModel_spAbundance.R](R/TryModel_spAbundance.R). Runs one chain; run command 5 times to run 5 chains.
    2. [run_UpdateModel.sb](run_UpdateModel.sb): Runs the code in [UpdateModel_spAbundance.R](R/UpdateModel_spAbundance.R). Updates the specified chain only. 
    3. [run_PostProcess.sb](run_PostProcess.sb): Runs the code in [PostProcess_MsAbund.R](R/PostProcess_MsAbund.R)

*Note that additional folders (e.g. 'Output', 'Output/1_runs') would be necessary to save output from the Post-process, Results, and Synthesis steps, as specified in the code.*

## Data

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
4. [TraitsButterflyCrossley.csv](https://github.com/zipkinlab/Leuenberger_etal_2025_PNAS/blob/main/Data/TraitsButterflyCrossley.csv): Trait data from [Crossley et al. 2025](https://onlinelibrary.wiley.com/doi/10.1111/gcb.15582), [GitHub repo](https://github.com/mcrossley3/NorthAmericanButterflies/tree/master), used to check a couple of overwintering stages (diapause).
5. [CountyKey.RData](Data/CountyKey.RData): County indicator and county name, used to keep track of which county is which during post-processing. Originally created in [TryModel_spAbundance.R](Code/R/TryModel_spAbundance.R), line 214 (commented out)

## Data availability statement:
Covariate data and data regarding butterfly names and groupings are available at this GitHub repository: [https://github.com/zipkinlab/Leuenberger_etal_2025_PNAS](https://github.com/zipkinlab/Leuenberger_etal_2025_PNAS). Butterfly data are proprietary and were obtained from the North American Butterfly Association (https://www.naba.org/), the Iowa Butterfly Survey Network (https://www.reimangardens.com/collections/insects/iowa-butterfly-survey-network/), the Illinois Butterfly Monitoring Network (https://bfly.org/), the Michigan Butterfly Network (https://michiganbutterfly.org/) and the Ohio Lepidopterists (http://www.ohiolepidopterists.org/). These data may be available upon reasonable request to LR and with permission from the aforementioned programs.

## Code availability statement: 
Code needed to run analyses (R scripts) is available at this GitHub repository: [https://github.com/zipkinlab/Leuenberger_etal_2025_PNAS](https://github.com/zipkinlab/Leuenberger_etal_2025_PNAS).
