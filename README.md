# Three decades of declines restructure butterfly communities in the Midwestern U.S.

### Wendy Leuenberger, [Jeffrey W. Doser](https://www.jeffdoser.com/), Michael W. Belitz, Leslie Ries, Nick M. Haddad, Wayne E. Thogmartin, and [Elise F. Zipkin](https://zipkinlab.org/)

### In review

### Please contact the first author for questions about the code or data: Wendy Leuenberger (leuenbe9@msu.edu)

-------------------------------------------------------------------------------

## Abstract:

Insects are declining worldwide (Wepprich et al. 2019, van Klink et al. 2020, Forister et al. 2021). These declines have been documented across taxonomic groups and are worrisome given ecosystem services provided by insects (Wilson 1987, Kawahara et al. 2021). Long-term data have illuminated butterfly declines across geographic regions (Warren et al. 2021, Edwards et al. In press). However, critical questions remain on how butterfly declines are distributed across species and functional groups, limiting effective conservation. Here we show unprecedented changes in butterfly biodiversity resulting from 32 years of species-levels declines throughout the Midwestern United States. No species increased over the three-decade study period and abundance declined across every functional group (e.g., rare, common, migratory, resident). Species richness declined across all but one functional group, with concomitant increases in evenness (e.g., abundance among species) in several groups resulting from steeper losses in abundance for common species as compared to rare species. Our results paint a bleaker picture than other butterfly studies (Breed et al. 2013, Halsch et al. 2021, Edwards et al. In press), likely due to our long time series of data and ability to include rare species. Such widespread declines undoubtedly affect other trophic levels and ecosystem services. Our results indicate that focusing risk assessments and management interventions only on rare species is insufficient given broad declines across species, which have fundamentally restructured butterfly communities in the region. As such, conservation efforts should shift focus to species assemblages and entire communities when possible (Belitz et al. In press).

-------------------------------------------------------------------------------

## Repository Directory

1. [Code](Code): Contains all code files used for the project

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

## Data availability statement:
Covariate data will be made available before publication and is available upon request during peer review. Butterfly data are proprietary and were obtained from the North American Butterfly Association (https://www.naba.org/), the Iowa Butterfly Survey Network (https://www.reimangardens.com/collections/insects/iowa-butterfly-survey-network/), the Illinois Butterfly Monitoring Network (https://bfly.org/), the Michigan Butterfly Network (https://michiganbutterfly.org/) and the Ohio Lepidopterists (http://www.ohiolepidopterists.org/). These data may be available upon reasonable request to LR and with permission from the aforementioned programs.

## Code availability statement: 
Code needed to run analyses (R scripts) is available at this GitHub repository: https://github.com/wleuenberger/ButterflyCommunityTrends.
