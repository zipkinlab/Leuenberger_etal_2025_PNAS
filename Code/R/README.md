## R Code and Workflow

1. Data preparation: 
    1. [FormatDataCh1.R](FormatDataCh1.R): R code used to take raw data and create cleaned objects ready for analyses. Output is used in TryModel_spAbundance.R 
    2. [SpeciesFunctionalTraits.R](SpeciesFunctionalTraits.R): R code used to link species with their traits. Used in Abundances.R, Diversity.R, Richness.R, and MultiMetricFigures.R
2. Model running: 
    1. [TryModel_spAbundance.R](TryModel_spAbundance.R): R code used to select a specific group of butterflies for analyses, perform some final data preparation, and begin the model fitting process. Models take a long time to finish, so these models are continued in UpdateModel_spAbundance.R
    2. [UpdateModel_spAbundance.R](UpdateModel_spAbundance.R): R code used to continue the model runs after TryModel_spAbundance.R.
3. Post processing:
    1. [PostProcess_MsAbund.R](PostProcess_MsAbund.R): Code that takes several model runs, puts them together, checks convergence, and does some initial abundance and richness processing. This code is run separately for each group that is used for modeling. These outputs aren't ready for interpretation; needs to be combined across groups in the next step.
4. Results: 
    1. [Abundances.R](Abundances.R): Calculates species-specific, group, and overall patterns in abundance. Includes code to make figure 3 (abundance trends). 
    2. [Richness.R](Richness.R): Calculates species-specific, group, and overall patterns in richness. 
    3. [Diversity.R](Diversity.R): Calculates species-specific, group, and overall patterns in evenness (Shannon diversity). 
5. Synthesis:
    1. [MultiMetricFigures.R](MultiMetricFigures.R): Takes output from Abundances.R, Richness.R, and Diversity.R to make composite figures (Figure 2) and tables. 
    2. [Map_Figure1.R](Map_Figure1.R): Takes survey data and creates a map of the number of surveys per location (Figure 1). 

*Note that additional folders (e.g. 'Output', 'Output/1_runs') would be necessary to save output from the Post-process, Results, and Synthesis steps, as specified in the code.*