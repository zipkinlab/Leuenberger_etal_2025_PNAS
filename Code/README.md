## R Code and Workflow

1. Data preparation: 
    1. [FormatDataCh1.R](R/FormatDataCh1.R): R code used to take raw data and create cleaned objects ready for analyses. Output is used in TryModel_spAbundance.R 
    2. [SpeciesFunctionalTraits.R](R/SpeciesFunctionalTraits.R): R code used to link species with their traits. Used in Abundances.R, Diversity.R, Richness.R, and MultiMetricFigures.R
2. Model running: 
    1. [TryModel_spAbundance.R](R/TryModel_spAbundance.R): R code used to select a specific group of butterflies for analyses, perform some final data preparation, and begin the model fitting process. Models take a long time to finish, so these models are continued in UpdateModel_spAbundance.R
    2. [UpdateModel_spAbundance.R](R/UpdateModel_spAbundance.R): R code used to continue the model runs after TryModel_spAbundance.R.
3. Post processing:
    1. [PostProcess_MsAbund.R](R/PostProcess_MsAbund.R): Code that takes several model runs, puts them together, checks convergence, and does some initial abundance and richness processing. This code is run separately for each group that is used for modeling. These outputs aren't ready for interpretation; needs to be combined across groups in the next step.
4. Results: 
    1. [Abundances.R](R/Abundances.R): Calculates species-specific, group, and overall patterns in abundance. Includes code to make figure 3 (abundance trends). 
    2. [Richness.R](R/Richness.R): Calculates species-specific, group, and overall patterns in richness. 
    3. [Diversity.R](R/Diversity.R): Calculates species-specific, group, and overall patterns in evenness (Shannon diversity). 
5. Synthesis:
    1. [MultiMetricFigures.R](R/MultiMetricFigures.R): Takes output from Abundances.R, Richness.R, and Diversity.R to make composite figures (Figure 2) and tables. 
    2. [Map_Figure1.R](R/Map_Figure1.R): Takes survey data and creates a map of the number of surveys per location (Figure 1). 
6. Batch files for super computer:
    Note: likely specific to MSU's ICER HPCC at the time of running. It might be a helpful starting place but wouldn't be ready to use as-is. 
    1. [run_BCSmodel_msAbund.sb](run_BCSmodel_msAbund.sb): Runs the code in [TryModel_spAbundance.R](R/TryModel_spAbundance.R). Runs one chain; run command 5 times to run 5 chains.
    2. [run_UpdateModel.sb](run_UpdateModel.sb): Runs the code in [UpdateModel_spAbundance.R](R/UpdateModel_spAbundance.R). Updates the specified chain only. 
    3. [run_PostProcess.sb](run_PostProcess.sb): Runs the code in [PostProcess_MsAbund.R](R/PostProcess_MsAbund.R)
    
*Note that additional folders (e.g. 'Output', 'Output/1_runs') would be necessary to save output from the Post-process, Results, and Synthesis steps, as specified in the code.*