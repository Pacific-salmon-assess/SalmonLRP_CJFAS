---
output:
  pdf_document: default
  html_document: default
---

# Projection-based reference points for Pacific salmon

Contacts: 
Carrie Holt, DFO (Carrie.Holt@dfo-mpo.gc.ca)
Kendra Holt, DFO (Kendra.Holt@dfo-mpo.gc.ca). 

Code developed by Kendra Holt and Carrie Holt, DFO. 

### Overview

This readme file provides steps to run projection-based reference points for West coast Vancouver Island (WCVI) Chinook salmon, as implemented for the manuscript "Adapting quantitative tools to support the Policy in an era of legislative and environmental change" for a CJFAS special issue on the Wild Salmon Policy, WSP. The files needed to run these analyses are primarily in this repository "SalmonLRP_CJFAS", but some inputs files are generated in a second repository, "Watershed-Area-Model" repository and are saved here. The "Watershed-Area-Model" repository contains files and code required to estimate benchmarks for WCVI Chinook, and is provided as an additional resource but is not needed to run this code.

The "SalmonLRP_CJFAS" repository is forked from the "SalmonLRP_RetroEval" repository which provides code used to estimate projection-based reference points that are documented in Holt, K. et al. (2023), including WCVI Chinook as a case study.

The analyses and results from this repository are reported in Holt, C. et al. (in review). Citations are listed below.

Code and associated files are organized into the following sub-folders: 


1. **Code**
 * Contains files/functions used to estimate projection LRPs. All files in this folder are intended to be generic code that can be applied to any stock aggregate. This code was developed for Holt, K. et al. (2023) to run for case studies on WCVI Chinook, Interior Fraser Coho, and South Coast Chum salmon, with if statements for individual SMUs were specificity is required.
 
 2. **WCVIChinookStudy**
 * Contains files specific to West Coast Vancouver Island Chinook, which was a case study in Holt, K. et al. (2023), and the focus of Brown et al. (2025). The Rproj files is run from this folder.
 This folder contains subfolders:
 
 **DataIn**: Folder of input data required for analyses described below in steps 3-7
 
 **DataOut**: Folder of output data derived from analyses described below in step 3-7.
 
 **Figures**: Folder of figures generated from analyses described below in steps 3-8.
 
 **RmdReports**: Folder of R Markdown reports generated for Holt, K. et al (2023) case study analysis on WCVI Chinook
 
 **samSimInputs**: Folder of inputs required to projection-based reference points, executed using samSim R package (see Step 6 below).
 
 **samSimOutputs**: Folder of outputs generated from estimation of projection-based reference points (see Step 6 below)



#### To run analyses, follow these steps:

**Step 1)** Open Rproj file in the 'SalmonLRP_WCVI_CK\WCVIChinookStudy' folder.

*File*: 'SalmonLRP_WCVI_CK\WCVIChinookStudy\WCVIChinookStudy.Rproj'

**Step 2)** Open 'runWCVIChinook_projLRP.r' file and run the 'Set-up:Libraries and Source Code' section

*File*: runWCVIChinook_projLRP.r (Section: Set-up:Libraries and Source Code)


**Step 3)** Generate a bubble plot of pairwise correlations in spawner abundances between inlets. 

*File*: CorrPlots.R 

*Inputs*:

'DataIn/Inlet_Sum.csv' for summed spawner abundances to the inlet level

'samSimInputs/CUPars.csv' for names of inlets

*Outputs*: 

'Figures/SpawnerCorrelation.png', a bubble plot of pairwise correlations in spawner abundances

*Required for*: Figure S1b in manuscript. Fig S1d is from  Brown et al. (In review)


**Step 4)** Manually complete input files of inlet-specific parameters for projections. All input parameters used in the manuscript are saved to this repository, but can be revised by the user.

*File*: samSimInput/CUPars.csv

*Inputs*: 

column manUnit = WCVIChinook (stock management unit, SMU)

column stk = (each inlet numbered sequentially)

column model = ricker (default stock-recruitment model assumption)

column minER = 0.05 (minimum exploitation rate assuming minimal unavoidable level of bycatch)

column usERscalar = 1 (scalar applied for additional US harvest, assumed default 1)

column canERscalar = 1 (scalar applied for additional CDN harvest, assumed default 1)

column tauCycAge = tau parameter of the multivariate logistic distribution of proportion of ages in recruitment by inlet (all inlets within a CU are identical), estimated from time-series of ages at maturity

column alpha = ln alpha parameter in Ricker model by inlet, provided from Step 3 above, by inlet (all inlets within a CU are identical)

column beta0 = beta parameter of the Ricker model by inlet, calculated from lnalpha/SREP, where lnalpha is generated in Step 3 and SREP from accessible watershed area model for each inlet (file of accessible watershed-area based benchmarks: 'DataIn/WCVI_SMSY_noEnh_wBC.csv').

column sigma = SD of Ricker residuals by inlet, provided from Step 3 above, by inlet (all inlets within a CU are identical)

column aveGenTime = 4 (average generation time in years)

column ageFirstRec = 2 (age at first recruitment)

column ageMaxRec = 4 (generation time used to set length of the initialization period [=ageMaxRec*10], and is set to 4 here)

columns meanRec2, meanRec3, meanRec4, meanRec5, meanRec6 = mean proportion of ages in recruitment, provided from Step 4 above, by inlet (all inlets are identical)

column Sinit = initial spawner abundances, set to spawners at equilibrium generated from accessible watershed area model, by inlet, predicted values from IWAM.r.

The following columns are not used:  domCycle, cvER, coef1, covarInit, mu_logCovar1, sig_logCovar1, min_logCovar, max_logCovar, larkAlpha, larkBeta0, larkBeta1, larkBeta2, larkBeta3, larkSigma, medianRec, lowQRec, highQRec, meanDBE, sdDBE, medMA, meanForecast, sdForecast 

*Output*: updated samSimInput/CUPars.csv

*Required for*: Generating random draws of Ricker model to be used in projections (Step 5 below) and running projections for projection-based reference points (Step 6 below)

**Step 5)** Generate a series of random draws for Ricker model to be used in the projections. Random draws have been generated and saved to the repository already to allow for comparison in outputs for different scenarios using the same series of random draws. Code provided for reference.

*File*: 'runWCVIChinook_projLRP.R' (Section 5) 

*Inputs*:

'samSimInputs/CUPars.csv' for inlet names, and productivity  (ln alpha) and sigma (SD of Ricker recruitment deviations) by inlet

'DataIn/WCVI_SMSY_noEnh_wBC.csv' for SREP estimates based on accessible watershed area model

*Outputs*:

'SamSimInputs/Ricker_mcmc.csv' file of random draws for Ricker parameters

'Figures/AlphaDensity.png', figure of density in ln alpha values from random draws

'Figures/SREPDensity.png', figure of density of SREP values from random draws

*Required for*: Running projection-based reference points (Step 6 below)

**Step 6)** Run projections to derive projection-based reference points. This code requires R package, samSim (branch LRP), which is sourced in the code below. See samSim repository (https://github.com/Pacific-salmon-assess/samSim/tree/LRP) for code, and Holt, K. et al. (2023) Appendix B for equations and documentation.

*File*: runWCVIChinook_projLRP.R file (See Sections 1-3). Section 3 of this R file sources the generic function, run_scenarioProj() in the Code subfolder, which generates projections.


*Inputs*:

'samSimInputs/CUPars.csv' (see Step 4 above)

'samSimInputs/Ricker_mcmc.csv' (see Step 5 above)

*Outputs*:

'SamSimOutputs/diagnostics/baseER/baseER_baseER_singleTrialFig.pdf' Example trajectory from a single Monte Carlo trial in sub-directory, labeled by scenarioName 

'SamSimOutputs/simData/projLRPDat_baseER.csv' Projection data to estimate projection based reference point, where file is labeled by the scenarioName

'SamSimOutputs/simData/projSpawnDat_baseER.csv' Projected spawner-recruit time-series by CU (or inlet), where file is labeled by the scenarioName

'SamSimOutputs/simData/baseER/' standard outputs from samSim, where the sub-directory is labelled by scenarioName. See repository, 'samSim' (LRP branch) for code, and Holt, K. et al. (2023) Appendix B for equations and documentation.

*Required for*: Generating projection-based reference points and figure of projections


**Step 7)** Plot projection results: binned aggregate abundances against the proportion of Monte Carlo trials in that bin where all component inlets were above their lower benchmark. From this plot, projection-based reference points can be identified at various probability levels. Overlain on this plot is the probability of individual inlets being above their lower benchmark.

*File*: runWCVIChinook_projLRP.r (see Section 6).

*Inputs*: 

'SamSimOutputs/simData/projLRPDat_baseER.csv', output from samSim that shows in which Monte Carlo trial and year all inlets were above their lower benchmarks.

'SamSimOutputs/simData/projCUBenchDat_baseER.csv', output file derived in step 10.

'SamSimInputs/CUPars.csv' for list of inlet names

*Outputs*:

'Figures/ProjectedLRPs/baseER-ProjLRPCurve-ALLp.png' for figure showing projection-based reference point for various probabilities of all inlets being above lower benchmarks, with individual inlet probabilities included.

'DataOut/ProjectedLRPs/ProjectedLRPsbaseER_ALLp.csv' for projection-based reference points at various probabilities

'DataOut/ProjectedLRPs/ProjectedLRP_databaseER_Allp.csv' for underlying binned data used to generate Figure of probabilities of all inlets being above lower benchmarks along gradient in aggregate abundances, above.

**Step 8)** Run sensitivity analyses for projection-based reference points, generated in Step 6. See documentation of sensitivity analyses in Holt, K. et al. (2023).

*File*: runWCVIChinook_projLRP.r (see Section 4)

*Inputs*: as in Step 6, for various alternative assumptions

*Outputs*: as in Step 6, for various alternative assumptions


#### Citations
Brown, N. et al. 2025. West Coast of Vancouver Island Natural-Origin Chinook Salmon (Oncorhynchus tshawytscha) Stock Assessment. CSAS Working Paper2025/nnn.

Holt, K.R., Holt, C.A., Warkentin, L., Wor, C., Davis, B., Arbeider, M., Bokvist, J., Crowley, S., Grant, S., Luedke, W., McHugh, D., Picco, C., and Van Will, P. 2023. Case Study Applications of LRP Estimation Methods to Pacific Salmon Stock Management Units. DFO Can. Sci. Advis. Sec. Res. Doc. 2023/010. iv+129p.


