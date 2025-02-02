# ==============================================================================
# Calculation of Projected Limit Reference Points for CJFAS special issue
# Carrie Holt, Last update: 2 Feb 2025
# ==============================================================================

# Projected LRPs represent the aggregate abundance that is associated with a
# specified probability that a required proportion of CUs will be above their
# lower benchmarks in projections
# --  E.g., the LRP may represent the aggregate abundance that is projected to
# have a 50% probability that all CUs (100% of CUs) will be above their lower
# benchmarks

# Projections are done using the LRP branch of the samSim model
# (https://github.com/Pacific-salmon-assess/samSim)
# The code in this file is divided into the following sections:
#     (1) Read-in WCVI Chinook data
#     (2) Create directories
#     (3) Run base projections (Using 4 different OM models at present)
#     (4) Run sensitivity analysis projections
#     (5) Code to create mcmcOut for Ricker pars
#     (6) Plot LRPs with various plevels

# ==============================================================================

# ==============================================================================
# Set-up:Libraries and Source Code
# ==============================================================================

# A polite helper for installing packages
# (adapted from Hadley Wickham's install script) -------------------------------

please_install <- function(pkgs, install_fun = install.packages) {
  if (length(pkgs) == 0) {
    return(invisible())
  }
  if (!interactive()) {
    stop("Please run in interactive session", call. = FALSE)
  }
  title <- paste0(
    "Ok to install these packges?\n",
    paste("* ", pkgs, collapse = "\n")
  )
  ok <- menu(c("Yes", "No"), title = title) == 1
  if (!ok) {
    return(invisible())
  }
  install_fun(pkgs)
}

# Do you have all the needed packages? ------------------------------------

pkgs <- c(
  "rsample", "tidyverse", "ggplot2", "gridExtra", "reshape2", "TMB", "viridis",
  "here", "zoo", "corrplot", "RColorBrewer"
)
have <- rownames(installed.packages())
needed <- setdiff(pkgs, have)

please_install(needed)


setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
wcviCKDir<-paste(rootDir,"/WCVIChinookStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("ProjLRP_Functions.r")
  source("plotFunctions.r")
  source("helperFunctions.r")
  source("get.mv.logistic.tau.r")
}
sourceAll()

# # Load TMB models for fitting Bayesian Ricker stock recruit models;
# #   outputs from these model fits are used to parameterize samSim
#
# compile("TMB_Files/SR_IndivRicker_Surv_noLRP.cpp")
# dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_noLRP"))
#
# compile("TMB_Files/SR_HierRicker_Surv_noLRP.cpp")
# dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv_noLRP"))
#
# compile("TMB_Files/SR_IndivRicker_SurvCap_noLRP.cpp")
# dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap_noLRP"))
#
# compile("TMB_Files/SR_HierRicker_SurvCap_noLRP.cpp")
# dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap_noLRP"))
#

# ======================================================================
#(1)  Read-in WCVI Chinook data and plot:
# =====================================================================
setwd(wcviCKDir)
remove.EnhStocks <- TRUE#FALSE#TRUE
CoreInd <- FALSE
AllExMH <- FALSE

# Data to estimate correlation matrix
if(!CoreInd & !AllExMH){
  if(remove.EnhStocks) {wcviCKSRDat <- read.csv("DataIn/Inlet_Sum.csv")}
  if(!remove.EnhStocks) {wcviCKSRDat <- read.csv("DataIn/Inlet_Sum_wEnh.csv")}
}
if(CoreInd){wcviCKSRDat <- read.csv("DataIn/Inlet_Sum_CoreInd.csv")}
if(AllExMH){wcviCKSRDat <- read.csv("DataIn/Inlet_Sum_AllExMH.csv")}

# Get this file from Watershed-Area-Model/DataOut

wcviCKSRDat$yr_num <- group_by(wcviCKSRDat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
wcviCKSRDat$CU_ID <- group_by(wcviCKSRDat, Inlet_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing

wcviInletsDF <- wcviCKSRDat %>% group_by(Inlet_Name) %>%
  mutate(genS=genSmooth(Spawners)) %>% ungroup() %>% rename(Year=BroodYear,
                                                            Stock=Inlet_Name,
                                                            SiteEsc=Spawners,
                                                            Value=genS)

#Value = generational geometric smoothed Value
wcviInletsDF <- wcviInletsDF %>% select(Year, Stock, SiteEsc, Value)

wcviCK_inlet <- unique(wcviCKSRDat$Inlet_Name)

# For plotting purposes only here: (these are old files, not used otherwise)
if(remove.EnhStocks) wcviCKbench <- read.csv("DataIn/wcviRPs_noEnh.csv")
if(!remove.EnhStocks) wcviCKbench <- read.csv("DataIn/wcviRPs_wEnh.csv")
# San Juan is still wrong in this file, for enhanced stocks

wcviCKbench <- wcviCKbench %>% filter(Stock %in% wcviCK_inlet) %>%
  select(Stock, SGEN)


statusFn <- function(x, LBM){
    if(is.na(x)) out <- NA
    if(!is.na(x) ) {
      if(x <= LBM) out <- "Red"
      if(x > LBM) out <- "Amber"#Could be green, but not relevant for this calc & plot
    }
  return(out)
}

# Calculate status for generational smoothed spawner abundance
wcviInletsDF <- wcviInletsDF %>% left_join(wcviCKbench, by = "Stock")
Status <- unlist(pmap(list(x=c(wcviInletsDF$Value), LBM=c(wcviInletsDF$SGEN)), statusFn))
wcviInletsDF <- wcviInletsDF %>% add_column(Status=Status)


# Plot timeseries with status
inletPlot <- ggplot(wcviInletsDF) +
  geom_path(aes(y=SiteEsc, x=Year), colour='black', alpha=0.5) +
  geom_point(aes(y=Value, x=Year, colour=Status)) +
  geom_hline(aes(yintercept=SGEN), colour="orange") +
  #geom_hline(aes(yintercept=UBM), colour="black", linetype=2) +
  # ylab("Échappées") +
  # xlab("Année") +
  scale_colour_manual(guide = NULL, breaks = c("None", "Amber", "Green", "Red"), values=c("gray", "gray", "gray","red")) +
  facet_wrap(~interaction(Stock), scales = "free_y") +
  theme_classic()
# ggsave("Figures/chinook-inlet-timeseries.png", plot=inletPlot, width = 6,
#        height = 4, units = "in")
# ggsave("Figures/chinook-inlet-timeseriesFR.png", plot=inletPlot, width = 6,
#        height = 4, units = "in")


# ======================================================================
# (2) Create directories
# =====================================================================


# Create output directories for Projected LRP outputs
figDir <- here(wcviCKDir, "Figures")
if (file.exists(figDir) == FALSE){
  dir.create(figDir)
}
figDir2 <- here(figDir, "ProjectedLRPs")
if (file.exists(figDir2) == FALSE){
  dir.create(figDir2)
}

projOutDir <- here(wcviCKDir, "DataOut")
if (file.exists(projOutDir) == FALSE){
  dir.create(projOutDir)
}
projOutDir2 <- here(projOutDir, "ProjectedLRPs")
if (file.exists(projOutDir2) == FALSE){
  dir.create(projOutDir2)
}

samSimOutDir <- here(wcviCKDir, "samSimOutputs")
if (file.exists(samSimOutDir) == FALSE){
  dir.create(samSimOutDir)
}

# ===================================================================
# (3) Run Base Projections
# ==================================================================

setwd(codeDir)
devtools::install_github("https://github.com/Pacific-salmon-assess/samSim", ref="LRP")
# remotes::install_github('https://github.com/Pacific-salmon-assess/samSim', ref='LRP')


# Create a correlation matrix from spawner time-series, as a proxy for
# correlation in recruitment residuals assuming no density dependence and
# constant harvest. Only used if recruitment time-series are missing
dum <- wcviCKSRDat %>% dplyr::select(CU_ID, BroodYear, Spawners)
dum <- dum %>% tidyr::pivot_wider(names_from=CU_ID,
                           values_from=Spawners) %>% dplyr::select (!BroodYear)
dum <- dum %>% drop_na()
corMat <- cor(dum)


#-------------------------------------------------------------------------------
# Create samSim input files for current scenario
# This takes a long time (~11hours)

scenarioName <- "baseER"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=50000, cvER = 0.085,
                                cvERSMU=0.17, recCorScalar=1, corMat=corMat,
                                agePpnConst=FALSE, annualcvERCU=FALSE,
                                biasCorrectProj=TRUE, ER=0.3)


# ===================================================================
# (4) Run Sensitivity Analyses on Projections
# ==================================================================

scenarioName <- "ER0"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0)


scenarioName <- "ER0.05"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.05)
scenarioName <- "ER0.10"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.10)

scenarioName <- "ER0.15"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.15)

scenarioName <- "ER0.2"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.2)

scenarioName <- "ER0.25"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.25)



scenarioName <- "ER0.35"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.35)

scenarioName <- "ER0.4"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.4)

# scenarioName <- "ER0.45"
#
# projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
#                                 scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,
#                                 TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
#                                 nMCMC=NULL, nProj=10000, cvER = 0.085, cvERSMU=0.17,
#                                 recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
#                                 annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.45)

scenarioName <- "alphaScalar0.75n50000"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=50000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3,
                                alphaScalar=0.75, SREPScalar=1)


scenarioName <- "alphaScalar1.25n50000"

projSpawners <-run_ScenarioProj(SRDat = NULL, BMmodel = NULL,
                                scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,
                                TMB_Inputs=NULL, outDir=wcviCKDir, runMCMC=T,
                                nMCMC=NULL, nProj=50000, cvER = 0.085, cvERSMU=0.17,
                                recCorScalar=1, corMat=corMat, agePpnConst=FALSE,
                                annualcvERCU=FALSE, biasCorrectProj=TRUE, ER=0.3,
                                alphaScalar=1.25, SREPScalar=1)




# ===================================================================
# (5) Code to create mcmcOut for SR parameters for WCVI CK
# ===================================================================


# The random draws are written to a file that is read-in by ProjRLRP_Functions.r
# These are generated and stored so that the same random numbers can be drawn
# with and without log-normal back-transformation bias adjustment

# Inputs:
# Choice of indicators:  base case with or without enhancement,
 # CoreInd (core indicators),
  # or AllEXMH (all indicators except major hatchery facilities)
# SREP: pulled from "WCVI_SMSY_noEnh_wBC.csv". This file originates from the
# github repository, Watershed-Area-Model,
# "Watershed-Area-Model/DataOut/WCVI_SMSY_noEnh_wBC.csv"
# uncertainty in ln(SREP) is assumed to be normal with 95% uncertainty intervals
# at LL and UL prediction intervals, provided by the integrated watershed area
# model, pulled from WCVI_SMSY_AllExMH.csv
# productivity (log alpha) and sigma (sd of Ricker residuals) are from
# "samSimInputs/CUPars.csv" in the SalmonLRP_wCVI_CK repository
# SD in productivity is assumed 0.5 as in the life-cycle model, i.e., relatively
# large uncertainty in productivity (W. Luedke pers. comm.)
# SD in sigma is assumed 0 (i.e., each MCMC draw has the same sigma)


createMCMCout <- TRUE
setwd(wcviCKDir)
alphaScalar <- 1
SREPScalar <- 1
evenPars <- FALSE#TRUE
remove.EnhStocks <- TRUE
CoreInd <- FALSE #Core 6 indicators only
AllExMH <- FALSE#TRUE # all except major hatchery facilities

# Only need to run once to create mcmcOut.csv file with a given assumed
# distribution of alpha and SREP
if(!CoreInd & !AllExMH) Inlet_Names <- read.csv(paste("samSimInputs/CUPars.csv"))$stkName
if(CoreInd) Inlet_Names <- read.csv(paste("samSimInputs/CUPars_CoreInd.csv"))$stkName
if(AllExMH) Inlet_Names <- read.csv(paste("samSimInputs/CUPars_AllExMH.csv"))$stkName

CU_inlet <- data.frame(Inlet_Names=Inlet_Names, CU_Names=NA)
CU_inlet[Inlet_Names=="Barkley",2] <- "WCVI South"#"Southwest_Vancouver_Island"
CU_inlet[Inlet_Names=="Clayoquot",2] <- "WCVI South"#"Southwest_Vancouver_Island"
CU_inlet[Inlet_Names=="Kyuquot",2] <- "WCVI Nootka & Kyuquot"#"Nootka_Kyuquot"
CU_inlet[Inlet_Names=="Nootka/Esperanza",2] <- "WCVI Nootka & Kyuquot"#"Nootka_Kyuquot"
CU_inlet[Inlet_Names=="Quatsino",2] <- "WCVI North"#"Northwest_Vancouver_Island"
CU_inlet[Inlet_Names=="San Juan",2] <- "WCVI South"#"Northwest_Vancouver_Island"

if(!remove.EnhStocks & !CoreInd & !AllExMH){
  Inlet_Names <- read.csv(paste("samSimInputs/CUPars.csv"))$stkName
  Inlet_Names <- c(Inlet_Names,"Nitinat", "San Juan")
  CU_inlet <- data.frame(Inlet_Names=Inlet_Names, CU_Names=NA)
  CU_inlet[Inlet_Names=="Barkley",2] <- "WCVI South"#"Southwest_Vancouver_Island"
  CU_inlet[Inlet_Names=="Clayoquot",2] <- "WCVI South"#"Southwest_Vancouver_Island"
  CU_inlet[Inlet_Names=="Kyuquot",2] <- "WCVI Nootka & Kyuquot"#"Nootka_Kyuquot"
  CU_inlet[Inlet_Names=="Nootka/Esperanza",2] <- "WCVI Nootka & Kyuquot"#"Nootka_Kyuquot"
  CU_inlet[Inlet_Names=="Nitinat",2] <- "WCVI South"#"Northwest_Vancouver_Island"
  CU_inlet[Inlet_Names=="San Juan",2] <- "WCVI South"#"Northwest_Vancouver_Island"
  CU_inlet[Inlet_Names=="Quatsino",2] <- "WCVI North"#"Northwest_Vancouver_Island"
}



if(createMCMCout){
  set.seed(1)
  nTrials <- 50000
  # Set up matrix of random numbers to use for generating alphas, so that
  # the same random numbers are used for Ricka estimates with bias correction
  # and without bias correction when using alpha to estimata beta (lnA/SREP)
  a_rand <- matrix(runif(nTrials*1.5*length(Inlet_Names)), nrow=nTrials*1.5,
                   ncol=length(Inlet_Names))

  # Pull SREP estimates from Watershed-Area model (see repository "Watershed-
  # Area-Model"). That model included a bias correction for back-transformation
  # from log-space
  if (CoreInd){
    SREP <- data.frame(read.csv(
      "DataIn/WCVI_SMSY_CoreInd.csv"))
  }
  if (AllExMH){
    SREP <- data.frame(read.csv(
      "DataIn/WCVI_SMSY_AllExMH.csv"))
  }

  if (!CoreInd & !AllExMH){
    if (remove.EnhStocks) SREP <- data.frame(read.csv(
      "DataIn/WCVI_SMSY_noEnh_wBC.csv"))
    if (!remove.EnhStocks) SREP <- data.frame(read.csv(
      "DataIn/WCVI_SMSY_wEnh_wBC.csv"))
  }



  #Get lnaplha
  if(!CoreInd & !AllExMH){
    # lnalpha from Diana Dobson's run reconstruction coded in R/TMB with bias correction
    lnalpha_inlet <- read.csv("samSimInputs/CUPars.csv") %>%
      dplyr::select(alpha,stkName) %>%
      rename(inlets=stkName, lnalpha=alpha)
    lnalpha_nBC_inlet <- read.csv("samSimInputs/CUPars_nBC.csv") %>%
      dplyr::select(alpha,stkName) %>%
      rename(inlets=stkName, lnalpha_nBC=alpha)
    lnalpha_inlet$lnalpha <- lnalpha_inlet$lnalpha * alphaScalar#1 #(value for life-stage model =1)
    lnalpha_nBC_inlet$lnalpha_nBC <- lnalpha_nBC_inlet$lnalpha_nBC * alphaScalar#1-0.5^2/2#
    if(!remove.EnhStocks){
      # Add two enhanced inlets, assuming same prod at SWVI
      lnalpha_inlet[6:7,1] <- lnalpha_inlet$lnalpha[2]
      lnalpha_inlet[6:7,2]  <- c("Nitinat", "San Juan")
      lnalpha_nBC_inlet[6:7,1] <- lnalpha_nBC_inlet$lnalpha[2]
      lnalpha_nBC_inlet[6:7,2] <- c("Nitinat", "San Juan")
    }
  }

  if(CoreInd){
    lnalpha_inlet <- read.csv("samSimInputs/CUPars_CoreInd.csv") %>%
      dplyr::select(alpha,stkName) %>%
      rename(inlets=stkName, lnalpha=alpha)
    lnalpha_nBC_inlet <- read.csv("samSimInputs/CUPars_nBC_CoreInd.csv") %>%
      dplyr::select(alpha,stkName) %>%
      rename(inlets=stkName, lnalpha_nBC=alpha)
    lnalpha_inlet$lnalpha <- lnalpha_inlet$lnalpha * alphaScalar#1 #(value for life-stage model =1)
    lnalpha_nBC_inlet$lnalpha_nBC <- lnalpha_nBC_inlet$lnalpha_nBC * alphaScalar#1-0.5^2/2#

  }
  if(AllExMH){
    lnalpha_inlet <- read.csv("samSimInputs/CUPars_AllExMH.csv") %>%
      dplyr::select(alpha,stkName) %>%
      rename(inlets=stkName, lnalpha=alpha)
    lnalpha_inlet$lnalpha <- lnalpha_inlet$lnalpha * alphaScalar#1 #(value for life-stage model =1)
    lnalpha_nBC_inlet <- read.csv("samSimInputs/CUPars_nBC_AllExMH.csv") %>%
      dplyr::select(alpha,stkName) %>%
      rename(inlets=stkName, lnalpha_nBC=alpha)
    lnalpha_nBC_inlet$lnalpha_nBC <- lnalpha_nBC_inlet$lnalpha_nBC * alphaScalar#1 #(value for life-stage model =1)

  }

    SREP <- SREP %>% filter(Stock %in% Inlet_Names) %>% filter(Param=="SREP") %>%
      dplyr::select(!c(X, Param)) %>% rename(SREP=Estimate, inlets=Stock)
    SREP <- SREP %>% mutate(SREP=SREP * SREPScalar) %>%
      mutate(LL=LL * SREPScalar) %>%
      mutate(UL=UL * SREPScalar)

    out <- SREP %>% left_join(lnalpha_inlet, by="inlets") %>%
      left_join(lnalpha_nBC_inlet, by="inlets")



  #Draw alpha value, then draw logSREP parameters,then calc beta for that draw
  # (lnalpha/SREP)
  for (i in 1:length(Inlet_Names)){
    meanSREP <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(SREP)
    logmeanSREP <- log(meanSREP)
    ULSREP <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(UL)
    logULSREP <- log(ULSREP)
    LLSREP <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(LL)
    logLLSREP <- log(LLSREP)
    sigSREP <- (logmeanSREP-logLLSREP)/1.96
    #sigSREP <- (logULSREP-logmeanSREP)/1.96 #Check should be same. Yes.
    # Without log-normal bias correction when sampling log-beta
    # rSREP <- exp(rnorm(nTrials*1.5, logmeanSREP,sigSREP))
    # With a log-normal bias correction when sampling log-beta
    rSREP <- exp(rnorm(nTrials*1.5, logmeanSREP - 0.5*sigSREP^2, sigSREP))

    if(!CoreInd & !AllExMH){
      rsig <-  read.csv(paste("samSimInputs/CUPars.csv")) %>%
        filter(stkName==Inlet_Names[i]) %>% dplyr::select(sigma,stk)
      if(!remove.EnhStocks){ # use same rsig for enhanced indicators as for Barkely (in SWVI)
        if(Inlet_Names[i] == "Nitinat" | Inlet_Names[i] == "San Juan") {
          rsig <- read.csv(paste("samSimInputs/CUPars.csv")) %>%
            filter(stkName == "Barkley") %>% dplyr::select(sigma,stk)
          rsig$stk <- i
        }
      }
    }
    if(CoreInd){
      rsig <- read.csv(paste("samSimInputs/CUPars_CoreInd.csv")) %>%
        filter(stkName==Inlet_Names[i]) %>% dplyr::select(sigma,stk)
    }
    if(AllExMH){
      rsig <- read.csv(paste("samSimInputs/CUPars_AllExMH.csv")) %>%
        filter(stkName==Inlet_Names[i]) %>% dplyr::select(sigma,stk)
    }

    meanlnalpha_nBC <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(lnalpha_nBC)#1-(rsig$sigma^2)/2# (life-stage model)
    meanlnalpha <- out %>% filter(inlets==Inlet_Names[i]) %>% pull(lnalpha)#1# (life-stage model)
    # ULlnalpha <- 2
    # LLlnalpha <- 0
    siglnalpha <- 0.5 # Assuming 95% CIs at 0 and 2, sig ~0.5.#0.25 (narrow)


    # Generate random lnalpha values using same random numbers with and withtout
    # bias correction (but diff for each CU or inlet)
    rlnalpha_nBC <- data.frame(a=qnorm(a_rand[,i], meanlnalpha_nBC, siglnalpha))
    rlnalpha <- data.frame(a=qnorm(a_rand[,i], meanlnalpha, siglnalpha))
    amin <- 0#(meanlnalpha - siglnalpha)# (narrow)- not implemented
    amax <- max(2,alphaScalar*2)#(meanlnalpha + siglnalpha)# (narrow)- not implemented



    # Create a dataframe of alpha (with BC), beta (from alpha w/out BC to
    #stabilize beta with and without BC)
    df <- data.frame( stk=rsig$stk, alpha=rlnalpha$a,
                      beta=rlnalpha_nBC$a/rSREP, sigma= rsig$sigma, SREP=rSREP,
                      stkName=Inlet_Names[i], alpha_nBC = rlnalpha_nBC$a )
    #Remove all rows with Ricker a greater or less than bounds
    df <- df %>% filter(alpha > amin & alpha < amax & alpha_nBC > amin &
                          alpha_nBC < amax) %>% slice(1:nTrials)

    if (i==1) mcmcOut <- df
    if (i>1) mcmcOut <- mcmcOut %>% add_row(df)

  }



  if(CoreInd){
    write.csv(mcmcOut, paste(wcviCKDir, "SamSimInputs","Ricker_mcmc_CoreInd.csv", sep="/"),#"Ricker_mcmc_narrow.csv",#_lifeStageModel
                                                   row.names=F)
  }

  if(AllExMH){
      write.csv(mcmcOut, paste(wcviCKDir, "SamSimInputs","Ricker_mcmc_AllExMH.csv", sep="/"),#"Ricker_mcmc_narrow.csv",#_lifeStageModel
                row.names=F)
  }


  if(!CoreInd & !AllExMH){
    if(remove.EnhStocks){
      if(alphaScalar==1&SREPScalar==1) write.csv(mcmcOut, paste(wcviCKDir, "SamSimInputs","Ricker_mcmc.csv", sep="/"),#"Ricker_mcmc_narrow.csv",#_lifeStageModel
                                                 row.names=F)
      # 15 March 2024. Created new MCMC outputs for updated benchmarks for San Juan, with Enh
      # Since San Juan is an ehcaned stock, no need to update LRPs wihtout enhancement

    }

  # if(remove.EnhStocks){
  #   if(alphaScalar==1&SREPScalar==1) write.csv(mcmcOut, paste(wcviCKDir, "SamSimInputs","Ricker_mcmc_wEnh_2024.csv", sep="/"),#"Ricker_mcmc_narrow.csv",#_lifeStageModel
  #                                              row.names=F)
  #   # 15 March 2024. Created new MCMC outputs for updated benchmarks for San Juan, with Enh
  #   # Since San Juan is an ehcaned stock, no need to update LRPs wihtout enhancement
  #
  # }

  if(alphaScalar!=1 | SREPScalar!=1) write.csv(mcmcOut, paste(wcviCKDir, "/SamSimInputs/Ricker_mcmc_alphaScalar",alphaScalar ,"_SREPScalar",SREPScalar,".csv", sep=""),
            row.names=F)
  }
  #plot of alpha and SREP density
  alphaDensity <- mcmcOut %>% ggplot(aes(alpha, colour=factor(stkName))) + #geom_histogram()
    geom_density() +theme(legend.title = element_blank())

  # Histogram might be better to see smoothed curves with 50,000 mcmc
  mcmcOut %>% ggplot(aes(alpha)) + geom_histogram() + facet_wrap(~factor(stkName))

  SREPDensity <- mcmcOut  %>%
    ggplot(aes(SREP, colour=factor(stkName), fill=factor(stkName))) +
    geom_density(alpha=0.1) +theme(legend.title = element_blank()) +xlim(0,30000)

  if(alphaScalar==1 & SREPScalar==1) {
    ggsave(paste(wcviCKDir,"/Figures/AlphaDensity.png",sep=""),#_lifeStageModel#_narrow
           plot = alphaDensity,
           width = 6, height = 4, units = "in")
    ggsave(paste(wcviCKDir,"/Figures/SREPDensity.png",sep=""),
           plot = SREPDensity,#"/Figures/AlphaDensity_narrow.png"
           width = 6, height = 4, units = "in")
  }

  if(alphaScalar!=1)  {
    ggsave(paste(wcviCKDir,"/Figures/AlphaDensity_alphaScalar",alphaScalar,
                 ".png", sep=""),
           plot = alphaDensity,
           width = 6, height = 4, units = "in")
  }
  if(SREPScalar!=1) {
    ggsave(paste(wcviCKDir,"/Figures/SREPDensity_SREPScalar",SREPScalar,
                 ".png", sep=""),
           plot = SREPDensity,
           width = 6, height = 4, units = "in")

  }

}

# sd((mcmcOut %>% filter(stkName=="Quatsino"))$beta)
# 0.0001692311

set.seed(1)
nInlets <- length(Inlet_Names)
rlnalpha_even <- data.frame(a=qnorm(runif(nTrials * nInlets), 1.5, 0.5))
rbeta_even <- data.frame(a=qnorm(runif(nTrials * nInlets), 1/3155, 0.0001692344))
rsigma <- rep(0.6821667,nTrials * nInlets)

mcmc_even <- data.frame(stk = rep(1:5, 1, each=nTrials), alpha_ = rlnalpha_even,
                        beta_ = rbeta_even,
                        sigma = rsigma,
                        SREP_ = 1/rbeta_even,
                        stkName = rep(Inlet_Names, 1, each=nTrials),
                        alpha_nBC_ = rlnalpha_even )
colnames(mcmc_even) <- c("stk", "alpha", "beta", "sigma", "SREP", "stkName", "alpha_nBC")

# write.csv(mcmc_even, paste(wcviCKDir, "SamSimInputs","Even_mcmc.csv", sep="/"),
#           row.names=F)
#
# not sure if I need alpha_nBC for mcmc_even??

#Look at Sgens by inlet
# test <- mcmcOut#
 ##test <- read.table("samSimInputs/Ricker_mcmc_alphaScalar0.75_SREPScalar1.csv")
# SGEN <- unlist(purrr::map2 (.x=test$alpha, .y=test$beta, .f=sGenSolver))
# SGEN <- data.frame(SGEN=SGEN, stkName=test$stkName)
# sGenDensity <- SGEN %>% ggplot(aes(SGEN)) + geom_density() +
#   xlim(0,3000) + facet_wrap(~stkName)
# ggsave(paste(wcviCKDir,"/Figures/sGenDensity.png", sep=""),
#        plot=sGenDensity, width=4, height=3, units="in")



# ===================================================================
# (6) Plot LRPs with various plevels
# ==================================================================

# Specify threshold to use when calculating LRP
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThresh<-c(0.50,0.66)# probability theshhold; the LRP is set as the aggregate abundance that has this
# probability that the propCUThreshold is met
plot.CUs <- FALSE

# Specify scenarios to calculate LRPs and make plots for.
# These scenarios will be looped over below with a LRP (and LRP plot) saved for each scenario
OMsToInclude<-c(
  "baseER")
  ## Uncomment various OMs for different sensitivity analyses, below
  # "ER0",
  # "ER0.05",
  # "ER0.10",
  # "ER0.15",
  # "ER0.2",
  # "ER0.25",
  # "baseER",
  # "ER0.35",
  # "ER0.4")
  # "alphaScalar0.75n50000",
  # "baseER",
  # "alphaScalar1.25n50000")


if(length(OMsToInclude)==1) OMsToIncludeName <- OMsToInclude[1]
if(length(OMsToInclude)==9) OMsToIncludeName <- "ERs"
if(length(OMsToInclude)==3) OMsToIncludeName <- "Alphas"

LRP <- NA

for (OM in 1:length(OMsToInclude)){

  # Loop over OM Scenarios
  for (i in 1:length(probThresh)) {

    # Read in samSim outputs for OM
    filename<-paste("projLRPDat_",OMsToInclude[OM],".csv",sep="")
    filenameCU<-paste("projCUBenchDat_",OMsToInclude[OM],".csv",sep="")
    projLRPDat<-read.csv(here::here(wcviCKDir, "SamSimOutputs", "simData",filename))
    if(plot.CUs) projCUBenchDat<-read.csv(here::here(wcviCKDir, "SamSimOutputs", "simData",filenameCU))

    # NEED TO CHANGE THIS WHEN I CHANGE OM TO INCLUDE
    CUpars <- read.csv(paste(wcviCKDir, "SamSimInputs/CUPars.csv",sep="/"))
    #######################################################

    if(plot.CUs) projCUBenchDat<-projCUBenchDat %>% filter(year > CUpars$ageMaxRec[1]*10)
    projLRPDat<-projLRPDat %>% filter(year > CUpars$ageMaxRec[1]*10)

    # Create bins for projected spawner abundances
    minBreak<-0
    maxBreak<-round(max(projLRPDat$sAg),digits=-2)
    binSize<-200 # Note: bin size is currently set here
    breaks<-seq(minBreak, maxBreak,by=binSize)

    # Set bin labels as the mid-point
    projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)))
    if(plot.CUs) projCUBenchDat$bins<-
      cut(projCUBenchDat$sAg,breaks=breaks,
          labels=as.character(rollmean(breaks,k=2)))

    tmp<-projLRPDat %>% group_by(bins) %>%
      summarise(nSims=(length(ppnCUsLowerBM)))

    if(plot.CUs){
      # Summarize nSims in each bin
        tmpCU<-projCUBenchDat %>% group_by(bins) %>%
          summarise(nSims=(length(X1)))
    }


    # Filter out bins with < 100 nSims
      tmp2<-projLRPDat %>% group_by(bins) %>%
        summarise(nSimsProp1=(length(ppnCUsLowerBM[ppnCUsLowerBM == propCUThresh]))) %>%
        add_column(nSims=tmp$nSims) %>% filter(nSims>=10)

    if(plot.CUs){

      # For CU level probability, note filter for low nSims is below
      if(length(CUpars$stk)==2){
        tmp2CU_x1<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX1=(sum(X1)))
        tmp2CU_x2<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX2=(sum(X2)))
        tmp2CU <- tmp2CU_x1 %>% left_join(tmp2CU_x2, by="bins") %>%
          add_column(nSims=tmpCU$nSims) %>% filter(nSims>=50)
      }
      if(length(CUpars$stk)==5){
        tmp2CU_x1<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX1=(sum(X1)))
        tmp2CU_x2<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX2=(sum(X2)))
        tmp2CU_x3<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX3=(sum(X3)))
        tmp2CU_x4<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX4=(sum(X4)))
        tmp2CU_x5<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX5=(sum(X5)))
        tmp2CU <- tmp2CU_x1 %>% left_join(tmp2CU_x2, by="bins") %>%
          left_join(tmp2CU_x3, by="bins") %>%
          left_join(tmp2CU_x4, by="bins") %>%
          left_join(tmp2CU_x5, by="bins") %>%
          add_column(nSims=tmpCU$nSims) %>% filter(nSims>=50)
      }
      if(length(CUpars$stk)==6){
        tmp2CU_x1<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX1=(sum(X1)))
        tmp2CU_x2<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX2=(sum(X2)))
        tmp2CU_x3<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX3=(sum(X3)))
        tmp2CU_x4<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX4=(sum(X4)))
        tmp2CU_x5<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX5=(sum(X5)))
        tmp2CU_x6<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX6=(sum(X6)))
        tmp2CU <- tmp2CU_x1 %>% left_join(tmp2CU_x2, by="bins") %>%
          left_join(tmp2CU_x3, by="bins") %>%
          left_join(tmp2CU_x4, by="bins") %>%
          left_join(tmp2CU_x5, by="bins") %>%
          left_join(tmp2CU_x6, by="bins") %>%
          add_column(nSims=tmpCU$nSims) %>% filter(nSims>=50)
      }
      if(length(CUpars$stk)==7){
        tmp2CU_x1<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX1=(sum(X1)))
        tmp2CU_x2<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX2=(sum(X2)))
        tmp2CU_x3<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX3=(sum(X3)))
        tmp2CU_x4<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX4=(sum(X4)))
        tmp2CU_x5<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX5=(sum(X5)))
        tmp2CU_x6<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX6=(sum(X6)))
        tmp2CU_x7<-projCUBenchDat %>% group_by(bins) %>% summarise(nSimsPropX7=(sum(X7)))
        tmp2CU <- tmp2CU_x1 %>% left_join(tmp2CU_x2, by="bins") %>%
          left_join(tmp2CU_x3, by="bins") %>%
          left_join(tmp2CU_x4, by="bins") %>%
          left_join(tmp2CU_x5, by="bins") %>%
          left_join(tmp2CU_x6, by="bins") %>%
          left_join(tmp2CU_x7, by="bins") %>%
          add_column(nSims=tmpCU$nSims) %>% filter(nSims>=50)
      }
    }

    # For each bin, calculate probability that required proportion of CUs above benchmark
    projLRPDat<-tmp2 %>% add_column(prob=tmp2$nSimsProp1/tmp2$nSims)

    if(plot.CUs){
        if(length(CUpars$stk)==2){
        projCUBenchDat<-tmp2CU %>%
          add_column(prob1=tmp2CU$nSimsPropX1/tmp2CU$nSims) %>%
          add_column(prob2=tmp2CU$nSimsPropX2/tmp2CU$nSims)
      }
      if(length(CUpars$stk)==5){
        projCUBenchDat<-tmp2CU %>%
          add_column(prob1=tmp2CU$nSimsPropX1/tmp2CU$nSims) %>%
          add_column(prob2=tmp2CU$nSimsPropX2/tmp2CU$nSims) %>%
          add_column(prob3=tmp2CU$nSimsPropX3/tmp2CU$nSims) %>%
          add_column(prob4=tmp2CU$nSimsPropX4/tmp2CU$nSims) %>%
          add_column(prob5=tmp2CU$nSimsPropX5/tmp2CU$nSims)
      }
      if(length(CUpars$stk)==6){
        projCUBenchDat<-tmp2CU %>%
          add_column(prob1=tmp2CU$nSimsPropX1/tmp2CU$nSims) %>%
          add_column(prob2=tmp2CU$nSimsPropX2/tmp2CU$nSims) %>%
          add_column(prob3=tmp2CU$nSimsPropX3/tmp2CU$nSims) %>%
          add_column(prob4=tmp2CU$nSimsPropX4/tmp2CU$nSims) %>%
          add_column(prob5=tmp2CU$nSimsPropX5/tmp2CU$nSims) %>%
          add_column(prob6=tmp2CU$nSimsPropX6/tmp2CU$nSims)
      }
      if(length(CUpars$stk)==7){
        projCUBenchDat<-tmp2CU %>%
          add_column(prob1=tmp2CU$nSimsPropX1/tmp2CU$nSims) %>%
          add_column(prob2=tmp2CU$nSimsPropX2/tmp2CU$nSims) %>%
          add_column(prob3=tmp2CU$nSimsPropX3/tmp2CU$nSims) %>%
          add_column(prob4=tmp2CU$nSimsPropX4/tmp2CU$nSims) %>%
          add_column(prob5=tmp2CU$nSimsPropX5/tmp2CU$nSims) %>%
          add_column(prob6=tmp2CU$nSimsPropX6/tmp2CU$nSims) %>%
          add_column(prob7=tmp2CU$nSimsPropX7/tmp2CU$nSims)
      }

    }

    # For each bin, calculate the difference between the threshold probability and the calculated probability
    tmp3 <- projLRPDat %>% filter(nSims>100)# Remove bins where there are very few nSims among LRP options
    min <- min(abs(probThresh[i]-tmp3$prob))

    projLRPDat$diff<-abs(probThresh[i]-projLRPDat$prob)

    # Save projection summaries used to create plots
    projLRPDat$OM.Name<-OMsToInclude[OM]
    if (i == 1) projLRPDat.plot<-projLRPDat
    if (i > 1) projLRPDat.plot<-rbind(projLRPDat.plot,projLRPDat)

    # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
    #LRP[i]<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))
    LRP[i]<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min]))

    # Create a table of LRP estimates to be saved for each OM model
    if (i ==1) {
      LRP_Ests<-data.frame(OMsToInclude[OM], probThresh[i], propCUThresh, LRP[i], binSize)
      names(LRP_Ests)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
    } else {
      tmp.df<-data.frame(OMsToInclude[OM], probThresh[i], propCUThresh, LRP[i], binSize)
      names(tmp.df)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
      LRP_Ests<-rbind(LRP_Ests,tmp.df)
    }

    if(i==1){# Plot projected LRP abundance relationship ===============================================================
      if(OM==1) {
        if(length(OMsToInclude)==1|length(OMsToInclude)==9) {
          plot.width <- 5
          plot.height <- 4
        }
        if(length(OMsToInclude)==3) {
          plot.width <- 5
          plot.height <- 1.5
        }

        png(paste(wcviCKDir,"/Figures/ProjectedLRPs/", OMsToIncludeName,
                  "-ProjLRPCurve-ALLp.png", sep=""), width=plot.width,
                  # "-ProjLRPCurve-ALLpFR.png", sep=""), width=plot.width,
            height=plot.height,
            units="in", res=300)#500
        if(length(OMsToInclude)==9) layout(matrix(c(1:9), 3, 3, byrow = TRUE))
        if(length(OMsToInclude)==3) layout(matrix(c(1:3), 1, 3, byrow = TRUE))

      }# End of if(OM==1) {


    if(length(OMsToInclude)==1){
      plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19,
           xlim=c(0, max( as.numeric(as.character(projLRPDat$bins)),
                          na.rm=T)*1.0 ),
           ylim=c(0,1),
           cex=0.5, cex.lab=1,#1.5,
           xlab="Aggregate Abundance", ylab="Pr (Inlets > Lower Benchmark)")
           # xlab="Aggregate Abundance", ylab="Pr (All inlets > Lower Benchmark)")
           # xlab="Abondance agrégée", ylab="Prob(tous les inlets) > PRI")
          yaxt <- "s"
      if(plot.CUs){
        points(as.numeric(as.character(projCUBenchDat$bins)),
               projCUBenchDat$prob1, pch=19,
               # col = palette(hcl.colors(7, "viridis"))[1],
               col = RColorBrewer::brewer.pal(n=6, name="Dark2")[1],
               cex=0.1) # or 'Classic Tableau'
        points(as.numeric(as.character(projCUBenchDat$bins)),
               projCUBenchDat$prob2, pch=19,
               # col=palette(hcl.colors(7, "viridis"))[2],
               col = RColorBrewer::brewer.pal(n=6, name="Dark2")[2],
               cex=0.1) # or 'Classic Tableau'
        if(length(CUpars$stk)>=5){
          points(as.numeric(as.character(projCUBenchDat$bins)),
                 projCUBenchDat$prob3, pch=19,
                 # col=palette(hcl.colors(7, "viridis"))[3],
                 col = RColorBrewer::brewer.pal(n=6, name="Dark2")[3],
                 cex=0.1) # or 'Classic Tableau'
          points(as.numeric(as.character(projCUBenchDat$bins)),
                 projCUBenchDat$prob4, pch=19,
                 # col=palette(hcl.colors(7, "viridis"))[4],
                 col = RColorBrewer::brewer.pal(n=6, name="Dark2")[4],
                 cex=0.1) # or 'Classic Tableau'
          points(as.numeric(as.character(projCUBenchDat$bins)),
                 projCUBenchDat$prob5, pch=19,
                 # col=palette(hcl.colors(7, "viridis"))[5],
                 col = RColorBrewer::brewer.pal(n=6, name="Dark2")[5],
                 cex=0.1) # or 'Classic Tableau'
        }
        if(length(CUpars$stk)>=6){
          points(as.numeric(as.character(projCUBenchDat$bins)),
                 projCUBenchDat$prob6, pch=19,
                 # col=palette(hcl.colors(7, "viridis"))[6],
                 col = RColorBrewer::brewer.pal(n=6, name="Dark2")[6],
                 cex=0.1) # or 'Classic Tableau'
        }
        if(length(CUpars$stk)==7){
          points(as.numeric(as.character(projCUBenchDat$bins)),
                 projCUBenchDat$prob7, pch=19,
                 col=palette(hcl.colors(7, "viridis"))[7],
                 cex=0.1) # or 'Classic Tableau'
        }
      loc.leg <- max(as.numeric(as.character(projCUBenchDat$bins)))*0.88
      # legend(x=loc.leg, y=0.4, legend = CUpars$stkName, pch=19, cex=0.5, col=palette(hcl.colors(7, "viridis")))
      legend(x=loc.leg, y=0.4, legend = c(CUpars$stkName, "All inlets"), pch=19,
             cex=0.5, col=c(RColorBrewer::brewer.pal(n=6, name="Dark2"), "black"),
             bty="n")
      }
    }# End of if(length(OMsToInclude)==1){

    if(length(OMsToInclude)==9){
      par(mar=c(2.8,3,0.6,0.6))
      xMax <- 50000
      if(OM %in% c(1,4,7)) yaxt <- "s"
      if(OM %in% c(2,3,5,6,8,9)) yaxt <- "n"
      if(OM<7){
        xaxt <- "n"#par(xaxt="n")
      }
      if(OM>=7){
        xaxt <- "s"#par(xaxt="s")
      }
    }# End of if(length(OMsToInclude)==9){
    if(length(OMsToInclude)==3){
        par(mar=c(2.8,2.5,1.1,1))
        xaxt <- "s"
        xMax <- 70000
        if(OM>1) yaxt <- "n"
        if(OM==1) yaxt <- "s"

      }

      if(length(OMsToInclude)>1){
        plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19,
             xlim=c(0, xMax ),
             ylim=c(0,1),
             cex=0.3, cex.lab=1,#1.5,
             xlab="", ylab="", xaxt=xaxt, yaxt=yaxt)
            aty <- seq(0, 1, 0.2)
            if(OM %in% c(2,3,5,6,8,9)) axis(side =2,  at=aty, labels = FALSE)


      }


      if(length(OMsToInclude)==9){
        if(OM<8){
          at1 <- seq(0, 50000, 10000)
          axis(side =1,  at=at1, labels = FALSE)
        }


        panel.title <- c("0%", "5%", "10%", "15%", "20%", "25%", "30%", "35%",
                         "40%")
                         # "45%")
        mtext(text=panel.title[OM], side=3, line=0, at=5000, cex=0.4)

        LRP_50 <- LRP_Ests$LRP[1]#(read.csv(paste(wcviCKDir,
                  #              "/DataOut/ProjectedLRPs/ProjectedLRPs",
                  #              OMsToInclude[OM], "_ALLp.csv", sep="") )%>%
                  # pull(LRP))[1]

        # LRP_66 <- (read.csv(paste(wcviCKDir,
        #                           "/DataOut/ProjectedLRPs/ProjectedLRPs",
        #                           OMsToInclude[OM], "_ALLp.csv", sep="") )%>%
        #              pull(LRP))[2]
        # text(x=35000, y=0.15, labels=paste("LRP(p=0.5)= ", LRP_50), cex=0.6)
        # if(OM<7)  text(x=35000, y=0.05, labels=paste("LRP(p=0.66)= ", LRP_66), cex=0.6)
        # text(x=35000, y=0.15, labels=paste("PRL(p=0,5)= ", LRP_50), cex=0.6)
        # if(OM<7)  text(x=35000, y=0.05, labels=paste("PRL(p=0,66)= ", LRP_66), cex=0.6)
        #text(x=35000, y=0.05, labels=paste("LRP(p=0.66)= ", LRP_66), cex=0.6)

        if(OM==4) {mtext("Probability of all inlets > lower benchmark", side=2,
        # if(OM==4) {mtext("Probabilité que tous les inlets > PRI", side=2,

                        line=1.8,at=0.5, cex=1) }
        if(OM==8) {mtext("Aggregate Abundance", side=1, line=1.8, at=25000,
        # if(OM==8) {mtext("Abondance agrégée", side=1, line=1.8, at=25000,

                         cex=0.7) }

      }# End of if(length(OMsToInclude)==9){


    if(length(OMsToInclude)==3){

      panel.title <- c("0.75 x productivity", "Current productivity",
                        "1.25 x productivity")
      # panel.title <- c("0.75 x productivité", "Productivité actuelle",
      #                   "1.5 x productivité")
      # panel.title <- c("cv = 0", "cv = 0.085",
      #                  "cv = 0.17")
      mtext(text=panel.title[OM], side=3, line=0.2, at=15000, cex=0.5)

      # LRP_50 <- (read.csv(paste(wcviCKDir,
      #                         "/DataOut/ProjectedLRPs/ProjectedLRPs",
      #                         OMsToInclude[OM], "_ALLp.csv", sep="") )%>%
      #             pull(LRP))[1]
      # LRP_66 <- (read.csv(paste(wcviCKDir,
      #                           "/DataOut/ProjectedLRPs/ProjectedLRPs",
      #                           OMsToInclude[OM], "_ALLp.csv", sep="") )%>%
      #             pull(LRP))[2]
      #
      # text(x=40000, y=0.15, labels=paste("LRP(p=0.5)= ", LRP_50), cex=0.4)#if (OM<3): alpha
      # text(x=40000, y=0.05, labels=paste("LRP(p=0.66)= ", LRP_66), cex=0.4)# if (OM>1): alpha
      # text(x=40000, y=0.15, labels=paste("PRL(p=0,5)= ", LRP_50), cex=0.4)#if (OM<3): alpha
      # text(x=40000, y=0.05, labels=paste("PRL(p=0,66)= ", LRP_66), cex=0.4)# if (OM>1): alpha

      if(OM==1) {
                mtext("Prob(all inlets)>lower benchmark", side=2,
                # mtext("Prob(inlets for core indicators)>lower benchmark", side=2,
                # mtext("Prob(tous les inlets) > PRI", side=2,

                       line=1.8,at=0.4, cex=0.55)
                  yaxt <- "s"}
      if(OM==2) {mtext("Aggregate Abundance", side=1, line=1.8, at=30000,
      # if(OM==2) {mtext("Abondance agrégée", side=1, line=1.8, at=30000,
                          cex=0.7) }
      if(OM>1)  yaxt <- "n"


    }# End of if(length(OMsToInclude)==3){

    }# End of if(i==1){

    if (length(OMsToInclude == 9)) lrp.lwd <- 1
    if (length(OMsToInclude != 9)) lrp.lwd <- 2
    abline(h=probThresh[i], lty=2, lwd=lrp.lwd)
    # if(OMsToInclude[OM]!="alphaScalar1.25n50000") {
      if (i==1) abline(v=LRP[i], col="#E69F00", lwd=lrp.lwd)
      # }# "orange" "#E69F00",
    if(OMsToInclude[OM]!="alphaScalar0.75n50000"&OM < 8) { if (i==2)
      abline(v=LRP[i], col="#56B4E9", lwd=lrp.lwd) }#viridis(4, alpha=0.3)[3] #"adjustcolor("#56B4E9", alpha.f = 0.5)
    # if(OMsToInclude[OM]!="alphaScalar0.75"&OM < 7) { if (i==3)
    #   abline(v=LRP[i], col="#009E73", lwd=lrp.lwd) }#viridis(4, alpha=0.3)[3] #"adjustcolor("#56B4E9", alpha.f = 0.5)
    # if(OMsToInclude[OM]!="alphaScalar0.75"&OM < 7) { if (i==4)
    #   abline(v=LRP[i], col="#D55E00", lwd=lrp.lwd) }#viridis(4, alpha=0.3)[3] #"adjustcolor("#56B4E9", alpha.f = 0.5)


    # abline(h=0.9, lty=2)
    # abline(h=0.99, lty=2)

    # if(i==3) abline(v=LRP[i], lwd=4, col=viridis(4, alpha=0.3)[2] )
    # if(i==4) abline(v=LRP[i], lwd=4, col=viridis(4, alpha=0.2)[1] )

    if(i==length(probThresh)) {
      if(OM==length(OMsToInclude)) {
        dev.off()
      }
    }

    # Option to plot histogram of nSims in each Agg Abundance Bin
    #barplot(height = projLRPDat$nSims,names.arg = projLRPDat$bins)

  }

  # # # Save LRPs for all OM scenarios

  # # UNCOMMENT THIS AFTER FINSIHING FRENCH TRANSLATIION!!
  write.csv(LRP_Ests, paste(projOutDir2, "/ProjectedLRPs",  OMsToInclude[OM],
                            "_ALLp.csv", sep=""), row.names=F)
  # Save LRP projection summaries used for calculating and plotting LRP (Optional)
  write.csv(projLRPDat.plot, paste(projOutDir2, "/ProjectedLRP_data", OMsToInclude[OM],
                                   "_Allp.csv", sep=""), row.names=F)
  # Save CU benchmark projection summaries used for calculating and plotting prob of CUs>LBM
  if(plot.CUs) write.csv(projCUBenchDat, paste(projOutDir2, "/projCUBench_data", OMsToInclude[OM],
                                   "_Allp.csv", sep=""), row.names=F)

}# End of for OM in 1:length(OMsToInclude)

