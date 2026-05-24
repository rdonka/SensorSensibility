#################### SENSOR SENSIBILITY ANALYSIS #########################
# This script contains code used to analyze the acute effects of morphine administration
# across escalating doses on DA activity assayed via in vivo fiber photometry. We compared
# three sensors: VTA DA cell body GCaMP6f, NAcLS dLight1.3b, and NAcLS GRABDA2h. Results 
# (Donka et al, 2026) are published at: 

# This script requires processed outputs of transient events and whole session stream
# quantification from the accompanying MATLAB script "SensorSensibility_SignalProcessing.m".
# Data were prepared using PASTa (v1.1.0), an opensource MATLAB toolbox for the processing
# of fiber photometry data and the detection and quantification of spontaneous transient events. 
# See Donka et al 2025 for analysis details:

#  VARIABLES
# SubjectID = subject number - as numeric
# Sex = string; M or F
# FiberPlacement = string; VTA or NAcLS
# FiberSide = string; L or R
# Virus = string; Sensor used for fiber photometry recordings. GCaMP6f, dLight1.3b, or GRABDA2h
# Rig = string; rig name used for data collection
# Power = numeric; power level used for data collection (30 to 50 uW)
# Phase = string; PREINJ or POSTINJ
# TreatNum = numeric; 
    # 1 = 2.5ml Saline
    # 2 = 2.5mg Morph
    # 3 = 5.0ml Saline
    # 4 = 5.0mg Morph
    # 5 = 7.5ml Saline
    # 6 = 7.5mg Morph
    # 7 = 10.0ml Saline
    # 8 = 10.0mg Morph
# Dose = numeric; 2.5, 5, 7.5, or 10
# Weight = numeric; subject body weight in grams
# InjType = string; S (saline) or M (morphine)
# CondNum = numeric
    # 0 = pre injection
    # 1 = post injection
# maxloc = timestamp of each peak
# maxval = actual value of the peak (absolute, not relative to BL)
# blstartloc = timestamp of pre-peak BL window start
# blendloc = timestamp of pre-peak BL window end
# blloc = timestamp of middle of pre-peak BL window
# blval = actual value of the pre-peak baseline (BL window mean)
# amp = amplitude of the peak in Z score (maxval - blval)
# quantheightval = value at half height
# risestartloc = timestamp of half height pre-peak
# risesamples = number of samples from rise start to peak
# risems = duration from rise start to peak in milliseconds
# riseslope = slope of transient rise from half-height to peak
# fallendloc = timestamp of half height post-peak
# fallsamples = number of samples from peak to fall half height
# fallms = duration from peak to fall half height in milliseconds
# fallslope = slope of transient fall from peak to half-height
# widthsamples = number of samples from rise half height to fall half height
# widthms = duration from rise half height to fall half height in milliseconds
# AUC = trapezoidal AUC of transient
# IEIsamples = samples between previous peak and current peak
# IEIms = milliseconds between previous peak and current peak
# IEIs = seconds between previous peak and current peak
# compoundeventnum = 0 or 1; 1 means the peak falls within the window of another detected peak
# AUCwindow = AUC of transient based on the window from BL end to post-peak (default set to 1500ms post peak)
# Bin_3min = Bin ID, with events binned into 3 minute time bins
# Phase = PREINJ or POSTINJ

# PREPARE ENVIRONMENT -----
rootgithubpath <- c('C:/Users/rmdon/Desktop/GitHub_MyRepositories/') # Path for Rachel's laptop to github repositories
rootboxpath <- c('C:/Users/rmdon/Box/') # Path for Rachel's laptop to Box drive

# Load libraries we will need
source(paste(rootgithubpath,"Analysis-FP/Morphine FP/RFunctions/Functions_LoadPackages.R",sep=''))

# Load additional functions
source(paste(rootgithubpath,"Analysis-FP/Morphine FP/RFunctions/Functions_SetPlotVariables.R",sep=''))
source(paste(rootgithubpath,"Analysis-FP/Morphine FP/RFunctions/Functions_quantifyRise.R",sep=''))
source(paste(rootgithubpath,"Analysis-FP/Morphine FP/RFunctions/Functions_formatTables.R",sep=''))


# Set working directory to where the data is saved
setwd("C:/Users/rmdon/Box/JRoit Lab/Publications/In Process/Morphine Fiber Photometry/Analysis/Processed Data Outputs/")

# Prepare output paths
figurepath <- paste(rootboxpath, 'JRoit Lab/Publications/In Process/Morphine Fiber Photometry/Figures/Figure Panels/',sep='')
analysispath <- paste(rootboxpath, 'JRoit Lab/Publications/In Process/Morphine Fiber Photometry/Analysis/Tables/',sep='')

# TRANSIENTS: PREPARE DATA -----
subjectexc <- c()

## Read in data and prep variables -----
pkdataz_raw <- read_csv('transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold.csv')

# Read in other AUCwindowsizes for supplementary comparisons
pkdataz_AUCwindow1000ms_raw <- read_csv('transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold_AUCwindow1000ms.csv')
pkdataz_AUCwindow2000ms_raw <- read_csv('transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold_AUCwindow2000ms.csv')
pkdataz_AUCwindow3000ms_raw <- read_csv('transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold_AUCwindow3000ms.csv')

pkdataz_AUCwindow1000ms_raw$AUCwindow1000ms <- pkdataz_AUCwindow1000ms_raw$AUCwindow
pkdataz_AUCwindow2000ms_raw$AUCwindow2000ms <- pkdataz_AUCwindow2000ms_raw$AUCwindow
pkdataz_AUCwindow3000ms_raw$AUCwindow3000ms <- pkdataz_AUCwindow3000ms_raw$AUCwindow

pkdataz_AUCwindow1000ms <- pkdataz_AUCwindow1000ms_raw %>% filter(!SubjectID %in% subjectexc) %>% select(SubjectID, Dose, InjType, maxloc, AUCwindow1000ms)
pkdataz_AUCwindow2000ms <- pkdataz_AUCwindow2000ms_raw %>% filter(!SubjectID %in% subjectexc) %>% select(SubjectID, Dose, InjType, maxloc, AUCwindow2000ms)
pkdataz_AUCwindow3000ms <- pkdataz_AUCwindow3000ms_raw %>% filter(!SubjectID %in% subjectexc) %>% select(SubjectID, Dose, InjType, maxloc, AUCwindow3000ms)

# Filter if there are excluded subjects or transients outside the 3 minute bins
pkdataz <- pkdataz_raw %>% filter(!SubjectID %in% subjectexc, !is.na(Bin_3min))

# Add AUCwindow2000ms and AUCwindow3000ms
pkdataz <- pkdataz %>% left_join(pkdataz_AUCwindow1000ms,
                                 by = c("SubjectID", "Dose", "InjType", "maxloc"))

pkdataz <- pkdataz %>% left_join(pkdataz_AUCwindow2000ms,
            by = c("SubjectID", "Dose", "InjType", "maxloc"))
            
pkdataz <- pkdataz %>% left_join(pkdataz_AUCwindow3000ms,
                                 by = c("SubjectID", "Dose", "InjType", "maxloc"))

# Add CondNum and Bin
pkdataz$Bin <- pkdataz$Bin_3min
pkdataz$CondNum <- if_else(pkdataz$Bin <= 4, 0, 1)

# Check variables
unique(pkdataz$SubjectID)
unique(pkdataz$Phase)
unique(pkdataz$Sex)
unique(pkdataz$FiberPlacement)
unique(pkdataz$Virus)
unique(pkdataz$Dose)
unique(pkdataz$TreatNum)
unique(pkdataz$CondNum)
unique(pkdataz$Bin)

# Factor variables
pkdataz$Sex <- factor(pkdataz$Sex, levels = c('F','M'))
pkdataz$Virus <- factor(pkdataz$Virus, levels = c('GCaMP6f','dLight1.3b','GRABDA2h'))
pkdataz$TreatNum <- factor(pkdataz$TreatNum)
pkdataz$CondNum <- factor(pkdataz$CondNum)
pkdataz$InjType <- factor(pkdataz$InjType, levels = c('S','M'))

pkdataz$DoseNum <- pkdataz$Dose
pkdataz$Dose <- factor(pkdataz$Dose, levels = c(2.5,5,7.5,10))

pkdataz$BinF <- factor(pkdataz$Bin)

# If fall is NA, set to 2000
pkdataz$fallms <- if_else(is.na(pkdataz$fallms), 2000, pkdataz$fallms)

## Find subject means -----
### Find subject values by condition number - overall session
subjectmeans <- pkdataz %>%
                  group_by(SubjectID, Sex, FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum) %>%
                  dplyr::summarize(freq = n_distinct(maxloc, na.rm=TRUE), 
                                   amp = mean(amp, na.rm=TRUE),
                                   risems = mean(risems, na.rm=TRUE),
                                   fallms = mean(fallms, na.rm=TRUE),
                                   widthms = mean(widthms, na.rm=TRUE),
                                   AUC = mean(AUC, na.rm=TRUE),
                                   AUCwindow = mean(AUCwindow, na.rm=TRUE),
                                   AUCwindow1000ms = mean(AUCwindow1000ms, na.rm=TRUE),
                                   AUCwindow2000ms = mean(AUCwindow2000ms, na.rm=TRUE),
                                   AUCwindow3000ms = mean(AUCwindow3000ms, na.rm=TRUE),
                                   IEIs = mean(IEIs, na.rm=TRUE),
                                   ncompoundevents = sum(compoundeventnum))
# Add pks per minute
subjectmeans$pkspermin <- if_else(subjectmeans$CondNum==0, subjectmeans$freq/12, subjectmeans$freq/60) # Add pks per minute

### Find subject values by bin - 3 minute bins
subjectmeansbin <- pkdataz %>%
  group_by(SubjectID, Sex, FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, Bin) %>%
  dplyr::summarize(freq = n_distinct(maxloc, na.rm=TRUE), 
                   amp = mean(amp, na.rm=TRUE),
                   risems = mean(risems, na.rm=TRUE),
                   fallms = mean(fallms, na.rm=TRUE),
                   widthms = mean(widthms, na.rm=TRUE),
                   AUC = mean(AUC, na.rm=TRUE),
                   AUCwindow = mean(AUCwindow, na.rm=TRUE),
                   AUCwindow1000ms = mean(AUCwindow1000ms, na.rm=TRUE),
                   AUCwindow2000ms = mean(AUCwindow2000ms, na.rm=TRUE),
                   AUCwindow3000ms = mean(AUCwindow3000ms, na.rm=TRUE),
                   IEIs = mean(IEIs, na.rm=TRUE),
                   ncompoundevents = sum(compoundeventnum)) %>%
  ungroup() %>%
  complete(Bin, nesting(SubjectID, Sex, FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum), fill = list(freq = 0))

subjectmeansbin$CondNum <- factor(if_else(subjectmeansbin$Bin <= 4, 0, 1))
subjectmeansbin$pkspermin <- subjectmeansbin$freq/3 # Add pks per minute
subjectmeansbin$BinF <- factor(subjectmeansbin$Bin)

subjectmeansbin$BinTime <- subjectmeansbin$Bin * 3
subjectmeansbin$BinTimeStart <- (subjectmeansbin$Bin*3)-3
subjectmeansbin$BinTimeEnd <- subjectmeansbin$BinTimeStart + 3
subjectmeansbin$BinTimeInj <- (subjectmeansbin$BinTimeStart-12) + 1.5

### Find subject change values -----
## Find subject percent change from saline control session by condition number
epsilon   <- 1e-3
variables <- c("freq", "pkspermin", "amp", "risems", "fallms", "widthms", "AUC", 
               "AUCwindow", "AUCwindow1000ms", "AUCwindow2000ms", "AUCwindow3000ms",
               "IEIs", "ncompoundevents")
keys <- c("SubjectID","Dose","CondNum")

setDT(subjectmeans)

# Build saline side
S_side <- subjectmeans[InjType == "S", c(keys, variables), with = FALSE]
setnames(S_side, variables, paste0(variables, "_S"))

# Join saline values onto all rows
subjectmeans <- S_side[subjectmeans, on = keys]

# Vectorized calculations only on M rows
for (v in variables) {
  Sname <- paste0(v, "_S")
  cname <- paste0(v, "C")
  pcname <- paste0(v, "PC")
  subjectmeans[InjType == "M", (cname)  := get(v) - get(Sname)]
  subjectmeans[InjType == "M", (pcname) := (get(v) - get(Sname)) / (get(Sname) + epsilon) * 100]
  subjectmeans[, (Sname) := NULL]
}

subjectmeans <- as.data.frame(subjectmeans)


## Find subject percent change from saline control session by condition number and bin - 3 minute bins
epsilon   <- 1e-3
variables <- c("freq", "pkspermin", "amp", "risems", "fallms", "widthms", "AUC", 
               "AUCwindow", "AUCwindow1000ms", "AUCwindow2000ms", "AUCwindow3000ms", "IEIs", "ncompoundevents")
binkeys <- c("SubjectID","Dose","CondNum","Bin")

setDT(subjectmeansbin)

# Build saline side
binS_side <- subjectmeansbin[InjType == "S", c(binkeys, variables), with = FALSE]
setnames(binS_side, variables, paste0(variables, "_S"))

# Join saline values onto all rows
subjectmeansbin <- binS_side[subjectmeansbin, on = binkeys]

# Vectorized calculations only on M rows
for (v in variables) {
  Sname <- paste0(v, "_S")
  cname <- paste0(v, "C")
  pcname <- paste0(v, "PC")
  subjectmeansbin[InjType == "M", (cname)  := get(v) - get(Sname)]
  subjectmeansbin[InjType == "M", (pcname) := (get(v) - get(Sname)) / (get(Sname) + epsilon) * 100]
  subjectmeansbin[, (Sname) := NULL]
}

subjectmeansbin <- as.data.frame(subjectmeansbin)

# Add Bin as numeric and re-scale POST bins
subjectmeansbin$BinNum <- as.numeric(subjectmeansbin$Bin)
subjectmeansbin$BinNumPOST <- subjectmeansbin$BinNum - 4


## Normalize variables -----
# Normalize variables to SNR within each sensor for between sensor comparisons
# SNR: Pooled baseline SD within each sensor (from saline rows)
subjectmeansbin_virus_sd <- subjectmeansbin %>%
  filter(InjType == "S") %>%
  group_by(Virus) %>%
  summarise(
    pkspermin_SD_S_pool = sd(pkspermin, na.rm = TRUE),
    amp_SD_S_pool        = sd(amp, na.rm = TRUE),
    AUCwindow_SD_S_pool  = sd(AUCwindow, na.rm = TRUE),
    risems_SD_S_pool     = sd(risems, na.rm = TRUE),
    fallms_SD_S_pool     = sd(fallms, na.rm = TRUE),
    .groups = "drop"
  )

subjectmeansbin <- subjectmeansbin %>%
  left_join(subjectmeansbin_virus_sd, by = "Virus") %>%
  mutate(
    pksperminSNR_virus = if_else(InjType == "M", pksperminC / (pkspermin_SD_S_pool + 1e-6), NA_real_),
    ampSNR_virus       = if_else(InjType == "M", ampC / (amp_SD_S_pool + 1e-6), NA_real_),
    AUCwindowSNR_virus = if_else(InjType == "M", AUCwindowC / (AUCwindow_SD_S_pool + 1e-6), NA_real_),
    risemsSNR_virus    = if_else(InjType == "M", risemsC / (risems_SD_S_pool + 1e-6), NA_real_),
    fallmsSNR_virus    = if_else(InjType == "M", fallmsC / (fallms_SD_S_pool + 1e-6), NA_real_)
  )

subjectmeansbin_subject_sd <- subjectmeansbin %>%
  filter(InjType == "S") %>%
  group_by(SubjectID, CondNum) %>%
  summarise(
    pkspermin_SD_S_subject = sd(pkspermin, na.rm = TRUE),
    amp_SD_S_subject        = sd(amp, na.rm = TRUE),
    AUCwindow_SD_S_subject  = sd(AUCwindow, na.rm = TRUE),
    risems_SD_S_subject     = sd(risems, na.rm = TRUE),
    fallms_SD_S_subject     = sd(fallms, na.rm = TRUE),
    .groups = "drop"
  )

subjectmeansbin <- subjectmeansbin %>%
  left_join(subjectmeansbin_subject_sd, by = c("SubjectID", "CondNum")) %>%
  mutate(
    pksperminSNR_subject = if_else(InjType == "M", pksperminC / (pkspermin_SD_S_subject + 1e-6), NA_real_),
    ampSNR_subject       = if_else(InjType == "M", ampC / (amp_SD_S_subject + 1e-6), NA_real_),
    AUCwindowSNR_subject = if_else(InjType == "M", AUCwindowC / (AUCwindow_SD_S_subject + 1e-6), NA_real_),
    risemsSNR_subject    = if_else(InjType == "M", risemsC / (risems_SD_S_subject + 1e-6), NA_real_),
    fallmsSNR_subject    = if_else(InjType == "M", fallmsC / (fallms_SD_S_subject + 1e-6), NA_real_)
  )


# subjectmeans
subjectmeans_virus_sd <- subjectmeans %>%
  filter(InjType == "S") %>%
  group_by(Virus) %>%
  summarise(
    pkspermin_SD_S_pool = sd(pkspermin, na.rm = TRUE),
    amp_SD_S_pool        = sd(amp, na.rm = TRUE),
    AUCwindow_SD_S_pool  = sd(AUCwindow, na.rm = TRUE),
    risems_SD_S_pool     = sd(risems, na.rm = TRUE),
    fallms_SD_S_pool     = sd(fallms, na.rm = TRUE),
    .groups = "drop"
  )

subjectmeans <- subjectmeans %>%
  left_join(subjectmeans_virus_sd, by = "Virus") %>%
  mutate(
    pksperminSNR_virus = if_else(InjType == "M", pksperminC / (pkspermin_SD_S_pool + 1e-6), NA_real_),
    ampSNR_virus       = if_else(InjType == "M", ampC / (amp_SD_S_pool + 1e-6), NA_real_),
    AUCwindowSNR_virus = if_else(InjType == "M", AUCwindowC / (AUCwindow_SD_S_pool + 1e-6), NA_real_),
    risemsSNR_virus    = if_else(InjType == "M", risemsC / (risems_SD_S_pool + 1e-6), NA_real_),
    fallmsSNR_virus    = if_else(InjType == "M", fallmsC / (fallms_SD_S_pool + 1e-6), NA_real_)
  )

subjectmeans_subject_sd <- subjectmeans %>%
  filter(InjType == "S") %>%
  group_by(SubjectID, CondNum) %>%
  summarise(
    pkspermin_SD_S_subject = sd(pkspermin, na.rm = TRUE),
    amp_SD_S_subject        = sd(amp, na.rm = TRUE),
    AUCwindow_SD_S_subject  = sd(AUCwindow, na.rm = TRUE),
    risems_SD_S_subject     = sd(risems, na.rm = TRUE),
    fallms_SD_S_subject     = sd(fallms, na.rm = TRUE),
    .groups = "drop"
  )

subjectmeans <- subjectmeans %>%
  left_join(subjectmeans_subject_sd, by = c("SubjectID", "CondNum")) %>%
  mutate(
    pksperminSNR_subject = if_else(InjType == "M", pksperminC / (pkspermin_SD_S_subject + 1e-6), NA_real_),
    ampSNR_subject       = if_else(InjType == "M", ampC / (amp_SD_S_subject + 1e-6), NA_real_),
    AUCwindowSNR_subject = if_else(InjType == "M", AUCwindowC / (AUCwindow_SD_S_subject + 1e-6), NA_real_),
    risemsSNR_subject    = if_else(InjType == "M", risemsC / (risems_SD_S_subject + 1e-6), NA_real_),
    fallmsSNR_subject    = if_else(InjType == "M", fallmsC / (fallms_SD_S_subject + 1e-6), NA_real_)
  )



## Find treatment means -----
## Mean by sensor and treatment number
treatmeans <- subjectmeans %>%
                group_by(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum) %>%
                dplyr::summarize(freq_mean = mean(freq, na.rm=TRUE), freq_sd = sd(freq, na.rm=TRUE), freq_n = length(freq), freq_se = freq_sd/(sqrt(freq_n)),
                                 freqC_mean = mean(freqC, na.rm=TRUE), freqC_sd = sd(freqC, na.rm=TRUE), freqC_n = length(freqC), freqC_se = freqC_sd/(sqrt(freqC_n)),
                                 freqPC_mean = mean(freqPC, na.rm=TRUE), freqPC_sd = sd(freqPC, na.rm=TRUE), freqPC_n = length(freqPC), freqPC_se = freqPC_sd/(sqrt(freqPC_n)),
                                 
                                 pkspermin_mean = mean(pkspermin, na.rm=TRUE), pkspermin_sd = sd(pkspermin, na.rm=TRUE), pkspermin_n = length(pkspermin), pkspermin_se = pkspermin_sd/(sqrt(pkspermin_n)),
                                 pksperminC_mean = mean(pksperminC, na.rm=TRUE), pksperminC_sd = sd(pksperminC, na.rm=TRUE), pksperminC_n = length(pksperminC), pksperminC_se = pksperminC_sd/(sqrt(pksperminC_n)),
                                 pksperminPC_mean = mean(pksperminPC, na.rm=TRUE), pksperminPC_sd = sd(pksperminPC, na.rm=TRUE), pksperminPC_n = length(pksperminPC), pksperminPC_se = pksperminPC_sd/(sqrt(pksperminPC_n)),
                                 
                                 pksperminSNR_virus_mean = mean(pksperminSNR_virus, na.rm=TRUE), pksperminSNR_virus_sd = sd(pksperminSNR_virus, na.rm=TRUE), pksperminSNR_virus_n = length(pksperminSNR_virus), pksperminSNR_virus_se = pksperminSNR_virus_sd/(sqrt(pksperminSNR_virus_n)),
                                 pksperminSNR_subject_mean = mean(pksperminSNR_subject, na.rm=TRUE), pksperminSNR_subject_sd = sd(pksperminSNR_subject, na.rm=TRUE), pksperminSNR_subject_n = length(pksperminSNR_subject), pksperminSNR_subject_se = pksperminSNR_subject_sd/(sqrt(pksperminSNR_subject_n)),
                                 
                                 amp_mean = mean(amp, na.rm=TRUE), amp_sd = sd(amp, na.rm=TRUE), amp_n = length(amp), amp_se = amp_sd/(sqrt(amp_n)),
                                 ampC_mean = mean(ampC, na.rm=TRUE), ampC_sd = sd(ampC, na.rm=TRUE), ampC_n = length(ampC), ampC_se = ampC_sd/(sqrt(ampC_n)),
                                 ampPC_mean = mean(ampPC, na.rm=TRUE), ampPC_sd = sd(ampPC, na.rm=TRUE), ampPC_n = length(ampPC), ampPC_se = ampPC_sd/(sqrt(ampPC_n)),
                                 
                                 ampSNR_virus_mean = mean(ampSNR_virus, na.rm=TRUE), ampSNR_virus_sd = sd(ampSNR_virus, na.rm=TRUE), ampSNR_virus_n = length(ampSNR_virus), ampSNR_virus_se = ampSNR_virus_sd/(sqrt(ampSNR_virus_n)),
                                 ampSNR_subject_mean = mean(ampSNR_subject, na.rm=TRUE), ampSNR_subject_sd = sd(ampSNR_subject, na.rm=TRUE), ampSNR_subject_n = length(ampSNR_subject), ampSNR_subject_se = ampSNR_subject_sd/(sqrt(ampSNR_subject_n)),
                                 
                                 risems_mean = mean(risems, na.rm=TRUE), risems_sd = sd(risems, na.rm=TRUE), risems_n = length(risems), risems_se = risems_sd/(sqrt(risems_n)),
                                 risemsC_mean = mean(risemsC, na.rm=TRUE), risemsC_sd = sd(risemsC, na.rm=TRUE), risemsC_n = length(risemsC), risemsC_se = risemsC_sd/(sqrt(risemsC_n)),
                                 risemsPC_mean = mean(risemsPC, na.rm=TRUE), risemsPC_sd = sd(risemsPC, na.rm=TRUE), risemsPC_n = length(risemsPC), risemsPC_se = risemsPC_sd/(sqrt(risemsPC_n)),

                                 risemsSNR_virus_mean = mean(risemsSNR_virus, na.rm=TRUE), risemsSNR_virus_sd = sd(risemsSNR_virus, na.rm=TRUE), risemsSNR_virus_n = length(risemsSNR_virus), risemsSNR_virus_se = risemsSNR_virus_sd/(sqrt(risemsSNR_virus_n)),
                                 risemsSNR_subject_mean = mean(risemsSNR_subject, na.rm=TRUE), risemsSNR_subject_sd = sd(risemsSNR_subject, na.rm=TRUE), risemsSNR_subject_n = length(risemsSNR_subject), risemsSNR_subject_se = risemsSNR_subject_sd/(sqrt(risemsSNR_subject_n)),
                                 
                                 fallms_mean = mean(fallms, na.rm=TRUE), fallms_sd = sd(fallms, na.rm=TRUE), fallms_n = length(fallms), fallms_se = fallms_sd/(sqrt(fallms_n)),
                                 fallmsC_mean = mean(fallmsC, na.rm=TRUE), fallmsC_sd = sd(fallmsC, na.rm=TRUE), fallmsC_n = length(fallmsC), fallmsC_se = fallmsC_sd/(sqrt(fallmsC_n)),
                                 fallmsPC_mean = mean(fallmsPC, na.rm=TRUE), fallmsPC_sd = sd(fallmsPC, na.rm=TRUE), fallmsPC_n = length(fallmsPC), fallmsPC_se = fallmsPC_sd/(sqrt(fallmsPC_n)),
                                 
                                 fallmsSNR_virus_mean = mean(fallmsSNR_virus, na.rm=TRUE), fallmsSNR_virus_sd = sd(fallmsSNR_virus, na.rm=TRUE), fallmsSNR_virus_n = length(fallmsSNR_virus), fallmsSNR_virus_se = fallmsSNR_virus_sd/(sqrt(fallmsSNR_virus_n)),
                                 fallmsSNR_subject_mean = mean(fallmsSNR_subject, na.rm=TRUE), fallmsSNR_subject_sd = sd(fallmsSNR_subject, na.rm=TRUE), fallmsSNR_subject_n = length(fallmsSNR_subject), fallmsSNR_subject_se = fallmsSNR_subject_sd/(sqrt(fallmsSNR_subject_n)),
                                 
                                 widthms_mean = mean(widthms, na.rm=TRUE), widthms_sd = sd(widthms, na.rm=TRUE), widthms_n = length(widthms), widthms_se = widthms_sd/(sqrt(widthms_n)),
                                 widthmsC_mean = mean(widthmsC, na.rm=TRUE), widthmsC_sd = sd(widthmsC, na.rm=TRUE), widthmsC_n = length(widthmsC), widthmsC_se = widthmsC_sd/(sqrt(widthmsC_n)),
                                 widthmsPC_mean = mean(widthmsPC, na.rm=TRUE), widthmsPC_sd = sd(widthmsPC, na.rm=TRUE), widthmsPC_n = length(widthmsPC), widthmsPC_se = widthmsPC_sd/(sqrt(widthmsPC_n)),
                                 
                                 AUC_mean = mean(AUC, na.rm=TRUE), AUC_sd = sd(AUC, na.rm=TRUE), AUC_n = length(AUC), AUC_se = AUC_sd/(sqrt(AUC_n)),
                                 AUCC_mean = mean(AUCC, na.rm=TRUE), AUCC_sd = sd(AUCC, na.rm=TRUE), AUCC_n = length(AUCC), AUCC_se = AUCC_sd/(sqrt(AUCC_n)),
                                 AUCPC_mean = mean(AUCPC, na.rm=TRUE), AUCPC_sd = sd(AUCPC, na.rm=TRUE), AUCPC_n = length(AUCPC), AUCPC_se = AUCPC_sd/(sqrt(AUCPC_n)),
 
                                 AUCwindow_mean = mean(AUCwindow, na.rm=TRUE), AUCwindow_sd = sd(AUCwindow, na.rm=TRUE), AUCwindow_n = length(AUCwindow), AUCwindow_se = AUCwindow_sd/(sqrt(AUCwindow_n)),
                                 AUCwindowC_mean = mean(AUCwindowC, na.rm=TRUE), AUCwindowC_sd = sd(AUCwindowC, na.rm=TRUE), AUCwindowC_n = length(AUCwindowC), AUCwindowC_se = AUCwindowC_sd/(sqrt(AUCwindowC_n)),
                                 AUCwindowPC_mean = mean(AUCwindowPC, na.rm=TRUE), AUCwindowPC_sd = sd(AUCwindowPC, na.rm=TRUE), AUCwindowPC_n = length(AUCwindowPC), AUCwindowPC_se = AUCwindowPC_sd/(sqrt(AUCwindowPC_n)),
                                 
                                 AUCwindowSNR_virus_mean = mean(AUCwindowSNR_virus, na.rm=TRUE), AUCwindowSNR_virus_sd = sd(AUCwindowSNR_virus, na.rm=TRUE), AUCwindowSNR_virus_n = length(AUCwindowSNR_virus), AUCwindowSNR_virus_se = AUCwindowSNR_virus_sd/(sqrt(AUCwindowSNR_virus_n)),
                                 AUCwindowSNR_subject_mean = mean(AUCwindowSNR_subject, na.rm=TRUE), AUCwindowSNR_subject_sd = sd(AUCwindowSNR_subject, na.rm=TRUE), AUCwindowSNR_subject_n = length(AUCwindowSNR_subject), AUCwindowSNR_subject_se = AUCwindowSNR_subject_sd/(sqrt(AUCwindowSNR_subject_n)),
                                 
                                 AUCwindow1000ms_mean = mean(AUCwindow1000ms, na.rm=TRUE), AUCwindow1000ms_sd = sd(AUCwindow1000ms, na.rm=TRUE), AUCwindow1000ms_n = length(AUCwindow1000ms), AUCwindow1000ms_se = AUCwindow1000ms_sd/(sqrt(AUCwindow1000ms_n)),
                                 AUCwindow1000msC_mean = mean(AUCwindow1000msC, na.rm=TRUE), AUCwindow1000msC_sd = sd(AUCwindow1000msC, na.rm=TRUE), AUCwindow1000msC_n = length(AUCwindow1000msC), AUCwindow1000msC_se = AUCwindow1000msC_sd/(sqrt(AUCwindow1000msC_n)),
                                 AUCwindow1000msPC_mean = mean(AUCwindow1000msPC, na.rm=TRUE), AUCwindow1000msPC_sd = sd(AUCwindow1000msPC, na.rm=TRUE), AUCwindow1000msPC_n = length(AUCwindow1000msPC), AUCwindow1000msPC_se = AUCwindow1000msPC_sd/(sqrt(AUCwindow1000msPC_n)),
                                 
                                 AUCwindow2000ms_mean = mean(AUCwindow2000ms, na.rm=TRUE), AUCwindow2000ms_sd = sd(AUCwindow2000ms, na.rm=TRUE), AUCwindow2000ms_n = length(AUCwindow2000ms), AUCwindow2000ms_se = AUCwindow2000ms_sd/(sqrt(AUCwindow2000ms_n)),
                                 AUCwindow2000msC_mean = mean(AUCwindow2000msC, na.rm=TRUE), AUCwindow2000msC_sd = sd(AUCwindow2000msC, na.rm=TRUE), AUCwindow2000msC_n = length(AUCwindow2000msC), AUCwindow2000msC_se = AUCwindow2000msC_sd/(sqrt(AUCwindow2000msC_n)),
                                 AUCwindow2000msPC_mean = mean(AUCwindow2000msPC, na.rm=TRUE), AUCwindow2000msPC_sd = sd(AUCwindow2000msPC, na.rm=TRUE), AUCwindow2000msPC_n = length(AUCwindow2000msPC), AUCwindow2000msPC_se = AUCwindow2000msPC_sd/(sqrt(AUCwindow2000msPC_n)),
                                 
                                 AUCwindow3000ms_mean = mean(AUCwindow3000ms, na.rm=TRUE), AUCwindow3000ms_sd = sd(AUCwindow3000ms, na.rm=TRUE), AUCwindow3000ms_n = length(AUCwindow3000ms), AUCwindow3000ms_se = AUCwindow3000ms_sd/(sqrt(AUCwindow3000ms_n)),
                                 AUCwindow3000msC_mean = mean(AUCwindow3000msC, na.rm=TRUE), AUCwindow3000msC_sd = sd(AUCwindow3000msC, na.rm=TRUE), AUCwindow3000msC_n = length(AUCwindow3000msC), AUCwindow3000msC_se = AUCwindow3000msC_sd/(sqrt(AUCwindow3000msC_n)),
                                 AUCwindow3000msPC_mean = mean(AUCwindow3000msPC, na.rm=TRUE), AUCwindow3000msPC_sd = sd(AUCwindow3000msPC, na.rm=TRUE), AUCwindow3000msPC_n = length(AUCwindow3000msPC), AUCwindow3000msPC_se = AUCwindow3000msPC_sd/(sqrt(AUCwindow3000msPC_n)),
                                 
                                 IEIs_mean = mean(IEIs, na.rm=TRUE), IEIs_sd = sd(IEIs, na.rm=TRUE), IEIs_n = length(IEIs), IEIs_se = IEIs_sd/(sqrt(IEIs_n)),
                                 IEIsC_mean = mean(IEIsC, na.rm=TRUE), IEIsC_sd = sd(IEIsC, na.rm=TRUE), IEIsC_n = length(IEIsC), IEIsC_se = IEIsC_sd/(sqrt(IEIsC_n)),
                                 IEIsPC_mean = mean(IEIsPC, na.rm=TRUE), IEIsPC_sd = sd(IEIsPC, na.rm=TRUE), IEIsPC_n = length(IEIsPC), IEIsPC_se = IEIsPC_sd/(sqrt(IEIsPC_n)),
                                 
                                 ncompoundevents_mean = mean(ncompoundevents, na.rm=TRUE), ncompoundevents_sd = sd(ncompoundevents, na.rm=TRUE), ncompoundevents_n = length(ncompoundevents), ncompoundevents_se = ncompoundevents_sd/(sqrt(ncompoundevents_n)),
                                 ncompoundeventsC_mean = mean(ncompoundeventsC, na.rm=TRUE), ncompoundeventsC_sd = sd(ncompoundeventsC, na.rm=TRUE), ncompoundeventsC_n = length(ncompoundeventsC), ncompoundeventsC_se = ncompoundeventsC_sd/(sqrt(ncompoundeventsC_n)),
                                 ncompoundeventsPC_mean = mean(ncompoundeventsPC, na.rm=TRUE), ncompoundeventsPC_sd = sd(ncompoundeventsPC, na.rm=TRUE), ncompoundeventsPC_n = length(ncompoundeventsPC), ncompoundeventsPC_se = ncompoundeventsPC_sd/(sqrt(ncompoundeventsPC_n)))

treatmeans


## Find bin means -----
### 3 minute bins
treatmeansbin <- subjectmeansbin %>%
  group_by(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, Bin, BinTimeInj, CondNum) %>%
  dplyr::summarize(freq_mean = mean(freq, na.rm=TRUE), freq_sd = sd(freq, na.rm=TRUE), freq_n = length(freq), freq_se = freq_sd/(sqrt(freq_n)),
                   freqC_mean = mean(freqC, na.rm=TRUE), freqC_sd = sd(freqC, na.rm=TRUE), freqC_n = length(freqC), freqC_se = freqC_sd/(sqrt(freqC_n)),
                   freqPC_mean = mean(freqPC, na.rm=TRUE), freqPC_sd = sd(freqPC, na.rm=TRUE), freqPC_n = length(freqPC), freqPC_se = freqPC_sd/(sqrt(freqPC_n)),
                   
                   pkspermin_mean = mean(pkspermin, na.rm=TRUE), pkspermin_sd = sd(pkspermin, na.rm=TRUE), pkspermin_n = length(pkspermin), pkspermin_se = pkspermin_sd/(sqrt(pkspermin_n)),
                   pksperminC_mean = mean(pksperminC, na.rm=TRUE), pksperminC_sd = sd(pksperminC, na.rm=TRUE), pksperminC_n = length(pksperminC), pksperminC_se = pksperminC_sd/(sqrt(pksperminC_n)),
                   pksperminPC_mean = mean(pksperminPC, na.rm=TRUE), pksperminPC_sd = sd(pksperminPC, na.rm=TRUE), pksperminPC_n = length(pksperminPC), pksperminPC_se = pksperminPC_sd/(sqrt(pksperminPC_n)),
                   
                   pksperminSNR_virus_mean = mean(pksperminSNR_virus, na.rm=TRUE), pksperminSNR_virus_sd = sd(pksperminSNR_virus, na.rm=TRUE), pksperminSNR_virus_n = length(pksperminSNR_virus), pksperminSNR_virus_se = pksperminSNR_virus_sd/(sqrt(pksperminSNR_virus_n)),
                   pksperminSNR_subject_mean = mean(pksperminSNR_subject, na.rm=TRUE), pksperminSNR_subject_sd = sd(pksperminSNR_subject, na.rm=TRUE), pksperminSNR_subject_n = length(pksperminSNR_subject), pksperminSNR_subject_se = pksperminSNR_subject_sd/(sqrt(pksperminSNR_subject_n)),
                   
                   amp_mean = mean(amp, na.rm=TRUE), amp_sd = sd(amp, na.rm=TRUE), amp_n = length(amp), amp_se = amp_sd/(sqrt(amp_n)),
                   ampC_mean = mean(ampC, na.rm=TRUE), ampC_sd = sd(ampC, na.rm=TRUE), ampC_n = length(ampC), ampC_se = ampC_sd/(sqrt(ampC_n)),
                   ampPC_mean = mean(ampPC, na.rm=TRUE), ampPC_sd = sd(ampPC, na.rm=TRUE), ampPC_n = length(ampPC), ampPC_se = ampPC_sd/(sqrt(ampPC_n)),
                   
                   ampSNR_virus_mean = mean(ampSNR_virus, na.rm=TRUE), ampSNR_virus_sd = sd(ampSNR_virus, na.rm=TRUE), ampSNR_virus_n = length(ampSNR_virus), ampSNR_virus_se = ampSNR_virus_sd/(sqrt(ampSNR_virus_n)),
                   ampSNR_subject_mean = mean(ampSNR_subject, na.rm=TRUE), ampSNR_subject_sd = sd(ampSNR_subject, na.rm=TRUE), ampSNR_subject_n = length(ampSNR_subject), ampSNR_subject_se = ampSNR_subject_sd/(sqrt(ampSNR_subject_n)),
                   
                   risems_mean = mean(risems, na.rm=TRUE), risems_sd = sd(risems, na.rm=TRUE), risems_n = length(risems), risems_se = risems_sd/(sqrt(risems_n)),
                   risemsC_mean = mean(risemsC, na.rm=TRUE), risemsC_sd = sd(risemsC, na.rm=TRUE), risemsC_n = length(risemsC), risemsC_se = risemsC_sd/(sqrt(risemsC_n)),
                   risemsPC_mean = mean(risemsPC, na.rm=TRUE), risemsPC_sd = sd(risemsPC, na.rm=TRUE), risemsPC_n = length(risemsPC), risemsPC_se = risemsPC_sd/(sqrt(risemsPC_n)),
                   
                   risemsSNR_virus_mean = mean(risemsSNR_virus, na.rm=TRUE), risemsSNR_virus_sd = sd(risemsSNR_virus, na.rm=TRUE), risemsSNR_virus_n = length(risemsSNR_virus), risemsSNR_virus_se = risemsSNR_virus_sd/(sqrt(risemsSNR_virus_n)),
                   risemsSNR_subject_mean = mean(risemsSNR_subject, na.rm=TRUE), risemsSNR_subject_sd = sd(risemsSNR_subject, na.rm=TRUE), risemsSNR_subject_n = length(risemsSNR_subject), risemsSNR_subject_se = risemsSNR_subject_sd/(sqrt(risemsSNR_subject_n)),
                   
                   fallms_mean = mean(fallms, na.rm=TRUE), fallms_sd = sd(fallms, na.rm=TRUE), fallms_n = length(fallms), fallms_se = fallms_sd/(sqrt(fallms_n)),
                   fallmsC_mean = mean(fallmsC, na.rm=TRUE), fallmsC_sd = sd(fallmsC, na.rm=TRUE), fallmsC_n = length(fallmsC), fallmsC_se = fallmsC_sd/(sqrt(fallmsC_n)),
                   fallmsPC_mean = mean(fallmsPC, na.rm=TRUE), fallmsPC_sd = sd(fallmsPC, na.rm=TRUE), fallmsPC_n = length(fallmsPC), fallmsPC_se = fallmsPC_sd/(sqrt(fallmsPC_n)),
                   
                   fallmsSNR_virus_mean = mean(fallmsSNR_virus, na.rm=TRUE), fallmsSNR_virus_sd = sd(fallmsSNR_virus, na.rm=TRUE), fallmsSNR_virus_n = length(fallmsSNR_virus), fallmsSNR_virus_se = fallmsSNR_virus_sd/(sqrt(fallmsSNR_virus_n)),
                   fallmsSNR_subject_mean = mean(fallmsSNR_subject, na.rm=TRUE), fallmsSNR_subject_sd = sd(fallmsSNR_subject, na.rm=TRUE), fallmsSNR_subject_n = length(fallmsSNR_subject), fallmsSNR_subject_se = fallmsSNR_subject_sd/(sqrt(fallmsSNR_subject_n)),
                   
                   widthms_mean = mean(widthms, na.rm=TRUE), widthms_sd = sd(widthms, na.rm=TRUE), widthms_n = length(widthms), widthms_se = widthms_sd/(sqrt(widthms_n)),
                   widthmsC_mean = mean(widthmsC, na.rm=TRUE), widthmsC_sd = sd(widthmsC, na.rm=TRUE), widthmsC_n = length(widthmsC), widthmsC_se = widthmsC_sd/(sqrt(widthmsC_n)),
                   widthmsPC_mean = mean(widthmsPC, na.rm=TRUE), widthmsPC_sd = sd(widthmsPC, na.rm=TRUE), widthmsPC_n = length(widthmsPC), widthmsPC_se = widthmsPC_sd/(sqrt(widthmsPC_n)),
                   
                   AUC_mean = mean(AUC, na.rm=TRUE), AUC_sd = sd(AUC, na.rm=TRUE), AUC_n = length(AUC), AUC_se = AUC_sd/(sqrt(AUC_n)),
                   AUCC_mean = mean(AUCC, na.rm=TRUE), AUCC_sd = sd(AUCC, na.rm=TRUE), AUCC_n = length(AUCC), AUCC_se = AUCC_sd/(sqrt(AUCC_n)),
                   AUCPC_mean = mean(AUCPC, na.rm=TRUE), AUCPC_sd = sd(AUCPC, na.rm=TRUE), AUCPC_n = length(AUCPC), AUCPC_se = AUCPC_sd/(sqrt(AUCPC_n)),
                   
                   AUCwindow_mean = mean(AUCwindow, na.rm=TRUE), AUCwindow_sd = sd(AUCwindow, na.rm=TRUE), AUCwindow_n = length(AUCwindow), AUCwindow_se = AUCwindow_sd/(sqrt(AUCwindow_n)),
                   AUCwindowC_mean = mean(AUCwindowC, na.rm=TRUE), AUCwindowC_sd = sd(AUCwindowC, na.rm=TRUE), AUCwindowC_n = length(AUCwindowC), AUCwindowC_se = AUCwindowC_sd/(sqrt(AUCwindowC_n)),
                   AUCwindowPC_mean = mean(AUCwindowPC, na.rm=TRUE), AUCwindowPC_sd = sd(AUCwindowPC, na.rm=TRUE), AUCwindowPC_n = length(AUCwindowPC), AUCwindowPC_se = AUCwindowPC_sd/(sqrt(AUCwindowPC_n)),
                   
                   AUCwindowSNR_virus_mean = mean(AUCwindowSNR_virus, na.rm=TRUE), AUCwindowSNR_virus_sd = sd(AUCwindowSNR_virus, na.rm=TRUE), AUCwindowSNR_virus_n = length(AUCwindowSNR_virus), AUCwindowSNR_virus_se = AUCwindowSNR_virus_sd/(sqrt(AUCwindowSNR_virus_n)),
                   AUCwindowSNR_subject_mean = mean(AUCwindowSNR_subject, na.rm=TRUE), AUCwindowSNR_subject_sd = sd(AUCwindowSNR_subject, na.rm=TRUE), AUCwindowSNR_subject_n = length(AUCwindowSNR_subject), AUCwindowSNR_subject_se = AUCwindowSNR_subject_sd/(sqrt(AUCwindowSNR_subject_n)),
                   
                   AUCwindow1000ms_mean = mean(AUCwindow1000ms, na.rm=TRUE), AUCwindow1000ms_sd = sd(AUCwindow1000ms, na.rm=TRUE), AUCwindow1000ms_n = length(AUCwindow1000ms), AUCwindow1000ms_se = AUCwindow1000ms_sd/(sqrt(AUCwindow1000ms_n)),
                   AUCwindow1000msC_mean = mean(AUCwindow1000msC, na.rm=TRUE), AUCwindow1000msC_sd = sd(AUCwindow1000msC, na.rm=TRUE), AUCwindow1000msC_n = length(AUCwindow1000msC), AUCwindow1000msC_se = AUCwindow1000msC_sd/(sqrt(AUCwindow1000msC_n)),
                   AUCwindow1000msPC_mean = mean(AUCwindow1000msPC, na.rm=TRUE), AUCwindow1000msPC_sd = sd(AUCwindow1000msPC, na.rm=TRUE), AUCwindow1000msPC_n = length(AUCwindow1000msPC), AUCwindow1000msPC_se = AUCwindow1000msPC_sd/(sqrt(AUCwindow1000msPC_n)),
                   
                   AUCwindow2000ms_mean = mean(AUCwindow2000ms, na.rm=TRUE), AUCwindow2000ms_sd = sd(AUCwindow2000ms, na.rm=TRUE), AUCwindow2000ms_n = length(AUCwindow2000ms), AUCwindow2000ms_se = AUCwindow2000ms_sd/(sqrt(AUCwindow2000ms_n)),
                   AUCwindow2000msC_mean = mean(AUCwindow2000msC, na.rm=TRUE), AUCwindow2000msC_sd = sd(AUCwindow2000msC, na.rm=TRUE), AUCwindow2000msC_n = length(AUCwindow2000msC), AUCwindow2000msC_se = AUCwindow2000msC_sd/(sqrt(AUCwindow2000msC_n)),
                   AUCwindow2000msPC_mean = mean(AUCwindow2000msPC, na.rm=TRUE), AUCwindow2000msPC_sd = sd(AUCwindow2000msPC, na.rm=TRUE), AUCwindow2000msPC_n = length(AUCwindow2000msPC), AUCwindow2000msPC_se = AUCwindow2000msPC_sd/(sqrt(AUCwindow2000msPC_n)),
                   
                   AUCwindow3000ms_mean = mean(AUCwindow3000ms, na.rm=TRUE), AUCwindow3000ms_sd = sd(AUCwindow3000ms, na.rm=TRUE), AUCwindow3000ms_n = length(AUCwindow3000ms), AUCwindow3000ms_se = AUCwindow3000ms_sd/(sqrt(AUCwindow3000ms_n)),
                   AUCwindow3000msC_mean = mean(AUCwindow3000msC, na.rm=TRUE), AUCwindow3000msC_sd = sd(AUCwindow3000msC, na.rm=TRUE), AUCwindow3000msC_n = length(AUCwindow3000msC), AUCwindow3000msC_se = AUCwindow3000msC_sd/(sqrt(AUCwindow3000msC_n)),
                   AUCwindow3000msPC_mean = mean(AUCwindow3000msPC, na.rm=TRUE), AUCwindow3000msPC_sd = sd(AUCwindow3000msPC, na.rm=TRUE), AUCwindow3000msPC_n = length(AUCwindow3000msPC), AUCwindow3000msPC_se = AUCwindow3000msPC_sd/(sqrt(AUCwindow3000msPC_n)),
                   
                   IEIs_mean = mean(IEIs, na.rm=TRUE), IEIs_sd = sd(IEIs, na.rm=TRUE), IEIs_n = length(IEIs), IEIs_se = IEIs_sd/(sqrt(IEIs_n)),
                   IEIsC_mean = mean(IEIsC, na.rm=TRUE), IEIsC_sd = sd(IEIsC, na.rm=TRUE), IEIsC_n = length(IEIsC), IEIsC_se = IEIsC_sd/(sqrt(IEIsC_n)),
                   IEIsPC_mean = mean(IEIsPC, na.rm=TRUE), IEIsPC_sd = sd(IEIsPC, na.rm=TRUE), IEIsPC_n = length(IEIsPC), IEIsPC_se = IEIsPC_sd/(sqrt(IEIsPC_n)),
                   
                   ncompoundevents_mean = mean(ncompoundevents, na.rm=TRUE), ncompoundevents_sd = sd(ncompoundevents, na.rm=TRUE), ncompoundevents_n = length(ncompoundevents), ncompoundevents_se = ncompoundevents_sd/(sqrt(ncompoundevents_n)),
                   ncompoundeventsC_mean = mean(ncompoundeventsC, na.rm=TRUE), ncompoundeventsC_sd = sd(ncompoundeventsC, na.rm=TRUE), ncompoundeventsC_n = length(ncompoundeventsC), ncompoundeventsC_se = ncompoundeventsC_sd/(sqrt(ncompoundeventsC_n)),
                   ncompoundeventsPC_mean = mean(ncompoundeventsPC, na.rm=TRUE), ncompoundeventsPC_sd = sd(ncompoundeventsPC, na.rm=TRUE), ncompoundeventsPC_n = length(ncompoundeventsPC), ncompoundeventsPC_se = ncompoundeventsPC_sd/(sqrt(ncompoundeventsPC_n)))

treatmeansbin



## Scale amp for plots -----
threshold <- 2.6

subjectmeans$ampscaled <- subjectmeans$amp - threshold
subjectmeansbin$ampscaled <- subjectmeansbin$amp - threshold
treatmeans$ampscaled_mean <- treatmeans$amp_mean - threshold
treatmeansbin$ampscaled_mean <- treatmeansbin$amp_mean - threshold


## Prepare AUC Window comparison data -----
treatmeans_AUCcomparison <- treatmeans %>% select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum,
                                                        AUCwindow_mean, AUCwindow1000ms_mean, AUCwindow2000ms_mean, AUCwindow3000ms_mean) %>%
  pivot_longer(cols = c(AUCwindow_mean, AUCwindow1000ms_mean, AUCwindow2000ms_mean, AUCwindow3000ms_mean), names_to = 'Window', values_to = 'AUC')

treatmeans_AUCcomparison$WindowSize <- if_else(treatmeans_AUCcomparison$Window == 'AUCwindow_mean', 1500,
                                                  if_else(treatmeans_AUCcomparison$Window == 'AUCwindow1000ms_mean', 1000,
                                                          if_else(treatmeans_AUCcomparison$Window == 'AUCwindow2000ms_mean', 2000,
                                                                  if_else(treatmeans_AUCcomparison$Window == 'AUCwindow3000ms_mean', 3000, NA))))

treatmeans_AUCcomparison$WindowSizeF <- factor(treatmeans_AUCcomparison$WindowSize)

treatmeans_AUCcomparison$VirusInjTypeDose <- paste(treatmeans_AUCcomparison$Virus,treatmeans_AUCcomparison$InjType, treatmeans_AUCcomparison$Dose, sep='_')

## Find n sizes -----
### VTA GCaMP6f
subjectmeans_GCaMP6f <- subjectmeans %>% filter(Virus=='GCaMP6f')
n_GCaMP6f_2.5_S <- subjectmeans_GCaMP6f %>% filter(Dose==2.5,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_GCaMP6f_2.5_M <- subjectmeans_GCaMP6f %>% filter(Dose==2.5,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_GCaMP6f_5_S <- subjectmeans_GCaMP6f %>% filter(Dose==5,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_GCaMP6f_5_M <- subjectmeans_GCaMP6f %>% filter(Dose==5,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_GCaMP6f_7.5_S <- subjectmeans_GCaMP6f %>% filter(Dose==7.5,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_GCaMP6f_7.5_M <- subjectmeans_GCaMP6f %>% filter(Dose==7.5,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_GCaMP6f_10_S <- subjectmeans_GCaMP6f %>% filter(Dose==10,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_GCaMP6f_10_M <- subjectmeans_GCaMP6f %>% filter(Dose==10,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_GCaMP6f_2.5_S
n_GCaMP6f_2.5_M

n_GCaMP6f_5_S
n_GCaMP6f_5_M

n_GCaMP6f_7.5_S
n_GCaMP6f_7.5_M

n_GCaMP6f_10_S
n_GCaMP6f_10_M

### NAcLS dLight1.3b
subjectmeans_dLight1.3b <- subjectmeans %>% filter(Virus=='dLight1.3b')
n_dLight1.3b_2.5_S <- subjectmeans_dLight1.3b %>% filter(Dose==2.5,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_dLight1.3b_2.5_M <- subjectmeans_dLight1.3b %>% filter(Dose==2.5,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_dLight1.3b_5_S <- subjectmeans_dLight1.3b %>% filter(Dose==5,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_dLight1.3b_5_M <- subjectmeans_dLight1.3b %>% filter(Dose==5,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_dLight1.3b_7.5_S <- subjectmeans_dLight1.3b %>% filter(Dose==7.5,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_dLight1.3b_7.5_M <- subjectmeans_dLight1.3b %>% filter(Dose==7.5,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_dLight1.3b_10_S <- subjectmeans_dLight1.3b %>% filter(Dose==10,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_dLight1.3b_10_M <- subjectmeans_dLight1.3b %>% filter(Dose==10,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_dLight1.3b_2.5_S
n_dLight1.3b_2.5_M

n_dLight1.3b_5_S
n_dLight1.3b_5_M

n_dLight1.3b_7.5_S
n_dLight1.3b_7.5_M

n_dLight1.3b_10_S
n_dLight1.3b_10_M

### NAcLS GRABDA2h
subjectmeans_GRABDA2h <- subjectmeans %>% filter(Virus=='GRABDA2h')
n_GRABDA2h_2.5_S <- subjectmeans_GRABDA2h %>% filter(Dose==2.5,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_GRABDA2h_2.5_M <- subjectmeans_GRABDA2h %>% filter(Dose==2.5,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_GRABDA2h_5_S <- subjectmeans_GRABDA2h %>% filter(Dose==5,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_GRABDA2h_5_M <- subjectmeans_GRABDA2h %>% filter(Dose==5,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_GRABDA2h_7.5_S <- subjectmeans_GRABDA2h %>% filter(Dose==7.5,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_GRABDA2h_7.5_M <- subjectmeans_GRABDA2h %>% filter(Dose==7.5,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_GRABDA2h_10_S <- subjectmeans_GRABDA2h %>% filter(Dose==10,InjType=='S') %>% summarise(n = n_distinct(SubjectID))
n_GRABDA2h_10_M <- subjectmeans_GRABDA2h %>% filter(Dose==10,InjType=='M') %>% summarise(n = n_distinct(SubjectID))

n_GRABDA2h_2.5_S
n_GRABDA2h_2.5_M

n_GRABDA2h_5_S
n_GRABDA2h_5_M

n_GRABDA2h_7.5_S
n_GRABDA2h_7.5_M

n_GRABDA2h_10_S
n_GRABDA2h_10_M

## Prepare Rise Analysis data -----
# Add dose centered
dosemean <- mean(unique(subjectmeansbin$DoseNum))
doselevels_c <- list(DoseNum_c = c(-1.25, -3.75, 1.25, 3.75))
subjectmeansbin$DoseNum_c <- subjectmeansbin$DoseNum - dosemean

subjectrise_pksperminC <- quantifyRise(subjectmeansbin, 'pksperminC', 12, 3)
subjectrise_ampC <- quantifyRise(subjectmeansbin, 'ampC', 12, 3)
subjectrise_AUCwindowC <- quantifyRise(subjectmeansbin, 'AUCwindowC', 12, 3)

subjectrise_pksperminC_GCaMP6f <-subjectrise_pksperminC %>% filter(Virus=='GCaMP6f'&InjType=='M')
subjectrise_pksperminC_dLight1.3b <-subjectrise_pksperminC %>% filter(Virus=='dLight1.3b'&InjType=='M')
subjectrise_pksperminC_GRABDA2h <-subjectrise_pksperminC %>% filter(Virus=='GRABDA2h'&InjType=='M')

subjectrise_ampC_GCaMP6f <-subjectrise_ampC %>% filter(Virus=='GCaMP6f'&InjType=='M')
subjectrise_ampC_dLight1.3b <-subjectrise_ampC %>% filter(Virus=='dLight1.3b'&InjType=='M')
subjectrise_ampC_GRABDA2h <-subjectrise_ampC %>% filter(Virus=='GRABDA2h'&InjType=='M')

subjectrise_AUCwindowC_GCaMP6f <-subjectrise_AUCwindowC %>% filter(Virus=='GCaMP6f'&InjType=='M')
subjectrise_AUCwindowC_dLight1.3b <-subjectrise_AUCwindowC %>% filter(Virus=='dLight1.3b'&InjType=='M')
subjectrise_AUCwindowC_GRABDA2h <-subjectrise_AUCwindowC %>% filter(Virus=='GRABDA2h'&InjType=='M')

## Prepare Saline Control Individual Sensor Means -----
### Read in data and prep variables -----
pkdatazS_raw <- read_csv('transientquantification_sigfiltz_normbaseline_injcropped_SDthreshold_S1Only.csv')

# Filter excluded subjects
pkdatazS <- pkdatazS_raw %>% filter(!SubjectID %in% subjectexc, !is.na(Bin_3min))

# Add CondNum and Bin
pkdatazS$Bin <- pkdatazS$Bin_3min
pkdatazS$CondNum <- if_else(pkdatazS$Bin <= 4, 0, 1)

# Check variables
unique(pkdatazS$SubjectID)
unique(pkdatazS$Phase)
unique(pkdatazS$Sex)
unique(pkdatazS$FiberPlacement)
unique(pkdatazS$Virus)
unique(pkdatazS$Dose)
unique(pkdatazS$TreatNum)
unique(pkdatazS$CondNum)
unique(pkdatazS$Bin)

# Factor variables
pkdatazS$Sex <- factor(pkdatazS$Sex, levels = c('F','M'))
pkdatazS$Virus <- factor(pkdatazS$Virus, levels = c('GCaMP6f','dLight1.3b','GRABDA2h'))
pkdatazS$CondNum <- factor(pkdatazS$CondNum)

pkdatazS$BinF <- factor(pkdatazS$Bin)

# If fall is NA, set to 2000
pkdatazS$fallms <- if_else(is.na(pkdatazS$fallms), 2000, pkdatazS$fallms)

# If slopes are Inf, set to NA
pkdatazS$riseslope <- if_else(!is.infinite(pkdatazS$riseslope), pkdatazS$riseslope, NA)
pkdatazS$fallslope <- if_else(!is.infinite(pkdatazS$fallslope), pkdatazS$fallslope, NA)


## Find subject means -----
### Find subject values by bin - 3 minute bins
subjectmeansbinSonly <- pkdatazS %>%
  group_by(SubjectID, Sex, FiberPlacement, Virus, Bin) %>%
  dplyr::summarize(freq = n_distinct(maxloc, na.rm=TRUE), 
                   amp = mean(amp, na.rm=TRUE),
                   risems = mean(risems, na.rm=TRUE),
                   riseslope = mean(riseslope, na.rm=TRUE),
                   fallms = mean(fallms, na.rm=TRUE),
                   fallslope = mean(fallslope, na.rm=TRUE),
                   widthms = mean(widthms, na.rm=TRUE),
                   AUC = mean(AUC, na.rm=TRUE),
                   AUCwindow = mean(AUCwindow, na.rm=TRUE),
                   IEIs = mean(IEIs, na.rm=TRUE),
                   ncompoundevents = sum(compoundeventnum)) %>%
  ungroup() %>%
  complete(Bin, nesting(SubjectID, Sex, FiberPlacement, Virus), fill = list(freq = 0))

subjectmeansbinSonly$CondNum <- factor(if_else(subjectmeansbinSonly$Bin <= 4, 0, 1))
subjectmeansbinSonly$pkspermin <- subjectmeansbinSonly$freq/3 # Add pks per minute
subjectmeansbinSonly$BinF <- factor(subjectmeansbinSonly$Bin)

subjectmeansbinSonly$BinTime <- subjectmeansbinSonly$Bin * 3
subjectmeansbinSonly$BinTimeStart <- (subjectmeansbinSonly$Bin*3)-3
subjectmeansbinSonly$BinTimeEnd <- subjectmeansbinSonly$BinTimeStart + 3
subjectmeansbinSonly$BinTimeInj <- (subjectmeansbinSonly$BinTimeStart-12) + 1.5


# Subject Mean by Sensor
salinesubjectmeans <- subjectmeansbinSonly %>%
  group_by(SubjectID, Sex, FiberPlacement, Virus) %>%
  dplyr::summarize(freq = mean(freq, na.rm=TRUE),
                   pkspermin = mean(pkspermin, na.rm=TRUE), 
                   amp = mean(amp, na.rm=TRUE), 
                   risems = mean(risems, na.rm=TRUE), 
                   riseslope = mean(riseslope, na.rm=TRUE), 
                   fallms = mean(fallms, na.rm=TRUE), 
                   fallslope = mean(fallslope, na.rm=TRUE), 
                   widthms = mean(widthms, na.rm=TRUE), 
                   AUC = mean(AUC, na.rm=TRUE),
                   AUCwindow = mean(AUCwindow, na.rm=TRUE),
                   IEIs = mean(IEIs, na.rm=TRUE), 
                   ncompoundevents = mean(ncompoundevents, na.rm=TRUE),
                   ncompoundevents_sum = sum(ncompoundevents, na.rm=TRUE),
  )

salinesubjectmeans
salinesubjectmeans$ampscaled <- salinesubjectmeans$amp - threshold

## Mean by sensor 
salinemeans <- salinesubjectmeans %>%
  group_by(FiberPlacement, Virus) %>%
  dplyr::summarize(freq_mean = mean(freq, na.rm=TRUE), freq_sd = sd(freq, na.rm=TRUE), freq_n = sum(!is.na(freq)), freq_se = freq_sd/(sqrt(freq_n)),
                   pkspermin_mean = mean(pkspermin, na.rm=TRUE), pkspermin_sd = sd(pkspermin, na.rm=TRUE), pkspermin_n = sum(!is.na(pkspermin)), pkspermin_se = pkspermin_sd/(sqrt(pkspermin_n)),
                   amp_mean = mean(amp, na.rm=TRUE), amp_sd = sd(amp, na.rm=TRUE), amp_n = sum(!is.na(amp)), amp_se = amp_sd/(sqrt(amp_n)),
                   risems_mean = mean(risems, na.rm=TRUE), risems_sd = sd(risems, na.rm=TRUE), risems_n = sum(!is.na(risems)), risems_se = risems_sd/(sqrt(risems_n)),
                   riseslope_mean = mean(riseslope, na.rm=TRUE), riseslope_sd = sd(riseslope, na.rm=TRUE), riseslope_n = sum(!is.na(riseslope)), riseslope_se = riseslope_sd/(sqrt(riseslope_n)),
                   fallms_mean = mean(fallms, na.rm=TRUE), fallms_sd = sd(fallms, na.rm=TRUE), fallms_n = sum(!is.na(fallms)), fallms_se = fallms_sd/(sqrt(fallms_n)),
                   fallslope_mean = mean(fallslope, na.rm=TRUE), fallslope_sd = sd(fallslope, na.rm=TRUE), fallslope_n = sum(!is.na(fallslope)), fallslope_se = fallslope_sd/(sqrt(fallslope_n)),
                   widthms_mean = mean(widthms, na.rm=TRUE), widthms_sd = sd(widthms, na.rm=TRUE), widthms_n = sum(!is.na(widthms)), widthms_se = widthms_sd/(sqrt(widthms_n)),
                   AUC_mean = mean(AUC, na.rm=TRUE), AUC_sd = sd(AUC, na.rm=TRUE), AUC_n = sum(!is.na(AUC)), AUC_se = AUC_sd/(sqrt(AUC_n)),
                   AUCwindow_mean = mean(AUCwindow, na.rm=TRUE), AUCwindow_sd = sd(AUCwindow, na.rm=TRUE), AUCwindow_n = sum(!is.na(AUCwindow)), AUCwindow_se = AUCwindow_sd/(sqrt(AUCwindow_n)),
                   IEIs_mean = mean(IEIs, na.rm=TRUE), IEIs_sd = sd(IEIs, na.rm=TRUE), IEIs_n = sum(!is.na(IEIs)), IEIs_se = IEIs_sd/(sqrt(IEIs_n)),
                   ncompoundevents_mean = mean(ncompoundevents, na.rm=TRUE), ncompoundevents_sd = sd(ncompoundevents, na.rm=TRUE), ncompoundevents_n = sum(!is.na(ncompoundevents)), ncompoundevents_se = ncompoundevents_sd/(sqrt(ncompoundevents_n)),
                   ncompoundevents_summean = mean(ncompoundevents, na.rm=TRUE), ncompoundevents_sumsd = sd(ncompoundevents, na.rm=TRUE), ncompoundevents_sumn = sum(!is.na(ncompoundevents)), ncompoundevents_sumse = ncompoundevents_sumsd/(sqrt(ncompoundevents_sumn)),
                   )

salinemeans
salinemeans$ampscaled_mean <- salinemeans$amp_mean - threshold

# WHOLE STREAM: PREPARE DATA -----
## Read in data  -----
streamdata_raw <- read_csv('WholeSessionValues_allstreams_180sBins.csv')
#view(AUCdata_raw)

# Filter excluded subjects
streamdata <- streamdata_raw %>% filter(!SubjectID %in% subjectexc, Stream == 'sigsubz_normbaseline_smooth_20hz_180')

# Add CondNum
streamdata$CondNum <- if_else(streamdata$Bin <= 4, 0, 1)
streamdata$Condition <- if_else(streamdata$CondNum == 0, 'PRE',if_else(streamdata$CondNum == 1, 'POST', NA))

# Factor variables
streamdata$Sex <- factor(streamdata$Sex, levels = c('F','M'))
streamdata$TreatNum <- factor(streamdata$TreatNum)
streamdata$CondNum <- factor(streamdata$CondNum)
streamdata$DoseNum <- streamdata$Dose
streamdata$Dose <- factor(streamdata$Dose, levels = c(2.5,5,7.5,10))
streamdata$Virus <- factor(streamdata$Virus, levels = c('GCaMP6f','dLight1.3b','GRABDA2h'))
streamdata$InjType <- factor(streamdata$InjType, levels = c('S','M'))
streamdata$BinF <- factor(streamdata$Bin)

dosemean <- mean(unique(streamdata$DoseNum))
streamdata$DoseNum_c <- streamdata$DoseNum - dosemean

streamdata$VirusDosecolor <- factor(paste(streamdata$Virus, streamdata$Dose, sep='_'), levels=c('GCaMP6f_2.5','GCaMP6f_5','GCaMP6f_7.5','GCaMP6f_10',
                                                                                                'dLight1.3b_2.5','dLight1.3b_5','dLight1.3b_7.5','dLight1.3b_10',
                                                                                                'GRABDA2h_2.5','GRABDA2h_5','GRABDA2h_7.5','GRABDA2h_10'))

streamdata$BinTime <- streamdata$Bin * 3
streamdata$BinTimeStart <- (streamdata$Bin*3)-3
streamdata$BinTimeEnd <- streamdata$BinTimeStart + 3
streamdata$BinTimeInj <- (streamdata$BinTimeStart-12) + 1.5

# Check variables
unique(streamdata$SubjectID)
unique(streamdata$Sex)
unique(streamdata$FiberPlacement)
unique(streamdata$Virus)
unique(streamdata$InjType)
unique(streamdata$TreatNum)
unique(streamdata$VirusDosecolor)
unique(streamdata$CondNum)
unique(streamdata$Dose)
unique(streamdata$DoseNum)
unique(streamdata$DoseNum_c)
unique(streamdata$Stream)
unique(streamdata$Bin)
unique(streamdata$BinF)
unique(streamdata$BinTimeInj)

## Find subject change values -----
epsilon   <- 1e-3
variables <- c('mean', 'median', 'sd', 'sdRatio','MAD','MADRatio','cv','binAUC','binAUCint')
keys <- c("SubjectID","Dose","Stream","CondNum","Bin")

setDT(streamdata)

# Build saline side
S_side <- streamdata[InjType == "S", c(keys, variables), with = FALSE]
setnames(S_side, variables, paste0(variables, "_S"))

# Join saline values onto all rows
streamdata <- S_side[streamdata, on = keys]

# Vectorized calculations only on M rows
for (v in variables) {
  Sname <- paste0(v, "_S")
  cname <- paste0(v, "C")
  pcname <- paste0(v, "PC")
  streamdata[InjType == "M", (cname)  := get(v) - get(Sname)]
  streamdata[InjType == "M", (pcname) := (get(v) - get(Sname)) / (get(Sname) + epsilon) * 100]
  streamdata[, (Sname) := NULL]
}

streamdata <- as.data.frame(streamdata)

## Find means -----
## Subject Overall Means
streamsubjectmeans <- streamdata %>%
  group_by(SubjectID, Sex, FiberPlacement, Virus, VirusDosecolor, Dose, DoseNum, DoseNum_c, TreatNum, InjType, Stream, Condition, CondNum) %>%
  dplyr::summarize(binmean = mean(mean, na.rm=TRUE), binmeanC = mean(meanC, na.rm=TRUE),
                   binmedian = mean(median, na.rm=TRUE), binmedianC = mean(medianC, na.rm=TRUE),
                   binsd = mean(sd, na.rm=TRUE), binsdC = mean(sdC, na.rm=TRUE),
                   binsdRatio = mean(sdRatio, na.rm=TRUE), binsdRatioC = mean(sdRatioC, na.rm=TRUE),
                   binMAD = mean(MAD, na.rm=TRUE), binMADC = mean(MADC, na.rm=TRUE),
                   binMADRatio = mean(MADRatio, na.rm=TRUE), binMADRatioC = mean(MADRatioC, na.rm=TRUE),
                   bincv = mean(cv, na.rm=TRUE), bincvC = mean(cvC, na.rm=TRUE),
                   binAUC = mean(binAUC, na.rm=TRUE), binAUCC = mean(binAUCC, na.rm=TRUE),
                   binAUCint = mean(binAUCint, na.rm=TRUE), binAUCintC = mean(binAUCintC, na.rm=TRUE),
                   binAUCsum = sum(binAUC, na.rm=TRUE), binAUCintsum = sum(binAUCint, na.rm=TRUE))

streamsubjectmeans

## Sensor and Dose Overall Means
streamoverallmeans <- streamsubjectmeans %>%
  group_by(FiberPlacement, Virus, VirusDosecolor, Dose, DoseNum, DoseNum_c, TreatNum, InjType, Stream, Condition, CondNum) %>%
  dplyr::summarize(binmean_mean = mean(binmean, na.rm=TRUE), binmean_sd = sd(binmean, na.rm=TRUE), 
                   binmean_n = length(binmean), binmean_se = binmean_sd/(sqrt(binmean_n)),
                   
                   binmeanC_mean = mean(binmeanC, na.rm=TRUE), binmeanC_sd = sd(binmeanC, na.rm=TRUE), 
                   binmeanC_n = length(binmeanC), binmeanC_se = binmeanC_sd/(sqrt(binmeanC_n)),
                   
                   binmedian_mean = mean(binmedian, na.rm=TRUE), binmedian_sd = sd(binmedian, na.rm=TRUE), 
                   binmedian_n = length(binmedian), binmedian_se = binmedian_sd/(sqrt(binmedian_n)),
                   
                   binmedianC_mean = mean(binmedianC, na.rm=TRUE), binmedianC_sd = sd(binmedianC, na.rm=TRUE), 
                   binmedianC_n = length(binmedianC), binmedianC_se = binmedianC_sd/(sqrt(binmedianC_n)),
                   
                   binsd_mean = mean(binsd, na.rm=TRUE), binsd_sd = sd(binsd, na.rm=TRUE), 
                   binsd_n = length(binsd), binsd_se = binsd_sd/(sqrt(binsd_n)),
                   
                   binsdC_mean = mean(binsdC, na.rm=TRUE), binsdC_sd = sd(binsdC, na.rm=TRUE), 
                   binsdC_n = length(binsdC), binsdC_se = binsdC_sd/(sqrt(binsdC_n)),
                   
                   binsdRatio_mean = mean(binsdRatio, na.rm=TRUE), binsdRatio_sd = sd(binsdRatio, na.rm=TRUE), 
                   binsdRatio_n = length(binsdRatio), binsdRatio_se = binsdRatio_sd/(sqrt(binsdRatio_n)),
                   
                   binsdRatioC_mean = mean(binsdRatioC, na.rm=TRUE), binsdRatioC_sd = sd(binsdRatioC, na.rm=TRUE), 
                   binsdRatioC_n = length(binsdRatioC), binsdRatioC_se = binsdRatioC_sd/(sqrt(binsdRatioC_n)),
                   
                   binMAD_mean = mean(binMAD, na.rm=TRUE), binMAD_sd = sd(binMAD, na.rm=TRUE), 
                   binMAD_n = length(binMAD), binMAD_se = binMAD_sd/(sqrt(binMAD_n)),
                   
                   binMADC_mean = mean(binMADC, na.rm=TRUE), binMADC_sd = sd(binMADC, na.rm=TRUE), 
                   binMADC_n = length(binMADC), binMADC_se = binMADC_sd/(sqrt(binMADC_n)),
                   
                   binMADRatio_mean = mean(binMADRatio, na.rm=TRUE), binMADRatio_sd = sd(binMADRatio, na.rm=TRUE), 
                   binMADRatio_n = length(binMADRatio), binMADRatio_se = binMADRatio_sd/(sqrt(binMADRatio_n)),
                   
                   binMADRatioC_mean = mean(binMADRatioC, na.rm=TRUE), binMADRatioC_sd = sd(binMADRatioC, na.rm=TRUE), 
                   binMADRatioC_n = length(binMADRatioC), binMADRatioC_se = binMADRatioC_sd/(sqrt(binMADRatioC_n)),
                   
                   bincv_mean = mean(bincv, na.rm=TRUE), bincv_sd = sd(bincv, na.rm=TRUE), 
                   bincv_n = length(bincv), bincv_se = bincv_sd/(sqrt(bincv_n)),
                   
                   bincvC_mean = mean(bincvC, na.rm=TRUE), bincvC_sd = sd(bincvC, na.rm=TRUE), 
                   bincvC_n = length(bincvC), bincvC_se = bincvC_sd/(sqrt(bincvC_n)),
                   
                   binAUC_mean = mean(binAUC, na.rm=TRUE), binAUC_sd = sd(binAUC, na.rm=TRUE), 
                   binAUC_n = length(binAUC), binAUC_se = binAUC_sd/(sqrt(binAUC_n)),
                   
                   binAUCC_mean = mean(binAUCC, na.rm=TRUE), binAUCC_sd = sd(binAUCC, na.rm=TRUE), 
                   binAUCC_n = length(binAUCC), binAUCC_se = binAUCC_sd/(sqrt(binAUCC_n)),
                   
                   binAUCint_mean = mean(binAUCint, na.rm=TRUE), binAUCint_sd = sd(binAUCint, na.rm=TRUE), 
                   binAUCint_n = length(binAUCint), binAUCint_se = binAUCint_sd/(sqrt(binAUCint_n)),
                   
                   binAUCintC_mean = mean(binAUCintC, na.rm=TRUE), binAUCintC_sd = sd(binAUCintC, na.rm=TRUE), 
                   binAUCintC_n = length(binAUCintC), binAUCintC_se = binAUCintC_sd/(sqrt(binAUCintC_n)),
                   
                   binAUCsum_mean = mean(binAUCsum, na.rm=TRUE), binAUCsum_sd = sd(binAUCsum, na.rm=TRUE), 
                   binAUCsum_n = length(binAUCsum), binAUCsum_se = binAUCsum_sd/(sqrt(binAUCsum_n)),
                   
                   binAUCintsum_mean = mean(binAUCintsum, na.rm=TRUE), binAUCintsum_sd = sd(binAUCintsum, na.rm=TRUE), 
                   binAUCintsum_n = length(binAUCintsum), binAUCintsum_se = binAUCintsum_sd/(sqrt(binAUCintsum_n)))

streamoverallmeans


## Sensor and Dose Bin Means
streambinmeans <- streamdata %>%
  group_by(FiberPlacement, Virus, VirusDosecolor, Dose, DoseNum, DoseNum_c, TreatNum, InjType, Stream, Condition, CondNum, Bin, BinF, BinTimeInj) %>%
  dplyr::summarize(binmean_mean = mean(mean, na.rm=TRUE), binmean_sd = sd(mean, na.rm=TRUE), 
                   binmean_n = length(mean), binmean_se = binmean_sd/(sqrt(binmean_n)),
                   
                   binmeanC_mean = mean(meanC, na.rm=TRUE), binmeanC_sd = sd(meanC, na.rm=TRUE), 
                   binmeanC_n = length(meanC), binmeanC_se = binmeanC_sd/(sqrt(binmeanC_n)),
                   
                   binmedian_mean = mean(median, na.rm=TRUE), binmedian_sd = sd(median, na.rm=TRUE), 
                   binmedian_n = length(median), binmedian_se = binmedian_sd/(sqrt(binmedian_n)),
                   
                   binmedianC_mean = mean(medianC, na.rm=TRUE), binmedianC_sd = sd(medianC, na.rm=TRUE), 
                   binmedianC_n = length(medianC), binmedianC_se = binmedianC_sd/(sqrt(binmedianC_n)),
                   
                   binsd_mean = mean(sd, na.rm=TRUE), binsd_sd = sd(sd, na.rm=TRUE), 
                   binsd_n = length(sd), binsd_se = binsd_sd/(sqrt(binsd_n)),
                   
                   binsdC_mean = mean(sdC, na.rm=TRUE), binsdC_sd = sd(sdC, na.rm=TRUE), 
                   binsdC_n = length(sdC), binsdC_se = binsdC_sd/(sqrt(binsdC_n)),
                   
                   binsdRatio_mean = mean(sdRatio, na.rm=TRUE), binsdRatio_sd = sd(sdRatio, na.rm=TRUE), 
                   binsdRatio_n = length(sdRatio), binsdRatio_se = binsdRatio_sd/(sqrt(binsdRatio_n)),
                   
                   binsdRatioC_mean = mean(sdRatioC, na.rm=TRUE), binsdRatioC_sd = sd(sdRatioC, na.rm=TRUE), 
                   binsdRatioC_n = length(sdRatioC), binsdRatioC_se = binsdRatioC_sd/(sqrt(binsdRatioC_n)),
                   
                   binMAD_mean = mean(MAD, na.rm=TRUE), binMAD_sd = sd(MAD, na.rm=TRUE), 
                   binMAD_n = length(MAD), binMAD_se = binMAD_sd/(sqrt(binMAD_n)),
                   
                   binMADC_mean = mean(MADC, na.rm=TRUE), binMADC_sd = sd(MADC, na.rm=TRUE), 
                   binMADC_n = length(MADC), binMADC_se = binMADC_sd/(sqrt(binMADC_n)),
                   
                   binMADRatio_mean = mean(MADRatio, na.rm=TRUE), binMADRatio_sd = sd(MADRatio, na.rm=TRUE), 
                   binMADRatio_n = length(MADRatio), binMADRatio_se = binMADRatio_sd/(sqrt(binMADRatio_n)),
                   
                   binMADRatioC_mean = mean(MADRatioC, na.rm=TRUE), binMADRatioC_sd = sd(MADRatioC, na.rm=TRUE), 
                   binMADRatioC_n = length(MADRatioC), binMADRatioC_se = binMADRatioC_sd/(sqrt(binMADRatioC_n)),
                   
                   bincv_mean = mean(cv, na.rm=TRUE), bincv_sd = sd(cv, na.rm=TRUE), 
                   bincv_n = length(cv), bincv_se = bincv_sd/(sqrt(bincv_n)),
                   
                   bincvC_mean = mean(cvC, na.rm=TRUE), bincvC_sd = sd(cvC, na.rm=TRUE), 
                   bincvC_n = length(cvC), bincvC_se = bincvC_sd/(sqrt(bincvC_n)),
                   
                   binAUC_mean = mean(binAUC, na.rm=TRUE), binAUC_sd = sd(binAUC, na.rm=TRUE), 
                   binAUC_n = length(binAUC), binAUC_se = binAUC_sd/(sqrt(binAUC_n)),
                   
                   binAUCC_mean = mean(binAUCC, na.rm=TRUE), binAUCC_sd = sd(binAUCC, na.rm=TRUE), 
                   binAUCC_n = length(binAUCC), binAUCC_se = binAUCC_sd/(sqrt(binAUCC_n)),
                   
                   binAUCint_mean = mean(binAUCint, na.rm=TRUE), binAUCint_sd = sd(binAUCint, na.rm=TRUE), 
                   binAUCint_n = length(binAUCint), binAUCint_se = binAUCint_sd/(sqrt(binAUCint_n)),
                   
                   binAUCintC_mean = mean(binAUCintC, na.rm=TRUE), binAUCintC_sd = sd(binAUCintC, na.rm=TRUE), 
                   binAUCintC_n = length(binAUCintC), binAUCintC_se = binAUCintC_sd/(sqrt(binAUCintC_n)))

streambinmeans

## PREPARE GRAPHING VARIABLES -----
treatcolors <- colorlists$treatcolors
treatcolors_raw <- colorlists$treatcolors_raw
overallcolors <- colorlists$overallcolors
sensorcompcolors <- colorlists$sensorcompcolors
sensordosecolors <- colorlists$sensordosecolors

injlinecolor <- '#575757'

sexshapes <- c('M'=0,'F'=1)

bintimexticks <- seq(-12,60,3)
bintimexlabels <- ifelse(bintimexticks %% 6 == 0, bintimexticks, "")

bintimexmin <- -12
bintimexmax <- 60

# FIGURE 1: BASELINE TRANSIENT DYNAMICS BY SENSOR -----
## STATS ------
### Frequency -----
#### Mixed effects model
freqmodel_BetweenSensors <- glmmTMB(pkspermin ~ Virus + BinF + (1|SubjectID), # random intercept by subject
                                            data = subjectmeansbinSonly, family = stats::gaussian,
                                            control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqmodel_BetweenSensors_aov <- Anova(freqmodel_BetweenSensors, type=2)
freqmodel_BetweenSensors_aov
summary(freqmodel_BetweenSensors)

freqmodel_BetweenSensors_r2 <- as.data.frame(r.squaredGLMM(freqmodel_BetweenSensors)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqmodel_BetweenSensors_r2

# Follow up tests - sensor pairs
freqmodel_BetweenSensors_sensorpairs <- emmeans(freqmodel_BetweenSensors, pairwise~Virus, adjust="none") # Pairs
freqmodel_BetweenSensors_sensorcontrasts <- as.data.frame(freqmodel_BetweenSensors_sensorpairs$contrasts) # Extract emmeans
freqmodel_BetweenSensors_sensoremmeans <- freqmodel_BetweenSensors_sensorpairs$emmeans # Extract emmeans
freqmodel_BetweenSensors_sensorpairs_es <- as.data.frame(eff_size(freqmodel_BetweenSensors_sensoremmeans, sigma = sigma(freqmodel_BetweenSensors), edf = df.residual(freqmodel_BetweenSensors))) # Effect Size
freqmodel_BetweenSensors_sensorpairs_es <- freqmodel_BetweenSensors_sensorpairs_es[, c("contrast", "effect.size")]

freqmodel_BetweenSensors_sensorcontrasts$d <- freqmodel_BetweenSensors_sensorpairs_es$effect.size
freqmodel_BetweenSensors_sensorcontrasts

# Format Tables
table_freqmodel_BetweenSensors <- format_regsummary_table(freqmodel_BetweenSensors, "ControlSensorComparison", "Frequency", "StandardModel")
table_freqmodel_BetweenSensors_aov <- format_anova_table(freqmodel_BetweenSensors_aov, "ControlSensorComparison", "Frequency", "StandardModel")
table_freqmodel_BetweenSensors_pairs_sensor <- format_emm_table(freqmodel_BetweenSensors_sensorcontrasts, "ControlSensorComparison", "Frequency", "Pairs", "Sensor")

#### Mixed effects AR1 model
freqmodel_BetweenSensors_ar1 <- glmmTMB(pkspermin ~ Virus + BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                                data = subjectmeansbinSonly, family = stats::gaussian, 
                                                control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqmodel_BetweenSensors_ar1_aov <- Anova(freqmodel_BetweenSensors_ar1, type=2)
freqmodel_BetweenSensors_ar1_aov
summary(freqmodel_BetweenSensors_ar1)

freqmodel_BetweenSensors_ar1_r2 <- as.data.frame(r.squaredGLMM(freqmodel_BetweenSensors_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqmodel_BetweenSensors_ar1_r2

freqmodel_BetweenSensors_ar1_AIC <- as.data.frame(AIC(freqmodel_BetweenSensors, freqmodel_BetweenSensors_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "freqmodel_BetweenSensors" ~ "StandardModel", ModelName == "freqmodel_BetweenSensors_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

freqmodel_BetweenSensors_ar1_AIC


# Follow up tests - sensor pairs
freqmodel_BetweenSensors_ar1_sensorpairs <- emmeans(freqmodel_BetweenSensors_ar1, pairwise~Virus, adjust="none") # Pairs
freqmodel_BetweenSensors_ar1_sensorcontrasts <- as.data.frame(freqmodel_BetweenSensors_ar1_sensorpairs$contrasts) # Extract emmeans
freqmodel_BetweenSensors_ar1_sensoremmeans <- freqmodel_BetweenSensors_ar1_sensorpairs$emmeans # Extract emmeans
freqmodel_BetweenSensors_ar1_sensorpairs_es <- as.data.frame(eff_size(freqmodel_BetweenSensors_ar1_sensoremmeans, sigma = sigma(freqmodel_BetweenSensors_ar1), edf = df.residual(freqmodel_BetweenSensors_ar1))) # Effect Size
freqmodel_BetweenSensors_ar1_sensorpairs_es <- freqmodel_BetweenSensors_ar1_sensorpairs_es[, c("contrast", "effect.size")]

freqmodel_BetweenSensors_ar1_sensorcontrasts$d <- freqmodel_BetweenSensors_ar1_sensorpairs_es$effect.size
freqmodel_BetweenSensors_ar1_sensorcontrasts

# Format Tables
table_freqmodel_BetweenSensors_ar1 <- format_regsummary_table(freqmodel_BetweenSensors_ar1, 'ControlSensorComparison', "Frequency", "AR1Model")
table_freqmodel_BetweenSensors_ar1_aov <- format_anova_table(freqmodel_BetweenSensors_ar1_aov, "ControlSensorComparison", "Frequency", "AR1Model") 
table_freqmodel_BetweenSensors_ar1_pairs_sensor <- format_emm_table(freqmodel_BetweenSensors_ar1_sensorcontrasts, "ControlSensorComparison", "Frequency", "Pairs", "Sensor")

### Amplitude -----
#### Mixed effects model
ampmodel_BetweenSensors <- glmmTMB(amp ~ Virus + BinF + (1|SubjectID), # random intercept by subject
                                           data = subjectmeansbinSonly, family = stats::gaussian,
                                           control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampmodel_BetweenSensors_aov <- Anova(ampmodel_BetweenSensors, type=2)
ampmodel_BetweenSensors_aov
summary(ampmodel_BetweenSensors)

ampmodel_BetweenSensors_r2 <- as.data.frame(r.squaredGLMM(ampmodel_BetweenSensors)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampmodel_BetweenSensors_r2

# Follow up tests - sensor pairs
ampmodel_BetweenSensors_sensorpairs <- emmeans(ampmodel_BetweenSensors, pairwise~Virus, adjust="none") # Pairs
ampmodel_BetweenSensors_sensorcontrasts <- as.data.frame(ampmodel_BetweenSensors_sensorpairs$contrasts) # Extract emmeans
ampmodel_BetweenSensors_sensoremmeans <- ampmodel_BetweenSensors_sensorpairs$emmeans # Extract emmeans
ampmodel_BetweenSensors_sensorpairs_es <- as.data.frame(eff_size(ampmodel_BetweenSensors_sensoremmeans, sigma = sigma(ampmodel_BetweenSensors), edf = df.residual(ampmodel_BetweenSensors))) # Effect Size
ampmodel_BetweenSensors_sensorpairs_es <- ampmodel_BetweenSensors_sensorpairs_es[, c("contrast", "effect.size")]

ampmodel_BetweenSensors_sensorcontrasts$d <- ampmodel_BetweenSensors_sensorpairs_es$effect.size
ampmodel_BetweenSensors_sensorcontrasts

# Format Tables
table_ampmodel_BetweenSensors <- format_regsummary_table(ampmodel_BetweenSensors, "ControlSensorComparison", "Amplitude", "StandardModel")
table_ampmodel_BetweenSensors_aov <- format_anova_table(ampmodel_BetweenSensors_aov, "ControlSensorComparison", "Amplitude", "StandardModel")
table_ampmodel_BetweenSensors_pairs_sensor <- format_emm_table(ampmodel_BetweenSensors_sensorcontrasts, "ControlSensorComparison", "Amplitude", "Pairs", "Sensor")

#### Mixed effects AR1 model
ampmodel_BetweenSensors_ar1 <- glmmTMB(amp ~ Virus + BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                               data = subjectmeansbinSonly, family = stats::gaussian, 
                                               control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampmodel_BetweenSensors_ar1_aov <- Anova(ampmodel_BetweenSensors_ar1, type=2)
ampmodel_BetweenSensors_ar1_aov
summary(ampmodel_BetweenSensors_ar1)

ampmodel_BetweenSensors_ar1_r2 <- as.data.frame(r.squaredGLMM(ampmodel_BetweenSensors_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampmodel_BetweenSensors_ar1_r2

ampmodel_BetweenSensors_ar1_AIC <- as.data.frame(AIC(ampmodel_BetweenSensors, ampmodel_BetweenSensors_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "ampmodel_BetweenSensors" ~ "StandardModel", ModelName == "ampmodel_BetweenSensors_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

ampmodel_BetweenSensors_ar1_AIC


# Follow up tests - sensor pairs
ampmodel_BetweenSensors_ar1_sensorpairs <- emmeans(ampmodel_BetweenSensors_ar1, pairwise~Virus, adjust="none") # Pairs
ampmodel_BetweenSensors_ar1_sensorcontrasts <- as.data.frame(ampmodel_BetweenSensors_ar1_sensorpairs$contrasts) # Extract emmeans
ampmodel_BetweenSensors_ar1_sensoremmeans <- ampmodel_BetweenSensors_ar1_sensorpairs$emmeans # Extract emmeans
ampmodel_BetweenSensors_ar1_sensorpairs_es <- as.data.frame(eff_size(ampmodel_BetweenSensors_ar1_sensoremmeans, sigma = sigma(ampmodel_BetweenSensors_ar1), edf = df.residual(ampmodel_BetweenSensors_ar1))) # Effect Size
ampmodel_BetweenSensors_ar1_sensorpairs_es <- ampmodel_BetweenSensors_ar1_sensorpairs_es[, c("contrast", "effect.size")]

ampmodel_BetweenSensors_ar1_sensorcontrasts$d <- ampmodel_BetweenSensors_ar1_sensorpairs_es$effect.size
ampmodel_BetweenSensors_ar1_sensorcontrasts

# Format Tables
table_ampmodel_BetweenSensors_ar1 <- format_regsummary_table(ampmodel_BetweenSensors_ar1, 'ControlSensorComparison', "Amplitude", "AR1Model")
table_ampmodel_BetweenSensors_ar1_aov <- format_anova_table(ampmodel_BetweenSensors_ar1_aov, "ControlSensorComparison", "Amplitude", "AR1Model") 
table_ampmodel_BetweenSensors_ar1_pairs_sensor <- format_emm_table(ampmodel_BetweenSensors_ar1_sensorcontrasts, "ControlSensorComparison", "Amplitude", "Pairs", "Sensor")

### AUCwindow -----
#### Mixed effects model
AUCwindowmodel_BetweenSensors <- glmmTMB(scale(AUCwindow) ~ Virus + BinF + (1|SubjectID), # random intercept by subject
                                                 data = subjectmeansbinSonly, family = stats::gaussian,
                                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowmodel_BetweenSensors_aov <- Anova(AUCwindowmodel_BetweenSensors, type=2)
AUCwindowmodel_BetweenSensors_aov
summary(AUCwindowmodel_BetweenSensors)

AUCwindowmodel_BetweenSensors_r2 <- as.data.frame(r.squaredGLMM(AUCwindowmodel_BetweenSensors)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowmodel_BetweenSensors_r2

# Follow up tests - sensor pairs
AUCwindowmodel_BetweenSensors_sensorpairs <- emmeans(AUCwindowmodel_BetweenSensors, pairwise~Virus, adjust="none") # Pairs
AUCwindowmodel_BetweenSensors_sensorcontrasts <- as.data.frame(AUCwindowmodel_BetweenSensors_sensorpairs$contrasts) # Extract emmeans
AUCwindowmodel_BetweenSensors_sensoremmeans <- AUCwindowmodel_BetweenSensors_sensorpairs$emmeans # Extract emmeans
AUCwindowmodel_BetweenSensors_sensorpairs_es <- as.data.frame(eff_size(AUCwindowmodel_BetweenSensors_sensoremmeans, sigma = sigma(AUCwindowmodel_BetweenSensors), edf = df.residual(AUCwindowmodel_BetweenSensors))) # Effect Size
AUCwindowmodel_BetweenSensors_sensorpairs_es <- AUCwindowmodel_BetweenSensors_sensorpairs_es[, c("contrast", "effect.size")]

AUCwindowmodel_BetweenSensors_sensorcontrasts$d <- AUCwindowmodel_BetweenSensors_sensorpairs_es$effect.size
AUCwindowmodel_BetweenSensors_sensorcontrasts

# Format Tables
table_AUCwindowmodel_BetweenSensors <- format_regsummary_table(AUCwindowmodel_BetweenSensors, "ControlSensorComparison", "AUCwindow", "StandardModel")
table_AUCwindowmodel_BetweenSensors_aov <- format_anova_table(AUCwindowmodel_BetweenSensors_aov, "ControlSensorComparison", "AUCwindow", "StandardModel")
table_AUCwindowmodel_BetweenSensors_pairs_sensor <- format_emm_table(AUCwindowmodel_BetweenSensors_sensorcontrasts, "ControlSensorComparison", "AUCwindow", "Pairs", "Sensor")

#### Mixed effects AR1 model
AUCwindowmodel_BetweenSensors_ar1 <- glmmTMB(scale(AUCwindow) ~ Virus + BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                                     data = subjectmeansbinSonly, family = stats::gaussian, 
                                                     control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowmodel_BetweenSensors_ar1_aov <- Anova(AUCwindowmodel_BetweenSensors_ar1, type=2)
AUCwindowmodel_BetweenSensors_ar1_aov
summary(AUCwindowmodel_BetweenSensors_ar1)

AUCwindowmodel_BetweenSensors_ar1_r2 <- as.data.frame(r.squaredGLMM(AUCwindowmodel_BetweenSensors_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowmodel_BetweenSensors_ar1_r2

AUCwindowmodel_BetweenSensors_ar1_AIC <- as.data.frame(AIC(AUCwindowmodel_BetweenSensors, AUCwindowmodel_BetweenSensors_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "AUCwindowmodel_BetweenSensors" ~ "StandardModel", ModelName == "AUCwindowmodel_BetweenSensors_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

AUCwindowmodel_BetweenSensors_ar1_AIC


# Follow up tests - sensor pairs
AUCwindowmodel_BetweenSensors_ar1_sensorpairs <- emmeans(AUCwindowmodel_BetweenSensors_ar1, pairwise~Virus, adjust="none") # Pairs
AUCwindowmodel_BetweenSensors_ar1_sensorcontrasts <- as.data.frame(AUCwindowmodel_BetweenSensors_ar1_sensorpairs$contrasts) # Extract emmeans
AUCwindowmodel_BetweenSensors_ar1_sensoremmeans <- AUCwindowmodel_BetweenSensors_ar1_sensorpairs$emmeans # Extract emmeans
AUCwindowmodel_BetweenSensors_ar1_sensorpairs_es <- as.data.frame(eff_size(AUCwindowmodel_BetweenSensors_ar1_sensoremmeans, sigma = sigma(AUCwindowmodel_BetweenSensors_ar1), edf = df.residual(AUCwindowmodel_BetweenSensors_ar1))) # Effect Size
AUCwindowmodel_BetweenSensors_ar1_sensorpairs_es <- AUCwindowmodel_BetweenSensors_ar1_sensorpairs_es[, c("contrast", "effect.size")]

AUCwindowmodel_BetweenSensors_ar1_sensorcontrasts$d <- AUCwindowmodel_BetweenSensors_ar1_sensorpairs_es$effect.size
AUCwindowmodel_BetweenSensors_ar1_sensorcontrasts

# Format Tables
table_AUCwindowmodel_BetweenSensors_ar1 <- format_regsummary_table(AUCwindowmodel_BetweenSensors_ar1, 'ControlSensorComparison', "AUCwindow", "AR1Model")
table_AUCwindowmodel_BetweenSensors_ar1_aov <- format_anova_table(AUCwindowmodel_BetweenSensors_ar1_aov, "ControlSensorComparison", "AUCwindow", "AR1Model") 
table_AUCwindowmodel_BetweenSensors_ar1_pairs_sensor <- format_emm_table(AUCwindowmodel_BetweenSensors_ar1_sensorcontrasts, "ControlSensorComparison", "AUCwindow", "Pairs", "Sensor")


### Rise Duration -----
#### Mixed effects model
risemsmodel_BetweenSensors <- glmmTMB(scale(risems) ~ Virus + BinF + (1|SubjectID), # random intercept by subject
                                              data = subjectmeansbinSonly, family = stats::gaussian,
                                              control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

risemsmodel_BetweenSensors_aov <- Anova(risemsmodel_BetweenSensors, type=2)
risemsmodel_BetweenSensors_aov
summary(risemsmodel_BetweenSensors)

risemsmodel_BetweenSensors_r2 <- as.data.frame(r.squaredGLMM(risemsmodel_BetweenSensors)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
risemsmodel_BetweenSensors_r2

# Follow up tests - sensor pairs
risemsmodel_BetweenSensors_sensorpairs <- emmeans(risemsmodel_BetweenSensors, pairwise~Virus, adjust="none") # Pairs
risemsmodel_BetweenSensors_sensorcontrasts <- as.data.frame(risemsmodel_BetweenSensors_sensorpairs$contrasts) # Extract emmeans
risemsmodel_BetweenSensors_sensoremmeans <- risemsmodel_BetweenSensors_sensorpairs$emmeans # Extract emmeans
risemsmodel_BetweenSensors_sensorpairs_es <- as.data.frame(eff_size(risemsmodel_BetweenSensors_sensoremmeans, sigma = sigma(risemsmodel_BetweenSensors), edf = df.residual(risemsmodel_BetweenSensors))) # Effect Size
risemsmodel_BetweenSensors_sensorpairs_es <- risemsmodel_BetweenSensors_sensorpairs_es[, c("contrast", "effect.size")]

risemsmodel_BetweenSensors_sensorcontrasts$d <- risemsmodel_BetweenSensors_sensorpairs_es$effect.size
risemsmodel_BetweenSensors_sensorcontrasts

# Format Tables
table_risemsmodel_BetweenSensors <- format_regsummary_table(risemsmodel_BetweenSensors, "ControlSensorComparison", "RiseDuration", "StandardModel")
table_risemsmodel_BetweenSensors_aov <- format_anova_table(risemsmodel_BetweenSensors_aov, "ControlSensorComparison", "RiseDuration", "StandardModel")
table_risemsmodel_BetweenSensors_pairs_sensor <- format_emm_table(risemsmodel_BetweenSensors_sensorcontrasts, "ControlSensorComparison", "RiseDuration", "Pairs", "Sensor")

#### Mixed effects AR1 model
risemsmodel_BetweenSensors_ar1 <- glmmTMB(scale(risems) ~ Virus + BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                                  data = subjectmeansbinSonly, family = stats::gaussian, 
                                                  control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

risemsmodel_BetweenSensors_ar1_aov <- Anova(risemsmodel_BetweenSensors_ar1, type=2)
risemsmodel_BetweenSensors_ar1_aov
summary(risemsmodel_BetweenSensors_ar1)

risemsmodel_BetweenSensors_ar1_r2 <- as.data.frame(r.squaredGLMM(risemsmodel_BetweenSensors_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
risemsmodel_BetweenSensors_ar1_r2

risemsmodel_BetweenSensors_ar1_AIC <- as.data.frame(AIC(risemsmodel_BetweenSensors, risemsmodel_BetweenSensors_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "risemsmodel_BetweenSensors" ~ "StandardModel", ModelName == "risemsmodel_BetweenSensors_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

risemsmodel_BetweenSensors_ar1_AIC


# Follow up tests - sensor pairs
risemsmodel_BetweenSensors_ar1_sensorpairs <- emmeans(risemsmodel_BetweenSensors_ar1, pairwise~Virus, adjust="none") # Pairs
risemsmodel_BetweenSensors_ar1_sensorcontrasts <- as.data.frame(risemsmodel_BetweenSensors_ar1_sensorpairs$contrasts) # Extract emmeans
risemsmodel_BetweenSensors_ar1_sensoremmeans <- risemsmodel_BetweenSensors_ar1_sensorpairs$emmeans # Extract emmeans
risemsmodel_BetweenSensors_ar1_sensorpairs_es <- as.data.frame(eff_size(risemsmodel_BetweenSensors_ar1_sensoremmeans, sigma = sigma(risemsmodel_BetweenSensors_ar1), edf = df.residual(risemsmodel_BetweenSensors_ar1))) # Effect Size
risemsmodel_BetweenSensors_ar1_sensorpairs_es <- risemsmodel_BetweenSensors_ar1_sensorpairs_es[, c("contrast", "effect.size")]

risemsmodel_BetweenSensors_ar1_sensorcontrasts$d <- risemsmodel_BetweenSensors_ar1_sensorpairs_es$effect.size
risemsmodel_BetweenSensors_ar1_sensorcontrasts

# Format Tables
table_risemsmodel_BetweenSensors_ar1 <- format_regsummary_table(risemsmodel_BetweenSensors_ar1, 'ControlSensorComparison', "RiseDuration", "AR1Model")
table_risemsmodel_BetweenSensors_ar1_aov <- format_anova_table(risemsmodel_BetweenSensors_ar1_aov, "ControlSensorComparison", "RiseDuration", "AR1Model") 
table_risemsmodel_BetweenSensors_ar1_pairs_sensor <- format_emm_table(risemsmodel_BetweenSensors_ar1_sensorcontrasts, "ControlSensorComparison", "RiseDuration", "Pairs", "Sensor")

### Rise Slope -----
#### Mixed effects model
riseslopemodel_BetweenSensors <- glmmTMB(riseslope ~ Virus + BinF + (1|SubjectID), # random intercept by subject
                                                 data = subset(subjectmeansbinSonly, !is.na(riseslope)), family = stats::gaussian,
                                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

riseslopemodel_BetweenSensors_aov <- Anova(riseslopemodel_BetweenSensors, type=2)
riseslopemodel_BetweenSensors_aov
summary(riseslopemodel_BetweenSensors)

riseslopemodel_BetweenSensors_r2 <- as.data.frame(r.squaredGLMM(riseslopemodel_BetweenSensors)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
riseslopemodel_BetweenSensors_r2

# Follow up tests - sensor pairs
riseslopemodel_BetweenSensors_sensorpairs <- emmeans(riseslopemodel_BetweenSensors, pairwise~Virus, adjust="none") # Pairs
riseslopemodel_BetweenSensors_sensorcontrasts <- as.data.frame(riseslopemodel_BetweenSensors_sensorpairs$contrasts) # Extract emmeans
riseslopemodel_BetweenSensors_sensoremmeans <- riseslopemodel_BetweenSensors_sensorpairs$emmeans # Extract emmeans
riseslopemodel_BetweenSensors_sensorpairs_es <- as.data.frame(eff_size(riseslopemodel_BetweenSensors_sensoremmeans, sigma = sigma(riseslopemodel_BetweenSensors), edf = df.residual(riseslopemodel_BetweenSensors))) # Effect Size
riseslopemodel_BetweenSensors_sensorpairs_es <- riseslopemodel_BetweenSensors_sensorpairs_es[, c("contrast", "effect.size")]

riseslopemodel_BetweenSensors_sensorcontrasts$d <- riseslopemodel_BetweenSensors_sensorpairs_es$effect.size
riseslopemodel_BetweenSensors_sensorcontrasts

# Format Tables
table_riseslopemodel_BetweenSensors <- format_regsummary_table(riseslopemodel_BetweenSensors, "ControlSensorComparison", "RiseSlope", "StandardModel")
table_riseslopemodel_BetweenSensors_aov <- format_anova_table(riseslopemodel_BetweenSensors_aov, "ControlSensorComparison", "RiseSlope", "StandardModel")
table_riseslopemodel_BetweenSensors_pairs_sensor <- format_emm_table(riseslopemodel_BetweenSensors_sensorcontrasts, "ControlSensorComparison", "RiseSlope", "Pairs", "Sensor")

#### Mixed effects AR1 model
riseslopemodel_BetweenSensors_ar1 <- glmmTMB(riseslope ~ Virus + BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                                     data = subset(subjectmeansbinSonly, !is.na(riseslope)), family = stats::gaussian, 
                                                     control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

riseslopemodel_BetweenSensors_ar1_aov <- Anova(riseslopemodel_BetweenSensors_ar1, type=2)
riseslopemodel_BetweenSensors_ar1_aov
summary(riseslopemodel_BetweenSensors_ar1)

riseslopemodel_BetweenSensors_ar1_r2 <- as.data.frame(r.squaredGLMM(riseslopemodel_BetweenSensors_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
riseslopemodel_BetweenSensors_ar1_r2

riseslopemodel_BetweenSensors_ar1_AIC <- as.data.frame(AIC(riseslopemodel_BetweenSensors, riseslopemodel_BetweenSensors_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "riseslopemodel_BetweenSensors" ~ "StandardModel", ModelName == "riseslopemodel_BetweenSensors_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

riseslopemodel_BetweenSensors_ar1_AIC


# Follow up tests - sensor pairs
riseslopemodel_BetweenSensors_ar1_sensorpairs <- emmeans(riseslopemodel_BetweenSensors_ar1, pairwise~Virus, adjust="none") # Pairs
riseslopemodel_BetweenSensors_ar1_sensorcontrasts <- as.data.frame(riseslopemodel_BetweenSensors_ar1_sensorpairs$contrasts) # Extract emmeans
riseslopemodel_BetweenSensors_ar1_sensoremmeans <- riseslopemodel_BetweenSensors_ar1_sensorpairs$emmeans # Extract emmeans
riseslopemodel_BetweenSensors_ar1_sensorpairs_es <- as.data.frame(eff_size(riseslopemodel_BetweenSensors_ar1_sensoremmeans, sigma = sigma(riseslopemodel_BetweenSensors_ar1), edf = df.residual(riseslopemodel_BetweenSensors_ar1))) # Effect Size
riseslopemodel_BetweenSensors_ar1_sensorpairs_es <- riseslopemodel_BetweenSensors_ar1_sensorpairs_es[, c("contrast", "effect.size")]

riseslopemodel_BetweenSensors_ar1_sensorcontrasts$d <- riseslopemodel_BetweenSensors_ar1_sensorpairs_es$effect.size
riseslopemodel_BetweenSensors_ar1_sensorcontrasts

# Format Tables
table_riseslopemodel_BetweenSensors_ar1 <- format_regsummary_table(riseslopemodel_BetweenSensors_ar1, 'ControlSensorComparison', "RiseSlope", "AR1Model")
table_riseslopemodel_BetweenSensors_ar1_aov <- format_anova_table(riseslopemodel_BetweenSensors_ar1_aov, "ControlSensorComparison", "RiseSlope", "AR1Model") 
table_riseslopemodel_BetweenSensors_ar1_pairs_sensor <- format_emm_table(riseslopemodel_BetweenSensors_ar1_sensorcontrasts, "ControlSensorComparison", "RiseSlope", "Pairs", "Sensor")


### Fall Duration -----
#### Mixed effects model
fallmsmodel_BetweenSensors <- glmmTMB(scale(fallms) ~ Virus + BinF + (1|SubjectID), # random intercept by subject
                                              data = subjectmeansbinSonly, family = stats::gaussian,
                                              control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

fallmsmodel_BetweenSensors_aov <- Anova(fallmsmodel_BetweenSensors, type=2)
fallmsmodel_BetweenSensors_aov
summary(fallmsmodel_BetweenSensors)

fallmsmodel_BetweenSensors_r2 <- as.data.frame(r.squaredGLMM(fallmsmodel_BetweenSensors)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
fallmsmodel_BetweenSensors_r2

# Follow up tests - sensor pairs
fallmsmodel_BetweenSensors_sensorpairs <- emmeans(fallmsmodel_BetweenSensors, pairwise~Virus, adjust="none") # Pairs
fallmsmodel_BetweenSensors_sensorcontrasts <- as.data.frame(fallmsmodel_BetweenSensors_sensorpairs$contrasts) # Extract emmeans
fallmsmodel_BetweenSensors_sensoremmeans <- fallmsmodel_BetweenSensors_sensorpairs$emmeans # Extract emmeans
fallmsmodel_BetweenSensors_sensorpairs_es <- as.data.frame(eff_size(fallmsmodel_BetweenSensors_sensoremmeans, sigma = sigma(fallmsmodel_BetweenSensors), edf = df.residual(fallmsmodel_BetweenSensors))) # Effect Size
fallmsmodel_BetweenSensors_sensorpairs_es <- fallmsmodel_BetweenSensors_sensorpairs_es[, c("contrast", "effect.size")]

fallmsmodel_BetweenSensors_sensorcontrasts$d <- fallmsmodel_BetweenSensors_sensorpairs_es$effect.size
fallmsmodel_BetweenSensors_sensorcontrasts

# Format Tables
table_fallmsmodel_BetweenSensors <- format_regsummary_table(fallmsmodel_BetweenSensors, "ControlSensorComparison", "FallDuration", "StandardModel")
table_fallmsmodel_BetweenSensors_aov <- format_anova_table(fallmsmodel_BetweenSensors_aov, "ControlSensorComparison", "FallDuration", "StandardModel")
table_fallmsmodel_BetweenSensors_pairs_sensor <- format_emm_table(fallmsmodel_BetweenSensors_sensorcontrasts, "ControlSensorComparison", "FallDuration", "Pairs", "Sensor")

#### Mixed effects AR1 model
fallmsmodel_BetweenSensors_ar1 <- glmmTMB(scale(fallms) ~ Virus + BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                                  data = subjectmeansbinSonly, family = stats::gaussian, 
                                                  control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

fallmsmodel_BetweenSensors_ar1_aov <- Anova(fallmsmodel_BetweenSensors_ar1, type=2)
fallmsmodel_BetweenSensors_ar1_aov
summary(fallmsmodel_BetweenSensors_ar1)

fallmsmodel_BetweenSensors_ar1_r2 <- as.data.frame(r.squaredGLMM(fallmsmodel_BetweenSensors_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
fallmsmodel_BetweenSensors_ar1_r2

fallmsmodel_BetweenSensors_ar1_AIC <- as.data.frame(AIC(fallmsmodel_BetweenSensors, fallmsmodel_BetweenSensors_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "fallmsmodel_BetweenSensors" ~ "StandardModel", ModelName == "fallmsmodel_BetweenSensors_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

fallmsmodel_BetweenSensors_ar1_AIC


# Follow up tests - sensor pairs
fallmsmodel_BetweenSensors_ar1_sensorpairs <- emmeans(fallmsmodel_BetweenSensors_ar1, pairwise~Virus, adjust="none") # Pairs
fallmsmodel_BetweenSensors_ar1_sensorcontrasts <- as.data.frame(fallmsmodel_BetweenSensors_ar1_sensorpairs$contrasts) # Extract emmeans
fallmsmodel_BetweenSensors_ar1_sensoremmeans <- fallmsmodel_BetweenSensors_ar1_sensorpairs$emmeans # Extract emmeans
fallmsmodel_BetweenSensors_ar1_sensorpairs_es <- as.data.frame(eff_size(fallmsmodel_BetweenSensors_ar1_sensoremmeans, sigma = sigma(fallmsmodel_BetweenSensors_ar1), edf = df.residual(fallmsmodel_BetweenSensors_ar1))) # Effect Size
fallmsmodel_BetweenSensors_ar1_sensorpairs_es <- fallmsmodel_BetweenSensors_ar1_sensorpairs_es[, c("contrast", "effect.size")]

fallmsmodel_BetweenSensors_ar1_sensorcontrasts$d <- fallmsmodel_BetweenSensors_ar1_sensorpairs_es$effect.size
fallmsmodel_BetweenSensors_ar1_sensorcontrasts

# Format Tables
table_fallmsmodel_BetweenSensors_ar1 <- format_regsummary_table(fallmsmodel_BetweenSensors_ar1, 'ControlSensorComparison', "FallDuration", "AR1Model")
table_fallmsmodel_BetweenSensors_ar1_aov <- format_anova_table(fallmsmodel_BetweenSensors_ar1_aov, "ControlSensorComparison", "FallDuration", "AR1Model") 
table_fallmsmodel_BetweenSensors_ar1_pairs_sensor <- format_emm_table(fallmsmodel_BetweenSensors_ar1_sensorcontrasts, "ControlSensorComparison", "FallDuration", "Pairs", "Sensor")

### Fall Slope -----
#### Mixed effects model
fallslopemodel_BetweenSensors <- glmmTMB(scale(fallslope) ~ Virus + BinF + (1|SubjectID), # random intercept by subject
                                                 data = subset(subjectmeansbinSonly, !is.na(fallslope)), family = stats::gaussian,
                                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

fallslopemodel_BetweenSensors_aov <- Anova(fallslopemodel_BetweenSensors, type=2)
fallslopemodel_BetweenSensors_aov
summary(fallslopemodel_BetweenSensors)

fallslopemodel_BetweenSensors_r2 <- as.data.frame(r.squaredGLMM(fallslopemodel_BetweenSensors)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
fallslopemodel_BetweenSensors_r2

# Follow up tests - sensor pairs
fallslopemodel_BetweenSensors_sensorpairs <- emmeans(fallslopemodel_BetweenSensors, pairwise~Virus, adjust="none") # Pairs
fallslopemodel_BetweenSensors_sensorcontrasts <- as.data.frame(fallslopemodel_BetweenSensors_sensorpairs$contrasts) # Extract emmeans
fallslopemodel_BetweenSensors_sensoremmeans <- fallslopemodel_BetweenSensors_sensorpairs$emmeans # Extract emmeans
fallslopemodel_BetweenSensors_sensorpairs_es <- as.data.frame(eff_size(fallslopemodel_BetweenSensors_sensoremmeans, sigma = sigma(fallslopemodel_BetweenSensors), edf = df.residual(fallslopemodel_BetweenSensors))) # Effect Size
fallslopemodel_BetweenSensors_sensorpairs_es <- fallslopemodel_BetweenSensors_sensorpairs_es[, c("contrast", "effect.size")]

fallslopemodel_BetweenSensors_sensorcontrasts$d <- fallslopemodel_BetweenSensors_sensorpairs_es$effect.size
fallslopemodel_BetweenSensors_sensorcontrasts

# Format Tables
table_fallslopemodel_BetweenSensors <- format_regsummary_table(fallslopemodel_BetweenSensors, "ControlSensorComparison", "FallSlope", "StandardModel")
table_fallslopemodel_BetweenSensors_aov <- format_anova_table(fallslopemodel_BetweenSensors_aov, "ControlSensorComparison", "FallSlope", "StandardModel")
table_fallslopemodel_BetweenSensors_pairs_sensor <- format_emm_table(fallslopemodel_BetweenSensors_sensorcontrasts, "ControlSensorComparison", "FallSlope", "Pairs", "Sensor")

#### Mixed effects AR1 model
fallslopemodel_BetweenSensors_ar1 <- glmmTMB(scale(fallslope) ~ Virus + BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                                     data = subset(subjectmeansbinSonly, !is.na(fallslope)), family = stats::gaussian, 
                                                     control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

fallslopemodel_BetweenSensors_ar1_aov <- Anova(fallslopemodel_BetweenSensors_ar1, type=2)
fallslopemodel_BetweenSensors_ar1_aov
summary(fallslopemodel_BetweenSensors_ar1)

fallslopemodel_BetweenSensors_ar1_r2 <- as.data.frame(r.squaredGLMM(fallslopemodel_BetweenSensors_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
fallslopemodel_BetweenSensors_ar1_r2

fallslopemodel_BetweenSensors_ar1_AIC <- as.data.frame(AIC(fallslopemodel_BetweenSensors, fallslopemodel_BetweenSensors_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "fallslopemodel_BetweenSensors" ~ "StandardModel", ModelName == "fallslopemodel_BetweenSensors_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

fallslopemodel_BetweenSensors_ar1_AIC


# Follow up tests - sensor pairs
fallslopemodel_BetweenSensors_ar1_sensorpairs <- emmeans(fallslopemodel_BetweenSensors_ar1, pairwise~Virus, adjust="none") # Pairs
fallslopemodel_BetweenSensors_ar1_sensorcontrasts <- as.data.frame(fallslopemodel_BetweenSensors_ar1_sensorpairs$contrasts) # Extract emmeans
fallslopemodel_BetweenSensors_ar1_sensoremmeans <- fallslopemodel_BetweenSensors_ar1_sensorpairs$emmeans # Extract emmeans
fallslopemodel_BetweenSensors_ar1_sensorpairs_es <- as.data.frame(eff_size(fallslopemodel_BetweenSensors_ar1_sensoremmeans, sigma = sigma(fallslopemodel_BetweenSensors_ar1), edf = df.residual(fallslopemodel_BetweenSensors_ar1))) # Effect Size
fallslopemodel_BetweenSensors_ar1_sensorpairs_es <- fallslopemodel_BetweenSensors_ar1_sensorpairs_es[, c("contrast", "effect.size")]

fallslopemodel_BetweenSensors_ar1_sensorcontrasts$d <- fallslopemodel_BetweenSensors_ar1_sensorpairs_es$effect.size
fallslopemodel_BetweenSensors_ar1_sensorcontrasts

# Format Tables
table_fallslopemodel_BetweenSensors_ar1 <- format_regsummary_table(fallslopemodel_BetweenSensors_ar1, 'ControlSensorComparison', "FallSlope", "AR1Model")
table_fallslopemodel_BetweenSensors_ar1_aov <- format_anova_table(fallslopemodel_BetweenSensors_ar1_aov, "ControlSensorComparison", "FallSlope", "AR1Model") 
table_fallslopemodel_BetweenSensors_ar1_pairs_sensor <- format_emm_table(fallslopemodel_BetweenSensors_ar1_sensorcontrasts, "ControlSensorComparison", "FallSlope", "Pairs", "Sensor")



### Compound Events -----
#### Mixed effects model
ncompoundeventsmodel_BetweenSensors <- glmmTMB(ncompoundevents ~ Virus + BinF + (1|SubjectID), # random intercept by subject
                                                       data = subjectmeansbinSonly, family = stats::gaussian,
                                                       control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ncompoundeventsmodel_BetweenSensors_aov <- Anova(ncompoundeventsmodel_BetweenSensors, type=2)
ncompoundeventsmodel_BetweenSensors_aov
summary(ncompoundeventsmodel_BetweenSensors)

ncompoundeventsmodel_BetweenSensors_r2 <- as.data.frame(r.squaredGLMM(ncompoundeventsmodel_BetweenSensors)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ncompoundeventsmodel_BetweenSensors_r2

# Follow up tests - sensor pairs
ncompoundeventsmodel_BetweenSensors_sensorpairs <- emmeans(ncompoundeventsmodel_BetweenSensors, pairwise~Virus, adjust="none") # Pairs
ncompoundeventsmodel_BetweenSensors_sensorcontrasts <- as.data.frame(ncompoundeventsmodel_BetweenSensors_sensorpairs$contrasts) # Extract emmeans
ncompoundeventsmodel_BetweenSensors_sensoremmeans <- ncompoundeventsmodel_BetweenSensors_sensorpairs$emmeans # Extract emmeans
ncompoundeventsmodel_BetweenSensors_sensorpairs_es <- as.data.frame(eff_size(ncompoundeventsmodel_BetweenSensors_sensoremmeans, sigma = sigma(ncompoundeventsmodel_BetweenSensors), edf = df.residual(ncompoundeventsmodel_BetweenSensors))) # Effect Size
ncompoundeventsmodel_BetweenSensors_sensorpairs_es <- ncompoundeventsmodel_BetweenSensors_sensorpairs_es[, c("contrast", "effect.size")]

ncompoundeventsmodel_BetweenSensors_sensorcontrasts$d <- ncompoundeventsmodel_BetweenSensors_sensorpairs_es$effect.size
ncompoundeventsmodel_BetweenSensors_sensorcontrasts

# Format Tables
table_ncompoundeventsmodel_BetweenSensors <- format_regsummary_table(ncompoundeventsmodel_BetweenSensors, "ControlSensorComparison", "nCompoundEvents", "StandardModel")
table_ncompoundeventsmodel_BetweenSensors_aov <- format_anova_table(ncompoundeventsmodel_BetweenSensors_aov, "ControlSensorComparison", "nCompoundEvents", "StandardModel")
table_ncompoundeventsmodel_BetweenSensors_pairs_sensor <- format_emm_table(ncompoundeventsmodel_BetweenSensors_sensorcontrasts, "ControlSensorComparison", "nCompoundEvents", "Pairs", "Sensor")

#### Mixed effects AR1 model
ncompoundeventsmodel_BetweenSensors_ar1 <- glmmTMB(ncompoundevents ~ Virus + BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                                           data = subjectmeansbinSonly, family = stats::gaussian, 
                                                           control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ncompoundeventsmodel_BetweenSensors_ar1_aov <- Anova(ncompoundeventsmodel_BetweenSensors_ar1, type=2)
ncompoundeventsmodel_BetweenSensors_ar1_aov
summary(ncompoundeventsmodel_BetweenSensors_ar1)

ncompoundeventsmodel_BetweenSensors_ar1_r2 <- as.data.frame(r.squaredGLMM(ncompoundeventsmodel_BetweenSensors_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ncompoundeventsmodel_BetweenSensors_ar1_r2

ncompoundeventsmodel_BetweenSensors_ar1_AIC <- as.data.frame(AIC(ncompoundeventsmodel_BetweenSensors, ncompoundeventsmodel_BetweenSensors_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "ncompoundeventsmodel_BetweenSensors" ~ "StandardModel", ModelName == "ncompoundeventsmodel_BetweenSensors_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

ncompoundeventsmodel_BetweenSensors_ar1_AIC


# Follow up tests - sensor pairs
ncompoundeventsmodel_BetweenSensors_ar1_sensorpairs <- emmeans(ncompoundeventsmodel_BetweenSensors_ar1, pairwise~Virus, adjust="none") # Pairs
ncompoundeventsmodel_BetweenSensors_ar1_sensorcontrasts <- as.data.frame(ncompoundeventsmodel_BetweenSensors_ar1_sensorpairs$contrasts) # Extract emmeans
ncompoundeventsmodel_BetweenSensors_ar1_sensoremmeans <- ncompoundeventsmodel_BetweenSensors_ar1_sensorpairs$emmeans # Extract emmeans
ncompoundeventsmodel_BetweenSensors_ar1_sensorpairs_es <- as.data.frame(eff_size(ncompoundeventsmodel_BetweenSensors_ar1_sensoremmeans, sigma = sigma(ncompoundeventsmodel_BetweenSensors_ar1), edf = df.residual(ncompoundeventsmodel_BetweenSensors_ar1))) # Effect Size
ncompoundeventsmodel_BetweenSensors_ar1_sensorpairs_es <- ncompoundeventsmodel_BetweenSensors_ar1_sensorpairs_es[, c("contrast", "effect.size")]

ncompoundeventsmodel_BetweenSensors_ar1_sensorcontrasts$d <- ncompoundeventsmodel_BetweenSensors_ar1_sensorpairs_es$effect.size
ncompoundeventsmodel_BetweenSensors_ar1_sensorcontrasts

# Format Tables
table_ncompoundeventsmodel_BetweenSensors_ar1 <- format_regsummary_table(ncompoundeventsmodel_BetweenSensors_ar1, 'ControlSensorComparison', "nCompoundEvents", "AR1Model")
table_ncompoundeventsmodel_BetweenSensors_ar1_aov <- format_anova_table(ncompoundeventsmodel_BetweenSensors_ar1_aov, "ControlSensorComparison", "nCompoundEvents", "AR1Model") 
table_ncompoundeventsmodel_BetweenSensors_ar1_pairs_sensor <- format_emm_table(ncompoundeventsmodel_BetweenSensors_ar1_sensorcontrasts, "ControlSensorComparison", "nCompoundEvents", "Pairs", "Sensor")


## EXPORT STATS -----
### Prepare Means -----
# Frequency
treatmeanstable_freq_BLSensorComp <- salinemeans %>%
  select(FiberPlacement, Virus, 
         pkspermin_mean, pkspermin_sd, pkspermin_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Amplitude
treatmeanstable_amp_BLSensorComp <- salinemeans %>%
  select(FiberPlacement, Virus, 
         amp_mean, amp_sd, amp_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# AUC window
treatmeanstable_AUCwindow_BLSensorComp <- salinemeans %>%
  select(FiberPlacement, Virus, 
         AUCwindow_mean, AUCwindow_sd, AUCwindow_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))


# Rise ms
treatmeanstable_risems_BLSensorComp <- salinemeans %>%
  select(FiberPlacement, Virus, 
         risems_mean, risems_sd, risems_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Rise slope
treatmeanstable_riseslope_BLSensorComp <- salinemeans %>%
  select(FiberPlacement, Virus, 
         riseslope_mean, riseslope_sd, riseslope_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))

# Fall ms
treatmeanstable_fallms_BLSensorComp <- salinemeans %>%
  select(FiberPlacement, Virus, 
         fallms_mean, fallms_sd, fallms_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Fall slope
treatmeanstable_fallslope_BLSensorComp <- salinemeans %>%
  select(FiberPlacement, Virus, 
         fallslope_mean, fallslope_sd, fallslope_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 6)))


### Frequency Export -----
wb_BLSensorComp_freq <- createWorkbook()
wb_BLSensorComp_freq_ar1 <- createWorkbook()

## Main results for BLSensorComp
# Standard Models
write_nice_sheet(wb_BLSensorComp_freq, "MeansRaw", treatmeanstable_freq_BLSensorComp)
write_nice_sheet(wb_BLSensorComp_freq, "ModelSummary", table_freqmodel_BetweenSensors)
write_nice_sheet(wb_BLSensorComp_freq, "ModelAnova", table_freqmodel_BetweenSensors_aov)
write_nice_sheet(wb_BLSensorComp_freq, "ModelR2", freqmodel_BetweenSensors_r2)
write_nice_sheet(wb_BLSensorComp_freq, "ModelAIC", freqmodel_BetweenSensors_ar1_AIC)
write_nice_sheet(wb_BLSensorComp_freq, "Pairs_Sensor", table_freqmodel_BetweenSensors_pairs_sensor)

saveWorkbook(wb_BLSensorComp_freq, paste(analysispath,"BLSensorComp_Frequency_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_BLSensorComp_freq_ar1, "MeansRaw", treatmeanstable_freq_BLSensorComp)
write_nice_sheet(wb_BLSensorComp_freq_ar1, "ModelSummary", table_freqmodel_BetweenSensors_ar1)
write_nice_sheet(wb_BLSensorComp_freq_ar1, "ModelAnova", table_freqmodel_BetweenSensors_ar1_aov)
write_nice_sheet(wb_BLSensorComp_freq_ar1, "ModelR2", freqmodel_BetweenSensors_ar1_r2)
write_nice_sheet(wb_BLSensorComp_freq_ar1, "ModelAIC", freqmodel_BetweenSensors_ar1_AIC)
write_nice_sheet(wb_BLSensorComp_freq_ar1, "Pairs_Sensor", table_freqmodel_BetweenSensors_ar1_pairs_sensor)

saveWorkbook(wb_BLSensorComp_freq_ar1, paste(analysispath,"BLSensorComp_Frequency_ar1models_results.xlsx",sep=''), overwrite = TRUE)

### Amplitude Export -----
wb_BLSensorComp_amp <- createWorkbook()

## Main results for BLSensorComp
# Standard Models
write_nice_sheet(wb_BLSensorComp_amp, "MeansRaw", treatmeanstable_amp_BLSensorComp)
write_nice_sheet(wb_BLSensorComp_amp, "ModelSummary", table_ampmodel_BetweenSensors)
write_nice_sheet(wb_BLSensorComp_amp, "ModelAnova", table_ampmodel_BetweenSensors_aov)
write_nice_sheet(wb_BLSensorComp_amp, "ModelR2", ampmodel_BetweenSensors_r2)
write_nice_sheet(wb_BLSensorComp_amp, "Pairs_Sensor", table_ampmodel_BetweenSensors_pairs_sensor)

saveWorkbook(wb_BLSensorComp_amp, paste(analysispath,"BLSensorComp_Amplitude_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

### AUCwindow Export -----
wb_BLSensorComp_AUCwindow <- createWorkbook()

## Main results for BLSensorComp
# Standard Models
write_nice_sheet(wb_BLSensorComp_AUCwindow, "MeansRaw", treatmeanstable_AUCwindow_BLSensorComp)
write_nice_sheet(wb_BLSensorComp_AUCwindow, "ModelSummary", table_AUCwindowmodel_BetweenSensors)
write_nice_sheet(wb_BLSensorComp_AUCwindow, "ModelAnova", table_AUCwindowmodel_BetweenSensors_aov)
write_nice_sheet(wb_BLSensorComp_AUCwindow, "ModelR2", AUCwindowmodel_BetweenSensors_r2)
write_nice_sheet(wb_BLSensorComp_AUCwindow, "Pairs_Sensor", table_AUCwindowmodel_BetweenSensors_pairs_sensor)

saveWorkbook(wb_BLSensorComp_AUCwindow, paste(analysispath,"BLSensorComp_AUCwindow_standardmodels_results.xlsx",sep=''), overwrite = TRUE)


### RiseDuration Export -----
wb_BLSensorComp_risems <- createWorkbook()
wb_BLSensorComp_risems_ar1 <- createWorkbook()

## Main results for BLSensorComp
# Standard Models
write_nice_sheet(wb_BLSensorComp_risems, "MeansRaw", treatmeanstable_risems_BLSensorComp)
write_nice_sheet(wb_BLSensorComp_risems, "ModelSummary", table_risemsmodel_BetweenSensors)
write_nice_sheet(wb_BLSensorComp_risems, "ModelAnova", table_risemsmodel_BetweenSensors_aov)
write_nice_sheet(wb_BLSensorComp_risems, "ModelR2", risemsmodel_BetweenSensors_r2)
write_nice_sheet(wb_BLSensorComp_risems, "Pairs_Sensor", table_risemsmodel_BetweenSensors_pairs_sensor)

saveWorkbook(wb_BLSensorComp_risems, paste(analysispath,"BLSensorComp_RiseDuration_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_BLSensorComp_risems_ar1, "MeansRaw", treatmeanstable_risems_BLSensorComp)
write_nice_sheet(wb_BLSensorComp_risems_ar1, "ModelSummary", table_risemsmodel_BetweenSensors_ar1)
write_nice_sheet(wb_BLSensorComp_risems_ar1, "ModelAnova", table_risemsmodel_BetweenSensors_ar1_aov)
write_nice_sheet(wb_BLSensorComp_risems_ar1, "ModelR2", risemsmodel_BetweenSensors_ar1_r2)
write_nice_sheet(wb_BLSensorComp_risems_ar1, "ModelAIC", risemsmodel_BetweenSensors_ar1_AIC)
write_nice_sheet(wb_BLSensorComp_risems_ar1, "Pairs_Sensor", table_risemsmodel_BetweenSensors_ar1_pairs_sensor)

saveWorkbook(wb_BLSensorComp_risems_ar1, paste(analysispath,"BLSensorComp_RiseDuration_ar1models_results.xlsx",sep=''), overwrite = TRUE)

### RiseSlope Export -----
wb_BLSensorComp_riseslope <- createWorkbook()

## Main results for BLSensorComp
# Standard Models
write_nice_sheet(wb_BLSensorComp_riseslope, "MeansRaw", treatmeanstable_riseslope_BLSensorComp)
write_nice_sheet(wb_BLSensorComp_riseslope, "ModelSummary", table_riseslopemodel_BetweenSensors)
write_nice_sheet(wb_BLSensorComp_riseslope, "ModelAnova", table_riseslopemodel_BetweenSensors_aov)
write_nice_sheet(wb_BLSensorComp_riseslope, "ModelR2", riseslopemodel_BetweenSensors_r2)
write_nice_sheet(wb_BLSensorComp_riseslope, "Pairs_Sensor", table_riseslopemodel_BetweenSensors_pairs_sensor)

saveWorkbook(wb_BLSensorComp_riseslope, paste(analysispath,"BLSensorComp_RiseSlope_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

### FallDuration Export -----
wb_BLSensorComp_fallms <- createWorkbook()

## Main results for BLSensorComp
# Standard Models
write_nice_sheet(wb_BLSensorComp_fallms, "MeansRaw", treatmeanstable_fallms_BLSensorComp)
write_nice_sheet(wb_BLSensorComp_fallms, "ModelSummary", table_fallmsmodel_BetweenSensors)
write_nice_sheet(wb_BLSensorComp_fallms, "ModelAnova", table_fallmsmodel_BetweenSensors_aov)
write_nice_sheet(wb_BLSensorComp_fallms, "ModelR2", fallmsmodel_BetweenSensors_r2)
write_nice_sheet(wb_BLSensorComp_fallms, "Pairs_Sensor", table_fallmsmodel_BetweenSensors_pairs_sensor)

saveWorkbook(wb_BLSensorComp_fallms, paste(analysispath,"BLSensorComp_FallDuration_standardmodels_results.xlsx",sep=''), overwrite = TRUE)


### FallSlope Export -----
wb_BLSensorComp_fallslope <- createWorkbook()

## Main results for BLSensorComp
# Standard Models
write_nice_sheet(wb_BLSensorComp_fallslope, "MeansRaw", treatmeanstable_fallslope_BLSensorComp)
write_nice_sheet(wb_BLSensorComp_fallslope, "ModelSummary", table_fallslopemodel_BetweenSensors)
write_nice_sheet(wb_BLSensorComp_fallslope, "ModelAnova", table_fallslopemodel_BetweenSensors_aov)
write_nice_sheet(wb_BLSensorComp_fallslope, "ModelR2", fallslopemodel_BetweenSensors_r2)
write_nice_sheet(wb_BLSensorComp_fallslope, "Pairs_Sensor", table_fallslopemodel_BetweenSensors_pairs_sensor)

saveWorkbook(wb_BLSensorComp_fallslope, paste(analysispath,"BLSensorComp_FallSlope_standardmodels_results.xlsx",sep=''), overwrite = TRUE)



## PLOTS -----
### Frequency -----
sensorcontrolcomp_freqbar_min <- 0
sensorcontrolcomp_freqbar_max <- 8
sensorcontrolcomp_freqbar_ticks <- 2
sensorcontrolcomp_freqbar_breaks <- seq(sensorcontrolcomp_freqbar_min, sensorcontrolcomp_freqbar_max, sensorcontrolcomp_freqbar_ticks)

# Bar Graph
sensorcontrolcomp_freqbar <- ggplot(salinemeans,aes(x=Virus, y=pkspermin_mean, fill=Virus)) + 
  geom_col(linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=salinemeans, show.legend=FALSE,
                aes(x=Virus, ymin=pkspermin_mean-pkspermin_se, ymax=pkspermin_mean+pkspermin_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  geom_jitter(data=salinesubjectmeans, aes(x=Virus, y=pkspermin, fill=Virus, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = 0), show.legend=FALSE) +
  scale_shape_manual(values = sexshapes) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcontrolcomp_freqbar_min, sensorcontrolcomp_freqbar_max), 
                     breaks=sensorcontrolcomp_freqbar_breaks) +
  xlab("Sensor") +
  ylab("Frequency (n/min)") + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(7.6),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(1.96), y.position=c(6.8),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.04), xmax = c(3), y.position=c(6.8),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcontrolcomp_freqbar

### Amplitude -----
sensorcontrolcomp_ampbar_min <- 0
sensorcontrolcomp_ampbar_max <- 1.8
sensorcontrolcomp_ampbar_ticks <- .4
sensorcontrolcomp_ampbar_breaks <- seq(sensorcontrolcomp_ampbar_min, sensorcontrolcomp_ampbar_max, sensorcontrolcomp_ampbar_ticks)
sensorcontrolcomp_ampbar_labels <- sensorcontrolcomp_ampbar_breaks + threshold

# Bar Graph
sensorcontrolcomp_ampbar <- ggplot(salinemeans,aes(x=Virus, y=ampscaled_mean, fill=Virus)) + 
  geom_col(linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=salinemeans, show.legend=FALSE,
                aes(x=Virus, ymin=ampscaled_mean-amp_se, ymax=ampscaled_mean+amp_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  geom_jitter(data=salinesubjectmeans, aes(x=Virus, y=ampscaled, fill=Virus, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = 0), show.legend=FALSE) +
  scale_shape_manual(values = sexshapes) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcontrolcomp_ampbar_min, sensorcontrolcomp_ampbar_max), 
                     breaks=sensorcontrolcomp_ampbar_breaks, labels=sensorcontrolcomp_ampbar_labels) +
  xlab("Sensor") +
  ylab("Amplitude (Z Score)") + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(1.72),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(1.96), y.position=c(1.58),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.04), xmax = c(3), y.position=c(1.58),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcontrolcomp_ampbar

### AUCwindow -----
sensorcontrolcomp_AUCwindowbar_min <- 0
sensorcontrolcomp_AUCwindowbar_max <- 3.2
sensorcontrolcomp_AUCwindowbar_ticks <- .8
sensorcontrolcomp_AUCwindowbar_breaks <- seq(sensorcontrolcomp_AUCwindowbar_min, sensorcontrolcomp_AUCwindowbar_max, sensorcontrolcomp_AUCwindowbar_ticks)

# Bar Graph
sensorcontrolcomp_AUCwindowbar <- ggplot(salinemeans,aes(x=Virus, y=AUCwindow_mean, fill=Virus)) + 
  geom_col(linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=salinemeans, show.legend=FALSE,
                aes(x=Virus, ymin=AUCwindow_mean-AUCwindow_se, ymax=AUCwindow_mean+AUCwindow_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  geom_jitter(data=salinesubjectmeans, aes(x=Virus, y=AUCwindow, fill=Virus, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = 0), show.legend=FALSE) +
  scale_shape_manual(values = sexshapes) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcontrolcomp_AUCwindowbar_min, sensorcontrolcomp_AUCwindowbar_max), 
                     breaks=sensorcontrolcomp_AUCwindowbar_breaks) +
  xlab("Sensor") +
  ylab("Transient AUC") + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(3.1),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(1.96), y.position=c(2.9),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.04), xmax = c(3), y.position=c(2.9),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcontrolcomp_AUCwindowbar

### Rise Duration -----
sensorcontrolcomp_risemsbar_min <- 0
sensorcontrolcomp_risemsbar_max <- 250
sensorcontrolcomp_risemsbar_ticks <- 60
sensorcontrolcomp_risemsbar_breaks <- seq(sensorcontrolcomp_risemsbar_min, sensorcontrolcomp_risemsbar_max, sensorcontrolcomp_risemsbar_ticks)

# Bar Graph
sensorcontrolcomp_risemsbar <- ggplot(salinemeans,aes(x=Virus, y=risems_mean, fill=Virus)) + 
  geom_col(linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=salinemeans, show.legend=FALSE,
                aes(x=Virus, ymin=risems_mean-risems_se, ymax=risems_mean+risems_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  geom_jitter(data=salinesubjectmeans, aes(x=Virus, y=risems, fill=Virus, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = 0), show.legend=FALSE) +
  scale_shape_manual(values = sexshapes) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcontrolcomp_risemsbar_min, sensorcontrolcomp_risemsbar_max), 
                     breaks=sensorcontrolcomp_risemsbar_breaks) +
  xlab("Sensor") +
  ylab("Rise Duration (ms)") + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(1), xmax = c(1.96), y.position=c(240),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.04), xmax = c(3), y.position=c(240),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcontrolcomp_risemsbar

### Rise Slope -----
sensorcontrolcomp_riseslopebar_min <- 0
sensorcontrolcomp_riseslopebar_max <- .016
sensorcontrolcomp_riseslopebar_ticks <- .004
sensorcontrolcomp_riseslopebar_breaks <- seq(sensorcontrolcomp_riseslopebar_min, sensorcontrolcomp_riseslopebar_max, sensorcontrolcomp_riseslopebar_ticks)

# Bar Graph
sensorcontrolcomp_riseslopebar <- ggplot(salinemeans,aes(x=Virus, y=riseslope_mean, fill=Virus)) + 
  geom_col(linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=salinemeans, show.legend=FALSE,
                aes(x=Virus, ymin=riseslope_mean-riseslope_se, ymax=riseslope_mean+riseslope_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  geom_jitter(data=salinesubjectmeans, aes(x=Virus, y=riseslope, fill=Virus, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = 0), show.legend=FALSE) +
  scale_shape_manual(values = sexshapes) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcontrolcomp_riseslopebar_min, sensorcontrolcomp_riseslopebar_max), 
                     breaks=sensorcontrolcomp_riseslopebar_breaks) +
  xlab("Sensor") +
  ylab("Rise Slope") + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(.0155),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1), xmax = c(1.96), y.position=c(.0142),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.04), xmax = c(3), y.position=c(.0142),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcontrolcomp_riseslopebar


### Fall Duration -----
sensorcontrolcomp_fallmsbar_min <- 0
sensorcontrolcomp_fallmsbar_max <- 1200
sensorcontrolcomp_fallmsbar_ticks <- 300
sensorcontrolcomp_fallmsbar_breaks <- seq(sensorcontrolcomp_fallmsbar_min, sensorcontrolcomp_fallmsbar_max, sensorcontrolcomp_fallmsbar_ticks)

# Bar Graph
sensorcontrolcomp_fallmsbar <- ggplot(salinemeans,aes(x=Virus, y=fallms_mean, fill=Virus)) + 
  geom_col(linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=salinemeans, show.legend=FALSE,
                aes(x=Virus, ymin=fallms_mean-fallms_se, ymax=fallms_mean+fallms_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  geom_jitter(data=salinesubjectmeans, aes(x=Virus, y=fallms, fill=Virus, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = 0), show.legend=FALSE) +
  scale_shape_manual(values = sexshapes) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcontrolcomp_fallmsbar_min, sensorcontrolcomp_fallmsbar_max), 
                     breaks=sensorcontrolcomp_fallmsbar_breaks) +
  xlab("Sensor") +
  ylab("Fall Duration (ms)") + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(1150),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.04), xmax = c(3), y.position=c(1050),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcontrolcomp_fallmsbar

### Fall Slope -----
sensorcontrolcomp_fallslopebar_min <- -.016
sensorcontrolcomp_fallslopebar_max <- 0
sensorcontrolcomp_fallslopebar_ticks <- .004
sensorcontrolcomp_fallslopebar_breaks <- seq(sensorcontrolcomp_fallslopebar_min, sensorcontrolcomp_fallslopebar_max, sensorcontrolcomp_fallslopebar_ticks)

# Bar Graph
sensorcontrolcomp_fallslopebar <- ggplot(salinemeans,aes(x=Virus, y=fallslope_mean, fill=Virus)) + 
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=myaxislinewidth) +
  geom_col(linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=salinemeans, show.legend=FALSE,
                aes(x=Virus, ymin=fallslope_mean-fallslope_se, ymax=fallslope_mean+fallslope_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  geom_jitter(data=salinesubjectmeans, aes(x=Virus, y=fallslope, fill=Virus, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = 0), show.legend=FALSE) +
  scale_shape_manual(values = sexshapes) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcontrolcomp_fallslopebar_min, sensorcontrolcomp_fallslopebar_max), 
                     breaks=sensorcontrolcomp_fallslopebar_breaks) +
  xlab("Sensor") +
  ylab("Fall Slope") + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(1), xmax = c(1.96), y.position=c(-.014),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.04), xmax = c(3), y.position=c(-.014),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcontrolcomp_fallslopebar


### N Compound Events -----
sensorcontrolcomp_ncompoundeventsbar_min <- 0
sensorcontrolcomp_ncompoundeventsbar_max <- 9
sensorcontrolcomp_ncompoundeventsbar_ticks <- 3
sensorcontrolcomp_ncompoundeventsbar_breaks <- seq(sensorcontrolcomp_ncompoundeventsbar_min, sensorcontrolcomp_ncompoundeventsbar_max, sensorcontrolcomp_ncompoundeventsbar_ticks)

# Bar Graph
sensorcontrolcomp_ncompoundeventsbar <- ggplot(salinemeans,aes(x=Virus, y=ncompoundevents_summean, fill=Virus)) + 
  geom_col(linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=salinemeans, show.legend=FALSE,
                aes(x=Virus, ymin=ncompoundevents_summean-ncompoundevents_sumse, ymax=ncompoundevents_summean+ncompoundevents_sumse), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  geom_jitter(data=salinesubjectmeans, aes(x=Virus, y=ncompoundevents_sum, fill=Virus, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = 0), show.legend=FALSE) +
  scale_shape_manual(values = sexshapes) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcontrolcomp_ncompoundeventsbar_min, sensorcontrolcomp_ncompoundeventsbar_max), 
                     breaks=sensorcontrolcomp_ncompoundeventsbar_breaks) +
  xlab("Sensor") +
  ylab("Compound Events (n)") + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(1), xmax = c(1.96), y.position=c(8),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2.04), xmax = c(3), y.position=c(8),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcontrolcomp_ncompoundeventsbar

## EXPORT FIGURE 1 -----
# BAR AND REG GRAPHS
emptypanel <- ggplot() + mytheme

suppfigure_ControlTransients_panels <- plot_grid(emptypanel, sensorcontrolcomp_freqbar, sensorcontrolcomp_ampbar, sensorcontrolcomp_AUCwindowbar,
                                                 sensorcontrolcomp_risemsbar, sensorcontrolcomp_riseslopebar, sensorcontrolcomp_fallmsbar, sensorcontrolcomp_fallslopebar,
                                                 ncol = 4, nrow=2, align = "hv",  # align horizontally & vertically
                                                 axis  = "tblr", hjust=-.5)

suppfigure_ControlTransients_panels

ggsave(filename = c("SupplementaryFigure_ControlTransientsSensorComparison.pdf"), path = figurepath, plot = suppfigure_ControlTransients_panels, device = "pdf", 
       width = figpanelwidthin, height = figpanelheightinbarline, dpi=figdpi, units="in")


# FIGURE 2 - VTA GCaMP6f TRANSIENTS -----
## STATS ------
### Frequency -----
#### Frequency Model -----
#### Mixed effects model
freqmodel_GCaMP6f <- glmmTMB(pkspermin ~ InjType*Dose*BinF + (1|SubjectID), # random intercept by subject
                             data = subset(subjectmeansbin, Virus=='GCaMP6f' & CondNum == 1), family = stats::gaussian,
                             control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqmodel_GCaMP6f_aov <- Anova(freqmodel_GCaMP6f, type=2)
freqmodel_GCaMP6f_aov
summary(freqmodel_GCaMP6f)

freqmodel_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(freqmodel_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqmodel_GCaMP6f_r2

# Follow up tests - joint
freqmodel_GCaMP6f_emmeans <-emmeans(freqmodel_GCaMP6f, pairwise~InjType|Dose, adjust="none") # Find means
freqmodel_GCaMP6f_jt_Dose <- joint_tests(freqmodel_GCaMP6f_emmeans, by = "Dose")
freqmodel_GCaMP6f_jt_InjType <- joint_tests(freqmodel_GCaMP6f_emmeans, by = "InjType")

freqmodel_GCaMP6f_jt_Dose
freqmodel_GCaMP6f_jt_InjType

# Follow up tests - treat pairs
freqmodel_GCaMP6f_treatpairs <- emmeans(freqmodel_GCaMP6f, pairwise~InjType|Dose, adjust="none") # Pairs
freqmodel_GCaMP6f_treatcontrasts <- as.data.frame(freqmodel_GCaMP6f_treatpairs$contrasts) # Extract emmeans
freqmodel_GCaMP6f_treatemmeans <- freqmodel_GCaMP6f_treatpairs$emmeans # Extract emmeans
freqmodel_GCaMP6f_treatpairs_es <- as.data.frame(eff_size(freqmodel_GCaMP6f_treatemmeans, sigma = sigma(freqmodel_GCaMP6f), edf = df.residual(freqmodel_GCaMP6f))) # Effect Size
freqmodel_GCaMP6f_treatpairs_es <- freqmodel_GCaMP6f_treatpairs_es[, c("Dose", "contrast", "effect.size")]

freqmodel_GCaMP6f_treatcontrasts$d <- freqmodel_GCaMP6f_treatpairs_es$effect.size
freqmodel_GCaMP6f_treatcontrasts

# Follow up tests - dose pairs
freqmodel_GCaMP6f_dosepairs <- emmeans(freqmodel_GCaMP6f, pairwise~Dose|InjType, adjust="none") # Pairs
freqmodel_GCaMP6f_dosecontrasts <- as.data.frame(freqmodel_GCaMP6f_dosepairs$contrasts) # Extract emmeans
freqmodel_GCaMP6f_doseemmeans <- freqmodel_GCaMP6f_dosepairs$emmeans # Extract emmeans
freqmodel_GCaMP6f_dosepairs_es <- as.data.frame(eff_size(freqmodel_GCaMP6f_doseemmeans, sigma = sigma(freqmodel_GCaMP6f), edf = df.residual(freqmodel_GCaMP6f))) # Effect Size
freqmodel_GCaMP6f_dosepairs_es <- freqmodel_GCaMP6f_dosepairs_es[, c("InjType", "contrast", "effect.size")]

freqmodel_GCaMP6f_dosecontrasts$d <- freqmodel_GCaMP6f_dosepairs_es$effect.size
freqmodel_GCaMP6f_dosecontrasts

# Format Tables
table_freqmodel_GCaMP6f <- format_regsummary_table(freqmodel_GCaMP6f, 'GCaMP6f', "Frequency", "StandardModel")
table_freqmodel_GCaMP6f_aov <- format_anova_table(freqmodel_GCaMP6f_aov, "GCaMP6f", "Frequency", "StandardModel")
table_freqmodel_GCaMP6f_joint_Dose <- format_emm_table(freqmodel_GCaMP6f_jt_Dose, "GCaMP6f", "Frequency", "Joint", "Dose")
table_freqmodel_GCaMP6f_joint_InjType <- format_emm_table(freqmodel_GCaMP6f_jt_InjType, "GCaMP6f", "Frequency", "Joint", "InjType")
table_freqmodel_GCaMP6f_pairs_Dose <- format_emm_table(freqmodel_GCaMP6f_dosecontrasts, "GCaMP6f", "Frequency", "Pairs", "Dose|InjType")
table_freqmodel_GCaMP6f_pairs_InjType <- format_emm_table(freqmodel_GCaMP6f_treatcontrasts, "GCaMP6f", "Frequency", "Pairs", "InjType|Dose")

#### Mixed effects AR1 model
freqmodel_GCaMP6f_ar1 <- glmmTMB(pkspermin ~ InjType*Dose*BinF + ar1(BinF + 0|SubjectID:InjType) + (1|SubjectID), # random intercept by subject
                                 data = subset(subjectmeansbin, Virus=='GCaMP6f' & CondNum == 1), family = stats::gaussian, 
                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqmodel_GCaMP6f_ar1_aov <- Anova(freqmodel_GCaMP6f_ar1, type=2)
freqmodel_GCaMP6f_ar1_aov
summary(freqmodel_GCaMP6f_ar1)

freqmodel_GCaMP6f_ar1_r2 <- as.data.frame(r.squaredGLMM(freqmodel_GCaMP6f_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqmodel_GCaMP6f_ar1_r2

freqmodel_GCaMP6f_ar1_AIC <- as.data.frame(AIC(freqmodel_GCaMP6f, freqmodel_GCaMP6f_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "freqmodel_GCaMP6f" ~ "StandardModel", ModelName == "freqmodel_GCaMP6f_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
           select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

freqmodel_GCaMP6f_ar1_AIC

# Follow up tests - joint
freqmodel_GCaMP6f_ar1_emmeans <-emmeans(freqmodel_GCaMP6f_ar1, pairwise~InjType|Dose, adjust="none") # Find means
freqmodel_GCaMP6f_ar1_jt_Dose <- joint_tests(freqmodel_GCaMP6f_ar1, by = "Dose")
freqmodel_GCaMP6f_ar1_jt_InjType <- joint_tests(freqmodel_GCaMP6f_ar1, by = "InjType")

freqmodel_GCaMP6f_ar1_jt_Dose
freqmodel_GCaMP6f_ar1_jt_InjType

# Follow up tests - treat pairs
freqmodel_GCaMP6f_ar1_treatpairs <- emmeans(freqmodel_GCaMP6f_ar1, pairwise~InjType|Dose, adjust="none") # Pairs
freqmodel_GCaMP6f_ar1_treatcontrasts <- as.data.frame(freqmodel_GCaMP6f_ar1_treatpairs$contrasts) # Extract emmeans
freqmodel_GCaMP6f_ar1_treatemmeans <- freqmodel_GCaMP6f_ar1_treatpairs$emmeans # Extract emmeans
freqmodel_GCaMP6f_ar1_treatpairs_es <- as.data.frame(eff_size(freqmodel_GCaMP6f_ar1_treatemmeans, sigma = sigma(freqmodel_GCaMP6f_ar1), edf = df.residual(freqmodel_GCaMP6f_ar1))) # Effect Size
freqmodel_GCaMP6f_ar1_treatpairs_es <- freqmodel_GCaMP6f_ar1_treatpairs_es[, c("Dose", "contrast", "effect.size")]

freqmodel_GCaMP6f_ar1_treatcontrasts$d <- freqmodel_GCaMP6f_ar1_treatpairs_es$effect.size
freqmodel_GCaMP6f_ar1_treatcontrasts

# Follow up tests - dose pairs
freqmodel_GCaMP6f_ar1_dosepairs <- emmeans(freqmodel_GCaMP6f_ar1, pairwise~Dose|InjType, adjust="none") # Pairs
freqmodel_GCaMP6f_ar1_dosecontrasts <- as.data.frame(freqmodel_GCaMP6f_ar1_dosepairs$contrasts) # Extract emmeans
freqmodel_GCaMP6f_ar1_doseemmeans <- freqmodel_GCaMP6f_ar1_dosepairs$emmeans # Extract emmeans
freqmodel_GCaMP6f_ar1_dosepairs_es <- as.data.frame(eff_size(freqmodel_GCaMP6f_ar1_doseemmeans, sigma = sigma(freqmodel_GCaMP6f_ar1), edf = df.residual(freqmodel_GCaMP6f_ar1))) # Effect Size
freqmodel_GCaMP6f_ar1_dosepairs_es <- freqmodel_GCaMP6f_ar1_dosepairs_es[, c("InjType", "contrast", "effect.size")]

freqmodel_GCaMP6f_ar1_dosecontrasts$d <- freqmodel_GCaMP6f_ar1_dosepairs_es$effect.size
freqmodel_GCaMP6f_ar1_dosecontrasts

# Format Tables
table_freqmodel_GCaMP6f_ar1 <- format_regsummary_table(freqmodel_GCaMP6f_ar1, 'GCaMP6f', "Frequency", "AR1Model")
table_freqmodel_GCaMP6f_ar1_aov <- format_anova_table(freqmodel_GCaMP6f_ar1_aov, "GCaMP6f", "Frequency", "AR1Model") 
table_freqmodel_GCaMP6f_ar1_joint_Dose <- format_emm_table(freqmodel_GCaMP6f_ar1_jt_Dose, "GCaMP6f", "Frequency", "Joint", "Dose")
table_freqmodel_GCaMP6f_ar1_joint_InjType <- format_emm_table(freqmodel_GCaMP6f_ar1_jt_InjType, "GCaMP6f", "Frequency", "Joint", "InjType")
table_freqmodel_GCaMP6f_ar1_pairs_Dose <- format_emm_table(freqmodel_GCaMP6f_ar1_dosecontrasts, "GCaMP6f", "Frequency", "Pairs", "Dose|InjType")
table_freqmodel_GCaMP6f_ar1_pairs_InjType <- format_emm_table(freqmodel_GCaMP6f_ar1_treatcontrasts, "GCaMP6f", "Frequency", "Pairs", "InjType|Dose")

#### Frequency Change Model -----
# Mixed effects model
freqCmodel_GCaMP6f <- glmmTMB(pksperminC ~ Dose*BinF + (1|SubjectID), # random intercept by subject
                              data = subset(subjectmeansbin, Virus=='GCaMP6f' & CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                              control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqCmodel_GCaMP6f_aov <- Anova(freqCmodel_GCaMP6f, type=2)
freqCmodel_GCaMP6f_aov
summary(freqCmodel_GCaMP6f)

freqCmodel_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(freqCmodel_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqCmodel_GCaMP6f_r2

# Follow up tests - joint
freqCmodel_GCaMP6f_emmeans <-emmeans(freqCmodel_GCaMP6f, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
freqCmodel_GCaMP6f_dosepairs <- emmeans(freqCmodel_GCaMP6f, pairwise~Dose, adjust="none") # Pairs
freqCmodel_GCaMP6f_dosecontrasts <- as.data.frame(freqCmodel_GCaMP6f_dosepairs$contrasts) # Extract emmeans
freqCmodel_GCaMP6f_doseemmeans <- freqCmodel_GCaMP6f_dosepairs$emmeans # Extract emmeans
freqCmodel_GCaMP6f_dosepairs_es <- as.data.frame(eff_size(freqCmodel_GCaMP6f_doseemmeans, sigma = sigma(freqCmodel_GCaMP6f), edf = df.residual(freqCmodel_GCaMP6f))) # Effect Size
freqCmodel_GCaMP6f_dosepairs_es <- freqCmodel_GCaMP6f_dosepairs_es[, c("contrast", "effect.size")]

freqCmodel_GCaMP6f_dosecontrasts$d <- freqCmodel_GCaMP6f_dosepairs_es$effect.size
freqCmodel_GCaMP6f_dosecontrasts

# Format Tables
table_freqCmodel_GCaMP6f <- format_regsummary_table(freqCmodel_GCaMP6f, 'GCaMP6f', "Frequency", "StandardModel")
table_freqCmodel_GCaMP6f_aov <- format_anova_table(freqCmodel_GCaMP6f_aov, "GCaMP6f", "Frequency", "StandardModel")
table_freqCmodel_GCaMP6f_pairs_Dose <- format_emm_table(freqCmodel_GCaMP6f_dosecontrasts, "GCaMP6f", "Frequency", "Pairs","Dose")

#### Mixed effects AR1 model
freqCmodel_GCaMP6f_ar1 <- glmmTMB(pksperminC ~ Dose*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                  data = subset(subjectmeansbin, Virus=='GCaMP6f' & CondNum == 1& InjType=='M'), family = stats::gaussian, 
                                  control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqCmodel_GCaMP6f_ar1_aov <- Anova(freqCmodel_GCaMP6f_ar1, type=2)
freqCmodel_GCaMP6f_ar1_aov
summary(freqCmodel_GCaMP6f_ar1)

freqCmodel_GCaMP6f_ar1_r2 <- as.data.frame(r.squaredGLMM(freqCmodel_GCaMP6f_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqCmodel_GCaMP6f_ar1_r2

freqCmodel_GCaMP6f_ar1_AIC <- as.data.frame(AIC(freqCmodel_GCaMP6f, freqCmodel_GCaMP6f_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "freqCmodel_GCaMP6f" ~ "StandardModel", ModelName == "freqCmodel_GCaMP6f_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

freqCmodel_GCaMP6f_ar1_AIC

# Follow up tests - joint
freqCmodel_GCaMP6f_ar1_emmeans <-emmeans(freqCmodel_GCaMP6f_ar1, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
freqCmodel_GCaMP6f_ar1_dosepairs <- emmeans(freqCmodel_GCaMP6f_ar1, pairwise~Dose, adjust="none") # Pairs
freqCmodel_GCaMP6f_ar1_dosecontrasts <- as.data.frame(freqCmodel_GCaMP6f_ar1_dosepairs$contrasts) # Extract emmeans
freqCmodel_GCaMP6f_ar1_doseemmeans <- freqCmodel_GCaMP6f_ar1_dosepairs$emmeans # Extract emmeans
freqCmodel_GCaMP6f_ar1_dosepairs_es <- as.data.frame(eff_size(freqCmodel_GCaMP6f_ar1_doseemmeans, sigma = sigma(freqCmodel_GCaMP6f_ar1), edf = df.residual(freqCmodel_GCaMP6f_ar1))) # Effect Size
freqCmodel_GCaMP6f_ar1_dosepairs_es <- freqCmodel_GCaMP6f_ar1_dosepairs_es[, c("contrast", "effect.size")]

freqCmodel_GCaMP6f_ar1_dosecontrasts$d <- freqCmodel_GCaMP6f_ar1_dosepairs_es$effect.size
freqCmodel_GCaMP6f_ar1_dosecontrasts

# Format Tables
table_freqCmodel_GCaMP6f_ar1 <- format_regsummary_table(freqCmodel_GCaMP6f_ar1, 'GCaMP6f', "Frequency", "AR1Model")
table_freqCmodel_GCaMP6f_ar1_aov <- format_anova_table(freqCmodel_GCaMP6f_ar1_aov, "GCaMP6f", "Frequency", "AR1Model") 
table_freqCmodel_GCaMP6f_ar1_pairs_Dose <- format_emm_table(freqCmodel_GCaMP6f_ar1_dosecontrasts, "GCaMP6f", "Frequency", "Pairs", "Dose")

#### Max Height Model -----
## Mixed effects model
MaxHeightmodel_freq_GCaMP6f <- glmmTMB(MaxHeight ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                       data = subjectrise_pksperminC_GCaMP6f, family = stats::gaussian)

MaxHeightmodel_freq_GCaMP6f_aov <- Anova(MaxHeightmodel_freq_GCaMP6f, type=2)
MaxHeightmodel_freq_GCaMP6f_aov
summary(MaxHeightmodel_freq_GCaMP6f)

MaxHeightmodel_freq_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(MaxHeightmodel_freq_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxHeightmodel_freq_GCaMP6f_r2

## Prepare fit data for plotting
MaxHeight_freq_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_GCaMP6f$DoseNum)))
MaxHeight_freq_GCaMP6f_regdata$DoseNum_c <- MaxHeight_freq_GCaMP6f_regdata$DoseNum - dosemean

MaxHeight_freq_GCaMP6f_regpred <- predict(MaxHeightmodel_freq_GCaMP6f, newdata = MaxHeight_freq_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxHeight_freq_GCaMP6f_regdata$fit <- MaxHeight_freq_GCaMP6f_regpred$fit
MaxHeight_freq_GCaMP6f_regdata$se <- MaxHeight_freq_GCaMP6f_regpred$se
MaxHeight_freq_GCaMP6f_regdata$lower <- MaxHeight_freq_GCaMP6f_regdata$fit - 1.96 * MaxHeight_freq_GCaMP6f_regdata$se
MaxHeight_freq_GCaMP6f_regdata$upper <- MaxHeight_freq_GCaMP6f_regdata$fit + 1.96 * MaxHeight_freq_GCaMP6f_regdata$se

MaxHeight_freq_GCaMP6f_regdata

# Format Tables
table_freqmodel_GCaMP6f_MaxHeight <- format_regsummary_table(MaxHeightmodel_freq_GCaMP6f, 'GCaMP6f', "Frequency", "MaxHeightModel")
table_freqmodel_GCaMP6f_MaxHeight_aov <- format_anova_table(MaxHeightmodel_freq_GCaMP6f_aov, "GCaMP6f", "Frequency", "MaxHeightModel") 
table_freqmodel_GCaMP6f_MaxHeight_reg <- format_reg_table(MaxHeight_freq_GCaMP6f_regdata, "GCaMP6f", "Frequency", "MaxHeightReg")

#### Max Duration Model -----
## Mixed effects model
MaxDurationmodel_freq_GCaMP6f <- glmmTMB(MaxDuration~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                         data = subjectrise_pksperminC_GCaMP6f, family = stats::gaussian)

MaxDurationmodel_freq_GCaMP6f_aov<- Anova(MaxDurationmodel_freq_GCaMP6f, type=2)
MaxDurationmodel_freq_GCaMP6f_aov
summary(MaxDurationmodel_freq_GCaMP6f)

MaxDurationmodel_freq_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(MaxDurationmodel_freq_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxDurationmodel_freq_GCaMP6f_r2

## Prepare fit data for plotting
MaxDuration_freq_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_GCaMP6f$DoseNum)))
MaxDuration_freq_GCaMP6f_regdata$DoseNum_c <- MaxDuration_freq_GCaMP6f_regdata$DoseNum - dosemean

MaxDuration_freq_GCaMP6f_regpred <- predict(MaxDurationmodel_freq_GCaMP6f, newdata = MaxDuration_freq_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxDuration_freq_GCaMP6f_regdata$fit <- MaxDuration_freq_GCaMP6f_regpred$fit
MaxDuration_freq_GCaMP6f_regdata$se <- MaxDuration_freq_GCaMP6f_regpred$se
MaxDuration_freq_GCaMP6f_regdata$lower <- MaxDuration_freq_GCaMP6f_regdata$fit - 1.96 * MaxDuration_freq_GCaMP6f_regdata$se
MaxDuration_freq_GCaMP6f_regdata$upper <- MaxDuration_freq_GCaMP6f_regdata$fit + 1.96 * MaxDuration_freq_GCaMP6f_regdata$se

MaxDuration_freq_GCaMP6f_regdata

# Format Tables
table_freqmodel_GCaMP6f_MaxDuration <- format_regsummary_table(MaxDurationmodel_freq_GCaMP6f, 'GCaMP6f', "Frequency", "MaxDurationModel")
table_freqmodel_GCaMP6f_MaxDuration_aov <- format_anova_table(MaxDurationmodel_freq_GCaMP6f_aov, "GCaMP6f", "Frequency", "MaxDurationModel") 
table_freqmodel_GCaMP6f_MaxDuration_reg <- format_reg_table(MaxDuration_freq_GCaMP6f_regdata, "GCaMP6f", "Frequency", "MaxDurationReg")

#### Rise Slope Model -----
## Mixed effects model
MaxSlopeLMmodel_freq_GCaMP6f <- glmmTMB(MaxSlopeLM ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                        data = subjectrise_pksperminC_GCaMP6f, family = stats::gaussian)

MaxSlopeLMmodel_freq_GCaMP6f_aov <- Anova(MaxSlopeLMmodel_freq_GCaMP6f, type=2)
MaxSlopeLMmodel_freq_GCaMP6f_aov
summary(MaxSlopeLMmodel_freq_GCaMP6f)

MaxSlopeLMmodel_freq_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(MaxSlopeLMmodel_freq_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxSlopeLMmodel_freq_GCaMP6f_r2

## Prepare fit data for plotting
MaxSlopeLM_freq_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_GCaMP6f$DoseNum)))
MaxSlopeLM_freq_GCaMP6f_regdata$DoseNum_c <- MaxSlopeLM_freq_GCaMP6f_regdata$DoseNum - dosemean

MaxSlopeLM_freq_GCaMP6f_regpred <- predict(MaxSlopeLMmodel_freq_GCaMP6f, newdata = MaxSlopeLM_freq_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxSlopeLM_freq_GCaMP6f_regdata$fit <- MaxSlopeLM_freq_GCaMP6f_regpred$fit
MaxSlopeLM_freq_GCaMP6f_regdata$se <- MaxSlopeLM_freq_GCaMP6f_regpred$se
MaxSlopeLM_freq_GCaMP6f_regdata$lower <- MaxSlopeLM_freq_GCaMP6f_regdata$fit - 1.96 * MaxSlopeLM_freq_GCaMP6f_regdata$se
MaxSlopeLM_freq_GCaMP6f_regdata$upper <- MaxSlopeLM_freq_GCaMP6f_regdata$fit + 1.96 * MaxSlopeLM_freq_GCaMP6f_regdata$se

MaxSlopeLM_freq_GCaMP6f_regdata

# Format Tables
table_freqmodel_GCaMP6f_MaxSlopeLM <- format_regsummary_table(MaxSlopeLMmodel_freq_GCaMP6f, 'GCaMP6f', "Frequency", "MaxSlopeLMModel")
table_freqmodel_GCaMP6f_MaxSlopeLM_aov <- format_anova_table(MaxSlopeLMmodel_freq_GCaMP6f_aov, "GCaMP6f", "Frequency", "MaxSlopeLMModel") 
table_freqmodel_GCaMP6f_MaxSlopeLM_reg <- format_reg_table(MaxSlopeLM_freq_GCaMP6f_regdata, "GCaMP6f", "Frequency", "MaxSlopeLMReg")

#### AUC Model -----
## Mixed effects model
AUCmodel_freq_GCaMP6f <- glmmTMB(AUC_POSTINJ ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                          data = subset(subjectrise_pksperminC_GCaMP6f, Virus=='GCaMP6f'&InjType=='M'), family = stats::gaussian)

AUCmodel_freq_GCaMP6f_aov <- Anova(AUCmodel_freq_GCaMP6f, type=2)
AUCmodel_freq_GCaMP6f_aov
summary(AUCmodel_freq_GCaMP6f)

AUCmodel_freq_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(AUCmodel_freq_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCmodel_freq_GCaMP6f_r2

## Prepare fit data for plotting
AUC_freq_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_GCaMP6f$DoseNum)))
AUC_freq_GCaMP6f_regdata$DoseNum_c <- AUC_freq_GCaMP6f_regdata$DoseNum - dosemean

AUC_freq_GCaMP6f_regpred <- predict(AUCmodel_freq_GCaMP6f, newdata = AUC_freq_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
AUC_freq_GCaMP6f_regdata$fit <- AUC_freq_GCaMP6f_regpred$fit
AUC_freq_GCaMP6f_regdata$se <- AUC_freq_GCaMP6f_regpred$se
AUC_freq_GCaMP6f_regdata$lower <- AUC_freq_GCaMP6f_regdata$fit - 1.96 * AUC_freq_GCaMP6f_regdata$se
AUC_freq_GCaMP6f_regdata$upper <- AUC_freq_GCaMP6f_regdata$fit + 1.96 * AUC_freq_GCaMP6f_regdata$se

AUC_freq_GCaMP6f_regdata

# Format Tables
table_freqmodel_GCaMP6f_AUC <- format_regsummary_table(AUCmodel_freq_GCaMP6f, 'GCaMP6f', "Frequency", "AUCModel")
table_freqmodel_GCaMP6f_AUC_aov <- format_anova_table(AUCmodel_freq_GCaMP6f_aov, "GCaMP6f", "Frequency", "AUCModel") 
table_freqmodel_GCaMP6f_AUC_reg <- format_reg_table(AUC_freq_GCaMP6f_regdata, "GCaMP6f", "Frequency", "AUCReg")

### Amplitude -----
#### Amplitude Model -----
# Mixed effects model
ampmodel_GCaMP6f <- glmmTMB(amp ~ InjType*Dose*BinF + (1|SubjectID), # random intercept by subject
                            data = subset(pkdataz, Virus=='GCaMP6f' & CondNum == 1), family = stats::gaussian, 
                            control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampmodel_GCaMP6f_aov <- Anova(ampmodel_GCaMP6f, type=2)
ampmodel_GCaMP6f_aov
summary(ampmodel_GCaMP6f)

ampmodel_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(ampmodel_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampmodel_GCaMP6f_r2

# Follow up tests - joint
ampmodel_GCaMP6f_emmeans <-emmeans(ampmodel_GCaMP6f, pairwise~InjType|Dose, adjust="none") # Find means
ampmodel_GCaMP6f_jt_Dose <- joint_tests(ampmodel_GCaMP6f_emmeans, by = "Dose")
ampmodel_GCaMP6f_jt_InjType <- joint_tests(ampmodel_GCaMP6f_emmeans, by = "InjType")

ampmodel_GCaMP6f_jt_Dose
ampmodel_GCaMP6f_jt_InjType

# Follow up tests - treat pairs
ampmodel_GCaMP6f_treatpairs <- emmeans(ampmodel_GCaMP6f, pairwise~InjType|Dose, adjust="none") # Pairs
ampmodel_GCaMP6f_treatcontrasts <- as.data.frame(ampmodel_GCaMP6f_treatpairs$contrasts) # Extract emmeans
ampmodel_GCaMP6f_treatemmeans <- ampmodel_GCaMP6f_treatpairs$emmeans # Extract emmeans
ampmodel_GCaMP6f_treatpairs_es <- as.data.frame(eff_size(ampmodel_GCaMP6f_treatemmeans, sigma = sigma(ampmodel_GCaMP6f), edf = df.residual(ampmodel_GCaMP6f))) # Effect Size
ampmodel_GCaMP6f_treatpairs_es <- ampmodel_GCaMP6f_treatpairs_es[, c("Dose", "contrast", "effect.size")]

ampmodel_GCaMP6f_treatcontrasts$d <- ampmodel_GCaMP6f_treatpairs_es$effect.size
ampmodel_GCaMP6f_treatcontrasts

# Follow up tests - dose pairs
ampmodel_GCaMP6f_dosepairs <- emmeans(ampmodel_GCaMP6f, pairwise~Dose|InjType, adjust="none") # Pairs
ampmodel_GCaMP6f_dosecontrasts <- as.data.frame(ampmodel_GCaMP6f_dosepairs$contrasts) # Extract emmeans
ampmodel_GCaMP6f_doseemmeans <- ampmodel_GCaMP6f_dosepairs$emmeans # Extract emmeans
ampmodel_GCaMP6f_dosepairs_es <- as.data.frame(eff_size(ampmodel_GCaMP6f_doseemmeans, sigma = sigma(ampmodel_GCaMP6f), edf = df.residual(ampmodel_GCaMP6f))) # Effect Size
ampmodel_GCaMP6f_dosepairs_es <- ampmodel_GCaMP6f_dosepairs_es[, c("InjType", "contrast", "effect.size")]

ampmodel_GCaMP6f_dosecontrasts$d <- ampmodel_GCaMP6f_dosepairs_es$effect.size
ampmodel_GCaMP6f_dosecontrasts

# Format Tables
table_ampmodel_GCaMP6f <- format_regsummary_table(ampmodel_GCaMP6f, 'GCaMP6f', "Amplitude", "StandardModel")
table_ampmodel_GCaMP6f_aov <- format_anova_table(ampmodel_GCaMP6f_aov, "GCaMP6f", "Amplitude", "StandardModel")
table_ampmodel_GCaMP6f_joint_Dose <- format_emm_table(ampmodel_GCaMP6f_jt_Dose, "GCaMP6f", "Amplitude", "Joint", "Dose")
table_ampmodel_GCaMP6f_joint_InjType <- format_emm_table(ampmodel_GCaMP6f_jt_InjType, "GCaMP6f", "Amplitude", "Joint", "InjType")
table_ampmodel_GCaMP6f_pairs_Dose <- format_emm_table(ampmodel_GCaMP6f_dosecontrasts, "GCaMP6f", "Amplitude", "Pairs", "Dose|InjType")
table_ampmodel_GCaMP6f_pairs_InjType <- format_emm_table(ampmodel_GCaMP6f_treatcontrasts, "GCaMP6f", "Amplitude", "Pairs", "InjType|Dose")

#### Mixed effects AR1 model
ampmodel_GCaMP6f_ar1 <- glmmTMB(amp ~ InjType*Dose*BinF + ar1(BinF + 0|SubjectID:InjType) + (1|SubjectID), # random intercept by subject
                                 data = subset(pkdataz, Virus=='GCaMP6f' & CondNum == 1), family = stats::gaussian, 
                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampmodel_GCaMP6f_ar1_aov <- Anova(ampmodel_GCaMP6f_ar1, type=2)
ampmodel_GCaMP6f_ar1_aov
summary(ampmodel_GCaMP6f_ar1)

ampmodel_GCaMP6f_ar1_r2 <- as.data.frame(r.squaredGLMM(ampmodel_GCaMP6f_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampmodel_GCaMP6f_ar1_r2

ampmodel_GCaMP6f_ar1_AIC <- as.data.frame(AIC(ampmodel_GCaMP6f, ampmodel_GCaMP6f_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "ampmodel_GCaMP6f" ~ "StandardModel", ModelName == "ampmodel_GCaMP6f_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

ampmodel_GCaMP6f_ar1_AIC

# Follow up tests - joint
ampmodel_GCaMP6f_ar1_emmeans <-emmeans(ampmodel_GCaMP6f_ar1, pairwise~InjType|Dose, adjust="none") # Find means
ampmodel_GCaMP6f_ar1_jt_Dose <- joint_tests(ampmodel_GCaMP6f_ar1, by = "Dose")
ampmodel_GCaMP6f_ar1_jt_InjType <- joint_tests(ampmodel_GCaMP6f_ar1, by = "InjType")

ampmodel_GCaMP6f_ar1_jt_Dose
ampmodel_GCaMP6f_ar1_jt_InjType

# Follow up tests - treat pairs
ampmodel_GCaMP6f_ar1_treatpairs <- emmeans(ampmodel_GCaMP6f_ar1, pairwise~InjType|Dose, adjust="none") # Pairs
ampmodel_GCaMP6f_ar1_treatcontrasts <- as.data.frame(ampmodel_GCaMP6f_ar1_treatpairs$contrasts) # Extract emmeans
ampmodel_GCaMP6f_ar1_treatemmeans <- ampmodel_GCaMP6f_ar1_treatpairs$emmeans # Extract emmeans
ampmodel_GCaMP6f_ar1_treatpairs_es <- as.data.frame(eff_size(ampmodel_GCaMP6f_ar1_treatemmeans, sigma = sigma(ampmodel_GCaMP6f_ar1), edf = df.residual(ampmodel_GCaMP6f_ar1))) # Effect Size
ampmodel_GCaMP6f_ar1_treatpairs_es <- ampmodel_GCaMP6f_ar1_treatpairs_es[, c("Dose", "contrast", "effect.size")]

ampmodel_GCaMP6f_ar1_treatcontrasts$d <- ampmodel_GCaMP6f_ar1_treatpairs_es$effect.size
ampmodel_GCaMP6f_ar1_treatcontrasts

# Follow up tests - dose pairs
ampmodel_GCaMP6f_ar1_dosepairs <- emmeans(ampmodel_GCaMP6f_ar1, pairwise~Dose|InjType, adjust="none") # Pairs
ampmodel_GCaMP6f_ar1_dosecontrasts <- as.data.frame(ampmodel_GCaMP6f_ar1_dosepairs$contrasts) # Extract emmeans
ampmodel_GCaMP6f_ar1_doseemmeans <- ampmodel_GCaMP6f_ar1_dosepairs$emmeans # Extract emmeans
ampmodel_GCaMP6f_ar1_dosepairs_es <- as.data.frame(eff_size(ampmodel_GCaMP6f_ar1_doseemmeans, sigma = sigma(ampmodel_GCaMP6f_ar1), edf = df.residual(ampmodel_GCaMP6f_ar1))) # Effect Size
ampmodel_GCaMP6f_ar1_dosepairs_es <- ampmodel_GCaMP6f_ar1_dosepairs_es[, c("InjType", "contrast", "effect.size")]

ampmodel_GCaMP6f_ar1_dosecontrasts$d <- ampmodel_GCaMP6f_ar1_dosepairs_es$effect.size
ampmodel_GCaMP6f_ar1_dosecontrasts

# Format Tables
table_ampmodel_GCaMP6f_ar1 <- format_regsummary_table(ampmodel_GCaMP6f_ar1, 'GCaMP6f', "Amplitude", "AR1Model")
table_ampmodel_GCaMP6f_ar1_aov <- format_anova_table(ampmodel_GCaMP6f_ar1_aov, "GCaMP6f", "Amplitude", "AR1Model") 
table_ampmodel_GCaMP6f_ar1_joint_Dose <- format_emm_table(ampmodel_GCaMP6f_ar1_jt_Dose, "GCaMP6f", "Amplitude", "Joint", "Dose")
table_ampmodel_GCaMP6f_ar1_joint_InjType <- format_emm_table(ampmodel_GCaMP6f_ar1_jt_InjType, "GCaMP6f", "Amplitude", "Joint", "InjType")
table_ampmodel_GCaMP6f_ar1_pairs_Dose <- format_emm_table(ampmodel_GCaMP6f_ar1_dosecontrasts, "GCaMP6f", "Amplitude", "Pairs", "Dose|InjType")
table_ampmodel_GCaMP6f_ar1_pairs_InjType <- format_emm_table(ampmodel_GCaMP6f_ar1_treatcontrasts, "GCaMP6f", "Amplitude", "Pairs", "InjType|Dose")

#### Amplitude Change Model -----
# Mixed effects model
ampCmodel_GCaMP6f <- glmmTMB(ampC ~ Dose*BinF + (1|SubjectID), # random intercept by subject
                                 data = subset(subjectmeansbin, Virus=='GCaMP6f' & CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampCmodel_GCaMP6f_aov <- Anova(ampCmodel_GCaMP6f, type=2)
ampCmodel_GCaMP6f_aov
summary(ampCmodel_GCaMP6f)

ampCmodel_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(ampCmodel_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampCmodel_GCaMP6f_r2

# Follow up tests - joint
ampCmodel_GCaMP6f_emmeans <-emmeans(ampCmodel_GCaMP6f, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
ampCmodel_GCaMP6f_dosepairs <- emmeans(ampCmodel_GCaMP6f, pairwise~Dose, adjust="none") # Pairs
ampCmodel_GCaMP6f_dosecontrasts <- as.data.frame(ampCmodel_GCaMP6f_dosepairs$contrasts) # Extract emmeans
ampCmodel_GCaMP6f_doseemmeans <- ampCmodel_GCaMP6f_dosepairs$emmeans # Extract emmeans
ampCmodel_GCaMP6f_dosepairs_es <- as.data.frame(eff_size(ampCmodel_GCaMP6f_doseemmeans, sigma = sigma(ampCmodel_GCaMP6f), edf = df.residual(ampCmodel_GCaMP6f))) # Effect Size
ampCmodel_GCaMP6f_dosepairs_es <- ampCmodel_GCaMP6f_dosepairs_es[, c("contrast", "effect.size")]

ampCmodel_GCaMP6f_dosecontrasts$d <- ampCmodel_GCaMP6f_dosepairs_es$effect.size
ampCmodel_GCaMP6f_dosecontrasts

# Format Tables
table_ampCmodel_GCaMP6f <- format_regsummary_table(ampCmodel_GCaMP6f, 'GCaMP6f', "Amplitude", "StandardModel")
table_ampCmodel_GCaMP6f_aov <- format_anova_table(ampCmodel_GCaMP6f_aov, "GCaMP6f", "Amplitude", "StandardModel")
table_ampCmodel_GCaMP6f_pairs_Dose <- format_emm_table(ampCmodel_GCaMP6f_dosecontrasts, "GCaMP6f", "Amplitude", "Pairs","Dose")

#### Mixed effects AR1 model
ampCmodel_GCaMP6f_ar1 <- glmmTMB(ampC ~ Dose*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                     data = subset(subjectmeansbin, Virus=='GCaMP6f' & CondNum == 1& InjType=='M'), family = stats::gaussian, 
                                     control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampCmodel_GCaMP6f_ar1_aov <- Anova(ampCmodel_GCaMP6f_ar1, type=2)
ampCmodel_GCaMP6f_ar1_aov
summary(ampCmodel_GCaMP6f_ar1)

ampCmodel_GCaMP6f_ar1_r2 <- as.data.frame(r.squaredGLMM(ampCmodel_GCaMP6f_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampCmodel_GCaMP6f_ar1_r2

ampCmodel_GCaMP6f_ar1_AIC <- as.data.frame(AIC(ampCmodel_GCaMP6f, ampCmodel_GCaMP6f_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "ampCmodel_GCaMP6f" ~ "StandardModel", ModelName == "ampCmodel_GCaMP6f_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

ampCmodel_GCaMP6f_ar1_AIC

# Follow up tests - joint
ampCmodel_GCaMP6f_ar1_emmeans <-emmeans(ampCmodel_GCaMP6f_ar1, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
ampCmodel_GCaMP6f_ar1_dosepairs <- emmeans(ampCmodel_GCaMP6f_ar1, pairwise~Dose, adjust="none") # Pairs
ampCmodel_GCaMP6f_ar1_dosecontrasts <- as.data.frame(ampCmodel_GCaMP6f_ar1_dosepairs$contrasts) # Extract emmeans
ampCmodel_GCaMP6f_ar1_doseemmeans <- ampCmodel_GCaMP6f_ar1_dosepairs$emmeans # Extract emmeans
ampCmodel_GCaMP6f_ar1_dosepairs_es <- as.data.frame(eff_size(ampCmodel_GCaMP6f_ar1_doseemmeans, sigma = sigma(ampCmodel_GCaMP6f_ar1), edf = df.residual(ampCmodel_GCaMP6f_ar1))) # Effect Size
ampCmodel_GCaMP6f_ar1_dosepairs_es <- ampCmodel_GCaMP6f_ar1_dosepairs_es[, c("contrast", "effect.size")]

ampCmodel_GCaMP6f_ar1_dosecontrasts$d <- ampCmodel_GCaMP6f_ar1_dosepairs_es$effect.size
ampCmodel_GCaMP6f_ar1_dosecontrasts

# Format Tables
table_ampCmodel_GCaMP6f_ar1 <- format_regsummary_table(ampCmodel_GCaMP6f_ar1, 'GCaMP6f', "Amplitude", "AR1Model")
table_ampCmodel_GCaMP6f_ar1_aov <- format_anova_table(ampCmodel_GCaMP6f_ar1_aov, "GCaMP6f", "Amplitude", "AR1Model") 
table_ampCmodel_GCaMP6f_ar1_pairs_Dose <- format_emm_table(ampCmodel_GCaMP6f_ar1_dosecontrasts, "GCaMP6f", "Amplitude", "Pairs", "Dose")


#### Max Height Model -----
## Mixed effects model
MaxHeightmodel_amp_GCaMP6f <- glmmTMB(MaxHeight ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                       data = subjectrise_ampC_GCaMP6f, family = stats::gaussian)

MaxHeightmodel_amp_GCaMP6f_aov <- Anova(MaxHeightmodel_amp_GCaMP6f, type=2)
MaxHeightmodel_amp_GCaMP6f_aov
summary(MaxHeightmodel_amp_GCaMP6f)

MaxHeightmodel_amp_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(MaxHeightmodel_amp_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxHeightmodel_amp_GCaMP6f_r2

## Prepare fit data for plotting
MaxHeight_amp_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_GCaMP6f$DoseNum)))
MaxHeight_amp_GCaMP6f_regdata$DoseNum_c <- MaxHeight_amp_GCaMP6f_regdata$DoseNum - dosemean

MaxHeight_amp_GCaMP6f_regpred <- predict(MaxHeightmodel_amp_GCaMP6f, newdata = MaxHeight_amp_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxHeight_amp_GCaMP6f_regdata$fit <- MaxHeight_amp_GCaMP6f_regpred$fit
MaxHeight_amp_GCaMP6f_regdata$se <- MaxHeight_amp_GCaMP6f_regpred$se
MaxHeight_amp_GCaMP6f_regdata$lower <- MaxHeight_amp_GCaMP6f_regdata$fit - 1.96 * MaxHeight_amp_GCaMP6f_regdata$se
MaxHeight_amp_GCaMP6f_regdata$upper <- MaxHeight_amp_GCaMP6f_regdata$fit + 1.96 * MaxHeight_amp_GCaMP6f_regdata$se

MaxHeight_amp_GCaMP6f_regdata

# Format Tables
table_ampmodel_GCaMP6f_MaxHeight <- format_regsummary_table(MaxHeightmodel_amp_GCaMP6f, 'GCaMP6f', "Amplitude", "MaxHeightModel")
table_ampmodel_GCaMP6f_MaxHeight_aov <- format_anova_table(MaxHeightmodel_amp_GCaMP6f_aov, "GCaMP6f", "Amplitude", "MaxHeightModel") 
table_ampmodel_GCaMP6f_MaxHeight_reg <- format_reg_table(MaxHeight_amp_GCaMP6f_regdata, "GCaMP6f", "Amplitude", "MaxHeightReg")

#### Max Duration Model -----
## Mixed effects model
MaxDurationmodel_amp_GCaMP6f <- glmmTMB(MaxDuration~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                         data = subjectrise_ampC_GCaMP6f, family = stats::gaussian)

MaxDurationmodel_amp_GCaMP6f_aov<- Anova(MaxDurationmodel_amp_GCaMP6f, type=2)
MaxDurationmodel_amp_GCaMP6f_aov
summary(MaxDurationmodel_amp_GCaMP6f)

MaxDurationmodel_amp_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(MaxDurationmodel_amp_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxDurationmodel_amp_GCaMP6f_r2

## Prepare fit data for plotting
MaxDuration_amp_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_GCaMP6f$DoseNum)))
MaxDuration_amp_GCaMP6f_regdata$DoseNum_c <- MaxDuration_amp_GCaMP6f_regdata$DoseNum - dosemean

MaxDuration_amp_GCaMP6f_regpred <- predict(MaxDurationmodel_amp_GCaMP6f, newdata = MaxDuration_amp_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxDuration_amp_GCaMP6f_regdata$fit <- MaxDuration_amp_GCaMP6f_regpred$fit
MaxDuration_amp_GCaMP6f_regdata$se <- MaxDuration_amp_GCaMP6f_regpred$se
MaxDuration_amp_GCaMP6f_regdata$lower <- MaxDuration_amp_GCaMP6f_regdata$fit - 1.96 * MaxDuration_amp_GCaMP6f_regdata$se
MaxDuration_amp_GCaMP6f_regdata$upper <- MaxDuration_amp_GCaMP6f_regdata$fit + 1.96 * MaxDuration_amp_GCaMP6f_regdata$se

MaxDuration_amp_GCaMP6f_regdata

# Format Tables
table_ampmodel_GCaMP6f_MaxDuration <- format_regsummary_table(MaxDurationmodel_amp_GCaMP6f, 'GCaMP6f', "Amplitude", "MaxDurationModel")
table_ampmodel_GCaMP6f_MaxDuration_aov <- format_anova_table(MaxDurationmodel_amp_GCaMP6f_aov, "GCaMP6f", "Amplitude", "MaxDurationModel") 
table_ampmodel_GCaMP6f_MaxDuration_reg <- format_reg_table(MaxDuration_amp_GCaMP6f_regdata, "GCaMP6f", "Amplitude", "MaxDurationReg")

#### Rise Slope Model -----
## Mixed effects model
MaxSlopeLMmodel_amp_GCaMP6f <- glmmTMB(MaxSlopeLM ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                        data = subjectrise_ampC_GCaMP6f, family = stats::gaussian)

MaxSlopeLMmodel_amp_GCaMP6f_aov <- Anova(MaxSlopeLMmodel_amp_GCaMP6f, type=2)
MaxSlopeLMmodel_amp_GCaMP6f_aov
summary(MaxSlopeLMmodel_amp_GCaMP6f)

MaxSlopeLMmodel_amp_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(MaxSlopeLMmodel_amp_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxSlopeLMmodel_amp_GCaMP6f_r2

## Prepare fit data for plotting
MaxSlopeLM_amp_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_GCaMP6f$DoseNum)))
MaxSlopeLM_amp_GCaMP6f_regdata$DoseNum_c <- MaxSlopeLM_amp_GCaMP6f_regdata$DoseNum - dosemean

MaxSlopeLM_amp_GCaMP6f_regpred <- predict(MaxSlopeLMmodel_amp_GCaMP6f, newdata = MaxSlopeLM_amp_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxSlopeLM_amp_GCaMP6f_regdata$fit <- MaxSlopeLM_amp_GCaMP6f_regpred$fit
MaxSlopeLM_amp_GCaMP6f_regdata$se <- MaxSlopeLM_amp_GCaMP6f_regpred$se
MaxSlopeLM_amp_GCaMP6f_regdata$lower <- MaxSlopeLM_amp_GCaMP6f_regdata$fit - 1.96 * MaxSlopeLM_amp_GCaMP6f_regdata$se
MaxSlopeLM_amp_GCaMP6f_regdata$upper <- MaxSlopeLM_amp_GCaMP6f_regdata$fit + 1.96 * MaxSlopeLM_amp_GCaMP6f_regdata$se

MaxSlopeLM_amp_GCaMP6f_regdata

# Format Tables
table_ampmodel_GCaMP6f_MaxSlopeLM <- format_regsummary_table(MaxSlopeLMmodel_amp_GCaMP6f, 'GCaMP6f', "Amplitude", "MaxSlopeLMModel")
table_ampmodel_GCaMP6f_MaxSlopeLM_aov <- format_anova_table(MaxSlopeLMmodel_amp_GCaMP6f_aov, "GCaMP6f", "Amplitude", "MaxSlopeLMModel") 
table_ampmodel_GCaMP6f_MaxSlopeLM_reg <- format_reg_table(MaxSlopeLM_amp_GCaMP6f_regdata, "GCaMP6f", "Amplitude", "MaxSlopeLMReg")

#### AUC Model -----
## Mixed effects model
AUCmodel_amp_GCaMP6f <- glmmTMB(AUC_POSTINJ ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                 data = subset(subjectrise_ampC_GCaMP6f, Virus=='GCaMP6f'&InjType=='M'), family = stats::gaussian)

AUCmodel_amp_GCaMP6f_aov <- Anova(AUCmodel_amp_GCaMP6f, type=2)
AUCmodel_amp_GCaMP6f_aov
summary(AUCmodel_amp_GCaMP6f)

AUCmodel_amp_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(AUCmodel_amp_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCmodel_amp_GCaMP6f_r2

## Prepare fit data for plotting
AUC_amp_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_GCaMP6f$DoseNum)))
AUC_amp_GCaMP6f_regdata$DoseNum_c <- AUC_amp_GCaMP6f_regdata$DoseNum - dosemean

AUC_amp_GCaMP6f_regpred <- predict(AUCmodel_amp_GCaMP6f, newdata = AUC_amp_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
AUC_amp_GCaMP6f_regdata$fit <- AUC_amp_GCaMP6f_regpred$fit
AUC_amp_GCaMP6f_regdata$se <- AUC_amp_GCaMP6f_regpred$se
AUC_amp_GCaMP6f_regdata$lower <- AUC_amp_GCaMP6f_regdata$fit - 1.96 * AUC_amp_GCaMP6f_regdata$se
AUC_amp_GCaMP6f_regdata$upper <- AUC_amp_GCaMP6f_regdata$fit + 1.96 * AUC_amp_GCaMP6f_regdata$se

AUC_amp_GCaMP6f_regdata

# Format Tables
table_ampmodel_GCaMP6f_AUC <- format_regsummary_table(AUCmodel_amp_GCaMP6f, 'GCaMP6f', "Amplitude", "AUCModel")
table_ampmodel_GCaMP6f_AUC_aov <- format_anova_table(AUCmodel_amp_GCaMP6f_aov, "GCaMP6f", "Amplitude", "AUCModel") 
table_ampmodel_GCaMP6f_AUC_reg <- format_reg_table(AUC_amp_GCaMP6f_regdata, "GCaMP6f", "Amplitude", "AUCReg")

### AUCwindow -----
#### AUCwindow Model -----
# Mixed effects model
AUCwindowmodel_GCaMP6f <- glmmTMB(AUCwindow ~ InjType*Dose*BinF + (1|SubjectID), # random intercept by subject
                            data = subset(pkdataz, Virus=='GCaMP6f' & CondNum == 1), family = stats::gaussian, 
                            control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowmodel_GCaMP6f_aov <- Anova(AUCwindowmodel_GCaMP6f, type=2)
AUCwindowmodel_GCaMP6f_aov
summary(AUCwindowmodel_GCaMP6f)

AUCwindowmodel_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(AUCwindowmodel_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowmodel_GCaMP6f_r2

# Follow up tests - joint
AUCwindowmodel_GCaMP6f_emmeans <-emmeans(AUCwindowmodel_GCaMP6f, pairwise~InjType|Dose, adjust="none") # Find means
AUCwindowmodel_GCaMP6f_jt_Dose <- joint_tests(AUCwindowmodel_GCaMP6f_emmeans, by = "Dose")
AUCwindowmodel_GCaMP6f_jt_InjType <- joint_tests(AUCwindowmodel_GCaMP6f_emmeans, by = "InjType")

AUCwindowmodel_GCaMP6f_jt_Dose
AUCwindowmodel_GCaMP6f_jt_InjType

# Follow up tests - treat pairs
AUCwindowmodel_GCaMP6f_treatpairs <- emmeans(AUCwindowmodel_GCaMP6f, pairwise~InjType|Dose, adjust="none") # Pairs
AUCwindowmodel_GCaMP6f_treatcontrasts <- as.data.frame(AUCwindowmodel_GCaMP6f_treatpairs$contrasts) # Extract emmeans
AUCwindowmodel_GCaMP6f_treatemmeans <- AUCwindowmodel_GCaMP6f_treatpairs$emmeans # Extract emmeans
AUCwindowmodel_GCaMP6f_treatpairs_es <- as.data.frame(eff_size(AUCwindowmodel_GCaMP6f_treatemmeans, sigma = sigma(AUCwindowmodel_GCaMP6f), edf = df.residual(AUCwindowmodel_GCaMP6f))) # Effect Size
AUCwindowmodel_GCaMP6f_treatpairs_es <- AUCwindowmodel_GCaMP6f_treatpairs_es[, c("Dose", "contrast", "effect.size")]

AUCwindowmodel_GCaMP6f_treatcontrasts$d <- AUCwindowmodel_GCaMP6f_treatpairs_es$effect.size
AUCwindowmodel_GCaMP6f_treatcontrasts

# Follow up tests - dose pairs
AUCwindowmodel_GCaMP6f_dosepairs <- emmeans(AUCwindowmodel_GCaMP6f, pairwise~Dose|InjType, adjust="none") # Pairs
AUCwindowmodel_GCaMP6f_dosecontrasts <- as.data.frame(AUCwindowmodel_GCaMP6f_dosepairs$contrasts) # Extract emmeans
AUCwindowmodel_GCaMP6f_doseemmeans <- AUCwindowmodel_GCaMP6f_dosepairs$emmeans # Extract emmeans
AUCwindowmodel_GCaMP6f_dosepairs_es <- as.data.frame(eff_size(AUCwindowmodel_GCaMP6f_doseemmeans, sigma = sigma(AUCwindowmodel_GCaMP6f), edf = df.residual(AUCwindowmodel_GCaMP6f))) # Effect Size
AUCwindowmodel_GCaMP6f_dosepairs_es <- AUCwindowmodel_GCaMP6f_dosepairs_es[, c("InjType", "contrast", "effect.size")]

AUCwindowmodel_GCaMP6f_dosecontrasts$d <- AUCwindowmodel_GCaMP6f_dosepairs_es$effect.size
AUCwindowmodel_GCaMP6f_dosecontrasts

# Format Tables
table_AUCwindowmodel_GCaMP6f <- format_regsummary_table(AUCwindowmodel_GCaMP6f, 'GCaMP6f', "AUCwindow", "StandardModel")
table_AUCwindowmodel_GCaMP6f_aov <- format_anova_table(AUCwindowmodel_GCaMP6f_aov, "GCaMP6f", "AUCwindow", "StandardModel")
table_AUCwindowmodel_GCaMP6f_joint_Dose <- format_emm_table(AUCwindowmodel_GCaMP6f_jt_Dose, "GCaMP6f", "AUCwindow", "Joint", "Dose")
table_AUCwindowmodel_GCaMP6f_joint_InjType <- format_emm_table(AUCwindowmodel_GCaMP6f_jt_InjType, "GCaMP6f", "AUCwindow", "Joint", "InjType")
table_AUCwindowmodel_GCaMP6f_pairs_Dose <- format_emm_table(AUCwindowmodel_GCaMP6f_dosecontrasts, "GCaMP6f", "AUCwindow", "Pairs", "Dose|InjType")
table_AUCwindowmodel_GCaMP6f_pairs_InjType <- format_emm_table(AUCwindowmodel_GCaMP6f_treatcontrasts, "GCaMP6f", "AUCwindow", "Pairs", "InjType|Dose")

#### Mixed effects AR1 model
AUCwindowmodel_GCaMP6f_ar1 <- glmmTMB(AUCwindow ~ InjType*Dose*BinF + ar1(BinF + 0|SubjectID:InjType) + (1|SubjectID), # random intercept by subject
                                data = subset(pkdataz, Virus=='GCaMP6f' & CondNum == 1), family = stats::gaussian, 
                                control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowmodel_GCaMP6f_ar1_aov <- Anova(AUCwindowmodel_GCaMP6f_ar1, type=2)
AUCwindowmodel_GCaMP6f_ar1_aov
summary(AUCwindowmodel_GCaMP6f_ar1)

AUCwindowmodel_GCaMP6f_ar1_r2 <- as.data.frame(r.squaredGLMM(AUCwindowmodel_GCaMP6f_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowmodel_GCaMP6f_ar1_r2

AUCwindowmodel_GCaMP6f_ar1_AIC <- as.data.frame(AIC(AUCwindowmodel_GCaMP6f, AUCwindowmodel_GCaMP6f_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "AUCwindowmodel_GCaMP6f" ~ "StandardModel", ModelName == "AUCwindowmodel_GCaMP6f_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

AUCwindowmodel_GCaMP6f_ar1_AIC

# Follow up tests - joint
AUCwindowmodel_GCaMP6f_ar1_emmeans <-emmeans(AUCwindowmodel_GCaMP6f_ar1, pairwise~InjType|Dose, adjust="none") # Find means
AUCwindowmodel_GCaMP6f_ar1_jt_Dose <- joint_tests(AUCwindowmodel_GCaMP6f_ar1, by = "Dose")
AUCwindowmodel_GCaMP6f_ar1_jt_InjType <- joint_tests(AUCwindowmodel_GCaMP6f_ar1, by = "InjType")

AUCwindowmodel_GCaMP6f_ar1_jt_Dose
AUCwindowmodel_GCaMP6f_ar1_jt_InjType

# Follow up tests - treat pairs
AUCwindowmodel_GCaMP6f_ar1_treatpairs <- emmeans(AUCwindowmodel_GCaMP6f_ar1, pairwise~InjType|Dose, adjust="none") # Pairs
AUCwindowmodel_GCaMP6f_ar1_treatcontrasts <- as.data.frame(AUCwindowmodel_GCaMP6f_ar1_treatpairs$contrasts) # Extract emmeans
AUCwindowmodel_GCaMP6f_ar1_treatemmeans <- AUCwindowmodel_GCaMP6f_ar1_treatpairs$emmeans # Extract emmeans
AUCwindowmodel_GCaMP6f_ar1_treatpairs_es <- as.data.frame(eff_size(AUCwindowmodel_GCaMP6f_ar1_treatemmeans, sigma = sigma(AUCwindowmodel_GCaMP6f_ar1), edf = df.residual(AUCwindowmodel_GCaMP6f_ar1))) # Effect Size
AUCwindowmodel_GCaMP6f_ar1_treatpairs_es <- AUCwindowmodel_GCaMP6f_ar1_treatpairs_es[, c("Dose", "contrast", "effect.size")]

AUCwindowmodel_GCaMP6f_ar1_treatcontrasts$d <- AUCwindowmodel_GCaMP6f_ar1_treatpairs_es$effect.size
AUCwindowmodel_GCaMP6f_ar1_treatcontrasts

# Follow up tests - dose pairs
AUCwindowmodel_GCaMP6f_ar1_dosepairs <- emmeans(AUCwindowmodel_GCaMP6f_ar1, pairwise~Dose|InjType, adjust="none") # Pairs
AUCwindowmodel_GCaMP6f_ar1_dosecontrasts <- as.data.frame(AUCwindowmodel_GCaMP6f_ar1_dosepairs$contrasts) # Extract emmeans
AUCwindowmodel_GCaMP6f_ar1_doseemmeans <- AUCwindowmodel_GCaMP6f_ar1_dosepairs$emmeans # Extract emmeans
AUCwindowmodel_GCaMP6f_ar1_dosepairs_es <- as.data.frame(eff_size(AUCwindowmodel_GCaMP6f_ar1_doseemmeans, sigma = sigma(AUCwindowmodel_GCaMP6f_ar1), edf = df.residual(AUCwindowmodel_GCaMP6f_ar1))) # Effect Size
AUCwindowmodel_GCaMP6f_ar1_dosepairs_es <- AUCwindowmodel_GCaMP6f_ar1_dosepairs_es[, c("InjType", "contrast", "effect.size")]

AUCwindowmodel_GCaMP6f_ar1_dosecontrasts$d <- AUCwindowmodel_GCaMP6f_ar1_dosepairs_es$effect.size
AUCwindowmodel_GCaMP6f_ar1_dosecontrasts

# Format Tables
table_AUCwindowmodel_GCaMP6f_ar1 <- format_regsummary_table(AUCwindowmodel_GCaMP6f_ar1, 'GCaMP6f', "AUCwindow", "AR1Model")
table_AUCwindowmodel_GCaMP6f_ar1_aov <- format_anova_table(AUCwindowmodel_GCaMP6f_ar1_aov, "GCaMP6f", "AUCwindow", "AR1Model") 
table_AUCwindowmodel_GCaMP6f_ar1_joint_Dose <- format_emm_table(AUCwindowmodel_GCaMP6f_ar1_jt_Dose, "GCaMP6f", "AUCwindow", "Joint", "Dose")
table_AUCwindowmodel_GCaMP6f_ar1_joint_InjType <- format_emm_table(AUCwindowmodel_GCaMP6f_ar1_jt_InjType, "GCaMP6f", "AUCwindow", "Joint", "InjType")
table_AUCwindowmodel_GCaMP6f_ar1_pairs_Dose <- format_emm_table(AUCwindowmodel_GCaMP6f_ar1_dosecontrasts, "GCaMP6f", "AUCwindow", "Pairs", "Dose|InjType")
table_AUCwindowmodel_GCaMP6f_ar1_pairs_InjType <- format_emm_table(AUCwindowmodel_GCaMP6f_ar1_treatcontrasts, "GCaMP6f", "AUCwindow", "Pairs", "InjType|Dose")

#### AUCwindow Change Model -----
# Mixed effects model
AUCwindowCmodel_GCaMP6f <- glmmTMB(AUCwindowC ~ Dose*BinF + (1|SubjectID), # random intercept by subject
                             data = subset(subjectmeansbin, Virus=='GCaMP6f' & CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                             control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowCmodel_GCaMP6f_aov <- Anova(AUCwindowCmodel_GCaMP6f, type=2)
AUCwindowCmodel_GCaMP6f_aov
summary(AUCwindowCmodel_GCaMP6f)

AUCwindowCmodel_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(AUCwindowCmodel_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowCmodel_GCaMP6f_r2

# Follow up tests - joint
AUCwindowCmodel_GCaMP6f_emmeans <-emmeans(AUCwindowCmodel_GCaMP6f, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
AUCwindowCmodel_GCaMP6f_dosepairs <- emmeans(AUCwindowCmodel_GCaMP6f, pairwise~Dose, adjust="none") # Pairs
AUCwindowCmodel_GCaMP6f_dosecontrasts <- as.data.frame(AUCwindowCmodel_GCaMP6f_dosepairs$contrasts) # Extract emmeans
AUCwindowCmodel_GCaMP6f_doseemmeans <- AUCwindowCmodel_GCaMP6f_dosepairs$emmeans # Extract emmeans
AUCwindowCmodel_GCaMP6f_dosepairs_es <- as.data.frame(eff_size(AUCwindowCmodel_GCaMP6f_doseemmeans, sigma = sigma(AUCwindowCmodel_GCaMP6f), edf = df.residual(AUCwindowCmodel_GCaMP6f))) # Effect Size
AUCwindowCmodel_GCaMP6f_dosepairs_es <- AUCwindowCmodel_GCaMP6f_dosepairs_es[, c("contrast", "effect.size")]

AUCwindowCmodel_GCaMP6f_dosecontrasts$d <- AUCwindowCmodel_GCaMP6f_dosepairs_es$effect.size
AUCwindowCmodel_GCaMP6f_dosecontrasts

# Format Tables
table_AUCwindowCmodel_GCaMP6f <- format_regsummary_table(AUCwindowCmodel_GCaMP6f, 'GCaMP6f', "Amplitude", "StandardModel")
table_AUCwindowCmodel_GCaMP6f_aov <- format_anova_table(AUCwindowCmodel_GCaMP6f_aov, "GCaMP6f", "Amplitude", "StandardModel")
table_AUCwindowCmodel_GCaMP6f_pairs_Dose <- format_emm_table(AUCwindowCmodel_GCaMP6f_dosecontrasts, "GCaMP6f", "Amplitude", "Pairs","Dose")

#### Mixed effects AR1 model
AUCwindowCmodel_GCaMP6f_ar1 <- glmmTMB(AUCwindowC ~ Dose*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                 data = subset(subjectmeansbin, Virus=='GCaMP6f' & CondNum == 1& InjType=='M'), family = stats::gaussian, 
                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowCmodel_GCaMP6f_ar1_aov <- Anova(AUCwindowCmodel_GCaMP6f_ar1, type=2)
AUCwindowCmodel_GCaMP6f_ar1_aov
summary(AUCwindowCmodel_GCaMP6f_ar1)

AUCwindowCmodel_GCaMP6f_ar1_r2 <- as.data.frame(r.squaredGLMM(AUCwindowCmodel_GCaMP6f_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowCmodel_GCaMP6f_ar1_r2

AUCwindowCmodel_GCaMP6f_ar1_AIC <- as.data.frame(AIC(AUCwindowCmodel_GCaMP6f, AUCwindowCmodel_GCaMP6f_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "AUCwindowCmodel_GCaMP6f" ~ "StandardModel", ModelName == "AUCwindowCmodel_GCaMP6f_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

AUCwindowCmodel_GCaMP6f_ar1_AIC

# Follow up tests - joint
AUCwindowCmodel_GCaMP6f_ar1_emmeans <-emmeans(AUCwindowCmodel_GCaMP6f_ar1, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
AUCwindowCmodel_GCaMP6f_ar1_dosepairs <- emmeans(AUCwindowCmodel_GCaMP6f_ar1, pairwise~Dose, adjust="none") # Pairs
AUCwindowCmodel_GCaMP6f_ar1_dosecontrasts <- as.data.frame(AUCwindowCmodel_GCaMP6f_ar1_dosepairs$contrasts) # Extract emmeans
AUCwindowCmodel_GCaMP6f_ar1_doseemmeans <- AUCwindowCmodel_GCaMP6f_ar1_dosepairs$emmeans # Extract emmeans
AUCwindowCmodel_GCaMP6f_ar1_dosepairs_es <- as.data.frame(eff_size(AUCwindowCmodel_GCaMP6f_ar1_doseemmeans, sigma = sigma(AUCwindowCmodel_GCaMP6f_ar1), edf = df.residual(AUCwindowCmodel_GCaMP6f_ar1))) # Effect Size
AUCwindowCmodel_GCaMP6f_ar1_dosepairs_es <- AUCwindowCmodel_GCaMP6f_ar1_dosepairs_es[, c("contrast", "effect.size")]

AUCwindowCmodel_GCaMP6f_ar1_dosecontrasts$d <- AUCwindowCmodel_GCaMP6f_ar1_dosepairs_es$effect.size
AUCwindowCmodel_GCaMP6f_ar1_dosecontrasts

# Format Tables
table_AUCwindowCmodel_GCaMP6f_ar1 <- format_regsummary_table(AUCwindowCmodel_GCaMP6f_ar1, 'GCaMP6f', "Amplitude", "AR1Model")
table_AUCwindowCmodel_GCaMP6f_ar1_aov <- format_anova_table(AUCwindowCmodel_GCaMP6f_ar1_aov, "GCaMP6f", "Amplitude", "AR1Model") 
table_AUCwindowCmodel_GCaMP6f_ar1_pairs_Dose <- format_emm_table(AUCwindowCmodel_GCaMP6f_ar1_dosecontrasts, "GCaMP6f", "Amplitude", "Pairs", "Dose")

#### Max Height Model -----
## Mixed effects model
MaxHeightmodel_AUCwindow_GCaMP6f <- glmmTMB(MaxHeight ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                      data = subjectrise_AUCwindowC_GCaMP6f, family = stats::gaussian)

MaxHeightmodel_AUCwindow_GCaMP6f_aov <- Anova(MaxHeightmodel_AUCwindow_GCaMP6f, type=2)
MaxHeightmodel_AUCwindow_GCaMP6f_aov
summary(MaxHeightmodel_AUCwindow_GCaMP6f)

MaxHeightmodel_AUCwindow_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(MaxHeightmodel_AUCwindow_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxHeightmodel_AUCwindow_GCaMP6f_r2

## Prepare fit data for plotting
MaxHeight_AUCwindow_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_GCaMP6f$DoseNum)))
MaxHeight_AUCwindow_GCaMP6f_regdata$DoseNum_c <- MaxHeight_AUCwindow_GCaMP6f_regdata$DoseNum - dosemean

MaxHeight_AUCwindow_GCaMP6f_regpred <- predict(MaxHeightmodel_AUCwindow_GCaMP6f, newdata = MaxHeight_AUCwindow_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxHeight_AUCwindow_GCaMP6f_regdata$fit <- MaxHeight_AUCwindow_GCaMP6f_regpred$fit
MaxHeight_AUCwindow_GCaMP6f_regdata$se <- MaxHeight_AUCwindow_GCaMP6f_regpred$se
MaxHeight_AUCwindow_GCaMP6f_regdata$lower <- MaxHeight_AUCwindow_GCaMP6f_regdata$fit - 1.96 * MaxHeight_AUCwindow_GCaMP6f_regdata$se
MaxHeight_AUCwindow_GCaMP6f_regdata$upper <- MaxHeight_AUCwindow_GCaMP6f_regdata$fit + 1.96 * MaxHeight_AUCwindow_GCaMP6f_regdata$se

MaxHeight_AUCwindow_GCaMP6f_regdata

# Format Tables
table_AUCwindowmodel_GCaMP6f_MaxHeight <- format_regsummary_table(MaxHeightmodel_AUCwindow_GCaMP6f, 'GCaMP6f', "AUCwindow", "MaxHeightModel")
table_AUCwindowmodel_GCaMP6f_MaxHeight_aov <- format_anova_table(MaxHeightmodel_AUCwindow_GCaMP6f_aov, "GCaMP6f", "AUCwindow", "MaxHeightModel") 
table_AUCwindowmodel_GCaMP6f_MaxHeight_reg <- format_reg_table(MaxHeight_AUCwindow_GCaMP6f_regdata, "GCaMP6f", "AUCwindow", "MaxHeightReg")

#### Max Duration Model -----
## Mixed effects model
MaxDurationmodel_AUCwindow_GCaMP6f <- glmmTMB(MaxDuration~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                        data = subjectrise_AUCwindowC_GCaMP6f, family = stats::gaussian)

MaxDurationmodel_AUCwindow_GCaMP6f_aov<- Anova(MaxDurationmodel_AUCwindow_GCaMP6f, type=2)
MaxDurationmodel_AUCwindow_GCaMP6f_aov
summary(MaxDurationmodel_AUCwindow_GCaMP6f)

MaxDurationmodel_AUCwindow_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(MaxDurationmodel_AUCwindow_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxDurationmodel_AUCwindow_GCaMP6f_r2

## Prepare fit data for plotting
MaxDuration_AUCwindow_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_GCaMP6f$DoseNum)))
MaxDuration_AUCwindow_GCaMP6f_regdata$DoseNum_c <- MaxDuration_AUCwindow_GCaMP6f_regdata$DoseNum - dosemean

MaxDuration_AUCwindow_GCaMP6f_regpred <- predict(MaxDurationmodel_AUCwindow_GCaMP6f, newdata = MaxDuration_AUCwindow_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxDuration_AUCwindow_GCaMP6f_regdata$fit <- MaxDuration_AUCwindow_GCaMP6f_regpred$fit
MaxDuration_AUCwindow_GCaMP6f_regdata$se <- MaxDuration_AUCwindow_GCaMP6f_regpred$se
MaxDuration_AUCwindow_GCaMP6f_regdata$lower <- MaxDuration_AUCwindow_GCaMP6f_regdata$fit - 1.96 * MaxDuration_AUCwindow_GCaMP6f_regdata$se
MaxDuration_AUCwindow_GCaMP6f_regdata$upper <- MaxDuration_AUCwindow_GCaMP6f_regdata$fit + 1.96 * MaxDuration_AUCwindow_GCaMP6f_regdata$se

MaxDuration_AUCwindow_GCaMP6f_regdata

# Format Tables
table_AUCwindowmodel_GCaMP6f_MaxDuration <- format_regsummary_table(MaxDurationmodel_AUCwindow_GCaMP6f, 'GCaMP6f', "AUCwindow", "MaxDurationModel")
table_AUCwindowmodel_GCaMP6f_MaxDuration_aov <- format_anova_table(MaxDurationmodel_AUCwindow_GCaMP6f_aov, "GCaMP6f", "AUCwindow", "MaxDurationModel") 
table_AUCwindowmodel_GCaMP6f_MaxDuration_reg <- format_reg_table(MaxDuration_AUCwindow_GCaMP6f_regdata, "GCaMP6f", "AUCwindow", "MaxDurationReg")

#### Rise Slope Model -----
## Mixed effects model
MaxSlopeLMmodel_AUCwindow_GCaMP6f <- glmmTMB(MaxSlopeLM ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                       data = subjectrise_AUCwindowC_GCaMP6f, family = stats::gaussian)

MaxSlopeLMmodel_AUCwindow_GCaMP6f_aov <- Anova(MaxSlopeLMmodel_AUCwindow_GCaMP6f, type=2)
MaxSlopeLMmodel_AUCwindow_GCaMP6f_aov
summary(MaxSlopeLMmodel_AUCwindow_GCaMP6f)

MaxSlopeLMmodel_AUCwindow_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(MaxSlopeLMmodel_AUCwindow_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxSlopeLMmodel_AUCwindow_GCaMP6f_r2

## Prepare fit data for plotting
MaxSlopeLM_AUCwindow_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_GCaMP6f$DoseNum)))
MaxSlopeLM_AUCwindow_GCaMP6f_regdata$DoseNum_c <- MaxSlopeLM_AUCwindow_GCaMP6f_regdata$DoseNum - dosemean

MaxSlopeLM_AUCwindow_GCaMP6f_regpred <- predict(MaxSlopeLMmodel_AUCwindow_GCaMP6f, newdata = MaxSlopeLM_AUCwindow_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxSlopeLM_AUCwindow_GCaMP6f_regdata$fit <- MaxSlopeLM_AUCwindow_GCaMP6f_regpred$fit
MaxSlopeLM_AUCwindow_GCaMP6f_regdata$se <- MaxSlopeLM_AUCwindow_GCaMP6f_regpred$se
MaxSlopeLM_AUCwindow_GCaMP6f_regdata$lower <- MaxSlopeLM_AUCwindow_GCaMP6f_regdata$fit - 1.96 * MaxSlopeLM_AUCwindow_GCaMP6f_regdata$se
MaxSlopeLM_AUCwindow_GCaMP6f_regdata$upper <- MaxSlopeLM_AUCwindow_GCaMP6f_regdata$fit + 1.96 * MaxSlopeLM_AUCwindow_GCaMP6f_regdata$se

MaxSlopeLM_AUCwindow_GCaMP6f_regdata

# Format Tables
table_AUCwindowmodel_GCaMP6f_MaxSlopeLM <- format_regsummary_table(MaxSlopeLMmodel_AUCwindow_GCaMP6f, 'GCaMP6f', "AUCwindow", "MaxSlopeLMModel")
table_AUCwindowmodel_GCaMP6f_MaxSlopeLM_aov <- format_anova_table(MaxSlopeLMmodel_AUCwindow_GCaMP6f_aov, "GCaMP6f", "AUCwindow", "MaxSlopeLMModel") 
table_AUCwindowmodel_GCaMP6f_MaxSlopeLM_reg <- format_reg_table(MaxSlopeLM_AUCwindow_GCaMP6f_regdata, "GCaMP6f", "AUCwindow", "MaxSlopeLMReg")

#### AUC Model -----
## Mixed effects model
AUCmodel_AUCwindow_GCaMP6f <- glmmTMB(AUC_POSTINJ ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                data = subset(subjectrise_AUCwindowC_GCaMP6f, Virus=='GCaMP6f'&InjType=='M'), family = stats::gaussian)

AUCmodel_AUCwindow_GCaMP6f_aov <- Anova(AUCmodel_AUCwindow_GCaMP6f, type=2)
AUCmodel_AUCwindow_GCaMP6f_aov
summary(AUCmodel_AUCwindow_GCaMP6f)

AUCmodel_AUCwindow_GCaMP6f_r2 <- as.data.frame(r.squaredGLMM(AUCmodel_AUCwindow_GCaMP6f)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCmodel_AUCwindow_GCaMP6f_r2

## Prepare fit data for plotting
AUC_AUCwindow_GCaMP6f_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_GCaMP6f$DoseNum)))
AUC_AUCwindow_GCaMP6f_regdata$DoseNum_c <- AUC_AUCwindow_GCaMP6f_regdata$DoseNum - dosemean

AUC_AUCwindow_GCaMP6f_regpred <- predict(AUCmodel_AUCwindow_GCaMP6f, newdata = AUC_AUCwindow_GCaMP6f_regdata, se.fit = TRUE, re.form = NA, type = "response")
AUC_AUCwindow_GCaMP6f_regdata$fit <- AUC_AUCwindow_GCaMP6f_regpred$fit
AUC_AUCwindow_GCaMP6f_regdata$se <- AUC_AUCwindow_GCaMP6f_regpred$se
AUC_AUCwindow_GCaMP6f_regdata$lower <- AUC_AUCwindow_GCaMP6f_regdata$fit - 1.96 * AUC_AUCwindow_GCaMP6f_regdata$se
AUC_AUCwindow_GCaMP6f_regdata$upper <- AUC_AUCwindow_GCaMP6f_regdata$fit + 1.96 * AUC_AUCwindow_GCaMP6f_regdata$se

AUC_AUCwindow_GCaMP6f_regdata

# Format Tables
table_AUCwindowmodel_GCaMP6f_AUC <- format_regsummary_table(AUCmodel_AUCwindow_GCaMP6f, 'GCaMP6f', "AUCwindow", "AUCModel")
table_AUCwindowmodel_GCaMP6f_AUC_aov <- format_anova_table(AUCmodel_AUCwindow_GCaMP6f_aov, "GCaMP6f", "AUCwindow", "AUCModel") 
table_AUCwindowmodel_GCaMP6f_AUC_reg <- format_reg_table(AUC_AUCwindow_GCaMP6f_regdata, "GCaMP6f", "AUCwindow", "AUCReg")


## EXPORT STATS -----
### Prepare Means -----
# Frequency
treatmeanstable_freq_GCaMP6f <- treatmeans %>% filter(Virus=='GCaMP6f') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         pkspermin_mean, pkspermin_sd, pkspermin_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansCtable_freq_GCaMP6f <- treatmeans %>% filter(InjType!='S', Virus=='GCaMP6f') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         pksperminC_mean, pksperminC_sd, pksperminC_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansSNRtable_freq_GCaMP6f <- treatmeans %>% filter(InjType!='S', Virus=='GCaMP6f') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         pksperminSNR_virus_mean, pksperminSNR_virus_sd, pksperminSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Amplitude
treatmeanstable_amp_GCaMP6f <- treatmeans %>% filter(Virus=='GCaMP6f') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         amp_mean, amp_sd, amp_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansCtable_amp_GCaMP6f <- treatmeans %>% filter(InjType!='S', Virus=='GCaMP6f') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         ampC_mean, ampC_sd, ampC_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansSNRtable_amp_GCaMP6f <- treatmeans %>% filter(InjType!='S', Virus=='GCaMP6f') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         ampSNR_virus_mean, ampSNR_virus_sd, ampSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# AUC Window
treatmeanstable_AUCwindow_GCaMP6f <- treatmeans %>% filter(Virus=='GCaMP6f') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         AUCwindow_mean, AUCwindow_sd, AUCwindow_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansCtable_AUCwindow_GCaMP6f <- treatmeans %>% filter(InjType!='S', Virus=='GCaMP6f') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         AUCwindowC_mean, AUCwindowC_sd, AUCwindowC_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansSNRtable_AUCwindow_GCaMP6f <- treatmeans %>% filter(InjType!='S', Virus=='GCaMP6f') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         AUCwindowSNR_virus_mean, AUCwindowSNR_virus_sd, AUCwindowSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))


### Frequency Export -----
wb_GCaMP6f_freq <- createWorkbook()
wb_GCaMP6f_freq_ar1 <- createWorkbook()
wb_GCaMP6f_freqC <- createWorkbook()
wb_GCaMP6f_freqC_ar1 <- createWorkbook()
wb_GCaMP6f_freq_reg <- createWorkbook()

## Main results for GCaMP6f
# Standard Models
write_nice_sheet(wb_GCaMP6f_freq, "MeansRaw", treatmeanstable_freq_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_freq, "ModelSummary", table_freqmodel_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_freq, "ModelAnova", table_freqmodel_GCaMP6f_aov)
write_nice_sheet(wb_GCaMP6f_freq, "ModelR2", freqmodel_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_freq, "Joint_InjType",  table_freqmodel_GCaMP6f_joint_InjType)
write_nice_sheet(wb_GCaMP6f_freq, "Joint_Dose", table_freqmodel_GCaMP6f_joint_Dose)
write_nice_sheet(wb_GCaMP6f_freq, "Pairs_InjType", table_freqmodel_GCaMP6f_pairs_InjType)
write_nice_sheet(wb_GCaMP6f_freq, "Pairs_Dose",  table_freqmodel_GCaMP6f_pairs_Dose)

saveWorkbook(wb_GCaMP6f_freq, paste(analysispath,"GCaMP6f_Frequency_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_GCaMP6f_freq_ar1, "MeansRaw", treatmeanstable_freq_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_freq_ar1, "ModelSummary", table_freqmodel_GCaMP6f_ar1)
write_nice_sheet(wb_GCaMP6f_freq_ar1, "ModelAnova", table_freqmodel_GCaMP6f_ar1_aov)
write_nice_sheet(wb_GCaMP6f_freq_ar1, "ModelR2", freqmodel_GCaMP6f_ar1_r2)
write_nice_sheet(wb_GCaMP6f_freq_ar1, "ModelAIC", freqmodel_GCaMP6f_ar1_AIC)
write_nice_sheet(wb_GCaMP6f_freq_ar1, "Joint_InjType",  table_freqmodel_GCaMP6f_ar1_joint_InjType)
write_nice_sheet(wb_GCaMP6f_freq_ar1, "Joint_Dose", table_freqmodel_GCaMP6f_ar1_joint_Dose)
write_nice_sheet(wb_GCaMP6f_freq_ar1, "Pairs_InjType", table_freqmodel_GCaMP6f_ar1_pairs_InjType)
write_nice_sheet(wb_GCaMP6f_freq_ar1, "Pairs_Dose",  table_freqmodel_GCaMP6f_ar1_pairs_Dose)

saveWorkbook(wb_GCaMP6f_freq_ar1, paste(analysispath,"GCaMP6f_Frequency_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Change Models
write_nice_sheet(wb_GCaMP6f_freqC, "MeansChange", treatmeanstable_freq_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_freqC, "ModelSummary", table_freqCmodel_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_freqC, "ModelAnova", table_freqCmodel_GCaMP6f_aov)
write_nice_sheet(wb_GCaMP6f_freqC, "ModelR2", freqCmodel_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_freqC, "Pairs_Dose",  table_freqCmodel_GCaMP6f_pairs_Dose)

saveWorkbook(wb_GCaMP6f_freqC, paste(analysispath,"GCaMP6f_FrequencyChange_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Change Models
write_nice_sheet(wb_GCaMP6f_freqC_ar1, "MeansChange", treatmeanstable_freq_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_freqC_ar1, "ModelSummary", table_freqCmodel_GCaMP6f_ar1)
write_nice_sheet(wb_GCaMP6f_freqC_ar1, "ModelAnova", table_freqCmodel_GCaMP6f_ar1_aov)
write_nice_sheet(wb_GCaMP6f_freqC_ar1, "ModelR2", freqCmodel_GCaMP6f_ar1_r2)
write_nice_sheet(wb_GCaMP6f_freqC_ar1, "ModelAIC", freqCmodel_GCaMP6f_ar1_AIC)
write_nice_sheet(wb_GCaMP6f_freqC_ar1, "Pairs_Dose",  table_freqCmodel_GCaMP6f_ar1_pairs_Dose)

saveWorkbook(wb_GCaMP6f_freqC_ar1, paste(analysispath,"GCaMP6f_FrequencyChange_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Regression Models
write_nice_sheet(wb_GCaMP6f_freq_reg, "MaxHeightSummary", table_freqmodel_GCaMP6f_MaxHeight)
write_nice_sheet(wb_GCaMP6f_freq_reg, "MaxHeightAnova", table_freqmodel_GCaMP6f_MaxHeight_aov)
write_nice_sheet(wb_GCaMP6f_freq_reg, "MaxHeightR2", MaxHeightmodel_freq_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_freq_reg, "MaxHeightReg",  table_freqmodel_GCaMP6f_MaxHeight_reg)
write_nice_sheet(wb_GCaMP6f_freq_reg, "MaxDurationSummary", table_freqmodel_GCaMP6f_MaxDuration)
write_nice_sheet(wb_GCaMP6f_freq_reg, "MaxDurationAnova", table_freqmodel_GCaMP6f_MaxDuration_aov)
write_nice_sheet(wb_GCaMP6f_freq_reg, "MaxDurationR2", MaxDurationmodel_freq_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_freq_reg, "MaxDurationReg", table_freqmodel_GCaMP6f_MaxDuration_reg)
write_nice_sheet(wb_GCaMP6f_freq_reg, "SlopeSummary", table_freqmodel_GCaMP6f_MaxSlopeLM)
write_nice_sheet(wb_GCaMP6f_freq_reg, "SlopeAnova", table_freqmodel_GCaMP6f_MaxSlopeLM_aov)
write_nice_sheet(wb_GCaMP6f_freq_reg, "SlopeR2", MaxSlopeLMmodel_freq_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_freq_reg, "SlopeReg", table_freqmodel_GCaMP6f_MaxSlopeLM_reg)
write_nice_sheet(wb_GCaMP6f_freq_reg, "AUCSummary", table_freqmodel_GCaMP6f_AUC)
write_nice_sheet(wb_GCaMP6f_freq_reg, "AUCAnova", table_freqmodel_GCaMP6f_AUC_aov)
write_nice_sheet(wb_GCaMP6f_freq_reg, "AUCR2", AUCmodel_freq_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_freq_reg, "AUCReg", table_freqmodel_GCaMP6f_AUC_reg)

saveWorkbook(wb_GCaMP6f_freq_reg, paste(analysispath,"GCaMP6f_Frequency_regressionmodels_results.xlsx",sep=''), overwrite = TRUE)

### Amplitude Export -----
wb_GCaMP6f_amp <- createWorkbook()
wb_GCaMP6f_amp_ar1 <- createWorkbook()
wb_GCaMP6f_ampC <- createWorkbook()
wb_GCaMP6f_ampC_ar1 <- createWorkbook()
wb_GCaMP6f_amp_reg <- createWorkbook()

## Main results for GCaMP6f
# Standard Models
write_nice_sheet(wb_GCaMP6f_amp, "MeansRaw", treatmeansCtable_amp_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_amp, "ModelSummary", table_ampmodel_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_amp, "ModelAnova", table_ampmodel_GCaMP6f_aov)
write_nice_sheet(wb_GCaMP6f_amp, "ModelR2", ampmodel_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_amp, "Joint_InjType",  table_ampmodel_GCaMP6f_joint_InjType)
write_nice_sheet(wb_GCaMP6f_amp, "Joint_Dose", table_ampmodel_GCaMP6f_joint_Dose)
write_nice_sheet(wb_GCaMP6f_amp, "Pairs_InjType", table_ampmodel_GCaMP6f_pairs_InjType)
write_nice_sheet(wb_GCaMP6f_amp, "Pairs_Dose",  table_ampmodel_GCaMP6f_pairs_Dose)

saveWorkbook(wb_GCaMP6f_amp, paste(analysispath,"GCaMP6f_Amplitude_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_GCaMP6f_amp_ar1, "MeansRaw", treatmeanstable_amp_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_amp_ar1, "ModelSummary", table_ampmodel_GCaMP6f_ar1)
write_nice_sheet(wb_GCaMP6f_amp_ar1, "ModelAnova", table_ampmodel_GCaMP6f_ar1_aov)
write_nice_sheet(wb_GCaMP6f_amp_ar1, "ModelR2", ampmodel_GCaMP6f_ar1_r2)
write_nice_sheet(wb_GCaMP6f_amp_ar1, "ModelAIC", ampmodel_GCaMP6f_ar1_AIC)
write_nice_sheet(wb_GCaMP6f_amp_ar1, "Joint_InjType",  table_ampmodel_GCaMP6f_ar1_joint_InjType)
write_nice_sheet(wb_GCaMP6f_amp_ar1, "Joint_Dose", table_ampmodel_GCaMP6f_ar1_joint_Dose)
write_nice_sheet(wb_GCaMP6f_amp_ar1, "Pairs_InjType", table_ampmodel_GCaMP6f_ar1_pairs_InjType)
write_nice_sheet(wb_GCaMP6f_amp_ar1, "Pairs_Dose",  table_ampmodel_GCaMP6f_ar1_pairs_Dose)

saveWorkbook(wb_GCaMP6f_amp_ar1, paste(analysispath,"GCaMP6f_Amplitude_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Change Models
write_nice_sheet(wb_GCaMP6f_ampC, "MeansChange", treatmeanstable_amp_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_ampC, "ModelSummary", table_ampCmodel_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_ampC, "ModelAnova", table_ampCmodel_GCaMP6f_aov)
write_nice_sheet(wb_GCaMP6f_ampC, "ModelR2", ampCmodel_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_ampC, "Pairs_Dose",  table_ampCmodel_GCaMP6f_pairs_Dose)

saveWorkbook(wb_GCaMP6f_ampC, paste(analysispath,"GCaMP6f_AmplitudeChange_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Change Models
write_nice_sheet(wb_GCaMP6f_ampC_ar1, "MeansChange", treatmeanstable_amp_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_ampC_ar1, "ModelSummary", table_ampCmodel_GCaMP6f_ar1)
write_nice_sheet(wb_GCaMP6f_ampC_ar1, "ModelAnova", table_ampCmodel_GCaMP6f_ar1_aov)
write_nice_sheet(wb_GCaMP6f_ampC_ar1, "ModelR2", ampCmodel_GCaMP6f_ar1_r2)
write_nice_sheet(wb_GCaMP6f_ampC_ar1, "ModelAIC", ampCmodel_GCaMP6f_ar1_AIC)
write_nice_sheet(wb_GCaMP6f_ampC_ar1, "Pairs_Dose",  table_ampCmodel_GCaMP6f_ar1_pairs_Dose)

saveWorkbook(wb_GCaMP6f_ampC_ar1, paste(analysispath,"GCaMP6f_AmplitudeChange_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Regression Models
write_nice_sheet(wb_GCaMP6f_amp_reg, "MaxHeightSummary", table_ampmodel_GCaMP6f_MaxHeight)
write_nice_sheet(wb_GCaMP6f_amp_reg, "MaxHeightAnova", table_ampmodel_GCaMP6f_MaxHeight_aov)
write_nice_sheet(wb_GCaMP6f_amp_reg, "MaxHeightR2", MaxHeightmodel_amp_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_amp_reg, "MaxHeightReg",  table_ampmodel_GCaMP6f_MaxHeight_reg)
write_nice_sheet(wb_GCaMP6f_amp_reg, "MaxDurationSummary", table_ampmodel_GCaMP6f_MaxDuration)
write_nice_sheet(wb_GCaMP6f_amp_reg, "MaxDurationAnova", table_ampmodel_GCaMP6f_MaxDuration_aov)
write_nice_sheet(wb_GCaMP6f_amp_reg, "MaxDurationR2", MaxDurationmodel_amp_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_amp_reg, "MaxDurationReg", table_ampmodel_GCaMP6f_MaxDuration_reg)
write_nice_sheet(wb_GCaMP6f_amp_reg, "SlopeSummary", table_ampmodel_GCaMP6f_MaxSlopeLM)
write_nice_sheet(wb_GCaMP6f_amp_reg, "SlopeAnova", table_ampmodel_GCaMP6f_MaxSlopeLM_aov)
write_nice_sheet(wb_GCaMP6f_amp_reg, "SlopeR2", MaxSlopeLMmodel_amp_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_amp_reg, "SlopeReg", table_ampmodel_GCaMP6f_MaxSlopeLM_reg)
write_nice_sheet(wb_GCaMP6f_amp_reg, "AUCSummary", table_ampmodel_GCaMP6f_AUC)
write_nice_sheet(wb_GCaMP6f_amp_reg, "AUCAnova", table_ampmodel_GCaMP6f_AUC_aov)
write_nice_sheet(wb_GCaMP6f_amp_reg, "AUCR2", AUCmodel_amp_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_amp_reg, "AUCReg", table_ampmodel_GCaMP6f_AUC_reg)

saveWorkbook(wb_GCaMP6f_amp_reg, paste(analysispath,"GCaMP6f_Amplitude_regressionmodels_results.xlsx",sep=''), overwrite = TRUE)


### AUCwindow Export -----
wb_GCaMP6f_AUCwindow <- createWorkbook()
wb_GCaMP6f_AUCwindow_ar1 <- createWorkbook()
wb_GCaMP6f_AUCwindowC <- createWorkbook()
wb_GCaMP6f_AUCwindowC_ar1 <- createWorkbook()
wb_GCaMP6f_AUCwindow_reg <- createWorkbook()

## Main results for GCaMP6f
# Standard Models
write_nice_sheet(wb_GCaMP6f_AUCwindow, "MeansRaw", treatmeanstable_AUCwindow_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_AUCwindow, "ModelSummary", table_AUCwindowmodel_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_AUCwindow, "ModelAnova", table_AUCwindowmodel_GCaMP6f_aov)
write_nice_sheet(wb_GCaMP6f_AUCwindow, "ModelR2", AUCwindowmodel_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_AUCwindow, "Joint_InjType",  table_AUCwindowmodel_GCaMP6f_joint_InjType)
write_nice_sheet(wb_GCaMP6f_AUCwindow, "Joint_Dose", table_AUCwindowmodel_GCaMP6f_joint_Dose)
write_nice_sheet(wb_GCaMP6f_AUCwindow, "Pairs_InjType", table_AUCwindowmodel_GCaMP6f_pairs_InjType)
write_nice_sheet(wb_GCaMP6f_AUCwindow, "Pairs_Dose",  table_AUCwindowmodel_GCaMP6f_pairs_Dose)

saveWorkbook(wb_GCaMP6f_AUCwindow, paste(analysispath,"GCaMP6f_AUCwindow_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_GCaMP6f_AUCwindow_ar1, "MeansRaw", treatmeanstable_AUCwindow_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_AUCwindow_ar1, "ModelSummary", table_AUCwindowmodel_GCaMP6f_ar1)
write_nice_sheet(wb_GCaMP6f_AUCwindow_ar1, "ModelAnova", table_AUCwindowmodel_GCaMP6f_ar1_aov)
write_nice_sheet(wb_GCaMP6f_AUCwindow_ar1, "ModelR2", AUCwindowmodel_GCaMP6f_ar1_r2)
write_nice_sheet(wb_GCaMP6f_AUCwindow_ar1, "ModelAIC", AUCwindowmodel_GCaMP6f_ar1_AIC)
write_nice_sheet(wb_GCaMP6f_AUCwindow_ar1, "Joint_InjType",  table_AUCwindowmodel_GCaMP6f_ar1_joint_InjType)
write_nice_sheet(wb_GCaMP6f_AUCwindow_ar1, "Joint_Dose", table_AUCwindowmodel_GCaMP6f_ar1_joint_Dose)
write_nice_sheet(wb_GCaMP6f_AUCwindow_ar1, "Pairs_InjType", table_AUCwindowmodel_GCaMP6f_ar1_pairs_InjType)
write_nice_sheet(wb_GCaMP6f_AUCwindow_ar1, "Pairs_Dose",  table_AUCwindowmodel_GCaMP6f_ar1_pairs_Dose)

saveWorkbook(wb_GCaMP6f_AUCwindow_ar1, paste(analysispath,"GCaMP6f_AUCwindow_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Change Models
write_nice_sheet(wb_GCaMP6f_AUCwindowC, "MeansChange", treatmeanstable_AUCwindow_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_AUCwindowC, "ModelSummary", table_AUCwindowCmodel_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_AUCwindowC, "ModelAnova", table_AUCwindowCmodel_GCaMP6f_aov)
write_nice_sheet(wb_GCaMP6f_AUCwindowC, "ModelR2", AUCwindowCmodel_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_AUCwindowC, "Pairs_Dose",  table_AUCwindowCmodel_GCaMP6f_pairs_Dose)

saveWorkbook(wb_GCaMP6f_AUCwindowC, paste(analysispath,"GCaMP6f_AUCwindowChange_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Change Models
write_nice_sheet(wb_GCaMP6f_AUCwindowC_ar1, "MeansChange", treatmeanstable_AUCwindow_GCaMP6f)
write_nice_sheet(wb_GCaMP6f_AUCwindowC_ar1, "ModelSummary", table_AUCwindowCmodel_GCaMP6f_ar1)
write_nice_sheet(wb_GCaMP6f_AUCwindowC_ar1, "ModelAnova", table_AUCwindowCmodel_GCaMP6f_ar1_aov)
write_nice_sheet(wb_GCaMP6f_AUCwindowC_ar1, "ModelR2", AUCwindowCmodel_GCaMP6f_ar1_r2)
write_nice_sheet(wb_GCaMP6f_AUCwindowC_ar1, "ModelAIC", AUCwindowCmodel_GCaMP6f_ar1_AIC)
write_nice_sheet(wb_GCaMP6f_AUCwindowC_ar1, "Pairs_Dose",  table_AUCwindowCmodel_GCaMP6f_ar1_pairs_Dose)

saveWorkbook(wb_GCaMP6f_AUCwindowC_ar1, paste(analysispath,"GCaMP6f_AUCwindowChange_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Regression Models
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "MaxHeightSummary", table_AUCwindowmodel_GCaMP6f_MaxHeight)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "MaxHeightAnova", table_AUCwindowmodel_GCaMP6f_MaxHeight_aov)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "MaxHeightR2", MaxHeightmodel_AUCwindow_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "MaxHeightReg",  table_AUCwindowmodel_GCaMP6f_MaxHeight_reg)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "MaxDurationSummary", table_AUCwindowmodel_GCaMP6f_MaxDuration)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "MaxDurationAnova", table_AUCwindowmodel_GCaMP6f_MaxDuration_aov)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "MaxDurationR2", MaxDurationmodel_AUCwindow_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "MaxDurationReg", table_AUCwindowmodel_GCaMP6f_MaxDuration_reg)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "SlopeSummary", table_AUCwindowmodel_GCaMP6f_MaxSlopeLM)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "SlopeAnova", table_AUCwindowmodel_GCaMP6f_MaxSlopeLM_aov)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "SlopeR2", MaxSlopeLMmodel_AUCwindow_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "SlopeReg", table_AUCwindowmodel_GCaMP6f_MaxSlopeLM_reg)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "AUCSummary", table_AUCwindowmodel_GCaMP6f_AUC)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "AUCAnova", table_AUCwindowmodel_GCaMP6f_AUC_aov)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "AUCR2", AUCmodel_AUCwindow_GCaMP6f_r2)
write_nice_sheet(wb_GCaMP6f_AUCwindow_reg, "AUCReg", table_AUCwindowmodel_GCaMP6f_AUC_reg)

saveWorkbook(wb_GCaMP6f_AUCwindow_reg, paste(analysispath,"GCaMP6f_AUCwindow_regressionmodels_results.xlsx",sep=''), overwrite = TRUE)


## PLOTS -----

### Frequency -----
# Bar Graph (Post-Inj Mean)
freqbarmin_GCaMP6f <- 0
freqbarmax_GCaMP6f <- 10.5
freqbarticks_GCaMP6f <- 2.5
freqbarbreaks_GCaMP6f <- seq(freqbarmin_GCaMP6f, freqbarmax_GCaMP6f, freqbarticks_GCaMP6f)

freqbar_GCaMP6f <- ggplot(subset(treatmeans, Virus == "GCaMP6f"&CondNum==1),aes(x=Dose, y=pkspermin_mean, fill=TreatNum, group=InjType)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_jitter(data=subset(subjectmeans,Virus == "GCaMP6f" & CondNum == 1), aes(x=Dose, y=pkspermin, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GCaMP6f"&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=pkspermin_mean-pkspermin_se, ymax=pkspermin_mean+pkspermin_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(freqbarmin_GCaMP6f, freqbarmax_GCaMP6f), breaks=freqbarbreaks_GCaMP6f) +
  xlab("Dose (mg/kg)") +
  ylab("Frequency (n/min)") + 
  mytheme +
  # STATS FOR AR1 MIXED MODEL
  geom_bracket(xmin=c(2.765), xmax = c(3.235), y.position=c(6.8),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(3.765), xmax = c(4.235), y.position=c(7.5),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

freqbar_GCaMP6f

# Bar graph (Change from Control)
freqCbarmin_GCaMP6f <- -3
freqCbarmax_GCaMP6f <- 9.5
freqCbarticks_GCaMP6f <- 3
freqCbarbreaks_GCaMP6f <- seq(-3,freqCbarmax_GCaMP6f,freqCbarticks_GCaMP6f)

freqCbar_GCaMP6f <- ggplot(subset(treatmeans, Virus == "GCaMP6f"&CondNum==1&InjType=='M'),aes(x=Dose, y=pksperminC_mean, fill=TreatNum, group=InjType)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(subjectmeans,Virus == "GCaMP6f" & CondNum == 1&InjType=='M'), aes(x=Dose, y=pksperminC, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GCaMP6f"&CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=pksperminC_mean-pksperminC_se, ymax=pksperminC_mean+pksperminC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(freqCbarmin_GCaMP6f, freqCbarmax_GCaMP6f), breaks=freqCbarbreaks_GCaMP6f)+
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" Frequency (n/min)")) + 
  mytheme +
  # AR1 and ALL AMP STATS
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(8.8),label.size=mystarsize, label = "",inherit.aes = FALSE, tip.length = c(.07, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(8.8),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2), xmax = c(3), y.position=c(7),label.size=mystarsize, label = "",inherit.aes = FALSE, tip.length = c(.07, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2), xmax = c(4), y.position=c(7),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3), xmax = c(4), y.position=c(5.4),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

freqCbar_GCaMP6f

# Line Graph (Change from Control)
freqlinemin_GCaMP6f <- -1.7
freqlinemax_GCaMP6f <- 3.3
freqlineticks_GCaMP6f <- 1
freqlinebreaks_GCaMP6f <- seq(-1, freqlinemax_GCaMP6f, freqlineticks_GCaMP6f)

freqline_GCaMP6f <- ggplot(subset(treatmeansbin, Virus=='GCaMP6f'&InjType!='S'), aes(x=BinTimeInj, y=pksperminC_mean, group=Dose, fill=Dose, color=Dose)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  scale_x_continuous(name="Minutes from Injection",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(freqlinemin_GCaMP6f, freqlinemax_GCaMP6f), breaks = freqlinebreaks_GCaMP6f)+
  scale_color_manual("legend", values = treatcolors) +
  mytheme +
  ylab(expression(Delta*" Frequency (n/min)")) + 
  mytheme

freqline_GCaMP6f

#### Temporal Dynamics Panels -----
#### Max Height
MaxHeight_freq_ymin_GCaMP6f <- -.8
MaxHeight_freq_ymax_GCaMP6f <- 15.2
MaxHeight_freq_yticks_GCaMP6f <- 5
MaxHeight_freq_ybreaks_GCaMP6f <- seq(0, MaxHeight_freq_ymax_GCaMP6f, MaxHeight_freq_yticks_GCaMP6f)

MaxHeight_freq_GCaMP6f <- ggplot(subset(subjectrise_pksperminC, Virus=='GCaMP6f'&InjType=='M'), aes(x=DoseNum, y=MaxHeight, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxHeight_freq_GCaMP6f_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxHeight_freq_GCaMP6f_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxHeight_freq_ymin_GCaMP6f, MaxHeight_freq_ymax_GCaMP6f), breaks = MaxHeight_freq_ybreaks_GCaMP6f) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Max "*Delta*" n")) + 
  annotate("text", x=10.6, y=6, label= "*", size = mystarsize*1.5, color = 'black')

MaxHeight_freq_GCaMP6f

#### Max Duration
MaxDuration_freq_ymin_GCaMP6f <- 0
MaxDuration_freq_ymax_GCaMP6f <- 63
MaxDuration_freq_yticks_GCaMP6f <- 20
MaxDuration_freq_ybreaks_GCaMP6f <- seq(MaxDuration_freq_ymin_GCaMP6f, MaxDuration_freq_ymax_GCaMP6f, MaxDuration_freq_yticks_GCaMP6f)

MaxDuration_freq_GCaMP6f <- ggplot(subset(subjectrise_pksperminC, Virus=='GCaMP6f'&InjType=='M'), aes(x=DoseNum, y=MaxDuration, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxDuration_freq_GCaMP6f_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxDuration_freq_GCaMP6f_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxDuration_freq_ymin_GCaMP6f, MaxDuration_freq_ymax_GCaMP6f), breaks = MaxDuration_freq_ybreaks_GCaMP6f) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Time to Max "*Delta*" n"))

MaxDuration_freq_GCaMP6f

#### Max Height Slope
MaxSlopeLM_freq_ymin_GCaMP6f <- -.1
MaxSlopeLM_freq_ymax_GCaMP6f <- 1
MaxSlopeLM_freq_yticks_GCaMP6f <- .25
MaxSlopeLM_freq_ybreaks_GCaMP6f <- seq(0, MaxSlopeLM_freq_ymax_GCaMP6f, MaxSlopeLM_freq_yticks_GCaMP6f)

MaxSlopeLM_freq_GCaMP6f <- ggplot(subset(subjectrise_pksperminC, Virus=='GCaMP6f'&InjType=='M'), aes(x=DoseNum, y=MaxSlopeLM, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxSlopeLM_freq_GCaMP6f_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxSlopeLM_freq_GCaMP6f_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxSlopeLM_freq_ymin_GCaMP6f, MaxSlopeLM_freq_ymax_GCaMP6f), breaks = MaxSlopeLM_freq_ybreaks_GCaMP6f)+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab('Time-Course Slope')

MaxSlopeLM_freq_GCaMP6f

#### AUC
AUC_freq_ymin_GCaMP6f <- -5
AUC_freq_ymax_GCaMP6f <- 405
AUC_freq_yticks_GCaMP6f <- 125
AUC_freq_ybreaks_GCaMP6f <- seq(0, AUC_freq_ymax_GCaMP6f, AUC_freq_yticks_GCaMP6f)

AUC_freq_GCaMP6f <- ggplot(subset(subjectrise_pksperminC, Virus=='GCaMP6f'&InjType=='M'), aes(x=DoseNum, y=AUC_POSTINJ, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = AUC_freq_GCaMP6f_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = AUC_freq_GCaMP6f_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(AUC_freq_ymin_GCaMP6f, AUC_freq_ymax_GCaMP6f), breaks = AUC_freq_ybreaks_GCaMP6f)+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab("Time-Course AUC") +
  annotate("text", x=10.6, y=135, label= "*", size = mystarsize*1.5, color = 'black')

AUC_freq_GCaMP6f


### Amplitude -----
# Bar Graph (Post-Injection Mean)
ampbarmin_GCaMP6f <- 0
ampbarmax_GCaMP6f <- 1.2
ampbarticks_GCaMP6f <- .3
ampbarybreaks_GCaMP6f <- seq(ampbarmin_GCaMP6f,ampbarmax_GCaMP6f,ampbarticks_GCaMP6f)
ampbarylabels_GCaMP6f <- ampbarybreaks_GCaMP6f + 2.6

ampbar_GCaMP6f <- ggplot(subset(treatmeans, Virus == "GCaMP6f"&CondNum==1),aes(x=Dose, y=ampscaled_mean, fill=TreatNum, group=InjType)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_jitter(data=subset(subjectmeans,Virus == "GCaMP6f" & CondNum == 1), aes(x=Dose, y=ampscaled, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GCaMP6f"&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=ampscaled_mean-amp_se, ymax=ampscaled_mean+amp_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(ampbarmin_GCaMP6f, ampbarmax_GCaMP6f), breaks=ampbarybreaks_GCaMP6f, labels=ampbarylabels_GCaMP6f)+
  xlab("Dose (mg/kg)") +
  ylab("Amplitude (z)") + 
  mytheme +
  # AR1 and ALL AMP STATS
  geom_bracket(xmin=c(0.765), xmax = c(1.235), y.position=c(1),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.765), xmax = c(4.235), y.position=c(1.12),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

ampbar_GCaMP6f

# Bar graph (Change from Control)
ampCbarmin_GCaMP6f <- -.55
ampCbarmax_GCaMP6f <- .8
ampCbarticks_GCaMP6f <- .4
ampCbarbreaks_GCaMP6f <- seq(-.4,ampCbarmax_GCaMP6f,ampCbarticks_GCaMP6f)

ampCbar_GCaMP6f <- ggplot(subset(treatmeans, Virus == "GCaMP6f"&CondNum==1&InjType=='M'),aes(x=Dose, y=ampC_mean, fill=TreatNum, group=InjType)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(subjectmeans,Virus == "GCaMP6f" & CondNum == 1&InjType=='M'), aes(x=Dose, y=ampC, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GCaMP6f"&CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=ampC_mean-ampC_se, ymax=ampC_mean+ampC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(ampCbarmin_GCaMP6f, ampCbarmax_GCaMP6f), breaks=ampCbarbreaks_GCaMP6f)+
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" Amplitude (z)")) + 
  mytheme +
  # AR1 and ALL AMP STATS
  geom_bracket(xmin=c(1), xmax = c(1.95), y.position=c(.6),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(.75),label.size=mystarsize, label = "",inherit.aes = FALSE, tip.length = c(.12, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(.75),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.05), xmax = c(4), y.position=c(.6),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

ampCbar_GCaMP6f

# Line Graph (Change from Control)
amplinemin_GCaMP6f <- -.6
amplinemax_GCaMP6f <- .6
amplineticks_GCaMP6f <- .3
amplinebreaks_GCaMP6f <- seq(amplinemin_GCaMP6f, amplinemax_GCaMP6f, amplineticks_GCaMP6f)

ampline_GCaMP6f <- ggplot(subset(treatmeansbin, Virus=='GCaMP6f'&InjType!='S'), aes(x=BinTimeInj, y=ampC_mean, group=Dose, fill=Dose, color=Dose)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  scale_x_continuous(name="Minutes from Injection",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(amplinemin_GCaMP6f, amplinemax_GCaMP6f), breaks = amplinebreaks_GCaMP6f)+
  scale_color_manual("legend", values = treatcolors) +
  scale_fill_manual("legend", values = treatcolors) +
  mytheme +
  ylab(expression(Delta*" Amplitude (z)")) + 
  mytheme

ampline_GCaMP6f

#### Temporal Dynamics Panels -----
# Max Height
MaxHeight_amp_ymin_GCaMP6f <- -.2
MaxHeight_amp_ymax_GCaMP6f <- 1.9
MaxHeight_amp_yticks_GCaMP6f <- .6
MaxHeight_amp_ybreaks_GCaMP6f <- seq(0, MaxHeight_amp_ymax_GCaMP6f, MaxHeight_amp_yticks_GCaMP6f)

MaxHeight_amp_GCaMP6f <- ggplot(subset(subjectrise_ampC, Virus=='GCaMP6f'&InjType=='M'), aes(x=DoseNum, y=MaxHeight, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxHeight_amp_GCaMP6f_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxHeight_amp_GCaMP6f_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxHeight_amp_ymin_GCaMP6f, MaxHeight_amp_ymax_GCaMP6f), breaks = MaxHeight_amp_ybreaks_GCaMP6f) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Max "*Delta*" z")) + 
  annotate("text", x=10.6, y=1.2, label= "**", size = mystarsize*1.5, color = 'black')

MaxHeight_amp_GCaMP6f

# Max Duration
MaxDuration_amp_ymin_GCaMP6f <- 0
MaxDuration_amp_ymax_GCaMP6f <- 63
MaxDuration_amp_yticks_GCaMP6f <- 20
MaxDuration_amp_ybreaks_GCaMP6f <- seq(MaxDuration_amp_ymin_GCaMP6f, MaxDuration_amp_ymax_GCaMP6f, MaxDuration_amp_yticks_GCaMP6f)

MaxDuration_amp_GCaMP6f <- ggplot(subset(subjectrise_ampC, Virus=='GCaMP6f'&InjType=='M'), aes(x=DoseNum, y=MaxDuration, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxDuration_amp_GCaMP6f_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxDuration_amp_GCaMP6f_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxDuration_amp_ymin_GCaMP6f, MaxDuration_amp_ymax_GCaMP6f), breaks = MaxDuration_amp_ybreaks_GCaMP6f) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Time to Max "*Delta*" z"))

MaxDuration_amp_GCaMP6f

# Rise Slope
MaxSlopeLM_amp_ymin_GCaMP6f <- -.02
MaxSlopeLM_amp_ymax_GCaMP6f <- .22
MaxSlopeLM_amp_yticks_GCaMP6f <- .07
MaxSlopeLM_amp_ybreaks_GCaMP6f <- seq(0, MaxSlopeLM_amp_ymax_GCaMP6f, MaxSlopeLM_amp_yticks_GCaMP6f)

MaxSlopeLM_amp_GCaMP6f <- ggplot(subset(subjectrise_ampC, Virus=='GCaMP6f'&InjType=='M'), aes(x=DoseNum, y=MaxSlopeLM, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxSlopeLM_amp_GCaMP6f_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxSlopeLM_amp_GCaMP6f_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxSlopeLM_amp_ymin_GCaMP6f, MaxSlopeLM_amp_ymax_GCaMP6f), breaks = MaxSlopeLM_amp_ybreaks_GCaMP6f)+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab('Time-Course Slope')+
  annotate("text", x=10.6, y=.069, label= "*", size = mystarsize*1.5, color = 'black')

MaxSlopeLM_amp_GCaMP6f

# AUC 
AUC_amp_ymin_GCaMP6f <- -1
AUC_amp_ymax_GCaMP6f <- 37
AUC_amp_yticks_GCaMP6f <- 12
AUC_amp_ybreaks_GCaMP6f <- seq(0, AUC_amp_ymax_GCaMP6f, AUC_amp_yticks_GCaMP6f)

AUC_amp_GCaMP6f <- ggplot(subset(subjectrise_ampC, Virus=='GCaMP6f'&InjType=='M'), aes(x=DoseNum, y=AUC_POSTINJ, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = AUC_amp_GCaMP6f_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = AUC_amp_GCaMP6f_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(AUC_amp_ymin_GCaMP6f, AUC_amp_ymax_GCaMP6f), breaks = AUC_amp_ybreaks_GCaMP6f)+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab("Time-Course AUC")

AUC_amp_GCaMP6f

### AUC Window -----
# Bar Graph (Post-Injection Mean)
AUCwindowbarmin_GCaMP6f <- 0
AUCwindowbarmax_GCaMP6f <- 3.2
AUCwindowbarticks_GCaMP6f <- .8
AUCwindowbarbreaks_GCaMP6f <- seq(AUCwindowbarmin_GCaMP6f,AUCwindowbarmax_GCaMP6f,AUCwindowbarticks_GCaMP6f)

AUCwindowbar_GCaMP6f <- ggplot(subset(treatmeans, Virus == "GCaMP6f"&CondNum==1),aes(x=Dose, y=AUCwindow_mean, fill=TreatNum, group=InjType)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_jitter(data=subset(subjectmeans,Virus == "GCaMP6f" & CondNum == 1), aes(x=Dose, y=AUCwindow, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GCaMP6f"&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=AUCwindow_mean-AUCwindow_se, ymax=AUCwindow_mean+AUCwindow_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(AUCwindowbarmin_GCaMP6f, AUCwindowbarmax_GCaMP6f), breaks=AUCwindowbarbreaks_GCaMP6f)+
  xlab("Dose (mg/kg)") +
  ylab("AUC per Transient") + 
  mytheme +
  # AR1 STATS
  geom_bracket(xmin=c(0.765), xmax = c(1.235), y.position=c(3),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1.765), xmax = c(2.235), y.position=c(3),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.765), xmax = c(4.235), y.position=c(3),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

AUCwindowbar_GCaMP6f

# Bar graph (Change from Control)
AUCwindowCbarmin_GCaMP6f <- -.6
AUCwindowCbarmax_GCaMP6f <- .62
AUCwindowCbarticks_GCaMP6f <- .3
AUCwindowCbarbreaks_GCaMP6f <- seq(AUCwindowCbarmin_GCaMP6f,AUCwindowCbarmax_GCaMP6f,AUCwindowCbarticks_GCaMP6f)

AUCwindowCbar_GCaMP6f <- ggplot(subset(treatmeans, Virus == "GCaMP6f"&CondNum==1&InjType=='M'),aes(x=Dose, y=AUCwindowC_mean, fill=TreatNum, group=InjType)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(subjectmeans,Virus == "GCaMP6f" & CondNum == 1&InjType=='M'), aes(x=Dose, y=AUCwindowC, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GCaMP6f"&CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=AUCwindowC_mean-AUCwindowC_se, ymax=AUCwindowC_mean+AUCwindowC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(AUCwindowCbarmin_GCaMP6f, AUCwindowCbarmax_GCaMP6f), breaks=AUCwindowCbarbreaks_GCaMP6f)+
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" AUC per Transient")) + 
  mytheme +
  # AR1 and ALL AUCwindow STATS
  geom_bracket(xmin=c(1), xmax = c(3), y.position=c(.58),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(4), y.position=c(.58),label.size=mystarsize, label = "",inherit.aes = FALSE, tip.length = c(.12, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2), xmax = c(3), y.position=c(.28),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2), xmax = c(4), y.position=c(.42),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.13, .03), size = .4, vjust = .5)

AUCwindowCbar_GCaMP6f

# Line Graph (Change from Control)
AUCwindowlineCmin_GCaMP6f <- -.8
AUCwindowlineCmax_GCaMP6f <- .8
AUCwindowlineCticks_GCaMP6f <- .4
AUCwindowlineCbreaks_GCaMP6f <- seq(AUCwindowlineCmin_GCaMP6f, AUCwindowlineCmax_GCaMP6f, AUCwindowlineCticks_GCaMP6f)

AUCwindowline_GCaMP6f <- ggplot(subset(treatmeansbin, Virus=='GCaMP6f'&InjType!='S'), aes(x=BinTimeInj, y=AUCwindowC_mean, group=Dose, fill=Dose, color=Dose)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  scale_x_continuous(name="Minutes from Injection",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(AUCwindowlineCmin_GCaMP6f, AUCwindowlineCmax_GCaMP6f), breaks = AUCwindowlineCbreaks_GCaMP6f)+
  scale_color_manual("legend", values = treatcolors) +
  scale_fill_manual("legend", values = treatcolors) +
  mytheme +
  ylab(expression(Delta*" AUC per Transient")) + 
  mytheme

AUCwindowline_GCaMP6f

### EXPORT FIGURE 2 -----
# BAR AND LINE GRAPHS
figure2panels_barline <- plot_grid(freqbar_GCaMP6f, freqCbar_GCaMP6f, freqline_GCaMP6f, 
                                   ampbar_GCaMP6f, ampCbar_GCaMP6f, ampline_GCaMP6f,
                                   AUCwindowbar_GCaMP6f, AUCwindowCbar_GCaMP6f, AUCwindowline_GCaMP6f,
                                   ncol = 3, nrow=3, align = "hv",  # align horizontally & vertically
                                   axis  = "tblr", hjust=-.5, rel_widths = c(1,.7,1.2,
                                                                            1,.7,1.2,
                                                                            1,.7,1.2))
figure2panels_barline

ggsave(filename = c("Figure2_VTAGCaMP6f_Bar_TransientPanels_BarLineGraphs.pdf"), path = figurepath, plot = figure2panels_barline, device = "pdf", 
       width = figpanelwidthin, height = figpanelheightinbarline, dpi=figdpi, units="in")


# TEMPORAL DYNAMICS GRAPHS
figure2panels_models <- plot_grid(MaxHeight_freq_GCaMP6f, MaxDuration_freq_GCaMP6f,
                                  MaxSlopeLM_freq_GCaMP6f, AUC_freq_GCaMP6f,
                                  
                                  MaxHeight_amp_GCaMP6f, MaxDuration_amp_GCaMP6f,
                                  MaxSlopeLM_amp_GCaMP6f, AUC_amp_GCaMP6f,
                                  
                                  ncol = 4, nrow=2, align = "hv",  # align horizontally & vertically
                                  axis  = "tblr", hjust=-.5)
figure2panels_models

ggsave(filename = c("Figure2_VTAGCaMP6f_Bar_TransientPanels_ModelGraphs.pdf"), path = figurepath, plot = figure2panels_models, device = "pdf", 
       width = figpanelwidthin, height = figpanelheightinmodel, dpi=figdpi, units="in")


# FIGURE 3 - NAcLS dLight1.3b -----
## STATS ------
### Frequency -----
#### Frequency Model -----
#### Mixed effects model
freqmodel_dLight1.3b <- glmmTMB(pkspermin ~ InjType*Dose*BinF + (1|SubjectID), # random intercept by subject
                             data = subset(subjectmeansbin, Virus=='dLight1.3b' & CondNum == 1), family = stats::gaussian,
                             control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqmodel_dLight1.3b_aov <- Anova(freqmodel_dLight1.3b, type=2)
freqmodel_dLight1.3b_aov
summary(freqmodel_dLight1.3b)

freqmodel_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(freqmodel_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqmodel_dLight1.3b_r2

# Follow up tests - joint
freqmodel_dLight1.3b_emmeans <-emmeans(freqmodel_dLight1.3b, pairwise~InjType|Dose, adjust="none") # Find means
freqmodel_dLight1.3b_jt_Dose <- joint_tests(freqmodel_dLight1.3b_emmeans, by = "Dose")
freqmodel_dLight1.3b_jt_InjType <- joint_tests(freqmodel_dLight1.3b_emmeans, by = "InjType")

freqmodel_dLight1.3b_jt_Dose
freqmodel_dLight1.3b_jt_InjType

# Follow up tests - treat pairs
freqmodel_dLight1.3b_treatpairs <- emmeans(freqmodel_dLight1.3b, pairwise~InjType|Dose, adjust="none") # Pairs
freqmodel_dLight1.3b_treatcontrasts <- as.data.frame(freqmodel_dLight1.3b_treatpairs$contrasts) # Extract emmeans
freqmodel_dLight1.3b_treatemmeans <- freqmodel_dLight1.3b_treatpairs$emmeans # Extract emmeans
freqmodel_dLight1.3b_treatpairs_es <- as.data.frame(eff_size(freqmodel_dLight1.3b_treatemmeans, sigma = sigma(freqmodel_dLight1.3b), edf = df.residual(freqmodel_dLight1.3b))) # Effect Size
freqmodel_dLight1.3b_treatpairs_es <- freqmodel_dLight1.3b_treatpairs_es[, c("Dose", "contrast", "effect.size")]

freqmodel_dLight1.3b_treatcontrasts$d <- freqmodel_dLight1.3b_treatpairs_es$effect.size
freqmodel_dLight1.3b_treatcontrasts

# Follow up tests - dose pairs
freqmodel_dLight1.3b_dosepairs <- emmeans(freqmodel_dLight1.3b, pairwise~Dose|InjType, adjust="none") # Pairs
freqmodel_dLight1.3b_dosecontrasts <- as.data.frame(freqmodel_dLight1.3b_dosepairs$contrasts) # Extract emmeans
freqmodel_dLight1.3b_doseemmeans <- freqmodel_dLight1.3b_dosepairs$emmeans # Extract emmeans
freqmodel_dLight1.3b_dosepairs_es <- as.data.frame(eff_size(freqmodel_dLight1.3b_doseemmeans, sigma = sigma(freqmodel_dLight1.3b), edf = df.residual(freqmodel_dLight1.3b))) # Effect Size
freqmodel_dLight1.3b_dosepairs_es <- freqmodel_dLight1.3b_dosepairs_es[, c("InjType", "contrast", "effect.size")]

freqmodel_dLight1.3b_dosecontrasts$d <- freqmodel_dLight1.3b_dosepairs_es$effect.size
freqmodel_dLight1.3b_dosecontrasts

# Format Tables
table_freqmodel_dLight1.3b <- format_regsummary_table(freqmodel_dLight1.3b, 'dLight1.3b', "Frequency", "StandardModel")
table_freqmodel_dLight1.3b_aov <- format_anova_table(freqmodel_dLight1.3b_aov, "dLight1.3b", "Frequency", "StandardModel")
table_freqmodel_dLight1.3b_joint_Dose <- format_emm_table(freqmodel_dLight1.3b_jt_Dose, "dLight1.3b", "Frequency", "Joint", "Dose")
table_freqmodel_dLight1.3b_joint_InjType <- format_emm_table(freqmodel_dLight1.3b_jt_InjType, "dLight1.3b", "Frequency", "Joint", "InjType")
table_freqmodel_dLight1.3b_pairs_Dose <- format_emm_table(freqmodel_dLight1.3b_dosecontrasts, "dLight1.3b", "Frequency", "Pairs", "Dose|InjType")
table_freqmodel_dLight1.3b_pairs_InjType <- format_emm_table(freqmodel_dLight1.3b_treatcontrasts, "dLight1.3b", "Frequency", "Pairs", "InjType|Dose")

#### Mixed effects AR1 model
freqmodel_dLight1.3b_ar1 <- glmmTMB(pkspermin ~ InjType*Dose*BinF + ar1(BinF + 0|SubjectID:InjType) + (1|SubjectID), # random intercept by subject
                                 data = subset(subjectmeansbin, Virus=='dLight1.3b' & CondNum == 1), family = stats::gaussian, 
                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqmodel_dLight1.3b_ar1_aov <- Anova(freqmodel_dLight1.3b_ar1, type=2)
freqmodel_dLight1.3b_ar1_aov
summary(freqmodel_dLight1.3b_ar1)

freqmodel_dLight1.3b_ar1_r2 <- as.data.frame(r.squaredGLMM(freqmodel_dLight1.3b_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqmodel_dLight1.3b_ar1_r2

freqmodel_dLight1.3b_ar1_AIC <- as.data.frame(AIC(freqmodel_dLight1.3b, freqmodel_dLight1.3b_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "freqmodel_dLight1.3b" ~ "StandardModel", ModelName == "freqmodel_dLight1.3b_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

freqmodel_dLight1.3b_ar1_AIC

# Follow up tests - joint
freqmodel_dLight1.3b_ar1_emmeans <-emmeans(freqmodel_dLight1.3b_ar1, pairwise~InjType|Dose, adjust="none") # Find means
freqmodel_dLight1.3b_ar1_jt_Dose <- joint_tests(freqmodel_dLight1.3b_ar1, by = "Dose")
freqmodel_dLight1.3b_ar1_jt_InjType <- joint_tests(freqmodel_dLight1.3b_ar1, by = "InjType")

freqmodel_dLight1.3b_ar1_jt_Dose
freqmodel_dLight1.3b_ar1_jt_InjType

# Follow up tests - treat pairs
freqmodel_dLight1.3b_ar1_treatpairs <- emmeans(freqmodel_dLight1.3b_ar1, pairwise~InjType|Dose, adjust="none") # Pairs
freqmodel_dLight1.3b_ar1_treatcontrasts <- as.data.frame(freqmodel_dLight1.3b_ar1_treatpairs$contrasts) # Extract emmeans
freqmodel_dLight1.3b_ar1_treatemmeans <- freqmodel_dLight1.3b_ar1_treatpairs$emmeans # Extract emmeans
freqmodel_dLight1.3b_ar1_treatpairs_es <- as.data.frame(eff_size(freqmodel_dLight1.3b_ar1_treatemmeans, sigma = sigma(freqmodel_dLight1.3b_ar1), edf = df.residual(freqmodel_dLight1.3b_ar1))) # Effect Size
freqmodel_dLight1.3b_ar1_treatpairs_es <- freqmodel_dLight1.3b_ar1_treatpairs_es[, c("Dose", "contrast", "effect.size")]

freqmodel_dLight1.3b_ar1_treatcontrasts$d <- freqmodel_dLight1.3b_ar1_treatpairs_es$effect.size
freqmodel_dLight1.3b_ar1_treatcontrasts

# Follow up tests - dose pairs
freqmodel_dLight1.3b_ar1_dosepairs <- emmeans(freqmodel_dLight1.3b_ar1, pairwise~Dose|InjType, adjust="none") # Pairs
freqmodel_dLight1.3b_ar1_dosecontrasts <- as.data.frame(freqmodel_dLight1.3b_ar1_dosepairs$contrasts) # Extract emmeans
freqmodel_dLight1.3b_ar1_doseemmeans <- freqmodel_dLight1.3b_ar1_dosepairs$emmeans # Extract emmeans
freqmodel_dLight1.3b_ar1_dosepairs_es <- as.data.frame(eff_size(freqmodel_dLight1.3b_ar1_doseemmeans, sigma = sigma(freqmodel_dLight1.3b_ar1), edf = df.residual(freqmodel_dLight1.3b_ar1))) # Effect Size
freqmodel_dLight1.3b_ar1_dosepairs_es <- freqmodel_dLight1.3b_ar1_dosepairs_es[, c("InjType", "contrast", "effect.size")]

freqmodel_dLight1.3b_ar1_dosecontrasts$d <- freqmodel_dLight1.3b_ar1_dosepairs_es$effect.size
freqmodel_dLight1.3b_ar1_dosecontrasts

# Format Tables
table_freqmodel_dLight1.3b_ar1 <- format_regsummary_table(freqmodel_dLight1.3b_ar1, 'dLight1.3b', "Frequency", "AR1Model")
table_freqmodel_dLight1.3b_ar1_aov <- format_anova_table(freqmodel_dLight1.3b_ar1_aov, "dLight1.3b", "Frequency", "AR1Model") 
table_freqmodel_dLight1.3b_ar1_joint_Dose <- format_emm_table(freqmodel_dLight1.3b_ar1_jt_Dose, "dLight1.3b", "Frequency", "Joint", "Dose")
table_freqmodel_dLight1.3b_ar1_joint_InjType <- format_emm_table(freqmodel_dLight1.3b_ar1_jt_InjType, "dLight1.3b", "Frequency", "Joint", "InjType")
table_freqmodel_dLight1.3b_ar1_pairs_Dose <- format_emm_table(freqmodel_dLight1.3b_ar1_dosecontrasts, "dLight1.3b", "Frequency", "Pairs", "Dose|InjType")
table_freqmodel_dLight1.3b_ar1_pairs_InjType <- format_emm_table(freqmodel_dLight1.3b_ar1_treatcontrasts, "dLight1.3b", "Frequency", "Pairs", "InjType|Dose")

#### Frequency Change Model -----
# Mixed effects model
freqCmodel_dLight1.3b <- glmmTMB(pksperminC ~ Dose*BinF + (1|SubjectID), # random intercept by subject
                              data = subset(subjectmeansbin, Virus=='dLight1.3b' & CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                              control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqCmodel_dLight1.3b_aov <- Anova(freqCmodel_dLight1.3b, type=2)
freqCmodel_dLight1.3b_aov
summary(freqCmodel_dLight1.3b)

freqCmodel_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(freqCmodel_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqCmodel_dLight1.3b_r2

# Follow up tests - joint
freqCmodel_dLight1.3b_emmeans <-emmeans(freqCmodel_dLight1.3b, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
freqCmodel_dLight1.3b_dosepairs <- emmeans(freqCmodel_dLight1.3b, pairwise~Dose, adjust="none") # Pairs
freqCmodel_dLight1.3b_dosecontrasts <- as.data.frame(freqCmodel_dLight1.3b_dosepairs$contrasts) # Extract emmeans
freqCmodel_dLight1.3b_doseemmeans <- freqCmodel_dLight1.3b_dosepairs$emmeans # Extract emmeans
freqCmodel_dLight1.3b_dosepairs_es <- as.data.frame(eff_size(freqCmodel_dLight1.3b_doseemmeans, sigma = sigma(freqCmodel_dLight1.3b), edf = df.residual(freqCmodel_dLight1.3b))) # Effect Size
freqCmodel_dLight1.3b_dosepairs_es <- freqCmodel_dLight1.3b_dosepairs_es[, c("contrast", "effect.size")]

freqCmodel_dLight1.3b_dosecontrasts$d <- freqCmodel_dLight1.3b_dosepairs_es$effect.size
freqCmodel_dLight1.3b_dosecontrasts

# Format Tables
table_freqCmodel_dLight1.3b <- format_regsummary_table(freqCmodel_dLight1.3b, 'dLight1.3b', "freqlitude", "StandardModel")
table_freqCmodel_dLight1.3b_aov <- format_anova_table(freqCmodel_dLight1.3b_aov, "dLight1.3b", "freqlitude", "StandardModel")
table_freqCmodel_dLight1.3b_pairs_Dose <- format_emm_table(freqCmodel_dLight1.3b_dosecontrasts, "dLight1.3b", "freqlitude", "Pairs","Dose")

#### Mixed effects AR1 model
freqCmodel_dLight1.3b_ar1 <- glmmTMB(pksperminC ~ Dose*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                  data = subset(subjectmeansbin, Virus=='dLight1.3b' & CondNum == 1& InjType=='M'), family = stats::gaussian, 
                                  control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqCmodel_dLight1.3b_ar1_aov <- Anova(freqCmodel_dLight1.3b_ar1, type=2)
freqCmodel_dLight1.3b_ar1_aov
summary(freqCmodel_dLight1.3b_ar1)

freqCmodel_dLight1.3b_ar1_r2 <- as.data.frame(r.squaredGLMM(freqCmodel_dLight1.3b_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqCmodel_dLight1.3b_ar1_r2

freqCmodel_dLight1.3b_ar1_AIC <- as.data.frame(AIC(freqCmodel_dLight1.3b, freqCmodel_dLight1.3b_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "freqCmodel_dLight1.3b" ~ "StandardModel", ModelName == "freqCmodel_dLight1.3b_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

freqCmodel_dLight1.3b_ar1_AIC

# Follow up tests - joint
freqCmodel_dLight1.3b_ar1_emmeans <-emmeans(freqCmodel_dLight1.3b_ar1, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
freqCmodel_dLight1.3b_ar1_dosepairs <- emmeans(freqCmodel_dLight1.3b_ar1, pairwise~Dose, adjust="none") # Pairs
freqCmodel_dLight1.3b_ar1_dosecontrasts <- as.data.frame(freqCmodel_dLight1.3b_ar1_dosepairs$contrasts) # Extract emmeans
freqCmodel_dLight1.3b_ar1_doseemmeans <- freqCmodel_dLight1.3b_ar1_dosepairs$emmeans # Extract emmeans
freqCmodel_dLight1.3b_ar1_dosepairs_es <- as.data.frame(eff_size(freqCmodel_dLight1.3b_ar1_doseemmeans, sigma = sigma(freqCmodel_dLight1.3b_ar1), edf = df.residual(freqCmodel_dLight1.3b_ar1))) # Effect Size
freqCmodel_dLight1.3b_ar1_dosepairs_es <- freqCmodel_dLight1.3b_ar1_dosepairs_es[, c("contrast", "effect.size")]

freqCmodel_dLight1.3b_ar1_dosecontrasts$d <- freqCmodel_dLight1.3b_ar1_dosepairs_es$effect.size
freqCmodel_dLight1.3b_ar1_dosecontrasts

# Format Tables
table_freqCmodel_dLight1.3b_ar1 <- format_regsummary_table(freqCmodel_dLight1.3b_ar1, 'dLight1.3b', "freqlitude", "AR1Model")
table_freqCmodel_dLight1.3b_ar1_aov <- format_anova_table(freqCmodel_dLight1.3b_ar1_aov, "dLight1.3b", "freqlitude", "AR1Model") 
table_freqCmodel_dLight1.3b_ar1_pairs_Dose <- format_emm_table(freqCmodel_dLight1.3b_ar1_dosecontrasts, "dLight1.3b", "freqlitude", "Pairs", "Dose")

#### Max Height Model -----
## Mixed effects model
MaxHeightmodel_freq_dLight1.3b <- glmmTMB(MaxHeight ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                       data = subjectrise_pksperminC_dLight1.3b, family = stats::gaussian)

MaxHeightmodel_freq_dLight1.3b_aov <- Anova(MaxHeightmodel_freq_dLight1.3b, type=2)
MaxHeightmodel_freq_dLight1.3b_aov
summary(MaxHeightmodel_freq_dLight1.3b)

MaxHeightmodel_freq_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(MaxHeightmodel_freq_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxHeightmodel_freq_dLight1.3b_r2

## Prepare fit data for plotting
MaxHeight_freq_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_dLight1.3b$DoseNum)))
MaxHeight_freq_dLight1.3b_regdata$DoseNum_c <- MaxHeight_freq_dLight1.3b_regdata$DoseNum - dosemean

MaxHeight_freq_dLight1.3b_regpred <- predict(MaxHeightmodel_freq_dLight1.3b, newdata = MaxHeight_freq_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxHeight_freq_dLight1.3b_regdata$fit <- MaxHeight_freq_dLight1.3b_regpred$fit
MaxHeight_freq_dLight1.3b_regdata$se <- MaxHeight_freq_dLight1.3b_regpred$se
MaxHeight_freq_dLight1.3b_regdata$lower <- MaxHeight_freq_dLight1.3b_regdata$fit - 1.96 * MaxHeight_freq_dLight1.3b_regdata$se
MaxHeight_freq_dLight1.3b_regdata$upper <- MaxHeight_freq_dLight1.3b_regdata$fit + 1.96 * MaxHeight_freq_dLight1.3b_regdata$se

MaxHeight_freq_dLight1.3b_regdata

# Format Tables
table_freqmodel_dLight1.3b_MaxHeight <- format_regsummary_table(MaxHeightmodel_freq_dLight1.3b, 'dLight1.3b', "Frequency", "MaxHeightModel")
table_freqmodel_dLight1.3b_MaxHeight_aov <- format_anova_table(MaxHeightmodel_freq_dLight1.3b_aov, "dLight1.3b", "Frequency", "MaxHeightModel") 
table_freqmodel_dLight1.3b_MaxHeight_reg <- format_reg_table(MaxHeight_freq_dLight1.3b_regdata, "dLight1.3b", "Frequency", "MaxHeightReg")

#### Max Duration Model -----
## Mixed effects model
MaxDurationmodel_freq_dLight1.3b <- glmmTMB(MaxDuration~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                         data = subjectrise_pksperminC_dLight1.3b, family = stats::gaussian)

MaxDurationmodel_freq_dLight1.3b_aov<- Anova(MaxDurationmodel_freq_dLight1.3b, type=2)
MaxDurationmodel_freq_dLight1.3b_aov
summary(MaxDurationmodel_freq_dLight1.3b)

MaxDurationmodel_freq_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(MaxDurationmodel_freq_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxDurationmodel_freq_dLight1.3b_r2

## Prepare fit data for plotting
MaxDuration_freq_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_dLight1.3b$DoseNum)))
MaxDuration_freq_dLight1.3b_regdata$DoseNum_c <- MaxDuration_freq_dLight1.3b_regdata$DoseNum - dosemean

MaxDuration_freq_dLight1.3b_regpred <- predict(MaxDurationmodel_freq_dLight1.3b, newdata = MaxDuration_freq_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxDuration_freq_dLight1.3b_regdata$fit <- MaxDuration_freq_dLight1.3b_regpred$fit
MaxDuration_freq_dLight1.3b_regdata$se <- MaxDuration_freq_dLight1.3b_regpred$se
MaxDuration_freq_dLight1.3b_regdata$lower <- MaxDuration_freq_dLight1.3b_regdata$fit - 1.96 * MaxDuration_freq_dLight1.3b_regdata$se
MaxDuration_freq_dLight1.3b_regdata$upper <- MaxDuration_freq_dLight1.3b_regdata$fit + 1.96 * MaxDuration_freq_dLight1.3b_regdata$se

MaxDuration_freq_dLight1.3b_regdata

# Format Tables
table_freqmodel_dLight1.3b_MaxDuration <- format_regsummary_table(MaxDurationmodel_freq_dLight1.3b, 'dLight1.3b', "Frequency", "MaxDurationModel")
table_freqmodel_dLight1.3b_MaxDuration_aov <- format_anova_table(MaxDurationmodel_freq_dLight1.3b_aov, "dLight1.3b", "Frequency", "MaxDurationModel") 
table_freqmodel_dLight1.3b_MaxDuration_reg <- format_reg_table(MaxDuration_freq_dLight1.3b_regdata, "dLight1.3b", "Frequency", "MaxDurationReg")

#### Rise Slope Model -----
## Mixed effects model
MaxSlopeLMmodel_freq_dLight1.3b <- glmmTMB(MaxSlopeLM ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                        data = subjectrise_pksperminC_dLight1.3b, family = stats::gaussian)

MaxSlopeLMmodel_freq_dLight1.3b_aov <- Anova(MaxSlopeLMmodel_freq_dLight1.3b, type=2)
MaxSlopeLMmodel_freq_dLight1.3b_aov
summary(MaxSlopeLMmodel_freq_dLight1.3b)

MaxSlopeLMmodel_freq_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(MaxSlopeLMmodel_freq_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxSlopeLMmodel_freq_dLight1.3b_r2

## Prepare fit data for plotting
MaxSlopeLM_freq_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_dLight1.3b$DoseNum)))
MaxSlopeLM_freq_dLight1.3b_regdata$DoseNum_c <- MaxSlopeLM_freq_dLight1.3b_regdata$DoseNum - dosemean

MaxSlopeLM_freq_dLight1.3b_regpred <- predict(MaxSlopeLMmodel_freq_dLight1.3b, newdata = MaxSlopeLM_freq_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxSlopeLM_freq_dLight1.3b_regdata$fit <- MaxSlopeLM_freq_dLight1.3b_regpred$fit
MaxSlopeLM_freq_dLight1.3b_regdata$se <- MaxSlopeLM_freq_dLight1.3b_regpred$se
MaxSlopeLM_freq_dLight1.3b_regdata$lower <- MaxSlopeLM_freq_dLight1.3b_regdata$fit - 1.96 * MaxSlopeLM_freq_dLight1.3b_regdata$se
MaxSlopeLM_freq_dLight1.3b_regdata$upper <- MaxSlopeLM_freq_dLight1.3b_regdata$fit + 1.96 * MaxSlopeLM_freq_dLight1.3b_regdata$se

MaxSlopeLM_freq_dLight1.3b_regdata

# Format Tables
table_freqmodel_dLight1.3b_MaxSlopeLM <- format_regsummary_table(MaxSlopeLMmodel_freq_dLight1.3b, 'dLight1.3b', "Frequency", "MaxSlopeLMModel")
table_freqmodel_dLight1.3b_MaxSlopeLM_aov <- format_anova_table(MaxSlopeLMmodel_freq_dLight1.3b_aov, "dLight1.3b", "Frequency", "MaxSlopeLMModel") 
table_freqmodel_dLight1.3b_MaxSlopeLM_reg <- format_reg_table(MaxSlopeLM_freq_dLight1.3b_regdata, "dLight1.3b", "Frequency", "MaxSlopeLMReg")

#### AUC Model -----
## Mixed effects model
AUCmodel_freq_dLight1.3b <- glmmTMB(AUC_POSTINJ ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                 data = subset(subjectrise_pksperminC_dLight1.3b, Virus=='dLight1.3b'&InjType=='M'), family = stats::gaussian)

AUCmodel_freq_dLight1.3b_aov <- Anova(AUCmodel_freq_dLight1.3b, type=2)
AUCmodel_freq_dLight1.3b_aov
summary(AUCmodel_freq_dLight1.3b)

AUCmodel_freq_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(AUCmodel_freq_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCmodel_freq_dLight1.3b_r2

## Prepare fit data for plotting
AUC_freq_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_dLight1.3b$DoseNum)))
AUC_freq_dLight1.3b_regdata$DoseNum_c <- AUC_freq_dLight1.3b_regdata$DoseNum - dosemean

AUC_freq_dLight1.3b_regpred <- predict(AUCmodel_freq_dLight1.3b, newdata = AUC_freq_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
AUC_freq_dLight1.3b_regdata$fit <- AUC_freq_dLight1.3b_regpred$fit
AUC_freq_dLight1.3b_regdata$se <- AUC_freq_dLight1.3b_regpred$se
AUC_freq_dLight1.3b_regdata$lower <- AUC_freq_dLight1.3b_regdata$fit - 1.96 * AUC_freq_dLight1.3b_regdata$se
AUC_freq_dLight1.3b_regdata$upper <- AUC_freq_dLight1.3b_regdata$fit + 1.96 * AUC_freq_dLight1.3b_regdata$se

AUC_freq_dLight1.3b_regdata

# Format Tables
table_freqmodel_dLight1.3b_AUC <- format_regsummary_table(AUCmodel_freq_dLight1.3b, 'dLight1.3b', "Frequency", "AUCModel")
table_freqmodel_dLight1.3b_AUC_aov <- format_anova_table(AUCmodel_freq_dLight1.3b_aov, "dLight1.3b", "Frequency", "AUCModel") 
table_freqmodel_dLight1.3b_AUC_reg <- format_reg_table(AUC_freq_dLight1.3b_regdata, "dLight1.3b", "Frequency", "AUCReg")

### Amplitude -----
#### Amplitude Model -----
# Mixed effects model
ampmodel_dLight1.3b <- glmmTMB(amp ~ InjType*Dose*BinF + (1|SubjectID), # random intercept by subject
                            data = subset(pkdataz, Virus=='dLight1.3b' & CondNum == 1), family = stats::gaussian, 
                            control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampmodel_dLight1.3b_aov <- Anova(ampmodel_dLight1.3b, type=2)
ampmodel_dLight1.3b_aov
summary(ampmodel_dLight1.3b)

ampmodel_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(ampmodel_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampmodel_dLight1.3b_r2

# Follow up tests - joint
ampmodel_dLight1.3b_emmeans <-emmeans(ampmodel_dLight1.3b, pairwise~InjType|Dose, adjust="none") # Find means
ampmodel_dLight1.3b_jt_Dose <- joint_tests(ampmodel_dLight1.3b_emmeans, by = "Dose")
ampmodel_dLight1.3b_jt_InjType <- joint_tests(ampmodel_dLight1.3b_emmeans, by = "InjType")

ampmodel_dLight1.3b_jt_Dose
ampmodel_dLight1.3b_jt_InjType

# Follow up tests - treat pairs
ampmodel_dLight1.3b_treatpairs <- emmeans(ampmodel_dLight1.3b, pairwise~InjType|Dose, adjust="none") # Pairs
ampmodel_dLight1.3b_treatcontrasts <- as.data.frame(ampmodel_dLight1.3b_treatpairs$contrasts) # Extract emmeans
ampmodel_dLight1.3b_treatemmeans <- ampmodel_dLight1.3b_treatpairs$emmeans # Extract emmeans
ampmodel_dLight1.3b_treatpairs_es <- as.data.frame(eff_size(ampmodel_dLight1.3b_treatemmeans, sigma = sigma(ampmodel_dLight1.3b), edf = df.residual(ampmodel_dLight1.3b))) # Effect Size
ampmodel_dLight1.3b_treatpairs_es <- ampmodel_dLight1.3b_treatpairs_es[, c("Dose", "contrast", "effect.size")]

ampmodel_dLight1.3b_treatcontrasts$d <- ampmodel_dLight1.3b_treatpairs_es$effect.size
ampmodel_dLight1.3b_treatcontrasts

# Follow up tests - dose pairs
ampmodel_dLight1.3b_dosepairs <- emmeans(ampmodel_dLight1.3b, pairwise~Dose|InjType, adjust="none") # Pairs
ampmodel_dLight1.3b_dosecontrasts <- as.data.frame(ampmodel_dLight1.3b_dosepairs$contrasts) # Extract emmeans
ampmodel_dLight1.3b_doseemmeans <- ampmodel_dLight1.3b_dosepairs$emmeans # Extract emmeans
ampmodel_dLight1.3b_dosepairs_es <- as.data.frame(eff_size(ampmodel_dLight1.3b_doseemmeans, sigma = sigma(ampmodel_dLight1.3b), edf = df.residual(ampmodel_dLight1.3b))) # Effect Size
ampmodel_dLight1.3b_dosepairs_es <- ampmodel_dLight1.3b_dosepairs_es[, c("InjType", "contrast", "effect.size")]

ampmodel_dLight1.3b_dosecontrasts$d <- ampmodel_dLight1.3b_dosepairs_es$effect.size
ampmodel_dLight1.3b_dosecontrasts

# Format Tables
table_ampmodel_dLight1.3b <- format_regsummary_table(ampmodel_dLight1.3b, 'dLight1.3b', "Amplitude", "StandardModel")
table_ampmodel_dLight1.3b_aov <- format_anova_table(ampmodel_dLight1.3b_aov, "dLight1.3b", "Amplitude", "StandardModel")
table_ampmodel_dLight1.3b_joint_Dose <- format_emm_table(ampmodel_dLight1.3b_jt_Dose, "dLight1.3b", "Amplitude", "Joint", "Dose")
table_ampmodel_dLight1.3b_joint_InjType <- format_emm_table(ampmodel_dLight1.3b_jt_InjType, "dLight1.3b", "Amplitude", "Joint", "InjType")
table_ampmodel_dLight1.3b_pairs_Dose <- format_emm_table(ampmodel_dLight1.3b_dosecontrasts, "dLight1.3b", "Amplitude", "Pairs", "Dose|InjType")
table_ampmodel_dLight1.3b_pairs_InjType <- format_emm_table(ampmodel_dLight1.3b_treatcontrasts, "dLight1.3b", "Amplitude", "Pairs", "InjType|Dose")

#### Mixed effects AR1 model
ampmodel_dLight1.3b_ar1 <- glmmTMB(amp ~ InjType*Dose*BinF + ar1(BinF + 0|SubjectID:InjType) + (1|SubjectID), # random intercept by subject
                                data = subset(pkdataz, Virus=='dLight1.3b' & CondNum == 1), family = stats::gaussian, 
                                control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampmodel_dLight1.3b_ar1_aov <- Anova(ampmodel_dLight1.3b_ar1, type=2)
ampmodel_dLight1.3b_ar1_aov
summary(ampmodel_dLight1.3b_ar1)

ampmodel_dLight1.3b_ar1_r2 <- as.data.frame(r.squaredGLMM(ampmodel_dLight1.3b_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampmodel_dLight1.3b_ar1_r2

ampmodel_dLight1.3b_ar1_AIC <- as.data.frame(AIC(ampmodel_dLight1.3b, ampmodel_dLight1.3b_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "ampmodel_dLight1.3b" ~ "StandardModel", ModelName == "ampmodel_dLight1.3b_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

ampmodel_dLight1.3b_ar1_AIC

# Follow up tests - joint
ampmodel_dLight1.3b_ar1_emmeans <-emmeans(ampmodel_dLight1.3b_ar1, pairwise~InjType|Dose, adjust="none") # Find means
ampmodel_dLight1.3b_ar1_jt_Dose <- joint_tests(ampmodel_dLight1.3b_ar1, by = "Dose")
ampmodel_dLight1.3b_ar1_jt_InjType <- joint_tests(ampmodel_dLight1.3b_ar1, by = "InjType")

ampmodel_dLight1.3b_ar1_jt_Dose
ampmodel_dLight1.3b_ar1_jt_InjType

# Follow up tests - treat pairs
ampmodel_dLight1.3b_ar1_treatpairs <- emmeans(ampmodel_dLight1.3b_ar1, pairwise~InjType|Dose, adjust="none") # Pairs
ampmodel_dLight1.3b_ar1_treatcontrasts <- as.data.frame(ampmodel_dLight1.3b_ar1_treatpairs$contrasts) # Extract emmeans
ampmodel_dLight1.3b_ar1_treatemmeans <- ampmodel_dLight1.3b_ar1_treatpairs$emmeans # Extract emmeans
ampmodel_dLight1.3b_ar1_treatpairs_es <- as.data.frame(eff_size(ampmodel_dLight1.3b_ar1_treatemmeans, sigma = sigma(ampmodel_dLight1.3b_ar1), edf = df.residual(ampmodel_dLight1.3b_ar1))) # Effect Size
ampmodel_dLight1.3b_ar1_treatpairs_es <- ampmodel_dLight1.3b_ar1_treatpairs_es[, c("Dose", "contrast", "effect.size")]

ampmodel_dLight1.3b_ar1_treatcontrasts$d <- ampmodel_dLight1.3b_ar1_treatpairs_es$effect.size
ampmodel_dLight1.3b_ar1_treatcontrasts

# Follow up tests - dose pairs
ampmodel_dLight1.3b_ar1_dosepairs <- emmeans(ampmodel_dLight1.3b_ar1, pairwise~Dose|InjType, adjust="none") # Pairs
ampmodel_dLight1.3b_ar1_dosecontrasts <- as.data.frame(ampmodel_dLight1.3b_ar1_dosepairs$contrasts) # Extract emmeans
ampmodel_dLight1.3b_ar1_doseemmeans <- ampmodel_dLight1.3b_ar1_dosepairs$emmeans # Extract emmeans
ampmodel_dLight1.3b_ar1_dosepairs_es <- as.data.frame(eff_size(ampmodel_dLight1.3b_ar1_doseemmeans, sigma = sigma(ampmodel_dLight1.3b_ar1), edf = df.residual(ampmodel_dLight1.3b_ar1))) # Effect Size
ampmodel_dLight1.3b_ar1_dosepairs_es <- ampmodel_dLight1.3b_ar1_dosepairs_es[, c("InjType", "contrast", "effect.size")]

ampmodel_dLight1.3b_ar1_dosecontrasts$d <- ampmodel_dLight1.3b_ar1_dosepairs_es$effect.size
ampmodel_dLight1.3b_ar1_dosecontrasts

# Format Tables
table_ampmodel_dLight1.3b_ar1 <- format_regsummary_table(ampmodel_dLight1.3b_ar1, 'dLight1.3b', "Amplitude", "AR1Model")
table_ampmodel_dLight1.3b_ar1_aov <- format_anova_table(ampmodel_dLight1.3b_ar1_aov, "dLight1.3b", "Amplitude", "AR1Model") 
table_ampmodel_dLight1.3b_ar1_joint_Dose <- format_emm_table(ampmodel_dLight1.3b_ar1_jt_Dose, "dLight1.3b", "Amplitude", "Joint", "Dose")
table_ampmodel_dLight1.3b_ar1_joint_InjType <- format_emm_table(ampmodel_dLight1.3b_ar1_jt_InjType, "dLight1.3b", "Amplitude", "Joint", "InjType")
table_ampmodel_dLight1.3b_ar1_pairs_Dose <- format_emm_table(ampmodel_dLight1.3b_ar1_dosecontrasts, "dLight1.3b", "Amplitude", "Pairs", "Dose|InjType")
table_ampmodel_dLight1.3b_ar1_pairs_InjType <- format_emm_table(ampmodel_dLight1.3b_ar1_treatcontrasts, "dLight1.3b", "Amplitude", "Pairs", "InjType|Dose")

#### Amplitude Change Model -----
# Mixed effects model
ampCmodel_dLight1.3b <- glmmTMB(ampC ~ Dose*BinF + (1|SubjectID), # random intercept by subject
                             data = subset(subjectmeansbin, Virus=='dLight1.3b' & CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                             control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampCmodel_dLight1.3b_aov <- Anova(ampCmodel_dLight1.3b, type=2)
ampCmodel_dLight1.3b_aov
summary(ampCmodel_dLight1.3b)

ampCmodel_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(ampCmodel_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampCmodel_dLight1.3b_r2

# Follow up tests - joint
ampCmodel_dLight1.3b_emmeans <-emmeans(ampCmodel_dLight1.3b, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
ampCmodel_dLight1.3b_dosepairs <- emmeans(ampCmodel_dLight1.3b, pairwise~Dose, adjust="none") # Pairs
ampCmodel_dLight1.3b_dosecontrasts <- as.data.frame(ampCmodel_dLight1.3b_dosepairs$contrasts) # Extract emmeans
ampCmodel_dLight1.3b_doseemmeans <- ampCmodel_dLight1.3b_dosepairs$emmeans # Extract emmeans
ampCmodel_dLight1.3b_dosepairs_es <- as.data.frame(eff_size(ampCmodel_dLight1.3b_doseemmeans, sigma = sigma(ampCmodel_dLight1.3b), edf = df.residual(ampCmodel_dLight1.3b))) # Effect Size
ampCmodel_dLight1.3b_dosepairs_es <- ampCmodel_dLight1.3b_dosepairs_es[, c("contrast", "effect.size")]

ampCmodel_dLight1.3b_dosecontrasts$d <- ampCmodel_dLight1.3b_dosepairs_es$effect.size
ampCmodel_dLight1.3b_dosecontrasts

# Format Tables
table_ampCmodel_dLight1.3b <- format_regsummary_table(ampCmodel_dLight1.3b, 'dLight1.3b', "Amplitude", "StandardModel")
table_ampCmodel_dLight1.3b_aov <- format_anova_table(ampCmodel_dLight1.3b_aov, "dLight1.3b", "Amplitude", "StandardModel")
table_ampCmodel_dLight1.3b_pairs_Dose <- format_emm_table(ampCmodel_dLight1.3b_dosecontrasts, "dLight1.3b", "Amplitude", "Pairs","Dose")

#### Mixed effects AR1 model
ampCmodel_dLight1.3b_ar1 <- glmmTMB(ampC ~ Dose*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                 data = subset(subjectmeansbin, Virus=='dLight1.3b' & CondNum == 1& InjType=='M'), family = stats::gaussian, 
                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampCmodel_dLight1.3b_ar1_aov <- Anova(ampCmodel_dLight1.3b_ar1, type=2)
ampCmodel_dLight1.3b_ar1_aov
summary(ampCmodel_dLight1.3b_ar1)

ampCmodel_dLight1.3b_ar1_r2 <- as.data.frame(r.squaredGLMM(ampCmodel_dLight1.3b_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampCmodel_dLight1.3b_ar1_r2

ampCmodel_dLight1.3b_ar1_AIC <- as.data.frame(AIC(ampCmodel_dLight1.3b, ampCmodel_dLight1.3b_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "ampCmodel_dLight1.3b" ~ "StandardModel", ModelName == "ampCmodel_dLight1.3b_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

ampCmodel_dLight1.3b_ar1_AIC

# Follow up tests - joint
ampCmodel_dLight1.3b_ar1_emmeans <-emmeans(ampCmodel_dLight1.3b_ar1, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
ampCmodel_dLight1.3b_ar1_dosepairs <- emmeans(ampCmodel_dLight1.3b_ar1, pairwise~Dose, adjust="none") # Pairs
ampCmodel_dLight1.3b_ar1_dosecontrasts <- as.data.frame(ampCmodel_dLight1.3b_ar1_dosepairs$contrasts) # Extract emmeans
ampCmodel_dLight1.3b_ar1_doseemmeans <- ampCmodel_dLight1.3b_ar1_dosepairs$emmeans # Extract emmeans
ampCmodel_dLight1.3b_ar1_dosepairs_es <- as.data.frame(eff_size(ampCmodel_dLight1.3b_ar1_doseemmeans, sigma = sigma(ampCmodel_dLight1.3b_ar1), edf = df.residual(ampCmodel_dLight1.3b_ar1))) # Effect Size
ampCmodel_dLight1.3b_ar1_dosepairs_es <- ampCmodel_dLight1.3b_ar1_dosepairs_es[, c("contrast", "effect.size")]

ampCmodel_dLight1.3b_ar1_dosecontrasts$d <- ampCmodel_dLight1.3b_ar1_dosepairs_es$effect.size
ampCmodel_dLight1.3b_ar1_dosecontrasts

# Format Tables
table_ampCmodel_dLight1.3b_ar1 <- format_regsummary_table(ampCmodel_dLight1.3b_ar1, 'dLight1.3b', "Amplitude", "AR1Model")
table_ampCmodel_dLight1.3b_ar1_aov <- format_anova_table(ampCmodel_dLight1.3b_ar1_aov, "dLight1.3b", "Amplitude", "AR1Model") 
table_ampCmodel_dLight1.3b_ar1_pairs_Dose <- format_emm_table(ampCmodel_dLight1.3b_ar1_dosecontrasts, "dLight1.3b", "Amplitude", "Pairs", "Dose")


#### Max Height Model -----
## Mixed effects model
MaxHeightmodel_amp_dLight1.3b <- glmmTMB(MaxHeight ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                      data = subjectrise_ampC_dLight1.3b, family = stats::gaussian)

MaxHeightmodel_amp_dLight1.3b_aov <- Anova(MaxHeightmodel_amp_dLight1.3b, type=2)
MaxHeightmodel_amp_dLight1.3b_aov
summary(MaxHeightmodel_amp_dLight1.3b)

MaxHeightmodel_amp_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(MaxHeightmodel_amp_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxHeightmodel_amp_dLight1.3b_r2

## Prepare fit data for plotting
MaxHeight_amp_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_dLight1.3b$DoseNum)))
MaxHeight_amp_dLight1.3b_regdata$DoseNum_c <- MaxHeight_amp_dLight1.3b_regdata$DoseNum - dosemean

MaxHeight_amp_dLight1.3b_regpred <- predict(MaxHeightmodel_amp_dLight1.3b, newdata = MaxHeight_amp_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxHeight_amp_dLight1.3b_regdata$fit <- MaxHeight_amp_dLight1.3b_regpred$fit
MaxHeight_amp_dLight1.3b_regdata$se <- MaxHeight_amp_dLight1.3b_regpred$se
MaxHeight_amp_dLight1.3b_regdata$lower <- MaxHeight_amp_dLight1.3b_regdata$fit - 1.96 * MaxHeight_amp_dLight1.3b_regdata$se
MaxHeight_amp_dLight1.3b_regdata$upper <- MaxHeight_amp_dLight1.3b_regdata$fit + 1.96 * MaxHeight_amp_dLight1.3b_regdata$se

MaxHeight_amp_dLight1.3b_regdata

# Format Tables
table_ampmodel_dLight1.3b_MaxHeight <- format_regsummary_table(MaxHeightmodel_amp_dLight1.3b, 'dLight1.3b', "Amplitude", "MaxHeightModel")
table_ampmodel_dLight1.3b_MaxHeight_aov <- format_anova_table(MaxHeightmodel_amp_dLight1.3b_aov, "dLight1.3b", "Amplitude", "MaxHeightModel") 
table_ampmodel_dLight1.3b_MaxHeight_reg <- format_reg_table(MaxHeight_amp_dLight1.3b_regdata, "dLight1.3b", "Amplitude", "MaxHeightReg")

#### Max Duration Model -----
## Mixed effects model
MaxDurationmodel_amp_dLight1.3b <- glmmTMB(MaxDuration~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                        data = subjectrise_ampC_dLight1.3b, family = stats::gaussian)

MaxDurationmodel_amp_dLight1.3b_aov<- Anova(MaxDurationmodel_amp_dLight1.3b, type=2)
MaxDurationmodel_amp_dLight1.3b_aov
summary(MaxDurationmodel_amp_dLight1.3b)

MaxDurationmodel_amp_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(MaxDurationmodel_amp_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxDurationmodel_amp_dLight1.3b_r2

## Prepare fit data for plotting
MaxDuration_amp_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_dLight1.3b$DoseNum)))
MaxDuration_amp_dLight1.3b_regdata$DoseNum_c <- MaxDuration_amp_dLight1.3b_regdata$DoseNum - dosemean

MaxDuration_amp_dLight1.3b_regpred <- predict(MaxDurationmodel_amp_dLight1.3b, newdata = MaxDuration_amp_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxDuration_amp_dLight1.3b_regdata$fit <- MaxDuration_amp_dLight1.3b_regpred$fit
MaxDuration_amp_dLight1.3b_regdata$se <- MaxDuration_amp_dLight1.3b_regpred$se
MaxDuration_amp_dLight1.3b_regdata$lower <- MaxDuration_amp_dLight1.3b_regdata$fit - 1.96 * MaxDuration_amp_dLight1.3b_regdata$se
MaxDuration_amp_dLight1.3b_regdata$upper <- MaxDuration_amp_dLight1.3b_regdata$fit + 1.96 * MaxDuration_amp_dLight1.3b_regdata$se

MaxDuration_amp_dLight1.3b_regdata

# Format Tables
table_ampmodel_dLight1.3b_MaxDuration <- format_regsummary_table(MaxDurationmodel_amp_dLight1.3b, 'dLight1.3b', "Amplitude", "MaxDurationModel")
table_ampmodel_dLight1.3b_MaxDuration_aov <- format_anova_table(MaxDurationmodel_amp_dLight1.3b_aov, "dLight1.3b", "Amplitude", "MaxDurationModel") 
table_ampmodel_dLight1.3b_MaxDuration_reg <- format_reg_table(MaxDuration_amp_dLight1.3b_regdata, "dLight1.3b", "Amplitude", "MaxDurationReg")

#### Rise Slope Model -----
## Mixed effects model
MaxSlopeLMmodel_amp_dLight1.3b <- glmmTMB(MaxSlopeLM ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                       data = subjectrise_ampC_dLight1.3b, family = stats::gaussian)

MaxSlopeLMmodel_amp_dLight1.3b_aov <- Anova(MaxSlopeLMmodel_amp_dLight1.3b, type=2)
MaxSlopeLMmodel_amp_dLight1.3b_aov
summary(MaxSlopeLMmodel_amp_dLight1.3b)

MaxSlopeLMmodel_amp_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(MaxSlopeLMmodel_amp_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxSlopeLMmodel_amp_dLight1.3b_r2

## Prepare fit data for plotting
MaxSlopeLM_amp_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_dLight1.3b$DoseNum)))
MaxSlopeLM_amp_dLight1.3b_regdata$DoseNum_c <- MaxSlopeLM_amp_dLight1.3b_regdata$DoseNum - dosemean

MaxSlopeLM_amp_dLight1.3b_regpred <- predict(MaxSlopeLMmodel_amp_dLight1.3b, newdata = MaxSlopeLM_amp_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxSlopeLM_amp_dLight1.3b_regdata$fit <- MaxSlopeLM_amp_dLight1.3b_regpred$fit
MaxSlopeLM_amp_dLight1.3b_regdata$se <- MaxSlopeLM_amp_dLight1.3b_regpred$se
MaxSlopeLM_amp_dLight1.3b_regdata$lower <- MaxSlopeLM_amp_dLight1.3b_regdata$fit - 1.96 * MaxSlopeLM_amp_dLight1.3b_regdata$se
MaxSlopeLM_amp_dLight1.3b_regdata$upper <- MaxSlopeLM_amp_dLight1.3b_regdata$fit + 1.96 * MaxSlopeLM_amp_dLight1.3b_regdata$se

MaxSlopeLM_amp_dLight1.3b_regdata

# Format Tables
table_ampmodel_dLight1.3b_MaxSlopeLM <- format_regsummary_table(MaxSlopeLMmodel_amp_dLight1.3b, 'dLight1.3b', "Amplitude", "MaxSlopeLMModel")
table_ampmodel_dLight1.3b_MaxSlopeLM_aov <- format_anova_table(MaxSlopeLMmodel_amp_dLight1.3b_aov, "dLight1.3b", "Amplitude", "MaxSlopeLMModel") 
table_ampmodel_dLight1.3b_MaxSlopeLM_reg <- format_reg_table(MaxSlopeLM_amp_dLight1.3b_regdata, "dLight1.3b", "Amplitude", "MaxSlopeLMReg")

#### AUC Model -----
## Mixed effects model
AUCmodel_amp_dLight1.3b <- glmmTMB(AUC_POSTINJ ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                data = subset(subjectrise_ampC_dLight1.3b, Virus=='dLight1.3b'&InjType=='M'), family = stats::gaussian)

AUCmodel_amp_dLight1.3b_aov <- Anova(AUCmodel_amp_dLight1.3b, type=2)
AUCmodel_amp_dLight1.3b_aov
summary(AUCmodel_amp_dLight1.3b)

AUCmodel_amp_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(AUCmodel_amp_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCmodel_amp_dLight1.3b_r2

## Prepare fit data for plotting
AUC_amp_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_dLight1.3b$DoseNum)))
AUC_amp_dLight1.3b_regdata$DoseNum_c <- AUC_amp_dLight1.3b_regdata$DoseNum - dosemean

AUC_amp_dLight1.3b_regpred <- predict(AUCmodel_amp_dLight1.3b, newdata = AUC_amp_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
AUC_amp_dLight1.3b_regdata$fit <- AUC_amp_dLight1.3b_regpred$fit
AUC_amp_dLight1.3b_regdata$se <- AUC_amp_dLight1.3b_regpred$se
AUC_amp_dLight1.3b_regdata$lower <- AUC_amp_dLight1.3b_regdata$fit - 1.96 * AUC_amp_dLight1.3b_regdata$se
AUC_amp_dLight1.3b_regdata$upper <- AUC_amp_dLight1.3b_regdata$fit + 1.96 * AUC_amp_dLight1.3b_regdata$se

AUC_amp_dLight1.3b_regdata

# Format Tables
table_ampmodel_dLight1.3b_AUC <- format_regsummary_table(AUCmodel_amp_dLight1.3b, 'dLight1.3b', "Amplitude", "AUCModel")
table_ampmodel_dLight1.3b_AUC_aov <- format_anova_table(AUCmodel_amp_dLight1.3b_aov, "dLight1.3b", "Amplitude", "AUCModel") 
table_ampmodel_dLight1.3b_AUC_reg <- format_reg_table(AUC_amp_dLight1.3b_regdata, "dLight1.3b", "Amplitude", "AUCReg")

### AUCwindow -----
#### AUCwindow Model -----
# Mixed effects model
AUCwindowmodel_dLight1.3b <- glmmTMB(AUCwindow ~ InjType*Dose*BinF + (1|SubjectID), # random intercept by subject
                                  data = subset(pkdataz, Virus=='dLight1.3b' & CondNum == 1), family = stats::gaussian, 
                                  control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowmodel_dLight1.3b_aov <- Anova(AUCwindowmodel_dLight1.3b, type=2)
AUCwindowmodel_dLight1.3b_aov
summary(AUCwindowmodel_dLight1.3b)

AUCwindowmodel_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(AUCwindowmodel_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowmodel_dLight1.3b_r2

# Follow up tests - joint
AUCwindowmodel_dLight1.3b_emmeans <-emmeans(AUCwindowmodel_dLight1.3b, pairwise~InjType|Dose, adjust="none") # Find means
AUCwindowmodel_dLight1.3b_jt_Dose <- joint_tests(AUCwindowmodel_dLight1.3b_emmeans, by = "Dose")
AUCwindowmodel_dLight1.3b_jt_InjType <- joint_tests(AUCwindowmodel_dLight1.3b_emmeans, by = "InjType")

AUCwindowmodel_dLight1.3b_jt_Dose
AUCwindowmodel_dLight1.3b_jt_InjType

# Follow up tests - treat pairs
AUCwindowmodel_dLight1.3b_treatpairs <- emmeans(AUCwindowmodel_dLight1.3b, pairwise~InjType|Dose, adjust="none") # Pairs
AUCwindowmodel_dLight1.3b_treatcontrasts <- as.data.frame(AUCwindowmodel_dLight1.3b_treatpairs$contrasts) # Extract emmeans
AUCwindowmodel_dLight1.3b_treatemmeans <- AUCwindowmodel_dLight1.3b_treatpairs$emmeans # Extract emmeans
AUCwindowmodel_dLight1.3b_treatpairs_es <- as.data.frame(eff_size(AUCwindowmodel_dLight1.3b_treatemmeans, sigma = sigma(AUCwindowmodel_dLight1.3b), edf = df.residual(AUCwindowmodel_dLight1.3b))) # Effect Size
AUCwindowmodel_dLight1.3b_treatpairs_es <- AUCwindowmodel_dLight1.3b_treatpairs_es[, c("Dose", "contrast", "effect.size")]

AUCwindowmodel_dLight1.3b_treatcontrasts$d <- AUCwindowmodel_dLight1.3b_treatpairs_es$effect.size
AUCwindowmodel_dLight1.3b_treatcontrasts

# Follow up tests - dose pairs
AUCwindowmodel_dLight1.3b_dosepairs <- emmeans(AUCwindowmodel_dLight1.3b, pairwise~Dose|InjType, adjust="none") # Pairs
AUCwindowmodel_dLight1.3b_dosecontrasts <- as.data.frame(AUCwindowmodel_dLight1.3b_dosepairs$contrasts) # Extract emmeans
AUCwindowmodel_dLight1.3b_doseemmeans <- AUCwindowmodel_dLight1.3b_dosepairs$emmeans # Extract emmeans
AUCwindowmodel_dLight1.3b_dosepairs_es <- as.data.frame(eff_size(AUCwindowmodel_dLight1.3b_doseemmeans, sigma = sigma(AUCwindowmodel_dLight1.3b), edf = df.residual(AUCwindowmodel_dLight1.3b))) # Effect Size
AUCwindowmodel_dLight1.3b_dosepairs_es <- AUCwindowmodel_dLight1.3b_dosepairs_es[, c("InjType", "contrast", "effect.size")]

AUCwindowmodel_dLight1.3b_dosecontrasts$d <- AUCwindowmodel_dLight1.3b_dosepairs_es$effect.size
AUCwindowmodel_dLight1.3b_dosecontrasts

# Format Tables
table_AUCwindowmodel_dLight1.3b <- format_regsummary_table(AUCwindowmodel_dLight1.3b, 'dLight1.3b', "AUCwindow", "StandardModel")
table_AUCwindowmodel_dLight1.3b_aov <- format_anova_table(AUCwindowmodel_dLight1.3b_aov, "dLight1.3b", "AUCwindow", "StandardModel")
table_AUCwindowmodel_dLight1.3b_joint_Dose <- format_emm_table(AUCwindowmodel_dLight1.3b_jt_Dose, "dLight1.3b", "AUCwindow", "Joint", "Dose")
table_AUCwindowmodel_dLight1.3b_joint_InjType <- format_emm_table(AUCwindowmodel_dLight1.3b_jt_InjType, "dLight1.3b", "AUCwindow", "Joint", "InjType")
table_AUCwindowmodel_dLight1.3b_pairs_Dose <- format_emm_table(AUCwindowmodel_dLight1.3b_dosecontrasts, "dLight1.3b", "AUCwindow", "Pairs", "Dose|InjType")
table_AUCwindowmodel_dLight1.3b_pairs_InjType <- format_emm_table(AUCwindowmodel_dLight1.3b_treatcontrasts, "dLight1.3b", "AUCwindow", "Pairs", "InjType|Dose")

#### Mixed effects AR1 model
AUCwindowmodel_dLight1.3b_ar1 <- glmmTMB(AUCwindow ~ InjType*Dose*BinF + ar1(BinF + 0|SubjectID:InjType) + (1|SubjectID), # random intercept by subject
                                      data = subset(pkdataz, Virus=='dLight1.3b' & CondNum == 1), family = stats::gaussian, 
                                      control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowmodel_dLight1.3b_ar1_aov <- Anova(AUCwindowmodel_dLight1.3b_ar1, type=2)
AUCwindowmodel_dLight1.3b_ar1_aov
summary(AUCwindowmodel_dLight1.3b_ar1)

AUCwindowmodel_dLight1.3b_ar1_r2 <- as.data.frame(r.squaredGLMM(AUCwindowmodel_dLight1.3b_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowmodel_dLight1.3b_ar1_r2

AUCwindowmodel_dLight1.3b_ar1_AIC <- as.data.frame(AIC(AUCwindowmodel_dLight1.3b, AUCwindowmodel_dLight1.3b_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "AUCwindowmodel_dLight1.3b" ~ "StandardModel", ModelName == "AUCwindowmodel_dLight1.3b_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

AUCwindowmodel_dLight1.3b_ar1_AIC

# Follow up tests - joint
AUCwindowmodel_dLight1.3b_ar1_emmeans <-emmeans(AUCwindowmodel_dLight1.3b_ar1, pairwise~InjType|Dose, adjust="none") # Find means
AUCwindowmodel_dLight1.3b_ar1_jt_Dose <- joint_tests(AUCwindowmodel_dLight1.3b_ar1, by = "Dose")
AUCwindowmodel_dLight1.3b_ar1_jt_InjType <- joint_tests(AUCwindowmodel_dLight1.3b_ar1, by = "InjType")

AUCwindowmodel_dLight1.3b_ar1_jt_Dose
AUCwindowmodel_dLight1.3b_ar1_jt_InjType

# Follow up tests - treat pairs
AUCwindowmodel_dLight1.3b_ar1_treatpairs <- emmeans(AUCwindowmodel_dLight1.3b_ar1, pairwise~InjType|Dose, adjust="none") # Pairs
AUCwindowmodel_dLight1.3b_ar1_treatcontrasts <- as.data.frame(AUCwindowmodel_dLight1.3b_ar1_treatpairs$contrasts) # Extract emmeans
AUCwindowmodel_dLight1.3b_ar1_treatemmeans <- AUCwindowmodel_dLight1.3b_ar1_treatpairs$emmeans # Extract emmeans
AUCwindowmodel_dLight1.3b_ar1_treatpairs_es <- as.data.frame(eff_size(AUCwindowmodel_dLight1.3b_ar1_treatemmeans, sigma = sigma(AUCwindowmodel_dLight1.3b_ar1), edf = df.residual(AUCwindowmodel_dLight1.3b_ar1))) # Effect Size
AUCwindowmodel_dLight1.3b_ar1_treatpairs_es <- AUCwindowmodel_dLight1.3b_ar1_treatpairs_es[, c("Dose", "contrast", "effect.size")]

AUCwindowmodel_dLight1.3b_ar1_treatcontrasts$d <- AUCwindowmodel_dLight1.3b_ar1_treatpairs_es$effect.size
AUCwindowmodel_dLight1.3b_ar1_treatcontrasts

# Follow up tests - dose pairs
AUCwindowmodel_dLight1.3b_ar1_dosepairs <- emmeans(AUCwindowmodel_dLight1.3b_ar1, pairwise~Dose|InjType, adjust="none") # Pairs
AUCwindowmodel_dLight1.3b_ar1_dosecontrasts <- as.data.frame(AUCwindowmodel_dLight1.3b_ar1_dosepairs$contrasts) # Extract emmeans
AUCwindowmodel_dLight1.3b_ar1_doseemmeans <- AUCwindowmodel_dLight1.3b_ar1_dosepairs$emmeans # Extract emmeans
AUCwindowmodel_dLight1.3b_ar1_dosepairs_es <- as.data.frame(eff_size(AUCwindowmodel_dLight1.3b_ar1_doseemmeans, sigma = sigma(AUCwindowmodel_dLight1.3b_ar1), edf = df.residual(AUCwindowmodel_dLight1.3b_ar1))) # Effect Size
AUCwindowmodel_dLight1.3b_ar1_dosepairs_es <- AUCwindowmodel_dLight1.3b_ar1_dosepairs_es[, c("InjType", "contrast", "effect.size")]

AUCwindowmodel_dLight1.3b_ar1_dosecontrasts$d <- AUCwindowmodel_dLight1.3b_ar1_dosepairs_es$effect.size
AUCwindowmodel_dLight1.3b_ar1_dosecontrasts

# Format Tables
table_AUCwindowmodel_dLight1.3b_ar1 <- format_regsummary_table(AUCwindowmodel_dLight1.3b_ar1, 'dLight1.3b', "AUCwindow", "AR1Model")
table_AUCwindowmodel_dLight1.3b_ar1_aov <- format_anova_table(AUCwindowmodel_dLight1.3b_ar1_aov, "dLight1.3b", "AUCwindow", "AR1Model") 
table_AUCwindowmodel_dLight1.3b_ar1_joint_Dose <- format_emm_table(AUCwindowmodel_dLight1.3b_ar1_jt_Dose, "dLight1.3b", "AUCwindow", "Joint", "Dose")
table_AUCwindowmodel_dLight1.3b_ar1_joint_InjType <- format_emm_table(AUCwindowmodel_dLight1.3b_ar1_jt_InjType, "dLight1.3b", "AUCwindow", "Joint", "InjType")
table_AUCwindowmodel_dLight1.3b_ar1_pairs_Dose <- format_emm_table(AUCwindowmodel_dLight1.3b_ar1_dosecontrasts, "dLight1.3b", "AUCwindow", "Pairs", "Dose|InjType")
table_AUCwindowmodel_dLight1.3b_ar1_pairs_InjType <- format_emm_table(AUCwindowmodel_dLight1.3b_ar1_treatcontrasts, "dLight1.3b", "AUCwindow", "Pairs", "InjType|Dose")

#### AUCwindow Change Model -----
# Mixed effects model
AUCwindowCmodel_dLight1.3b <- glmmTMB(AUCwindowC ~ Dose*BinF + (1|SubjectID), # random intercept by subject
                                   data = subset(subjectmeansbin, Virus=='dLight1.3b' & CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                   control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowCmodel_dLight1.3b_aov <- Anova(AUCwindowCmodel_dLight1.3b, type=2)
AUCwindowCmodel_dLight1.3b_aov
summary(AUCwindowCmodel_dLight1.3b)

AUCwindowCmodel_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(AUCwindowCmodel_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowCmodel_dLight1.3b_r2

# Follow up tests - joint
AUCwindowCmodel_dLight1.3b_emmeans <-emmeans(AUCwindowCmodel_dLight1.3b, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
AUCwindowCmodel_dLight1.3b_dosepairs <- emmeans(AUCwindowCmodel_dLight1.3b, pairwise~Dose, adjust="none") # Pairs
AUCwindowCmodel_dLight1.3b_dosecontrasts <- as.data.frame(AUCwindowCmodel_dLight1.3b_dosepairs$contrasts) # Extract emmeans
AUCwindowCmodel_dLight1.3b_doseemmeans <- AUCwindowCmodel_dLight1.3b_dosepairs$emmeans # Extract emmeans
AUCwindowCmodel_dLight1.3b_dosepairs_es <- as.data.frame(eff_size(AUCwindowCmodel_dLight1.3b_doseemmeans, sigma = sigma(AUCwindowCmodel_dLight1.3b), edf = df.residual(AUCwindowCmodel_dLight1.3b))) # Effect Size
AUCwindowCmodel_dLight1.3b_dosepairs_es <- AUCwindowCmodel_dLight1.3b_dosepairs_es[, c("contrast", "effect.size")]

AUCwindowCmodel_dLight1.3b_dosecontrasts$d <- AUCwindowCmodel_dLight1.3b_dosepairs_es$effect.size
AUCwindowCmodel_dLight1.3b_dosecontrasts

# Format Tables
table_AUCwindowCmodel_dLight1.3b <- format_regsummary_table(AUCwindowCmodel_dLight1.3b, 'dLight1.3b', "Amplitude", "StandardModel")
table_AUCwindowCmodel_dLight1.3b_aov <- format_anova_table(AUCwindowCmodel_dLight1.3b_aov, "dLight1.3b", "Amplitude", "StandardModel")
table_AUCwindowCmodel_dLight1.3b_pairs_Dose <- format_emm_table(AUCwindowCmodel_dLight1.3b_dosecontrasts, "dLight1.3b", "Amplitude", "Pairs","Dose")

#### Mixed effects AR1 model
AUCwindowCmodel_dLight1.3b_ar1 <- glmmTMB(AUCwindowC ~ Dose*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                       data = subset(subjectmeansbin, Virus=='dLight1.3b' & CondNum == 1& InjType=='M'), family = stats::gaussian, 
                                       control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowCmodel_dLight1.3b_ar1_aov <- Anova(AUCwindowCmodel_dLight1.3b_ar1, type=2)
AUCwindowCmodel_dLight1.3b_ar1_aov
summary(AUCwindowCmodel_dLight1.3b_ar1)

AUCwindowCmodel_dLight1.3b_ar1_r2 <- as.data.frame(r.squaredGLMM(AUCwindowCmodel_dLight1.3b_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowCmodel_dLight1.3b_ar1_r2

AUCwindowCmodel_dLight1.3b_ar1_AIC <- as.data.frame(AIC(AUCwindowCmodel_dLight1.3b, AUCwindowCmodel_dLight1.3b_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "AUCwindowCmodel_dLight1.3b" ~ "StandardModel", ModelName == "AUCwindowCmodel_dLight1.3b_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

AUCwindowCmodel_dLight1.3b_ar1_AIC

# Follow up tests - joint
AUCwindowCmodel_dLight1.3b_ar1_emmeans <-emmeans(AUCwindowCmodel_dLight1.3b_ar1, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
AUCwindowCmodel_dLight1.3b_ar1_dosepairs <- emmeans(AUCwindowCmodel_dLight1.3b_ar1, pairwise~Dose, adjust="none") # Pairs
AUCwindowCmodel_dLight1.3b_ar1_dosecontrasts <- as.data.frame(AUCwindowCmodel_dLight1.3b_ar1_dosepairs$contrasts) # Extract emmeans
AUCwindowCmodel_dLight1.3b_ar1_doseemmeans <- AUCwindowCmodel_dLight1.3b_ar1_dosepairs$emmeans # Extract emmeans
AUCwindowCmodel_dLight1.3b_ar1_dosepairs_es <- as.data.frame(eff_size(AUCwindowCmodel_dLight1.3b_ar1_doseemmeans, sigma = sigma(AUCwindowCmodel_dLight1.3b_ar1), edf = df.residual(AUCwindowCmodel_dLight1.3b_ar1))) # Effect Size
AUCwindowCmodel_dLight1.3b_ar1_dosepairs_es <- AUCwindowCmodel_dLight1.3b_ar1_dosepairs_es[, c("contrast", "effect.size")]

AUCwindowCmodel_dLight1.3b_ar1_dosecontrasts$d <- AUCwindowCmodel_dLight1.3b_ar1_dosepairs_es$effect.size
AUCwindowCmodel_dLight1.3b_ar1_dosecontrasts

# Format Tables
table_AUCwindowCmodel_dLight1.3b_ar1 <- format_regsummary_table(AUCwindowCmodel_dLight1.3b_ar1, 'dLight1.3b', "Amplitude", "AR1Model")
table_AUCwindowCmodel_dLight1.3b_ar1_aov <- format_anova_table(AUCwindowCmodel_dLight1.3b_ar1_aov, "dLight1.3b", "Amplitude", "AR1Model") 
table_AUCwindowCmodel_dLight1.3b_ar1_pairs_Dose <- format_emm_table(AUCwindowCmodel_dLight1.3b_ar1_dosecontrasts, "dLight1.3b", "Amplitude", "Pairs", "Dose")

#### Max Height Model -----
## Mixed effects model
MaxHeightmodel_AUCwindow_dLight1.3b <- glmmTMB(MaxHeight ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                            data = subjectrise_AUCwindowC_dLight1.3b, family = stats::gaussian)

MaxHeightmodel_AUCwindow_dLight1.3b_aov <- Anova(MaxHeightmodel_AUCwindow_dLight1.3b, type=2)
MaxHeightmodel_AUCwindow_dLight1.3b_aov
summary(MaxHeightmodel_AUCwindow_dLight1.3b)

MaxHeightmodel_AUCwindow_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(MaxHeightmodel_AUCwindow_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxHeightmodel_AUCwindow_dLight1.3b_r2

## Prepare fit data for plotting
MaxHeight_AUCwindow_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_dLight1.3b$DoseNum)))
MaxHeight_AUCwindow_dLight1.3b_regdata$DoseNum_c <- MaxHeight_AUCwindow_dLight1.3b_regdata$DoseNum - dosemean

MaxHeight_AUCwindow_dLight1.3b_regpred <- predict(MaxHeightmodel_AUCwindow_dLight1.3b, newdata = MaxHeight_AUCwindow_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxHeight_AUCwindow_dLight1.3b_regdata$fit <- MaxHeight_AUCwindow_dLight1.3b_regpred$fit
MaxHeight_AUCwindow_dLight1.3b_regdata$se <- MaxHeight_AUCwindow_dLight1.3b_regpred$se
MaxHeight_AUCwindow_dLight1.3b_regdata$lower <- MaxHeight_AUCwindow_dLight1.3b_regdata$fit - 1.96 * MaxHeight_AUCwindow_dLight1.3b_regdata$se
MaxHeight_AUCwindow_dLight1.3b_regdata$upper <- MaxHeight_AUCwindow_dLight1.3b_regdata$fit + 1.96 * MaxHeight_AUCwindow_dLight1.3b_regdata$se

MaxHeight_AUCwindow_dLight1.3b_regdata

# Format Tables
table_AUCwindowmodel_dLight1.3b_MaxHeight <- format_regsummary_table(MaxHeightmodel_AUCwindow_dLight1.3b, 'dLight1.3b', "AUCwindow", "MaxHeightModel")
table_AUCwindowmodel_dLight1.3b_MaxHeight_aov <- format_anova_table(MaxHeightmodel_AUCwindow_dLight1.3b_aov, "dLight1.3b", "AUCwindow", "MaxHeightModel") 
table_AUCwindowmodel_dLight1.3b_MaxHeight_reg <- format_reg_table(MaxHeight_AUCwindow_dLight1.3b_regdata, "dLight1.3b", "AUCwindow", "MaxHeightReg")

#### Max Duration Model -----
## Mixed effects model
MaxDurationmodel_AUCwindow_dLight1.3b <- glmmTMB(MaxDuration~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                              data = subjectrise_AUCwindowC_dLight1.3b, family = stats::gaussian)

MaxDurationmodel_AUCwindow_dLight1.3b_aov<- Anova(MaxDurationmodel_AUCwindow_dLight1.3b, type=2)
MaxDurationmodel_AUCwindow_dLight1.3b_aov
summary(MaxDurationmodel_AUCwindow_dLight1.3b)

MaxDurationmodel_AUCwindow_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(MaxDurationmodel_AUCwindow_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxDurationmodel_AUCwindow_dLight1.3b_r2

## Prepare fit data for plotting
MaxDuration_AUCwindow_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_dLight1.3b$DoseNum)))
MaxDuration_AUCwindow_dLight1.3b_regdata$DoseNum_c <- MaxDuration_AUCwindow_dLight1.3b_regdata$DoseNum - dosemean

MaxDuration_AUCwindow_dLight1.3b_regpred <- predict(MaxDurationmodel_AUCwindow_dLight1.3b, newdata = MaxDuration_AUCwindow_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxDuration_AUCwindow_dLight1.3b_regdata$fit <- MaxDuration_AUCwindow_dLight1.3b_regpred$fit
MaxDuration_AUCwindow_dLight1.3b_regdata$se <- MaxDuration_AUCwindow_dLight1.3b_regpred$se
MaxDuration_AUCwindow_dLight1.3b_regdata$lower <- MaxDuration_AUCwindow_dLight1.3b_regdata$fit - 1.96 * MaxDuration_AUCwindow_dLight1.3b_regdata$se
MaxDuration_AUCwindow_dLight1.3b_regdata$upper <- MaxDuration_AUCwindow_dLight1.3b_regdata$fit + 1.96 * MaxDuration_AUCwindow_dLight1.3b_regdata$se

MaxDuration_AUCwindow_dLight1.3b_regdata

# Format Tables
table_AUCwindowmodel_dLight1.3b_MaxDuration <- format_regsummary_table(MaxDurationmodel_AUCwindow_dLight1.3b, 'dLight1.3b', "AUCwindow", "MaxDurationModel")
table_AUCwindowmodel_dLight1.3b_MaxDuration_aov <- format_anova_table(MaxDurationmodel_AUCwindow_dLight1.3b_aov, "dLight1.3b", "AUCwindow", "MaxDurationModel") 
table_AUCwindowmodel_dLight1.3b_MaxDuration_reg <- format_reg_table(MaxDuration_AUCwindow_dLight1.3b_regdata, "dLight1.3b", "AUCwindow", "MaxDurationReg")

#### Rise Slope Model -----
## Mixed effects model
MaxSlopeLMmodel_AUCwindow_dLight1.3b <- glmmTMB(MaxSlopeLM ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                             data = subjectrise_AUCwindowC_dLight1.3b, family = stats::gaussian)

MaxSlopeLMmodel_AUCwindow_dLight1.3b_aov <- Anova(MaxSlopeLMmodel_AUCwindow_dLight1.3b, type=2)
MaxSlopeLMmodel_AUCwindow_dLight1.3b_aov
summary(MaxSlopeLMmodel_AUCwindow_dLight1.3b)

MaxSlopeLMmodel_AUCwindow_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(MaxSlopeLMmodel_AUCwindow_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxSlopeLMmodel_AUCwindow_dLight1.3b_r2

## Prepare fit data for plotting
MaxSlopeLM_AUCwindow_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_dLight1.3b$DoseNum)))
MaxSlopeLM_AUCwindow_dLight1.3b_regdata$DoseNum_c <- MaxSlopeLM_AUCwindow_dLight1.3b_regdata$DoseNum - dosemean

MaxSlopeLM_AUCwindow_dLight1.3b_regpred <- predict(MaxSlopeLMmodel_AUCwindow_dLight1.3b, newdata = MaxSlopeLM_AUCwindow_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxSlopeLM_AUCwindow_dLight1.3b_regdata$fit <- MaxSlopeLM_AUCwindow_dLight1.3b_regpred$fit
MaxSlopeLM_AUCwindow_dLight1.3b_regdata$se <- MaxSlopeLM_AUCwindow_dLight1.3b_regpred$se
MaxSlopeLM_AUCwindow_dLight1.3b_regdata$lower <- MaxSlopeLM_AUCwindow_dLight1.3b_regdata$fit - 1.96 * MaxSlopeLM_AUCwindow_dLight1.3b_regdata$se
MaxSlopeLM_AUCwindow_dLight1.3b_regdata$upper <- MaxSlopeLM_AUCwindow_dLight1.3b_regdata$fit + 1.96 * MaxSlopeLM_AUCwindow_dLight1.3b_regdata$se

MaxSlopeLM_AUCwindow_dLight1.3b_regdata

# Format Tables
table_AUCwindowmodel_dLight1.3b_MaxSlopeLM <- format_regsummary_table(MaxSlopeLMmodel_AUCwindow_dLight1.3b, 'dLight1.3b', "AUCwindow", "MaxSlopeLMModel")
table_AUCwindowmodel_dLight1.3b_MaxSlopeLM_aov <- format_anova_table(MaxSlopeLMmodel_AUCwindow_dLight1.3b_aov, "dLight1.3b", "AUCwindow", "MaxSlopeLMModel") 
table_AUCwindowmodel_dLight1.3b_MaxSlopeLM_reg <- format_reg_table(MaxSlopeLM_AUCwindow_dLight1.3b_regdata, "dLight1.3b", "AUCwindow", "MaxSlopeLMReg")

#### AUC Model -----
## Mixed effects model
AUCmodel_AUCwindow_dLight1.3b <- glmmTMB(AUC_POSTINJ ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                      data = subset(subjectrise_AUCwindowC_dLight1.3b, Virus=='dLight1.3b'&InjType=='M'), family = stats::gaussian)

AUCmodel_AUCwindow_dLight1.3b_aov <- Anova(AUCmodel_AUCwindow_dLight1.3b, type=2)
AUCmodel_AUCwindow_dLight1.3b_aov
summary(AUCmodel_AUCwindow_dLight1.3b)

AUCmodel_AUCwindow_dLight1.3b_r2 <- as.data.frame(r.squaredGLMM(AUCmodel_AUCwindow_dLight1.3b)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCmodel_AUCwindow_dLight1.3b_r2

## Prepare fit data for plotting
AUC_AUCwindow_dLight1.3b_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_dLight1.3b$DoseNum)))
AUC_AUCwindow_dLight1.3b_regdata$DoseNum_c <- AUC_AUCwindow_dLight1.3b_regdata$DoseNum - dosemean

AUC_AUCwindow_dLight1.3b_regpred <- predict(AUCmodel_AUCwindow_dLight1.3b, newdata = AUC_AUCwindow_dLight1.3b_regdata, se.fit = TRUE, re.form = NA, type = "response")
AUC_AUCwindow_dLight1.3b_regdata$fit <- AUC_AUCwindow_dLight1.3b_regpred$fit
AUC_AUCwindow_dLight1.3b_regdata$se <- AUC_AUCwindow_dLight1.3b_regpred$se
AUC_AUCwindow_dLight1.3b_regdata$lower <- AUC_AUCwindow_dLight1.3b_regdata$fit - 1.96 * AUC_AUCwindow_dLight1.3b_regdata$se
AUC_AUCwindow_dLight1.3b_regdata$upper <- AUC_AUCwindow_dLight1.3b_regdata$fit + 1.96 * AUC_AUCwindow_dLight1.3b_regdata$se

AUC_AUCwindow_dLight1.3b_regdata

# Format Tables
table_AUCwindowmodel_dLight1.3b_AUC <- format_regsummary_table(AUCmodel_AUCwindow_dLight1.3b, 'dLight1.3b', "AUCwindow", "AUCModel")
table_AUCwindowmodel_dLight1.3b_AUC_aov <- format_anova_table(AUCmodel_AUCwindow_dLight1.3b_aov, "dLight1.3b", "AUCwindow", "AUCModel") 
table_AUCwindowmodel_dLight1.3b_AUC_reg <- format_reg_table(AUC_AUCwindow_dLight1.3b_regdata, "dLight1.3b", "AUCwindow", "AUCReg")

## EXPORT STATS -----
### Prepare Means -----
# Frequency
treatmeanstable_freq_dLight1.3b <- treatmeans %>% filter(Virus=='dLight1.3b') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         pkspermin_mean, pkspermin_sd, pkspermin_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansCtable_freq_dLight1.3b <- treatmeans %>% filter(InjType!='S', Virus=='dLight1.3b') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         pksperminC_mean, pksperminC_sd, pksperminC_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansSNRtable_freq_dLight1.3b <- treatmeans %>% filter(InjType!='S', Virus=='dLight1.3b') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         pksperminSNR_virus_mean, pksperminSNR_virus_sd, pksperminSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Amplitude
treatmeanstable_amp_dLight1.3b <- treatmeans %>% filter(Virus=='dLight1.3b') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         amp_mean, amp_sd, amp_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansCtable_amp_dLight1.3b <- treatmeans %>% filter(InjType!='S', Virus=='dLight1.3b') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         ampC_mean, ampC_sd, ampC_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansSNRtable_amp_dLight1.3b <- treatmeans %>% filter(InjType!='S', Virus=='dLight1.3b') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         ampSNR_virus_mean, ampSNR_virus_sd, ampSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# AUC Window
treatmeanstable_AUCwindow_dLight1.3b <- treatmeans %>% filter(Virus=='dLight1.3b') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         AUCwindow_mean, AUCwindow_sd, AUCwindow_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansCtable_AUCwindow_dLight1.3b <- treatmeans %>% filter(InjType!='S', Virus=='dLight1.3b') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         AUCwindowC_mean, AUCwindowC_sd, AUCwindowC_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansSNRtable_AUCwindow_dLight1.3b <- treatmeans %>% filter(InjType!='S', Virus=='dLight1.3b') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         AUCwindowSNR_virus_mean, AUCwindowSNR_virus_sd, AUCwindowSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))


### Frequency Export -----
wb_dLight1.3b_freq <- createWorkbook()
wb_dLight1.3b_freq_ar1 <- createWorkbook()
wb_dLight1.3b_freqC <- createWorkbook()
wb_dLight1.3b_freqC_ar1 <- createWorkbook()
wb_dLight1.3b_freq_reg <- createWorkbook()

## Main results for dLight1.3b
# Standard Models
write_nice_sheet(wb_dLight1.3b_freq, "MeansRaw", treatmeanstable_freq_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_freq, "ModelSummary", table_freqmodel_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_freq, "ModelAnova", table_freqmodel_dLight1.3b_aov)
write_nice_sheet(wb_dLight1.3b_freq, "ModelR2", freqmodel_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_freq, "Joint_InjType",  table_freqmodel_dLight1.3b_joint_InjType)
write_nice_sheet(wb_dLight1.3b_freq, "Joint_Dose", table_freqmodel_dLight1.3b_joint_Dose)
write_nice_sheet(wb_dLight1.3b_freq, "Pairs_InjType", table_freqmodel_dLight1.3b_pairs_InjType)
write_nice_sheet(wb_dLight1.3b_freq, "Pairs_Dose",  table_freqmodel_dLight1.3b_pairs_Dose)

saveWorkbook(wb_dLight1.3b_freq, paste(analysispath,"dLight1.3b_Frequency_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_dLight1.3b_freq_ar1, "MeansRaw", treatmeanstable_freq_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_freq_ar1, "ModelSummary", table_freqmodel_dLight1.3b_ar1)
write_nice_sheet(wb_dLight1.3b_freq_ar1, "ModelAnova", table_freqmodel_dLight1.3b_ar1_aov)
write_nice_sheet(wb_dLight1.3b_freq_ar1, "ModelR2", freqmodel_dLight1.3b_ar1_r2)
write_nice_sheet(wb_dLight1.3b_freq_ar1, "ModelAIC", freqmodel_dLight1.3b_ar1_AIC)
write_nice_sheet(wb_dLight1.3b_freq_ar1, "Joint_InjType",  table_freqmodel_dLight1.3b_ar1_joint_InjType)
write_nice_sheet(wb_dLight1.3b_freq_ar1, "Joint_Dose", table_freqmodel_dLight1.3b_ar1_joint_Dose)
write_nice_sheet(wb_dLight1.3b_freq_ar1, "Pairs_InjType", table_freqmodel_dLight1.3b_ar1_pairs_InjType)
write_nice_sheet(wb_dLight1.3b_freq_ar1, "Pairs_Dose",  table_freqmodel_dLight1.3b_ar1_pairs_Dose)

saveWorkbook(wb_dLight1.3b_freq_ar1, paste(analysispath,"dLight1.3b_Frequency_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Change Models
write_nice_sheet(wb_dLight1.3b_freqC, "MeansChange", treatmeanstable_freq_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_freqC, "ModelSummary", table_freqCmodel_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_freqC, "ModelAnova", table_freqCmodel_dLight1.3b_aov)
write_nice_sheet(wb_dLight1.3b_freqC, "ModelR2", freqCmodel_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_freqC, "Pairs_Dose",  table_freqCmodel_dLight1.3b_pairs_Dose)

saveWorkbook(wb_dLight1.3b_freqC, paste(analysispath,"dLight1.3b_FrequencyChange_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Change Models
write_nice_sheet(wb_dLight1.3b_freqC_ar1, "MeansChange", treatmeanstable_freq_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_freqC_ar1, "ModelSummary", table_freqCmodel_dLight1.3b_ar1)
write_nice_sheet(wb_dLight1.3b_freqC_ar1, "ModelAnova", table_freqCmodel_dLight1.3b_ar1_aov)
write_nice_sheet(wb_dLight1.3b_freqC_ar1, "ModelR2", freqCmodel_dLight1.3b_ar1_r2)
write_nice_sheet(wb_dLight1.3b_freqC_ar1, "ModelAIC", freqCmodel_dLight1.3b_ar1_AIC)
write_nice_sheet(wb_dLight1.3b_freqC_ar1, "Pairs_Dose",  table_freqCmodel_dLight1.3b_ar1_pairs_Dose)

saveWorkbook(wb_dLight1.3b_freqC_ar1, paste(analysispath,"dLight1.3b_FrequencyChange_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Regression Models
write_nice_sheet(wb_dLight1.3b_freq_reg, "MaxHeightSummary", table_freqmodel_dLight1.3b_MaxHeight)
write_nice_sheet(wb_dLight1.3b_freq_reg, "MaxHeightAnova", table_freqmodel_dLight1.3b_MaxHeight_aov)
write_nice_sheet(wb_dLight1.3b_freq_reg, "MaxHeightR2", MaxHeightmodel_freq_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_freq_reg, "MaxHeightReg",  table_freqmodel_dLight1.3b_MaxHeight_reg)
write_nice_sheet(wb_dLight1.3b_freq_reg, "MaxDurationSummary", table_freqmodel_dLight1.3b_MaxDuration)
write_nice_sheet(wb_dLight1.3b_freq_reg, "MaxDurationAnova", table_freqmodel_dLight1.3b_MaxDuration_aov)
write_nice_sheet(wb_dLight1.3b_freq_reg, "MaxDurationR2", MaxDurationmodel_freq_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_freq_reg, "MaxDurationReg", table_freqmodel_dLight1.3b_MaxDuration_reg)
write_nice_sheet(wb_dLight1.3b_freq_reg, "SlopeSummary", table_freqmodel_dLight1.3b_MaxSlopeLM)
write_nice_sheet(wb_dLight1.3b_freq_reg, "SlopeAnova", table_freqmodel_dLight1.3b_MaxSlopeLM_aov)
write_nice_sheet(wb_dLight1.3b_freq_reg, "SlopeR2", MaxSlopeLMmodel_freq_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_freq_reg, "SlopeReg", table_freqmodel_dLight1.3b_MaxSlopeLM_reg)
write_nice_sheet(wb_dLight1.3b_freq_reg, "AUCSummary", table_freqmodel_dLight1.3b_AUC)
write_nice_sheet(wb_dLight1.3b_freq_reg, "AUCAnova", table_freqmodel_dLight1.3b_AUC_aov)
write_nice_sheet(wb_dLight1.3b_freq_reg, "AUCR2", AUCmodel_freq_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_freq_reg, "AUCReg", table_freqmodel_dLight1.3b_AUC_reg)

saveWorkbook(wb_dLight1.3b_freq_reg, paste(analysispath,"dLight1.3b_Frequency_regressionmodels_results.xlsx",sep=''), overwrite = TRUE)

### Amplitude Export -----
wb_dLight1.3b_amp <- createWorkbook()
wb_dLight1.3b_amp_ar1 <- createWorkbook()
wb_dLight1.3b_ampC <- createWorkbook()
wb_dLight1.3b_ampC_ar1 <- createWorkbook()
wb_dLight1.3b_amp_reg <- createWorkbook()

## Main results for dLight1.3b
# Standard Models
write_nice_sheet(wb_dLight1.3b_amp, "MeansRaw", treatmeansCtable_amp_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_amp, "ModelSummary", table_ampmodel_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_amp, "ModelAnova", table_ampmodel_dLight1.3b_aov)
write_nice_sheet(wb_dLight1.3b_amp, "ModelR2", ampmodel_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_amp, "Joint_InjType",  table_ampmodel_dLight1.3b_joint_InjType)
write_nice_sheet(wb_dLight1.3b_amp, "Joint_Dose", table_ampmodel_dLight1.3b_joint_Dose)
write_nice_sheet(wb_dLight1.3b_amp, "Pairs_InjType", table_ampmodel_dLight1.3b_pairs_InjType)
write_nice_sheet(wb_dLight1.3b_amp, "Pairs_Dose",  table_ampmodel_dLight1.3b_pairs_Dose)

saveWorkbook(wb_dLight1.3b_amp, paste(analysispath,"dLight1.3b_Amplitude_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_dLight1.3b_amp_ar1, "MeansRaw", treatmeanstable_amp_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_amp_ar1, "ModelSummary", table_ampmodel_dLight1.3b_ar1)
write_nice_sheet(wb_dLight1.3b_amp_ar1, "ModelAnova", table_ampmodel_dLight1.3b_ar1_aov)
write_nice_sheet(wb_dLight1.3b_amp_ar1, "ModelR2", ampmodel_dLight1.3b_ar1_r2)
write_nice_sheet(wb_dLight1.3b_amp_ar1, "ModelAIC", ampmodel_dLight1.3b_ar1_AIC)
write_nice_sheet(wb_dLight1.3b_amp_ar1, "Joint_InjType",  table_ampmodel_dLight1.3b_ar1_joint_InjType)
write_nice_sheet(wb_dLight1.3b_amp_ar1, "Joint_Dose", table_ampmodel_dLight1.3b_ar1_joint_Dose)
write_nice_sheet(wb_dLight1.3b_amp_ar1, "Pairs_InjType", table_ampmodel_dLight1.3b_ar1_pairs_InjType)
write_nice_sheet(wb_dLight1.3b_amp_ar1, "Pairs_Dose",  table_ampmodel_dLight1.3b_ar1_pairs_Dose)

saveWorkbook(wb_dLight1.3b_amp_ar1, paste(analysispath,"dLight1.3b_Amplitude_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Change Models
write_nice_sheet(wb_dLight1.3b_ampC, "MeansChange", treatmeanstable_amp_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_ampC, "ModelSummary", table_ampCmodel_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_ampC, "ModelAnova", table_ampCmodel_dLight1.3b_aov)
write_nice_sheet(wb_dLight1.3b_ampC, "ModelR2", ampCmodel_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_ampC, "Pairs_Dose",  table_ampCmodel_dLight1.3b_pairs_Dose)

saveWorkbook(wb_dLight1.3b_ampC, paste(analysispath,"dLight1.3b_AmplitudeChange_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Change Models
write_nice_sheet(wb_dLight1.3b_ampC_ar1, "MeansChange", treatmeanstable_amp_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_ampC_ar1, "ModelSummary", table_ampCmodel_dLight1.3b_ar1)
write_nice_sheet(wb_dLight1.3b_ampC_ar1, "ModelAnova", table_ampCmodel_dLight1.3b_ar1_aov)
write_nice_sheet(wb_dLight1.3b_ampC_ar1, "ModelR2", ampCmodel_dLight1.3b_ar1_r2)
write_nice_sheet(wb_dLight1.3b_ampC_ar1, "ModelAIC", ampCmodel_dLight1.3b_ar1_AIC)
write_nice_sheet(wb_dLight1.3b_ampC_ar1, "Pairs_Dose",  table_ampCmodel_dLight1.3b_ar1_pairs_Dose)

saveWorkbook(wb_dLight1.3b_ampC_ar1, paste(analysispath,"dLight1.3b_AmplitudeChange_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Regression Models
write_nice_sheet(wb_dLight1.3b_amp_reg, "MaxHeightSummary", table_ampmodel_dLight1.3b_MaxHeight)
write_nice_sheet(wb_dLight1.3b_amp_reg, "MaxHeightAnova", table_ampmodel_dLight1.3b_MaxHeight_aov)
write_nice_sheet(wb_dLight1.3b_amp_reg, "MaxHeightR2", MaxHeightmodel_amp_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_amp_reg, "MaxHeightReg",  table_ampmodel_dLight1.3b_MaxHeight_reg)
write_nice_sheet(wb_dLight1.3b_amp_reg, "MaxDurationSummary", table_ampmodel_dLight1.3b_MaxDuration)
write_nice_sheet(wb_dLight1.3b_amp_reg, "MaxDurationAnova", table_ampmodel_dLight1.3b_MaxDuration_aov)
write_nice_sheet(wb_dLight1.3b_amp_reg, "MaxDurationR2", MaxDurationmodel_amp_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_amp_reg, "MaxDurationReg", table_ampmodel_dLight1.3b_MaxDuration_reg)
write_nice_sheet(wb_dLight1.3b_amp_reg, "SlopeSummary", table_ampmodel_dLight1.3b_MaxSlopeLM)
write_nice_sheet(wb_dLight1.3b_amp_reg, "SlopeAnova", table_ampmodel_dLight1.3b_MaxSlopeLM_aov)
write_nice_sheet(wb_dLight1.3b_amp_reg, "SlopeR2", MaxSlopeLMmodel_amp_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_amp_reg, "SlopeReg", table_ampmodel_dLight1.3b_MaxSlopeLM_reg)
write_nice_sheet(wb_dLight1.3b_amp_reg, "AUCSummary", table_ampmodel_dLight1.3b_AUC)
write_nice_sheet(wb_dLight1.3b_amp_reg, "AUCAnova", table_ampmodel_dLight1.3b_AUC_aov)
write_nice_sheet(wb_dLight1.3b_amp_reg, "AUCR2", AUCmodel_amp_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_amp_reg, "AUCReg", table_ampmodel_dLight1.3b_AUC_reg)

saveWorkbook(wb_dLight1.3b_amp_reg, paste(analysispath,"dLight1.3b_Amplitude_regressionmodels_results.xlsx",sep=''), overwrite = TRUE)


### AUCwindow Export -----
wb_dLight1.3b_AUCwindow <- createWorkbook()
wb_dLight1.3b_AUCwindow_ar1 <- createWorkbook()
wb_dLight1.3b_AUCwindowC <- createWorkbook()
wb_dLight1.3b_AUCwindowC_ar1 <- createWorkbook()
wb_dLight1.3b_AUCwindow_reg <- createWorkbook()

## Main results for dLight1.3b
# Standard Models
write_nice_sheet(wb_dLight1.3b_AUCwindow, "MeansRaw", treatmeanstable_AUCwindow_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_AUCwindow, "ModelSummary", table_AUCwindowmodel_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_AUCwindow, "ModelAnova", table_AUCwindowmodel_dLight1.3b_aov)
write_nice_sheet(wb_dLight1.3b_AUCwindow, "ModelR2", AUCwindowmodel_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_AUCwindow, "Joint_InjType",  table_AUCwindowmodel_dLight1.3b_joint_InjType)
write_nice_sheet(wb_dLight1.3b_AUCwindow, "Joint_Dose", table_AUCwindowmodel_dLight1.3b_joint_Dose)
write_nice_sheet(wb_dLight1.3b_AUCwindow, "Pairs_InjType", table_AUCwindowmodel_dLight1.3b_pairs_InjType)
write_nice_sheet(wb_dLight1.3b_AUCwindow, "Pairs_Dose",  table_AUCwindowmodel_dLight1.3b_pairs_Dose)

saveWorkbook(wb_dLight1.3b_AUCwindow, paste(analysispath,"dLight1.3b_AUCwindow_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_dLight1.3b_AUCwindow_ar1, "MeansRaw", treatmeanstable_AUCwindow_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_AUCwindow_ar1, "ModelSummary", table_AUCwindowmodel_dLight1.3b_ar1)
write_nice_sheet(wb_dLight1.3b_AUCwindow_ar1, "ModelAnova", table_AUCwindowmodel_dLight1.3b_ar1_aov)
write_nice_sheet(wb_dLight1.3b_AUCwindow_ar1, "ModelR2", AUCwindowmodel_dLight1.3b_ar1_r2)
write_nice_sheet(wb_dLight1.3b_AUCwindow_ar1, "ModelAIC", AUCwindowmodel_dLight1.3b_ar1_AIC)
write_nice_sheet(wb_dLight1.3b_AUCwindow_ar1, "Joint_InjType",  table_AUCwindowmodel_dLight1.3b_ar1_joint_InjType)
write_nice_sheet(wb_dLight1.3b_AUCwindow_ar1, "Joint_Dose", table_AUCwindowmodel_dLight1.3b_ar1_joint_Dose)
write_nice_sheet(wb_dLight1.3b_AUCwindow_ar1, "Pairs_InjType", table_AUCwindowmodel_dLight1.3b_ar1_pairs_InjType)
write_nice_sheet(wb_dLight1.3b_AUCwindow_ar1, "Pairs_Dose",  table_AUCwindowmodel_dLight1.3b_ar1_pairs_Dose)

saveWorkbook(wb_dLight1.3b_AUCwindow_ar1, paste(analysispath,"dLight1.3b_AUCwindow_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Change Models
write_nice_sheet(wb_dLight1.3b_AUCwindowC, "MeansChange", treatmeanstable_AUCwindow_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_AUCwindowC, "ModelSummary", table_AUCwindowCmodel_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_AUCwindowC, "ModelAnova", table_AUCwindowCmodel_dLight1.3b_aov)
write_nice_sheet(wb_dLight1.3b_AUCwindowC, "ModelR2", AUCwindowCmodel_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_AUCwindowC, "Pairs_Dose",  table_AUCwindowCmodel_dLight1.3b_pairs_Dose)

saveWorkbook(wb_dLight1.3b_AUCwindowC, paste(analysispath,"dLight1.3b_AUCwindowChange_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Change Models
write_nice_sheet(wb_dLight1.3b_AUCwindowC_ar1, "MeansChange", treatmeanstable_AUCwindow_dLight1.3b)
write_nice_sheet(wb_dLight1.3b_AUCwindowC_ar1, "ModelSummary", table_AUCwindowCmodel_dLight1.3b_ar1)
write_nice_sheet(wb_dLight1.3b_AUCwindowC_ar1, "ModelAnova", table_AUCwindowCmodel_dLight1.3b_ar1_aov)
write_nice_sheet(wb_dLight1.3b_AUCwindowC_ar1, "ModelR2", AUCwindowCmodel_dLight1.3b_ar1_r2)
write_nice_sheet(wb_dLight1.3b_AUCwindowC_ar1, "ModelAIC", AUCwindowCmodel_dLight1.3b_ar1_AIC)
write_nice_sheet(wb_dLight1.3b_AUCwindowC_ar1, "Pairs_Dose",  table_AUCwindowCmodel_dLight1.3b_ar1_pairs_Dose)

saveWorkbook(wb_dLight1.3b_AUCwindowC_ar1, paste(analysispath,"dLight1.3b_AUCwindowChange_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Regression Models
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "MaxHeightSummary", table_AUCwindowmodel_dLight1.3b_MaxHeight)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "MaxHeightAnova", table_AUCwindowmodel_dLight1.3b_MaxHeight_aov)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "MaxHeightR2", MaxHeightmodel_AUCwindow_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "MaxHeightReg",  table_AUCwindowmodel_dLight1.3b_MaxHeight_reg)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "MaxDurationSummary", table_AUCwindowmodel_dLight1.3b_MaxDuration)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "MaxDurationAnova", table_AUCwindowmodel_dLight1.3b_MaxDuration_aov)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "MaxDurationR2", MaxDurationmodel_AUCwindow_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "MaxDurationReg", table_AUCwindowmodel_dLight1.3b_MaxDuration_reg)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "SlopeSummary", table_AUCwindowmodel_dLight1.3b_MaxSlopeLM)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "SlopeAnova", table_AUCwindowmodel_dLight1.3b_MaxSlopeLM_aov)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "SlopeR2", MaxSlopeLMmodel_AUCwindow_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "SlopeReg", table_AUCwindowmodel_dLight1.3b_MaxSlopeLM_reg)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "AUCSummary", table_AUCwindowmodel_dLight1.3b_AUC)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "AUCAnova", table_AUCwindowmodel_dLight1.3b_AUC_aov)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "AUCR2", AUCmodel_AUCwindow_dLight1.3b_r2)
write_nice_sheet(wb_dLight1.3b_AUCwindow_reg, "AUCReg", table_AUCwindowmodel_dLight1.3b_AUC_reg)

saveWorkbook(wb_dLight1.3b_AUCwindow_reg, paste(analysispath,"dLight1.3b_AUCwindow_regressionmodels_results.xlsx",sep=''), overwrite = TRUE)


## PLOTS -----
### Frequency -----
# Bar Graph (Post-Injection Mean)
freqbarmin_dLight1.3b <- 0
freqbarmax_dLight1.3b <- 9.2
freqbarticks_dLight1.3b <- 2
freqbarbreaks_dLight1.3b <- seq(freqbarmin_dLight1.3b,freqbarmax_dLight1.3b,freqbarticks_dLight1.3b)

freqbar_dLight1.3b <- ggplot(subset(treatmeans, Virus == "dLight1.3b"&CondNum==1),aes(x=Dose, y=pkspermin_mean, fill=TreatNum, group=InjType)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_jitter(data=subset(subjectmeans,Virus == "dLight1.3b" & CondNum == 1), aes(x=Dose, y=pkspermin, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "dLight1.3b"&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=pkspermin_mean-pkspermin_se, ymax=pkspermin_mean+pkspermin_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(freqbarmin_dLight1.3b, freqbarmax_dLight1.3b), 
                     breaks=freqbarbreaks_dLight1.3b) +
  xlab("Dose (mg/kg)") +
  ylab("Frequency (n/min)") + 
  mytheme +
  # STATS FOR AR1 MIXED MODEL
  geom_bracket(xmin=c(0.765), xmax = c(1.235), y.position=c(8),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(3.765), xmax = c(4.235), y.position=c(8),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

freqbar_dLight1.3b

# Bar graph (Change from Control)
freqCbarmin_dLight1.3b <- -4
freqCbarmax_dLight1.3b <- 8
freqCbarticks_dLight1.3b <- 4
freqCbarbreaks_dLight1.3b <- seq(freqCbarmin_dLight1.3b,freqCbarmax_dLight1.3b,freqCbarticks_dLight1.3b)

freqCbar_dLight1.3b <- ggplot(subset(treatmeans, Virus == "dLight1.3b"&CondNum==1&InjType=='M'),aes(x=Dose, y=pksperminC_mean, fill=TreatNum, group=InjType)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(subjectmeans,Virus == "dLight1.3b" & CondNum == 1&InjType=='M'), aes(x=Dose, y=pksperminC, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "dLight1.3b"&CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=pksperminC_mean-pksperminC_se, ymax=pksperminC_mean+pksperminC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(freqCbarmin_dLight1.3b, freqCbarmax_dLight1.3b), breaks=freqCbarbreaks_dLight1.3b)+
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" Frequency (n/min)")) + 
  mytheme +
  # AR1 and ALL AMP STATS
  geom_bracket(xmin=c(1), xmax = c(2.9), y.position=c(7.4),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2), xmax = c(4), y.position=c(5.9),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.1), xmax = c(4), y.position=c(7.4),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

freqCbar_dLight1.3b

# Line Graph (Change from Control)
freqlineCmin_dLight1.3b <- -1.8
freqlineCmax_dLight1.3b <- 3.6
freqlineCticks_dLight1.3b <- 1.5
freqlineCbreaks_dLight1.3b <- seq(-1.5, freqlineCmax_dLight1.3b, freqlineCticks_dLight1.3b)

freqline_dLight1.3b <- ggplot(subset(treatmeansbin, Virus=='dLight1.3b'&InjType!='S'), aes(x=BinTimeInj, y=pksperminC_mean, group=Dose, fill=Dose, color=Dose)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_continuous(name="Minutes from Injection",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(freqlineCmin_dLight1.3b, freqlineCmax_dLight1.3b), 
                     breaks = freqlineCbreaks_dLight1.3b)+
  scale_color_manual("legend", values = treatcolors) +
  mytheme +
  ylab(expression(Delta*" Frequency (n/min)")) + 
  mytheme

freqline_dLight1.3b

#### Temporal Dynamics Panels -----
# Max Height
MaxHeight_freq_ymin_dLight1.3b <- -.35
MaxHeight_freq_ymax_dLight1.3b <- 13.5
MaxHeight_freq_ybreaks_dLight1.3b <- 4

MaxHeight_freq_dLight1.3b <- ggplot(subset(subjectrise_pksperminC, Virus=='dLight1.3b'&InjType=='M'), aes(x=DoseNum, y=MaxHeight, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxHeight_freq_dLight1.3b_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxHeight_freq_dLight1.3b_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxHeight_freq_ymin_dLight1.3b, MaxHeight_freq_ymax_dLight1.3b), 
                     breaks = seq(0, MaxHeight_freq_ymax_dLight1.3b, MaxHeight_freq_ybreaks_dLight1.3b)) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Max "*Delta*" n"))

MaxHeight_freq_dLight1.3b

# Max Duration
MaxDuration_freq_ymin_dLight1.3b <- 0
MaxDuration_freq_ymax_dLight1.3b <- 63
MaxDuration_freq_ybreaks_dLight1.3b <- 20

MaxDuration_freq_dLight1.3b <- ggplot(subset(subjectrise_pksperminC, Virus=='dLight1.3b'&InjType=='M'), aes(x=DoseNum, y=MaxDuration, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxDuration_freq_dLight1.3b_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxDuration_freq_dLight1.3b_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxDuration_freq_ymin_dLight1.3b, MaxDuration_freq_ymax_dLight1.3b), 
                     breaks = seq(MaxDuration_freq_ymin_dLight1.3b, MaxDuration_freq_ymax_dLight1.3b, MaxDuration_freq_ybreaks_dLight1.3b)) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Time to Max "*Delta*" n"))

MaxDuration_freq_dLight1.3b


# Max Height Slope
MaxSlopeLM_freq_ymin_dLight1.3b <- -.1
MaxSlopeLM_freq_ymax_dLight1.3b <- .55
MaxSlopeLM_freq_ybreaks_dLight1.3b <- .15

MaxSlopeLM_freq_dLight1.3b <- ggplot(subset(subjectrise_pksperminC, Virus=='dLight1.3b'&InjType=='M'), aes(x=DoseNum, y=MaxSlopeLM, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxSlopeLM_freq_dLight1.3b_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxSlopeLM_freq_dLight1.3b_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxSlopeLM_freq_ymin_dLight1.3b, MaxSlopeLM_freq_ymax_dLight1.3b), 
                     breaks = seq(0, MaxSlopeLM_freq_ymax_dLight1.3b, MaxSlopeLM_freq_ybreaks_dLight1.3b))+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab('Time-Course Slope')

MaxSlopeLM_freq_dLight1.3b

# AUC POST INJ
AUC_freq_ymin <- -5
AUC_freq_ymax <- 405
AUC_freq_ybreaks <- 125

AUC_freq_dLight1.3b <- ggplot(subset(subjectrise_pksperminC, Virus=='dLight1.3b'&InjType=='M'), aes(x=DoseNum, y=AUC_POSTINJ, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = AUC_freq_dLight1.3b_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = AUC_freq_dLight1.3b_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(AUC_freq_ymin, AUC_freq_ymax), breaks = seq(0, AUC_freq_ymax, AUC_freq_ybreaks))+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab("Time-Course AUC")

AUC_freq_dLight1.3b

### Amplitude -----
# Bar graph (Post-Injection Mean)
ampbarmin_dLight1.3b <- 0
ampbarmax_dLight1.3b <- 2
ampbarticks_dLight1.3b <- .5
ampbarybreaks_dLight1.3b <- seq(ampbarmin_dLight1.3b,ampbarmax_dLight1.3b,ampbarticks_dLight1.3b)
ampbarylabels_dLight1.3b <- ampbarybreaks_dLight1.3b + 2.6

ampbar_dLight1.3b <- ggplot(subset(treatmeans, Virus == "dLight1.3b"&CondNum==1),aes(x=Dose, y=ampscaled_mean, fill=TreatNum, group=InjType)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_jitter(data=subset(subjectmeans,Virus == "dLight1.3b" & CondNum == 1), aes(x=Dose, y=ampscaled, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "dLight1.3b"&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=ampscaled_mean-amp_se, ymax=ampscaled_mean+amp_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(ampbarmin_dLight1.3b, ampbarmax_dLight1.3b), breaks=ampbarybreaks_dLight1.3b,
                     labels=ampbarylabels_dLight1.3b)+
  xlab("Dose (mg/kg)") +
  ylab("Amplitude (z)") + 
  mytheme +
  # AR1 and ALL AMP STATS
  geom_bracket(xmin=c(0.765), xmax = c(1.235), y.position=c(1.92),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1.765), xmax = c(2.235), y.position=c(1.92),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.765), xmax = c(3.235), y.position=c(1.88),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.765), xmax = c(4.235), y.position=c(1.92),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

ampbar_dLight1.3b

# Bar graph (Change from Control)
ampCbarmin_dLight1.3b <- -.25
ampCbarmax_dLight1.3b <- .6
ampCbarticks_dLight1.3b <- .2
ampCbarbreaks_dLight1.3b <- seq(-.2,ampCbarmax_dLight1.3b,ampCbarticks_dLight1.3b)

ampCbar_dLight1.3b <- ggplot(subset(treatmeans, Virus == "dLight1.3b"&CondNum==1&InjType=='M'),aes(x=Dose, y=ampC_mean, fill=TreatNum, group=InjType)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(subjectmeans,Virus == "dLight1.3b" & CondNum == 1&InjType=='M'), aes(x=Dose, y=ampC, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "dLight1.3b"&CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=ampC_mean-ampC_se, ymax=ampC_mean+ampC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(ampCbarmin_dLight1.3b, ampCbarmax_dLight1.3b), breaks=ampCbarbreaks_dLight1.3b)+
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" Amplitude (z)")) + 
  mytheme

ampCbar_dLight1.3b

# Line Graph (Change from Control)
amplineCmin <- -.6
amplineCmax <- .65
amplineCticks <- .3

ampline_dLight1.3b <- ggplot(subset(treatmeansbin, Virus=='dLight1.3b'&InjType!='S'), aes(x=BinTimeInj, y=ampC_mean, group=Dose, fill=Dose, color=Dose)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  scale_x_continuous(name="Minutes from Injection",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(amplineCmin, amplineCmax), breaks = seq(amplineCmin, amplineCmax, amplineCticks))+
  scale_color_manual("legend", values = treatcolors) +
  scale_fill_manual("legend", values = treatcolors) +
  mytheme +
  ylab(expression(Delta*" Amplitude (z)")) + 
  mytheme

ampline_dLight1.3b

#### Temporal Dynamics Panels -----
# Max Height
MaxHeight_amp_ymin_dLight1.3b <- 0
MaxHeight_amp_ymax_dLight1.3b <- 3
MaxHeight_amp_ybreaks_dLight1.3b <- 1

MaxHeight_amp_dLight1.3b <- ggplot(subset(subjectrise_ampC, Virus=='dLight1.3b'&InjType=='M'), aes(x=DoseNum, y=MaxHeight, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxHeight_amp_dLight1.3b_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxHeight_amp_dLight1.3b_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxHeight_amp_ymin_dLight1.3b, MaxHeight_amp_ymax_dLight1.3b), 
                     breaks = seq(MaxHeight_amp_ymin_dLight1.3b, MaxHeight_amp_ymax_dLight1.3b, MaxHeight_amp_ybreaks_dLight1.3b)) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Max "*Delta*" z"))

MaxHeight_amp_dLight1.3b

# Max Duration
MaxDuration_amp_ymin_dLight1.3b <- 0
MaxDuration_amp_ymax_dLight1.3b <- 60
MaxDuration_amp_ybreaks_dLight1.3b <- 20

MaxDuration_amp_dLight1.3b <- ggplot(subset(subjectrise_ampC, Virus=='dLight1.3b'&InjType=='M'), aes(x=DoseNum, y=MaxDuration, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxDuration_amp_dLight1.3b_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxDuration_amp_dLight1.3b_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxDuration_amp_ymin_dLight1.3b, MaxDuration_amp_ymax_dLight1.3b), 
                     breaks = seq(MaxDuration_amp_ymin_dLight1.3b, MaxDuration_amp_ymax_dLight1.3b, MaxDuration_amp_ybreaks_dLight1.3b)) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Time to Max "*Delta*" z"))

MaxDuration_amp_dLight1.3b

# Max Height Slope
MaxSlopeLM_amp_ymin_dLight1.3b <- -.03
MaxSlopeLM_amp_ymax_dLight1.3b <- .38
MaxSlopeLM_amp_ybreaks_dLight1.3b <- .12

MaxSlopeLM_amp_dLight1.3b <- ggplot(subset(subjectrise_ampC, Virus=='dLight1.3b'&InjType=='M'), aes(x=DoseNum, y=MaxSlopeLM, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxSlopeLM_amp_dLight1.3b_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxSlopeLM_amp_dLight1.3b_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxSlopeLM_amp_ymin_dLight1.3b, MaxSlopeLM_amp_ymax_dLight1.3b), 
                     breaks = seq(0, MaxSlopeLM_amp_ymax_dLight1.3b, MaxSlopeLM_amp_ybreaks_dLight1.3b))+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab('Time-Course Slope')

MaxSlopeLM_amp_dLight1.3b

# AUC POST INJ
AUC_amp_ymin_dLight1.3b <- 0
AUC_amp_ymax_dLight1.3b <- 56
AUC_amp_ybreaks_dLight1.3b <- 18

AUC_amp_dLight1.3b <- ggplot(subset(subjectrise_ampC, Virus=='dLight1.3b'&InjType=='M'), aes(x=DoseNum, y=AUC_POSTINJ, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = AUC_amp_dLight1.3b_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = AUC_amp_dLight1.3b_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(AUC_amp_ymin_dLight1.3b, AUC_amp_ymax_dLight1.3b), 
                     breaks = seq(0, AUC_amp_ymax_dLight1.3b, AUC_amp_ybreaks_dLight1.3b))+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab("Time-Course AUC")

AUC_amp_dLight1.3b

### AUC Window -----
# Bar graph (Post-Injection Mean)
AUCwindowbarmin <- 0
AUCwindowbarmax <- 3.2
AUCwindowbarticks <- .8
AUCwindowbarbreaks <- seq(AUCwindowbarmin,AUCwindowbarmax,AUCwindowbarticks)

# Bar graph
AUCwindowbar_dLight1.3b <- ggplot(subset(treatmeans, Virus == "dLight1.3b"&CondNum==1),aes(x=Dose, y=AUCwindow_mean, fill=TreatNum, group=InjType)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_jitter(data=subset(subjectmeans,Virus == "dLight1.3b" & CondNum == 1), aes(x=Dose, y=AUCwindow, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "dLight1.3b"&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=AUCwindow_mean-AUCwindow_se, ymax=AUCwindow_mean+AUCwindow_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(AUCwindowbarmin, AUCwindowbarmax), breaks=AUCwindowbarbreaks)+
  xlab("Dose (mg/kg)") +
  ylab("AUC per Transient") + 
  mytheme

AUCwindowbar_dLight1.3b

# Bar graph (Change from Control)
AUCwindowCbarmin_dLight1.3b <- -.4
AUCwindowCbarmax_dLight1.3b <- .8
AUCwindowCbarticks_dLight1.3b <- .4
AUCwindowCbarbreaks_dLight1.3b <- seq(AUCwindowCbarmin_dLight1.3b,AUCwindowCbarmax_dLight1.3b,AUCwindowCbarticks_dLight1.3b)

AUCwindowCbar_dLight1.3b <- ggplot(subset(treatmeans, Virus == "dLight1.3b"&CondNum==1&InjType=='M'),aes(x=Dose, y=AUCwindowC_mean, fill=TreatNum, group=InjType)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(subjectmeans,Virus == "dLight1.3b" & CondNum == 1&InjType=='M'), aes(x=Dose, y=AUCwindowC, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "dLight1.3b"&CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=AUCwindowC_mean-AUCwindowC_se, ymax=AUCwindowC_mean+AUCwindowC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(AUCwindowCbarmin_dLight1.3b, AUCwindowCbarmax_dLight1.3b), breaks=AUCwindowCbarbreaks_dLight1.3b)+
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" AUC per Transient")) + 
  mytheme

AUCwindowCbar_dLight1.3b

# Line Graph (Change from Control)
AUCwindowlineCmin <- -.6
AUCwindowlineCmax <- .6
AUCwindowlineCticks <- .3
AUCwindowlineCbreaks <- seq(AUCwindowlineCmin, AUCwindowlineCmax, AUCwindowlineCticks)

AUCwindowline_dLight1.3b <- ggplot(subset(treatmeansbin, Virus=='dLight1.3b'&InjType!='S'), aes(x=BinTimeInj, y=AUCwindowC_mean, group=Dose, fill=Dose, color=Dose)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  scale_x_continuous(name="Minutes from Injection",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(AUCwindowlineCmin, AUCwindowlineCmax), breaks = AUCwindowlineCbreaks)+
  scale_color_manual("legend", values = treatcolors) +
  scale_fill_manual("legend", values = treatcolors) +
  mytheme +
  ylab(expression(Delta*" AUC per Transient")) + 
  mytheme

AUCwindowline_dLight1.3b

### EXPORT FIGURE 3 -----
# BAR AND LINE GRAPHS
figure3panels_barline <- plot_grid(freqbar_dLight1.3b, freqCbar_dLight1.3b, freqline_dLight1.3b, 
                                   ampbar_dLight1.3b, ampCbar_dLight1.3b, ampline_dLight1.3b,
                                   AUCwindowbar_dLight1.3b, AUCwindowCbar_dLight1.3b, AUCwindowline_dLight1.3b,
                                   ncol = 3, nrow=3, align = "hv",  # align horizontally & vertically
                                   axis  = "tblr", hjust=-.5, rel_widths = c(1,.8,1.2,
                                                                             1,.8,1.2,
                                                                             1,.8,1.2))
figure3panels_barline


ggsave(filename = c("Figure3_NAcLSdLight1.3b_Bar_TransientPanels_BarLineGraphs.pdf"), path = figurepath, plot = figure3panels_barline, device = "pdf", 
       width = figpanelwidthin, height = figpanelheightinbarline, dpi=figdpi, units="in")


# TEMPORAL DYNAMICS GRAPHS
figure3panels_models <- plot_grid(MaxHeight_freq_dLight1.3b, MaxDuration_freq_dLight1.3b,
                                  MaxSlopeLM_freq_dLight1.3b, AUC_freq_dLight1.3b,
                                  
                                  MaxHeight_amp_dLight1.3b, MaxDuration_amp_dLight1.3b,
                                  MaxSlopeLM_amp_dLight1.3b, AUC_amp_dLight1.3b,
                                  
                                  ncol = 4, nrow=2, align = "hv",  # align horizontally & vertically
                                  axis  = "tblr", hjust=-.5)
figure3panels_models

ggsave(filename = c("Figure3_NAcLSdLight1.3b_Bar_TransientPanels_ModelGraphs.pdf"), path = figurepath, plot = figure3panels_models, device = "pdf", 
       width = figpanelwidthin, height = figpanelheightinmodel, dpi=figdpi, units="in")


# FIGURE 4 - NAcLS GRABDA2h -----
## STATS ------
### Frequency -----
#### Frequency Model -----
#### Mixed effects model
freqmodel_GRABDA2h <- glmmTMB(pkspermin ~ InjType*Dose*BinF + (1|SubjectID), # random intercept by subject
                                data = subset(subjectmeansbin, Virus=='GRABDA2h' & CondNum == 1), family = stats::gaussian,
                                control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqmodel_GRABDA2h_aov <- Anova(freqmodel_GRABDA2h, type=2)
freqmodel_GRABDA2h_aov
summary(freqmodel_GRABDA2h)

freqmodel_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(freqmodel_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqmodel_GRABDA2h_r2

# Follow up tests - joint
freqmodel_GRABDA2h_emmeans <-emmeans(freqmodel_GRABDA2h, pairwise~InjType|Dose|BinF, adjust="none") # Find means
freqmodel_GRABDA2h_jt_Dose <- joint_tests(freqmodel_GRABDA2h_emmeans, by = "Dose")
freqmodel_GRABDA2h_jt_InjType <- joint_tests(freqmodel_GRABDA2h_emmeans, by = "InjType")

freqmodel_GRABDA2h_jt_Dose
freqmodel_GRABDA2h_jt_InjType

# Follow up tests - treat pairs
freqmodel_GRABDA2h_treatpairs <- emmeans(freqmodel_GRABDA2h, pairwise~InjType|Dose, adjust="none") # Pairs
freqmodel_GRABDA2h_treatcontrasts <- as.data.frame(freqmodel_GRABDA2h_treatpairs$contrasts) # Extract emmeans
freqmodel_GRABDA2h_treatemmeans <- freqmodel_GRABDA2h_treatpairs$emmeans # Extract emmeans
freqmodel_GRABDA2h_treatpairs_es <- as.data.frame(eff_size(freqmodel_GRABDA2h_treatemmeans, sigma = sigma(freqmodel_GRABDA2h), edf = df.residual(freqmodel_GRABDA2h))) # Effect Size
freqmodel_GRABDA2h_treatpairs_es <- freqmodel_GRABDA2h_treatpairs_es[, c("Dose", "contrast", "effect.size")]

freqmodel_GRABDA2h_treatcontrasts$d <- freqmodel_GRABDA2h_treatpairs_es$effect.size
freqmodel_GRABDA2h_treatcontrasts

# Follow up tests - dose pairs
freqmodel_GRABDA2h_dosepairs <- emmeans(freqmodel_GRABDA2h, pairwise~Dose|InjType, adjust="none") # Pairs
freqmodel_GRABDA2h_dosecontrasts <- as.data.frame(freqmodel_GRABDA2h_dosepairs$contrasts) # Extract emmeans
freqmodel_GRABDA2h_doseemmeans <- freqmodel_GRABDA2h_dosepairs$emmeans # Extract emmeans
freqmodel_GRABDA2h_dosepairs_es <- as.data.frame(eff_size(freqmodel_GRABDA2h_doseemmeans, sigma = sigma(freqmodel_GRABDA2h), edf = df.residual(freqmodel_GRABDA2h))) # Effect Size
freqmodel_GRABDA2h_dosepairs_es <- freqmodel_GRABDA2h_dosepairs_es[, c("InjType", "contrast", "effect.size")]

freqmodel_GRABDA2h_dosecontrasts$d <- freqmodel_GRABDA2h_dosepairs_es$effect.size
freqmodel_GRABDA2h_dosecontrasts

# Format Tables
table_freqmodel_GRABDA2h <- format_regsummary_table(freqmodel_GRABDA2h, 'GRABDA2h', "Frequency", "StandardModel")
table_freqmodel_GRABDA2h_aov <- format_anova_table(freqmodel_GRABDA2h_aov, "GRABDA2h", "Frequency", "StandardModel")
table_freqmodel_GRABDA2h_joint_Dose <- format_emm_table(freqmodel_GRABDA2h_jt_Dose, "GRABDA2h", "Frequency", "Joint", "Dose")
table_freqmodel_GRABDA2h_joint_InjType <- format_emm_table(freqmodel_GRABDA2h_jt_InjType, "GRABDA2h", "Frequency", "Joint", "InjType")
table_freqmodel_GRABDA2h_pairs_Dose <- format_emm_table(freqmodel_GRABDA2h_dosecontrasts, "GRABDA2h", "Frequency", "Pairs", "Dose|InjType")
table_freqmodel_GRABDA2h_pairs_InjType <- format_emm_table(freqmodel_GRABDA2h_treatcontrasts, "GRABDA2h", "Frequency", "Pairs", "InjType|Dose")

#### Mixed effects AR1 model
freqmodel_GRABDA2h_ar1 <- glmmTMB(pkspermin ~ InjType*Dose*BinF + ar1(BinF + 0|SubjectID:InjType) + (1|SubjectID), # random intercept by subject
                                    data = subset(subjectmeansbin, Virus=='GRABDA2h' & CondNum == 1), family = stats::gaussian, 
                                    control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqmodel_GRABDA2h_ar1_aov <- Anova(freqmodel_GRABDA2h_ar1, type=2)
freqmodel_GRABDA2h_ar1_aov
summary(freqmodel_GRABDA2h_ar1)

freqmodel_GRABDA2h_ar1_r2 <- as.data.frame(r.squaredGLMM(freqmodel_GRABDA2h_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqmodel_GRABDA2h_ar1_r2

freqmodel_GRABDA2h_ar1_AIC <- as.data.frame(AIC(freqmodel_GRABDA2h, freqmodel_GRABDA2h_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "freqmodel_GRABDA2h" ~ "StandardModel", ModelName == "freqmodel_GRABDA2h_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

freqmodel_GRABDA2h_ar1_AIC

# Follow up tests - joint
freqmodel_GRABDA2h_ar1_emmeans <-emmeans(freqmodel_GRABDA2h_ar1, pairwise~InjType|Dose, adjust="none") # Find means
freqmodel_GRABDA2h_ar1_jt_Dose <- joint_tests(freqmodel_GRABDA2h_ar1, by = "Dose")
freqmodel_GRABDA2h_ar1_jt_InjType <- joint_tests(freqmodel_GRABDA2h_ar1, by = "InjType")

freqmodel_GRABDA2h_ar1_jt_Dose
freqmodel_GRABDA2h_ar1_jt_InjType

# Follow up tests - treat pairs
freqmodel_GRABDA2h_ar1_treatpairs <- emmeans(freqmodel_GRABDA2h_ar1, pairwise~InjType|Dose, adjust="none") # Pairs
freqmodel_GRABDA2h_ar1_treatcontrasts <- as.data.frame(freqmodel_GRABDA2h_ar1_treatpairs$contrasts) # Extract emmeans
freqmodel_GRABDA2h_ar1_treatemmeans <- freqmodel_GRABDA2h_ar1_treatpairs$emmeans # Extract emmeans
freqmodel_GRABDA2h_ar1_treatpairs_es <- as.data.frame(eff_size(freqmodel_GRABDA2h_ar1_treatemmeans, sigma = sigma(freqmodel_GRABDA2h_ar1), edf = df.residual(freqmodel_GRABDA2h_ar1))) # Effect Size
freqmodel_GRABDA2h_ar1_treatpairs_es <- freqmodel_GRABDA2h_ar1_treatpairs_es[, c("Dose", "contrast", "effect.size")]

freqmodel_GRABDA2h_ar1_treatcontrasts$d <- freqmodel_GRABDA2h_ar1_treatpairs_es$effect.size
freqmodel_GRABDA2h_ar1_treatcontrasts

# Follow up tests - dose pairs
freqmodel_GRABDA2h_ar1_dosepairs <- emmeans(freqmodel_GRABDA2h_ar1, pairwise~Dose|InjType, adjust="none") # Pairs
freqmodel_GRABDA2h_ar1_dosecontrasts <- as.data.frame(freqmodel_GRABDA2h_ar1_dosepairs$contrasts) # Extract emmeans
freqmodel_GRABDA2h_ar1_doseemmeans <- freqmodel_GRABDA2h_ar1_dosepairs$emmeans # Extract emmeans
freqmodel_GRABDA2h_ar1_dosepairs_es <- as.data.frame(eff_size(freqmodel_GRABDA2h_ar1_doseemmeans, sigma = sigma(freqmodel_GRABDA2h_ar1), edf = df.residual(freqmodel_GRABDA2h_ar1))) # Effect Size
freqmodel_GRABDA2h_ar1_dosepairs_es <- freqmodel_GRABDA2h_ar1_dosepairs_es[, c("InjType", "contrast", "effect.size")]

freqmodel_GRABDA2h_ar1_dosecontrasts$d <- freqmodel_GRABDA2h_ar1_dosepairs_es$effect.size
freqmodel_GRABDA2h_ar1_dosecontrasts

# Format Tables
table_freqmodel_GRABDA2h_ar1 <- format_regsummary_table(freqmodel_GRABDA2h_ar1, 'GRABDA2h', "Frequency", "AR1Model")
table_freqmodel_GRABDA2h_ar1_aov <- format_anova_table(freqmodel_GRABDA2h_ar1_aov, "GRABDA2h", "Frequency", "AR1Model") 
table_freqmodel_GRABDA2h_ar1_joint_Dose <- format_emm_table(freqmodel_GRABDA2h_ar1_jt_Dose, "GRABDA2h", "Frequency", "Joint", "Dose")
table_freqmodel_GRABDA2h_ar1_joint_InjType <- format_emm_table(freqmodel_GRABDA2h_ar1_jt_InjType, "GRABDA2h", "Frequency", "Joint", "InjType")
table_freqmodel_GRABDA2h_ar1_pairs_Dose <- format_emm_table(freqmodel_GRABDA2h_ar1_dosecontrasts, "GRABDA2h", "Frequency", "Pairs", "Dose|InjType")
table_freqmodel_GRABDA2h_ar1_pairs_InjType <- format_emm_table(freqmodel_GRABDA2h_ar1_treatcontrasts, "GRABDA2h", "Frequency", "Pairs", "InjType|Dose")

#### Frequency Change Model -----
# Mixed effects model
freqCmodel_GRABDA2h <- glmmTMB(pksperminC ~ Dose*BinF + (1|SubjectID), # random intercept by subject
                                 data = subset(subjectmeansbin, Virus=='GRABDA2h' & CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                 control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqCmodel_GRABDA2h_aov <- Anova(freqCmodel_GRABDA2h, type=2)
freqCmodel_GRABDA2h_aov
summary(freqCmodel_GRABDA2h)

freqCmodel_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(freqCmodel_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqCmodel_GRABDA2h_r2

# Follow up tests - joint
freqCmodel_GRABDA2h_emmeans <-emmeans(freqCmodel_GRABDA2h, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
freqCmodel_GRABDA2h_dosepairs <- emmeans(freqCmodel_GRABDA2h, pairwise~Dose, adjust="none") # Pairs
freqCmodel_GRABDA2h_dosecontrasts <- as.data.frame(freqCmodel_GRABDA2h_dosepairs$contrasts) # Extract emmeans
freqCmodel_GRABDA2h_doseemmeans <- freqCmodel_GRABDA2h_dosepairs$emmeans # Extract emmeans
freqCmodel_GRABDA2h_dosepairs_es <- as.data.frame(eff_size(freqCmodel_GRABDA2h_doseemmeans, sigma = sigma(freqCmodel_GRABDA2h), edf = df.residual(freqCmodel_GRABDA2h))) # Effect Size
freqCmodel_GRABDA2h_dosepairs_es <- freqCmodel_GRABDA2h_dosepairs_es[, c("contrast", "effect.size")]

freqCmodel_GRABDA2h_dosecontrasts$d <- freqCmodel_GRABDA2h_dosepairs_es$effect.size
freqCmodel_GRABDA2h_dosecontrasts

# Format Tables
table_freqCmodel_GRABDA2h <- format_regsummary_table(freqCmodel_GRABDA2h, 'GRABDA2h', "Frequency", "StandardModel")
table_freqCmodel_GRABDA2h_aov <- format_anova_table(freqCmodel_GRABDA2h_aov, "GRABDA2h", "Frequency", "StandardModel")
table_freqCmodel_GRABDA2h_pairs_Dose <- format_emm_table(freqCmodel_GRABDA2h_dosecontrasts, "GRABDA2h", "Frequency", "Pairs","Dose")

#### Mixed effects AR1 model
freqCmodel_GRABDA2h_ar1 <- glmmTMB(pksperminC ~ Dose*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                     data = subset(subjectmeansbin, Virus=='GRABDA2h' & CondNum == 1& InjType=='M'), family = stats::gaussian, 
                                     control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqCmodel_GRABDA2h_ar1_aov <- Anova(freqCmodel_GRABDA2h_ar1, type=2)
freqCmodel_GRABDA2h_ar1_aov
summary(freqCmodel_GRABDA2h_ar1)

freqCmodel_GRABDA2h_ar1_r2 <- as.data.frame(r.squaredGLMM(freqCmodel_GRABDA2h_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqCmodel_GRABDA2h_ar1_r2

freqCmodel_GRABDA2h_ar1_AIC <- as.data.frame(AIC(freqCmodel_GRABDA2h, freqCmodel_GRABDA2h_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "freqCmodel_GRABDA2h" ~ "StandardModel", ModelName == "freqCmodel_GRABDA2h_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

freqCmodel_GRABDA2h_ar1_AIC

# Follow up tests - joint
freqCmodel_GRABDA2h_ar1_emmeans <-emmeans(freqCmodel_GRABDA2h_ar1, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
freqCmodel_GRABDA2h_ar1_dosepairs <- emmeans(freqCmodel_GRABDA2h_ar1, pairwise~Dose, adjust="none") # Pairs
freqCmodel_GRABDA2h_ar1_dosecontrasts <- as.data.frame(freqCmodel_GRABDA2h_ar1_dosepairs$contrasts) # Extract emmeans
freqCmodel_GRABDA2h_ar1_doseemmeans <- freqCmodel_GRABDA2h_ar1_dosepairs$emmeans # Extract emmeans
freqCmodel_GRABDA2h_ar1_dosepairs_es <- as.data.frame(eff_size(freqCmodel_GRABDA2h_ar1_doseemmeans, sigma = sigma(freqCmodel_GRABDA2h_ar1), edf = df.residual(freqCmodel_GRABDA2h_ar1))) # Effect Size
freqCmodel_GRABDA2h_ar1_dosepairs_es <- freqCmodel_GRABDA2h_ar1_dosepairs_es[, c("contrast", "effect.size")]

freqCmodel_GRABDA2h_ar1_dosecontrasts$d <- freqCmodel_GRABDA2h_ar1_dosepairs_es$effect.size
freqCmodel_GRABDA2h_ar1_dosecontrasts

# Format Tables
table_freqCmodel_GRABDA2h_ar1 <- format_regsummary_table(freqCmodel_GRABDA2h_ar1, 'GRABDA2h', "Frequency", "AR1Model")
table_freqCmodel_GRABDA2h_ar1_aov <- format_anova_table(freqCmodel_GRABDA2h_ar1_aov, "GRABDA2h", "Frequency", "AR1Model") 
table_freqCmodel_GRABDA2h_ar1_pairs_Dose <- format_emm_table(freqCmodel_GRABDA2h_ar1_dosecontrasts, "GRABDA2h", "Frequency", "Pairs", "Dose")

#### Max Height Model -----
## Mixed effects model
MaxHeightmodel_freq_GRABDA2h <- glmmTMB(MaxHeight ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                          data = subjectrise_pksperminC_GRABDA2h, family = stats::gaussian)

MaxHeightmodel_freq_GRABDA2h_aov <- Anova(MaxHeightmodel_freq_GRABDA2h, type=2)
MaxHeightmodel_freq_GRABDA2h_aov
summary(MaxHeightmodel_freq_GRABDA2h)

MaxHeightmodel_freq_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(MaxHeightmodel_freq_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxHeightmodel_freq_GRABDA2h_r2

## Prepare fit data for plotting
MaxHeight_freq_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_GRABDA2h$DoseNum)))
MaxHeight_freq_GRABDA2h_regdata$DoseNum_c <- MaxHeight_freq_GRABDA2h_regdata$DoseNum - dosemean

MaxHeight_freq_GRABDA2h_regpred <- predict(MaxHeightmodel_freq_GRABDA2h, newdata = MaxHeight_freq_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxHeight_freq_GRABDA2h_regdata$fit <- MaxHeight_freq_GRABDA2h_regpred$fit
MaxHeight_freq_GRABDA2h_regdata$se <- MaxHeight_freq_GRABDA2h_regpred$se
MaxHeight_freq_GRABDA2h_regdata$lower <- MaxHeight_freq_GRABDA2h_regdata$fit - 1.96 * MaxHeight_freq_GRABDA2h_regdata$se
MaxHeight_freq_GRABDA2h_regdata$upper <- MaxHeight_freq_GRABDA2h_regdata$fit + 1.96 * MaxHeight_freq_GRABDA2h_regdata$se

MaxHeight_freq_GRABDA2h_regdata

# Format Tables
table_freqmodel_GRABDA2h_MaxHeight <- format_regsummary_table(MaxHeightmodel_freq_GRABDA2h, 'GRABDA2h', "Frequency", "MaxHeightModel")
table_freqmodel_GRABDA2h_MaxHeight_aov <- format_anova_table(MaxHeightmodel_freq_GRABDA2h_aov, "GRABDA2h", "Frequency", "MaxHeightModel") 
table_freqmodel_GRABDA2h_MaxHeight_reg <- format_reg_table(MaxHeight_freq_GRABDA2h_regdata, "GRABDA2h", "Frequency", "MaxHeightReg")

#### Max Duration Model -----
## Mixed effects model
MaxDurationmodel_freq_GRABDA2h <- glmmTMB(MaxDuration~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                            data = subjectrise_pksperminC_GRABDA2h, family = stats::gaussian)

MaxDurationmodel_freq_GRABDA2h_aov<- Anova(MaxDurationmodel_freq_GRABDA2h, type=2)
MaxDurationmodel_freq_GRABDA2h_aov
summary(MaxDurationmodel_freq_GRABDA2h)

MaxDurationmodel_freq_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(MaxDurationmodel_freq_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxDurationmodel_freq_GRABDA2h_r2

## Prepare fit data for plotting
MaxDuration_freq_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_GRABDA2h$DoseNum)))
MaxDuration_freq_GRABDA2h_regdata$DoseNum_c <- MaxDuration_freq_GRABDA2h_regdata$DoseNum - dosemean

MaxDuration_freq_GRABDA2h_regpred <- predict(MaxDurationmodel_freq_GRABDA2h, newdata = MaxDuration_freq_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxDuration_freq_GRABDA2h_regdata$fit <- MaxDuration_freq_GRABDA2h_regpred$fit
MaxDuration_freq_GRABDA2h_regdata$se <- MaxDuration_freq_GRABDA2h_regpred$se
MaxDuration_freq_GRABDA2h_regdata$lower <- MaxDuration_freq_GRABDA2h_regdata$fit - 1.96 * MaxDuration_freq_GRABDA2h_regdata$se
MaxDuration_freq_GRABDA2h_regdata$upper <- MaxDuration_freq_GRABDA2h_regdata$fit + 1.96 * MaxDuration_freq_GRABDA2h_regdata$se

MaxDuration_freq_GRABDA2h_regdata

# Format Tables
table_freqmodel_GRABDA2h_MaxDuration <- format_regsummary_table(MaxDurationmodel_freq_GRABDA2h, 'GRABDA2h', "Frequency", "MaxDurationModel")
table_freqmodel_GRABDA2h_MaxDuration_aov <- format_anova_table(MaxDurationmodel_freq_GRABDA2h_aov, "GRABDA2h", "Frequency", "MaxDurationModel") 
table_freqmodel_GRABDA2h_MaxDuration_reg <- format_reg_table(MaxDuration_freq_GRABDA2h_regdata, "GRABDA2h", "Frequency", "MaxDurationReg")

#### Rise Slope Model -----
## Mixed effects model
MaxSlopeLMmodel_freq_GRABDA2h <- glmmTMB(MaxSlopeLM ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                           data = subjectrise_pksperminC_GRABDA2h, family = stats::gaussian)

MaxSlopeLMmodel_freq_GRABDA2h_aov <- Anova(MaxSlopeLMmodel_freq_GRABDA2h, type=2)
MaxSlopeLMmodel_freq_GRABDA2h_aov
summary(MaxSlopeLMmodel_freq_GRABDA2h)

MaxSlopeLMmodel_freq_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(MaxSlopeLMmodel_freq_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxSlopeLMmodel_freq_GRABDA2h_r2

## Prepare fit data for plotting
MaxSlopeLM_freq_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_GRABDA2h$DoseNum)))
MaxSlopeLM_freq_GRABDA2h_regdata$DoseNum_c <- MaxSlopeLM_freq_GRABDA2h_regdata$DoseNum - dosemean

MaxSlopeLM_freq_GRABDA2h_regpred <- predict(MaxSlopeLMmodel_freq_GRABDA2h, newdata = MaxSlopeLM_freq_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxSlopeLM_freq_GRABDA2h_regdata$fit <- MaxSlopeLM_freq_GRABDA2h_regpred$fit
MaxSlopeLM_freq_GRABDA2h_regdata$se <- MaxSlopeLM_freq_GRABDA2h_regpred$se
MaxSlopeLM_freq_GRABDA2h_regdata$lower <- MaxSlopeLM_freq_GRABDA2h_regdata$fit - 1.96 * MaxSlopeLM_freq_GRABDA2h_regdata$se
MaxSlopeLM_freq_GRABDA2h_regdata$upper <- MaxSlopeLM_freq_GRABDA2h_regdata$fit + 1.96 * MaxSlopeLM_freq_GRABDA2h_regdata$se

MaxSlopeLM_freq_GRABDA2h_regdata

# Format Tables
table_freqmodel_GRABDA2h_MaxSlopeLM <- format_regsummary_table(MaxSlopeLMmodel_freq_GRABDA2h, 'GRABDA2h', "Frequency", "MaxSlopeLMModel")
table_freqmodel_GRABDA2h_MaxSlopeLM_aov <- format_anova_table(MaxSlopeLMmodel_freq_GRABDA2h_aov, "GRABDA2h", "Frequency", "MaxSlopeLMModel") 
table_freqmodel_GRABDA2h_MaxSlopeLM_reg <- format_reg_table(MaxSlopeLM_freq_GRABDA2h_regdata, "GRABDA2h", "Frequency", "MaxSlopeLMReg")

#### AUC Model -----
## Mixed effects model
AUCmodel_freq_GRABDA2h <- glmmTMB(AUC_POSTINJ ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                    data = subset(subjectrise_pksperminC_GRABDA2h, Virus=='GRABDA2h'&InjType=='M'), family = stats::gaussian)

AUCmodel_freq_GRABDA2h_aov <- Anova(AUCmodel_freq_GRABDA2h, type=2)
AUCmodel_freq_GRABDA2h_aov
summary(AUCmodel_freq_GRABDA2h)

AUCmodel_freq_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(AUCmodel_freq_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCmodel_freq_GRABDA2h_r2

## Prepare fit data for plotting
AUC_freq_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_pksperminC_GRABDA2h$DoseNum)))
AUC_freq_GRABDA2h_regdata$DoseNum_c <- AUC_freq_GRABDA2h_regdata$DoseNum - dosemean

AUC_freq_GRABDA2h_regpred <- predict(AUCmodel_freq_GRABDA2h, newdata = AUC_freq_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
AUC_freq_GRABDA2h_regdata$fit <- AUC_freq_GRABDA2h_regpred$fit
AUC_freq_GRABDA2h_regdata$se <- AUC_freq_GRABDA2h_regpred$se
AUC_freq_GRABDA2h_regdata$lower <- AUC_freq_GRABDA2h_regdata$fit - 1.96 * AUC_freq_GRABDA2h_regdata$se
AUC_freq_GRABDA2h_regdata$upper <- AUC_freq_GRABDA2h_regdata$fit + 1.96 * AUC_freq_GRABDA2h_regdata$se

AUC_freq_GRABDA2h_regdata

# Format Tables
table_freqmodel_GRABDA2h_AUC <- format_regsummary_table(AUCmodel_freq_GRABDA2h, 'GRABDA2h', "Frequency", "AUCModel")
table_freqmodel_GRABDA2h_AUC_aov <- format_anova_table(AUCmodel_freq_GRABDA2h_aov, "GRABDA2h", "Frequency", "AUCModel") 
table_freqmodel_GRABDA2h_AUC_reg <- format_reg_table(AUC_freq_GRABDA2h_regdata, "GRABDA2h", "Frequency", "AUCReg")

### Amplitude -----
#### Amplitude Model -----
# Mixed effects model
ampmodel_GRABDA2h <- glmmTMB(amp ~ InjType*Dose*BinF + (1|SubjectID), # random intercept by subject
                               data = subset(pkdataz, Virus=='GRABDA2h' & CondNum == 1), family = stats::gaussian, 
                               control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampmodel_GRABDA2h_aov <- Anova(ampmodel_GRABDA2h, type=2)
ampmodel_GRABDA2h_aov
summary(ampmodel_GRABDA2h)

ampmodel_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(ampmodel_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampmodel_GRABDA2h_r2

# Follow up tests - joint
ampmodel_GRABDA2h_emmeans <-emmeans(ampmodel_GRABDA2h, pairwise~InjType|Dose|BinF, adjust="none") # Find means
ampmodel_GRABDA2h_jt_Dose <- joint_tests(ampmodel_GRABDA2h_emmeans, by = "Dose")
ampmodel_GRABDA2h_jt_InjType <- joint_tests(ampmodel_GRABDA2h_emmeans, by = "InjType")

ampmodel_GRABDA2h_jt_Dose
ampmodel_GRABDA2h_jt_InjType

# Follow up tests - treat pairs
ampmodel_GRABDA2h_treatpairs <- emmeans(ampmodel_GRABDA2h, pairwise~InjType|Dose, adjust="none") # Pairs
ampmodel_GRABDA2h_treatcontrasts <- as.data.frame(ampmodel_GRABDA2h_treatpairs$contrasts) # Extract emmeans
ampmodel_GRABDA2h_treatemmeans <- ampmodel_GRABDA2h_treatpairs$emmeans # Extract emmeans
ampmodel_GRABDA2h_treatpairs_es <- as.data.frame(eff_size(ampmodel_GRABDA2h_treatemmeans, sigma = sigma(ampmodel_GRABDA2h), edf = df.residual(ampmodel_GRABDA2h))) # Effect Size
ampmodel_GRABDA2h_treatpairs_es <- ampmodel_GRABDA2h_treatpairs_es[, c("Dose", "contrast", "effect.size")]

ampmodel_GRABDA2h_treatcontrasts$d <- ampmodel_GRABDA2h_treatpairs_es$effect.size
ampmodel_GRABDA2h_treatcontrasts

# Follow up tests - dose pairs
ampmodel_GRABDA2h_dosepairs <- emmeans(ampmodel_GRABDA2h, pairwise~Dose|InjType, adjust="none") # Pairs
ampmodel_GRABDA2h_dosecontrasts <- as.data.frame(ampmodel_GRABDA2h_dosepairs$contrasts) # Extract emmeans
ampmodel_GRABDA2h_doseemmeans <- ampmodel_GRABDA2h_dosepairs$emmeans # Extract emmeans
ampmodel_GRABDA2h_dosepairs_es <- as.data.frame(eff_size(ampmodel_GRABDA2h_doseemmeans, sigma = sigma(ampmodel_GRABDA2h), edf = df.residual(ampmodel_GRABDA2h))) # Effect Size
ampmodel_GRABDA2h_dosepairs_es <- ampmodel_GRABDA2h_dosepairs_es[, c("InjType", "contrast", "effect.size")]

ampmodel_GRABDA2h_dosecontrasts$d <- ampmodel_GRABDA2h_dosepairs_es$effect.size
ampmodel_GRABDA2h_dosecontrasts

# Format Tables
table_ampmodel_GRABDA2h <- format_regsummary_table(ampmodel_GRABDA2h, 'GRABDA2h', "Amplitude", "StandardModel")
table_ampmodel_GRABDA2h_aov <- format_anova_table(ampmodel_GRABDA2h_aov, "GRABDA2h", "Amplitude", "StandardModel")
table_ampmodel_GRABDA2h_joint_Dose <- format_emm_table(ampmodel_GRABDA2h_jt_Dose, "GRABDA2h", "Amplitude", "Joint", "Dose")
table_ampmodel_GRABDA2h_joint_InjType <- format_emm_table(ampmodel_GRABDA2h_jt_InjType, "GRABDA2h", "Amplitude", "Joint", "InjType")
table_ampmodel_GRABDA2h_pairs_Dose <- format_emm_table(ampmodel_GRABDA2h_dosecontrasts, "GRABDA2h", "Amplitude", "Pairs", "Dose|InjType")
table_ampmodel_GRABDA2h_pairs_InjType <- format_emm_table(ampmodel_GRABDA2h_treatcontrasts, "GRABDA2h", "Amplitude", "Pairs", "InjType|Dose")

#### Mixed effects AR1 model
ampmodel_GRABDA2h_ar1 <- glmmTMB(amp ~ InjType*Dose*BinF + ar1(BinF + 0|SubjectID:InjType) + (1|SubjectID), # random intercept by subject
                                   data = subset(pkdataz, Virus=='GRABDA2h' & CondNum == 1), family = stats::gaussian, 
                                   control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampmodel_GRABDA2h_ar1_aov <- Anova(ampmodel_GRABDA2h_ar1, type=2)
ampmodel_GRABDA2h_ar1_aov
summary(ampmodel_GRABDA2h_ar1)

ampmodel_GRABDA2h_ar1_r2 <- as.data.frame(r.squaredGLMM(ampmodel_GRABDA2h_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampmodel_GRABDA2h_ar1_r2

ampmodel_GRABDA2h_ar1_AIC <- as.data.frame(AIC(ampmodel_GRABDA2h, ampmodel_GRABDA2h_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "ampmodel_GRABDA2h" ~ "StandardModel", ModelName == "ampmodel_GRABDA2h_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

ampmodel_GRABDA2h_ar1_AIC

# Follow up tests - joint
ampmodel_GRABDA2h_ar1_emmeans <-emmeans(ampmodel_GRABDA2h_ar1, pairwise~InjType|Dose, adjust="none") # Find means
ampmodel_GRABDA2h_ar1_jt_Dose <- joint_tests(ampmodel_GRABDA2h_ar1, by = "Dose")
ampmodel_GRABDA2h_ar1_jt_InjType <- joint_tests(ampmodel_GRABDA2h_ar1, by = "InjType")

ampmodel_GRABDA2h_ar1_jt_Dose
ampmodel_GRABDA2h_ar1_jt_InjType

# Follow up tests - treat pairs
ampmodel_GRABDA2h_ar1_treatpairs <- emmeans(ampmodel_GRABDA2h_ar1, pairwise~InjType|Dose, adjust="none") # Pairs
ampmodel_GRABDA2h_ar1_treatcontrasts <- as.data.frame(ampmodel_GRABDA2h_ar1_treatpairs$contrasts) # Extract emmeans
ampmodel_GRABDA2h_ar1_treatemmeans <- ampmodel_GRABDA2h_ar1_treatpairs$emmeans # Extract emmeans
ampmodel_GRABDA2h_ar1_treatpairs_es <- as.data.frame(eff_size(ampmodel_GRABDA2h_ar1_treatemmeans, sigma = sigma(ampmodel_GRABDA2h_ar1), edf = df.residual(ampmodel_GRABDA2h_ar1))) # Effect Size
ampmodel_GRABDA2h_ar1_treatpairs_es <- ampmodel_GRABDA2h_ar1_treatpairs_es[, c("Dose", "contrast", "effect.size")]

ampmodel_GRABDA2h_ar1_treatcontrasts$d <- ampmodel_GRABDA2h_ar1_treatpairs_es$effect.size
ampmodel_GRABDA2h_ar1_treatcontrasts

# Follow up tests - dose pairs
ampmodel_GRABDA2h_ar1_dosepairs <- emmeans(ampmodel_GRABDA2h_ar1, pairwise~Dose|InjType, adjust="none") # Pairs
ampmodel_GRABDA2h_ar1_dosecontrasts <- as.data.frame(ampmodel_GRABDA2h_ar1_dosepairs$contrasts) # Extract emmeans
ampmodel_GRABDA2h_ar1_doseemmeans <- ampmodel_GRABDA2h_ar1_dosepairs$emmeans # Extract emmeans
ampmodel_GRABDA2h_ar1_dosepairs_es <- as.data.frame(eff_size(ampmodel_GRABDA2h_ar1_doseemmeans, sigma = sigma(ampmodel_GRABDA2h_ar1), edf = df.residual(ampmodel_GRABDA2h_ar1))) # Effect Size
ampmodel_GRABDA2h_ar1_dosepairs_es <- ampmodel_GRABDA2h_ar1_dosepairs_es[, c("InjType", "contrast", "effect.size")]

ampmodel_GRABDA2h_ar1_dosecontrasts$d <- ampmodel_GRABDA2h_ar1_dosepairs_es$effect.size
ampmodel_GRABDA2h_ar1_dosecontrasts

# Format Tables
table_ampmodel_GRABDA2h_ar1 <- format_regsummary_table(ampmodel_GRABDA2h_ar1, 'GRABDA2h', "Amplitude", "AR1Model")
table_ampmodel_GRABDA2h_ar1_aov <- format_anova_table(ampmodel_GRABDA2h_ar1_aov, "GRABDA2h", "Amplitude", "AR1Model") 
table_ampmodel_GRABDA2h_ar1_joint_Dose <- format_emm_table(ampmodel_GRABDA2h_ar1_jt_Dose, "GRABDA2h", "Amplitude", "Joint", "Dose")
table_ampmodel_GRABDA2h_ar1_joint_InjType <- format_emm_table(ampmodel_GRABDA2h_ar1_jt_InjType, "GRABDA2h", "Amplitude", "Joint", "InjType")
table_ampmodel_GRABDA2h_ar1_pairs_Dose <- format_emm_table(ampmodel_GRABDA2h_ar1_dosecontrasts, "GRABDA2h", "Amplitude", "Pairs", "Dose|InjType")
table_ampmodel_GRABDA2h_ar1_pairs_InjType <- format_emm_table(ampmodel_GRABDA2h_ar1_treatcontrasts, "GRABDA2h", "Amplitude", "Pairs", "InjType|Dose")

#### Amplitude Change Model -----
# Mixed effects model
ampCmodel_GRABDA2h <- glmmTMB(ampC ~ Dose*BinF + (1|SubjectID), # random intercept by subject
                                data = subset(subjectmeansbin, Virus=='GRABDA2h' & CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampCmodel_GRABDA2h_aov <- Anova(ampCmodel_GRABDA2h, type=2)
ampCmodel_GRABDA2h_aov
summary(ampCmodel_GRABDA2h)

ampCmodel_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(ampCmodel_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampCmodel_GRABDA2h_r2

# Follow up tests - joint
ampCmodel_GRABDA2h_emmeans <-emmeans(ampCmodel_GRABDA2h, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
ampCmodel_GRABDA2h_dosepairs <- emmeans(ampCmodel_GRABDA2h, pairwise~Dose, adjust="none") # Pairs
ampCmodel_GRABDA2h_dosecontrasts <- as.data.frame(ampCmodel_GRABDA2h_dosepairs$contrasts) # Extract emmeans
ampCmodel_GRABDA2h_doseemmeans <- ampCmodel_GRABDA2h_dosepairs$emmeans # Extract emmeans
ampCmodel_GRABDA2h_dosepairs_es <- as.data.frame(eff_size(ampCmodel_GRABDA2h_doseemmeans, sigma = sigma(ampCmodel_GRABDA2h), edf = df.residual(ampCmodel_GRABDA2h))) # Effect Size
ampCmodel_GRABDA2h_dosepairs_es <- ampCmodel_GRABDA2h_dosepairs_es[, c("contrast", "effect.size")]

ampCmodel_GRABDA2h_dosecontrasts$d <- ampCmodel_GRABDA2h_dosepairs_es$effect.size
ampCmodel_GRABDA2h_dosecontrasts

# Format Tables
table_ampCmodel_GRABDA2h <- format_regsummary_table(ampCmodel_GRABDA2h, 'GRABDA2h', "Amplitude", "StandardModel")
table_ampCmodel_GRABDA2h_aov <- format_anova_table(ampCmodel_GRABDA2h_aov, "GRABDA2h", "Amplitude", "StandardModel")
table_ampCmodel_GRABDA2h_pairs_Dose <- format_emm_table(ampCmodel_GRABDA2h_dosecontrasts, "GRABDA2h", "Amplitude", "Pairs","Dose")

#### Mixed effects AR1 model
ampCmodel_GRABDA2h_ar1 <- glmmTMB(ampC ~ Dose*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                    data = subset(subjectmeansbin, Virus=='GRABDA2h' & CondNum == 1& InjType=='M'), family = stats::gaussian, 
                                    control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampCmodel_GRABDA2h_ar1_aov <- Anova(ampCmodel_GRABDA2h_ar1, type=2)
ampCmodel_GRABDA2h_ar1_aov
summary(ampCmodel_GRABDA2h_ar1)

ampCmodel_GRABDA2h_ar1_r2 <- as.data.frame(r.squaredGLMM(ampCmodel_GRABDA2h_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampCmodel_GRABDA2h_ar1_r2

ampCmodel_GRABDA2h_ar1_AIC <- as.data.frame(AIC(ampCmodel_GRABDA2h, ampCmodel_GRABDA2h_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "ampCmodel_GRABDA2h" ~ "StandardModel", ModelName == "ampCmodel_GRABDA2h_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

ampCmodel_GRABDA2h_ar1_AIC

# Follow up tests - joint
ampCmodel_GRABDA2h_ar1_emmeans <-emmeans(ampCmodel_GRABDA2h_ar1, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
ampCmodel_GRABDA2h_ar1_dosepairs <- emmeans(ampCmodel_GRABDA2h_ar1, pairwise~Dose, adjust="none") # Pairs
ampCmodel_GRABDA2h_ar1_dosecontrasts <- as.data.frame(ampCmodel_GRABDA2h_ar1_dosepairs$contrasts) # Extract emmeans
ampCmodel_GRABDA2h_ar1_doseemmeans <- ampCmodel_GRABDA2h_ar1_dosepairs$emmeans # Extract emmeans
ampCmodel_GRABDA2h_ar1_dosepairs_es <- as.data.frame(eff_size(ampCmodel_GRABDA2h_ar1_doseemmeans, sigma = sigma(ampCmodel_GRABDA2h_ar1), edf = df.residual(ampCmodel_GRABDA2h_ar1))) # Effect Size
ampCmodel_GRABDA2h_ar1_dosepairs_es <- ampCmodel_GRABDA2h_ar1_dosepairs_es[, c("contrast", "effect.size")]

ampCmodel_GRABDA2h_ar1_dosecontrasts$d <- ampCmodel_GRABDA2h_ar1_dosepairs_es$effect.size
ampCmodel_GRABDA2h_ar1_dosecontrasts

# Format Tables
table_ampCmodel_GRABDA2h_ar1 <- format_regsummary_table(ampCmodel_GRABDA2h_ar1, 'GRABDA2h', "Amplitude", "AR1Model")
table_ampCmodel_GRABDA2h_ar1_aov <- format_anova_table(ampCmodel_GRABDA2h_ar1_aov, "GRABDA2h", "Amplitude", "AR1Model") 
table_ampCmodel_GRABDA2h_ar1_pairs_Dose <- format_emm_table(ampCmodel_GRABDA2h_ar1_dosecontrasts, "GRABDA2h", "Amplitude", "Pairs", "Dose")


#### Max Height Model -----
## Mixed effects model
MaxHeightmodel_amp_GRABDA2h <- glmmTMB(MaxHeight ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                         data = subjectrise_ampC_GRABDA2h, family = stats::gaussian)

MaxHeightmodel_amp_GRABDA2h_aov <- Anova(MaxHeightmodel_amp_GRABDA2h, type=2)
MaxHeightmodel_amp_GRABDA2h_aov
summary(MaxHeightmodel_amp_GRABDA2h)

MaxHeightmodel_amp_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(MaxHeightmodel_amp_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxHeightmodel_amp_GRABDA2h_r2

## Prepare fit data for plotting
MaxHeight_amp_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_GRABDA2h$DoseNum)))
MaxHeight_amp_GRABDA2h_regdata$DoseNum_c <- MaxHeight_amp_GRABDA2h_regdata$DoseNum - dosemean

MaxHeight_amp_GRABDA2h_regpred <- predict(MaxHeightmodel_amp_GRABDA2h, newdata = MaxHeight_amp_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxHeight_amp_GRABDA2h_regdata$fit <- MaxHeight_amp_GRABDA2h_regpred$fit
MaxHeight_amp_GRABDA2h_regdata$se <- MaxHeight_amp_GRABDA2h_regpred$se
MaxHeight_amp_GRABDA2h_regdata$lower <- MaxHeight_amp_GRABDA2h_regdata$fit - 1.96 * MaxHeight_amp_GRABDA2h_regdata$se
MaxHeight_amp_GRABDA2h_regdata$upper <- MaxHeight_amp_GRABDA2h_regdata$fit + 1.96 * MaxHeight_amp_GRABDA2h_regdata$se

MaxHeight_amp_GRABDA2h_regdata

# Format Tables
table_ampmodel_GRABDA2h_MaxHeight <- format_regsummary_table(MaxHeightmodel_amp_GRABDA2h, 'GRABDA2h', "Amplitude", "MaxHeightModel")
table_ampmodel_GRABDA2h_MaxHeight_aov <- format_anova_table(MaxHeightmodel_amp_GRABDA2h_aov, "GRABDA2h", "Amplitude", "MaxHeightModel") 
table_ampmodel_GRABDA2h_MaxHeight_reg <- format_reg_table(MaxHeight_amp_GRABDA2h_regdata, "GRABDA2h", "Amplitude", "MaxHeightReg")

#### Max Duration Model -----
## Mixed effects model
MaxDurationmodel_amp_GRABDA2h <- glmmTMB(MaxDuration~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                           data = subjectrise_ampC_GRABDA2h, family = stats::gaussian)

MaxDurationmodel_amp_GRABDA2h_aov<- Anova(MaxDurationmodel_amp_GRABDA2h, type=2)
MaxDurationmodel_amp_GRABDA2h_aov
summary(MaxDurationmodel_amp_GRABDA2h)

MaxDurationmodel_amp_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(MaxDurationmodel_amp_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxDurationmodel_amp_GRABDA2h_r2

## Prepare fit data for plotting
MaxDuration_amp_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_GRABDA2h$DoseNum)))
MaxDuration_amp_GRABDA2h_regdata$DoseNum_c <- MaxDuration_amp_GRABDA2h_regdata$DoseNum - dosemean

MaxDuration_amp_GRABDA2h_regpred <- predict(MaxDurationmodel_amp_GRABDA2h, newdata = MaxDuration_amp_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxDuration_amp_GRABDA2h_regdata$fit <- MaxDuration_amp_GRABDA2h_regpred$fit
MaxDuration_amp_GRABDA2h_regdata$se <- MaxDuration_amp_GRABDA2h_regpred$se
MaxDuration_amp_GRABDA2h_regdata$lower <- MaxDuration_amp_GRABDA2h_regdata$fit - 1.96 * MaxDuration_amp_GRABDA2h_regdata$se
MaxDuration_amp_GRABDA2h_regdata$upper <- MaxDuration_amp_GRABDA2h_regdata$fit + 1.96 * MaxDuration_amp_GRABDA2h_regdata$se

MaxDuration_amp_GRABDA2h_regdata

# Format Tables
table_ampmodel_GRABDA2h_MaxDuration <- format_regsummary_table(MaxDurationmodel_amp_GRABDA2h, 'GRABDA2h', "Amplitude", "MaxDurationModel")
table_ampmodel_GRABDA2h_MaxDuration_aov <- format_anova_table(MaxDurationmodel_amp_GRABDA2h_aov, "GRABDA2h", "Amplitude", "MaxDurationModel") 
table_ampmodel_GRABDA2h_MaxDuration_reg <- format_reg_table(MaxDuration_amp_GRABDA2h_regdata, "GRABDA2h", "Amplitude", "MaxDurationReg")

#### Rise Slope Model -----
## Mixed effects model
MaxSlopeLMmodel_amp_GRABDA2h <- glmmTMB(MaxSlopeLM ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                          data = subjectrise_ampC_GRABDA2h, family = stats::gaussian)

MaxSlopeLMmodel_amp_GRABDA2h_aov <- Anova(MaxSlopeLMmodel_amp_GRABDA2h, type=2)
MaxSlopeLMmodel_amp_GRABDA2h_aov
summary(MaxSlopeLMmodel_amp_GRABDA2h)

MaxSlopeLMmodel_amp_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(MaxSlopeLMmodel_amp_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxSlopeLMmodel_amp_GRABDA2h_r2

## Prepare fit data for plotting
MaxSlopeLM_amp_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_GRABDA2h$DoseNum)))
MaxSlopeLM_amp_GRABDA2h_regdata$DoseNum_c <- MaxSlopeLM_amp_GRABDA2h_regdata$DoseNum - dosemean

MaxSlopeLM_amp_GRABDA2h_regpred <- predict(MaxSlopeLMmodel_amp_GRABDA2h, newdata = MaxSlopeLM_amp_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxSlopeLM_amp_GRABDA2h_regdata$fit <- MaxSlopeLM_amp_GRABDA2h_regpred$fit
MaxSlopeLM_amp_GRABDA2h_regdata$se <- MaxSlopeLM_amp_GRABDA2h_regpred$se
MaxSlopeLM_amp_GRABDA2h_regdata$lower <- MaxSlopeLM_amp_GRABDA2h_regdata$fit - 1.96 * MaxSlopeLM_amp_GRABDA2h_regdata$se
MaxSlopeLM_amp_GRABDA2h_regdata$upper <- MaxSlopeLM_amp_GRABDA2h_regdata$fit + 1.96 * MaxSlopeLM_amp_GRABDA2h_regdata$se

MaxSlopeLM_amp_GRABDA2h_regdata

# Format Tables
table_ampmodel_GRABDA2h_MaxSlopeLM <- format_regsummary_table(MaxSlopeLMmodel_amp_GRABDA2h, 'GRABDA2h', "Amplitude", "MaxSlopeLMModel")
table_ampmodel_GRABDA2h_MaxSlopeLM_aov <- format_anova_table(MaxSlopeLMmodel_amp_GRABDA2h_aov, "GRABDA2h", "Amplitude", "MaxSlopeLMModel") 
table_ampmodel_GRABDA2h_MaxSlopeLM_reg <- format_reg_table(MaxSlopeLM_amp_GRABDA2h_regdata, "GRABDA2h", "Amplitude", "MaxSlopeLMReg")

#### AUC Model -----
## Mixed effects model
AUCmodel_amp_GRABDA2h <- glmmTMB(AUC_POSTINJ ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                   data = subset(subjectrise_ampC_GRABDA2h, Virus=='GRABDA2h'&InjType=='M'), family = stats::gaussian)

AUCmodel_amp_GRABDA2h_aov <- Anova(AUCmodel_amp_GRABDA2h, type=2)
AUCmodel_amp_GRABDA2h_aov
summary(AUCmodel_amp_GRABDA2h)

AUCmodel_amp_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(AUCmodel_amp_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCmodel_amp_GRABDA2h_r2

## Prepare fit data for plotting
AUC_amp_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_ampC_GRABDA2h$DoseNum)))
AUC_amp_GRABDA2h_regdata$DoseNum_c <- AUC_amp_GRABDA2h_regdata$DoseNum - dosemean

AUC_amp_GRABDA2h_regpred <- predict(AUCmodel_amp_GRABDA2h, newdata = AUC_amp_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
AUC_amp_GRABDA2h_regdata$fit <- AUC_amp_GRABDA2h_regpred$fit
AUC_amp_GRABDA2h_regdata$se <- AUC_amp_GRABDA2h_regpred$se
AUC_amp_GRABDA2h_regdata$lower <- AUC_amp_GRABDA2h_regdata$fit - 1.96 * AUC_amp_GRABDA2h_regdata$se
AUC_amp_GRABDA2h_regdata$upper <- AUC_amp_GRABDA2h_regdata$fit + 1.96 * AUC_amp_GRABDA2h_regdata$se

AUC_amp_GRABDA2h_regdata

# Format Tables
table_ampmodel_GRABDA2h_AUC <- format_regsummary_table(AUCmodel_amp_GRABDA2h, 'GRABDA2h', "Amplitude", "AUCModel")
table_ampmodel_GRABDA2h_AUC_aov <- format_anova_table(AUCmodel_amp_GRABDA2h_aov, "GRABDA2h", "Amplitude", "AUCModel") 
table_ampmodel_GRABDA2h_AUC_reg <- format_reg_table(AUC_amp_GRABDA2h_regdata, "GRABDA2h", "Amplitude", "AUCReg")

### AUCwindow -----
#### AUCwindow Model -----
# Mixed effects model
AUCwindowmodel_GRABDA2h <- glmmTMB(AUCwindow ~ InjType*Dose*BinF + (1|SubjectID), # random intercept by subject
                                     data = subset(pkdataz, Virus=='GRABDA2h' & CondNum == 1), family = stats::gaussian, 
                                     control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowmodel_GRABDA2h_aov <- Anova(AUCwindowmodel_GRABDA2h, type=2)
AUCwindowmodel_GRABDA2h_aov
summary(AUCwindowmodel_GRABDA2h)

AUCwindowmodel_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(AUCwindowmodel_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowmodel_GRABDA2h_r2

# Follow up tests - joint
AUCwindowmodel_GRABDA2h_emmeans <-emmeans(AUCwindowmodel_GRABDA2h, pairwise~InjType|Dose|BinF, adjust="none") # Find means
AUCwindowmodel_GRABDA2h_jt_Dose <- joint_tests(AUCwindowmodel_GRABDA2h_emmeans, by = "Dose")
AUCwindowmodel_GRABDA2h_jt_InjType <- joint_tests(AUCwindowmodel_GRABDA2h_emmeans, by = "InjType")

AUCwindowmodel_GRABDA2h_jt_Dose
AUCwindowmodel_GRABDA2h_jt_InjType

# Follow up tests - treat pairs
AUCwindowmodel_GRABDA2h_treatpairs <- emmeans(AUCwindowmodel_GRABDA2h, pairwise~InjType|Dose, adjust="none") # Pairs
AUCwindowmodel_GRABDA2h_treatcontrasts <- as.data.frame(AUCwindowmodel_GRABDA2h_treatpairs$contrasts) # Extract emmeans
AUCwindowmodel_GRABDA2h_treatemmeans <- AUCwindowmodel_GRABDA2h_treatpairs$emmeans # Extract emmeans
AUCwindowmodel_GRABDA2h_treatpairs_es <- as.data.frame(eff_size(AUCwindowmodel_GRABDA2h_treatemmeans, sigma = sigma(AUCwindowmodel_GRABDA2h), edf = df.residual(AUCwindowmodel_GRABDA2h))) # Effect Size
AUCwindowmodel_GRABDA2h_treatpairs_es <- AUCwindowmodel_GRABDA2h_treatpairs_es[, c("Dose", "contrast", "effect.size")]

AUCwindowmodel_GRABDA2h_treatcontrasts$d <- AUCwindowmodel_GRABDA2h_treatpairs_es$effect.size
AUCwindowmodel_GRABDA2h_treatcontrasts

# Follow up tests - dose pairs
AUCwindowmodel_GRABDA2h_dosepairs <- emmeans(AUCwindowmodel_GRABDA2h, pairwise~Dose|InjType, adjust="none") # Pairs
AUCwindowmodel_GRABDA2h_dosecontrasts <- as.data.frame(AUCwindowmodel_GRABDA2h_dosepairs$contrasts) # Extract emmeans
AUCwindowmodel_GRABDA2h_doseemmeans <- AUCwindowmodel_GRABDA2h_dosepairs$emmeans # Extract emmeans
AUCwindowmodel_GRABDA2h_dosepairs_es <- as.data.frame(eff_size(AUCwindowmodel_GRABDA2h_doseemmeans, sigma = sigma(AUCwindowmodel_GRABDA2h), edf = df.residual(AUCwindowmodel_GRABDA2h))) # Effect Size
AUCwindowmodel_GRABDA2h_dosepairs_es <- AUCwindowmodel_GRABDA2h_dosepairs_es[, c("InjType", "contrast", "effect.size")]

AUCwindowmodel_GRABDA2h_dosecontrasts$d <- AUCwindowmodel_GRABDA2h_dosepairs_es$effect.size
AUCwindowmodel_GRABDA2h_dosecontrasts

# Format Tables
table_AUCwindowmodel_GRABDA2h <- format_regsummary_table(AUCwindowmodel_GRABDA2h, 'GRABDA2h', "AUCwindow", "StandardModel")
table_AUCwindowmodel_GRABDA2h_aov <- format_anova_table(AUCwindowmodel_GRABDA2h_aov, "GRABDA2h", "AUCwindow", "StandardModel")
table_AUCwindowmodel_GRABDA2h_joint_Dose <- format_emm_table(AUCwindowmodel_GRABDA2h_jt_Dose, "GRABDA2h", "AUCwindow", "Joint", "Dose")
table_AUCwindowmodel_GRABDA2h_joint_InjType <- format_emm_table(AUCwindowmodel_GRABDA2h_jt_InjType, "GRABDA2h", "AUCwindow", "Joint", "InjType")
table_AUCwindowmodel_GRABDA2h_pairs_Dose <- format_emm_table(AUCwindowmodel_GRABDA2h_dosecontrasts, "GRABDA2h", "AUCwindow", "Pairs", "Dose|InjType")
table_AUCwindowmodel_GRABDA2h_pairs_InjType <- format_emm_table(AUCwindowmodel_GRABDA2h_treatcontrasts, "GRABDA2h", "AUCwindow", "Pairs", "InjType|Dose")

#### Mixed effects AR1 model
AUCwindowmodel_GRABDA2h_ar1 <- glmmTMB(AUCwindow ~ InjType*Dose*BinF + ar1(BinF + 0|SubjectID:InjType) + (1|SubjectID), # random intercept by subject
                                         data = subset(pkdataz, Virus=='GRABDA2h' & CondNum == 1), family = stats::gaussian, 
                                         control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowmodel_GRABDA2h_ar1_aov <- Anova(AUCwindowmodel_GRABDA2h_ar1, type=2)
AUCwindowmodel_GRABDA2h_ar1_aov
summary(AUCwindowmodel_GRABDA2h_ar1)

AUCwindowmodel_GRABDA2h_ar1_r2 <- as.data.frame(r.squaredGLMM(AUCwindowmodel_GRABDA2h_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowmodel_GRABDA2h_ar1_r2

AUCwindowmodel_GRABDA2h_ar1_AIC <- as.data.frame(AIC(AUCwindowmodel_GRABDA2h, AUCwindowmodel_GRABDA2h_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "AUCwindowmodel_GRABDA2h" ~ "StandardModel", ModelName == "AUCwindowmodel_GRABDA2h_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

AUCwindowmodel_GRABDA2h_ar1_AIC

# Follow up tests - joint
AUCwindowmodel_GRABDA2h_ar1_emmeans <-emmeans(AUCwindowmodel_GRABDA2h_ar1, pairwise~InjType|Dose, adjust="none") # Find means
AUCwindowmodel_GRABDA2h_ar1_jt_Dose <- joint_tests(AUCwindowmodel_GRABDA2h_ar1, by = "Dose")
AUCwindowmodel_GRABDA2h_ar1_jt_InjType <- joint_tests(AUCwindowmodel_GRABDA2h_ar1, by = "InjType")

AUCwindowmodel_GRABDA2h_ar1_jt_Dose
AUCwindowmodel_GRABDA2h_ar1_jt_InjType

# Follow up tests - treat pairs
AUCwindowmodel_GRABDA2h_ar1_treatpairs <- emmeans(AUCwindowmodel_GRABDA2h_ar1, pairwise~InjType|Dose, adjust="none") # Pairs
AUCwindowmodel_GRABDA2h_ar1_treatcontrasts <- as.data.frame(AUCwindowmodel_GRABDA2h_ar1_treatpairs$contrasts) # Extract emmeans
AUCwindowmodel_GRABDA2h_ar1_treatemmeans <- AUCwindowmodel_GRABDA2h_ar1_treatpairs$emmeans # Extract emmeans
AUCwindowmodel_GRABDA2h_ar1_treatpairs_es <- as.data.frame(eff_size(AUCwindowmodel_GRABDA2h_ar1_treatemmeans, sigma = sigma(AUCwindowmodel_GRABDA2h_ar1), edf = df.residual(AUCwindowmodel_GRABDA2h_ar1))) # Effect Size
AUCwindowmodel_GRABDA2h_ar1_treatpairs_es <- AUCwindowmodel_GRABDA2h_ar1_treatpairs_es[, c("Dose", "contrast", "effect.size")]

AUCwindowmodel_GRABDA2h_ar1_treatcontrasts$d <- AUCwindowmodel_GRABDA2h_ar1_treatpairs_es$effect.size
AUCwindowmodel_GRABDA2h_ar1_treatcontrasts

# Follow up tests - dose pairs
AUCwindowmodel_GRABDA2h_ar1_dosepairs <- emmeans(AUCwindowmodel_GRABDA2h_ar1, pairwise~Dose|InjType, adjust="none") # Pairs
AUCwindowmodel_GRABDA2h_ar1_dosecontrasts <- as.data.frame(AUCwindowmodel_GRABDA2h_ar1_dosepairs$contrasts) # Extract emmeans
AUCwindowmodel_GRABDA2h_ar1_doseemmeans <- AUCwindowmodel_GRABDA2h_ar1_dosepairs$emmeans # Extract emmeans
AUCwindowmodel_GRABDA2h_ar1_dosepairs_es <- as.data.frame(eff_size(AUCwindowmodel_GRABDA2h_ar1_doseemmeans, sigma = sigma(AUCwindowmodel_GRABDA2h_ar1), edf = df.residual(AUCwindowmodel_GRABDA2h_ar1))) # Effect Size
AUCwindowmodel_GRABDA2h_ar1_dosepairs_es <- AUCwindowmodel_GRABDA2h_ar1_dosepairs_es[, c("InjType", "contrast", "effect.size")]

AUCwindowmodel_GRABDA2h_ar1_dosecontrasts$d <- AUCwindowmodel_GRABDA2h_ar1_dosepairs_es$effect.size
AUCwindowmodel_GRABDA2h_ar1_dosecontrasts

# Format Tables
table_AUCwindowmodel_GRABDA2h_ar1 <- format_regsummary_table(AUCwindowmodel_GRABDA2h_ar1, 'GRABDA2h', "AUCwindow", "AR1Model")
table_AUCwindowmodel_GRABDA2h_ar1_aov <- format_anova_table(AUCwindowmodel_GRABDA2h_ar1_aov, "GRABDA2h", "AUCwindow", "AR1Model") 
table_AUCwindowmodel_GRABDA2h_ar1_joint_Dose <- format_emm_table(AUCwindowmodel_GRABDA2h_ar1_jt_Dose, "GRABDA2h", "AUCwindow", "Joint", "Dose")
table_AUCwindowmodel_GRABDA2h_ar1_joint_InjType <- format_emm_table(AUCwindowmodel_GRABDA2h_ar1_jt_InjType, "GRABDA2h", "AUCwindow", "Joint", "InjType")
table_AUCwindowmodel_GRABDA2h_ar1_pairs_Dose <- format_emm_table(AUCwindowmodel_GRABDA2h_ar1_dosecontrasts, "GRABDA2h", "AUCwindow", "Pairs", "Dose|InjType")
table_AUCwindowmodel_GRABDA2h_ar1_pairs_InjType <- format_emm_table(AUCwindowmodel_GRABDA2h_ar1_treatcontrasts, "GRABDA2h", "AUCwindow", "Pairs", "InjType|Dose")

#### AUCwindow Change Model -----
# Mixed effects model
AUCwindowCmodel_GRABDA2h <- glmmTMB(AUCwindowC ~ Dose*BinF + (1|SubjectID), # random intercept by subject
                                      data = subset(subjectmeansbin, Virus=='GRABDA2h' & CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                      control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowCmodel_GRABDA2h_aov <- Anova(AUCwindowCmodel_GRABDA2h, type=2)
AUCwindowCmodel_GRABDA2h_aov
summary(AUCwindowCmodel_GRABDA2h)

AUCwindowCmodel_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(AUCwindowCmodel_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowCmodel_GRABDA2h_r2

# Follow up tests - joint
AUCwindowCmodel_GRABDA2h_emmeans <-emmeans(AUCwindowCmodel_GRABDA2h, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
AUCwindowCmodel_GRABDA2h_dosepairs <- emmeans(AUCwindowCmodel_GRABDA2h, pairwise~Dose, adjust="none") # Pairs
AUCwindowCmodel_GRABDA2h_dosecontrasts <- as.data.frame(AUCwindowCmodel_GRABDA2h_dosepairs$contrasts) # Extract emmeans
AUCwindowCmodel_GRABDA2h_doseemmeans <- AUCwindowCmodel_GRABDA2h_dosepairs$emmeans # Extract emmeans
AUCwindowCmodel_GRABDA2h_dosepairs_es <- as.data.frame(eff_size(AUCwindowCmodel_GRABDA2h_doseemmeans, sigma = sigma(AUCwindowCmodel_GRABDA2h), edf = df.residual(AUCwindowCmodel_GRABDA2h))) # Effect Size
AUCwindowCmodel_GRABDA2h_dosepairs_es <- AUCwindowCmodel_GRABDA2h_dosepairs_es[, c("contrast", "effect.size")]

AUCwindowCmodel_GRABDA2h_dosecontrasts$d <- AUCwindowCmodel_GRABDA2h_dosepairs_es$effect.size
AUCwindowCmodel_GRABDA2h_dosecontrasts

# Format Tables
table_AUCwindowCmodel_GRABDA2h <- format_regsummary_table(AUCwindowCmodel_GRABDA2h, 'GRABDA2h', "Amplitude", "StandardModel")
table_AUCwindowCmodel_GRABDA2h_aov <- format_anova_table(AUCwindowCmodel_GRABDA2h_aov, "GRABDA2h", "Amplitude", "StandardModel")
table_AUCwindowCmodel_GRABDA2h_pairs_Dose <- format_emm_table(AUCwindowCmodel_GRABDA2h_dosecontrasts, "GRABDA2h", "Amplitude", "Pairs","Dose")

#### Mixed effects AR1 model
AUCwindowCmodel_GRABDA2h_ar1 <- glmmTMB(AUCwindowC ~ Dose*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                          data = subset(subjectmeansbin, Virus=='GRABDA2h' & CondNum == 1& InjType=='M'), family = stats::gaussian, 
                                          control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowCmodel_GRABDA2h_ar1_aov <- Anova(AUCwindowCmodel_GRABDA2h_ar1, type=2)
AUCwindowCmodel_GRABDA2h_ar1_aov
summary(AUCwindowCmodel_GRABDA2h_ar1)

AUCwindowCmodel_GRABDA2h_ar1_r2 <- as.data.frame(r.squaredGLMM(AUCwindowCmodel_GRABDA2h_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowCmodel_GRABDA2h_ar1_r2

AUCwindowCmodel_GRABDA2h_ar1_AIC <- as.data.frame(AIC(AUCwindowCmodel_GRABDA2h, AUCwindowCmodel_GRABDA2h_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "AUCwindowCmodel_GRABDA2h" ~ "StandardModel", ModelName == "AUCwindowCmodel_GRABDA2h_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

AUCwindowCmodel_GRABDA2h_ar1_AIC

# Follow up tests - joint
AUCwindowCmodel_GRABDA2h_ar1_emmeans <-emmeans(AUCwindowCmodel_GRABDA2h_ar1, pairwise~Dose, adjust="none") # Find means

# Follow up tests - dose pairs
AUCwindowCmodel_GRABDA2h_ar1_dosepairs <- emmeans(AUCwindowCmodel_GRABDA2h_ar1, pairwise~Dose, adjust="none") # Pairs
AUCwindowCmodel_GRABDA2h_ar1_dosecontrasts <- as.data.frame(AUCwindowCmodel_GRABDA2h_ar1_dosepairs$contrasts) # Extract emmeans
AUCwindowCmodel_GRABDA2h_ar1_doseemmeans <- AUCwindowCmodel_GRABDA2h_ar1_dosepairs$emmeans # Extract emmeans
AUCwindowCmodel_GRABDA2h_ar1_dosepairs_es <- as.data.frame(eff_size(AUCwindowCmodel_GRABDA2h_ar1_doseemmeans, sigma = sigma(AUCwindowCmodel_GRABDA2h_ar1), edf = df.residual(AUCwindowCmodel_GRABDA2h_ar1))) # Effect Size
AUCwindowCmodel_GRABDA2h_ar1_dosepairs_es <- AUCwindowCmodel_GRABDA2h_ar1_dosepairs_es[, c("contrast", "effect.size")]

AUCwindowCmodel_GRABDA2h_ar1_dosecontrasts$d <- AUCwindowCmodel_GRABDA2h_ar1_dosepairs_es$effect.size
AUCwindowCmodel_GRABDA2h_ar1_dosecontrasts

# Format Tables
table_AUCwindowCmodel_GRABDA2h_ar1 <- format_regsummary_table(AUCwindowCmodel_GRABDA2h_ar1, 'GRABDA2h', "Amplitude", "AR1Model")
table_AUCwindowCmodel_GRABDA2h_ar1_aov <- format_anova_table(AUCwindowCmodel_GRABDA2h_ar1_aov, "GRABDA2h", "Amplitude", "AR1Model") 
table_AUCwindowCmodel_GRABDA2h_ar1_pairs_Dose <- format_emm_table(AUCwindowCmodel_GRABDA2h_ar1_dosecontrasts, "GRABDA2h", "Amplitude", "Pairs", "Dose")

#### Max Height Model -----
## Mixed effects model
MaxHeightmodel_AUCwindow_GRABDA2h <- glmmTMB(MaxHeight ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                               data = subjectrise_AUCwindowC_GRABDA2h, family = stats::gaussian)

MaxHeightmodel_AUCwindow_GRABDA2h_aov <- Anova(MaxHeightmodel_AUCwindow_GRABDA2h, type=2)
MaxHeightmodel_AUCwindow_GRABDA2h_aov
summary(MaxHeightmodel_AUCwindow_GRABDA2h)

MaxHeightmodel_AUCwindow_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(MaxHeightmodel_AUCwindow_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxHeightmodel_AUCwindow_GRABDA2h_r2

## Prepare fit data for plotting
MaxHeight_AUCwindow_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_GRABDA2h$DoseNum)))
MaxHeight_AUCwindow_GRABDA2h_regdata$DoseNum_c <- MaxHeight_AUCwindow_GRABDA2h_regdata$DoseNum - dosemean

MaxHeight_AUCwindow_GRABDA2h_regpred <- predict(MaxHeightmodel_AUCwindow_GRABDA2h, newdata = MaxHeight_AUCwindow_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxHeight_AUCwindow_GRABDA2h_regdata$fit <- MaxHeight_AUCwindow_GRABDA2h_regpred$fit
MaxHeight_AUCwindow_GRABDA2h_regdata$se <- MaxHeight_AUCwindow_GRABDA2h_regpred$se
MaxHeight_AUCwindow_GRABDA2h_regdata$lower <- MaxHeight_AUCwindow_GRABDA2h_regdata$fit - 1.96 * MaxHeight_AUCwindow_GRABDA2h_regdata$se
MaxHeight_AUCwindow_GRABDA2h_regdata$upper <- MaxHeight_AUCwindow_GRABDA2h_regdata$fit + 1.96 * MaxHeight_AUCwindow_GRABDA2h_regdata$se

MaxHeight_AUCwindow_GRABDA2h_regdata

# Format Tables
table_AUCwindowmodel_GRABDA2h_MaxHeight <- format_regsummary_table(MaxHeightmodel_AUCwindow_GRABDA2h, 'GRABDA2h', "AUCwindow", "MaxHeightModel")
table_AUCwindowmodel_GRABDA2h_MaxHeight_aov <- format_anova_table(MaxHeightmodel_AUCwindow_GRABDA2h_aov, "GRABDA2h", "AUCwindow", "MaxHeightModel") 
table_AUCwindowmodel_GRABDA2h_MaxHeight_reg <- format_reg_table(MaxHeight_AUCwindow_GRABDA2h_regdata, "GRABDA2h", "AUCwindow", "MaxHeightReg")

#### Max Duration Model -----
## Mixed effects model
MaxDurationmodel_AUCwindow_GRABDA2h <- glmmTMB(MaxDuration~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                                 data = subjectrise_AUCwindowC_GRABDA2h, family = stats::gaussian)

MaxDurationmodel_AUCwindow_GRABDA2h_aov<- Anova(MaxDurationmodel_AUCwindow_GRABDA2h, type=2)
MaxDurationmodel_AUCwindow_GRABDA2h_aov
summary(MaxDurationmodel_AUCwindow_GRABDA2h)

MaxDurationmodel_AUCwindow_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(MaxDurationmodel_AUCwindow_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxDurationmodel_AUCwindow_GRABDA2h_r2

## Prepare fit data for plotting
MaxDuration_AUCwindow_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_GRABDA2h$DoseNum)))
MaxDuration_AUCwindow_GRABDA2h_regdata$DoseNum_c <- MaxDuration_AUCwindow_GRABDA2h_regdata$DoseNum - dosemean

MaxDuration_AUCwindow_GRABDA2h_regpred <- predict(MaxDurationmodel_AUCwindow_GRABDA2h, newdata = MaxDuration_AUCwindow_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxDuration_AUCwindow_GRABDA2h_regdata$fit <- MaxDuration_AUCwindow_GRABDA2h_regpred$fit
MaxDuration_AUCwindow_GRABDA2h_regdata$se <- MaxDuration_AUCwindow_GRABDA2h_regpred$se
MaxDuration_AUCwindow_GRABDA2h_regdata$lower <- MaxDuration_AUCwindow_GRABDA2h_regdata$fit - 1.96 * MaxDuration_AUCwindow_GRABDA2h_regdata$se
MaxDuration_AUCwindow_GRABDA2h_regdata$upper <- MaxDuration_AUCwindow_GRABDA2h_regdata$fit + 1.96 * MaxDuration_AUCwindow_GRABDA2h_regdata$se

MaxDuration_AUCwindow_GRABDA2h_regdata

# Format Tables
table_AUCwindowmodel_GRABDA2h_MaxDuration <- format_regsummary_table(MaxDurationmodel_AUCwindow_GRABDA2h, 'GRABDA2h', "AUCwindow", "MaxDurationModel")
table_AUCwindowmodel_GRABDA2h_MaxDuration_aov <- format_anova_table(MaxDurationmodel_AUCwindow_GRABDA2h_aov, "GRABDA2h", "AUCwindow", "MaxDurationModel") 
table_AUCwindowmodel_GRABDA2h_MaxDuration_reg <- format_reg_table(MaxDuration_AUCwindow_GRABDA2h_regdata, "GRABDA2h", "AUCwindow", "MaxDurationReg")

#### Rise Slope Model -----
## Mixed effects model
MaxSlopeLMmodel_AUCwindow_GRABDA2h <- glmmTMB(MaxSlopeLM ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                                data = subjectrise_AUCwindowC_GRABDA2h, family = stats::gaussian)

MaxSlopeLMmodel_AUCwindow_GRABDA2h_aov <- Anova(MaxSlopeLMmodel_AUCwindow_GRABDA2h, type=2)
MaxSlopeLMmodel_AUCwindow_GRABDA2h_aov
summary(MaxSlopeLMmodel_AUCwindow_GRABDA2h)

MaxSlopeLMmodel_AUCwindow_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(MaxSlopeLMmodel_AUCwindow_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
MaxSlopeLMmodel_AUCwindow_GRABDA2h_r2

## Prepare fit data for plotting
MaxSlopeLM_AUCwindow_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_GRABDA2h$DoseNum)))
MaxSlopeLM_AUCwindow_GRABDA2h_regdata$DoseNum_c <- MaxSlopeLM_AUCwindow_GRABDA2h_regdata$DoseNum - dosemean

MaxSlopeLM_AUCwindow_GRABDA2h_regpred <- predict(MaxSlopeLMmodel_AUCwindow_GRABDA2h, newdata = MaxSlopeLM_AUCwindow_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
MaxSlopeLM_AUCwindow_GRABDA2h_regdata$fit <- MaxSlopeLM_AUCwindow_GRABDA2h_regpred$fit
MaxSlopeLM_AUCwindow_GRABDA2h_regdata$se <- MaxSlopeLM_AUCwindow_GRABDA2h_regpred$se
MaxSlopeLM_AUCwindow_GRABDA2h_regdata$lower <- MaxSlopeLM_AUCwindow_GRABDA2h_regdata$fit - 1.96 * MaxSlopeLM_AUCwindow_GRABDA2h_regdata$se
MaxSlopeLM_AUCwindow_GRABDA2h_regdata$upper <- MaxSlopeLM_AUCwindow_GRABDA2h_regdata$fit + 1.96 * MaxSlopeLM_AUCwindow_GRABDA2h_regdata$se

MaxSlopeLM_AUCwindow_GRABDA2h_regdata

# Format Tables
table_AUCwindowmodel_GRABDA2h_MaxSlopeLM <- format_regsummary_table(MaxSlopeLMmodel_AUCwindow_GRABDA2h, 'GRABDA2h', "AUCwindow", "MaxSlopeLMModel")
table_AUCwindowmodel_GRABDA2h_MaxSlopeLM_aov <- format_anova_table(MaxSlopeLMmodel_AUCwindow_GRABDA2h_aov, "GRABDA2h", "AUCwindow", "MaxSlopeLMModel") 
table_AUCwindowmodel_GRABDA2h_MaxSlopeLM_reg <- format_reg_table(MaxSlopeLM_AUCwindow_GRABDA2h_regdata, "GRABDA2h", "AUCwindow", "MaxSlopeLMReg")

#### AUC Model -----
## Mixed effects model
AUCmodel_AUCwindow_GRABDA2h <- glmmTMB(AUC_POSTINJ ~ DoseNum_c + (1|SubjectID), # random intercept by subject
                                         data = subset(subjectrise_AUCwindowC_GRABDA2h, Virus=='GRABDA2h'&InjType=='M'), family = stats::gaussian)

AUCmodel_AUCwindow_GRABDA2h_aov <- Anova(AUCmodel_AUCwindow_GRABDA2h, type=2)
AUCmodel_AUCwindow_GRABDA2h_aov
summary(AUCmodel_AUCwindow_GRABDA2h)

AUCmodel_AUCwindow_GRABDA2h_r2 <- as.data.frame(r.squaredGLMM(AUCmodel_AUCwindow_GRABDA2h)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCmodel_AUCwindow_GRABDA2h_r2

## Prepare fit data for plotting
AUC_AUCwindow_GRABDA2h_regdata <- data.frame(DoseNum = sort(unique(subjectrise_AUCwindowC_GRABDA2h$DoseNum)))
AUC_AUCwindow_GRABDA2h_regdata$DoseNum_c <- AUC_AUCwindow_GRABDA2h_regdata$DoseNum - dosemean

AUC_AUCwindow_GRABDA2h_regpred <- predict(AUCmodel_AUCwindow_GRABDA2h, newdata = AUC_AUCwindow_GRABDA2h_regdata, se.fit = TRUE, re.form = NA, type = "response")
AUC_AUCwindow_GRABDA2h_regdata$fit <- AUC_AUCwindow_GRABDA2h_regpred$fit
AUC_AUCwindow_GRABDA2h_regdata$se <- AUC_AUCwindow_GRABDA2h_regpred$se
AUC_AUCwindow_GRABDA2h_regdata$lower <- AUC_AUCwindow_GRABDA2h_regdata$fit - 1.96 * AUC_AUCwindow_GRABDA2h_regdata$se
AUC_AUCwindow_GRABDA2h_regdata$upper <- AUC_AUCwindow_GRABDA2h_regdata$fit + 1.96 * AUC_AUCwindow_GRABDA2h_regdata$se

AUC_AUCwindow_GRABDA2h_regdata

# Format Tables
table_AUCwindowmodel_GRABDA2h_AUC <- format_regsummary_table(AUCmodel_AUCwindow_GRABDA2h, 'GRABDA2h', "AUCwindow", "AUCModel")
table_AUCwindowmodel_GRABDA2h_AUC_aov <- format_anova_table(AUCmodel_AUCwindow_GRABDA2h_aov, "GRABDA2h", "AUCwindow", "AUCModel") 
table_AUCwindowmodel_GRABDA2h_AUC_reg <- format_reg_table(AUC_AUCwindow_GRABDA2h_regdata, "GRABDA2h", "AUCwindow", "AUCReg")

## EXPORT STATS -----
### Prepare Means -----
# Frequency
treatmeanstable_freq_GRABDA2h <- treatmeans %>% filter(Virus=='GRABDA2h') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         pkspermin_mean, pkspermin_sd, pkspermin_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansCtable_freq_GRABDA2h <- treatmeans %>% filter(InjType!='S', Virus=='GRABDA2h') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         pksperminC_mean, pksperminC_sd, pksperminC_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansSNRtable_freq_GRABDA2h <- treatmeans %>% filter(InjType!='S', Virus=='GRABDA2h') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         pksperminSNR_virus_mean, pksperminSNR_virus_sd, pksperminSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Amplitude
treatmeanstable_amp_GRABDA2h <- treatmeans %>% filter(Virus=='GRABDA2h') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         amp_mean, amp_sd, amp_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansCtable_amp_GRABDA2h <- treatmeans %>% filter(InjType!='S', Virus=='GRABDA2h') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         ampC_mean, ampC_sd, ampC_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansSNRtable_amp_GRABDA2h <- treatmeans %>% filter(InjType!='S', Virus=='GRABDA2h') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         ampSNR_virus_mean, ampSNR_virus_sd, ampSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# AUC Window
treatmeanstable_AUCwindow_GRABDA2h <- treatmeans %>% filter(Virus=='GRABDA2h') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         AUCwindow_mean, AUCwindow_sd, AUCwindow_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansCtable_AUCwindow_GRABDA2h <- treatmeans %>% filter(InjType!='S', Virus=='GRABDA2h') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         AUCwindowC_mean, AUCwindowC_sd, AUCwindowC_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

treatmeansSNRtable_AUCwindow_GRABDA2h <- treatmeans %>% filter(InjType!='S', Virus=='GRABDA2h') %>%
  select(FiberPlacement, Virus, InjType, Dose, DoseNum, TreatNum, CondNum, 
         AUCwindowSNR_virus_mean, AUCwindowSNR_virus_sd, AUCwindowSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))


### Frequency Export -----
wb_GRABDA2h_freq <- createWorkbook()
wb_GRABDA2h_freq_ar1 <- createWorkbook()
wb_GRABDA2h_freqC <- createWorkbook()
wb_GRABDA2h_freqC_ar1 <- createWorkbook()
wb_GRABDA2h_freq_reg <- createWorkbook()

## Main results for GRABDA2h
# Standard Models
write_nice_sheet(wb_GRABDA2h_freq, "MeansRaw", treatmeanstable_freq_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_freq, "ModelSummary", table_freqmodel_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_freq, "ModelAnova", table_freqmodel_GRABDA2h_aov)
write_nice_sheet(wb_GRABDA2h_freq, "ModelR2", freqmodel_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_freq, "Joint_InjType",  table_freqmodel_GRABDA2h_joint_InjType)
write_nice_sheet(wb_GRABDA2h_freq, "Joint_Dose", table_freqmodel_GRABDA2h_joint_Dose)
write_nice_sheet(wb_GRABDA2h_freq, "Pairs_InjType", table_freqmodel_GRABDA2h_pairs_InjType)
write_nice_sheet(wb_GRABDA2h_freq, "Pairs_Dose",  table_freqmodel_GRABDA2h_pairs_Dose)

saveWorkbook(wb_GRABDA2h_freq, paste(analysispath,"GRABDA2h_Frequency_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_GRABDA2h_freq_ar1, "MeansRaw", treatmeanstable_freq_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_freq_ar1, "ModelSummary", table_freqmodel_GRABDA2h_ar1)
write_nice_sheet(wb_GRABDA2h_freq_ar1, "ModelAnova", table_freqmodel_GRABDA2h_ar1_aov)
write_nice_sheet(wb_GRABDA2h_freq_ar1, "ModelR2", freqmodel_GRABDA2h_ar1_r2)
write_nice_sheet(wb_GRABDA2h_freq_ar1, "ModelAIC", freqmodel_GRABDA2h_ar1_AIC)
write_nice_sheet(wb_GRABDA2h_freq_ar1, "Joint_InjType",  table_freqmodel_GRABDA2h_ar1_joint_InjType)
write_nice_sheet(wb_GRABDA2h_freq_ar1, "Joint_Dose", table_freqmodel_GRABDA2h_ar1_joint_Dose)
write_nice_sheet(wb_GRABDA2h_freq_ar1, "Pairs_InjType", table_freqmodel_GRABDA2h_ar1_pairs_InjType)
write_nice_sheet(wb_GRABDA2h_freq_ar1, "Pairs_Dose",  table_freqmodel_GRABDA2h_ar1_pairs_Dose)

saveWorkbook(wb_GRABDA2h_freq_ar1, paste(analysispath,"GRABDA2h_Frequency_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Change Models
write_nice_sheet(wb_GRABDA2h_freqC, "MeansChange", treatmeanstable_freq_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_freqC, "ModelSummary", table_freqCmodel_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_freqC, "ModelAnova", table_freqCmodel_GRABDA2h_aov)
write_nice_sheet(wb_GRABDA2h_freqC, "ModelR2", freqCmodel_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_freqC, "Pairs_Dose",  table_freqCmodel_GRABDA2h_pairs_Dose)

saveWorkbook(wb_GRABDA2h_freqC, paste(analysispath,"GRABDA2h_FrequencyChange_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Change Models
write_nice_sheet(wb_GRABDA2h_freqC_ar1, "MeansChange", treatmeanstable_freq_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_freqC_ar1, "ModelSummary", table_freqCmodel_GRABDA2h_ar1)
write_nice_sheet(wb_GRABDA2h_freqC_ar1, "ModelAnova", table_freqCmodel_GRABDA2h_ar1_aov)
write_nice_sheet(wb_GRABDA2h_freqC_ar1, "ModelR2", freqCmodel_GRABDA2h_ar1_r2)
write_nice_sheet(wb_GRABDA2h_freqC_ar1, "ModelAIC", freqCmodel_GRABDA2h_ar1_AIC)
write_nice_sheet(wb_GRABDA2h_freqC_ar1, "Pairs_Dose",  table_freqCmodel_GRABDA2h_ar1_pairs_Dose)

saveWorkbook(wb_GRABDA2h_freqC_ar1, paste(analysispath,"GRABDA2h_FrequencyChange_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Regression Models
write_nice_sheet(wb_GRABDA2h_freq_reg, "MaxHeightSummary", table_freqmodel_GRABDA2h_MaxHeight)
write_nice_sheet(wb_GRABDA2h_freq_reg, "MaxHeightAnova", table_freqmodel_GRABDA2h_MaxHeight_aov)
write_nice_sheet(wb_GRABDA2h_freq_reg, "MaxHeightR2", MaxHeightmodel_freq_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_freq_reg, "MaxHeightReg",  table_freqmodel_GRABDA2h_MaxHeight_reg)
write_nice_sheet(wb_GRABDA2h_freq_reg, "MaxDurationSummary", table_freqmodel_GRABDA2h_MaxDuration)
write_nice_sheet(wb_GRABDA2h_freq_reg, "MaxDurationAnova", table_freqmodel_GRABDA2h_MaxDuration_aov)
write_nice_sheet(wb_GRABDA2h_freq_reg, "MaxDurationR2", MaxDurationmodel_freq_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_freq_reg, "MaxDurationReg", table_freqmodel_GRABDA2h_MaxDuration_reg)
write_nice_sheet(wb_GRABDA2h_freq_reg, "SlopeSummary", table_freqmodel_GRABDA2h_MaxSlopeLM)
write_nice_sheet(wb_GRABDA2h_freq_reg, "SlopeAnova", table_freqmodel_GRABDA2h_MaxSlopeLM_aov)
write_nice_sheet(wb_GRABDA2h_freq_reg, "SlopeR2", MaxSlopeLMmodel_freq_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_freq_reg, "SlopeReg", table_freqmodel_GRABDA2h_MaxSlopeLM_reg)
write_nice_sheet(wb_GRABDA2h_freq_reg, "AUCSummary", table_freqmodel_GRABDA2h_AUC)
write_nice_sheet(wb_GRABDA2h_freq_reg, "AUCAnova", table_freqmodel_GRABDA2h_AUC_aov)
write_nice_sheet(wb_GRABDA2h_freq_reg, "AUCR2", AUCmodel_freq_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_freq_reg, "AUCReg", table_freqmodel_GRABDA2h_AUC_reg)

saveWorkbook(wb_GRABDA2h_freq_reg, paste(analysispath,"GRABDA2h_Frequency_regressionmodels_results.xlsx",sep=''), overwrite = TRUE)

### Amplitude Export -----
wb_GRABDA2h_amp <- createWorkbook()
wb_GRABDA2h_amp_ar1 <- createWorkbook()
wb_GRABDA2h_ampC <- createWorkbook()
wb_GRABDA2h_ampC_ar1 <- createWorkbook()
wb_GRABDA2h_amp_reg <- createWorkbook()

## Main results for GRABDA2h
# Standard Models
write_nice_sheet(wb_GRABDA2h_amp, "MeansRaw", treatmeansCtable_amp_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_amp, "ModelSummary", table_ampmodel_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_amp, "ModelAnova", table_ampmodel_GRABDA2h_aov)
write_nice_sheet(wb_GRABDA2h_amp, "ModelR2", ampmodel_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_amp, "Joint_InjType",  table_ampmodel_GRABDA2h_joint_InjType)
write_nice_sheet(wb_GRABDA2h_amp, "Joint_Dose", table_ampmodel_GRABDA2h_joint_Dose)
write_nice_sheet(wb_GRABDA2h_amp, "Pairs_InjType", table_ampmodel_GRABDA2h_pairs_InjType)
write_nice_sheet(wb_GRABDA2h_amp, "Pairs_Dose",  table_ampmodel_GRABDA2h_pairs_Dose)

saveWorkbook(wb_GRABDA2h_amp, paste(analysispath,"GRABDA2h_Amplitude_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_GRABDA2h_amp_ar1, "MeansRaw", treatmeanstable_amp_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_amp_ar1, "ModelSummary", table_ampmodel_GRABDA2h_ar1)
write_nice_sheet(wb_GRABDA2h_amp_ar1, "ModelAnova", table_ampmodel_GRABDA2h_ar1_aov)
write_nice_sheet(wb_GRABDA2h_amp_ar1, "ModelR2", ampmodel_GRABDA2h_ar1_r2)
write_nice_sheet(wb_GRABDA2h_amp_ar1, "ModelAIC", ampmodel_GRABDA2h_ar1_AIC)
write_nice_sheet(wb_GRABDA2h_amp_ar1, "Joint_InjType",  table_ampmodel_GRABDA2h_ar1_joint_InjType)
write_nice_sheet(wb_GRABDA2h_amp_ar1, "Joint_Dose", table_ampmodel_GRABDA2h_ar1_joint_Dose)
write_nice_sheet(wb_GRABDA2h_amp_ar1, "Pairs_InjType", table_ampmodel_GRABDA2h_ar1_pairs_InjType)
write_nice_sheet(wb_GRABDA2h_amp_ar1, "Pairs_Dose",  table_ampmodel_GRABDA2h_ar1_pairs_Dose)

saveWorkbook(wb_GRABDA2h_amp_ar1, paste(analysispath,"GRABDA2h_Amplitude_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Change Models
write_nice_sheet(wb_GRABDA2h_ampC, "MeansChange", treatmeanstable_amp_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_ampC, "ModelSummary", table_ampCmodel_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_ampC, "ModelAnova", table_ampCmodel_GRABDA2h_aov)
write_nice_sheet(wb_GRABDA2h_ampC, "ModelR2", ampCmodel_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_ampC, "Pairs_Dose",  table_ampCmodel_GRABDA2h_pairs_Dose)

saveWorkbook(wb_GRABDA2h_ampC, paste(analysispath,"GRABDA2h_AmplitudeChange_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Change Models
write_nice_sheet(wb_GRABDA2h_ampC_ar1, "MeansChange", treatmeanstable_amp_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_ampC_ar1, "ModelSummary", table_ampCmodel_GRABDA2h_ar1)
write_nice_sheet(wb_GRABDA2h_ampC_ar1, "ModelAnova", table_ampCmodel_GRABDA2h_ar1_aov)
write_nice_sheet(wb_GRABDA2h_ampC_ar1, "ModelR2", ampCmodel_GRABDA2h_ar1_r2)
write_nice_sheet(wb_GRABDA2h_ampC_ar1, "ModelAIC", ampCmodel_GRABDA2h_ar1_AIC)
write_nice_sheet(wb_GRABDA2h_ampC_ar1, "Pairs_Dose",  table_ampCmodel_GRABDA2h_ar1_pairs_Dose)

saveWorkbook(wb_GRABDA2h_ampC_ar1, paste(analysispath,"GRABDA2h_AmplitudeChange_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Regression Models
write_nice_sheet(wb_GRABDA2h_amp_reg, "MaxHeightSummary", table_ampmodel_GRABDA2h_MaxHeight)
write_nice_sheet(wb_GRABDA2h_amp_reg, "MaxHeightAnova", table_ampmodel_GRABDA2h_MaxHeight_aov)
write_nice_sheet(wb_GRABDA2h_amp_reg, "MaxHeightR2", MaxHeightmodel_amp_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_amp_reg, "MaxHeightReg",  table_ampmodel_GRABDA2h_MaxHeight_reg)
write_nice_sheet(wb_GRABDA2h_amp_reg, "MaxDurationSummary", table_ampmodel_GRABDA2h_MaxDuration)
write_nice_sheet(wb_GRABDA2h_amp_reg, "MaxDurationAnova", table_ampmodel_GRABDA2h_MaxDuration_aov)
write_nice_sheet(wb_GRABDA2h_amp_reg, "MaxDurationR2", MaxDurationmodel_amp_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_amp_reg, "MaxDurationReg", table_ampmodel_GRABDA2h_MaxDuration_reg)
write_nice_sheet(wb_GRABDA2h_amp_reg, "SlopeSummary", table_ampmodel_GRABDA2h_MaxSlopeLM)
write_nice_sheet(wb_GRABDA2h_amp_reg, "SlopeAnova", table_ampmodel_GRABDA2h_MaxSlopeLM_aov)
write_nice_sheet(wb_GRABDA2h_amp_reg, "SlopeR2", MaxSlopeLMmodel_amp_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_amp_reg, "SlopeReg", table_ampmodel_GRABDA2h_MaxSlopeLM_reg)
write_nice_sheet(wb_GRABDA2h_amp_reg, "AUCSummary", table_ampmodel_GRABDA2h_AUC)
write_nice_sheet(wb_GRABDA2h_amp_reg, "AUCAnova", table_ampmodel_GRABDA2h_AUC_aov)
write_nice_sheet(wb_GRABDA2h_amp_reg, "AUCR2", AUCmodel_amp_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_amp_reg, "AUCReg", table_ampmodel_GRABDA2h_AUC_reg)

saveWorkbook(wb_GRABDA2h_amp_reg, paste(analysispath,"GRABDA2h_Amplitude_regressionmodels_results.xlsx",sep=''), overwrite = TRUE)


### AUCwindow Export -----
wb_GRABDA2h_AUCwindow <- createWorkbook()
wb_GRABDA2h_AUCwindow_ar1 <- createWorkbook()
wb_GRABDA2h_AUCwindowC <- createWorkbook()
wb_GRABDA2h_AUCwindowC_ar1 <- createWorkbook()
wb_GRABDA2h_AUCwindow_reg <- createWorkbook()

## Main results for GRABDA2h
# Standard Models
write_nice_sheet(wb_GRABDA2h_AUCwindow, "MeansRaw", treatmeanstable_AUCwindow_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_AUCwindow, "ModelSummary", table_AUCwindowmodel_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_AUCwindow, "ModelAnova", table_AUCwindowmodel_GRABDA2h_aov)
write_nice_sheet(wb_GRABDA2h_AUCwindow, "ModelR2", AUCwindowmodel_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_AUCwindow, "Joint_InjType",  table_AUCwindowmodel_GRABDA2h_joint_InjType)
write_nice_sheet(wb_GRABDA2h_AUCwindow, "Joint_Dose", table_AUCwindowmodel_GRABDA2h_joint_Dose)
write_nice_sheet(wb_GRABDA2h_AUCwindow, "Pairs_InjType", table_AUCwindowmodel_GRABDA2h_pairs_InjType)
write_nice_sheet(wb_GRABDA2h_AUCwindow, "Pairs_Dose",  table_AUCwindowmodel_GRABDA2h_pairs_Dose)

saveWorkbook(wb_GRABDA2h_AUCwindow, paste(analysispath,"GRABDA2h_AUCwindow_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Models
write_nice_sheet(wb_GRABDA2h_AUCwindow_ar1, "MeansRaw", treatmeanstable_AUCwindow_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_AUCwindow_ar1, "ModelSummary", table_AUCwindowmodel_GRABDA2h_ar1)
write_nice_sheet(wb_GRABDA2h_AUCwindow_ar1, "ModelAnova", table_AUCwindowmodel_GRABDA2h_ar1_aov)
write_nice_sheet(wb_GRABDA2h_AUCwindow_ar1, "ModelR2", AUCwindowmodel_GRABDA2h_ar1_r2)
write_nice_sheet(wb_GRABDA2h_AUCwindow_ar1, "ModelAIC", AUCwindowmodel_GRABDA2h_ar1_AIC)
write_nice_sheet(wb_GRABDA2h_AUCwindow_ar1, "Joint_InjType",  table_AUCwindowmodel_GRABDA2h_ar1_joint_InjType)
write_nice_sheet(wb_GRABDA2h_AUCwindow_ar1, "Joint_Dose", table_AUCwindowmodel_GRABDA2h_ar1_joint_Dose)
write_nice_sheet(wb_GRABDA2h_AUCwindow_ar1, "Pairs_InjType", table_AUCwindowmodel_GRABDA2h_ar1_pairs_InjType)
write_nice_sheet(wb_GRABDA2h_AUCwindow_ar1, "Pairs_Dose",  table_AUCwindowmodel_GRABDA2h_ar1_pairs_Dose)

saveWorkbook(wb_GRABDA2h_AUCwindow_ar1, paste(analysispath,"GRABDA2h_AUCwindow_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Change Models
write_nice_sheet(wb_GRABDA2h_AUCwindowC, "MeansChange", treatmeanstable_AUCwindow_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_AUCwindowC, "ModelSummary", table_AUCwindowCmodel_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_AUCwindowC, "ModelAnova", table_AUCwindowCmodel_GRABDA2h_aov)
write_nice_sheet(wb_GRABDA2h_AUCwindowC, "ModelR2", AUCwindowCmodel_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_AUCwindowC, "Pairs_Dose",  table_AUCwindowCmodel_GRABDA2h_pairs_Dose)

saveWorkbook(wb_GRABDA2h_AUCwindowC, paste(analysispath,"GRABDA2h_AUCwindowChange_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# ar1 Change Models
write_nice_sheet(wb_GRABDA2h_AUCwindowC_ar1, "MeansChange", treatmeanstable_AUCwindow_GRABDA2h)
write_nice_sheet(wb_GRABDA2h_AUCwindowC_ar1, "ModelSummary", table_AUCwindowCmodel_GRABDA2h_ar1)
write_nice_sheet(wb_GRABDA2h_AUCwindowC_ar1, "ModelAnova", table_AUCwindowCmodel_GRABDA2h_ar1_aov)
write_nice_sheet(wb_GRABDA2h_AUCwindowC_ar1, "ModelR2", AUCwindowCmodel_GRABDA2h_ar1_r2)
write_nice_sheet(wb_GRABDA2h_AUCwindowC_ar1, "ModelAIC", AUCwindowCmodel_GRABDA2h_ar1_AIC)
write_nice_sheet(wb_GRABDA2h_AUCwindowC_ar1, "Pairs_Dose",  table_AUCwindowCmodel_GRABDA2h_ar1_pairs_Dose)

saveWorkbook(wb_GRABDA2h_AUCwindowC_ar1, paste(analysispath,"GRABDA2h_AUCwindowChange_ar1models_results.xlsx",sep=''), overwrite = TRUE)

# Regression Models
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "MaxHeightSummary", table_AUCwindowmodel_GRABDA2h_MaxHeight)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "MaxHeightAnova", table_AUCwindowmodel_GRABDA2h_MaxHeight_aov)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "MaxHeightR2", MaxHeightmodel_AUCwindow_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "MaxHeightReg",  table_AUCwindowmodel_GRABDA2h_MaxHeight_reg)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "MaxDurationSummary", table_AUCwindowmodel_GRABDA2h_MaxDuration)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "MaxDurationAnova", table_AUCwindowmodel_GRABDA2h_MaxDuration_aov)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "MaxDurationR2", MaxDurationmodel_AUCwindow_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "MaxDurationReg", table_AUCwindowmodel_GRABDA2h_MaxDuration_reg)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "SlopeSummary", table_AUCwindowmodel_GRABDA2h_MaxSlopeLM)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "SlopeAnova", table_AUCwindowmodel_GRABDA2h_MaxSlopeLM_aov)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "SlopeR2", MaxSlopeLMmodel_AUCwindow_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "SlopeReg", table_AUCwindowmodel_GRABDA2h_MaxSlopeLM_reg)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "AUCSummary", table_AUCwindowmodel_GRABDA2h_AUC)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "AUCAnova", table_AUCwindowmodel_GRABDA2h_AUC_aov)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "AUCR2", AUCmodel_AUCwindow_GRABDA2h_r2)
write_nice_sheet(wb_GRABDA2h_AUCwindow_reg, "AUCReg", table_AUCwindowmodel_GRABDA2h_AUC_reg)

saveWorkbook(wb_GRABDA2h_AUCwindow_reg, paste(analysispath,"GRABDA2h_AUCwindow_regressionmodels_results.xlsx",sep=''), overwrite = TRUE)



## PLOTS -----
### Frequency -----
freqbarmin_GRABDA2h <- 0
freqbarmax_GRABDA2h <- 2.65
freqbarticks_GRABDA2h <- .6
freqbarbreaks_GRABDA2h <- seq(freqbarmin_GRABDA2h,freqbarmax_GRABDA2h,freqbarticks_GRABDA2h)

# Bar Graph
freqbar_GRABDA2h <- ggplot(subset(treatmeans, Virus == "GRABDA2h"&CondNum==1),aes(x=Dose, y=pkspermin_mean, fill=TreatNum, group=InjType)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_jitter(data=subset(subjectmeans,Virus == "GRABDA2h" & CondNum == 1), aes(x=Dose, y=pkspermin, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GRABDA2h"&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=pkspermin_mean-pkspermin_se, ymax=pkspermin_mean+pkspermin_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(freqbarmin_GRABDA2h, freqbarmax_GRABDA2h), 
                     breaks=freqbarbreaks_GRABDA2h) +
  xlab("Dose (mg/kg)") +
  ylab("Frequency (n/min)") + 
  mytheme +
  # STATS FOR AR1 MIXED MODEL
  geom_bracket(xmin=c(0.765), xmax = c(1.235), y.position=c(1.6),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1.765), xmax = c(2.235), y.position=c(2.05),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.765), xmax = c(3.235), y.position=c(1.6),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.765), xmax = c(4.235), y.position=c(2.05),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

freqbar_GRABDA2h

# Bar graph (Change from Control)
freqCbarmin_GRABDA2h <- -2
freqCbarmax_GRABDA2h <- 0
freqCbarticks_GRABDA2h <- .5
freqCbarbreaks_GRABDA2h <- seq(freqCbarmin_GRABDA2h,freqCbarmax_GRABDA2h,freqCbarticks_GRABDA2h)

freqCbar_GRABDA2h <- ggplot(subset(treatmeans, Virus == "GRABDA2h"&CondNum==1&InjType=='M'),aes(x=Dose, y=pksperminC_mean, fill=TreatNum, group=InjType)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(subjectmeans,Virus == "GRABDA2h" & CondNum == 1&InjType=='M'), aes(x=Dose, y=pksperminC, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GRABDA2h"&CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=pksperminC_mean-pksperminC_se, ymax=pksperminC_mean+pksperminC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(freqCbarmin_GRABDA2h, freqCbarmax_GRABDA2h), breaks=freqCbarbreaks_GRABDA2h)+
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" Frequency (n/min)")) + 
  mytheme

freqCbar_GRABDA2h

# Line Graph (Change)
freqlineCmin_GRABDA2h <- -1.5
freqlineCmax_GRABDA2h <- 1.5
freqlineCticks_GRABDA2h <- .75
freqlineCbreaks_GRABDA2h <- seq(freqlineCmin_GRABDA2h, freqlineCmax_GRABDA2h, freqlineCticks_GRABDA2h)

freqline_GRABDA2h <- ggplot(subset(treatmeansbin, Virus=='GRABDA2h'&InjType!='S'), aes(x=BinTimeInj, y=pksperminC_mean, group=Dose, fill=Dose, color=Dose)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  scale_x_continuous(name="Minutes from Injection",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(freqlineCmin_GRABDA2h, freqlineCmax_GRABDA2h), 
                     breaks = freqlineCbreaks_GRABDA2h)+
  scale_color_manual("legend", values = treatcolors) +
  mytheme +
  ylab(expression(Delta*" Frequency (n/min)")) + 
  mytheme

freqline_GRABDA2h

#### Temporal Dynamics Panels -----
# Max Height
MaxHeight_freq_ymin_GRABDA2h <- -.2
MaxHeight_freq_ymax_GRABDA2h <- 2.8
MaxHeight_freq_ybreaks_GRABDA2h <- .8

MaxHeight_freq_GRABDA2h <- ggplot(subset(subjectrise_pksperminC, Virus=='GRABDA2h'&InjType=='M'), aes(x=DoseNum, y=MaxHeight, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxHeight_freq_GRABDA2h_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxHeight_freq_GRABDA2h_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxHeight_freq_ymin_GRABDA2h, MaxHeight_freq_ymax_GRABDA2h), 
                     breaks = seq(0, MaxHeight_freq_ymax_GRABDA2h, MaxHeight_freq_ybreaks_GRABDA2h)) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Max "*Delta*" n"))

MaxHeight_freq_GRABDA2h

# Max Duration
MaxDuration_freq_ymin_GRABDA2h <- 0
MaxDuration_freq_ymax_GRABDA2h <- 62
MaxDuration_freq_ybreaks_GRABDA2h <- 20

MaxDuration_freq_GRABDA2h <- ggplot(subset(subjectrise_pksperminC, Virus=='GRABDA2h'&InjType=='M'), aes(x=DoseNum, y=MaxDuration, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxDuration_freq_GRABDA2h_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxDuration_freq_GRABDA2h_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxDuration_freq_ymin_GRABDA2h, MaxDuration_freq_ymax_GRABDA2h), 
                     breaks = seq(MaxDuration_freq_ymin_GRABDA2h, MaxDuration_freq_ymax_GRABDA2h, MaxDuration_freq_ybreaks_GRABDA2h)) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Time to Max "*Delta*" n"))

MaxDuration_freq_GRABDA2h

# Max Height Slope
MaxSlopeLM_freq_ymin_GRABDA2h <- -.1
MaxSlopeLM_freq_ymax_GRABDA2h <- .2
MaxSlopeLM_freq_ybreaks_GRABDA2h <- .1

MaxSlopeLM_freq_GRABDA2h <- ggplot(subset(subjectrise_pksperminC, Virus=='GRABDA2h'&InjType=='M'), aes(x=DoseNum, y=MaxSlopeLM, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxSlopeLM_freq_GRABDA2h_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxSlopeLM_freq_GRABDA2h_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxSlopeLM_freq_ymin_GRABDA2h, MaxSlopeLM_freq_ymax_GRABDA2h), 
                     breaks = seq(MaxSlopeLM_freq_ymin_GRABDA2h, MaxSlopeLM_freq_ymax_GRABDA2h, MaxSlopeLM_freq_ybreaks_GRABDA2h))+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab('Time-Course Slope')

MaxSlopeLM_freq_GRABDA2h

# AUC POST INJ
AUC_freq_ymin_GRABDA2h <- -4
AUC_freq_ymax_GRABDA2h <- 45
AUC_freq_ybreaks_GRABDA2h <- 15

AUC_freq_GRABDA2h <- ggplot(subset(subjectrise_pksperminC, Virus=='GRABDA2h'&InjType=='M'), aes(x=DoseNum, y=AUC_POSTINJ, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = AUC_freq_GRABDA2h_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = AUC_freq_GRABDA2h_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(AUC_freq_ymin_GRABDA2h, AUC_freq_ymax_GRABDA2h), 
                     breaks = seq(0, AUC_freq_ymax_GRABDA2h, AUC_freq_ybreaks_GRABDA2h))+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab("Time-Course AUC")

AUC_freq_GRABDA2h


### Amplitude -----
ampbarmin_GRABDA2h <- 0
ampbarmax_GRABDA2h <- 1.2
ampbarticks_GRABDA2h <- .3
ampbarybreaks_GRABDA2h <- seq(ampbarmin_GRABDA2h,ampbarmax_GRABDA2h,ampbarticks_GRABDA2h)
ampbarylabels_GRABDA2h <- ampbarybreaks_GRABDA2h + 2.6

# Bar graph
ampbar_GRABDA2h <- ggplot(subset(treatmeans, Virus == "GRABDA2h"&CondNum==1),aes(x=Dose, y=ampscaled_mean, fill=TreatNum, group=InjType)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_jitter(data=subset(subjectmeans,Virus == "GRABDA2h" & CondNum == 1), aes(x=Dose, y=ampscaled, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GRABDA2h"&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=ampscaled_mean-amp_se, ymax=ampscaled_mean+amp_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(ampbarmin_GRABDA2h, ampbarmax_GRABDA2h), 
                     breaks=ampbarybreaks_GRABDA2h, labels=ampbarylabels_GRABDA2h)+
  xlab("Dose (mg/kg)") +
  ylab("Amplitude (z)") + 
  mytheme +
  # AR1 and ALL AMP STATS
  geom_bracket(xmin=c(0.765), xmax = c(1.235), y.position=c(.85),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1.765), xmax = c(2.235), y.position=c(.9),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.765), xmax = c(3.235), y.position=c(.9),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.765), xmax = c(4.235), y.position=c(.75),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

ampbar_GRABDA2h

# Bar graph (Change from Control)
ampCbarmin_GRABDA2h <- -.9
ampCbarmax_GRABDA2h <- .3
ampCbarticks_GRABDA2h <- .3
ampCbarbreaks_GRABDA2h <- round(seq(ampCbarmin_GRABDA2h,ampCbarmax_GRABDA2h,ampCbarticks_GRABDA2h),1)

ampCbar_GRABDA2h <- ggplot(subset(treatmeans, Virus == "GRABDA2h"&CondNum==1&InjType=='M'),aes(x=Dose, y=ampC_mean, fill=TreatNum, group=InjType)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(subjectmeans,Virus == "GRABDA2h" & CondNum == 1&InjType=='M'), aes(x=Dose, y=ampC, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GRABDA2h"&CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=ampC_mean-ampC_se, ymax=ampC_mean+ampC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(ampCbarmin_GRABDA2h, ampCbarmax_GRABDA2h), breaks=ampCbarbreaks_GRABDA2h)+
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" Amplitude (z)")) + 
  mytheme

ampCbar_GRABDA2h

# Line Graph (Change from S)
amplineCmin_GRABDA2h <- -1.2
amplineCmax_GRABDA2h <- 0.4
amplineCticks_GRABDA2h <- 0.4

ampline_GRABDA2h <- ggplot(subset(treatmeansbin, Virus=='GRABDA2h'&InjType!='S'), aes(x=BinTimeInj, y=ampC_mean, group=Dose, fill=Dose, color=Dose)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  scale_x_continuous(name="Minutes from Injection",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(amplineCmin_GRABDA2h, amplineCmax_GRABDA2h), 
                     breaks = seq(amplineCmin_GRABDA2h, amplineCmax_GRABDA2h, amplineCticks_GRABDA2h), labels = c('-1.2','-0.8','-0.4','0','0.4'))+
  scale_color_manual("legend", values = treatcolors) +
  scale_fill_manual("legend", values = treatcolors) +
  mytheme +
  ylab(expression(Delta*" Amplitude (z)")) + 
  mytheme

ampline_GRABDA2h

#### Temporal Dynamics Panels -----
# Max Height
MaxHeight_amp_ymin_GRABDA2h <- -.25
MaxHeight_amp_ymax_GRABDA2h <- 1.5
MaxHeight_amp_ybreaks_GRABDA2h <- .5

MaxHeight_amp_GRABDA2h <- ggplot(subset(subjectrise_ampC, Virus=='GRABDA2h'&InjType=='M'), aes(x=DoseNum, y=MaxHeight, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxHeight_amp_GRABDA2h_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxHeight_amp_GRABDA2h_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxHeight_amp_ymin_GRABDA2h, MaxHeight_amp_ymax_GRABDA2h), 
                     breaks = seq(0, MaxHeight_amp_ymax_GRABDA2h, MaxHeight_amp_ybreaks_GRABDA2h)) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Max "*Delta*" z"))

MaxHeight_amp_GRABDA2h

# Max Duration
MaxDuration_amp_ymin_GRABDA2h <- 0
MaxDuration_amp_ymax_GRABDA2h <- 62
MaxDuration_amp_ybreaks_GRABDA2h <- 20

MaxDuration_amp_GRABDA2h <- ggplot(subset(subjectrise_ampC, Virus=='GRABDA2h'&InjType=='M'), aes(x=DoseNum, y=MaxDuration, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxDuration_amp_GRABDA2h_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxDuration_amp_GRABDA2h_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxDuration_amp_ymin_GRABDA2h, MaxDuration_amp_ymax_GRABDA2h), 
                     breaks = seq(MaxDuration_amp_ymin_GRABDA2h, MaxDuration_amp_ymax_GRABDA2h, MaxDuration_amp_ybreaks_GRABDA2h)) + 
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab(expression("Time to Max "*Delta*" z"))

MaxDuration_amp_GRABDA2h

# Max Height Slope
MaxSlopeLM_amp_ymin_GRABDA2h <- -.12
MaxSlopeLM_amp_ymax_GRABDA2h <- .28
MaxSlopeLM_amp_ybreaks_GRABDA2h <- .12

MaxSlopeLM_amp_GRABDA2h <- ggplot(subset(subjectrise_ampC, Virus=='GRABDA2h'&InjType=='M'), aes(x=DoseNum, y=MaxSlopeLM, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = MaxSlopeLM_amp_GRABDA2h_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = MaxSlopeLM_amp_GRABDA2h_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(MaxSlopeLM_amp_ymin_GRABDA2h, MaxSlopeLM_amp_ymax_GRABDA2h), 
                     breaks = seq(MaxSlopeLM_amp_ymin_GRABDA2h, MaxSlopeLM_amp_ymax_GRABDA2h, MaxSlopeLM_amp_ybreaks_GRABDA2h))+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab('Time-Course Slope')

MaxSlopeLM_amp_GRABDA2h

# AUC POST INJ
AUC_amp_ymin_GRABDA2h <- -5
AUC_amp_ymax_GRABDA2h <- 36
AUC_amp_ybreaks_GRABDA2h <- 12

AUC_amp_GRABDA2h <- ggplot(subset(subjectrise_ampC, Virus=='GRABDA2h'&InjType=='M'), aes(x=DoseNum, y=AUC_POSTINJ, Group=SubjectID, color=factor(DoseNum), shape=Sex)) +
  geom_ribbon(data = AUC_amp_GRABDA2h_regdata, aes(x = DoseNum, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, fill = "grey70") +
  geom_jitter(size=mysubjectpointsize, position = position_jitterdodge(jitter.width = 1, jitter.height = 0), show.legend=FALSE) +
  geom_line(data = AUC_amp_GRABDA2h_regdata,  aes(x = DoseNum, y = fit, group = 1), linewidth = 1, inherit.aes = FALSE)+
  scale_x_continuous(name='Dose',limits = c(2, 10.7), breaks = c(2.5, 5, 7.5, 10), labels = c(2.5, 5, 7.5, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(AUC_amp_ymin_GRABDA2h, AUC_amp_ymax_GRABDA2h), 
                     breaks = seq(0, AUC_amp_ymax_GRABDA2h, AUC_amp_ybreaks_GRABDA2h))+
  scale_color_manual("legend", values = treatcolors) +
  scale_shape_manual(values = sexshapes) +
  mytheme +
  ylab("Time-Course AUC") + 
  annotate("text", x=10.6, y=0, label= "*", size = mystarsize*1.5, color = 'black')

AUC_amp_GRABDA2h

### AUC Window -----
AUCwindowbarmin <- 0
AUCwindowbarmax <- 2.4
AUCwindowbarticks <- .6
AUCwindowbarbreaks <- seq(AUCwindowbarmin,AUCwindowbarmax,AUCwindowbarticks)

# Bar graph
AUCwindowbar_GRABDA2h <- ggplot(subset(treatmeans, Virus == "GRABDA2h"&CondNum==1),aes(x=Dose, y=AUCwindow_mean, fill=TreatNum, group=InjType)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_jitter(data=subset(subjectmeans,Virus == "GRABDA2h" & CondNum == 1), aes(x=Dose, y=AUCwindow, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GRABDA2h"&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=AUCwindow_mean-AUCwindow_se, ymax=AUCwindow_mean+AUCwindow_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(AUCwindowbarmin, AUCwindowbarmax), breaks=AUCwindowbarbreaks)+
  xlab("Dose (mg/kg)") +
  ylab("AUC per Transient") + 
  mytheme +
  geom_bracket(xmin=c(0.765), xmax = c(1.235), y.position=c(2.32),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

AUCwindowbar_GRABDA2h

# Bar graph (Change from Control)
AUCwindowCbarmin_GRABDA2h <- -.6
AUCwindowCbarmax_GRABDA2h <- .3
AUCwindowCbarticks_GRABDA2h <- .3
AUCwindowCbarbreaks_GRABDA2h <- seq(AUCwindowCbarmin_GRABDA2h,AUCwindowCbarmax_GRABDA2h,AUCwindowCbarticks_GRABDA2h)

AUCwindowCbar_GRABDA2h <- ggplot(subset(treatmeans, Virus == "GRABDA2h"&CondNum==1&InjType=='M'),aes(x=Dose, y=AUCwindowC_mean, fill=TreatNum, group=InjType)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(subjectmeans,Virus == "GRABDA2h" & CondNum == 1&InjType=='M'), aes(x=Dose, y=AUCwindowC, fill=InjType, group=InjType, shape=Sex), color='#4c4c4c', size=mysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  geom_errorbar(data=subset(treatmeans, Virus == "GRABDA2h"&CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = InjType, ymin=AUCwindowC_mean-AUCwindowC_se, ymax=AUCwindowC_mean+AUCwindowC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = treatcolors_raw) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(AUCwindowCbarmin_GRABDA2h, AUCwindowCbarmax_GRABDA2h), breaks=AUCwindowCbarbreaks_GRABDA2h)+
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" AUC per Transient")) + 
  mytheme

AUCwindowCbar_GRABDA2h

# Line Graph (Change from S)
AUCwindowlineCmin <- -.8
AUCwindowlineCmax <- .8
AUCwindowlineCticks <- .4
AUCwindowlineCbreaks <- seq(AUCwindowlineCmin, AUCwindowlineCmax, AUCwindowlineCticks)

AUCwindowline_GRABDA2h <- ggplot(subset(treatmeansbin, Virus=='GRABDA2h'&InjType!='S'), aes(x=BinTimeInj, y=AUCwindowC_mean, group=Dose, fill=Dose, color=Dose)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mypointsize, show.legend = FALSE) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  scale_x_continuous(name="Minutes from Injection",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(AUCwindowlineCmin, AUCwindowlineCmax), breaks = AUCwindowlineCbreaks)+
  scale_color_manual("legend", values = treatcolors) +
  scale_fill_manual("legend", values = treatcolors) +
  mytheme +
  ylab(expression(Delta*" AUC per Transient")) + 
  mytheme

AUCwindowline_GRABDA2h

### EXPORT FIGURE 4 -----
# BAR AND LINE GRAPHS
figure4panels_barline <- plot_grid(freqbar_GRABDA2h, freqCbar_GRABDA2h, freqline_GRABDA2h, 
                                   ampbar_GRABDA2h, ampCbar_GRABDA2h, ampline_GRABDA2h,
                                   AUCwindowbar_GRABDA2h, AUCwindowCbar_GRABDA2h, AUCwindowline_GRABDA2h,
                                   ncol = 3, nrow=3, align = "hv",  # align horizontally & vertically
                                   axis  = "tblr", hjust=-.5, rel_widths = c(1,.8,1.2,
                                                                             1,.8,1.2,
                                                                             1,.8,1.2))
figure4panels_barline

ggsave(filename = c("Figure4_NAcLSGRABDA2h_Bar_TransientPanels_BarLineGraphs.pdf"), path = figurepath, plot = figure4panels_barline, device = "pdf", 
       width = figpanelwidthin, height = figpanelheightinbarline, dpi=figdpi, units="in")


# TEMPORAL DYNAMICS GRAPHS
figure4panels_models <- plot_grid(MaxHeight_freq_GRABDA2h, MaxDuration_freq_GRABDA2h,
                                  MaxSlopeLM_freq_GRABDA2h, AUC_freq_GRABDA2h,
                                  
                                  MaxHeight_amp_GRABDA2h, MaxDuration_amp_GRABDA2h,
                                  MaxSlopeLM_amp_GRABDA2h, AUC_amp_GRABDA2h,
                                  
                                  ncol = 4, nrow=2, align = "hv",  # align horizontally & vertically
                                  axis  = "tblr", hjust=-.5)
figure4panels_models


ggsave(filename = c("Figure4_NAcLSGRABDA2h_Bar_TransientPanels_ModelGraphs.pdf"), path = figurepath, plot = figure4panels_models, device = "pdf", 
       width = figpanelwidthin, height = figpanelheightinmodel, dpi=figdpi, units="in")




# FIGURE 5 - SENSOR TRANSIENT COMPARISON -----
## STATS ------
viruses <- c('GCaMP6f','dLight1.3b','GRABDA2h')
DoseNum_c <- c(-3.75, -1.25, 1.25, 3.75)
bins <- subjectmeansbin %>% filter(CondNum == 1) %>% distinct(BinF) %>% pull(BinF)

doselevels_c <- list(DoseNum_c = c(-3.75, -1.25, 1.25, 3.75))
subjectmeansbin$DoseNum_c <- subjectmeansbin$DoseNum - mean(unique(subjectmeansbin$DoseNum))

### Frequency SNR -----
# Mixed effects model
freqSNRmodel_sensorcomp <- glmmTMB(pksperminSNR_virus ~ DoseNum_c*Virus*BinF + (1|SubjectID), # random intercept by subject
                               data = subset(subjectmeansbin, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                               control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqSNRmodel_sensorcomp_aov <- Anova(freqSNRmodel_sensorcomp, type=2)
freqSNRmodel_sensorcomp_aov
summary(freqSNRmodel_sensorcomp)

freqSNRmodel_sensorcomp_r2 <- as.data.frame(r.squaredGLMM(freqSNRmodel_sensorcomp)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqSNRmodel_sensorcomp_r2

# Follow up tests - joint
freqSNRmodel_sensorcomp_emmeans <-emmeans(freqSNRmodel_sensorcomp, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
freqSNRmodel_sensorcomp_jt_Dose <- joint_tests(freqSNRmodel_sensorcomp_emmeans, by = "DoseNum_c")
freqSNRmodel_sensorcomp_jt_Virus <- joint_tests(freqSNRmodel_sensorcomp_emmeans, by = "Virus")

freqSNRmodel_sensorcomp_jt_Dose
freqSNRmodel_sensorcomp_jt_Virus

# Follow up tests - virus pairs
freqSNRmodel_sensorcomp_viruspairs <- emmeans(freqSNRmodel_sensorcomp, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
freqSNRmodel_sensorcomp_viruscontrasts <- as.data.frame(freqSNRmodel_sensorcomp_viruspairs$contrasts) # Extract emmeans
freqSNRmodel_sensorcomp_virusemmeans <- freqSNRmodel_sensorcomp_viruspairs$emmeans # Extract emmeans
freqSNRmodel_sensorcomp_viruspairs_es <- as.data.frame(eff_size(freqSNRmodel_sensorcomp_virusemmeans, sigma = sigma(freqSNRmodel_sensorcomp), edf = df.residual(freqSNRmodel_sensorcomp))) # Effect Size
freqSNRmodel_sensorcomp_viruspairs_es <- freqSNRmodel_sensorcomp_viruspairs_es[, c("contrast", "effect.size")]

freqSNRmodel_sensorcomp_viruscontrasts$d <- freqSNRmodel_sensorcomp_viruspairs_es$effect.size
freqSNRmodel_sensorcomp_viruscontrasts

# Format Tables
table_freqSNRmodel_sensorcomp <- format_regsummary_table(freqSNRmodel_sensorcomp, 'SensorComp', "FrequencySNR", "StandardModel")
table_freqSNRmodel_sensorcomp_aov <- format_anova_table(freqSNRmodel_sensorcomp_aov, "SensorComp", "FrequencySNR", "StandardModel")
table_freqSNRmodel_sensorcomp_joint_Dose <- format_emm_table(freqSNRmodel_sensorcomp_jt_Dose, "SensorComp", "FrequencySNR", "Joint", "Dose")
table_freqSNRmodel_sensorcomp_joint_Virus <- format_emm_table(freqSNRmodel_sensorcomp_jt_Virus, "SensorComp", "FrequencySNR", "Joint", "Virus")
table_freqSNRmodel_sensorcomp_pairs_Virus <- format_emm_table(freqSNRmodel_sensorcomp_viruscontrasts, "SensorComp", "FrequencySNR", "Pairs","Virus")

#### Mixed effects AR1 model
freqSNRmodel_sensorcomp_ar1 <- glmmTMB(pksperminSNR_virus ~ DoseNum_c*Virus*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                   data = subset(subjectmeansbin, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                   control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

freqSNRmodel_sensorcomp_aov_ar1 <- Anova(freqSNRmodel_sensorcomp_ar1, type=2)
freqSNRmodel_sensorcomp_aov_ar1
summary(freqSNRmodel_sensorcomp_ar1)

freqSNRmodel_sensorcomp_ar1_r2 <- as.data.frame(r.squaredGLMM(freqSNRmodel_sensorcomp_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
freqSNRmodel_sensorcomp_ar1_r2

freqSNRmodel_sensorcomp_ar1_AIC <- as.data.frame(AIC(freqSNRmodel_sensorcomp, freqSNRmodel_sensorcomp_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "freqSNRmodel_sensorcomp" ~ "StandardModel", ModelName == "freqSNRmodel_sensorcomp_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

freqSNRmodel_sensorcomp_ar1_AIC

# Follow up tests - joint
freqSNRmodel_sensorcomp_emmeans_ar1 <-emmeans(freqSNRmodel_sensorcomp_ar1, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
freqSNRmodel_sensorcomp_jt_Dose_ar1 <- joint_tests(freqSNRmodel_sensorcomp_emmeans_ar1, by = "DoseNum_c")
freqSNRmodel_sensorcomp_jt_Virus_ar1 <- joint_tests(freqSNRmodel_sensorcomp_emmeans_ar1, by = "Virus")

freqSNRmodel_sensorcomp_jt_Dose_ar1
freqSNRmodel_sensorcomp_jt_Virus_ar1

# Follow up tests - virus pairs
freqSNRmodel_sensorcomp_viruspairs_ar1 <- emmeans(freqSNRmodel_sensorcomp_ar1, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
freqSNRmodel_sensorcomp_viruscontrasts_ar1 <- as.data.frame(freqSNRmodel_sensorcomp_viruspairs_ar1$contrasts) # Extract emmeans
freqSNRmodel_sensorcomp_virusemmeans_ar1 <- freqSNRmodel_sensorcomp_viruspairs_ar1$emmeans # Extract emmeans
freqSNRmodel_sensorcomp_viruspairs_es_ar1 <- as.data.frame(eff_size(freqSNRmodel_sensorcomp_virusemmeans_ar1, sigma = sigma(freqSNRmodel_sensorcomp_ar1), edf = df.residual(freqSNRmodel_sensorcomp_ar1))) # Effect Size
freqSNRmodel_sensorcomp_viruspairs_es_ar1 <- freqSNRmodel_sensorcomp_viruspairs_es_ar1[, c("contrast", "effect.size")]

freqSNRmodel_sensorcomp_viruscontrasts_ar1$d <- freqSNRmodel_sensorcomp_viruspairs_es_ar1$effect.size
freqSNRmodel_sensorcomp_viruscontrasts_ar1


## Prepare fit data for plotting
freqmodel_sensorcomp_regdata <- expand.grid(Virus = viruses, BinF = bins, DoseNum_c = doselevels_c$DoseNum_c)
freqmodel_sensorcomp_regpred <- predict(freqSNRmodel_sensorcomp_ar1, newdata = freqmodel_sensorcomp_regdata, se.fit = TRUE, re.form = NA, type = "response")

freqmodel_sensorcomp_regdata$fit <- freqmodel_sensorcomp_regpred$fit
freqmodel_sensorcomp_regdata$se <- freqmodel_sensorcomp_regpred$se

freqmodel_sensorcomp_regdata <- freqmodel_sensorcomp_regdata %>%
  group_by(DoseNum_c, Virus) %>% summarise(fit = mean(fit), se = mean(se), .groups = 'drop')

freqmodel_sensorcomp_regdata$lower <- freqmodel_sensorcomp_regdata$fit - 1.96 * freqmodel_sensorcomp_regdata$se
freqmodel_sensorcomp_regdata$upper <- freqmodel_sensorcomp_regdata$fit + 1.96 * freqmodel_sensorcomp_regdata$se

freqmodel_sensorcomp_regdata$Dose <- freqmodel_sensorcomp_regdata$DoseNum_c + dosemean
freqmodel_sensorcomp_regdata

# Compare slopes
freqmodel_sensorcomp_regslopes <- emtrends(freqSNRmodel_sensorcomp_ar1, ~ Virus, var = "DoseNum_c")

freqmodel_sensorcomp_regslopes_pairs <- pairs(freqmodel_sensorcomp_regslopes)
freqmodel_sensorcomp_regslopes_pairs

freqmodel_sensorcomp_regslopes <- freqmodel_sensorcomp_regslopes %>% as.data.frame() %>% rename(Sensor = Virus, DoseTrend = DoseNum_c.trend)

# Format Tables
table_freqSNRmodel_sensorcomp_ar1 <- format_regsummary_table(freqSNRmodel_sensorcomp_ar1, 'SensorComp', "FrequencySNR", "AR1Model")
table_freqSNRmodel_sensorcomp_aov_ar1 <- format_anova_table(freqSNRmodel_sensorcomp_aov_ar1, "SensorComp", "FrequencySNR", "AR1Model")
table_freqSNRmodel_sensorcomp_joint_Dose_ar1 <- format_emm_table(freqSNRmodel_sensorcomp_jt_Dose_ar1, "SensorComp", "FrequencySNR", "Joint", "Dose")
table_freqSNRmodel_sensorcomp_joint_Virus_ar1 <- format_emm_table(freqSNRmodel_sensorcomp_jt_Virus_ar1, "SensorComp", "FrequencySNR", "Joint", "Virus")
table_freqSNRmodel_sensorcomp_pairs_Virus_ar1 <- format_emm_table(freqSNRmodel_sensorcomp_viruscontrasts_ar1, "SensorComp", "FrequencySNR", "Pairs","Virus")
table_freqSNRmodel_sensorcomp_virusreg <- format_regmultigroup_table(freqmodel_sensorcomp_regdata, "FrequencySNR", "DoseReg")
table_freqSNRmodel_sensorcomp_virusregslopes <- format_slopes_table(freqmodel_sensorcomp_regslopes, "FrequencySNR", "SimpleSlope", "Dose")
table_freqSNRmodel_sensorcomp_virusregslopepairs <- format_emm_table(freqmodel_sensorcomp_regslopes_pairs, "SensorComp", "FrequencySNR", "Pairs", "Virus")


### Amplitude SNR -----
# Mixed effects model
ampSNRmodel_sensorcomp <- glmmTMB(ampSNR_virus ~ DoseNum_c*Virus*BinF + (1|SubjectID), # random intercept by subject
                                   data = subset(subjectmeansbin, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                   control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampSNRmodel_sensorcomp_aov <- Anova(ampSNRmodel_sensorcomp, type=2)
ampSNRmodel_sensorcomp_aov
summary(ampSNRmodel_sensorcomp)

ampSNRmodel_sensorcomp_r2 <- as.data.frame(r.squaredGLMM(ampSNRmodel_sensorcomp)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampSNRmodel_sensorcomp_r2

# Follow up tests - joint
ampSNRmodel_sensorcomp_emmeans <-emmeans(ampSNRmodel_sensorcomp, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
ampSNRmodel_sensorcomp_jt_Dose <- joint_tests(ampSNRmodel_sensorcomp_emmeans, by = "DoseNum_c")
ampSNRmodel_sensorcomp_jt_Virus <- joint_tests(ampSNRmodel_sensorcomp_emmeans, by = "Virus")

ampSNRmodel_sensorcomp_jt_Dose
ampSNRmodel_sensorcomp_jt_Virus

# Follow up tests - virus pairs
ampSNRmodel_sensorcomp_viruspairs <- emmeans(ampSNRmodel_sensorcomp, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
ampSNRmodel_sensorcomp_viruscontrasts <- as.data.frame(ampSNRmodel_sensorcomp_viruspairs$contrasts) # Extract emmeans
ampSNRmodel_sensorcomp_virusemmeans <- ampSNRmodel_sensorcomp_viruspairs$emmeans # Extract emmeans
ampSNRmodel_sensorcomp_viruspairs_es <- as.data.frame(eff_size(ampSNRmodel_sensorcomp_virusemmeans, sigma = sigma(ampSNRmodel_sensorcomp), edf = df.residual(ampSNRmodel_sensorcomp))) # Effect Size
ampSNRmodel_sensorcomp_viruspairs_es <- ampSNRmodel_sensorcomp_viruspairs_es[, c("contrast", "effect.size")]

ampSNRmodel_sensorcomp_viruscontrasts$d <- ampSNRmodel_sensorcomp_viruspairs_es$effect.size
ampSNRmodel_sensorcomp_viruscontrasts

## Prepare fit data for plotting
ampmodel_sensorcomp_regdata <- expand.grid(Virus = viruses, BinF = bins, DoseNum_c = doselevels_c$DoseNum_c)
ampmodel_sensorcomp_regpred <- predict(ampSNRmodel_sensorcomp, newdata = ampmodel_sensorcomp_regdata, se.fit = TRUE, re.form = NA, type = "response")

ampmodel_sensorcomp_regdata$fit <- ampmodel_sensorcomp_regpred$fit
ampmodel_sensorcomp_regdata$se <- ampmodel_sensorcomp_regpred$se

ampmodel_sensorcomp_regdata <- ampmodel_sensorcomp_regdata %>%
  group_by(DoseNum_c, Virus) %>% summarise(fit = mean(fit), se = mean(se), .groups = 'drop')

ampmodel_sensorcomp_regdata$lower <- ampmodel_sensorcomp_regdata$fit - 1.96 * ampmodel_sensorcomp_regdata$se
ampmodel_sensorcomp_regdata$upper <- ampmodel_sensorcomp_regdata$fit + 1.96 * ampmodel_sensorcomp_regdata$se

ampmodel_sensorcomp_regdata$Dose <- ampmodel_sensorcomp_regdata$DoseNum_c + dosemean
ampmodel_sensorcomp_regdata

# Compare slopes
ampmodel_sensorcomp_regslopes <- emtrends(ampSNRmodel_sensorcomp, ~ Virus, var = "DoseNum_c")
ampmodel_sensorcomp_regslopes_pairs <- pairs(ampmodel_sensorcomp_regslopes)
ampmodel_sensorcomp_regslopes_pairs

ampmodel_sensorcomp_regslopes <- ampmodel_sensorcomp_regslopes %>% as.data.frame() %>% rename(Sensor = Virus, DoseTrend = DoseNum_c.trend)

# Format Tables
table_ampSNRmodel_sensorcomp <- format_regsummary_table(ampSNRmodel_sensorcomp, 'SensorComp', "AmplitudeSNR", "StandardModel")
table_ampSNRmodel_sensorcomp_aov <- format_anova_table(ampSNRmodel_sensorcomp_aov, "SensorComp", "AmplitudeSNR", "StandardModel")
table_ampSNRmodel_sensorcomp_joint_Dose <- format_emm_table(ampSNRmodel_sensorcomp_jt_Dose, "SensorComp", "AmplitudeSNR", "Joint", "Dose")
table_ampSNRmodel_sensorcomp_joint_Virus <- format_emm_table(ampSNRmodel_sensorcomp_jt_Virus, "SensorComp", "AmplitudeSNR", "Joint", "Virus")
table_ampSNRmodel_sensorcomp_pairs_Virus <- format_emm_table(ampSNRmodel_sensorcomp_viruscontrasts, "SensorComp", "AmplitudeSNR", "Pairs","Virus")

table_ampSNRmodel_sensorcomp_virusreg <- format_regmultigroup_table(ampmodel_sensorcomp_regdata, "AmplitudeSNR", "DoseReg")
table_ampSNRmodel_sensorcomp_virusregslopepairs <- format_emm_table(ampmodel_sensorcomp_regslopes_pairs, "SensorComp", "FrequencySNR", "Pairs", "Virus")
table_ampSNRmodel_sensorcomp_virusregslopes <- format_slopes_table(ampmodel_sensorcomp_regslopes, "AmplitudeSNR", "SimpleSlope", "Dose")


#### Mixed effects AR1 model
ampSNRmodel_sensorcomp_ar1 <- glmmTMB(ampSNR_virus ~ DoseNum_c*Virus*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                       data = subset(subjectmeansbin, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                       control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

ampSNRmodel_sensorcomp_aov_ar1 <- Anova(ampSNRmodel_sensorcomp_ar1, type=2)
ampSNRmodel_sensorcomp_aov_ar1
summary(ampSNRmodel_sensorcomp_ar1)

ampSNRmodel_sensorcomp_ar1_r2 <- as.data.frame(r.squaredGLMM(ampSNRmodel_sensorcomp_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
ampSNRmodel_sensorcomp_ar1_r2

ampSNRmodel_sensorcomp_ar1_AIC <- as.data.frame(AIC(ampSNRmodel_sensorcomp, ampSNRmodel_sensorcomp_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "ampSNRmodel_sensorcomp" ~ "StandardModel", ModelName == "ampSNRmodel_sensorcomp_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

ampSNRmodel_sensorcomp_ar1_AIC

# Follow up tests - joint
ampSNRmodel_sensorcomp_emmeans_ar1 <-emmeans(ampSNRmodel_sensorcomp_ar1, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
ampSNRmodel_sensorcomp_jt_Dose_ar1 <- joint_tests(ampSNRmodel_sensorcomp_emmeans_ar1, by = "DoseNum_c")
ampSNRmodel_sensorcomp_jt_Virus_ar1 <- joint_tests(ampSNRmodel_sensorcomp_emmeans_ar1, by = "Virus")

ampSNRmodel_sensorcomp_jt_Dose_ar1
ampSNRmodel_sensorcomp_jt_Virus_ar1

# Follow up tests - virus pairs
ampSNRmodel_sensorcomp_viruspairs_ar1 <- emmeans(ampSNRmodel_sensorcomp_ar1, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
ampSNRmodel_sensorcomp_viruscontrasts_ar1 <- as.data.frame(ampSNRmodel_sensorcomp_viruspairs_ar1$contrasts) # Extract emmeans
ampSNRmodel_sensorcomp_virusemmeans_ar1 <- ampSNRmodel_sensorcomp_viruspairs_ar1$emmeans # Extract emmeans
ampSNRmodel_sensorcomp_viruspairs_es_ar1 <- as.data.frame(eff_size(ampSNRmodel_sensorcomp_virusemmeans_ar1, sigma = sigma(ampSNRmodel_sensorcomp_ar1), edf = df.residual(ampSNRmodel_sensorcomp_ar1))) # Effect Size
ampSNRmodel_sensorcomp_viruspairs_es_ar1 <- ampSNRmodel_sensorcomp_viruspairs_es_ar1[, c("contrast", "effect.size")]

ampSNRmodel_sensorcomp_viruscontrasts_ar1$d <- ampSNRmodel_sensorcomp_viruspairs_es_ar1$effect.size
ampSNRmodel_sensorcomp_viruscontrasts_ar1


# Format Tables
table_ampSNRmodel_sensorcomp_ar1 <- format_regsummary_table(ampSNRmodel_sensorcomp_ar1, 'SensorComp', "AmplitudeSNR", "AR1Model")
table_ampSNRmodel_sensorcomp_aov_ar1 <- format_anova_table(ampSNRmodel_sensorcomp_aov_ar1, "SensorComp", "AmplitudeSNR", "AR1Model")
table_ampSNRmodel_sensorcomp_joint_Dose_ar1 <- format_emm_table(ampSNRmodel_sensorcomp_jt_Dose_ar1, "SensorComp", "AmplitudeSNR", "Joint", "Dose")
table_ampSNRmodel_sensorcomp_joint_Virus_ar1 <- format_emm_table(ampSNRmodel_sensorcomp_jt_Virus_ar1, "SensorComp", "AmplitudeSNR", "Joint", "Virus")
table_ampSNRmodel_sensorcomp_pairs_Virus_ar1 <- format_emm_table(ampSNRmodel_sensorcomp_viruscontrasts_ar1, "SensorComp", "AmplitudeSNR", "Pairs","Virus")


### AUCwindow SNR -----
# Mixed effects model
AUCwindowSNRmodel_sensorcomp <- glmmTMB(AUCwindowSNR_virus ~ DoseNum_c*Virus*BinF + (1|SubjectID), # random intercept by subject
                                  data = subset(subjectmeansbin, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                  control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowSNRmodel_sensorcomp_aov <- Anova(AUCwindowSNRmodel_sensorcomp, type=2)
AUCwindowSNRmodel_sensorcomp_aov
summary(AUCwindowSNRmodel_sensorcomp)

AUCwindowSNRmodel_sensorcomp_r2 <- as.data.frame(r.squaredGLMM(AUCwindowSNRmodel_sensorcomp)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowSNRmodel_sensorcomp_r2

# Follow up tests - joint
AUCwindowSNRmodel_sensorcomp_emmeans <-emmeans(AUCwindowSNRmodel_sensorcomp, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
AUCwindowSNRmodel_sensorcomp_jt_Dose <- joint_tests(AUCwindowSNRmodel_sensorcomp_emmeans, by = "DoseNum_c")
AUCwindowSNRmodel_sensorcomp_jt_Virus <- joint_tests(AUCwindowSNRmodel_sensorcomp_emmeans, by = "Virus")

AUCwindowSNRmodel_sensorcomp_jt_Dose
AUCwindowSNRmodel_sensorcomp_jt_Virus

# Follow up tests - virus pairs
AUCwindowSNRmodel_sensorcomp_viruspairs <- emmeans(AUCwindowSNRmodel_sensorcomp, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
AUCwindowSNRmodel_sensorcomp_viruscontrasts <- as.data.frame(AUCwindowSNRmodel_sensorcomp_viruspairs$contrasts) # Extract emmeans
AUCwindowSNRmodel_sensorcomp_virusemmeans <- AUCwindowSNRmodel_sensorcomp_viruspairs$emmeans # Extract emmeans
AUCwindowSNRmodel_sensorcomp_viruspairs_es <- as.data.frame(eff_size(AUCwindowSNRmodel_sensorcomp_virusemmeans, sigma = sigma(AUCwindowSNRmodel_sensorcomp), edf = df.residual(AUCwindowSNRmodel_sensorcomp))) # Effect Size
AUCwindowSNRmodel_sensorcomp_viruspairs_es <- AUCwindowSNRmodel_sensorcomp_viruspairs_es[, c("contrast", "effect.size")]

AUCwindowSNRmodel_sensorcomp_viruscontrasts$d <- AUCwindowSNRmodel_sensorcomp_viruspairs_es$effect.size
AUCwindowSNRmodel_sensorcomp_viruscontrasts

## Prepare fit data for plotting
AUCwindowmodel_sensorcomp_regdata <- expand.grid(Virus = viruses, BinF = bins, DoseNum_c = doselevels_c$DoseNum_c)
AUCwindowmodel_sensorcomp_regpred <- predict(AUCwindowSNRmodel_sensorcomp, newdata = AUCwindowmodel_sensorcomp_regdata, se.fit = TRUE, re.form = NA, type = "response")

AUCwindowmodel_sensorcomp_regdata$fit <- AUCwindowmodel_sensorcomp_regpred$fit
AUCwindowmodel_sensorcomp_regdata$se <- AUCwindowmodel_sensorcomp_regpred$se

AUCwindowmodel_sensorcomp_regdata <- AUCwindowmodel_sensorcomp_regdata %>%
  group_by(DoseNum_c, Virus) %>% summarise(fit = mean(fit), se = mean(se), .groups = 'drop')

AUCwindowmodel_sensorcomp_regdata$lower <- AUCwindowmodel_sensorcomp_regdata$fit - 1.96 * AUCwindowmodel_sensorcomp_regdata$se
AUCwindowmodel_sensorcomp_regdata$upper <- AUCwindowmodel_sensorcomp_regdata$fit + 1.96 * AUCwindowmodel_sensorcomp_regdata$se

AUCwindowmodel_sensorcomp_regdata$Dose <- AUCwindowmodel_sensorcomp_regdata$DoseNum_c + dosemean
AUCwindowmodel_sensorcomp_regdata

# Compare slopes
AUCwindowmodel_sensorcomp_regslopes <- emtrends(AUCwindowSNRmodel_sensorcomp, ~ Virus, var = "DoseNum_c")
AUCwindowmodel_sensorcomp_regslopes_pairs <- pairs(AUCwindowmodel_sensorcomp_regslopes)
AUCwindowmodel_sensorcomp_regslopes_pairs

AUCwindowmodel_sensorcomp_regslopes <- AUCwindowmodel_sensorcomp_regslopes %>% as.data.frame() %>% rename(Sensor = Virus, DoseTrend = DoseNum_c.trend)

# Format Tables
table_AUCwindowSNRmodel_sensorcomp <- format_regsummary_table(AUCwindowSNRmodel_sensorcomp, 'SensorComp', "AUCwindowSNR", "StandardModel")
table_AUCwindowSNRmodel_sensorcomp_aov <- format_anova_table(AUCwindowSNRmodel_sensorcomp_aov, "SensorComp", "AUCwindowSNR", "StandardModel")
table_AUCwindowSNRmodel_sensorcomp_joint_Dose <- format_emm_table(AUCwindowSNRmodel_sensorcomp_jt_Dose, "SensorComp", "AUCwindowSNR", "Joint", "Dose")
table_AUCwindowSNRmodel_sensorcomp_joint_Virus <- format_emm_table(AUCwindowSNRmodel_sensorcomp_jt_Virus, "SensorComp", "AUCwindowSNR", "Joint", "Virus")
table_AUCwindowSNRmodel_sensorcomp_pairs_Virus <- format_emm_table(AUCwindowSNRmodel_sensorcomp_viruscontrasts, "SensorComp", "AUCwindowSNR", "Pairs","Virus")
table_AUCwindowSNRmodel_sensorcomp_virusreg <- format_regmultigroup_table(AUCwindowmodel_sensorcomp_regdata, "AUCwindowSNR", "DoseReg")
table_AUCwindowSNRmodel_sensorcomp_virusregslopepairs <- format_emm_table(AUCwindowmodel_sensorcomp_regslopes_pairs, "SensorComp", "FrequencySNR", "Pairs", "Virus")
table_AUCwindowSNRmodel_sensorcomp_virusregslopes <- format_slopes_table(AUCwindowmodel_sensorcomp_regslopes, "AmplitudeSNR", "SimpleSlope", "Dose")


#### Mixed effects AR1 model
AUCwindowSNRmodel_sensorcomp_ar1 <- glmmTMB(AUCwindowSNR_virus ~ DoseNum_c*Virus*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                      data = subset(subjectmeansbin, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                      control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

AUCwindowSNRmodel_sensorcomp_aov_ar1 <- Anova(AUCwindowSNRmodel_sensorcomp_ar1, type=2)
AUCwindowSNRmodel_sensorcomp_aov_ar1
summary(AUCwindowSNRmodel_sensorcomp_ar1)

AUCwindowSNRmodel_sensorcomp_ar1_r2 <- as.data.frame(r.squaredGLMM(AUCwindowSNRmodel_sensorcomp_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
AUCwindowSNRmodel_sensorcomp_ar1_r2

AUCwindowSNRmodel_sensorcomp_ar1_AIC <- as.data.frame(AIC(AUCwindowSNRmodel_sensorcomp, AUCwindowSNRmodel_sensorcomp_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "AUCwindowSNRmodel_sensorcomp" ~ "StandardModel", ModelName == "AUCwindowSNRmodel_sensorcomp_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

AUCwindowSNRmodel_sensorcomp_ar1_AIC

# Follow up tests - joint
AUCwindowSNRmodel_sensorcomp_emmeans_ar1 <-emmeans(AUCwindowSNRmodel_sensorcomp_ar1, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
AUCwindowSNRmodel_sensorcomp_jt_Dose_ar1 <- joint_tests(AUCwindowSNRmodel_sensorcomp_emmeans_ar1, by = "DoseNum_c")
AUCwindowSNRmodel_sensorcomp_jt_Virus_ar1 <- joint_tests(AUCwindowSNRmodel_sensorcomp_emmeans_ar1, by = "Virus")

AUCwindowSNRmodel_sensorcomp_jt_Dose_ar1
AUCwindowSNRmodel_sensorcomp_jt_Virus_ar1

# Follow up tests - virus pairs
AUCwindowSNRmodel_sensorcomp_viruspairs_ar1 <- emmeans(AUCwindowSNRmodel_sensorcomp_ar1, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
AUCwindowSNRmodel_sensorcomp_viruscontrasts_ar1 <- as.data.frame(AUCwindowSNRmodel_sensorcomp_viruspairs_ar1$contrasts) # Extract emmeans
AUCwindowSNRmodel_sensorcomp_virusemmeans_ar1 <- AUCwindowSNRmodel_sensorcomp_viruspairs_ar1$emmeans # Extract emmeans
AUCwindowSNRmodel_sensorcomp_viruspairs_es_ar1 <- as.data.frame(eff_size(AUCwindowSNRmodel_sensorcomp_virusemmeans_ar1, sigma = sigma(AUCwindowSNRmodel_sensorcomp_ar1), edf = df.residual(AUCwindowSNRmodel_sensorcomp_ar1))) # Effect Size
AUCwindowSNRmodel_sensorcomp_viruspairs_es_ar1 <- AUCwindowSNRmodel_sensorcomp_viruspairs_es_ar1[, c("contrast", "effect.size")]

AUCwindowSNRmodel_sensorcomp_viruscontrasts_ar1$d <- AUCwindowSNRmodel_sensorcomp_viruspairs_es_ar1$effect.size
AUCwindowSNRmodel_sensorcomp_viruscontrasts_ar1


# Format Tables
table_AUCwindowSNRmodel_sensorcomp_ar1 <- format_regsummary_table(AUCwindowSNRmodel_sensorcomp_ar1, 'SensorComp', "AUCwindowSNR", "AR1Model")
table_AUCwindowSNRmodel_sensorcomp_aov_ar1 <- format_anova_table(AUCwindowSNRmodel_sensorcomp_aov_ar1, "SensorComp", "AUCwindowSNR", "AR1Model")
table_AUCwindowSNRmodel_sensorcomp_joint_Dose_ar1 <- format_emm_table(AUCwindowSNRmodel_sensorcomp_jt_Dose_ar1, "SensorComp", "AUCwindowSNR", "Joint", "Dose")
table_AUCwindowSNRmodel_sensorcomp_joint_Virus_ar1 <- format_emm_table(AUCwindowSNRmodel_sensorcomp_jt_Virus_ar1, "SensorComp", "AUCwindowSNR", "Joint", "Virus")
table_AUCwindowSNRmodel_sensorcomp_pairs_Virus_ar1 <- format_emm_table(AUCwindowSNRmodel_sensorcomp_viruscontrasts_ar1, "SensorComp", "AUCwindowSNR", "Pairs","Virus")


## EXPORT STATS -----
### Prepare Means -----
# Frequency
sensorcompmeanstable_freqSNR <- treatmeans %>% filter(InjType=='M') %>%
  select(FiberPlacement, Virus, Dose, DoseNum, TreatNum, InjType, CondNum, 
         pksperminSNR_virus_mean, pksperminSNR_virus_sd, pksperminSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% ungroup() %>% select(!c(InjType, DoseNum, TreatNum))


# Amplitude
sensorcompmeanstable_ampSNR <- treatmeans %>% filter(InjType=='M') %>%
  select(FiberPlacement, Virus, Dose, DoseNum, TreatNum, InjType, CondNum, 
         ampSNR_virus_mean, ampSNR_virus_sd, ampSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% ungroup() %>% select(!c(InjType, DoseNum, TreatNum))


# AUC Window
sensorcompmeanstable_AUCwindowSNR <- treatmeans %>% filter(InjType=='M') %>%
  select(FiberPlacement, Virus, Dose, DoseNum, TreatNum, InjType, CondNum, 
         AUCwindowSNR_virus_mean, AUCwindowSNR_virus_sd, AUCwindowSNR_virus_se) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>% ungroup() %>% select(!c(InjType, DoseNum, TreatNum))


### Frequency Export -----
wb_sensorcomp_freq <- createWorkbook()
wb_sensorcomp_freq_ar1 <- createWorkbook()

# Standard Models
write_nice_sheet(wb_sensorcomp_freq, "Means", sensorcompmeanstable_freqSNR)
write_nice_sheet(wb_sensorcomp_freq, "ModelSummary", table_freqSNRmodel_sensorcomp)
write_nice_sheet(wb_sensorcomp_freq, "ModelAnova", table_freqSNRmodel_sensorcomp_aov)
write_nice_sheet(wb_sensorcomp_freq, "ModelR2", freqSNRmodel_sensorcomp_r2)
write_nice_sheet(wb_sensorcomp_freq, "ModelAIC", freqSNRmodel_sensorcomp_ar1_AIC)
write_nice_sheet(wb_sensorcomp_freq, "Joint_Dose",  table_freqSNRmodel_sensorcomp_joint_Dose)
write_nice_sheet(wb_sensorcomp_freq, "Joint_Virus", table_freqSNRmodel_sensorcomp_joint_Virus)
write_nice_sheet(wb_sensorcomp_freq, "Pairs_Virus",  table_freqSNRmodel_sensorcomp_pairs_Virus)

saveWorkbook(wb_sensorcomp_freq, paste(analysispath,"SensorComparison_FrequencySNR_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# AR1 Models
write_nice_sheet(wb_sensorcomp_freq_ar1, "Means", sensorcompmeanstable_freqSNR)
write_nice_sheet(wb_sensorcomp_freq_ar1, "ModelSummary", table_freqSNRmodel_sensorcomp_ar1)
write_nice_sheet(wb_sensorcomp_freq_ar1, "ModelAnova", table_freqSNRmodel_sensorcomp_aov_ar1)
write_nice_sheet(wb_sensorcomp_freq_ar1, "ModelR2", freqSNRmodel_sensorcomp_ar1_r2)
write_nice_sheet(wb_sensorcomp_freq_ar1, "ModelAIC", freqSNRmodel_sensorcomp_ar1_AIC)
write_nice_sheet(wb_sensorcomp_freq_ar1, "Joint_Dose",  table_freqSNRmodel_sensorcomp_joint_Dose_ar1)
write_nice_sheet(wb_sensorcomp_freq_ar1, "Joint_Virus", table_freqSNRmodel_sensorcomp_joint_Virus_ar1)
write_nice_sheet(wb_sensorcomp_freq_ar1, "Pairs_Virus",  table_freqSNRmodel_sensorcomp_pairs_Virus_ar1)
write_nice_sheet(wb_sensorcomp_freq_ar1, "Slopes_Virus",  table_freqSNRmodel_sensorcomp_virusregslopes)
write_nice_sheet(wb_sensorcomp_freq_ar1, "SlopePairs_Virus",  table_freqSNRmodel_sensorcomp_virusregslopepairs)

saveWorkbook(wb_sensorcomp_freq_ar1, paste(analysispath,"SensorComparison_FrequencySNR_ar1models_results.xlsx",sep=''), overwrite = TRUE)

### Amplitude Export -----
wb_sensorcomp_amp <- createWorkbook()
wb_sensorcomp_amp_ar1 <- createWorkbook()

# Standard Models
write_nice_sheet(wb_sensorcomp_amp, "Means", sensorcompmeanstable_ampSNR)
write_nice_sheet(wb_sensorcomp_amp, "ModelSummary", table_ampSNRmodel_sensorcomp)
write_nice_sheet(wb_sensorcomp_amp, "ModelAnova", table_ampSNRmodel_sensorcomp_aov)
write_nice_sheet(wb_sensorcomp_amp, "ModelR2", ampSNRmodel_sensorcomp_r2)
write_nice_sheet(wb_sensorcomp_amp, "ModelAIC", ampSNRmodel_sensorcomp_ar1_AIC)
write_nice_sheet(wb_sensorcomp_amp, "Joint_Dose",  table_ampSNRmodel_sensorcomp_joint_Dose)
write_nice_sheet(wb_sensorcomp_amp, "Joint_Virus", table_ampSNRmodel_sensorcomp_joint_Virus)
write_nice_sheet(wb_sensorcomp_amp, "Pairs_Virus",  table_ampSNRmodel_sensorcomp_pairs_Virus)
write_nice_sheet(wb_sensorcomp_amp, "Slopes_Virus",  table_ampSNRmodel_sensorcomp_virusregslopes)
write_nice_sheet(wb_sensorcomp_amp, "SlopePairs_Virus",  table_ampSNRmodel_sensorcomp_virusregslopepairs)

saveWorkbook(wb_sensorcomp_amp, paste(analysispath,"SensorComparison_AmplitudeSNR_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# AR1 Models
write_nice_sheet(wb_sensorcomp_amp_ar1, "Means", sensorcompmeanstable_ampSNR)
write_nice_sheet(wb_sensorcomp_amp_ar1, "ModelSummary", table_ampSNRmodel_sensorcomp_ar1)
write_nice_sheet(wb_sensorcomp_amp_ar1, "ModelAnova", table_ampSNRmodel_sensorcomp_aov_ar1)
write_nice_sheet(wb_sensorcomp_amp_ar1, "ModelR2", ampSNRmodel_sensorcomp_ar1_r2)
write_nice_sheet(wb_sensorcomp_amp_ar1, "ModelAIC", ampSNRmodel_sensorcomp_ar1_AIC)
write_nice_sheet(wb_sensorcomp_amp_ar1, "Joint_Dose",  table_ampSNRmodel_sensorcomp_joint_Dose_ar1)
write_nice_sheet(wb_sensorcomp_amp_ar1, "Joint_Virus", table_ampSNRmodel_sensorcomp_joint_Virus_ar1)
write_nice_sheet(wb_sensorcomp_amp_ar1, "Pairs_Virus",  table_ampSNRmodel_sensorcomp_pairs_Virus_ar1)

saveWorkbook(wb_sensorcomp_amp_ar1, paste(analysispath,"SensorComparison_AmplitudeSNR_ar1models_results.xlsx",sep=''), overwrite = TRUE)

### AUCwindow Export -----
wb_sensorcomp_AUCwindow <- createWorkbook()
wb_sensorcomp_AUCwindow_ar1 <- createWorkbook()

# Standard Models
write_nice_sheet(wb_sensorcomp_AUCwindow, "Means", sensorcompmeanstable_AUCwindowSNR)
write_nice_sheet(wb_sensorcomp_AUCwindow, "ModelSummary", table_AUCwindowSNRmodel_sensorcomp)
write_nice_sheet(wb_sensorcomp_AUCwindow, "ModelAnova", table_AUCwindowSNRmodel_sensorcomp_aov)
write_nice_sheet(wb_sensorcomp_AUCwindow, "ModelR2", AUCwindowSNRmodel_sensorcomp_r2)
write_nice_sheet(wb_sensorcomp_AUCwindow, "ModelAIC", AUCwindowSNRmodel_sensorcomp_ar1_AIC)
write_nice_sheet(wb_sensorcomp_AUCwindow, "Joint_Dose",  table_AUCwindowSNRmodel_sensorcomp_joint_Dose)
write_nice_sheet(wb_sensorcomp_AUCwindow, "Joint_Virus", table_AUCwindowSNRmodel_sensorcomp_joint_Virus)
write_nice_sheet(wb_sensorcomp_AUCwindow, "Pairs_Virus",  table_AUCwindowSNRmodel_sensorcomp_pairs_Virus)
write_nice_sheet(wb_sensorcomp_AUCwindow, "Slopes_Virus",  table_AUCwindowSNRmodel_sensorcomp_virusregslopes)
write_nice_sheet(wb_sensorcomp_AUCwindow, "SlopePairs_Virus",  table_AUCwindowSNRmodel_sensorcomp_virusregslopepairs)

saveWorkbook(wb_sensorcomp_AUCwindow, paste(analysispath,"SensorComparison_AUCwindowSNR_standardmodels_results.xlsx",sep=''), overwrite = TRUE)

# AR1 Models
write_nice_sheet(wb_sensorcomp_AUCwindow_ar1, "Means", sensorcompmeanstable_AUCwindowSNR)
write_nice_sheet(wb_sensorcomp_AUCwindow_ar1, "ModelSummary", table_AUCwindowSNRmodel_sensorcomp_ar1)
write_nice_sheet(wb_sensorcomp_AUCwindow_ar1, "ModelAnova", table_AUCwindowSNRmodel_sensorcomp_aov_ar1)
write_nice_sheet(wb_sensorcomp_AUCwindow_ar1, "ModelR2", AUCwindowSNRmodel_sensorcomp_ar1_r2)
write_nice_sheet(wb_sensorcomp_AUCwindow_ar1, "ModelAIC", AUCwindowSNRmodel_sensorcomp_ar1_AIC)
write_nice_sheet(wb_sensorcomp_AUCwindow_ar1, "Joint_Dose",  table_AUCwindowSNRmodel_sensorcomp_joint_Dose_ar1)
write_nice_sheet(wb_sensorcomp_AUCwindow_ar1, "Joint_Virus", table_AUCwindowSNRmodel_sensorcomp_joint_Virus_ar1)
write_nice_sheet(wb_sensorcomp_AUCwindow_ar1, "Pairs_Virus",  table_AUCwindowSNRmodel_sensorcomp_pairs_Virus_ar1)

saveWorkbook(wb_sensorcomp_AUCwindow_ar1, paste(analysispath,"SensorComparison_AUCwindowSNR_ar1models_results.xlsx",sep=''), overwrite = TRUE)



## PLOTS -----
treatmeans$VirusDose <- paste(treatmeans$Virus, treatmeans$Dose, sep='_')
subjectmeans$VirusDose <- paste(subjectmeans$Virus, subjectmeans$Dose, sep='_')

reglinexbreaks <- c(2.5,5,7.5,10)
reglinexlabels <- c('2.5','5','7.5','10')

# SYMBOL KEY
# *** = sig from GCaMP6f
# ### = sig from dLight1.3b

### Frequency SNR -----
# SNR Virus
sensorcomp_freqbar_min <- -1.2
sensorcomp_freqbar_max <- 3.6
sensorcomp_freqbar_ticks <- 1.2
sensorcomp_freqbar_breaks <- seq(sensorcomp_freqbar_min, sensorcomp_freqbar_max, sensorcomp_freqbar_ticks)

# Bar Graph
sensorcomp_freqbar <- ggplot(subset(treatmeans, CondNum==1&InjType=='M'),aes(x=Dose, y=pksperminSNR_virus_mean, group=Virus, fill=VirusDose)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_errorbar(data=subset(treatmeans, CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = Virus, ymin=pksperminSNR_virus_mean-pksperminSNR_virus_se, ymax=pksperminSNR_virus_mean+pksperminSNR_virus_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = sensorcompcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcomp_freqbar_min, sensorcomp_freqbar_max), 
                     breaks=sensorcomp_freqbar_breaks) +
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" Frequency (SNR Scaled)")) + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(1), xmax = c(1.315), y.position=c(1.4),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(2), xmax = c(2.315), y.position=c(1.4),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(1.685), xmax = c(2.315), y.position=c(1.8),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.685), xmax = c(3.315), y.position=c(2.9),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3), xmax = c(3.315), y.position=c(2.5),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.685), xmax = c(4.315), y.position=c(3.4),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(4), xmax = c(4.315), y.position=c(3),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcomp_freqbar

# SNR Virus Regression
sensorcomp_freqreg_min <- -1.8
sensorcomp_freqreg_max <- 2.2
sensorcomp_freqreg_ticks <- .9
sensorcomp_freqreg_breaks <- seq(sensorcomp_freqreg_min, sensorcomp_freqreg_max, sensorcomp_freqreg_ticks)

sensorcomp_freqreg <- ggplot(table_freqSNRmodel_sensorcomp_virusreg,aes(x=Dose, y=fit, group=Sensor, color=Sensor)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_ribbon(data = table_freqSNRmodel_sensorcomp_virusreg, aes(x = Dose, ymin = CIlower, ymax = CIupper, fill=Sensor), inherit.aes = FALSE, alpha = myregalpha) +
  geom_line(size = mylinewidth) +
  scale_color_manual(values = overallcolors) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcomp_freqreg_min, sensorcomp_freqreg_max), 
                     breaks=sensorcomp_freqreg_breaks) +
  scale_x_continuous(expand=c(0,0.3), name='Dose (mg/kg)', limits=c(2.5,10.3), breaks = reglinexbreaks, labels = reglinexlabels) + 
  ylab(expression(Delta*" Frequency (SNR Scaled)")) + 
  mytheme +
  geom_bracket(xmin=c(3), xmax = c(4), y.position=c(2.1),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(3), xmax = c(5), y.position=c(1.7),label.size=mystarsize-2, label = "###",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .1)

sensorcomp_freqreg

### Amplitude SNR -----
# SNR Virus
sensorcomp_ampbar_min <- -1.6
sensorcomp_ampbar_max <- 1.6
sensorcomp_ampbar_ticks <- .8
sensorcomp_ampbar_breaks <- seq(sensorcomp_ampbar_min, sensorcomp_ampbar_max, sensorcomp_ampbar_ticks)

# Bar Graph
sensorcomp_ampbar <- ggplot(subset(treatmeans, CondNum==1&InjType=='M'),aes(x=Dose, y=ampSNR_virus_mean, group=Virus, fill=VirusDose)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_errorbar(data=subset(treatmeans, CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = Virus, ymin=ampSNR_virus_mean-ampSNR_virus_se, ymax=ampSNR_virus_mean+ampSNR_virus_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = sensorcompcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcomp_ampbar_min, sensorcomp_ampbar_max), 
                     breaks=sensorcomp_ampbar_breaks) +
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" Amplitude (SNR Scaled)")) + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(0.685), xmax = c(.96), y.position=c(1.1),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1.04), xmax = c(1.315), y.position=c(1.1),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1.685), xmax = c(1.96), y.position=c(1.1),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(2.04), xmax = c(2.315), y.position=c(1.1),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  
  geom_bracket(xmin=c(2.685), xmax = c(3.315), y.position=c(1.5),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3), xmax = c(3.315), y.position=c(1.25),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(3.685), xmax = c(4.315), y.position=c(1.5),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)+
  geom_bracket(xmin=c(4), xmax = c(4.315), y.position=c(1.25),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcomp_ampbar

# SNR Virus Regression
sensorcomp_ampreg_min <- -2
sensorcomp_ampreg_max <- 1.2
sensorcomp_ampreg_ticks <- 1
sensorcomp_ampreg_breaks <- seq(sensorcomp_ampreg_min, sensorcomp_ampreg_max, sensorcomp_ampreg_ticks)

sensorcomp_ampreg <- ggplot(table_ampSNRmodel_sensorcomp_virusreg,aes(x=Dose, y=fit, group=Sensor, color=Sensor)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_ribbon(data = table_ampSNRmodel_sensorcomp_virusreg, aes(x = Dose, ymin = CIlower, ymax = CIupper, fill=Sensor), inherit.aes = FALSE, alpha = myregalpha) +
  geom_line(size = mylinewidth) +
  scale_color_manual(values = overallcolors) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcomp_ampreg_min, sensorcomp_ampreg_max), 
                     breaks=sensorcomp_ampreg_breaks) +
  scale_x_continuous(expand=c(0,0.3), name='Dose (mg/kg)', limits=c(2.5,10.3), breaks = reglinexbreaks, labels = reglinexlabels) + 
  ylab(expression(Delta*" Amplitude (SNR Scaled)")) + 
  mytheme +
  geom_bracket(xmin=c(3), xmax = c(3.4), y.position=c(-1.6),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(4), xmax = c(5.5), y.position=c(-1.6),label.size=mystarsize-2, label = "##",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .1)

sensorcomp_ampreg


### AUCwindow SNR -----
# SNR Virus
sensorcomp_AUCwindowbar_min <- -1.6
sensorcomp_AUCwindowbar_max <- 1.6
sensorcomp_AUCwindowbar_ticks <- .8
sensorcomp_AUCwindowbar_breaks <- seq(sensorcomp_AUCwindowbar_min, sensorcomp_AUCwindowbar_max, sensorcomp_AUCwindowbar_ticks)

# Bar Graph
sensorcomp_AUCwindowbar <- ggplot(subset(treatmeans, CondNum==1&InjType=='M'),aes(x=Dose, y=AUCwindowSNR_virus_mean, group=Virus, fill=VirusDose)) + 
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=mycolwidth, show.legend=FALSE) + 
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_errorbar(data=subset(treatmeans, CondNum==1&InjType=='M'), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = Virus, ymin=AUCwindowSNR_virus_mean-AUCwindowSNR_virus_se, ymax=AUCwindowSNR_virus_mean+AUCwindowSNR_virus_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth) +
  scale_fill_manual(values = sensorcompcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcomp_AUCwindowbar_min, sensorcomp_AUCwindowbar_max), 
                     breaks=sensorcomp_AUCwindowbar_breaks) +
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" AUC (SNR Scaled)")) + 
  mytheme +
  # STATS
  geom_bracket(xmin=c(0.685), xmax = c(1), y.position=c(1),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(1.685), xmax = c(2), y.position=c(1),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcomp_AUCwindowbar

# SNR Virus Regression
sensorcomp_AUCwindowreg_min <- -1.5
sensorcomp_AUCwindowreg_max <- 1.4
sensorcomp_AUCwindowreg_ticks <- .7
sensorcomp_AUCwindowreg_breaks <- seq(-1.4, sensorcomp_AUCwindowreg_max, sensorcomp_AUCwindowreg_ticks)

sensorcomp_AUCwindowreg <- ggplot(table_AUCwindowSNRmodel_sensorcomp_virusreg,aes(x=Dose, y=fit, group=Sensor, color=Sensor)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_ribbon(data = table_AUCwindowSNRmodel_sensorcomp_virusreg, aes(x = Dose, ymin = CIlower, ymax = CIupper, fill=Sensor), inherit.aes = FALSE, alpha = myregalpha) +
  geom_line(size = mylinewidth) +
  scale_color_manual(values = overallcolors) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(sensorcomp_AUCwindowreg_min, sensorcomp_AUCwindowreg_max), 
                     breaks=sensorcomp_AUCwindowreg_breaks) +
  scale_x_continuous(expand=c(0,0.3), name='Dose (mg/kg)', limits=c(2.5,10.3), breaks = reglinexbreaks, labels = reglinexlabels) + 
  ylab(expression(Delta*" AUC (SNR Scaled)")) + 
  mytheme +
  geom_bracket(xmin=c(3), xmax = c(3.4), y.position=c(1.2),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5)

sensorcomp_AUCwindowreg

### EXPORT FIGURE 5 -----
# BAR AND REG GRAPHS
figure5panels_barline <- plot_grid(sensorcomp_freqbar, sensorcomp_freqreg,
                                   sensorcomp_ampbar, sensorcomp_ampreg,
                                   sensorcomp_AUCwindowbar, sensorcomp_AUCwindowreg,
                                   ncol = 2, nrow=3, align = "hv",  # align horizontally & vertically
                                   axis  = "tblr", hjust=-.5, rel_widths = c(1,.6,
                                                                             1,.6,
                                                                             1,.6))
figure5panels_barline

ggsave(filename = c("Figure5_TransientsSensorComparison.pdf"), path = figurepath, plot = figure5panels_barline, device = "pdf", 
       width = figpanelwidthin, height = figpanelheightinbarline, dpi=figdpi, units="in")


# FIGURE 6 - WHOLE STREAM COMPARISON -----
## STATS ------
viruses <- c('GCaMP6f','dLight1.3b','GRABDA2h')
DoseNum_c <- c(-3.75, -1.25, 1.25, 3.75)
bins <- streambinmeans %>% filter(CondNum == 1) %>% distinct(BinF) %>% pull(BinF)

doselevels_c <- list(DoseNum_c = c(-3.75, -1.25, 1.25, 3.75))

### Mean (Change) -----
# Mixed effects model
streammeanCmodel_sigsub <- glmmTMB(meanC ~ DoseNum_c*Virus*BinF + (1|SubjectID), # random intercept by subject
                                   data = subset(streamdata, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                   control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

streammeanCmodel_sigsub_aov <- Anova(streammeanCmodel_sigsub, type=2)
streammeanCmodel_sigsub_aov
summary(streammeanCmodel_sigsub)

streammeanCmodel_sigsub_r2 <- as.data.frame(r.squaredGLMM(streammeanCmodel_sigsub)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
streammeanCmodel_sigsub_r2

# Follow up tests - joint
streammeanCmodel_sigsub_emmeans <-emmeans(streammeanCmodel_sigsub, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
streammeanCmodel_sigsub_jt_Dose <- joint_tests(streammeanCmodel_sigsub_emmeans, by = "DoseNum_c")
streammeanCmodel_sigsub_jt_Virus <- joint_tests(streammeanCmodel_sigsub_emmeans, by = "Virus")

streammeanCmodel_sigsub_jt_Dose
streammeanCmodel_sigsub_jt_Virus

# Follow up tests - virus pairs
streammeanCmodel_sigsub_viruspairs <- emmeans(streammeanCmodel_sigsub, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
streammeanCmodel_sigsub_viruscontrasts <- as.data.frame(streammeanCmodel_sigsub_viruspairs$contrasts) # Extract emmeans
streammeanCmodel_sigsub_virusemmeans <- streammeanCmodel_sigsub_viruspairs$emmeans # Extract emmeans
streammeanCmodel_sigsub_viruspairs_es <- as.data.frame(eff_size(streammeanCmodel_sigsub_virusemmeans, sigma = sigma(streammeanCmodel_sigsub), edf = df.residual(streammeanCmodel_sigsub))) # Effect Size
streammeanCmodel_sigsub_viruspairs_es <- streammeanCmodel_sigsub_viruspairs_es[, c("contrast", "effect.size")]

streammeanCmodel_sigsub_viruscontrasts$d <- streammeanCmodel_sigsub_viruspairs_es$effect.size
streammeanCmodel_sigsub_viruscontrasts

# Format Tables
table_streammeanCmodel_sigsub <- format_regsummary_table(streammeanCmodel_sigsub, 'SensorComp', "MeanC", "StandardModel")
table_streammeanCmodel_sigsub_aov <- format_anova_table(streammeanCmodel_sigsub_aov, "SensorComp", "MeanC", "StandardModel")
table_streammeanCmodel_sigsub_joint_Dose <- format_emm_table(streammeanCmodel_sigsub_jt_Dose, "SensorComp", "MeanC", "Joint", "Dose")
table_streammeanCmodel_sigsub_joint_Virus <- format_emm_table(streammeanCmodel_sigsub_jt_Virus, "SensorComp", "MeanC", "Joint", "Virus")
table_streammeanCmodel_sigsub_pairs_Virus <- format_emm_table(streammeanCmodel_sigsub_viruscontrasts, "SensorComp", "MeanC", "Pairs","Virus")

#### Mixed effects AR1 model
streammeanCmodel_sigsub_ar1 <- glmmTMB(meanC ~ DoseNum_c*Virus*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                       data = subset(streamdata, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                       control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

streammeanCmodel_sigsub_aov_ar1 <- Anova(streammeanCmodel_sigsub_ar1, type=2)
streammeanCmodel_sigsub_aov_ar1
summary(streammeanCmodel_sigsub_ar1)

streammeanCmodel_sigsub_ar1_r2 <- as.data.frame(r.squaredGLMM(streammeanCmodel_sigsub_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
streammeanCmodel_sigsub_ar1_r2

streammeanCmodel_sigsub_ar1_AIC <- as.data.frame(AIC(streammeanCmodel_sigsub, streammeanCmodel_sigsub_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "streammeanCmodel_sigsub" ~ "StandardModel", ModelName == "streammeanCmodel_sigsub_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

streammeanCmodel_sigsub_ar1_AIC

# Follow up tests - joint
streammeanCmodel_sigsub_emmeans_ar1 <-emmeans(streammeanCmodel_sigsub_ar1, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
streammeanCmodel_sigsub_jt_Dose_ar1 <- joint_tests(streammeanCmodel_sigsub_emmeans_ar1, by = "DoseNum_c")
streammeanCmodel_sigsub_jt_Virus_ar1 <- joint_tests(streammeanCmodel_sigsub_emmeans_ar1, by = "Virus")

streammeanCmodel_sigsub_jt_Dose_ar1
streammeanCmodel_sigsub_jt_Virus_ar1

# Follow up tests - virus pairs
streammeanCmodel_sigsub_viruspairs_ar1 <- emmeans(streammeanCmodel_sigsub_ar1, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
streammeanCmodel_sigsub_viruscontrasts_ar1 <- as.data.frame(streammeanCmodel_sigsub_viruspairs_ar1$contrasts) # Extract emmeans
streammeanCmodel_sigsub_virusemmeans_ar1 <- streammeanCmodel_sigsub_viruspairs_ar1$emmeans # Extract emmeans
streammeanCmodel_sigsub_viruspairs_es_ar1 <- as.data.frame(eff_size(streammeanCmodel_sigsub_virusemmeans_ar1, sigma = sigma(streammeanCmodel_sigsub_ar1), edf = df.residual(streammeanCmodel_sigsub_ar1))) # Effect Size
streammeanCmodel_sigsub_viruspairs_es_ar1 <- streammeanCmodel_sigsub_viruspairs_es_ar1[, c("contrast", "effect.size")]

streammeanCmodel_sigsub_viruscontrasts_ar1$d <- streammeanCmodel_sigsub_viruspairs_es_ar1$effect.size
streammeanCmodel_sigsub_viruscontrasts_ar1


## Prepare fit data for plotting
streammeanCmodel_sigsub_regdata <- expand.grid(Virus = viruses, BinF = bins, DoseNum_c = doselevels_c$DoseNum_c)
streammeanCmodel_sigsub_regpred <- predict(streammeanCmodel_sigsub_ar1, newdata = streammeanCmodel_sigsub_regdata, se.fit = TRUE, re.form = NA, type = "response")

streammeanCmodel_sigsub_regdata$fit <- streammeanCmodel_sigsub_regpred$fit
streammeanCmodel_sigsub_regdata$se <- streammeanCmodel_sigsub_regpred$se

streammeanCmodel_sigsub_regdata <- streammeanCmodel_sigsub_regdata %>%
  group_by(DoseNum_c, Virus) %>% summarise(fit = mean(fit), se = mean(se), .groups = 'drop')

streammeanCmodel_sigsub_regdata$lower <- streammeanCmodel_sigsub_regdata$fit - 1.96 * streammeanCmodel_sigsub_regdata$se
streammeanCmodel_sigsub_regdata$upper <- streammeanCmodel_sigsub_regdata$fit + 1.96 * streammeanCmodel_sigsub_regdata$se

streammeanCmodel_sigsub_regdata$Dose <- streammeanCmodel_sigsub_regdata$DoseNum_c + dosemean
streammeanCmodel_sigsub_regdata

# Compare slopes
streammeanCmodel_sigsub_regslopes <- emtrends(streammeanCmodel_sigsub_ar1, ~ Virus, var = "DoseNum_c")

streammeanCmodel_sigsub_regslopes_pairs <- pairs(streammeanCmodel_sigsub_regslopes)
streammeanCmodel_sigsub_regslopes_pairs

streammeanCmodel_sigsub_regslopes <- streammeanCmodel_sigsub_regslopes %>% as.data.frame() %>% rename(Sensor = Virus, DoseTrend = DoseNum_c.trend)

# Format Tables
table_streammeanCmodel_sigsub_ar1 <- format_regsummary_table(streammeanCmodel_sigsub_ar1, 'SensorComp', "MeanC", "AR1Model")
table_streammeanCmodel_sigsub_aov_ar1 <- format_anova_table(streammeanCmodel_sigsub_aov_ar1, "SensorComp", "MeanC", "AR1Model")
table_streammeanCmodel_sigsub_joint_Dose_ar1 <- format_emm_table(streammeanCmodel_sigsub_jt_Dose_ar1, "SensorComp", "MeanC", "Joint", "Dose")
table_streammeanCmodel_sigsub_joint_Virus_ar1 <- format_emm_table(streammeanCmodel_sigsub_jt_Virus_ar1, "SensorComp", "MeanC", "Joint", "Virus")
table_streammeanCmodel_sigsub_pairs_Virus_ar1 <- format_emm_table(streammeanCmodel_sigsub_viruscontrasts_ar1, "SensorComp", "MeanC", "Pairs","Virus")
table_streammeanCmodel_sigsub_virusreg <- format_regmultigroup_table(streammeanCmodel_sigsub_regdata, "MeanC", "DoseReg")
table_streammeanCmodel_sigsub_virusregslopes <- format_slopes_table(streammeanCmodel_sigsub_regslopes, "MeanC", "SimpleSlope", "Dose")
table_streammeanCmodel_sigsub_virusregslopepairs <- format_emm_table(streammeanCmodel_sigsub_regslopes_pairs, "SensorComp", "MeanC", "Pairs", "Virus")



### MAD Ratio (Change) -----
# Mixed effects model
streamMADRatioCmodel_sigsub <- glmmTMB(MADRatioC ~ DoseNum_c*Virus*BinF + (1|SubjectID), # random intercept by subject
                                       data = subset(streamdata, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                       control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

streamMADRatioCmodel_sigsub_aov <- Anova(streamMADRatioCmodel_sigsub, type=2)
streamMADRatioCmodel_sigsub_aov
summary(streamMADRatioCmodel_sigsub)

streamMADRatioCmodel_sigsub_r2 <- as.data.frame(r.squaredGLMM(streamMADRatioCmodel_sigsub)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
streamMADRatioCmodel_sigsub_r2

# Follow up tests - joint
streamMADRatioCmodel_sigsub_emmeans <-emmeans(streamMADRatioCmodel_sigsub, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
streamMADRatioCmodel_sigsub_jt_Dose <- joint_tests(streamMADRatioCmodel_sigsub_emmeans, by = "DoseNum_c")
streamMADRatioCmodel_sigsub_jt_Virus <- joint_tests(streamMADRatioCmodel_sigsub_emmeans, by = "Virus")

streamMADRatioCmodel_sigsub_jt_Dose
streamMADRatioCmodel_sigsub_jt_Virus

# Follow up tests - virus pairs
streamMADRatioCmodel_sigsub_viruspairs <- emmeans(streamMADRatioCmodel_sigsub, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
streamMADRatioCmodel_sigsub_viruscontrasts <- as.data.frame(streamMADRatioCmodel_sigsub_viruspairs$contrasts) # Extract emmeans
streamMADRatioCmodel_sigsub_virusemmeans <- streamMADRatioCmodel_sigsub_viruspairs$emmeans # Extract emmeans
streamMADRatioCmodel_sigsub_viruspairs_es <- as.data.frame(eff_size(streamMADRatioCmodel_sigsub_virusemmeans, sigma = sigma(streamMADRatioCmodel_sigsub), edf = df.residual(streamMADRatioCmodel_sigsub))) # Effect Size
streamMADRatioCmodel_sigsub_viruspairs_es <- streamMADRatioCmodel_sigsub_viruspairs_es[, c("contrast", "effect.size")]

streamMADRatioCmodel_sigsub_viruscontrasts$d <- streamMADRatioCmodel_sigsub_viruspairs_es$effect.size
streamMADRatioCmodel_sigsub_viruscontrasts

# Format Tables
table_streamMADRatioCmodel_sigsub <- format_regsummary_table(streamMADRatioCmodel_sigsub, 'SensorComp', "MADRatioC", "StandardModel")
table_streamMADRatioCmodel_sigsub_aov <- format_anova_table(streamMADRatioCmodel_sigsub_aov, "SensorComp", "MADRatioC", "StandardModel")
table_streamMADRatioCmodel_sigsub_joint_Dose <- format_emm_table(streamMADRatioCmodel_sigsub_jt_Dose, "SensorComp", "MADRatioC", "Joint", "Dose")
table_streamMADRatioCmodel_sigsub_joint_Virus <- format_emm_table(streamMADRatioCmodel_sigsub_jt_Virus, "SensorComp", "MADRatioC", "Joint", "Virus")
table_streamMADRatioCmodel_sigsub_pairs_Virus <- format_emm_table(streamMADRatioCmodel_sigsub_viruscontrasts, "SensorComp", "MADRatioC", "Pairs","Virus")

#### Mixed effects AR1 model
streamMADRatioCmodel_sigsub_ar1 <- glmmTMB(MADRatioC ~ DoseNum_c*Virus*BinF + ar1(BinF + 0|SubjectID) + (1|SubjectID), # random intercept by subject
                                           data = subset(streamdata, CondNum == 1 & InjType=='M'), family = stats::gaussian, 
                                           control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))

streamMADRatioCmodel_sigsub_aov_ar1 <- Anova(streamMADRatioCmodel_sigsub_ar1, type=2)
streamMADRatioCmodel_sigsub_aov_ar1
summary(streamMADRatioCmodel_sigsub_ar1)

streamMADRatioCmodel_sigsub_ar1_r2 <- as.data.frame(r.squaredGLMM(streamMADRatioCmodel_sigsub_ar1)) %>% rename(R2_Marginal = R2m, R2_Conditional = R2c)
streamMADRatioCmodel_sigsub_ar1_r2

streamMADRatioCmodel_sigsub_ar1_AIC <- as.data.frame(AIC(streamMADRatioCmodel_sigsub, streamMADRatioCmodel_sigsub_ar1)) %>% rownames_to_column(var = "ModelName") %>% 
  mutate(Model = case_when(ModelName == "streamMADRatioCmodel_sigsub" ~ "StandardModel", ModelName == "streamMADRatioCmodel_sigsub_ar1" ~ "AR1Model", TRUE ~ ModelName)) %>%
  select(!ModelName) %>% relocate(Model) %>% add_row(Model = "DeltaAIC", AIC = diff(.$AIC))

streamMADRatioCmodel_sigsub_ar1_AIC

# Follow up tests - joint
streamMADRatioCmodel_sigsub_emmeans_ar1 <-emmeans(streamMADRatioCmodel_sigsub_ar1, specs = "Virus", at = doselevels_c, by = "DoseNum_c", type  = "response",  adjust="none") # Find means
streamMADRatioCmodel_sigsub_jt_Dose_ar1 <- joint_tests(streamMADRatioCmodel_sigsub_emmeans_ar1, by = "DoseNum_c")
streamMADRatioCmodel_sigsub_jt_Virus_ar1 <- joint_tests(streamMADRatioCmodel_sigsub_emmeans_ar1, by = "Virus")

streamMADRatioCmodel_sigsub_jt_Dose_ar1
streamMADRatioCmodel_sigsub_jt_Virus_ar1

# Follow up tests - virus pairs
streamMADRatioCmodel_sigsub_viruspairs_ar1 <- emmeans(streamMADRatioCmodel_sigsub_ar1, pairwise~Virus, by = "DoseNum_c", at = doselevels_c, adjust="none") # Pairs
streamMADRatioCmodel_sigsub_viruscontrasts_ar1 <- as.data.frame(streamMADRatioCmodel_sigsub_viruspairs_ar1$contrasts) # Extract emmeans
streamMADRatioCmodel_sigsub_virusemmeans_ar1 <- streamMADRatioCmodel_sigsub_viruspairs_ar1$emmeans # Extract emmeans
streamMADRatioCmodel_sigsub_viruspairs_es_ar1 <- as.data.frame(eff_size(streamMADRatioCmodel_sigsub_virusemmeans_ar1, sigma = sigma(streamMADRatioCmodel_sigsub_ar1), edf = df.residual(streamMADRatioCmodel_sigsub_ar1))) # Effect Size
streamMADRatioCmodel_sigsub_viruspairs_es_ar1 <- streamMADRatioCmodel_sigsub_viruspairs_es_ar1[, c("contrast", "effect.size")]

streamMADRatioCmodel_sigsub_viruscontrasts_ar1$d <- streamMADRatioCmodel_sigsub_viruspairs_es_ar1$effect.size
streamMADRatioCmodel_sigsub_viruscontrasts_ar1

## Prepare fit data for plotting
streamMADRatioCmodel_sigsub_regdata <- expand.grid(Virus = viruses, BinF = bins, DoseNum_c = doselevels_c$DoseNum_c)
streamMADRatioCmodel_sigsub_regpred <- predict(streamMADRatioCmodel_sigsub_ar1, newdata = streamMADRatioCmodel_sigsub_regdata, se.fit = TRUE, re.form = NA, type = "response")

streamMADRatioCmodel_sigsub_regdata$fit <- streamMADRatioCmodel_sigsub_regpred$fit
streamMADRatioCmodel_sigsub_regdata$se <- streamMADRatioCmodel_sigsub_regpred$se

streamMADRatioCmodel_sigsub_regdata <- streamMADRatioCmodel_sigsub_regdata %>%
  group_by(DoseNum_c, Virus) %>% summarise(fit = mean(fit), se = mean(se), .groups = 'drop')

streamMADRatioCmodel_sigsub_regdata$lower <- streamMADRatioCmodel_sigsub_regdata$fit - 1.96 * streamMADRatioCmodel_sigsub_regdata$se
streamMADRatioCmodel_sigsub_regdata$upper <- streamMADRatioCmodel_sigsub_regdata$fit + 1.96 * streamMADRatioCmodel_sigsub_regdata$se

streamMADRatioCmodel_sigsub_regdata$Dose <- streamMADRatioCmodel_sigsub_regdata$DoseNum_c + dosemean
streamMADRatioCmodel_sigsub_regdata

# Compare slopes
streamMADRatioCmodel_sigsub_regslopes <- emtrends(streamMADRatioCmodel_sigsub_ar1, ~ Virus, var = "DoseNum_c")

streamMADRatioCmodel_sigsub_regslopes_pairs <- pairs(streamMADRatioCmodel_sigsub_regslopes)
streamMADRatioCmodel_sigsub_regslopes_pairs

streamMADRatioCmodel_sigsub_regslopes <- streamMADRatioCmodel_sigsub_regslopes %>% as.data.frame() %>% rename(Sensor = Virus, DoseTrend = DoseNum_c.trend)

# Format Tables
table_streamMADRatioCmodel_sigsub_ar1 <- format_regsummary_table(streamMADRatioCmodel_sigsub_ar1, 'SensorComp', "MADRatioC", "AR1Model")
table_streamMADRatioCmodel_sigsub_aov_ar1 <- format_anova_table(streamMADRatioCmodel_sigsub_aov_ar1, "SensorComp", "MADRatioC", "AR1Model")
table_streamMADRatioCmodel_sigsub_joint_Dose_ar1 <- format_emm_table(streamMADRatioCmodel_sigsub_jt_Dose_ar1, "SensorComp", "MADRatioC", "Joint", "Dose")
table_streamMADRatioCmodel_sigsub_joint_Virus_ar1 <- format_emm_table(streamMADRatioCmodel_sigsub_jt_Virus_ar1, "SensorComp", "MADRatioC", "Joint", "Virus")
table_streamMADRatioCmodel_sigsub_pairs_Virus_ar1 <- format_emm_table(streamMADRatioCmodel_sigsub_viruscontrasts_ar1, "SensorComp", "MADRatioC", "Pairs","Virus")
table_streamMADRatioCmodel_sigsub_virusreg <- format_regmultigroup_table(streamMADRatioCmodel_sigsub_regdata, "MADRatioC", "DoseReg")
table_streamMADRatioCmodel_sigsub_virusregslopes <- format_slopes_table(streamMADRatioCmodel_sigsub_regslopes, "MADRatioC", "SimpleSlope", "Dose")
table_streamMADRatioCmodel_sigsub_virusregslopepairs <- format_emm_table(streamMADRatioCmodel_sigsub_regslopes_pairs, "SensorComp", "MADRatioC", "Pairs", "Virus")

## EXPORT STATS -----
### Prepare Means -----
# Stream Mean
treatmeanstable_mean <- streamoverallmeans %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  select(FiberPlacement, 
         Sensor = Virus, 
         InjType, 
         Dose,
         Condition,
         n = binmean_n,
         Mean = binmean_mean, 
         SD = binmean_sd, 
         SE = binmean_se)

# Stream MeanC
treatmeanstable_meanC <- streamoverallmeans %>% filter(CondNum == 1, InjType == 'M') %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  select(FiberPlacement, 
         Sensor = Virus, 
         InjType, 
         Dose,
         n = binmeanC_n,
         Mean = binmeanC_mean, 
         SD = binmeanC_sd, 
         SE = binmeanC_se)

# Stream MAD
treatmeanstable_MADRatio <- streamoverallmeans %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  select(FiberPlacement, 
         Sensor = Virus, 
         InjType, 
         Dose,
         Condition,
         n = binMADRatio_n,
         Mean = binMADRatio_mean, 
         SD = binMADRatio_sd, 
         SE = binMADRatio_se)

# Stream MADC
treatmeanstable_MADRatioC <- streamoverallmeans %>% filter(CondNum == 1, InjType == 'M') %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  select(FiberPlacement, 
         Sensor = Virus, 
         InjType, 
         Dose,
         n = binMADRatioC_n,
         Mean = binMADRatioC_mean, 
         SD = binMADRatioC_sd, 
         SE = binMADRatioC_se)

### Mean Export -----
wb_wholestream_meanC_ar1 <- createWorkbook()

# ar1 Models
write_nice_sheet(wb_wholestream_meanC_ar1, "MeansRaw", treatmeanstable_mean)
write_nice_sheet(wb_wholestream_meanC_ar1, "MeansChange", treatmeanstable_meanC)
write_nice_sheet(wb_wholestream_meanC_ar1, "ModelSummary", table_streammeanCmodel_sigsub_ar1)
write_nice_sheet(wb_wholestream_meanC_ar1, "ModelAnova", table_streammeanCmodel_sigsub_aov_ar1)
write_nice_sheet(wb_wholestream_meanC_ar1, "ModelR2", streammeanCmodel_sigsub_ar1_r2)
write_nice_sheet(wb_wholestream_meanC_ar1, "ModelAIC", streammeanCmodel_sigsub_ar1_AIC)
write_nice_sheet(wb_wholestream_meanC_ar1, "Joint_Dose", table_streammeanCmodel_sigsub_joint_Dose_ar1)
write_nice_sheet(wb_wholestream_meanC_ar1, "Joint_Virus",  table_streammeanCmodel_sigsub_joint_Virus_ar1)
write_nice_sheet(wb_wholestream_meanC_ar1, "Pairs_Virus", table_streammeanCmodel_sigsub_pairs_Virus_ar1)
write_nice_sheet(wb_wholestream_meanC_ar1, "Reg_Slopes",  table_streammeanCmodel_sigsub_virusregslopes)
write_nice_sheet(wb_wholestream_meanC_ar1, "Reg_Pairs",  table_streammeanCmodel_sigsub_virusregslopepairs)

saveWorkbook(wb_wholestream_meanC_ar1, paste(analysispath,"WholeStream_Mean_ar1models_results.xlsx",sep=''), overwrite = TRUE)

### Mean Export -----
wb_wholestream_MADRatioC_ar1 <- createWorkbook()

# ar1 Models
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "MeansRaw", treatmeanstable_MADRatio)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "MeansChange", treatmeanstable_MADRatioC)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "ModelSummary", table_streamMADRatioCmodel_sigsub_ar1)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "ModelAnova", table_streamMADRatioCmodel_sigsub_aov_ar1)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "ModelR2", streamMADRatioCmodel_sigsub_ar1_r2)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "ModelAIC", streamMADRatioCmodel_sigsub_ar1_AIC)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "Joint_Dose", table_streamMADRatioCmodel_sigsub_joint_Dose_ar1)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "Joint_Virus",  table_streamMADRatioCmodel_sigsub_joint_Virus_ar1)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "Pairs_Virus", table_streamMADRatioCmodel_sigsub_pairs_Virus_ar1)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "Reg_Slopes",  table_streamMADRatioCmodel_sigsub_virusregslopes)
write_nice_sheet(wb_wholestream_MADRatioC_ar1, "Reg_Pairs",  table_streamMADRatioCmodel_sigsub_virusregslopepairs)

saveWorkbook(wb_wholestream_MADRatioC_ar1, paste(analysispath,"WholeStream_MADRatio_ar1models_results.xlsx",sep=''), overwrite = TRUE)


## PLOTS -----
## Set up graphing variables
treatcolors <- colorlists$treatcolors
treatcolors_raw <- colorlists$treatcolors_raw
overallcolors <- colorlists$overallcolors
sensorcompcolors <- colorlists$sensorcompcolors
sensordosecolors <- colorlists$sensordosecolors

injlinecolor <- '#575757'

sexshapes <- c('M'=0,'F'=1)

bintimexticks <- seq(-12,60,3)
bintimexlabels <- ifelse(bintimexticks %% 6 == 0, bintimexticks, "")

bintimexmin <- -12
bintimexmax <- 60

reglinexbreaks <- c(2.5,5,7.5,10)
reglinexlabels <- c('2.5','5','7.5','10')

### Mean (Change) -----
## Bar Graph (Change)
meanCbar_ymin <- -3.1
meanCbar_ymax <- 3.5
meanCbar_yticks <- 1.5
meanCbar_breaks <- seq(-3, meanCbar_ymax, meanCbar_yticks)

meanCbar <- ggplot(subset(streamoverallmeans,InjType=='M'&CondNum==1), aes(x=Dose, y=binmeanC_mean, group=VirusDosecolor, fill=VirusDosecolor, color=VirusDosecolor)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=myskinnycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=subset(streamoverallmeans,InjType=='M'&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = VirusDosecolor, ymin=binmeanC_mean-binmeanC_se, ymax=binmeanC_mean+binmeanC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth, color='#000000') +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(streamsubjectmeans,InjType=='M'&CondNum==1), aes(x=Dose, y=binmeanC, fill=VirusDosecolor, group=VirusDosecolor, shape=Sex), color='#4c4c4c', size=myskinnysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  scale_fill_manual(values = sensorcompcolors) +
  scale_color_manual(values = sensorcompcolors) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(meanCbar_ymin, meanCbar_ymax), breaks=meanCbar_breaks) +
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" Mean")) + 
  mytheme +
  geom_bracket(xmin=c(0.685), xmax = c(1.313), y.position=c(2.9),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = mybracketlinewidth, vjust = .5)+
  geom_bracket(xmin=c(1), xmax = c(1.313), y.position=c(2.2),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = mybracketlinewidth, vjust = .5)+
  geom_bracket(xmin=c(4), xmax = c(4.313), y.position=c(2.2),label.size=mystarsize, label = "*",inherit.aes = FALSE, tip.length = c(.03, .03), size = mybracketlinewidth, vjust = .5)

meanCbar

## Line Graph 
meanCline_ymin <- -1.5
meanCline_ymax <- 1.6
meanCline_yticks <- .75
meanCline_ybreaks <- seq(meanCline_ymin, meanCline_ymax, meanCline_yticks)

meanCline <- ggplot(subset(streambinmeans,InjType=='M'), aes(x=BinTimeInj, y=binmeanC_mean, group=VirusDosecolor, fill=VirusDosecolor, color=VirusDosecolor)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mysmallpointsize, show.legend = FALSE) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  scale_x_continuous(name="Time from Injection (min)",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(meanCline_ymin, meanCline_ymax), breaks = meanCline_ybreaks)+
  scale_color_manual("legend", values = sensorcompcolors) +
  scale_fill_manual("legend", values = sensorcompcolors) +
  mytheme +
  ylab(expression(Delta*" Mean")) + 
  mytheme

meanCline

# Regression
meanCreg_ymin <- -1.5
meanCreg_ymax <- 1.5
meanCreg_yticks <- .75
meanCreg_ybreaks <- seq(meanCreg_ymin,meanCreg_ymax,meanCreg_yticks)

meanCreg <- ggplot(table_streammeanCmodel_sigsub_virusreg, aes(x=Dose, y=fit, group=Sensor, color=Sensor)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_ribbon(data = table_streammeanCmodel_sigsub_virusreg, aes(x = Dose, ymin = CIlower, ymax = CIupper, fill=Sensor), inherit.aes = FALSE, alpha = myregalpha) +
  geom_line(size = mylinewidth) +
  scale_color_manual(values = overallcolors) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(meanCreg_ymin, meanCreg_ymax), breaks=meanCreg_ybreaks) +
  scale_x_continuous(expand=c(0,0.3), name='Dose (mg/kg)', limits=c(2.5,10.3), breaks = reglinexbreaks, labels = reglinexlabels) + 
  ylab(expression(Delta*" Mean")) + 
  mytheme +
  geom_bracket(xmin=c(5.4), xmax = c(6.5), y.position=c(-1),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .5) +
  geom_bracket(xmin=c(7), xmax = c(7.5), y.position=c(-1),label.size=mystarsize-2, label = "###",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .1)+
  geom_bracket(xmin=c(8), xmax = c(9.8), y.position=c(-1),label.size=mystarsize-2, label = "&&&",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .1)

meanCreg


### MAD Ratio (Change) -----
## Bar Graph (Change)
MADRatioCbar_ymin <- -.6
MADRatioCbar_ymax <- .6
MADRatioCbar_yticks <- .3
MADRatioCbar_breaks <- seq(MADRatioCbar_ymin, MADRatioCbar_ymax, MADRatioCbar_yticks)

MADRatioCbar <- ggplot(subset(streamoverallmeans,InjType=='M'&CondNum==1), aes(x=Dose, y=binMADRatioC_mean, group=VirusDosecolor, fill=VirusDosecolor, color=VirusDosecolor)) +
  geom_col(position = position_dodge(mydodgewidth), linewidth = mycollinewidth, width=myskinnycolwidth, show.legend=FALSE) + 
  geom_errorbar(data=subset(streamoverallmeans,InjType=='M'&CondNum==1), show.legend=FALSE, position=position_dodge(mydodgewidth),
                aes(x=Dose, group = VirusDosecolor, ymin=binMADRatioC_mean-binMADRatioC_se, ymax=binMADRatioC_mean+binMADRatioC_se), width=myerrorbarwidth, linewidth=myerrorbarlinewidth, color='#000000') +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_jitter(data=subset(streamsubjectmeans,InjType=='M'&CondNum==1), aes(x=Dose, y=binMADRatioC, fill=VirusDosecolor, group=VirusDosecolor, shape=Sex), color='#4c4c4c', size=myskinnysubjectpointsize, 
              position = position_jitterdodge(jitter.width = 0.4, jitter.height = 0, dodge.width  = mydodgewidth), show.legend=FALSE) +
  scale_fill_manual(values = sensorcompcolors) +
  scale_color_manual(values = sensorcompcolors) +
  scale_shape_manual(values = sexshapes) +
  scale_y_continuous(expand=c(0,0),limits = c(MADRatioCbar_ymin, MADRatioCbar_ymax), breaks=MADRatioCbar_breaks) +
  xlab("Dose (mg/kg)") +
  ylab(expression(Delta*" MAD Ratio")) + 
  mytheme +
  geom_bracket(xmin=c(1.685), xmax = c(2.313), y.position=c(.44),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = mybracketlinewidth, vjust = .5)+
  geom_bracket(xmin=c(2), xmax = c(2.313), y.position=c(.3),label.size=mystarsize, label = "**",inherit.aes = FALSE, tip.length = c(.03, .03), size = mybracketlinewidth, vjust = .5)+
  geom_bracket(xmin=c(2.685), xmax = c(3.313), y.position=c(.44),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = mybracketlinewidth, vjust = .5)+
  geom_bracket(xmin=c(3), xmax = c(3.313), y.position=c(.3),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = mybracketlinewidth, vjust = .5)+
  geom_bracket(xmin=c(3.685), xmax = c(4.313), y.position=c(.54),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = mybracketlinewidth, vjust = .5)+
  geom_bracket(xmin=c(4), xmax = c(4.313), y.position=c(.4),label.size=mystarsize, label = "***",inherit.aes = FALSE, tip.length = c(.03, .03), size = mybracketlinewidth, vjust = .5)

MADRatioCbar

## Line Graph 
MADRatioCline_ymin <- -.4
MADRatioCline_ymax <- .4
MADRatioCline_yticks <- .2
MADRatioCline_ybreaks <- seq(-.4, MADRatioCline_ymax, MADRatioCline_yticks)

MADRatioCline <- ggplot(subset(streambinmeans,InjType=='M'), aes(x=BinTimeInj, y=binMADRatioC_mean, group=VirusDosecolor, fill=VirusDosecolor, color=VirusDosecolor)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_line(stat="identity", size = mylinewidth, show.legend = FALSE) +
  geom_point(stat="identity", size = mysmallpointsize, show.legend = FALSE) +
  geom_vline(xintercept=0, linetype='dashed', col = 'black', linewidth = mylinewidth) +
  scale_x_continuous(name="Time from Injection (min)",expand=c(0,0),limits = c(bintimexmin, bintimexmax), breaks=bintimexticks, labels=bintimexlabels)+
  scale_y_continuous(expand = c(0, 0), limits = c(MADRatioCline_ymin, MADRatioCline_ymax), breaks = MADRatioCline_ybreaks)+
  scale_color_manual("legend", values = sensorcompcolors) +
  scale_fill_manual("legend", values = sensorcompcolors) +
  mytheme +
  ylab(expression(Delta*" MAD Ratio")) + 
  mytheme

MADRatioCline

# Regression
MADRatioCreg_ymin <- -.4
MADRatioCreg_ymax <- .4
MADRatioCreg_yticks <- .2
MADRatioCreg_ybreaks <- seq(MADRatioCreg_ymin,MADRatioCreg_ymax,MADRatioCreg_yticks)

MADRatioCreg <- ggplot(table_streamMADRatioCmodel_sigsub_virusreg, aes(x=Dose, y=fit, group=Sensor, color=Sensor)) +
  geom_hline(yintercept=0, linetype='solid', col = 'black', linewidth=.4) +
  geom_ribbon(data = table_streamMADRatioCmodel_sigsub_virusreg, aes(x = Dose, ymin = CIlower, ymax = CIupper, fill=Sensor), inherit.aes = FALSE, alpha = myregalpha) +
  geom_line(size = mylinewidth) +
  scale_color_manual(values = overallcolors) +
  scale_fill_manual(values = overallcolors) +
  scale_y_continuous(expand=c(0,0),limits = c(MADRatioCreg_ymin, MADRatioCreg_ymax), breaks=MADRatioCreg_ybreaks) +
  scale_x_continuous(expand=c(0,0.3), name='Dose (mg/kg)', limits=c(2.5,10.3), breaks = reglinexbreaks, labels = reglinexlabels) + 
  ylab(expression(Delta*" MAD Ratio")) + 
  mytheme +
  geom_bracket(xmin=c(1), xmax = c(5), y.position=c(.3),label.size=mystarsize-2, label = "###",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .1)+
  geom_bracket(xmin=c(5), xmax = c(8), y.position=c(.2),label.size=mystarsize-2, label = "&&&",inherit.aes = FALSE, tip.length = c(.03, .03), size = .4, vjust = .1)

MADRatioCreg


### EXPORT FIGURE 6 -----
# BAR AND REG GRAPHS
figure6panels <- plot_grid(meanCbar, meanCline, meanCreg,
                           MADRatioCbar, MADRatioCline, MADRatioCreg,
                           ncol = 3, nrow=2, align = "hv",  # align horizontally & vertically
                           axis  = "tblr", hjust=-.5, rel_widths = c(1,1,.6,
                                                                     1,1,.6,
                                                                     1,1,.6,
                                                                     1,1,.6))
figure6panels

ggsave(filename = c("Figure6_WholeStreamSensorComparison.pdf"), path = figurepath, plot = figure6panels, device = "pdf", 
       width = figpanelwidthin+.15, height = figpanelheightin/1.5, dpi=figdpi, units="in")

