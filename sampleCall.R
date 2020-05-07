
source('supportFunctionsEpimutation.R')
source('analyzePopulation.R')
source('epimutationNormalize.R')
source('epimut.R')
source('commonR.R')

resultFolderValue <- "~/experiments/IPOT/results"
sampleFolderValue <- "~/experiments/IPOT/samples"
logFolderValue <- "~/experiments/IPOT/logs"

epimut(
  slidingWindowSize = 11, # sliding windows size to define if the locus is subject to an hyper/hypo methylation
  resultFolder = resultFolderValue,
  sampleFolder = sampleFolderValue,
  logFolder = logFolderValue,
  popStudyDefinition  = data.frame("LABELS" = c("STUDY","STUDY2"), "SELECTORS" = c("Case_Ctrl == 1", "Case_Ctrl == 2")), # LABELS and SELCTOR are used as the appera into the sample sheet
  popControlSelector = "Case_Ctrl == 0", ## define which subselection of the sample sheet is the population control
  methodNormalization = "BMIQ", # BMIQ, SWAN, PBC,
  useMValue = TRUE, # FALSE means use BETA Value, TRUE means use MVALUE
  nIterations = 1, # number of repetitions of experiment
  bonferroniThreshold = 0.05 #threshold to define whihc pValue accept for lesions definition
)
