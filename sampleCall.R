# source('~/smarties/microarray/epigenetics/epimutation_analysis/epimut.R')

# GEOgse <- "GSE104812"
# resultFolderValue <- "~/experiments/test/results"
# sampleFolderValue <- "~/experiments/test/samples"

# epimut(
#   slidingWindowSize = 11,
#   resultFolder = resultFolderValue,
#   sampleFolder = sampleFolderValue,
#   popReferenceSelector = "Case_Ctrl == 1",
#   popCaseSelector = "Case_Ctrl == 0",
#   popControlSelector = "Case_Ctrl == 1 | Case_Ctrl == 0" 
# )

source('~/smarties/microarray/epigenetics/epimutation_analysis/supportFunctionsEpimutation.R')
source('~/smarties/microarray/epigenetics/epimutation_analysis/analyzePopulation.R')
source('~/smarties/microarray/epigenetics/epimutation_analysis/epimutationNormalize.R')
source('~/smarties/microarray/epigenetics/epimutation_analysis/epimut.R')
source('~/smarties/microarray/commonR.R')

# resultFolderValue <- "~/experiments/IPOT/results"
# sampleFolderValue <- "~/experiments/IPOT/samples"
# logFolderValue <- "~/experiments/IPOT/logs"

resultFolderValue <- "~/experiments/IPOT/results"
sampleFolderValue <- "~/experiments/IPOT/samples"
logFolderValue <- "~/experiments/IPOT/logs"

# savedData <- paste0(logFolderValue, "/data.Rda")
# if(file.exists(savedData)) localResult <- readRDS(file = savedData)

epimut(
  slidingWindowSize = 11,
  resultFolder = resultFolderValue,
  sampleFolder = sampleFolderValue,
  logFolder = logFolderValue,
  popStudyDefinition  = data.frame("LABELS" = c("STUDY","STUDY2"), "SELECTORS" = c("Case_Ctrl == 1", "Case_Ctrl == 2")),
  popControlSelector = "Case_Ctrl == 0",
  methodNormalization = "BMIQ", # BMIQ, SWAN, PBC,
  useMValue = TRUE, # FALSE means use BETA Value, TRUE means use MVALUE
  nIterations = 1, # number of repetitions of experiment
  bonferroniThreshold = 0.05 #threshold to define whihc pValue accept for lesions definition
)

# epimut(
#   slidingWindowSize = 11,
#   resultFolder = resultFolderValue,
#   sampleFolder = sampleFolderValue,
#   logFolder = logFolderValue,
#   popReferenceSelector = "Sample_Group == 'POPREF' ",
#   popCaseSelector = "Sample_Group == 'CASI' & Case_Ctrl == 0",
#   popControlSelector = "Sample_Group == 'CASI' & Case_Ctrl == 1"
# )


# source('~/smarties/microarray/commonR.R')
# syncSampleSheetFromFolder(renameFileStartingWithNumber = FALSE,
#                           sampleSheet = paste(sampleFolderValue, "/", "METHYLATION_giu_samples_NOV_2017.csv", sep = ""),
#                           workingFolder = sampleFolderValue,
#                           columnSeparator = ","
#                           )

# buildSampleSheetFromFolder(workingFolder =  sampleFolderValue, renameFileStartingWithNumber = FALSE)


# source('~/smarties/microarray/commonR.R')
# linkSampleFilesToExperimentFolder(
#     sampleFolder = "/Volumes/Data/idat/",
#     experimentFolder = "~/experiments/IPOT/samples",
#     sampleSheet = "~/experiments/IPOT/METHYLATION_giu_samples_NOV_2017.csv"
# )


# source('~/smarties/microarray/commonR.R')
# createSampleFilesToExperimentFolder(
#     sampleFolder = "/Volumes/Data/idat/",
#     experimentFolder = "~/experiments/IPOT/samples",
#     sampleSheet = "~/experiments/IPOT/METHYLATION_giu_samples_NOV_2017.csv"
# )