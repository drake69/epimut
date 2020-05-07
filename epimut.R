source('supportFunctionsEpimutation.R')
source('analyzePopulation.R')
source('epimutationNormalize.R'
)

### EPIMUT ############################################################################################################################################################################

epimut <-
  function(sampleFolder,
           resultFolder,
           logFolder,
           slidingWindowSize ,
           methodNormalization,
           popControlSelector,
           popStudyDefinition,
           useMValue = FALSE,
           nIterations = 1, 
           bonferroniThreshold = 0.05) {
    
    list.of.packages <-
      c(
        "plyr",
        "zoo",
        "dplyr",
        "foreach",
        "R.utils",
        "parallel",
        "doParallel",
        "gtools",
        "openxlsx",
        "stringi",
        "reshape2",
        "ChAMP",
        "checkr"
      )
    
    new.packages <-
      list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
    if (length(new.packages))
      install.packages(new.packages)
    
    library(plyr)
    
    ### Check the Sample Sheet has the expected column name ##################
    sampleSheets <-
      list.files(path = sampleFolder, pattern = "*.csv")
    if (length(sampleSheets) > 1)
    {
      stop(" More than one sample sheet ", Sys.time())
    }
    
    if (length(sampleSheets) == 0)
    {
      stop(" At least one sample sheet is needed.", Sys.time())
    }
    
    sampleSheet <-
      read.table(
        paste0(sampleFolder, "/", sampleSheets[1], sep = ""),
        header = TRUE,
        sep = ","
      )
    needColumn <- c("Sample_Name", "Sentrix_ID", "Sentrix_Position")
    missedColumns <-
      needColumn[!(needColumn  %in% colnames(sampleSheet))]
    if (length(missedColumns) > 0)
    {
      stop("File:" ,
           sampleSheets[1]  ,
           " Lost following columns ",
           missedColumns,
           " ",
           Sys.time())
    }
    
    ### Normalize ##########################################################################
    
    if (logFolder != "" && !dir.exists(logFolder)) {
      dir.create(logFolder)
    }
    
    #############  repeat 10 times the calculations
    for (j in 1:nIterations)
    {
      # main result folder check
      system(paste0("rm -rf ", paste0(resultFolder, "_", j, sep = ""), sep = ""))
      dir.create(paste0(resultFolder, "_", j, sep = ""))
      
      savedData <- paste(logFolder, "/data.Rda", sep = "")
      if (file.exists(savedData))
      {
        localResult <- readRDS(file = savedData)
      }
      else
      {
        localResult <-
          epimut_normalize_champ(
            sampleFolder = sampleFolder,
            resultFolder = paste0(resultFolder,"_", j),
            methodNormalization = methodNormalization,
            useMvalue = useMValue
          )
        #saveRDS(localResult, file = savedData)
      }
      
      normalizedData =  localResult$normalizedData
      probeFeatures = localResult$probeFeatures
      sampleSheet = localResult$sampleSheet
      
      ### Reference Study ##########################################################################
      
      populations <-
        rbind(
          popStudyDefinition,
          data.frame("LABELS" = "CONTROL", "SELECTORS" = popControlSelector)
        )
      
      # control population
      controlPopulationSampleSheet <-
        subset(sampleSheet, eval(parse(text = popControlSelector)))
      controlPopulationMatrix <-
        data.frame("PROBE" = row.names(normalizedData), normalizedData[, controlPopulationSampleSheet$Sample_Name])
      
      if (empty(controlPopulationMatrix) |
          dim(controlPopulationMatrix)[2] < 2)
      {
        message("Empty normalizedData ", Sys.time())
        stop("Empty normalizedData ")
      }
      
      populationControlRangeBetaValues <-
        rangeBetaValuePerProbeAsNormalizedDistribution(controlPopulationMatrix, iqrTimes = 3)
      # populationControlRangeBetaValues <- rangeBetaValuePerCHROverCHRAsNormalizedDistribution(controlPopulationMatrix, iqrTimes = 3, probeFeatures)
      
      
      for (i in 1:dim(populations)[1])
      {
        if (is.null(populations[i, "SELECTORS"]))
        {
          next
        }
        
        popSelector <- as.character(populations[i, "SELECTORS"])
        popName <- as.character(populations[i, "LABELS"])
        
        #browser()
        populationSampleSheet <-
          subset(sampleSheet, eval(parse(text = popSelector)))
        populationMatrix <-
          data.frame("PROBE" = row.names(normalizedData), normalizedData[, populationSampleSheet$Sample_Name])
        
        if (empty(populationMatrix) | dim(populationMatrix)[2] < 2)
        {
          message("Population with Selector ",
                  popSelector ,
                  " is empty ",
                  Sys.time())
          next
        }
        
        analizePopulation(
          populationMatrix = populationMatrix,
          slidingWindowSize = slidingWindowSize,
          resultFolder = paste0(resultFolder, "_", j, sep = ""),
          logFolder = logFolder,
          betaSuperiorThresholds = populationControlRangeBetaValues$betaSuperiorThresholds,
          betaInferiorThresholds = populationControlRangeBetaValues$betaInferiorThresholds,
          sampleSheet = populationSampleSheet,
          probeFeatures = probeFeatures,
          betaMeans = populationControlRangeBetaValues$betaMedianValues,
          populationName = popName,
          bonferroniThreshold = bonferroniThreshold
        )
        
        rm(populationSampleSheet)
        rm(populationMatrix)
      }
      
      rm(populationControlRangeBetaValues)
    }
  }
