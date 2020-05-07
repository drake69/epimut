source("supportFunctionsEpimutation.R")

### analyzeSingleSample ############################################################################################################################################################################

analyzeSingleSample <- function(values, slidingWindowSize, resultFolder, thresholds, comparison, sampleName, probeFeatures, means, subFileExtension, bonferroniThreshold = 0.05) {

  start_time_single_sample <- Sys.time()
  message(sampleName," ", "... Sample analysis warmingUP ", Sys.time() )
  library(zoo, quietly = TRUE)
  result <- ""
  colnames(values) <- "VALUE"

  ### get probeFeatures ################################################################################################
  
  probePositions <- data.frame("PROBE" = probeFeatures$Probe,"CHR" = probeFeatures$CHR, "START" = probeFeatures$START, "END" = probeFeatures$START )
  colnames(probePositions) <- c("PROBE" ,"CHR" , "START" , "END" )
  
  message(sampleName," ","Sample analysis WarmedUP ...", Sys.time())
  message(sampleName," ","Start sample analyze ", Sys.time())
  
  ### get probesOverThreshold ################################################################################################
  
  mutation <- as.numeric(comparison(values,  thresholds))
  ##browser()
  message(sampleName," ","Got probesOverThreshold ", Sys.time() )

  ### get mutationAnnotatedSorted ################################################################################################
  if (!test_match_order(row.names(mutation), probePositions$Probe))
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  
  mutationAnnotated <- data.frame( as.data.frame(probePositions),"MUTATIONS" = mutation)
  mutationAnnotatedSorted <- sortByCHRandSTART(mutationAnnotated)
  probePositions <- sortByCHRandSTART(probePositions)
  
  if (!test_match_order(row.names(mutationAnnotatedSorted), probePositions$Probe))
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  
  #browser()
  if (!test_match_order(mutationAnnotatedSorted$Probe, row.names(mutationAnnotatedSorted)))
    stop("Mutation annotation sorted is not coherent with probe informations order!")  
  
  result <- c( "mutationCount" = sum(mutationAnnotatedSorted$MUTATIONS), "lesionCount" = 0, "probesCount" = 0)
  
  mutationAnnotatedSortedToSave <- subset(mutationAnnotatedSorted, MUTATIONS == 1)[,c("CHR","START","END")]
  message(sampleName," ","Got mutationAnnotatedSorted ", Sys.time() )
  
  dumpSampleAsBedFile(
    dataToDump = mutationAnnotatedSortedToSave,
    fileExtension = paste(".", subFileExtension , ".MUTATIONS.bed",sep = ""),
    resultFolder =  resultFolder,
    resultSubFolder =  paste( "MUTATIONS", subFileExtension,sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR","START","END","SAMPLENAME")
  )

  ### get lesion #################################################################################################

  #browser()
  chromosomes <- unique(attributes(mutationAnnotatedSorted$CHR)$levels)
  lesionWeighted <- data.frame()
  for (chrome in chromosomes)
  {
    if (chrome == "" | chrome == "X" | chrome == "Y")
      browser()
    #browser(condition = chrome == "X")
    message(sampleName, " Working on Chromosome ", chrome, " ", Sys.time())
    mutationAnnotatedSortedTemp <-    subset(mutationAnnotatedSorted,CHR == chrome)
    
    if (plyr::empty(mutationAnnotatedSortedTemp))
      browser()
    probePositionsTemp <- subset(probePositions, CHR == chrome)
    lesionWeightedTemp <- getLesions(mutationAnnotatedSorted = mutationAnnotatedSortedTemp, slidingWindowSize = slidingWindowSize , sampleName = sampleName, probePositions =  probePositionsTemp, bonferroniThreshold = bonferroniThreshold)
    if (sum(lesionWeightedTemp$LESIONS) > sum(mutationAnnotatedSortedTemp$MUTATIONS))
      browser()
    lesionWeighted <- rbind(lesionWeighted,lesionWeightedTemp)
  }
  
  result["lesionCount"] <- dim(lesionWeighted)[1]
  result["probesCount"] <- dim(probeFeatures)[1]
  # if (result["lesionCount"] > dim(mutationAnnotatedSortedToSave)[1])
  # {
  #   #browser()
  #   lesionWeighted <- getLesions(mutationAnnotatedSorted = mutationAnnotatedSorted, slidingWindowSize = slidingWindowSize , sampleName = sampleName, probePositions =  probePositions)
  #   result["lesionCount"] <- dim(lesionWeighted)[1]
  # }
  

  dumpSampleAsBedFile(
    dataToDump = lesionWeighted,
    fileExtension = paste0(".", subFileExtension , ".LESIONS.bed"),
    resultFolder =  resultFolder,
    resultSubFolder =  paste( "LESIONS", subFileExtension, sep = "_"),
    sampleName = sampleName,
    multipleFileColNames = c("CHR","START","END","SAMPLENAME")
  )
  
  end_time_single_sample <- Sys.time()
  time_taken <- end_time_single_sample - start_time_single_sample
  message(sampleName," ","Completed sample ", time_taken)
  return(result)
  # rm(list = ls())
}



### analyzeDeltaSingleSample ############################################################################################################################################################################

deltaSingleSample <- function(values,  resultFolder, highThresholds,lowThresholds, sampleName, probeFeatures, means, subFileExtension) {
  
  message(sampleName," ", "... Sample analysis warmingUP ", Sys.time() )
  library(zoo, quietly = TRUE)
  
  colnames(values) <- "VALUE"
  
  ### get probeFeatures ################################################################################################
  
  probePositions <- data.frame("PROBE" = probeFeatures$Probe,"CHR" = probeFeatures$CHR, "START" = probeFeatures$START, "END" = probeFeatures$START )
  
  message(sampleName," ","Sample analysis WarmedUP ...", Sys.time())
  message(sampleName," ","Start sample analyze ", Sys.time())
  
  ### get probesOverThreshold ################################################################################################
  
  mutationAbove <- values > highThresholds
  mutationBelow <- values < lowThresholds
  mutation <- (mutationBelow + mutationAbove) > 0
  colnames(mutation) <- "MUTATIONS"
  
  message(sampleName," ","Got outliers ", Sys.time() )
  
  ### get deltas #########################################################
  
  deltas <- values - means  
  colnames(deltas) <- "DELTA"
  
  if (!test_match_order(row.names(mutation), probePositions$Probe))
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  if (!test_match_order(row.names(deltas), probePositions$Probe))
    stop("Wrong order matching Probes and Mutation!", Sys.time())

  deltasAnnotated <- data.frame( as.data.frame(probePositions), deltas, "MUTATIONS" = mutation)
  
  deltasAnnotatedSorted <- sortByCHRandSTART(deltasAnnotated)

  deltasAnnotatedSorted <- subset(deltasAnnotatedSorted, MUTATIONS == 1)[,c("CHR","START","END","DELTA")]
  
  dumpSampleAsBedFile(
    dataToDump = deltasAnnotatedSorted,
    fileExtension =  ".DELTAS.METHYLATION.bedgraph",
    resultFolder =  resultFolder,
    resultSubFolder = "DELTAS_METHYLATION" ,
    sampleName = sampleName,
    multipleFileColNames = c("CHR","START","END","DELTA","SAMPLENAME")
  )
}


