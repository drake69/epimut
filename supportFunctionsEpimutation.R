### rangeBetaValuePerProbeAsNormalizedDistribution ############################################################################################################################################################################

rangeBetaValuePerProbeAsNormalizedDistribution <- function(populationMatrix, iqrTimes = 1) {
  
  library(parallel)
  library(doParallel)
  
  computationCluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE) - 1) 
  registerDoParallel(computationCluster)
  
  populationMatrixDim <- dim(populationMatrix)
  betaValues <- populationMatrix[, 2:populationMatrixDim[2]]
  row.names(betaValues) <- populationMatrix[,1]

  betaQ1Values <- parApply(computationCluster, betaValues, 1, quantile, 0.25)  
  betaQ3Values <- parApply(computationCluster, betaValues, 1, quantile, 0.75)  

  betaMedianValues <- parApply(computationCluster, betaValues, 1, median)
  betaValuesIQR <- parApply(computationCluster, betaValues, 1, IQR)

  if (!test_match_order(row.names(betaValues), row.names(betaMedianValues)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  
  if (!test_match_order(row.names(betaValues), row.names(betaValuesIQR)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  
  betaInferiorThresholds <- (betaQ1Values - (iqrTimes * betaValuesIQR))
  row.names(betaInferiorThresholds) <- row.names(betaMedianValues)
  
  result <-
    list(
      betaInferiorThresholds = betaInferiorThresholds,
      betaSuperiorThresholds = (betaQ3Values + (iqrTimes * betaValuesIQR)),
      betaMedianValues = betaMedianValues
    )

  stopCluster(computationCluster)
  return(result)
}

### rangeBetaValuePerGeneOverGeneAsNormalizedDistribution ############################################################################################################################################################################

rangeBetaValuePerGeneOverGeneAsNormalizedDistribution <- function(populationMatrix, iqrTimes = 1, probeFeatures) {
  
  library(parallel)
  library(doParallel)
  # library(fitdistrplus)
  
  computationCluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE) - 1) 
  registerDoParallel(computationCluster)
  
  betaValues <- populationMatrix[,-(1)]

  betaValuesMergedComplete <- data.frame("CHR" = probeFeatures$CHR, "gene" = probeFeatures$gene, "START" = probeFeatures$START, "PROBE" = probeFeatures$Probe, betaValues)
  betaValuesMergedNoGene <- subset(betaValuesMergedComplete, is.na(betaValuesMergedComplete$gene) | is.null(betaValuesMergedComplete$gene) | betaValuesMergedComplete$gene == "")

  ## by gene median #
  betaValuesMergedGene <- subset(betaValuesMergedComplete, (!is.na(betaValuesMergedComplete$gene) & !is.null(betaValuesMergedComplete$gene) & betaValuesMergedComplete$gene != ""))
  

  # browser()
  betaValuesNoGeneMatrix <- betaValuesMergedNoGene[,-(1:4)]
  
  betaMedianValuesNoGene <- parApply(computationCluster, betaValuesNoGeneMatrix, 1, median)
  betaMedianValuesNoGene <- data.frame(betaMedianValuesNoGene)
  colnames(betaMedianValuesNoGene) <- "median"
  
  betaIQRValuesNoGene <- parApply(computationCluster, betaValuesNoGeneMatrix, 1, IQR)
  betaIQRValuesNoGene <- data.frame(betaIQRValuesNoGene)
  colnames(betaIQRValuesNoGene) <- "iqr"
  
  # browser()
  
  genes <- as.character(unique(attributes(betaValuesMergedGene$gene)$levels))
  genes <- subset( genes, genes != "")
  betaMedianIQRValues <- foreach(i = 1:length(genes), .combine = rbind) %dopar%
  # for (i in 1:length(genes))  
  {
    geneSelected <- genes[i]
    betaValuesGene <- subset(betaValuesMergedGene, gene == geneSelected)
    probes <- as.character(betaValuesGene[,"PROBE"])
    betaValuesGeneAsMatrix <- as.matrix(betaValuesGene[,-(1:4)])
    
    medianValues <- rep(median(betaValuesGeneAsMatrix), dim(betaValuesGene)[1])
    iqrValues <- rep(IQR(betaValuesGeneAsMatrix), dim(betaValuesGene)[1])
    data.frame(row.names = probes, "median" = medianValues, "iqr" = iqrValues)
  }
  
  # browser()

  betaMedianValues <- data.frame(row.names = rownames(betaMedianIQRValues), betaMedianIQRValues[,"median"])
  colnames(betaMedianValues) <- "median"

  betaMedianValues <- rbind(betaMedianValues, betaMedianValuesNoGene)
  betaMedianValues <- betaMedianValues[
    order(row.names(betaMedianValues), decreasing = FALSE),
    ]
  if (!test_match_order(row.names(betaValues), row.names(betaMedianValues)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())

  betaIQRValues <- data.frame(row.names = rownames(betaMedianIQRValues), betaMedianIQRValues[,"iqr"])
  colnames(betaIQRValues) <- "iqr"
  
  betaIQRValues <- rbind(betaIQRValues, betaIQRValuesNoGene)
  betaIQRValues <- betaIQRValues[
    order(row.names(betaIQRValues), decreasing = FALSE),
    ]
  if (!test_match_order(row.names(betaValues), row.names(betaIQRValues)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  
  # browser()
  
  betaInferiorThresholds <- (betaMedianValues - (iqrTimes * betaIQRValues))
  row.names(betaInferiorThresholds) <- row.names(betaMedianValues)
  
  result <-
    list(
      betaInferiorThresholds = betaInferiorThresholds,
      betaSuperiorThresholds = (betaMedianValues + (iqrTimes * betaIQRValues)),
      betaMedianValues = betaMedianValues
    )
  
  stopCluster(computationCluster)
  return(result)
}

### rangeBetaValuePerCHROverCHRAsNormalizedDistribution ############################################################################################################################################################################

rangeBetaValuePerCHROverCHRAsNormalizedDistribution <- function(populationMatrix, iqrTimes = 1, probeFeatures) {
  
  library(parallel)
  library(doParallel)
  # library(fitdistrplus)
  # library(gdata)
  # library(e1071)
  
  computationCluster <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE) - 1) 
  registerDoParallel(computationCluster)
  
  betaValues <- populationMatrix[,-(1)]
  
  betaValuesMergedComplete <- data.frame("gene" = probeFeatures$gene, "CHR" = probeFeatures$CHR, "START" = probeFeatures$START, "PROBE" = probeFeatures$Probe, betaValues)
  betaValuesMergedNoCHR <- subset(betaValuesMergedComplete, is.na(betaValuesMergedComplete$CHR) | is.null(betaValuesMergedComplete$CHR) | betaValuesMergedComplete$CHR == "")
  
  ## by CHR median #
  betaValuesMergedCHR <- subset(betaValuesMergedComplete, (!is.na(betaValuesMergedComplete$CHR) & !is.null(betaValuesMergedComplete$CHR) & betaValuesMergedComplete$CHR != ""))
  
  # browser()
  betaValuesNoCHRMatrix <- betaValuesMergedNoCHR[,-(1:4)]
  
  betaMedianValuesNoCHR <- parApply(computationCluster, betaValuesNoCHRMatrix, 1, median)
  betaMedianValuesNoCHR <- data.frame(betaMedianValuesNoCHR)
  colnames(betaMedianValuesNoCHR) <- "median"
  
  betaIQRValuesNoCHR <- parApply(computationCluster, betaValuesNoCHRMatrix, 1, IQR)
  betaIQRValuesNoCHR <- data.frame(betaIQRValuesNoCHR)
  colnames(betaIQRValuesNoCHR) <- "iqr"
  
  browser()
  
  chromosomes <- as.character(unique(attributes(betaValuesMergedCHR$CHR)$levels))
  chromosomes <- subset( chromosomes, chromosomes != "")
  # betaMedianIQRValues <- foreach(i = 1:length(chromosomes), .combine = rbind) %dopar%
    for (i in 1:length(chromosomes))
    {
      geneSelected <- chromosomes[i]
      betaValuesCHR <- subset(betaValuesMergedCHR, CHR == geneSelected)
      probes <- as.character(betaValuesCHR[,"PROBE"])
      betaValuesCHRAsMatrix <- as.matrix(betaValuesCHR[,-(1:4)], ncol = 1)
      betaValuesCHRAsMatrix <- unmatrix(betaValuesCHRAsMatrix)
      # print(e1071::skewness(betaValuesCHRAsMatrix))
      # hist(betaValuesCHRAsMatrix)
      # fw <- fitdistrplus::descdist(betaValuesCHRAsMatrix)
      # print(summary(fw))
      medianValues <- rep(median(betaValuesCHRAsMatrix), dim(betaValuesCHR)[1])
      iqrValues <- rep(IQR(betaValuesCHRAsMatrix), dim(betaValuesCHR)[1])
      data.frame(row.names = probes, "median" = medianValues, "iqr" = iqrValues)
    }
  
  browser()
  
  betaMedianValues <- data.frame(row.names = rownames(betaMedianIQRValues), betaMedianIQRValues[,"median"])
  colnames(betaMedianValues) <- "median"
  
  betaMedianValues <- rbind(betaMedianValues, betaMedianValuesNoCHR)
  betaMedianValues <- betaMedianValues[
    order(row.names(betaMedianValues), decreasing = FALSE),
    ]
  if (!test_match_order(row.names(betaValues), row.names(betaMedianValues)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  
  betaIQRValues <- data.frame(row.names = rownames(betaMedianIQRValues), betaMedianIQRValues[,"iqr"])
  colnames(betaIQRValues) <- "iqr"
  
  betaIQRValues <- rbind(betaIQRValues, betaIQRValuesNoCHR)
  betaIQRValues <- betaIQRValues[
    order(row.names(betaIQRValues), decreasing = FALSE),
    ]
  if (!test_match_order(row.names(betaValues), row.names(betaIQRValues)))
    stop("Wrong order matching Probes and Mutation!", Sys.time())
  
  # browser()
  
  betaInferiorThresholds <- (betaMedianValues - (iqrTimes * betaIQRValues))
  row.names(betaInferiorThresholds) <- row.names(betaMedianValues)
  
  result <-
    list(
      betaInferiorThresholds = betaInferiorThresholds,
      betaSuperiorThresholds = (betaMedianValues + (iqrTimes * betaIQRValues)),
      betaMedianValues = betaMedianValues
    )
  
  stopCluster(computationCluster)
  return(result)
}

### sortByCHRandStart ############################################################################################################################################################################

sortByCHRandSTART <- function(dataframe)
{
  library(gtools)
  dataframe$CHR <- droplevels(dataframe$CHR)
  
  dataframe$CHR <- factor(dataframe$CHR, levels = mixedsort(attributes(dataframe$CHR)$levels))
  dataframe <- dataframe[
       with(dataframe, order(CHR, START)),
     ]
  return(dataframe)
}

### sortPivot ############################################################################################################################################################################


sortPivot <- function(pivotTable)
{
  library(gtools)

  # browser()
  chr <- colnames(pivotTable)[2:length(colnames(pivotTable))]
  chr <- (mixedsort(chr))
  pivotTable <- data.frame(pivotTable[,1], pivotTable[,mixedsort(chr)])
  colnames(pivotTable)[1] <- "SAMPLENAME"
  pivotTable <- pivotTable[
    with(pivotTable, mixedsort(SAMPLENAME)),
    ]
  return(pivotTable)

}

### getLesions ############################################################################################################################################################################

getLesions <- function(mutationAnnotatedSorted,slidingWindowSize, sampleName, probePositions, bonferroniThreshold = 0.05)
{
  library(checkr)
  
  colnames_required <- c("CHR","START","END")
  check_colnames( probePositions, colnames_required, error = TRUE)
  check_colnames( mutationAnnotatedSorted, "MUTATIONS", error = TRUE)
  
  mutationAnnotatedSortedLocal <- sortByCHRandSTART(mutationAnnotatedSorted)
  probePositions <- sortByCHRandSTART(probePositions)

  if (!test_match_order(row.names(mutationAnnotatedSortedLocal), probePositions$Probe))
    stop("Wrong order matching Probes and Mutation!", Sys.time())

  mutationAnnotatedSortedLocal$CHR <- droplevels(mutationAnnotatedSortedLocal$CHR)
  
  # if (length(unique(attributes(mutationAnnotatedSortedLocal$CHR)$levels)) != 1)
  #   browser()
  
  # print(unique(attributes(mutationAnnotatedSortedLocal$CHR)$levels))
  
  if (length(mutationAnnotatedSortedLocal$CHR) <= slidingWindowSize )
    browser()

  missedWindowLength <- ((slidingWindowSize - 1) / 2)
  
  #browser()
  mutationAnnotatedSortedWindowed <- zoo::zoo(x =  mutationAnnotatedSortedLocal[,"MUTATIONS"])
  message(sampleName," ","Got mutationAnnotatedSortedWindowed ", Sys.time() )
  
  sumIntoCentreWindow <- function(x) {
    zoo::rollapply(x, width = slidingWindowSize, by = 1, FUN = sum, align = "center")
  }
  
  ## check burden into the window (sum)
  mutationAnnotatedSortedWindowedSum <- sumIntoCentreWindow(mutationAnnotatedSortedWindowed)
  
  message(sampleName," ","Got mutationAnnotatedSortedWindowedSum ", Sys.time() )
  
  # x	vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
  # m the number of white balls in the urn.
  # n the number of black balls in the urn.
  # k the number of balls drawn from the urn.
  # p probability, it must be between 0 and 1.
  # dhyper(x, m, n, k, log = FALSE)
  
  # slidingWindowSizeVector <- rep(slidingWindowSize, length(mutationAnnotatedSortedWindowedSum))
  # 
  # if (length(mutationAnnotatedSortedWindowedSum) != length(sum(mutationAnnotatedSortedWindowedSum))
  #    | length(slidingWindowSizeVector) != length( (length(mutationAnnotatedSortedWindowed) - sum(mutationAnnotatedSortedWindowedSum))))
  #   browser()

  # white balls == mutated
  
  x <- mutationAnnotatedSortedWindowedSum # white balls drawm / mutated probes every drawn
  m <- sum(mutationAnnotatedSortedWindowedSum) # white balls total / mutated probes total
  n <- (length(mutationAnnotatedSortedWindowed) - sum(mutationAnnotatedSortedWindowed)) # black balls total /non mutated probes total
  k <- slidingWindowSize # number of balls drawn / number of probes drawn

  lesionpValue <- dhyper(x , m , n , k)

  #lesionpValue <- phyper(mutationAnnotatedSortedWindowedSum, sum(mutationAnnotatedSortedWindowedSum), (length(mutationAnnotatedSortedWindowed) - sum(mutationAnnotatedSortedWindowedSum)), slidingWindowSize)
  # lesionpValue <-   round( lesionpValue , 10)
  
  lesionpValue <- coredata(lesionpValue)
  
  startingValue <- lesionpValue[1:missedWindowLength]
  startingValue <- 1
  
  endingValue <- lesionpValue[(length(lesionpValue) - missedWindowLength):length(lesionpValue)]
  endingValue <- 1
  
  message(sampleName," ","Got lesionpValue ", Sys.time() )
  ### remove missed pValue from hypergeometric function
  lesionpValue[is.nan(lesionpValue)] <- 1
  message(sampleName," ","Replaces NaN  pValue in lesionpValue ", Sys.time() )
  
  ## correction by Bonferroni
  lesionWeighted <-   (lesionpValue < (bonferroniThreshold/length(lesionpValue)))
  
  missedValue <- rep(FALSE, (slidingWindowSize - 1) / 2)
  row.names(missedValue) <- row.names(startingValue)
  lesionWeighted <- append(missedValue, lesionWeighted)
  row.names(missedValue) <- row.names(endingValue)
  lesionWeighted <- append(lesionWeighted, missedValue)

  if (dim(probePositions)[1] != length(lesionWeighted))
    browser()
  
  lesionWeighted <- data.frame(as.data.frame(probePositions),"LESIONS" = lesionWeighted)
  
  lesionWeighted <- sortByCHRandSTART(lesionWeighted)
  lesionWeighted <- subset(lesionWeighted, LESIONS == 1)[,c("CHR","START","END")]
  
  if (dim(lesionWeighted)[1] > dim(mutationAnnotatedSortedLocal)[1])
    browser()  
  
  if (dim(lesionWeighted)[1] > dim(subset(mutationAnnotatedSortedLocal, MUTATIONS == 1))[1])
  {
    
    # ext <- function(i)
    # {
    #   # x	vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls.
    #   # m the number of white balls in the urn.
    #   # n the number of black balls in the urn.
    #   # k the number of balls drawn from the urn.
    # 
    #   d <- mutationAnnotatedSortedWindowedSum >= i
    #   
    #   x <- as.numeric(d == 1) # white balls drawm / mutated probes every drawn
    #   m <- sum(d == 1) # white balls total / mutated probes total
    #   n <- (length(d) - sum(d == 1)) # black balls total /non mutated probes total (all balls less the white balls)
    #   k <- 1 # number of balls drawn / number of probes drawn
    #   
    #   lesionpValue_1 <- dhyper(x , m , n , k)
    #   lesionpValue_1[is.nan(lesionpValue_1)] <- 1
    #   lesionWeighted_1 <-   (lesionpValue_1 < (0.05/length(lesionpValue_1)))
    #   print(sum(lesionWeighted_1))
    # }
    # 
    # browser()
    # ext(1)
    # ext(2)
    # ext(3)
    # ext(4)
    # ext(5)
    # ext(6)
    # ext(7)
    # ext(8)
    # ext(9)
  }
  
  return(lesionWeighted)
}


### dumpSampleAsBedFile ############################################################################################################################################################################

dumpSampleAsBedFile <- function(dataToDump, fileExtension, resultFolder, resultSubFolder, sampleName, multipleFileColNames) {
  
  library(plyr)
  library(stringi)
  
  if (resultSubFolder != "" && !dir.exists(paste0(resultFolder, "/", resultSubFolder, "/", sep = ""))) {
    dir.create(paste0(resultFolder, "/", resultSubFolder, "/", sep = ""))
  }
  
  if (!empty(dataToDump) && !startsWith(x = as.character(dataToDump[1,"CHR"]),prefix = "CHR"))
  {
    chr <- rep(x = "chr", dim(dataToDump)[1])
    chr <- paste0(chr, dataToDump[,"CHR"], sep = "")
    dataToDump[,"CHR"] <- chr
  }
  
  filePath <- paste0(resultFolder, "/", resultSubFolder, "/", sampleName , fileExtension, sep = "")
  write.table(dataToDump, file = filePath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  if (!empty(dataToDump))
  {
    message(sampleName," ","saved file ", filePath, " ", Sys.time())
    
    oldFilePath <- paste0(resultFolder, "/", resultSubFolder, "/", "MULTIPLE", fileExtension, sep = "")
    
    ## needed random extension to avoid concurrency problem to access same file during parallel operations
    tempExtension <- stri_rand_strings(1, 10)
    filePath <- paste0(resultFolder, "/", resultSubFolder, "/", "MULTIPLE", fileExtension, tempExtension, sep = "")
    sampleNames <- rep(sampleName, dim(dataToDump)[1])
    dataToDump <- data.frame(dataToDump, sampleNames)
    
    colnames(dataToDump) <- multipleFileColNames
    
    resultMultiple <- dataToDump
    if (!file.exists(filePath))
    {
      resultMultiple <- dataToDump
    }
    else
    {
      # oldFile <- read.table(file = filePath,sep = "\t", header = FALSE, quote = "")
      # colnames(oldFile) <- multipleFileColNames
      # resultMultiple <- rbind(oldFile, dataToDump)
    }
    
    write.table(
      resultMultiple,
      file = filePath,
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
    
    ## append to multiple file the single sample file though system, afterwards remove it!
    system(paste0("cat ", filePath, " >> ", oldFilePath, sep = ""))
    system(paste0("rm ", filePath, sep = ""))
    message(sampleName," ","saved files ", fileExtension, Sys.time())
  }
  
}

########  add cell with content ################################################################################################

addCellToDataFrame <- function(dataFrame, colSelection, cellValueSelection, colname, cellValue)
{
  #browser()
  if (!colname %in% colnames(dataFrame)) dataFrame[,colname] <- ""
  dataFrame[ dataFrame[,colSelection] == cellValueSelection ,colname] <- cellValue
  return(dataFrame)
}

##### createPivotResultFromMultipleBed ####################################################################################################################################################################################

createPivotResultFromMultipleBed <- function(resultFolder, anomalyLabel ,figureLable, probeFeatures)
{
  library(reshape2)
  library(plyr)
  
  # browser()

  souceFolder <- paste(resultFolder,"/", anomalyLabel,"_",  figureLable, "/", sep = "")
  fileName <- paste(souceFolder,"/", "MULTIPLE",".",figureLable,".", anomalyLabel,".bed", sep = "")
  sourceData <- read.table(fileName,sep = "\t", blank.lines.skip = TRUE, fill = TRUE)
  colnames(sourceData) <- c("CHR","START","END","SAMPLENAME")

  probeFeatures$CHR <- as.factor(paste0("chr", probeFeatures$CHR))
  probeFeatures <- probeFeatures[(probeFeatures$CHR %in% unique((sourceData$CHR))),]
  droplevels(probeFeatures$CHR)
  droplevels(sourceData$CHR)
  
  sourceDatabyCHR <- plyr::count(df = sourceData, vars = c("SAMPLENAME","CHR"))
  finalResult <- dcast(data = sourceDatabyCHR, SAMPLENAME ~ CHR, value.var = "freq")
  finalResult[is.na(finalResult)] <- 0
  finalResult[finalResult == ""] <- 0
  
  resultByCHR <- finalResult
  rm(finalResult)
  
  sourceData$CHR <- as.factor(sourceData$CHR)
  sourceData$START <- as.integer(sourceData$START)

  sourceData <- dplyr::left_join(sourceData, probeFeatures, by = c('CHR','START'))
  sourceData[is.na(sourceData)] <- ""
  sourceData$CHR <- paste0( sourceData$CHR, 
                            ifelse(
                              test = is.na(sourceData$gene) | is.null(sourceData$gene) | sourceData$gene == "", 
                              yes = paste0(":", sourceData$START), 
                              no = "_" 
                              ), 
                            sourceData$gene, 
                            sep = "")

  sourceData$CHR <- as.factor(sourceData$CHR)
  sourceData <- plyr::count(df = sourceData, vars = c("CHR","SAMPLENAME"))
  finalResult <- dcast(data = sourceData, CHR ~ SAMPLENAME, value.var = "freq")
  finalResult[is.na(finalResult)] <- 0
  finalResult[finalResult == ""] <- 0
  
  resultByGene <- finalResult
  
  result <- list("byCHR" = resultByCHR, "byGene" = resultByGene)
  
  return(result)  
  # destFolder <- paste( resultFolder,"/", sep = "")
  # fileName <- paste(destFolder,"/", "results.xls", sep = "")
  # sheetName <- paste(anomalyLabel,figureLable, sep = "_")ì
  # sheet <- list( sheetName = finalResult)
  # setNames(sheet, sheetName)
  # openxlsx::write.xlsx(x =  sheet, file =  fileName, asTable = TRUE )
  
}

############   createSummaryExcelFromCumulativeBedFile ################################################################################################
createSummaryExcelFromCumulativeBedFile <- function(resultFolder, probeFeaturesOriginal, sampleSheet)
{
  fileName <- paste0(resultFolder,"/", "results.xlsx", sep = "")
  
  hyperLesionsResult <- createPivotResultFromMultipleBed(resultFolder = resultFolder, anomalyLabel = "LESIONS", figureLable = "HYPER", probeFeatures =  probeFeaturesOriginal)
  hyperLesions <- hyperLesionsResult$byCHR
  hyperLesions[hyperLesions == 0] <- ""
  hyperLesions <- sortPivot(pivotTable =  hyperLesions)
  
  
  hypoLesionsResult <- createPivotResultFromMultipleBed(resultFolder = resultFolder, anomalyLabel = "LESIONS", figureLable = "HYPO", probeFeatures =  probeFeaturesOriginal)
  hypoLesions <- hypoLesionsResult$byCHR
  hypoLesions[hypoLesions == 0] <- ""
  hypoLesions <- sortPivot(hypoLesions)
  
  hyperMutationsResult <- createPivotResultFromMultipleBed(resultFolder = resultFolder, anomalyLabel = "MUTATIONS", figureLable = "HYPER", probeFeatures =  probeFeaturesOriginal)
  hyperMutations <- hyperMutationsResult$byCHR
  hyperMutations[hyperMutations == 0] <- ""
  hyperMutations <- sortPivot(hyperMutations)
  
  hypoMutationsResult <- createPivotResultFromMultipleBed(resultFolder = resultFolder, anomalyLabel = "MUTATIONS", figureLable = "HYPO", probeFeatures =  probeFeaturesOriginal)
  hypoMutations <- hypoMutationsResult$byCHR
  hypoMutations[hypoMutations == 0] <- ""
  hypoMutations <-  sortPivot(hypoMutations)
  
  # browser()
  
  # hypoLesions[is.numeric(hypoLesions)] <- hypoLesions[is.numeric(hypoLesions)] * -1
  # lesionsComparison <- dplyr::bind_rows(hyperLesions, hypoLesions)
  # lesionsComparison[is.na(lesionsComparison)] <- ""
  # lesionsComparison[lesionsComparison == 0] <- ""
  # lesionsComparison <- sortPivot(lesionsComparison)
  # 
  # hyperLesionsResult$byGene[is.numeric(hyperLesionsResult$byGene)] <- hyperLesionsResult$byGene[is.numeric(hyperLesionsResult$byGene)] * -1
  # lesionsComparisonByGene <- dplyr::bind_rows(hyperLesionsResult$byGene, hypoLesionsResult$byGene)
  # lesionsComparisonByGene[is.na(lesionsComparisonByGene)] <- ""
  # lesionsComparisonByGene[lesionsComparisonByGene == 0] <- ""
  
  hyperMutationsResultByGene <- hyperMutationsResult$byGene
  hyperMutationsResultByGene[hyperMutationsResultByGene == 0] <- ""
  hyperLesionsResult$byGene[hyperLesionsResult$byGene == 0] <- ""
  hypoMutationsResult$byGene[hypoMutationsResult$byGene == 0] <- ""
  hypoLesionsResult$byGene[hypoLesionsResult$byGene == 0 ] <- ""
  
  sheets <- list( "SUMMARY" = sampleSheet, 
                  # "LESION_COMPARISON" = lesionsComparison, 
                  "HYPER_MUTATIONS" = hyperMutations, 
                  "HYPER_LESIONS" = hyperLesions, 
                  "HYPO_MUTATIONS" = hypoMutations, 
                  "HYPO_LESIONS" = hypoLesions,
                  # "LESION_COMPARISON_ByGene" = lesionsComparisonByGene, 
                  "HYPER_MUTATIONS_ByGene" = hyperMutationsResultByGene, 
                  "HYPER_LESIONS_ByGene" = hyperLesionsResult$byGene, 
                  "HYPO_MUTATIONS_ByGene" = hypoMutationsResult$byGene, 
                  "HYPO_LESIONS_ByGene" = hypoLesionsResult$byGene
  )
  openxlsx::write.xlsx(x =  sheets, file =  fileName, asTable = TRUE )
  
}