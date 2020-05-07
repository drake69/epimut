########## buildSampleSheetFromGEO ############################################################################################################################################

buildSampleSheetFromGEO <-
  function(GEOgse, workingFolder, downloadFiles = FALSE) {
    library(GEOquery)
    library(R.utils)
    
    gse <- getGEO(GEOgse, GSEMatrix = TRUE)
    
    phenoDataName <-
      paste(GEOgse, "_series_matrix.txt.gz", sep = "")
    
    samplesheet <- gse[[phenoDataName]]@phenoData@data
    
    final_samplesheet <-
      data.frame(
        "Sample_Name" = samplesheet$supplementary_file,
        "Sample_Well" = "",
        "Sample_Plate" = "",
        "Sample_Group" = "",
        "Pool_ID" = "",
        "Sentrix_ID" = samplesheet$supplementary_file,
        "Sentrix_Position" = samplesheet$supplementary_file,
        "Case_Ctrl" = "0"
      )
    
    tempData <- final_samplesheet$Sample_Name
    
    tempData <-
      gsub("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/", "", tempData)
    tempData <-  gsub("suppl/", "", tempData)
    tempData <-  gsub("_Grn.idat.gz", " ", tempData)
    tempData <-  gsub("/", "_", tempData)
    tempData <-  gsub(" ", "", tempData)
    tempData <-  noquote(strsplit(tempData, "_"))
    
    library(dplyr)
    
    tempDataItemLength <- length(tempData[[1]])
    if (tempDataItemLength == 5)
    {
      final_samplesheet$Sample_Name <-
        lapply(tempData, '[[', 2) %>% unlist()
      final_samplesheet$Sentrix_ID <-
        paste(lapply(tempData, '[[', 3), lapply(tempData, '[[', 4), sep = "_") %>% unlist()
      final_samplesheet$Sentrix_Position <-
        lapply(tempData, '[[', 5) %>% unlist()
    } else {
      final_samplesheet$Sample_Name <-
        lapply(tempData, '[[', 2) %>% unlist()
      final_samplesheet$Sentrix_ID <-
        lapply(tempData, '[[', 3) %>% unlist()
      final_samplesheet$Sentrix_Position <-
        lapply(tempData, '[[', 4) %>% unlist()
    }
    
    if ("disease state:ch1" %in% names(final_samplesheet))
    {
      final_samplesheet$Pool_ID <-
        samplesheet$`disease state:ch1` %>% unlist()
    }
    final_samplesheet$Sample_Well <-
      samplesheet$`gender:ch1` %>% unlist()
    final_samplesheet$Sample_Group <-
      samplesheet$`tissue:ch1` %>% unlist()
    
    write.table(
      final_samplesheet,
      paste(workingFolder, "/", "final_samplesheet.csv", sep = ""),
      sep = ",",
      row.names = FALSE,
      quote = FALSE
    )
    
    
    if (downloadFiles) {
      for (fileName in samplesheet$supplementary_file)
      {
        localFileName <-  noquote(strsplit(fileName, "/"))
        localFileNameDim <- length(localFileName[[1]])
        localFileName <- localFileName[[1]][localFileNameDim]
        localFileName <-
          paste(workingFolder, "/", localFileName, sep = "")
        
        download.file(
          url = fileName,
          destfile =  localFileName,
          method =  "internal",
          # method = "wininet",
          quiet = FALSE,
          mode = "wb",
          cacheOK = TRUE,
          extra = getOption("download.file.extra"),
          headers = NULL
        )
        
        GEOquery::gunzip(localFileName, overwrite = TRUE)
      }
      for (fileName in samplesheet$supplementary_file.1)
      {
        localFileName <-  noquote(strsplit(fileName, "/"))
        localFileNameDim <- length(localFileName[[1]])
        localFileName <- localFileName[[1]][localFileNameDim]
        localFileName <-
          paste(workingFolder, "/", localFileName, sep = "")
        
        download.file(
          url = fileName,
          destfile =  localFileName,
          method =  "internal",
          quiet = FALSE,
          mode = "wb",
          cacheOK = TRUE,
          extra = getOption("download.file.extra"),
          headers = NULL
        )
        GEOquery::gunzip(localFileName, overwrite = TRUE)
      }
    }
  }

########## buildSampleSheetFromFolder ############################################################################################################################################

buildSampleSheetFromFolder <-
  function(workingFolder,
           renameFileStartingWithNumber = FALSE) {
    file.names <- dir(workingFolder, pattern = ".idat")
    
    #browser()
    for (i in 1:length(file.names)) {
      if (renameFileStartingWithNumber)
      {
        # m <- regexpr("d+", file.names[i], perl=TRUE)
        from <- file.names[i]
        if (startsWith(from, "0") ||
            startsWith(from, "1") ||
            startsWith(from, "2") ||
            startsWith(from, "3") ||
            startsWith(from, "4") ||
            startsWith(from, "5") ||
            startsWith(from, "6") ||
            startsWith(from, "7") ||
            startsWith(from, "8") ||
            startsWith(from, "9"))
        {
          file.names[i] <- paste("x", file.names[i], sep = "")
          # to <- file.names[i]
          # file.rename(from, to)
        }
      }
      
      file.names[i] <- gsub("_Grn.idat", "", file.names[i])
      file.names[i] <- gsub("_Red.idat", "", file.names[i])
    }
    file.names = unique(file.names)
    
    library(dplyr)
    
    final_samplesheet <-
      data.frame(
        "Sample_Name" = file.names ,
        "Sample_Well" = "",
        "Sample_Plate" = "",
        "Sample_Group" = "",
        "Pool_ID" = "",
        "Sentrix_ID" = noquote(strsplit(file.names, "_")),
        "Sentrix_Position" = noquote(strsplit(file.names, "_")),
        "Case" = "0"
      )
    
    final_samplesheet$Sample_Name <-
      final_samplesheet$Sample_Name %>% unlist()
    final_samplesheet$Sentrix_ID <-
      lapply(final_samplesheet$Sentrix_ID, '[[', 1) %>% unlist()
    
    for (i in 1:length(file.names)) {
      file.names[[i]] <- noquote(strsplit(file.names[[i]], "_"))
      fileName <- file.names[[i]]
      fileName <- fileName[[1]][-c(1)]
      file.names[[i]] <- paste(fileName, collapse = "_")
    }
    
    final_samplesheet$Sentrix_Position <- file.names %>% unlist()
    
    write.table(
      final_samplesheet,
      paste(workingFolder, "/final_samplesheet.csv", sep = ""),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }


########## syncSampleSheetFromFolder ############################################################################################################################################

syncSampleSheetFromFolder <-
  function(workingFolder,
           renameFileStartingWithNumber = FALSE,
           sampleSheet,
           columnSeparator
           ) {
    file.names <- dir(workingFolder, pattern = ".idat")
    
    for (i in 1:length(file.names)) {
      if (renameFileStartingWithNumber)
      {
        # m <- regexpr("d+", file.names[i], perl=TRUE)
        from <- file.names[i]
        if (startsWith(from, "0") ||
            startsWith(from, "1") ||
            startsWith(from, "2") ||
            startsWith(from, "3") ||
            startsWith(from, "4") ||
            startsWith(from, "5") ||
            startsWith(from, "6") ||
            startsWith(from, "7") ||
            startsWith(from, "8") ||
            startsWith(from, "9"))
        {
          file.names[i] <- paste("x", file.names[i], sep = "")
          # to <- file.names[i]
          # file.rename(from, to)
        }
      }
      
      file.names[i] <- gsub("_Grn.idat", "", file.names[i])
      file.names[i] <- gsub("_Red.idat", "", file.names[i])
    }
    file.names = unique(file.names)
    
    library(dplyr)
    
    sample_sheet_to_sync <- data.frame()
    if (file.exists(sampleSheet))
      sample_sheet_to_sync <- read.table(file =  sampleSheet, header = TRUE, sep = columnSeparator, fill = TRUE )
    
    final_samplesheet <-
      data.frame(
        "Sample_Name" = file.names ,
        "Sample_Well" = "",
        "Sample_Plate" = "",
        "Sample_Group" = "",
        "Pool_ID" = "",
        "Sentrix_ID" = noquote(strsplit(file.names, "_")),
        "Sentrix_Position" = noquote(strsplit(file.names, "_"))
      )
    
    final_samplesheet$Sample_Name <-
      final_samplesheet$Sample_Name %>% unlist()
    final_samplesheet$Sentrix_ID <-
      lapply(final_samplesheet$Sentrix_ID, '[[', 1) %>% unlist()
    
    for (i in 1:length(file.names)) {
      file.names[[i]] <- noquote(strsplit(file.names[[i]], "_"))
      fileName <- file.names[[i]]
      fileName <- fileName[[1]][-c(1)]
      file.names[[i]] <- paste(fileName, collapse = "_")
    }
    
    final_samplesheet$Sentrix_Position <- file.names %>% unlist()
    
    #browser()
    for (i in 1:dim(final_samplesheet)[1])
    {
      if (dim(subset(sample_sheet_to_sync,Sentrix_ID == final_samplesheet$Sentrix_ID[i] & Sentrix_Position == final_samplesheet$Sentrix_Position[i] ))[1] == 0)
      {
        new_row <- sample_sheet_to_sync[nrow(sample_sheet_to_sync) + 1,] 
        new_row$Sentrix_ID <- final_samplesheet$Sentrix_ID[i]
        new_row$Sentrix_Position <- final_samplesheet$Sentrix_Position[i]
        new_row$Sample_Name <- final_samplesheet$Sample_Name[i]
        sample_sheet_to_sync <- rbind(sample_sheet_to_sync, new_row)
      }
    }

    # sample_sheet_to_sync$check <- ifelse(
    #   is.na(
    #     match(
    #       paste0(sample_sheet_to_sync$Sentrix_ID, sample_sheet_to_sync$Sentrix_Position),
    #       paste0(final_samplesheet$Sentrix_ID, final_samplesheet$Sentrix_Position)
    #       )
    #     ),
    #   "No",
    #   "Yes"
    #   )
    
    # test <- ifelse(
    #   is.na(
    #     match(
    #       paste0(final_samplesheet$Sentrix_ID, final_samplesheet$Sentrix_Position),
    #       paste0(sample_sheet_to_sync$Sentrix_ID, sample_sheet_to_sync$Sentrix_Position)
    #     )
    #   ),
    #   paste0(sample_sheet_to_sync$Sentrix_ID, sample_sheet_to_sync$Sentrix_Position),
    #   "Yes"
    # )
    # 
    # sample_sheet_to_sync <- subset(sample_sheet_to_sync,check == "No")

    write.table(
      sample_sheet_to_sync,
      paste(workingFolder, "/final_samplesheet.csv", sep = ""),
      sep = columnSeparator,
      row.names = FALSE,
      quote = FALSE
    )
    
  }

##### linkSampleFilesToExperimentFolder ################################################################################################################################################# 
linkSampleFilesToExperimentFolder <- function(sampleFolder, experimentFolder, sampleSheet)
{
  library(R.utils)

    
  if (experimentFolder != "" && !dir.exists(experimentFolder)) {
    dir.create(experimentFolder)
  }
  
  system(paste0("rm -rf ", experimentFolder, "/*", sep = ""))
  
  createLink(paste0(experimentFolder, "/", basename(sampleSheet), sep = ""), sampleSheet, overwrite = TRUE)
  
  sampleSheet <- paste0(experimentFolder, "/", basename(sampleSheet))
  
  sample_sheet <- read.table(file =  sampleSheet, header = TRUE, sep = ",", fill = TRUE )
  
  
  for (i in 1:dim(sample_sheet)[1])
  {
    linkNamePrefix <- paste0(experimentFolder,"/", sample_sheet$Sentrix_ID[i],"_", sample_sheet$Sentrix_Position[i], sep = "")
    sampleNamePrefix <- paste0(sampleFolder,"/", sample_sheet$Sentrix_ID[i],"_", sample_sheet$Sentrix_Position[i], sep = "")
    
    linkNameRed <- paste0(linkNamePrefix,"_Red.idat", sep = "")
    sampleNameRed <- paste0(sampleNamePrefix,"_Red.idat", sep = "")
    if (!file.exists(sampleNameRed))
    {
      message("Missed sample file ", sampleNameRed)
    }
    createLink( linkNameRed, sampleNameRed )
    
    linkNameGrn <- paste0(linkNamePrefix,"_Grn.idat", sep = "")
    sampleNameGrn <- paste0(sampleNamePrefix,"_Grn.idat", sep = "")
    if (!file.exists(sampleNameGrn))
    {
      message("Missed sample file ", sampleNameGrn)
    }
    createLink( linkNameGrn, sampleNameGrn )
  }  
}


##### createSampleFilesToExperimentFolder ################################################################################################################################################# 
createSampleFilesToExperimentFolder <- function(sampleFolder, experimentFolder, sampleSheet)
{
  library(R.utils)
  
  
  if (experimentFolder != "" && !dir.exists(experimentFolder)) {
    dir.create(experimentFolder)
  }
  
  system(paste0("rm -rf ", experimentFolder, "/*", sep = ""))
  
  file.copy( sampleSheet, paste0(experimentFolder, "/", basename(sampleSheet), sep = ""))
  
  sampleSheet <- paste0(experimentFolder, "/", basename(sampleSheet))
  
  sample_sheet <- read.table(file =  sampleSheet, header = TRUE, sep = ",", fill = TRUE )
  
  
  for (i in 1:dim(sample_sheet)[1])
  {
    destinationNamePrefix <- paste0(experimentFolder,"/", sample_sheet$Sentrix_ID[i],"_", sample_sheet$Sentrix_Position[i], sep = "")
    sourceNamePrefix <- paste0(sampleFolder,"/", sample_sheet$Sentrix_ID[i],"_", sample_sheet$Sentrix_Position[i], sep = "")
    
    destinationNameRed <- paste0(destinationNamePrefix,"_Red.idat", sep = "")
    sourceNameRed <- paste0(sourceNamePrefix,"_Red.idat", sep = "")
    if (!file.exists(sourceNameRed))
    {
      message("Missed sample file ", sourceNameRed)
    }
    file.copy(sourceNameRed, destinationNameRed )
    
    destinationNameGrn <- paste0(destinationNamePrefix,"_Grn.idat", sep = "")
    sourceNameGrn <- paste0(sourceNamePrefix,"_Grn.idat", sep = "")
    if (!file.exists(sourceNameGrn))
    {
      message("Missed sample file ", sourceNameGrn)
    }
    file.copy( sourceNameGrn, destinationNameGrn )
  }  
}



test_match_order <- function(x,y) {
  
  if (all(x == y)) return(TRUE)
    # print('Perfect match in same order')
  
  if (!all(x == y) && all(sort(x) == sort(y)))  return(FALSE)
    # print('Perfect match in wrong order')
  
  if (!all(x == y) && !all(sort(x) == sort(y)))  return(FALSE)
    # print('No match')
}
