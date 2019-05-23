function (name, reps = 20, writeTrees = T) 
{
  options(warn = -1)
  inFileName <- paste0(name, ".xml")
  inFile <- readLines(inFileName)
  ver2 <- grep(pattern = "version=\"2.", x = inFile, value = F)
  ver1 <- grep(pattern = "version=\"1.", x = inFile, value = F)
  ver <- length(ver1) + length(ver2)
  if (ver == 1) {
    matchLines <- grep(pattern = "date value=", x = inFile, 
                       value = T)
    if (length(matchLines) == 0) {
      stop("No dates found, check BEAST input files")
    }
    matchLinesPosition <- grep(pattern = "date value=", x = inFile, 
                               value = F)
    matchFileName <- grep(pattern = "fileName", x = inFile, 
                          value = T)
    matchFileNamePosition <- grep(pattern = "fileName", x = inFile, 
                                  value = F)
    for (i in 1:reps) {
      newFile <- inFile
      randLines <- sample(matchLines)
      newFile[matchLinesPosition] <- randLines
      log = paste0("\\.log")
      matchLog <- grep(pattern = log, x = inFile, value = T)
      matchLogPosition <- grep(pattern = log, x = inFile, 
                               value = F)
      logRep <- paste0("\\.Rep", i, log)
      if (length(matchLogPosition) != 0) {
        newFile[matchLogPosition] <- gsub(log, logRep, 
                                          matchLog)
      }
      trees = paste0("\\.trees")
      matchTrees <- grep(pattern = trees, x = inFile, value = T)
      matchTreesPosition <- grep(pattern = trees, x = inFile, 
                                 value = F)
      treesRep <- paste0("\\.Rep", i, trees)
      if (length(matchTreesPosition) != 0) {
        newFile[matchTreesPosition] <- gsub(trees, treesRep, 
                                            matchTrees)
      }
      csv = paste0("\\.csv")
      matchCsv <- grep(pattern = csv, x = inFile, value = T)
      matchCsvPosition <- grep(pattern = csv, x = inFile, 
                               value = F)
      csvRep <- paste0("\\.Rep", i, csv)
      if (length(matchCsvPosition) != 0) {
        newFile[matchCsvPosition] <- gsub(csv, csvRep, 
                                          matchCsv)
      }
      ops = paste0("\\.ops")
      matchOps <- grep(pattern = ops, x = inFile, value = T)
      matchOpsPosition <- grep(pattern = ops, x = inFile, 
                               value = F)
      opsRep <- paste0("\\.Rep", i, ops)
      if (length(matchOpsPosition) != 0) {
        newFile[matchOpsPosition] <- gsub(ops, opsRep, 
                                          matchOps)
      }
      if (writeTrees == F) {
        stop("writeTrees function is not longer available")
      }
      out <- paste0(name, ".Rep", i, ".xml")
      cat(newFile, file = out, sep = "\n")
    }
    cat("Replicates done:", i, "\n")
  }
 ######################"
  
  inFile <- readLines("/home/kadurand/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/Tip-Dating/TipDatingBeast/13paucamulti_Inv")
   if (ver == 2) {
    numberTaxa <- length(grep("taxon=", inFile))
    line <- grep(pattern = "traitname=\"date|traitname='date", 
                 x = inFile)
    line <- line 
    if (length(line) == 0) {
      stop("No date info found, check BEAST input file")
    }
    datePositions = c()
    repeat {
      if (length(grep("value=", inFile[line])) > 0) 
        line <- line 
      if (length(grep("alignment", inFile[line])) > 0) 
        break
      if (length(grep("=", inFile[line])) > 0) {
        datePositions <- c(datePositions, line)
      }
      line <- line 
    }
    numberDates <- numberTaxa
    dateLines <- inFile[datePositions]
    dateLines <- trimws(dateLines)
    date1 <- unlist(strsplit(dateLines, ","))
    date <- unlist(strsplit(date1, "="))
    dateHap <- date[c(T, F)]
    length(dateHap)->dateHaptot
    dateHaptot1<-dateHaptot  + 1
    dateHap <- dateHap[(dateHaptot1-numberDates) : dateHaptot]
    dateValues <- date[c(F, T)]
    dateValues <-dateValues[(dateHaptot1-numberDates) : dateHaptot]
    lastLine <- length(grep("<taxa", dateValues))
    if (lastLine == 1) {
      lastDate <- tail(dateValues, 2)
      lastDate <- unlist(strsplit(lastDate, ""))
      lastDate <- head(lastDate, 1)
      dateValues <- head(dateValues, numberTaxa - 1)
      dateValues <- c(dateValues, lastDate)
    }
    dateValues <- gsub(",$", "", dateValues)
    
    
    dateValues<- c("2009","1983","2009","2004","2012","2002","2009","1999","2016","1997","2006","2015","1992","2015","2010","2015","2015","2006","2009","2003","1997","2009","2015","2013","2015","2010","2015","2015","2009","1997","2010","2009","1994","1996","2015","2015")
    for (i in 1:reps) {
      newFile <- inFile
      dateValues <- sample(dateValues)
      newDate <- paste0("\t\t\t", dateHap, "=", dateValues)
      newFile[datePositions] <- paste0(newDate, ",")
      datePositions[numberDates]
      if (lastLine == 1) {
        newFile[(datePositions[numberDates])] <- paste0(newDate[numberDates], 
                                                        "\t\t\t\t<taxa id=", date[numberDates * 2 + 
                                                                                    1], "=", date[numberDates * 2 + 2])
      }
      log = paste0("\\.log")
      matchLog <- grep(pattern = log, x = inFile, value = T)
      matchLogPosition <- grep(pattern = log, x = inFile, 
                               value = F)
      logRep <- paste0("\\.Rep", i, log)
      if (length(matchLogPosition) != 0) {
        newFile[matchLogPosition] <- gsub(log, logRep, 
                                          matchLog)
      }
      trees = paste0("\\.trees")
      matchTrees <- grep(pattern = trees, x = inFile, value = T)
      matchTreesPosition <- grep(pattern = trees, x = inFile, 
                                 value = F)
      treesRep <- paste0("\\.Rep", i, trees)
      if (length(matchTreesPosition) != 0) {
        newFile[matchTreesPosition] <- gsub(trees, treesRep, 
                                            matchTrees)
      }
      csv = paste0("\\.csv")
      matchCsv <- grep(pattern = csv, x = inFile, value = T)
      matchCsvPosition <- grep(pattern = csv, x = inFile, 
                               value = F)
      csvRep <- paste0("\\.Rep", i, csv)
      if (length(matchCsvPosition) != 0) {
        newFile[matchCsvPosition] <- gsub(csv, csvRep, 
                                          matchCsv)
      }
      if (writeTrees == F) {
        stop("writeTrees function is not longer available")
      }
      out <- paste0(name, ".Rep", i, ".xml")
      cat(newFile, file = out, sep = "\n")
    }
    cat("Replicates done:", i, "\n")
  }
  if (ver != 1 & ver != 2) {
    stop("Error, check BEAST input file")
  }
}