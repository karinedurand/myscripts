#!/usr/bin/Rscript
source("/home/dmerda/extraction_tract_recombinant/CFoutputAnalysis_Functions.R") # mettre le chemin de CFoutputAnalysis_Functions.R
multiFasta <- readLines("/home/dmerda/extraction_tract_recombinant/extraction_mauve/align_clonalframeML.xmfa") # mettre le chemin de ton xmfa
alignResults <- construct.alignment(multiFasta)
alignment <- alignResults$alignment
gap <- alignResults$gap
blocksPos <- alignResults$blocksPos
outputName <- "/home/dmerda/extraction_tract_recombinant/extraction_mauve/result1ABmauve_111215" # mettre le chemin de l'output CF
results <- extract.info(output, gap = gap, thetaFixed = TRUE, burnin = 5e+04)
strainsNames <- results$strainsNames
consTree <- results$consTree
recEvents <- results$recEvents
subEvents <- results$subEvents
refSites <- results$refSites
chain <- results$chain
recomList <- construct.recombList(recEvents, refSites, blocksPos, thresholdLow = 0.5, thresholdHigh = 0.95)
recBlocks <- construct.recBlocks(recomList, strainsNames)
head(recBlocks)
write.table(recBlocks, file = "/home/dmerda/extraction_tract_recombinant/extraction_mauve/tableau_position.txt") # mettre le chemin pour le fichier de sortie + mettre un nom 
extract.tracts(alignment, recBlocks, consTree, "/home/dmerda/extraction_tract_recombinant/extraction_mauve/recombinant_tracts_mauve.fasta")  # mettre le chemin pour le fichier de sortie + mettre un nom 
