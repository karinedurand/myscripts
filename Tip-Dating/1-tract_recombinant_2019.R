#!/usr/bin/Rscript
source("/home/kadurand/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/Tip-Dating/CFoutputAnalysis_Functions.R") # mettre le chemin de CFoutputAnalysis_Functions.R
multiFasta <- readLines("/home/kadurand/partage_windows/Xylella/analyses_genomiques/Tip-dating/clonalframe/parsnpcut.xmfa") # mettre le chemin de ton xmfa
alignResults <- construct.alignment(multiFasta)
alignment <- alignResults$alignment
gap <- alignResults$gap
blocksPos <- alignResults$blocksPos
outputName <- "/home/kadurand/partage_windows/Xylella/analyses_genomiques/Tip-dating/clonalframe/resultpaucamulti1" # mettre le chemin de l'output CF
output <- scan(file = outputName, what = "")
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
write.table(recBlocks, file = "/home/kadurand/partage_windows/Xylella/analyses_genomiques/Tip-dating/clonalframe/tableau_position.txt") # mettre le chemin pour le fichier de sortie + mettre un nom 
extract.tracts(alignment, recBlocks, consTree, "/home/kadurand/partage_windows/Xylella/analyses_genomiques/Tip-dating/clonalframe/recombinant_tracts.fasta")  # mettre le chemin pour le fichier de sortie + mettre un nom 
write.phylip(alignment, consTree)
plot(consTree)
calcul.thetaBranch(subEvents, consTree)
listOfOrthologousClusters <- clusters
