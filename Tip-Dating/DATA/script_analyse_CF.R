source("/home/dmerda/Bureau/Documents/analyses_genomique/analyse_X.arboricola_sens_strict/clonalFrame_AB/CFoutputAnalysis_Functions.R")
multiFasta <- readLines("/home/dmerda/Bureau/Documents/Nico_Xylella/coregenome_renomme.xmfa")
alignResults <- construct.alignment(multiFasta)
alignment <- alignResults$alignment
gap <- alignResults$gap
blocksPos <- alignResults$blocksPos
outputName <- "/home/dmerda/Bureau/Documents/Nico_Xylella/result3_CF_Xf30_260116"
output <- scan(file = outputName, what = "")
results <- extract.info(output, gap = gap, thetaFixed = FALSE, burnin = 1e+04)
strainsNames <- results$StrainsNames
consTree <- results$consTree
recEvents <- results$recEvents
subEvents <- results$subEvents
chain <- results$chain
par(mfrow = c(nvar(chain), 2))
plot(chain, ask = FALSE, auto.layout = FALSE)
par(mfrow = c(1, 1))
crosscorr.plot(chain)
