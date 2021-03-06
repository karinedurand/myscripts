# R SCRIPT#1 to run a single ABC-RF model choice analysis starting from the reference table generated by the DIYABC
v2.1.0 package (.bin format) and from a single dataset of observed summary statistics
# This script includes some code to obtain some detailed numerical outputs and several graphical illustrations of analysis
results
# Date: 01/09/2016 version 1.0.
# Licence: GPL2
# Authors: Arnaud Estoup, Jean-Michel Marin, Julien Foucaud, Alex Dehne-Garcia, and Antoine Fraimout
# Note: the example below is based on a DIYABC reference table named "reftable.bin" including 4 scenarios (models)
#and a total of 40000 simulated datasets summarized with 130 summary statistics
# Preliminary steps:

# STEP1 - Install the R package abcrf (Pudlo P, Marin JM, Estoup A, Cornuet JM, Gautier M, Robert CP. 2016. Reliable
#ABC model choice via random forests. Bioinformatics. 32: 859-866.)
# Load and install automatically the package "abcrf" from the CRAN by clicking in the corresponding option of a R_Gui
#or R_Studio console (option "Packages" in RGui and "Intall Packages" in RStudio")
install.packages("abcrf")
# or if you got the abcrf source file then write: install.packages(“abcrf_1.1.tar.gz”, repos = NULL, type="source")
library(abcrf)
library("abcrf", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")

# STEP2 - Go into your working directory which may be the DIYABC project directory.
setwd("/pathway_for_my_working_directory/")
# At minimum the chosen working directory should contain the following DIYABC files:
# the reference table "reftable.bin" (generated with DIYABC v2.1.0), the file RNG_state_0000.bin and the datafile
#"statobs.txt" (the latter including the summary statistics of the observed dataset)
# Add in the working directory the DIYABC v2.1.0 core executable file which can be downloaded at
#http://www1.montpellier.inra.fr/CBGP/diyabc/ for different operating systems (cf. downloadable file diyabc_core-2.1.0-linWinOsXExe.zip).
# The name of the DIYABC 2.1 core executable files is:
# diyabc_core-2.1.0-win.exe for windows
# diyabc_core-2.1.0-OsX for mac
# diyabc_core-2.1.0-linux-i386 or diyabc_core-2.1.0-linux-x64 or diyabc_core-2.1.0-OldLinux-x64 (or the executable file
#obtained through your own compilation) for linux 28

# STEP 3 - Convert the reftable.bin file into a reftable.txt file (cf. reference table file with a txt format) that will be ready
#to use by the present R script
# To do that, open a terminal, go to your working directory and write the following instruction:
# In a linux OS write: ./ diyabc_core-2.1.0-linux-i386 -p ./ -x
# In a windows OS write: diyabc_core-2.1.0-win.exe -p ./ -x
# In a mac OS write: ./ diyabc_core-2.1.0-OsX -p ./ -x
# Definition of the key parameters of the RF analysis:
# Number of compared scenarios (models) in the reference table
nscenarios=4
# Number of summary statistics in the reference table
nSS=130
# Size of the reference table (i.e. number of simulated dataset) that will be used to do the RF analysis
Nref= 40000
# Number of trees in the random forest (default = 500)
ntrees_in_forest=500
# Reading of the DIYABC datafile "statobs.txt" (the latter including the summary statistics of the observed dataset)
# Note: Usually the file statobs.txt includes a single row corresponding to a single vector of observed summary statistics
# see the R script # 3 for the treatment of a statobs.txt file including several rows corresponding to several vectors of observed summary statistics
statobs=read.table("statobs.txt",header=TRUE)
# Procedure to read the file reftable.txt and store data in a correct way for RF analysis
# Note: The procedure is a bit complex as it allows dealing with DIYABC reftables in which the number of parameters may be different for the different scenarios
# T1 = starting time of the analysis (including reading the file reftable.txt)
T1<-Sys.time()
ncol=max(count.fields("reftable.txt",skip=2*nscenarios+3))
toto=read.table("reftable.txt",skip=2*nscenarios+3,header=FALSE,fill=TRUE,row.names=NULL,col.names=paste("V",1:
ncol,sep=""),nrows=Nref)
toto=toto[order(toto[,1]),]
rm(ncol)
import.reftable=function(nb.SS=nSS)
{
scenarios=unique(toto[,1])
sortie=data.frame()
param=list()
long=rep(0,length(scenarios))
for(scen in scenarios)
{
lignes.scen=which(toto[,1]==scen)
nb.col=ncol(toto[lignes.scen[1],])-sum(is.na(toto[lignes.scen[1], ]))
extrait=toto[lignes.scen,1:nb.col]
param0=extrait[, 2:(nb.col-nb.SS)]
long[scen]=ncol(param0)
extrait=extrait[,c(1, (nb.col - nb.SS + 1):nb.col)]
colnames(extrait)=1:(nb.SS+1)
colnames(param0)=1:(ncol(param0))
rownames(extrait)=rownames(param0)=NULL
sortie=rbind(sortie,extrait)
param[[scen]]=param0
}
max.param.dim=max(long)
for(scen in scenarios)
{
if(ncol(param[[scen]])<max.param.dim)
param[[scen]]=cbind(param[[scen]],matrix(NA,nrow=nrow(param[[scen]]),ncol=max.param.dim-ncol(param[[scen]])))
}
mes.param=data.frame()
for(scen in scenarios)
{
colnames(param[[scen]])=1:max.param.dim
rownames(param[[scen]])=NULL
mes.param=rbind(mes.param, as.matrix(param[[scen]]))
}
return(list(ss=sortie, param=mes.param))
}
tmp=import.reftable()
indi=sample(1:nrow(tmp$ss),Nref)
reftable=tmp$ss[indi,]
param=tmp$param[indi,]
# Further tuning of the structure of reftable
colnames(reftable) <- c("modindex", colnames(statobs))
modindex <- as.factor(reftable[1:Nref,1])
sumsta <- reftable[1:Nref,-1]


### WE CAN NOW START RANDOM FOREST COMPUTATIONS #########
# Computation related to the forest classification
mc.rf <- abcrf(modindex,sumsta,paral=TRUE,ntree=ntrees_in_forest)
## Computation to identify the best model for the observed dataset and estimate its posterior probability
# The function "predict" provides the ID of the selected (i.e. best) model (scenario), the votes for each compared scenarios (over ntrees_in_forest),
# knowing that the best scenario is the one with the highest vote value, and the value of the posterior probability of the best model
predict(mc.rf,statobs,paral=TRUE,ntree=ntrees_in_forest)
# Duration of the analysis (including reading the file reftable.txt)
T2<-Sys.time()
duration = difftime(T2, T1)
duration

### CODE TO VISUALIZE RESULTS AND VARIOUS GRAPHICAL ILLUSTRATIONS################
# Visualize numerical results for the global prior error rates and the matrix of confusion (i.e. prior error rates detailed for each scenario)
# Estimation from 10000 pseudo-observed datasets from the reference table
mc.rf
# Vizualize and save in pdf files several useful illustrative graphics
# Graphic providing prior error rates for forests with different number of trees (pdf file name = ""error_vs_ntree.pdf")
# e.g. Fig. 3 in Pudlo et al. 2016
err.rf <- err.abcrf(mc.rf)
x11()
plot(err.rf)
pdf(file="error_vs_ntree.pdf",h=18,w=18)
plot(err.rf)
dev.off()
# Two graphics/figures providing (i) LDA projections of the reference table for the different scenarios plus the observed
#dataset (cf. black star in the figure)
# (file named graph_lda.pdf) and (ii) the contributions of the 30 most important statistics to the RF (file named
#graph_varImpPlot.pdf)
# e.g. Fig. S6 and Fig. S7 in Pudlo et al. 2016
plot(mc.rf,statobs,pdf=TRUE,n.var=30)




#Script #2: this second R script allows to running successively several ABC-RF model choice analyses starting from the same DIYABC v2.1.0 reference table and this for single dataset of observed summary statistics
# R SCRIPT#2 to run successively SEVERAL ABC-RF model choice analyses starting from the SAME reference table generated by the DIYABC v2.1.0 package (.bin format) and from a single dataset of observed summary statistics
# Date: 01/09/2016 version 1.0.
# Licence: GPL2
# Authors: Arnaud Estoup, Jean-Michel Marin, Julien Foucaud, Alex Dehne-Garcia, and Antoine Fraimout
# Note: the example below is based on a DIYABC reference table named "reftable.bin" including 4 scenarios (models)
#and a total of 40000 simulated datasets summarized with 130 summary statistics
# One wants 3 successive ABC-RF model choice analyses from a given reference table: the first one with 20000
#simulated datasets, the second one with 30000 simulated datasets and third one with 40000 simulated datasets.
# Preliminary steps:
# STEP1 - Install the R package abcrf (Pudlo P, Marin JM, Estoup A, Cornuet JM, Gautier M, Robert CP. 2016. Reliable
#ABC model choice via random forests. Bioinformatics. 32: 859-866.)
# Load and install automatically the package "abcrf" from the CRAN by clicking in the corresponding option of a R_Gui
#or R_Studio console (option "Packages" in RGui and "Intall Packages" in RStudio")
install.packages("abcrf")
# or if you got the abcrf source files then write: install.packages(“abcrf_1.1.tar.gz”, repos = NULL, type="source")
library(abcrf)

# STEP2 - Go into your working directory which may be the DIYABC project directory.
setwd("/pathway_for_my_working_directory/")
# At minimum the chosen working directory should contain the following DIYABC files:
# the reference table "reftable.bin" (generated with DIYABC v2.1.0), the file RNG_state_0000.bin and the datafile
"statobs.txt" (the latter including the summary statistics of the observed dataset)
# Add in the working directory the DIYABC v2.1.0 core executable file which can be downloaded at
#http://www1.montpellier.inra.fr/CBGP/diyabc/ for different operating systems (cf. downloadable file diyabc_core-2.1.0-linWinOsXExe.zip).
# The name of the DIYABC 2.1 core executable files is:
# diyabc_core-2.1.0-win.exe for windows
# diyabc_core-2.1.0-OsX for mac
# diyabc_core-2.1.0-linux-i386 or diyabc_core-2.1.0-linux-x64 or diyabc_core-2.1.0-OldLinux-x64 (or the executable file
#obtained through your own compilation) for linux
# STEP 3 - Convert the reftable.bin file into a reftable.txt file (cf. reference table file with a txt format) that will be ready
#to use by the present R script
# To do that, open a terminal, go to your working directory and write the following instruction:
# In a linux OS write: ./ diyabc_core-2.1.0-linux-i386 -p ./ -x
# In a windows OS write: diyabc_core-2.1.0-win.exe -p ./ -x
# In a mac OS write: ./ diyabc_core-2.1.0-OsX -p ./ -x
# Definition of the key parameters of the RF analysis:
# Number of compared scenarios (models) in the reference table
nscenarios=4
# Number of summary statistics in the reference table
nSS=130
# Size of the reference table (i.e. number of simulated dataset) that will be used to do the successive RF analyses
Nrefmax=40000
# Vector defining the 3 successive ABC-RF model choice analyses from a given reference table a size Nrefmax: the third
#one with 20000 simulated datasets,
# the second one with 30000 simulated datasets and first one with 40000 simulated datasets.
Nref= c(20000,30000,40000)
# Number of trees in the random forest (default = 500)
ntrees_in_forest=500
# Name and structure of the output file that will summarize numerical results for all successive RF analyses
# Outputs for each analysis are presented as following: "Nref best_scenario posterior_probability vote_scen1 vote_scen2
#vote_scen3 vote_scen4 prior_error_rate
output_file_name <- "output_RF.txt"
output <- data.frame(matrix(ncol=nscenarios+4))
options(digits=4)

# Reading of the DIYABC datafile "statobs.txt" (the latter including the summary statistics of the observed dataset)
# Note: Usually the file statobs.txt includes a single row corresponding to a single vector of observed summary statistics
# see the R script # 3 for the treatment of a statobs.txt file including several rows corresponding to several vectors of
#observed summary statistics
statobs=read.table("statobs.txt",header=TRUE)
# Procedure to read the file reftable.txt and store data in a correct way for RF analyses
# Note: The procedure is a bit complex as it allows dealing with DIYABC reftables in which the number of parameters
#may be different for the different scenarios
# T1 = starting time of the analysis (including reading the file reftable.txt)
T1<-Sys.time()
ncol=max(count.fields("reftable.txt",skip=2*nscenarios+3))
toto=read.table("reftable.txt",skip=2*nscenarios+3,header=FALSE,fill=TRUE,row.names=NULL,col.names=paste("V",1:
ncol,sep=""),nrows=Nrefmax)
toto=toto[order(toto[,1]),]
rm(ncol)
import.reftable=function(nb.SS=nSS)
{
scenarios=unique(toto[,1])
sortie=data.frame()
param=list()
long=rep(0,length(scenarios))
for(scen in scenarios)
{
lignes.scen=which(toto[,1]==scen)
nb.col=ncol(toto[lignes.scen[1],])-sum(is.na(toto[lignes.scen[1], ]))
extrait=toto[lignes.scen,1:nb.col]
param0=extrait[, 2:(nb.col-nb.SS)]
long[scen]=ncol(param0)
extrait=extrait[,c(1, (nb.col - nb.SS + 1):nb.col)]
colnames(extrait)=1:(nb.SS+1)
colnames(param0)=1:(ncol(param0))
rownames(extrait)=rownames(param0)=NULL
sortie=rbind(sortie,extrait)
param[[scen]]=param0
}
max.param.dim=max(long)
for(scen in scenarios)
{
if(ncol(param[[scen]])<max.param.dim)
param[[scen]]=cbind(param[[scen]],matrix(NA,nrow=nrow(param[[scen]]),ncol=max.param.dim-ncol(param[[scen]])))
}
mes.param=data.frame()
for(scen in scenarios)
{
colnames(param[[scen]])=1:max.param.dim
rownames(param[[scen]])=NULL
mes.param=rbind(mes.param, as.matrix(param[[scen]]))
}
return(list(ss=sortie, param=mes.param))
}
tmp=import.reftable()
########## STARTING SERIAL ABCRF COMPUTATIONS #################
for (i in 1:length(Nref))
{
indi=sample(1:nrow(tmp$ss),Nref[i])
reftable=tmp$ss[indi,]
param=tmp$param[indi,]
colnames(reftable) <- c("modindex", colnames(statobs))
modindex <- as.factor(reftable[1:Nref[i],1])
sumsta <- reftable[1:Nref[i],-1]
mc.rf <- abcrf(modindex,sumsta,paral=TRUE,ntrees_in_forest)
pred.rf <- predict(mc.rf,statobs,paral=TRUE,ntrees_in_forest)
output <- rbind(output, c(Nref[i], as.numeric(pred.rf[1]),
round(as.numeric(pred.rf[length(pred.rf)]),4),
unlist(pred.rf[2:(length(pred.rf)-1)]),
round(as.numeric(mc.rf$prior.err),4)))
}
output <- output[-1,]
colnames(output) <- c("Nref", "Best_scenario", "Posterior_probability", rownames(as.data.frame(unlist(pred.rf[2:
(length(pred.rf)-1)]))), "Prior_error_rate")
write.table(output, output_file_name, quote=FALSE, row.names=FALSE)
# Duration of the analysis (including reading the file reftable.txt)
T2<-Sys.time()
duration = difftime(T2,T1)
duration
rm(list=ls())
# END OF COMPUTATION: Open and edit the output file (output_file_name as defined above) which contains the
#numerical results of the ABC RF analyses: for each analysis outputs are presented as following: Nref Statobs_row_ID
#Best_scnario Posterior_probability vote_scen1 vote_scen2 vote_scen3 vote_scen4 Prior_error_rate





#Script #3: this third R script allows to running successively several ABC-RF model choice analyses starting from
#the same DIYABC v2.1.0 reference table and this for a data file including several observed datasets.
# R SCRIPT#3 to run successively SEVERAL ABC-RF model choice analyses starting from the SAME reference table
#generated by the DIYABC v2.1.0 package (.bin format)
# and this for a data file including SEVERAL rows corresponding to several vectors of observed summary statistics.
# Date: 01/09/2016 version 1.0.
# Licence: GPL2
# Authors: Arnaud Estoup, Jean-Michel Marin, Julien Foucaud, Alex Dehne-Garcia, and Antoine Fraimout
# Note: the example below is based on a DIYABC reference table named "reftable.bin" including 4 scenarios (models)
#and a total of 40000 simulated datasets summarized with 130 summary statistics
# One wants 6 successive ABC-RF model choice analyses from a given reference table: the first one with 20000
#simulated datasets, the second one with 30000 simulated datasets and third one with 40000 simulated datasets;
# each of the three analyses are processed on 3 vectors of observed summary statistic (in the statobs.txt data file); hence a
#total of 6 RF analyses.
# Notice that the (observed) data file may include much larger number of vectors of observed summary statistic without
#increasing too much the RF computational cost.
# Preliminary steps:

# STEP1 - Install the R package abcrf (Pudlo P, Marin JM, Estoup A, Cornuet JM, Gautier M, Robert CP. 2016. Reliable
#ABC model choice via random forests. Bioinformatics. 32: 859-866.)
# Load and install automatically the package "abcrf" from the CRAN by clicking in the corresponding option of a R_Gui
#or R_Studio console (option "Packages" in RGui and "Intall Packages" in RStudio")
install.packages("abcrf")
# or if you got the abcrf sources file then write: install.packages(“abcrf_1.1.tar.gz”, repos = NULL, type="source")
library(abcrf)
# STEP2 - Go into your working directory which may be the DIYABC project directory.
setwd("/pathway_for_my_working_directory/")
# At minimum the chosen working directory should contain the following DIYABC files:
# the reference table "reftable.bin" (generated with DIYABC v2.1.0), the file RNG_state_0000.bin and the datafile
#"statobs.txt" (the latter including the summary statistics of the observed dataset)
# Add in the working directory the DIYABC v2.1.0 core executable file which can be downloaded at
#http://www1.montpellier.inra.fr/CBGP/diyabc/ for different operating systems (cf. downloadable file diyabc_core-2.1.0-linWinOsXExe.zip).
# The name of the DIYABC 2.1 core executable files is:
# diyabc_core-2.1.0-win.exe for windows
# diyabc_core-2.1.0-OsX for mac
# diyabc_core-2.1.0-linux-i386 or diyabc_core-2.1.0-linux-x64 or diyabc_core-2.1.0-OldLinux-x64 (or the executable file
#obtained through your own compilation) for linux
# STEP 3 - Convert the reftable.bin file into a reftable.txt file (cf. reference table file with a txt format) that will be ready
#to use by the present R script
# To do that, open a terminal, go to your working directory and write the following instruction:
# In a linux OS write: ./ diyabc_core-2.1.0-linux-i386 -p ./ -x
# In a windows OS write: diyabc_core-2.1.0-win.exe -p ./ -x
# In a mac OS write: ./ diyabc_core-2.1.0-OsX -p ./ -x
# Definition of the key parameters of the RF analysis:
# Number of compared scenarios (models) in the reference table
nscenarios=4
# Number of summary statistics in the reference table
nSS=130
# Size of the reference table (i.e. number of simulated dataset) that will be used to do the successive RF analyses
Nrefmax=40000
# Vector defining the 3 successive ABC-RF model choice analyses from a given reference table a size Nrefmax: the third
#one with 20000 simulated datasets,
# the second one with 30000 simulated datasets and first one with 40000 simulated datasets.
Nref= c(20000,30000,40000)
# Number of trees in the random forest (default = 500)
ntrees_in_forest=500
# Name and structure of the output file that will summarize numerical results for all successive RF analyses
# Outputs for each analysis are presented as following: Nref Statobs_row_ID Best_scenario Posterior_probability
#vote_scen1 vote_scen2 vote_scen3 vote_scen4 Prior_error_rate
output_file_name <- "output_RF.txt"
output <- data.frame(matrix(ncol=nscenarios+5))
options(digits=4)
# Reading of the DIYABC datafile "statobs.txt" (the latter including the summary statistics of the observed dataset(s))
# Here statobs.txt file includes 3 rows corresponding to 3 vectors of observed summary statistics that will be analysed serially.
# You can add a text identifier at the beginning of each vectors of observed summary statistics of the DIYABC datafile
#initially named "statobs.txt" and now renamed "statobs_multi.txt" (for instance "dataobs_case_x"): this additional text
#will not hinder the analyses
statobs=read.table("statobs_multi.txt",header=TRUE)
# Procedure to read the file reftable.txt and store data in a correct way for RF analyses
# Note: The procedure is a bit complex as it allows dealing with DIYABC reftables in which the number of parameters
#may be different for the different scenarios
# T1 = starting time of the analysis (including reading the file reftable.txt)
T1<-Sys.time()
ncol=max(count.fields("reftable.txt",skip=2*nscenarios+3))
toto=read.table("reftable.txt",skip=2*nscenarios+3,header=FALSE,fill=TRUE,row.names=NULL,col.names=paste("V",1:
ncol,sep=""),nrows=Nrefmax)
toto=toto[order(toto[,1]),]
rm(ncol)
import.reftable=function(nb.SS=nSS)
{
scenarios=unique(toto[,1])
sortie=data.frame()
param=list()
long=rep(0,length(scenarios))
for(scen in scenarios)
{
lignes.scen=which(toto[,1]==scen)
nb.col=ncol(toto[lignes.scen[1],])-sum(is.na(toto[lignes.scen[1], ]))
extrait=toto[lignes.scen,1:nb.col]
param0=extrait[, 2:(nb.col-nb.SS)]
long[scen]=ncol(param0)
extrait=extrait[,c(1, (nb.col - nb.SS + 1):nb.col)]
colnames(extrait)=1:(nb.SS+1)
colnames(param0)=1:(ncol(param0))
rownames(extrait)=rownames(param0)=NULL
sortie=rbind(sortie,extrait)
param[[scen]]=param0
}
max.param.dim=max(long)
for(scen in scenarios)
{
if(ncol(param[[scen]])<max.param.dim)
param[[scen]]=cbind(param[[scen]],matrix(NA,nrow=nrow(param[[scen]]),ncol=max.param.dim-ncol(param[[scen]])))
}
mes.param=data.frame()
for(scen in scenarios)
{
colnames(param[[scen]])=1:max.param.dim
rownames(param[[scen]])=NULL
mes.param=rbind(mes.param, as.matrix(param[[scen]]))
}
return(list(ss=sortie, param=mes.param))
}
tmp=import.reftable()
########## STARTING SERIAL ABCRF COMPUTATIONS #################
for (i in 1:length(Nref)) {
for (j in 1:dim(statobs)[1]) {
indi=sample(1:nrow(tmp$ss),Nref[i])
reftable=tmp$ss[indi,]
param=tmp$param[indi,]
colnames(reftable) <- c("modindex", colnames(statobs))
modindex <- as.factor(reftable[1:Nref[i],1])
sumsta <- reftable[1:Nref[i],-1]
mc.rf <- abcrf(modindex,sumsta,paral=TRUE,ntrees_in_forest)
pred.rf <- predict(mc.rf,statobs[j,],paral=TRUE,ntrees_in_forest)
output <- rbind(output, c(Nref[i], j, as.numeric(pred.rf[1]),
round(as.numeric(pred.rf[length(pred.rf)]),4),
unlist(pred.rf[2:(length(pred.rf)-1)]),
round(as.numeric(mc.rf$prior.err),4)))
}
}
output <- output[-1,]
colnames(output) <- c("Nref", "Statobs_row_ID", "Best_scenario", "Posterior_probability",
rownames(as.data.frame(unlist(pred.rf[2:(length(pred.rf)-1)]))), "Prior_error_rate")
write.table(output, output_file_name, quote=FALSE, row.names=FALSE)

# Duration of the analysis (including reading the file reftable.txt)
T2<-Sys.time()
duration = difftime(T2,T1)
duration
rm(list=ls())
# END OF COMPUTATION: Open and edit the output file (output_file_name as defined above) which contains the
#numerical results of the ABC RF analyses: for each analysis outputs are presented as following: Nref Statobs_row_ID
#Best_scenario Posterior_probability vote_scen1 vote_scen2 vote_scen3 vote_scen4 Prior_error_rate