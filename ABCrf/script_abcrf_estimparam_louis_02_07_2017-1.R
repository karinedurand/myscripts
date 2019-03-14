#Installing and loading the R package abcrf
#install.packages("abcrf") # To install the abcrf package (version 1.6)
library(abcrf) # To load the package.

#######Reading data: option 1 - using a .bin and .text files obtained using DIYABC
#We assume the reference table is recorded within the reftable.bin file and its corresponding header in the
#header.txt file. The function readRefTable is used to recover the data.
data <- readRefTable(filename = "reftableRF.bin", header = "headerRF.txt")
# data is a list containing the scenarios (or models) indices, the matrix with the
# parameters, the summary statistics and other informations.
# We are here interested in the simulated data of the scenario 1.
index1 <- data$scenarios == 1 # To store the model 1 indexes.
# We then create a data frame composed of the parameter of interest poi and the
# summary statistics of the scenario 1.
# PARAM = ra7
data.poi <- data.frame(poi = data$params[index1, "ra7"], sumsta = data$stats[index1, ])
# subdataset size (cf. number of simulations in the reference table) to process the analysis
subdataset <- 500
data.poi <- data.poi[1:subdataset,]
#data.poi <- data.frame(poi = data.poi$poi[1:subdataset],  sumsta = data.poi$sumsta[1:subdataset,] ) 

#######Training a random forest
#The random forest of the ABC-RF method is built thanks to the regAbcrf function, its principle arguments being a
#R formula and the corresponding data frame as training dataset. Additional arguments are available, e.g. the number
#of trees, the number of covariates randomly considered at each split. See the regAbcrf help for further details.
model.poi <- regAbcrf(formula = poi~., data = data.poi, paral = TRUE)
# The used formula means that we are interested in explaining the parameter poi thanks to
# all the remaining columns of data.poi (i.e. all the summary statistics).
# The paral argument determine if parallel computing will be activated or not.

######Graphical representations of to access the performance of the method
#The evolution of the out-of-bag mean squared error depending on the number of tree can be easily represented with
#the err.regAbcrf function (e.g. Figure 6 of the main text).
errorOOB <- err.regAbcrf(object = model.poi, training = data.poi, paral = TRUE)

#######The contribution of summary statistics in ABC-RF estimation for the parameter of interest can be retrieved with the
#plot function applied to an object resulting from regAbcrf.
plot(x = model.poi, n.var = 25)
# The contributions of the 25 most important summary statistics are represented
# (e.g. Figure S5).

#######Making predictions
#Finally, given a data frame obs.poi containing the summary statistics you want the predictions of posterior quantities
#of interest for a given dataset (usually the observed dataset). When using DIYABC, note that the summary statistics
#of the observed dataset are recorded in a file name statobs.txt. The predict method can be used for this purpose.
#The column names need to be the same than those in the summary statistics of data.poi.
# Reading the observed dataset with
obs.poi <- read.table("statobsRF.txt", skip=2)
# If obs.poi is not a dataframe or the column names do not match,
# you can use the following lines:
obs.poi <- as.data.frame(obs.poi)
colnames(obs.poi) <- colnames(data.poi[ ,-1])
# Prediction is complete by
pred.obsPoi <- predict(object = model.poi, obs = obs.poi, training = data.poi,
quantiles = c(0.025,0.975), paral = TRUE)
# The 2.5 and 97.5 order quantiles are computed by specifying quantiles = c(0.025,0.975).
# pred.obsPoi is a list containing predictions of interest.
# Posterior mean can be retrieved by
pred.obsPoi$expectation
# Posterior variance by
pred.obsPoi$variance
# Posterior quantiles by
pred.obsPoi$quantiles

