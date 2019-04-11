#script to perform ACP on abc data

#source('00-scripts/rscript/cv4abc.R')
if("data.table" %in% rownames(installed.packages()) == FALSE) 
{install.packages("data.table", repos="https://cloud.r-project.org") 
    print("installing packages data.table ..." ) }
if("ade4" %in% rownames(installed.packages()) == FALSE) 
{install.packages("ade4", repos="https://cloud.r-project.org") 
    print("installing packages ade4 ..." ) }

library(data.table)
library(ade4)
library(readr)
setwd("~/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/")

target=as.numeric(read.table("00.data/OBS.ABC.stat.txt",skip=2,h=F))

#load simulations
IM1=fread("zcat 00.data/im.homom.homon.ABC.stat.txt.gz"    )
IM2=fread("zcat 00.data/im.heterom.homon.ABC.stat.txt.gz"  )
IM3=fread("zcat 00.data/im.homom.heteron.ABC.stat.txt.gz"  )
IM4=fread("zcat 00.data/im.heterom.heteron.ABC.stat.txt.gz")
AM1=fread("zcat 00.data/am.homom.homon.ABC.stat.txt.gz"    )
AM2=fread("zcat 00.data/am.heterom.homon.ABC.stat.txt.gz"  )
AM3=fread("zcat 00.data/am.homom.heteron.ABC.stat.txt.gz"  )
AM4=fread("zcat 00.data/am.heterom.heteron.ABC.stat.txt.gz")
SC1=fread("zcat 00.data/sc.homom.homon.ABC.stat.txt.gz"    )
SC2=fread("zcat 00.data/sc.heterom.homon.ABC.stat.txt.gz"  )
SC3=fread("zcat 00.data/sc.homom.heteron.ABC.stat.txt.gz"  )
SC4=fread("zcat 00.data/sc.heterom.heteron.ABC.stat.txt.gz")
SI1=fread("zcat 00.data/si.homom.homon.ABC.stat.txt.gz"    )
SI3=fread("zcat 00.data/si.homom.heteron.ABC.stat.txt.gz"  )

#replace missing by means
f <- function(x){
    m <- mean(x, na.rm = TRUE)
    x[is.na(x)] <- m
    x
}

SI1 <- apply(SI1, 2, f)
SI3 <- apply(SI3, 2, f)
SC1 <- apply(SC1, 2, f)
SC2 <- apply(SC2, 2, f)
SC3 <- apply(SC3, 2, f)
SC4 <- apply(SC4, 2, f)
IM1 <- apply(IM1, 2, f)
IM2 <- apply(IM2, 2, f)
IM3 <- apply(IM3, 2, f)
IM4 <- apply(IM4, 2, f)
AM1 <- apply(AM1, 2, f)
AM2 <- apply(AM2, 2, f)
AM3 <- apply(AM3, 2, f)
AM4 <- apply(AM4, 2, f)

#min number of sim
nlinesFul=min(nrow(IM1), nrow(IM2), nrow(IM3), nrow(IM4),
              nrow(SI1), nrow(SI3),
              nrow(AM1), nrow(AM2), nrow(AM3), nrow(AM4),
              nrow(SC1), nrow(SC2), nrow(SC3), nrow(SC4))

#keep only usefull stats 

IM2b=IM2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
IM3b=IM3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
IM4b=IM4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC1b=SC1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC2b=SC2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC3b=SC3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SC4b=SC4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM1b=AM1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM2b=AM2[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM3b=AM3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
AM4b=AM4[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SI1b=SI1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
SI3b=SI3[c(1:nlinesFul),-c(1:3,12,13,18:25)]
IM1b=IM1[c(1:nlinesFul),-c(1:3,12,13,18:25)]

#observed stats
target=target[-c(1:3,12,13,18:25)]
obs=matrix(rep(target[1:19],1), byrow=T, nrow=1)

###IMPORTDATA#############################################################################################################
###############################################################DATA FASTSIMBAC#######################################
sumfastalnSI <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/fastSimBac_linux/fastsimbac_sur_alignement/statistiquefastSimbac_1712736bp/summary_fastSI", col_names = FALSE)
sumfastalnAM <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/fastSimBac_linux/fastsimbac_sur_alignement/statistiquefastSimbac_1712736bp/summary_fastAM", col_names = FALSE)
sumfastalnSC <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/fastSimBac_linux/fastsimbac_sur_alignement/statistiquefastSimbac_1712736bp/summary_fastSC", col_names = FALSE)
sumfastalnIM <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/fastSimBac_linux/fastsimbac_sur_alignement/statistiquefastSimbac_1712736bp/summary_fastIM", col_names = FALSE)
obsfast <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/fastSimBac_linux/fastsimbac_sur_alignement/statistiquefastSimbac_1712736bp/obs.txt",  col_names = FALSE)

###############################################################DATA MS#######################################

statmsSI <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/msms/statsmsms/priors9/SIms9", col_names = FALSE)
statmsAM <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/msms/statsmsms/priors9/AMms9",  col_names = FALSE)
statmsSC <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/msms/statsmsms/priors9/SCms9",  col_names = FALSE)
statmsIM <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/msms/statsmsms/priors9/IMms9", col_names = FALSE)

obs997 <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/msms/statsmsms/ABSstat997obs.txt",  col_names = FALSE)


View(summaryIM)
View(summarySC)
View(summarystatAM)
View(SUMMARYsi)
View(obs997)
###############################################################################################################################
### ENLEVER LES na#############################################################################################################
f <- function(x){
  m <- mean(x, na.rm = TRUE)
  x[is.na(x)] <- m
  x
}
###############################################################DATA FASTSIMBAC#######################################
sumfastalnSI<- apply(sumfastalnSI, 2, f)
sumfastalnSC<- apply(sumfastalnSC, 2, f)
sumfastalnAM<- apply(sumfastalnAM, 2, f)
sumfastalnIM<- apply(sumfastalnIM, 2, f)
###############################################################DATA MS#######################################
statmsAM <-apply(statmsAM , 2, f)
statmsIM<-apply(statmsIM , 2, f)
statmsSC<-apply(statmsSC , 2, f)
statmsSI<-apply(statmsSI , 2, f)

###############################################################DATA FASTSIMBAC#######################################
#VIRER DES COLONNES SANS RIEN
sumfastalnSI=sumfastalnSI[c(1:nlinesFul),c(1:39)]
sumfastalnSC=sumfastalnSC[c(1:nlinesFul),c(1:39)]
sumfastalnAM=sumfastalnAM[c(1:nlinesFul),c(1:39)]
sumfastalnIM=sumfastalnIM[c(1:nlinesFul),c(1:39)]
obsfast=obsfast[c(1:1),c(1:39)]


###############################################################DATA MS#######################################
###############################################################DATA FASTSIMBAC#######################################

#min number of sim
nlinesFul=min(nrow(summaryIM), nrow(summarySC), nrow(summarystatAM), nrow(SUMMARYsi))
###############################################################DATA FASTSIMBAC#######################################
nlinesFul=min(nrow(sumfastalnAM),nrow(sumfastalnSC),nrow(sumfastalnIM),nrow(sumfastalnSI))
###############################################################DATA ms#######################################
nlinesFul=min(nrow(statmsSI ),nrow(statmsSC ),nrow(statmsAM),nrow(statmsIM))
nlinesFul=min(nrow(statmsIM),nrow(statmsSI),nrow(statmsSC))

reduce = nlinesFul/1 #reduction factor to not keep all sims

x <- as.factor(c(rep("summaryIM",reduce),
        rep("summarySC",reduce),
       # rep("summarystatAM",reduce),
      #  rep("SUMMARYsi",reduce), 
        rep("obs",1))) 
###############################################################DATA FASTSIMBAC#######################################
x <- as.factor(c(rep("sumfastalnAM",reduce),
                 rep("sumfastalnIM",reduce),
                  rep("sumfastalnSC",reduce),
                  rep("sumfastalnSI",reduce), 
                 rep("obsfast",1))) 
View(x) 
###############################################################DATA ms#######################################
x_ms <- as.factor(c(rep("statmsIM",reduce), rep("statmsAM",reduce),rep("statmsSI",reduce),rep("statmsSC",reduce),
                 rep("obs997",1))) 
View(x_ms) 
#z=rbind(SI1b[1:reduce,],
#    IM4b[1:reduce,],
#    AM4b[1:reduce,],
#    SC4b[1:reduce,])
###############################################################DATA FASTSIMBAC#######################################
z=rbind(summaryIM[1:reduce,],
        summarySC[1:reduce,],
        summarystatAM[1:reduce,],
        SUMMARYsi[1:reduce,])
###############################################################DATA MS#######################################
z_ms=rbind(statmsIM[1:reduce,],
           statmsAM[1:reduce,],
           statmsSI[1:reduce,],
           statmsSC[1:reduce,]
        #, sumfastalnIM[1:reduce,],sumfastalnSC[1:reduce,],sumfastalnSI[1:reduce,]
        )

View(z_ms)
###############################################################DATA FASTSIMBAC#######################################
datafast=(as.data.frame(rbind(z,obsfast)))
###############################################################DATA MS#######################################
datams=(as.data.frame(rbind(z_ms,obs997)))
View(datams)
couleur <- c("red","blue","black","pink","orange")

acp1 <- dudi.pca(datafast,scannf=F,nf=2)
###############################################################DATA MS#######################################
acpms <- dudi.pca(datams,scannf=F,nf=2)
acpms <- dudi.pca(datams,scannf=F,nf=4)
matrice<-data.frame(acpms$li)
write.table(matrice,sep=" ", file="matrice")
matrice<- read.table("~/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/matrice.csv", sep=";", col.names = F,row.names=1)
ggplot(matrice, aes(x=V2,y=V3, shape=V1, color=V1))+geom_point()+scale_color_manual(values=c('grey','pink', 'red',"black","blue"))

ggplot(matrice, aes(x=V4,y=V5, shape=V1, color=V1))+geom_point()+scale_color_manual(values=c('blue','black', 'red',"pink","orange"))

pdf(file="acp_ade4_2.pdf")
s.class(acp1$li,x,col=couleur) # avec des ellipses de confiance
acpms$li->ESSAI
s.class(acpms$li,x_ms,col=couleur) # avec des ellipses de confiance
dev.off()

s.label(acpms$li,
        xax = 1,     # Dimension 1
        yax = 2) 
scatter(acpms,
        posieig = "none", # Cacher le scree plot
        clab.row = 0      # Caché l'annotation de slignes
)
s.class(acpms$li,x)

############################dISTRIBUTION STATS##############################################################"
SI=datams[c(2001:3000),c(1:73)]
IM=datams[c(1:1000),c(1:73)]

ggplot(SI, aes(x=X2)) + geom_histogram()+geom_vline(xintercept=134.5930, color="red")
ggplot(datams, aes(x=X4)) + geom_histogram()+geom_vline(xintercept=5.264530, color="red")
ggplot(IM, aes(x=X10)) + geom_histogram()+geom_vline(xintercept=0.002188940, color="red")#THETA_0
ggplot(datams, aes(x=X20)) + geom_histogram()+geom_vline(xintercept=106.63300, color="red")
ggplot(datams, aes(x=X26)) + geom_histogram()+geom_vline(xintercept=9.86335e-04, color="red")

############### Autre méthode ##########################################
#Autre methode
library(hexbin)
library(laeken)
library('grid')

PCA_stats  <- princomp(z, scale=T, center=T)   
#Observed data
obs=matrix(rep(target[1:19],1), byrow=T, nrow=1)
colnames(obs) <-c("V4","V5","V6","V7","V8","V9","V10",
    "V11","V14","V15","V16","V17",
     "V26","V27" ,"V28","V29","V30","V31","V32")
     
obs <- data.frame(obs)

PCA_target <- predict(PCA_stats, obs)
x1=as.factor(c(rep("SI",nlinesFul),
    rep("IM",nlinesFul),
    rep("AM",nlinesFul),
    rep("SC",nlinesFul)))

nmod=length(table(x1))
theindex=names(table(x1))
couleur <- c("red","blue","green","black","yellow")

pdf(file="PCA_3.pdf", width=11.7, height=8.3)
for (pci in 1:3){ #x_coords
  for (pcj in (pci+1):4){  # y_coords
    hbin<-hexbin(PCA_stats$scores[,pci], 
          PCA_stats$scores[,pcj],
          xbins=100,
          xlab=paste("PC",pci),
          ylab=paste("PC",pcj)) #,
          #col=couleur)
    pp<-plot(hbin,legend=FALSE, main="PCA on summary statistics")
    pushHexport(pp$plot.vp)
    grid.points(PCA_target[pci],PCA_target[pcj],pch="*",gp=gpar(cex=2,col="red"))
  }
}
dev.off ( which=dev.cur() )



#NE PAS RAJOUTER DE FACTEUR
scan(file = "~/partage_windows/Xylella/analyses_genomiques/ABC/msms/statsmsms/im.csv",what=character(),sep="\t" )->summaryIM$X1
scan(file = "~/partage_windows/Xylella/analyses_genomiques/ABC/msms/statsmsms/AM.csv",what=character(),sep="\t" )->summarystatAM$X1
scan(file = "~/partage_windows/Xylella/analyses_genomiques/ABC/msms/statsmsms/SC.csv",what=character(),sep="\t" )->summarySC$X1
scan(file = "~/partage_windows/Xylella/analyses_genomiques/ABC/msms/statsmsms/si.csv",what=character(),sep="\t" )->SUMMARYsi$X1
