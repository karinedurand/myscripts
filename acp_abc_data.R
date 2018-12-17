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
IM1b=IM1[c(1:nlinesFul),-c(1:3,12,13,18:25)]
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

#observed stats
target=target[-c(1:3,12,13,18:25)]
obs=matrix(rep(target[1:19],1), byrow=T, nrow=1)

reduce = nlinesFul/500 #reduction factor to not keep all sims

x <- as.factor(c(rep("SI",reduce),
        rep("IM",reduce),
        rep("AM",reduce),
        rep("SC",reduce), 
        rep("obs",1))) 
        
z=rbind(SI1b[1:reduce,],
    IM4b[1:reduce,],
    AM4b[1:reduce,],
    SC4b[1:reduce,])

data=(as.data.frame(rbind(z,obs)))
couleur <- c("red","blue","green","black","yellow","darkred")

acp1 <- dudi.pca(data,scannf=F,nf=2)

pdf(file="acp_ade4_2.pdf")
s.class(acp1$li,x,col=couleur) # avec des ellipses de confiance
dev.off()


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
