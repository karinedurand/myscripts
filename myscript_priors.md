# scripts pour générer les priors du modèle IMex

---
title: "Prior_production_IMex"
author: "karine Durand"
date: "9 novembre 2018"
output:
  pdf_document: default
  html_document: default
---

## Production des priors IMex


```{r}
library(stats)
library("KScorrect", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
######partie locus 
#variables locus
#-L=taille du gene
#-t=theta
#-r=rho
#-delta=taille du track recombinant
#boucle de 1000000 iterations(1000000 tirage demographique)
demo<-NULL
locus<-NULL
tbs<-NULL
#####TIRER un prior locus dans une distribution uniforme de bornes 
L<-round(runif(1368,100,2500))#bound_taille du gene[100-2500]
t<-runif(1368,0, 0.0003 )#bound_theta=[0-0.0003]bornes pour 13pauca_multiplex
r<-runif(1368,0,0.0006)##bound_rho=[0,0.0006]
delta<-round(runif(1368,10, 1000))#bound=[10-1000]
#print(L,t,r,delta)
m_locus=matrix(c(L,t,r,delta),ncol=4)
m_locus=as.data.frame(m_locus)

for (i in 1:10000){#tirage des priors demographiques
  #variables demographique modéle SI
  ##Param_demo (5) = Ts  N1, N2,  M12, M21
  Ts<-rlunif(1,100,10000000)#bound=[100-10000000]
  N1<-rlunif(1,100,100000)#bound=[100-100000]
  N2<-rlunif(1,100,100000)#bound=[100-100000]
  M12<-runif(1,0.01,6)#bound=[0.01-6]
  M21<-runif(1,0.01,6)#bound=[0.01-6]
  Te<-rlunif(1,10,10000000)#)#bound=[0-100]
  ex<-runif(1,0.0001,10000)#)#bound=[0-1]

  #print( Ts  N1, N2,  M12, M21)
m_demo=matrix(c(Ts,N1,N2,M12,M21,Te,ex),ncol=7)
m_demo=as.data.frame(m_demo)
locus<-cbind(m_locus,m_demo)

path <- "/home//kadurand/partage_windows/Xylella/analyses_genomiques/ABC/fastSimBac_linux/Priors_IMex_10000/locus" 
write.table(locus,file= paste(path,i, sep="-"),col.names=FALSE,row.names =FALSE)
}
```

