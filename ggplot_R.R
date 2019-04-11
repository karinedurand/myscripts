obs997 <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/obs/997_ortho13paucamulti_clean/997obs.txt",  col_names = T)
obs1169 <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/obs/1169Nrecomb_17multiplex6st53/obs1169.txt",  col_names = T)

obs1209 <- read_table2("~/partage_windows/Xylella/analyses_genomiques/ABC/obs/1209_PAUCA_nonrecomb/ABC1209obs.txt",  col_names = T)
rbind(obs997,obs1169,obs1209)->OBS
as.factor(c("CVCmulti","ST53multi","CVCST53"))->facteur
dataOBS <- data.frame(facteur, OBS[,2:42])
col.name<-names(dataOBS)
#####
as.data.frame(names(dataOBS))->mom.col
for ( i in 2:42) { 
  #y=mom.col[i]
  ggplot(dataOBS, aes(x=facteur, y=mom.col[i,] ))+ 
           geom_boxplot()
}

ggplot(dataOBS, aes(x=facteur, y=print(mom.col[2,]))) + 
  geom_boxplot()

print(mom.col[2,])

Ap<-ggplot(dataOBS) +  geom_boxplot(aes(x=facteur, y=bialsites_std))
Bp<-ggplot(dataOBS) + geom_boxplot(aes(x=facteur, y=bialsites_avg))
Cp<-ggplot(dataOBS) + geom_boxplot(aes(x=facteur, y=thetaA_avg))
Dp<-ggplot(dataOBS) + geom_boxplot(aes(x=facteur, y=thetaB_avg))
Ep<-ggplot(dataOBS) + geom_boxplot(aes(x=facteur, y=FST_avg))
Fp<-ggplot(dataOBS) + geom_boxplot(aes(x=facteur, y=netdivAB_avg))
Gp<-ggplot(dataOBS) + geom_boxplot(aes(x=facteur, y=Gmax_avg))
Hp<-ggplot(dataOBS) + geom_boxplot(aes(x=facteur, y=sf_avg))
grid.arrange(Ap, Bp,Cp,Dp,Ep,Fp,Gp,Hp, ncol=2, nrow = 4)
+
  geom_boxplot(aes(x=facteur, y=sf_avg))+
  geom_boxplot(aes(x=facteur, y=sf_std))
