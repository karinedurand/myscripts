library(ape)

ptm <- proc.time()

## importe un fasta
dna <- read.dna(file = "/home/kadurand/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/Tip-Dating/DATA/alignementavecR.fasta", format="fasta", as.character = "TRUE", as.matrix="TRUE")
## find SNPS
# pour chaque colonne, renvoi 1 si il y a au moins un A et 0 si il n'y en a pas - idem pour t,c,g.
isa=as.numeric(apply(dna=="a"|dna=="A",2,sum)>0)
ist=as.numeric(apply(dna=="t"|dna=="T",2,sum)>0)
isc=as.numeric(apply(dna=="c"|dna=="C",2,sum)>0)
isg=as.numeric(apply(dna=="g"|dna=="G",2,sum)>0)
# somme. Si la valeur est>1, cela veut dire qu'il y a un SNP
issnp=isa+ist+isc+isg
# save snp pos
pos_var=(1:length(issnp))[issnp>1]
snp=dna[,issnp>1]
dim(snp)
# write fasta
#write.dna(snp,"snp.fa",format="fasta",nbcol=-1,colsep="")
for (i in 1:nrow(snp)){
  new=snp[i,]
  cat(paste(">",rownames(snp)[i],sep="",collapse=""), "\n", file="/home/kadurand/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/Tip-Dating/DATA/snppaucamulti.fa", append=T)
  cat(paste(new,sep="",collapse=""), "\n", file="/home/kadurand/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/Tip-Dating/DATA/snppaucamulti.fa", append=T)
}
proc.time() - ptm


## read the SNP file and filter out R sites

SNPfile <- read.dna(file = "/home/kadurand/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/Tip-Dating/DATA/snppaucamulti_AVECGAP.fa", format="fasta", as.character = "TRUE", as.matrix="TRUE")
## find and exclude recombiniting positions
#faire la somme par colonne (2,sum) des gap
countR=apply(SNPfile =="-"|SNPfile =="-",2,sum)
pos_recomb=(1:length(countR))[countR>=1]
pos_non_recomb=(1:length(countR))[countR<1]
SNP_no_recombi=SNPfile[,pos_non_recomb]
## write fasta
for (i in 1:nrow(SNP_no_recombi)){
  new=SNP_no_recombi[i,]
  cat(paste(">",rownames(SNP_no_recombi)[i],sep="",collapse=""), "\n", file="/home/kadurand/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/Tip-Dating/DATA/SNPnorecomb.fasta", append=T)
  cat(paste(new,sep="",collapse=""), "\n", file="/home/kadurand/partage_windows/Xylella/analyses_genomiques/ABC/myscripts/Tip-Dating/DATA/SNPnorecomb.fasta", append=T)
}
