
####### extract info from "tableau pos recombi"....
Dir <- "/home/dmerda/Bureau/Documents/analyses_genomique/analyse_X.arboricola_sens_strict/extract_recombinant_tract/extraction_mauve"
setwd(Dir)
info <- scan("position_groupeA.txt",sep="\n",what="raw")
trim_names <- function(X,sep="_",pos=1,fix=TRUE){strsplit(X,split=sep,fixed=fix)[[1]][pos]}

cluster <- gsub(">cluster","",as.character(sapply(info,trim_names,pos=1,sep="|",fix=TRUE)))
start <- gsub(":.*","",as.character(sapply(info,trim_names,pos=2,sep="|",fix=TRUE)))
stop <- gsub(".*:","",as.character(sapply(info,trim_names,pos=2,sep="|",fix=TRUE)))
strain <- gsub(" ","",gsub("strain","",as.character(sapply(info,trim_names,pos=3,sep="|",fix=TRUE))))

info2 <- data.frame(cluster=as.numeric(cluster),start=as.numeric(start),stop=as.numeric(stop),strain=as.character(strain))

id <- as.character(unique(info2[,4]))

#### split sequence

test <- scan("align_clonalframeML.xmfa",what="raw",sep="\n")
pos <- sort(c(grep("=",test),grep(">",test)))

# this fct writes a fasta file for each strain by concatenating each cluster....

for(i in 1:length(id)){
  posi <- grep(paste(">",id[i],sep=""),test)
  beg <- posi+1
  end <- sapply(beg,function(X,Y){Y[Y>X][1]},pos) -1
  write(paste(">",id[i],sep=""),paste(id[i],".fasta",sep=""),append=TRUE)
  for(j in 1:length(beg)){
    sub <- test[beg[j]:end[j]]
    write(sub,paste(id[i],".fasta",sep=""),append=TRUE)
  }
}

####### lecture + modif

# This fct writes a new fasta for each strain with "R" being recombining sites

write_fasta <- function(s,id,outfile){
  write(paste(">",id,sep=""),outfile)
  lseq <- length(s)
  nseq <- floor(lseq/80)
  if(nseq==0) vseq <- c(0,lseq)
  if(nseq>0)
  {
    if(floor(lseq/80)==lseq/80) vseq <- c(0,(1:nseq)*80)
    if(floor(lseq/80)!=lseq/80) vseq <- c(0,(1:nseq)*80,lseq)
  }
  for(j in 1:(length(vseq)-1)){
    write(paste(s[(vseq[j]+1):vseq[j+1]],collapse=""),outfile,append=TRUE)
  }
}

for(i in 1:length(id)){
  s <- unlist(strsplit(paste(scan(paste(id[i],".fasta",sep=""),sep="\n",what="raw",skip=1),collapse=""),"")[[1]])
  sub <- info2[info2[,4]==id[i],]
  for(j in 1:nrow(sub)){
    print(paste(i,j,sep="_"))
    a <- sub[j,2]
    b <- sub[j,3]
    s[a:b] <- "R"
  }
  write_fasta(s,id[i],paste(id[i],"_MODIFIED.fasta",sep=""))
}
