require(ape)
require(coda)
## library(plotrix)


construct.alignment <- function(multiFasta, ref.strain = NULL)
{	
    ## returns a matrix "alignement" row=strains col=aligned Sites
    ## arg= multiFasta (obtained with readLines(fichier XMFA.fasta))
    ## attention à strainRef!
    ## Save one fasta per bloc
    cut <- grep("=", multiFasta)
    block1 <- multiFasta[1:(cut[1]-1)]
    write(file = "block1.fasta", block1)
    for (i in 2:length(cut)){
        assign(paste("block",i,sep = ""), multiFasta[(cut[i-1]+1) : (cut[i]-1)])
        write(file = paste("block",i,".fasta",sep = ""), get(paste("block",i,sep = "")))
    }
    
    ## Read .fasta and create a matrix per block
    for (i in 1:length(cut)){
        assign(paste("block",i,"Fasta",sep = ""), as.character(read.dna(file = paste("block",i,".fasta",sep = ""), format = "fasta")))
        unlink(paste("block",i,".fasta",sep = ""))
    }
    
    ## Complete Alignment
    listBlocks <- paste0("block", 1:length(cut), "Fasta")
    blockLengths <- unlist(lapply(listBlocks, function(x) ncol(get(x))))
    blocksPos <- cbind(1 + c(0, cumsum(blockLengths[-length(blockLengths)])),
                       cumsum(blockLengths))
    alignment <- matrix(NA, nrow = nrow(block1Fasta), ncol = blocksPos[nrow(blocksPos), 2])
    for (i in 1:length(cut)) {
        fasta <- get(paste("block",i,"Fasta",sep = ""))
        fasta <- fasta[order(rownames(fasta)),]            
        alignment[ , blocksPos[i, 1]:blocksPos[i, 2]] <- fasta
    }
    rownames(alignment) <- sort(rownames(fasta))
    if (!is.null(ref.strain)) {
        if (!(ref.strain %in% rownames(fasta)) ) { stop("Invalid reference strain") }
    } else {
          ref.strain <- rownames(block1Fasta)[1]
      }
    index <- which(sort(rownames(fasta)) == ref.strain)
    alignment[c(1, index), ] <- alignment[c(index, 1), ]
    rownames(alignment)[c(1, index)] <- rownames(alignment)[c(index, 1)]
    if (alignment[1, 1] == "-") {
        warning(paste0("The default reference strain (", sort(rownames(fasta))[1], ") starts with a gap, this may lead to further problems, you can choose a different one using \"ref.strain = \"."))
    }
    gap <- grep("-", alignment[1, ])        
    rownames(blocksPos) <- sub("_tcoffee.fasta","",sub("# ", "", multiFasta[grep("#", multiFasta)]))
    return(alignResults=list(alignment=alignment, blocksPos=blocksPos, gap=gap))
}




construct.positionsMatrix <- function(alignment, multiFasta) {
	# on crée une matrice contenant pour chaque site de l'alignment sa position sur les génomes des différentes souches
	info <- multiFasta[grep("[>=]", multiFasta)]
	cut <- grep("=", info) 
	cut <- c(0,cut)
	
	# on recupere les positions des blocs (sur le pseudogenome) pour chaque souche 
	blocksArray <- array(NA, dim = c(dim(alignment)[1], 2, (length(cut)-1)))   # lignes=strains, col=start/end, z=block
	for (b in 1:(length(cut)-1)) {
		for (i in 1:dim(alignment)[1]){
		chaine <- info[(cut[b]+1):(cut[b+1]-1)][i]
		chaine <- strsplit(chaine, split=" ")[[1]][1]
		chaine <- gsub(paste(">", i, ":", sep=""), "", chaine)
		blocksArray[i,1,b] <- strsplit(chaine, split = "-")[[1]][1]
		blocksArray[i,2,b] <- strsplit(chaine, split = "-")[[1]][2]
		}
	}
	
	# creation de la positionsMatrix
	positionsMatrix <- matrix(NA, nrow=dim(alignment)[1], ncol=dim(alignment)[2])
	for (s in 1:dim(alignment)[1]) {
		block <- 1
		startBlock <- as.numeric(blocksArray[s,1,block]) 
		endBlock <- as.numeric(blocksArray[s,2,block])
		gap <- 0
		i <- 0
		for (p in 1:dim(alignment)[2]) {
			if (alignment[s,p] == "-") { # si on est dans un gap
				positionsMatrix[s,p] <- 0
				gap <- gap + 1
			} else {
				if (i< (endBlock - startBlock + 1)) { # si on a pas fini le block
					positionsMatrix[s,p] <- startBlock + i
					i <- i+1
				} else { # si on démarre un nouveau block
					i <- 1
					block <- block + 1
					startBlock <- as.numeric(blocksArray[s,1,block])
					endBlock <- as.numeric(blocksArray[s,2,block])
					positionsMatrix[s,p] <- startBlock
				}
			}	
		}
	}
	return(positionsMatrix)
}

polymorphic.sites <- function(alignment) {
	#returns a vector storing polymorphic sites positions of a given alignment 
	# arg= alignment (obtained with construct.alignment)
	polSites <- c()
	for (i in 1:dim(alignment)[2]) {
		alleles <- unique(alignment[,i])
		if (length(alleles[alleles != '-']) > 1) { #unique: donne un vecteur sans doublons
			polSites <- c(polSites,i)
		}
	}
	return(polSites)
}



triallelic.sites <- function(alignment) {
	# returns a vector "triallelicSites" storing positions of triallelic sites from an alignment (obtained with construct.alignment)
	# arg= alignment (obtained with construct.alignment)
	triallelicSites <- c()
	for (i in 1:dim(alignment)[2]) {
		if (length(unique(alignment[,i])) > 2) { #plus de 3 valeurs différentes
			triallelicSites <- c(triallelicSites,i)
		}
	}
	return(triallelicSites)
}

quadraallelic.sites <- function(alignment) {
	# returns a vector "triallelicSites" storing positions of triallelic sites from an alignment (obtained with construct.alignment)
	# arg= alignment (obtained with construct.alignment)
	quatreSites <- c()
	for (i in 1:dim(alignment)[2]) {
		if (length(unique(alignment[,i])) > 3) { #plus de 3 valeurs différentes
			quatreSites <- c(quatreSites,i)
		}
	}
	return(quatreSites)
}




sites.properties <- function(tree, alignment) {
	# returns a list of 3 vectors describing the sequence
	# arg = alignment (obtained with construct.alignment), tree
	# polSites = positions of polymorphic sites (missing data is counted as an allele)
	# gapSites = positions of gaps
	# observedState = for each nucleotide 0=non polymorph, 1=compatible (with the given tree), 2=non compatible, 3=missing data
	part <- prop.part(tree)
	
	for (i in 1:length(part)) {
		vect <- part[i]
		if (!any(is.element(part[[i]],1))) { #si la partition ne contient pas la souche 1
			part[[i]] <- setdiff(seq(length(tree$tip.label)),part[[i]]) #on prend le complementaire
		}
	}
	
	strainOrder <- match(tree$tip.label,rownames(alignment))
	orderedAlignment <- alignment[strainOrder,] #rearrangement de l'alignement pour que les souches soient dans le même ordre que dans l'arbre
	
	
	polSites <- c()
	gapSites <- c()
	observedState <- rep(NA, len=dim(alignment)[2]) #observedState 0=non polymorph, 1=compatible, 2=non compatible, 3=non renseigné
	
	for (p in 1:length(observedState)) {
		allels <- unique(orderedAlignment[,p]) # vect contenant les differents aucleotides du site (unique: donne un vecteur sans doublons)
			
		if (length(allels) == 1) { #site non polymorph
			observedState[p] <- 0
		} else {
			polSites <- c(polSites,p)
			if (any(is.element(allels,"-"))) {
				gapSites <- c(gapSites,p)
				observedState[p] <- 3
			} else { #si on est a un site polymorph sans gap
				if (length(allels) == 2) { #biallelic
					observedState[p] <- 2
					partition <- grep(allels[1], orderedAlignment[,p])
					if (length(partition)==1|length(partition)==(length(tree$tip.label)-1)) { # si le SNP est souche spécifique
						observedState[p] <- 1
					} else { #si le SNP n'est pas souche spé
						if (!any(is.element(partition,1))) { #si la partition ne contient pas la souche 1
							partition <- seq(length(tree$tip.label))[-partition] #on prend le complementaire
						} else {}
						for (i in 1:length(part)){
							if ( length(part[[i]]) == length(partition) && all(part[[i]]==partition) ) {
								observedState[p] <- 1
								break
							}
						}
					}
				} else { # on est dans un site polymorph sans gap mais pas biallelic
					vec <- rep(FALSE, len=length(allels))
					for (a in 1:length(allels)) { #pour chaque allele on verifie la compatibilité
						partition <- grep(allels[a], orderedAlignment[,p])
						if (length(partition)==1|length(partition)==(length(tree$tip.label)-1)) { # si le SNP est souche spécifique
							vec[a] <- TRUE
						} else {
							if (!any(is.element(partition,1))) { #si la partition ne contient pas la souche 1
								partition <- seq(length(tree$tip.label))[-partition] #on prend le complementaire
							} else {}
							for (i in 1:length(part)){
								if ( length(part[[i]]) == length(partition) && all(part[[i]]==partition) )  {
									vec[a] <- TRUE
									break
								} else {}
							}
						}
					}
					if (all(vec)) {
						observedState[p] <- 1
					} else { 
						observedState[p] <- 2
					}
				}
			}
		}
	}
	sitesProperties=list(polSites=polSites, gapSites=gapSites,observedState=observedState)
	return(sitesProperties)
}



extract.info <- function(output, gap, thetaFixed=FALSE, burnin=0, completeGenome=FALSE) {
	# returns a list of data extracted from clonal frame output
	# arg= output (scan (output clonal frame))
	
	# output organisation
	titlesNumbers <- list(constree = 1, names = 2, events = 3, consinfo = 4, mcmc = 5, phy = 6, poly = 7, ll = 8, blocks = 9, theta = 10, nu = 11, delta = 12, R = 13) #recuperation : titlesNumbers$R=13
	titlesPositions <- grep("#", output)

		#Recup Names
		namesMin <- titlesPositions[titlesNumbers$names]+1
		namesMax <- titlesPositions[titlesNumbers$names+1]-1
		strainsNames <- output[namesMin:namesMax]
		
		# Recup constree
		consTreeNewick <- output[2] # il manque un ; à la fin
		consTreeNewick <- paste(consTreeNewick, ";", sep="")
		consTree <- read.tree(text = consTreeNewick)
		consTree$tip.label <- strainsNames[as.numeric(consTree$tip.label)]
		
		#Recup blocks
		blocksMin <- titlesPositions[titlesNumbers$blocks]+1
		blocksMax <- titlesPositions[titlesNumbers$blocks+1]-1
		blocks <- output[blocksMin:blocksMax] #blocks = vector containing the number of reference sites in each blocks
		#Blocks Positions
		blocks <- cumsum(as.numeric(blocks))	
			
		#Recup poly (position of referenced sites)
		polyMin <- titlesPositions[titlesNumbers$poly]+1
		polyMax <- titlesPositions[titlesNumbers$poly+1]-1
		poly <- output[polyMin:polyMax]
		refSites <- as.integer(poly)
			#si genome complet
		if (completeGenome==TRUE) {
			positionMatrix <- get(paste("positionMatrixGC_",pop, sep=""))
			gap <- get(paste("gapGC_",pop, sep=""))
			blocksPos <- get(paste("blocksPosGC_",pop, sep=""))
			refSites[which(refSites!=0)] <- match(refSites[which(refSites!=0)],positionMatrix[1,])
			refSites[which(refSites==1510)] <- 0
			refSites <- modif.refSites.compGenome(refSites, blocks, gap, blocksPos)
		} else { refSites <- modif.refSites(refSites,blocks,gap) }

		#Recup events
		eventsMin <- titlesPositions[titlesNumbers$events]+1
		eventsMax <- titlesPositions[titlesNumbers$events+1]-1
		events <- output[eventsMin:eventsMax]
		#Events matrix construction size=N*2M N=nodes number, M=reference sites
		nodesNumber <- length(consTree$tip.label) + consTree$Nnode #Nnode=internal nodes only
		events <- matrix(as.numeric(events), nrow = nodesNumber, ncol = 2*length(poly), byrow = T)
		recEvents <- events[,2*(1:length(poly))-1]
		subEvents <- events[,2*(1:length(poly))]
		colnames(recEvents)<-as.character(refSites)
		colnames(subEvents)<-as.character(refSites)
		recEvents <- order.events(strainsNames, consTree, recEvents)
		subEvents <- order.events(strainsNames, consTree, subEvents)



		#Recup consinfo
		consinfoMin <- titlesPositions[titlesNumbers$consinfo]+1
		consinfoMax <- titlesPositions[titlesNumbers$consinfo+1]-1
		consinfo <- output[consinfoMin:consinfoMax]
		
		#Recup mcmc
		mcmcMin <- titlesPositions[titlesNumbers$mcmc]+1
		mcmcMax <- titlesPositions[titlesNumbers$mcmc+1]-1
		mcmc <- as.numeric(output[mcmcMin:mcmcMax])
		thining <- as.numeric(mcmc[3])
		#burnin <- as.numeric(mcmc[2])
		iterations <- as.numeric(mcmc[1])
		if (burnin >= iterations) {
                    stop("The requested burnin is higher than the number of iterations.")
                }
		
		#Recup phy
		phyMin <- titlesPositions[titlesNumbers$phy]+1
		phyMax <- titlesPositions[titlesNumbers$phy+1]-1
		phy <- output[phyMin:phyMax]
		
		#Recup ll
		llMin <- titlesPositions[titlesNumbers$ll]+1
		llMax <- titlesPositions[titlesNumbers$ll+1]-1
		ll <- as.numeric(output[llMin:llMax])

		
		#Recup theta
		thetaMin <- titlesPositions[titlesNumbers$theta]+1
		thetaMax <- titlesPositions[titlesNumbers$theta+1]-1
		theta <- as.numeric(output[thetaMin:thetaMax])
				
		#Recup nu
		nuMin <- titlesPositions[titlesNumbers$nu]+1
		nuMax <- titlesPositions[titlesNumbers$nu+1]-1
		nu <- as.numeric(output[nuMin:nuMax])
				
		#Recup delta
		deltaMin <- titlesPositions[titlesNumbers$delta]+1
		deltaMax <- titlesPositions[titlesNumbers$delta+1]-1
		delta <- as.numeric(output[deltaMin:deltaMax])
				
		#Recup R
		RMin <- titlesPositions[titlesNumbers$R]+1
		RMax <- titlesPositions[titlesNumbers$R+1]-1
		R <- as.numeric(output[RMin:RMax])
		
		
		#Construction of MCMC
		if (burnin==0){
			llmcmc <- mcmc(ll, thin=thining, end=iterations)
			Rmcmc <- mcmc(R, thin=thining, end=iterations)
			deltamcmc <- mcmc(delta, thin=thining, end=iterations)
			numcmc <- mcmc(nu, thin=thining, end=iterations)
			thetamcmc <- mcmc(theta, thin=thining, end=iterations)
		} else {
			llmcmc <- mcmc(ll[(burnin/thining):(iterations/thining)], thin=thining, start=burnin, end=iterations)
			Rmcmc <- mcmc(R[(burnin/thining):(iterations/thining)], thin=thining, start=burnin, end=iterations)
			deltamcmc <- mcmc(delta[(burnin/thining):(iterations/thining)], thin=thining, start=burnin, end=iterations)
			numcmc <- mcmc(nu[(burnin/thining):(iterations/thining)], thin=thining, start=burnin, end=iterations)
			thetamcmc <- mcmc(theta[(burnin/thining):(iterations/thining)], thin=thining, start=burnin, end=iterations)
		}
		

		if (thetaFixed==TRUE) {
			dataChain <- matrix(c(ll,nu,delta,R), ncol=4)
			colnames(dataChain) <- c("ll","nu","delta","R")
		} else {
			dataChain <- matrix(c(ll,theta,nu,delta,R), ncol=5)
			colnames(dataChain) <- c("ll","theta","nu","delta","R")
		}
		if (burnin==0){
			chain <- mcmc(data=dataChain, thin=thining, end=iterations)
		} else {
			chain <- mcmc(data=dataChain[(burnin/thining):(iterations/thining), ], thin=thining, start=burnin, end=iterations)
		}

	listResults <- list(strainsNames=strainsNames, consTree=consTree, refSites=refSites, recEvents=recEvents, subEvents=subEvents, consinfo=consinfo, mcmc=mcmc, phy=phy, llmcmc=llmcmc, blocks=blocks, thetamcmc=thetamcmc, numcmc=numcmc, deltamcmc=deltamcmc, Rmcmc=Rmcmc, chain=chain)
	return(listResults)
	
}


order.events <- function(strainsNames, consTree, matEvents) {
	## On réordonne les lignes de recEvents et subEvents dans le bon ordre
	nbStrains <- length(strainsNames)
	cfNodeNames <- c(strainsNames, as.character((nbStrains+1):(2*nbStrains-1))) ## noms dans Clonal Frame
	rNodeNames <- c(consTree$tip.label, consTree$node.label) ## noms (Clonal Frame) des noeuds (dans l'ordre de R)
	cf2rNodeOrder <- match(rNodeNames, cfNodeNames)
	matEvents <- matEvents[cf2rNodeOrder, ]
	return(matEvents)
}



modif.refSites <- function(refSites,blocks,gap) {
    ## return a vector storing the true positions of referenced sites 
    ## used in extract.info
    refSites <- 1+refSites
    refSitesNew<- rep(NA, len = length(refSites))
    it <- c(0,0)
    b <- 1
    first <- 1
    for (p in 1:length(refSites)) {
        if (p > blocks[b]) { #on change de bloc
            b <- b+1
            it[1] <- it[1]+1
            first <- p
        } else {
          }
        if (refSites[p] != 0) { # on est pas dans un gap
            refSitesNew[p] <- refSites[p] - it[1] + it[2]
        } else {	#on est dans un gap
              it[2] <- it[2]+1
              if (first != p) { #on n'est pas en premiere position du bloc
                  refSitesNew[p] <- refSitesNew[first] + gap[it[2]] - 1 #si gap=3 il faut rajouter que 2 par rapport au 1er
              } else{
                    refSitesNew[p] <- refSitesNew[p-1] + 1
                }	
          }
    }
    return(refSitesNew)
}


modif.refSites.compGenome <- function(refSites,blocks,gap, blocksPos) {
	# return a vector storing the true positions of referenced sites 
	# used in extract.info
	refSitesNew<- rep(NA, len = length(refSites))
	it <- 0
	b <- 1
	first <- 1
		for (p in 1:length(refSites)) {
			if (p > blocks[b]) { #on change de bloc
				b <- b+1
				first <- p
			} else {
			}
			if (refSites[p] != 0) { # on est pas dans un gap
				refSitesNew[p] <- refSites[p]
			} else {	#on est dans un gap
				it <- it+1
				if (first != p) { #on n'est pas en premiere position du bloc
					refSitesNew[p] <- blocksPos[b,1] + gap[it] - 1 #si gap=3 il faut rajouter que 2 par rapport au 1er
				} else{
					refSitesNew[p] <- blocksPos[b,1]
				}	
			}
		}
	return(refSitesNew)
}


blocks.position <- function(blocks,refSites) {
	# Blocks Positions (start, end)
	blocksPos <- matrix(NA, nrow = length(blocks), ncol = 2)
		blocksPos[,2] <- refSites[blocks]
		blocksPos[1,1] <- 1
		blocksPos[-c(1),1] <- refSites[blocks[-length(blocks)]]+1
	return(blocksPos)
}


HMM.pairSeq <- function(pairNames, alignment) {
	pairSeq <- alignment[pairNames,]
	otherNames <- rownames(alignment)[-match(pairNames,rownames(alignment))]
	otherSeq <- alignment[otherNames,]
	
	HMMpair <- rep(NA, len=dim(alignment)[2])
	#0=non polymorphic 1=unique polymorphism 2=repeated polymorphism 3=missing data
	
	for (p in 1:dim(alignment)[2]) {
		if (pairSeq[1,p] == "-" || pairSeq[2,p] == "-") {
			HMMpair[p] <- 3
		} else {
			if (pairSeq[1,p] ==  pairSeq[2,p]) {HMMpair[p] <- 0}
			else {
				vec <- c(FALSE, FALSE)
				other <- unique(otherSeq[,p]) #regarde les différents alleles présentés par les autres souches
				if (length(other) == 1 && other != "-") {HMMpair[p] <- 1}
				else {
					if  (all(is.element(c(pairSeq[1,p],pairSeq[2,p]), other))) {HMMpair[p] <- 2}
					else {HMMpair[p] <- 3}
				}
			}
		}
	}
	return(HMMpair)
}


write.HMM <- function(HMM, fastaFile, blocksPos) {
	for (cluster in rownames(blocksPos)) {
	seqBeg <- blocksPos[cluster, 1]
	seqEnd <- blocksPos[cluster, 2]
	clusSeq <- c("A", "G", "C", "N")[1+HMM[seqBeg:seqEnd]]
	if (cluster == rownames(blocksPos)[1]) { 
	write(file = fastaFile, paste(">", cluster), append = FALSE)
	} else {
	write(file = fastaFile, paste(">", cluster), append = TRUE)
	}
	write(file = fastaFile, clusSeq, ncolumns = 120, append = TRUE, sep = "")
	}
}

write.alignment <- function(alignment, fastaFile) {
	for (strain in rownames(alignment)) {
		if (strain == rownames(alignment)[1]) { 
		write(file = fastaFile, paste(">", strain), append = FALSE)
		} else {
		write(file = fastaFile, paste(">", strain), append = TRUE)
		}
		write(file = fastaFile, toupper(alignment[strain,], ncolumns = 120), append = TRUE, sep = "")
	}
}

write.xmfa <- function(alignment, fastaFile) {
	for (strain in rownames(alignment)) {
		if (strain == rownames(alignment)[1]) { 
		write(file = fastaFile, paste(">", strain), append = TRUE)
		} else {
		write(file = fastaFile, paste(">", strain), append = TRUE)
		}
		write(file = fastaFile, toupper(alignment[strain,]),ncolumns = 120, append = TRUE, sep = "")
	}
	write(file = fastaFile, "=", append=TRUE)
}

write.phylip <- function(alignment, phylipFile) {
	write(file = phylipFile, dim(alignment), append = FALSE)
	for (strain in rownames(alignment)) {
		str <- paste(strain, paste(toupper(alignment[strain,]), collapse="",sep=""), sep="\t")
		write(file = phylipFile, str, append = TRUE)
	}
}

cluster.checking <- function(alignment, blocksPos) {
	# detects the clusters that contain a stretch of close polymorphic sites 
	allels <- rep(0, len=dim(alignment)[2])
	for (p in 1:dim(alignment)[2]) {
		position <- unique(alignment[,p]) 
		if (length(position)>1 && any(is.element(position,"-"))) position <- position[-match("-",position)]
		if (length(position)>1) allels[p] <- 1
	}
	
	clusterToCheck <- c()	
	for (p in seq(from=1, to=dim(alignment)[2]-15)){
		sum <- sum(allels[p:(p+15)])
		if (sum > 10 ) {
			clusterToAdd <- intersect(which(blocksPos[,2] > p),which(blocksPos[,1] < p))
			clusterToAdd <- rownames(blocksPos)[clusterToAdd]
			clusterToCheck <- unique(c(clusterToCheck,clusterToAdd))
		}
	}
	return(clusterToCheck)
}


extract.sites <- function(matEvents, nodes = "all", threshold = 0.5) {
	if (nodes == "all") {probSiteEvents <- colSums(matEvents)} 
	else { 
		if (length(nodes) == 1) {probSiteEvents <- matEvents[nodes,]}
		else {probSiteEvents <- colSums(matEvents[nodes,])}
	}
	return(colnames(matEvents)[probSiteEvents > threshold])
}

## Vérifier la vitesse de la fonction
construct.recombList <- function(recEvents, refSites, blocksPos, thresholdLow=0.5, thresholdHigh=0.95) {
	# Create a list of n matrix (n=nb of nodes in the tree)
	# Matrix i contains the clusters that are affected by a recombination at the node i (and the position of the tract)
	recombList <- list()
        if (is.null(rownames(blocksPos))) {
            rownames(blocksPos) <- paste0("cluster", 1:nrow(blocksPos))
        }
	for (node in seq(dim(recEvents)[1])) {
		print(paste("parsing node", node, "out of", nrow(recEvents)))
		recSites <- extract.sites(recEvents, nodes=node, threshold=thresholdLow)
		start <- 0
		end <- 0
		b <- 1
		matNode <- matrix()
		for (p in recSites) {
			refPos <- which(refSites == p)
			if ((refPos == 1) || 
			!is.element(refSites[refPos-1],recSites)) {
				start <- p
			} else { 
				if ((refPos == length(refSites)) || 
				!is.element(refSites[refPos+1],recSites)) {
					end <- p
					## print(c(p, start, end))
					if (max(recEvents[node, which(refSites==start):which(refSites==end)])>thresholdHigh) {
						clusterToAdd <- ((blocksPos[,2] >= as.numeric(start)) &
										(blocksPos[,1] <= as.numeric(end)))
                                                clusterToAdd <- rownames(blocksPos)[clusterToAdd]
						matrixToAdd <- blocksPos[clusterToAdd, , drop = F]
						matrixToAdd[1, 1] <- start
						matrixToAdd[nrow(matrixToAdd) , 2] <- end
						## Allouer matNode dès le départ
						if (b==1) { matNode <- matrixToAdd ; b <- b+1 }
						else { matNode <- rbind(matNode,matrixToAdd) }
					}
				}
			}
		}
		recombList[[node]] <- matNode
	}
		return(recombList)
}



construct.recBlocks <- function(recombList, strainsNames) {
    ## returns a matrix containing information about recombined tracts detected
    ## ncol= nb of recombination events
    ## rows = cluster, start position, end position, tracts size, node
    recBlocksCF <- matrix(nrow=5)
    for (i in (1:length(recombList))[-(length(strainsNames)+1)]) {
        ## skip empty blocks
        if (!is.na(recombList[[i]][1,1])) {
            blocksToAdd <- rownames(recombList[[i]])
            m <- matrix(ncol=length(blocksToAdd), nrow=5)
            m[1,] <- blocksToAdd
            m[5,] <- i
            tailleTracts <- c()
            startPos <- c()
            endPos <- c()
            for (r in 1:dim(recombList[[i]])[1]) {
                tailleToAdd <- as.numeric(recombList[[i]][r,2])-as.numeric(recombList[[i]][r,1])
                tailleTracts <- c(tailleTracts, tailleToAdd)
                startPos <- c(startPos,as.numeric(recombList[[i]][r,1]))
                endPos <- c(endPos,as.numeric(recombList[[i]][r,2]))
            }
            m[2,] <- startPos
            m[3,] <- endPos
            m[4,] <- tailleTracts
            recBlocksCF <- cbind(recBlocksCF, m)
        }
    }
    recBlocksCF <- recBlocksCF[,2:dim(recBlocksCF)[2]]
    rownames(recBlocksCF) <- c("cluster","start position","end position","tracts size","node")
    recBlocksCF <- data.frame(cluster = recBlocksCF[1, ],
                              start = as.integer(recBlocksCF[2, ]),
                              end = as.integer(recBlocksCF[3, ]),
                              size = as.integer(recBlocksCF[4, ]),
                              node = as.integer(recBlocksCF[5, ])
                              )
    return(recBlocksCF)
}

extract.tracts <- function(alignment, recBlocks, consTree, file) {
    if (nrow(recBlocks) == 0) {
        stop("No recombinant blocks in recBlocks")
    }
    recBlocks <- subset(recBlocks, node <= nrow(alignment))
    recBlocks$strain <- consTree$tip.label[recBlocks$node]
    sequence <- function(i) {
       seq <- toupper(alignment[recBlocks$strain[i], recBlocks$start[i]:recBlocks$end[i]])
       return(paste(seq[seq != "-"], collapse = ""))
    }
    seqs <- unlist(lapply(1:nrow(recBlocks), sequence))
    labels <- paste0(">", recBlocks$cluster, "|", recBlocks$start, ":",
                     recBlocks$end, "|", "strain", recBlocks$strain)
    for (i in 1:nrow(recBlocks)) {
        if (i == 1) {
            cat(labels[i], "\n", file = file)
        } else {
            cat(labels[i], "\n", file = file, append = TRUE)
        }
        cat(seqs[i], "\n", file = file, append = TRUE)
    }
}


calcul.branchLength <- function(consTree) {
	nbStrains <- length(consTree$tip.label)
	nbNodes <- consTree$Nnode + nbStrains
	branchLength <- rep(NA, len =  nbNodes)
	for (n in seq(nbNodes)) {
		if (n != (nbStrains+1)) {
			pos <- which(consTree$edge[,2]==n)
			branchLength[n] <- consTree$edge.length[pos]
		} else {
			branchLength[n] <- 0
		}
	}
	return(branchLength)
}

node2edge <- function(phy, values, default = NA) {
    ntaxa <- length(phy$tip.label)
    ## pad values to have total length equal to number of nodes and tips
    ## by adding NA values (corresponding to tips) in front of it
    edge.values <- c(rep(default, ntaxa), values)[phy$edge[ , 2]]
    return(edge.values)
}

construct.number.recevents <- function(consTree, recBlocks) {
    nbin <- consTree$Nnode + length(consTree$tip.label)
    ## Number of recombination events by node, named by node names (ape format)
    table <- rep(0, length = nbin)
    names(table) <- c(consTree$tip.label, consTree$node.label)
    ## Convert node number (CF format) to node name (ape format)
    recBlocks$name <- as.character(recBlocks$node)
    recBlocks$name[recBlocks$node <= length(consTree$tip.label)] <- consTree$tip.label[recBlocks$node[recBlocks$node <= length(consTree$tip.label)]]
    counts <- tapply(recBlocks$node, recBlocks$name, length)
    table[names(counts)] <- counts
    ## Return number of recombination events on edges
    return(table[consTree$edge[ , 2]])
}


calcul.thetaBranch <- function(subEvents, consTree,
							 branchLength=calcul.branchLength(consTree)){
	#branchLength <- c()
	thetaBranch <- c()
	nbStrains <- length(consTree$tip.label)
	for (n in seq(consTree$Nnode + nbStrains)) {
		if (n != (nbStrains+1)) {
			pos <- which(consTree$edge[,2]==n)
			#branchLength[n] <- consTree$edge.length[pos]
			branchPol <- sum(subEvents[n,])
			thetaBranch[n] <- 2*branchPol/branchLength[n]
		} else {
			#branchLength[n] <- 0
			thetaBranch[n] <- NA
		}
	}
	return(thetaBranch)
}

calcul.thetaBoot <- function(alignment, bootstrap=20, thetaBoot=c()) {
    ## returns a vector containing x=bootstrap values of thetaWat calculated on
    ## bootstrap alignments
    for (i in 1:bootstrap) {
        if (i %% 10 == 0) {
            print(paste("Boostrap iteration", i))
        }
        sample <- sample(dim(alignment)[2], replace = TRUE)
        polSitesBoot <- polymorphic.sites(alignment[,sample])
        thetaBoot <- c(thetaBoot, length(polSitesBoot)/sum(1/(1:(dim(alignment)[1]-1))))
    }
    return(thetaBoot)
}


extract.clusterInfo <- function () {
	## Reads annotations of CDS
	tmp <- readLines("correspondance_file_CDS.txt")
	tmp <- sub("\t\t\t\t\t\t=\t", "|", tmp)
	tmp <- strsplit(tmp, "|", fixed = TRUE)
	tmp <- lapply(tmp, function(x) if (length(x) == 5) c(x, "") else  x)
	
	correspondanceOfCDS <- matrix(unlist(tmp), byrow = TRUE,
								nrow = length(tmp), ncol = 6)
	colnames(correspondanceOfCDS) <- c("gene_ID", "locus_tag", "gene_tag",
									"product", "gene_name", "EC_number")
	correspondanceOfCDS <- cbind(correspondanceOfCDS,
								sub("#[0-9]+", "", correspondanceOfCDS[, 3]))
	colnames(correspondanceOfCDS)[7] <- "strain"

	# Reads (orthologous) clusters content and store them in a matrix
	# clusterGeneContent (nrow=clusters ncol=strains)
	# each CDS is coded with its gene_ID, which is also its line number in
	# correspondanceOfCDS
	
	listOfOrthologousClusters <- clusters
	clusterGeneContent <- matrix(0, nrow = length(listOfOrthologousClusters), ncol = length(strainsNames))
	colnames(clusterGeneContent) <- strainsNames
	rownames(clusterGeneContent) <- sub("_tcoffee.fasta", "", listOfOrthologousClusters)
	
	for (cluster in listOfOrthologousClusters) {
		print(paste("Treating cluster", cluster))
		clus <- sub("_tcoffee", "", cluster)
		clus <- paste(clus, ".fasta", sep="")
		clusCont <- readLines(paste("/home/mig/cblin/10StrainsLacto_Delbru/All_Clusters/", clus, sep = ""))
		clusGeneContent <- clusCont[grep(">", clusCont)]
		clusGeneContent <- strsplit(clusGeneContent, "|", fixed = TRUE)
		index <- unlist(lapply(clusGeneContent, function(x) grep(x[2], correspondanceOfCDS[, "locus_tag"])))
		subsetOfIndexes <- match(strainsNames, gsub(" ", "", correspondanceOfCDS[index, "strain"]))
		clusterGeneContent[sub("_tcoffee.fasta", "", cluster),] <- index[subsetOfIndexes]
	}
	
	# Reads classifications (COG and agmial) and create a matrix storing this information row= cluster_X col= COG, agmial classification
	COG <- readLines("/home/mig/cblin/10StrainsLacto_Delbru/Classification/Ldb-ATCC11842-COGclassification.csv")
	agClass <- readLines("/home/mig/cblin/10StrainsLacto_Delbru/Classification/Ldb-ATCC11842-functionnalclassification.csv")

	clusterGeneFunction <- matrix(0, nrow = dim(clusterGeneContent)[1], ncol = 2)
	colnames(clusterGeneFunction) <- c("COG", "Agmial_Classification")
	rownames(clusterGeneFunction) <- rownames(clusterGeneContent)

	for (cluster in rownames(clusterGeneContent)) {
		#récuperation gene_tag ex LB#54
		id <- clusterGeneContent[cluster,"LB"]
		geneTag <- correspondanceOfCDS[which(correspondanceOfCDS[,1]==id),"gene_tag"]
		geneTag <- gsub(" ", "", geneTag)
		#recuperation classification COG
		COGgene <- COG[grep(paste(geneTag,";", sep=""), COG)]
		COGgene <- strsplit(COGgene, split="[", fixed=TRUE)[[1]][2]
		COGgene <- gsub("]", "",COGgene, fixed=TRUE) 
		COGgene <- gsub(" +", " ",COGgene) # retire les espaces multiples
		clusterGeneFunction[cluster,"COG"] <- COGgene
		#recuperation classification Agmial
		agClassGene <- agClass[grep(paste(geneTag,";", sep=""), agClass)]
		agClassGene <- strsplit(agClassGene, split=";", fixed=TRUE)[[1]][3]
		agClassGene <- gsub(" +", " ",agClassGene)
		clusterGeneFunction[cluster,"Agmial_Classification"] <- agClassGene
	}

	

	return(clusterInfo=list(correspondanceOfCDS=correspondanceOfCDS, clusterGeneContent=clusterGeneContent, clusterGeneFunction= clusterGeneFunction))
}
	
	
write.clustersCSV <- function(clusters,pop, fileName="clusters.csv") {
	# csv file (delimitation=semicolon) containing the positions of rec event for each cluster of the list, and information about the gene
	recBlocksCF <- get(paste("recBlocksCF",pop, sep=""))
	blocksPos <- get(paste("blocksPos",pop, sep=""))
	
	strToWr <- character(0)
	for (cluster in clusters) {
		strToWr <- c(strToWr, paste(cluster, recBlocksCF[5,which(recBlocksCF[1,]==cluster)], paste(recBlocksCF[2,which(recBlocksCF[1,]==cluster)],recBlocksCF[3,which(recBlocksCF[1,]==cluster)], sep=" - "), 
		paste( as.numeric(recBlocksCF[2,which(recBlocksCF[1,]==cluster)])-blocksPos10LB[cluster,1],
			as.numeric(recBlocksCF[3,which(recBlocksCF[1,]==cluster)])-blocksPos10LB[cluster,1], sep=" - "),
		sep=";"))
	clusCont <- clusterGeneContent[cluster, ]
	strToWr <- c(strToWr,
				apply(correspondanceOfCDS[clusCont, ,drop=F], 1, function(x) paste(x, collapse = ";")), "")
	}
	cat(strToWr, file = fileName, sep = "\n")
}


construct.recBlocksShow <- function(vit, blocksPos) {
	# returns a matrix containing information about recombined tracts detected by show (vit=0 recombiné)
	# ncol= nb of recombination events
	# rows = cluster, start position, end position, tracts size
	recSites <- which(vit==1)
	startPos <- c()
	endPos <- c()
	clusters <- c()
	for (p in 1:length(recSites)) {
		if ((recSites[p]==1) || (vit[recSites[p]-1]==0)) {#1er site du tract 
			startPos <- c(startPos,p)
			clusterToAdd <- intersect(which(blocksPos[,2] >= recSites[p]),which(blocksPos[,1] <= recSites[p]))
			clusterToAdd <- rownames(blocksPos)[clusterToAdd]
			clusters <- c(clusters, clusterToAdd)
		} else {
			if ((recSites[p]==dim(alignment)[2]) || (vit[recSites[p]+1]==0)) { #dernier site
			endPos <- c(endPos, p)
			}
		}
	}
	recBlocksShow <- matrix(nrow=4, ncol=length(startPos))
	recBlocksShow[1,] <- clusters
	recBlocksShow[2,] <- startPos
	recBlocksShow[3,] <- endPos
	recBlocksShow[4,] <- endPos-startPos
	rownames(recBlocksShow) <- c("cluster","start position","end position","tracts size")
	return(recBlocksShow)
}
	
	
	
multhist.modif <- function (x, beside = TRUE, freq = NULL, probability = !freq, plot.it = TRUE, labels, colToAdd,...) {
	#Meme fonction que multhist mais permet de changer les labels des x
    hist.args <- formals(hist.default)
    args <- list(...)
    hargs <- names(args)[names(args) %in% names(hist.args)]
    hist.args[hargs] <- args[hargs]
    barplot.args <- formals(barplot.default)
    bargs <- names(args)[names(args) %in% names(barplot.args)]
    barplot.args[bargs] <- args[bargs]
    barplot.args$beside <- beside
    barplot.args$... <- barplot.args$inside <- NULL
    allhist <- hist(unlist(x), breaks = hist.args$breaks, plot = FALSE)
   # if (!"names.arg" %in% bargs) {
       # barplot.args$names.arg <- signif(allhist$mids, 2)
    #}
	barplot.args$names.arg <- labels
    if (is.null(freq)) {
        freq <- if (!missing(probability))
            !as.logical(probability)
        else TRUE
    }
    if (freq)
        comp <- "counts"
    else comp <- "density"
    combhist <- t(sapply(x, function(z) hist(z, breaks = allhist$breaks,
        plot = FALSE)[[comp]]))
	combhist <- cbind(combhist, colToAdd)
    if (plot.it)
        do.call("barplot", c(list(combhist), barplot.args))
    invisible(list(breaks = allhist$breaks, out = combhist))
}	


row.names <- function(blocksPos) {	
	blocksNames <- c()	
	for (i in 1:dim(blocksPos)[1]){
	blocksNames <- c(blocksNames, paste("block_", i, sep=""))
	}
	rownames(blocksPos) <- blocksNames
	return (blocksPos)
}
	
	

construct.COGModif <- function(clusterGeneFunction) {
	# matrix containing the names of clusters (column 1) and their COG function (column 2). If a cluster has two functions there are two rows
	COGModif <- matrix()
	cluster <- rownames(clusterGeneFunction)[1]
	COGgene <- clusterGeneFunction[cluster,"COG"]
	COGModif <- matrix(c(rep(cluster, length=length(COGgene[[1]])), COGgene[[1]]), ncol=2)
	for (cluster in rownames(clusterGeneFunction)[-1]) {
		COGgene <- clusterGeneFunction[cluster,"COG"]
		COGgene <- strsplit( COGgene, split=" / ")
		matToAdd <- matrix(c(rep(cluster, length=length(COGgene[[1]])), COGgene[[1]]), ncol=2)
		COGModif <- rbind(COGModif, matToAdd)
	}
	return(COGModif)
}

construct.lengthListByCOG <- function(COGmodif,blocksPos) {
	lengthListByCOG <- list()
	i<-0
	for (funct in COGcat) {
		i <- i+1
		print(funct)
		clusters <- COGmodif[which(COGmodif[,2]== funct),1]
		lengthListByCOG[[i]] <- blocksPos[clusters,2]-blocksPos[clusters,1]
	 }
	 return(lengthListByCOG)
}
	
	
