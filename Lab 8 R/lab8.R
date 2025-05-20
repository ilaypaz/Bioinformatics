#Installing the package dendextend
install.packages("phangorn")

#Installing the package phytools
install.packages("phytools")

#Installing the package taxize
install.packages("conditionz")
library(conditionz)
download.file(url = "https://cran.r-project.org/src/contrib/Archive/taxize/taxize_0.9.100.1.tar.gz", destfile = "taxize_0.9.100.1.tar.gz")
install.packages(pkgs="taxize_0.9.100.1.tar.gz",  type="source",  repos=NULL)
install.packages("usethis")
library(usethis)
install.packages("remotes")
remotes::install_github("ropensci/bold")

remotes::install_github("ropensci/taxize")

#lapply is part of the apply family of functions, one of those functions that act as a for loop
lapply(c("ape","seqinr","Biostrings", "tidyverse","msa", "phangorn",
         "phytools","stringr"),
       library, character.only = TRUE)

genomes <- readDNAStringSet(
  filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/SARS_MERS_coronavirus.raw_sequence.fasta")
spikes <- readAAStringSet(
  filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/spike.fa")
#CREATING MSA's for these two fasta's 
genomesMSA <- msa(genomes, method = "ClustalW")#MSA for genome
spikesMSA <- msa(spikes, method = "ClustalW")#MSA for spikes
genomesPHY <- as.phyDat(msaConvert(genomesMSA, "seqinr::alignment"),type = "DNA")#unit conversions the mSA's
spikesPHY <- as.phyDat(msaConvert(spikesMSA, "seqinr::alignment"),type = "AA")#unit conversion the MSA's 
#First step is to generate a distance matrix using a substitution model.
#The function dist.ml from phangorn can use 17 different AA substitution models. 
#You will use a fairly simple one. 
genomedml <-dist.ml(genomesPHY, model = "JC69")#distance matrix genome
spikesdml <-dist.ml(spikesPHY, model = "JTT")#distance matrix spikes

#Next, let's generate a neighborjoining tree
#Note you will use NJ (capitalized) from phangorn rather than nj from ape to generate the tree
genomesnj <- phangorn::NJ(genomedml)#create the tree for genome
spikesnj <- phangorn::NJ(spikesdml)#create spikes tree
genomesnj<-midpoint(genomesnj)## rooting for consistency
spikesnj <- midpoint(spikesnj)## rooting for consistency

#Note that the tree is rooted
plot(genomesnj,main="Rooted Neighboring Tree of covid genome")#NJ tree genome
plot(spikesnj,main="Rooted Neighboring Tree of covid protein spikes")#NJ spike tree

#BOOTSTRAPPING THE GENOME
genomesNJBS <-bootstrap.phyDat(genomesPHY,FUN = function(x)NJ(dist.ml(x, model = "JC69")),bs=500)#bootstrap genome
genomesNJBS<- midpoint(genomesNJBS)## rooting for consistency
#BOOTSTRAPPING PROTEIN SPIKES
spikesNJBS <-bootstrap.phyDat(spikesPHY,FUN = function(x)NJ(dist.ml(x, model = "JTT")),bs=500)#bootstrap spikes
spikesNJBS<- midpoint(spikesNJBS)## rooting for consistency

#Note that values smaller than 50% are not printed.
#plotting the NJ bootstrap tree
genomesnj$tip.label#preview the label names
genomesNJBS$tip.label#preview the label names
plotBS(genomesnj,genomesNJBS,"phylogram", main="Rooted Bootstrap Neighbor Joining Tree of covid genomes")#plotting bootstrap genome nj tree
plotBS(spikesnj,spikesNJBS,"phylogram", main="Rooted Bootstrap Neighbor Joining Tree of covid protein spikes")#ployting bootstrap spikes nj tree

##optimization....
genomesnj.pml <- pml(genomesnj, genomesPHY, model = "JC69", k = 4, inv = 0.2)#optimization
spikesnj.pml <- pml(spikesnj, spikesPHY, model = "JTT", k = 4, inv = 0.2)#optimization#optimization
genomesnj.pml <- optim.pml(genomesnj.pml,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)#optimization
spikesnj.pml <- optim.pml(spikesnj.pml,optNni=TRUE,optBf=TRUE,optQ=TRUE,optInv=TRUE,optGamma=TRUE,optEdge=TRUE)#optimization
#optimization and bootstrapping
genomesnj.pml.bs <- bootstrap.pml(genomesnj.pml,bs=100,trees=TRUE,optNni=TRUE)#optimization
spikesnj.pml.bs <- bootstrap.pml(spikesnj.pml,bs=100,trees=TRUE,optNni=TRUE)#optimization
genomesnj.pml.bs <- midpoint(genomesnj.pml.bs) #rooting for consistency
spikesnj.pml.bs <- midpoint(spikesnj.pml.bs)# rooting for consistency

#PLOTTINGT ML TREE
plotBS(genomesnj.pml$tree, genomesnj.pml.bs, type = "phylogram", main="Rooted Bootstrap MAXIMUM LIKELIHOOD Tree of COVID GENOME")#bootstap covid genome ML TREE
plotBS(spikesnj.pml$tree, spikesnj.pml.bs, type = "phylogram", main="Rooted Bootstrap MAXIMUM LIKELIHOOD Tree of COVID PROTEIN SPIKES")# bootstrap covid spikes ML TREE
bootstrap.pml(genomesnj.pml)#calculating ML bootstrap value for genome
bootstrap.pml(spikesnj.pml)#calculating ML bootstrap value for spikes
bootstrap.phyDat(genomesPHY,FUN=function(x)NJ(dist.hamming(x)))#calculating NJ genome bootstrap value
bootstrap.phyDat(spikesPHY,FUN=function(x)NJ(dist.hamming(x)))#calculating NJ spikes bootstrap value
#the  NJ values obtained are 100 phylogenetic trees for 
#ML  genome value is 113368.7 , spikes ML value is 9242.266

#Making Backups
genomesnj2 <- genomesnj#making backups

spikesnj2<-  spikesnj#making backups

#making backups
genomesnj.pml2 <- genomesnj.pml#making backups

spikesnj.pml2 <- spikesnj.pml#making backups


sp1 <- genomesnj2$tip.label#retrieve tip labels from NJ TREE OF genome
sp2 <- genomesnj.pml2$tree$tip.label#retrieve tip labels from ML TREE OF genome

sp1b <- spikesnj2$tip.label#retrieve tip labels from NJ TREE OF SPIKES
sp2b <- spikesnj.pml2$tree$tip.label#retrieve tip labels from ML TREE OF SPIKES

uids <- taxize::get_uid(sp1) #retrieve UID's 

uids2 <- taxize::get_uid(sp2)#retrieve UID's 

uidsb <- taxize::get_uid(sp1b)#retrieve UID's 

uids2b <- taxize::get_uid(sp2b)#retrieve UID's 

cnames <- taxize::sci2comm(uids, "ncbi")#retrieve common name from UIDS
cnames2 <- taxize::sci2comm(uids2, "ncbi")#retrieve common name from UIDS

cnamesb <- taxize::sci2comm(uidsb, "ncbi")#retrieve common name from UIDS
cnames2b <- taxize::sci2comm(uids2b, "ncbi")#retrieve common name from UIDS

a <- unlist(lapply(cnames, "[",1), use.names = F) #narrows down common name and puts it into vector

a2 <- unlist(lapply(cnames2, "[",1), use.names = F)#narrows down common name and puts it into vector

b <- unlist(lapply(cnamesb, "[",1), use.names = F)#narrows down common name and puts it into vector

b2 <- unlist(lapply(cnames2b, "[",1), use.names = F)#narrows down common name and puts it into vector

which(is.na(a))#checking to see if there are any NA's in the vector (no) 
which(is.na(a2))#checking to see if there are any NA's in the vector (no) 
which(is.na(b))#checking to see if there are any NA's in the vector (no) 
which(is.na(b2))#checking to see if there are any NA's in the vector (no) 

a[which(is.na(a))] <- sp1[which(is.na(a))] #Identifies the indices of NA values (none) 
a2[which(is.na(a2))] <- sp2[which(is.na(a2))]#Identifies the indices of NA values (none) 
b[which(is.na(b))] <- sp1b[which(is.na(b))]#Identifies the indices of NA values (none) 
b2[which(is.na(b2))] <- sp2b[which(is.na(b2))]#Identifies the indices of NA values (none) 

a #staring at species names
a2#staring at species names
b#staring at species names
b2#staring at species names

genomesnj2$tip.label <- a#assign labels to the NJ tree
genomesnj.pml2$tree$tip.label <- a2#assign labels to the ML tree
spikesnj2$tip.label <- b#assign labels to the NJ tree
spikesnj.pml2$tree$tip.label <- b2#assign labels to the ML tree
obj <- cophylo(genomesnj2, genomesnj.pml2$tree)#conctructing object for tanglegrams gennomes
objb <- cophylo(spikesnj2, spikesnj.pml2$tree)#conctructing object for tanglegrams spikes
objML<-cophylo(genomesnj.pml2$tree, spikesnj.pml2$tree)#conctructing object for tanglegrams ML genomes vs ML proteins.

#Plotting the trees
par(mfrow=c(1,1))
#mar sets the margins of the plot
plot(obj,mar=c(.1,.1,2,.1))
#Title is self explanatory
title("NJ genome                                                                   ML genome")

#Plotting the trees
par(mfrow=c(1,1))
#mar sets the margins of the plot
plot(objb,mar=c(.1,.1,2,.1))
#Title is self explanatory
title("NJ spikes                                                                   ML spikes")

#Plotting the trees
par(mfrow=c(1,1))
#mar sets the margins of the plot
plot(objML,mar=c(.1,.1,2,.1))
#Title is self explanatory
title("ML genome                                                                   ML spikes")



  
