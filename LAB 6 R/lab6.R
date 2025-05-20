#Installing the package dendextend
install.packages("dendextend")
#Installing the package ggmsa
BiocManager::install("ggmsa",force=TRUE)

#Installing the package mas from bioconductor
BiocManager::install("msa", force=TRUE)
install.packages("tidyverse")
install.packages("vegan")
install.packages("seqinr")
install.packages("ggmsa")
install.packages("ggplot2",force=true) 
install.packages("stringr")
install.packages("tidyverse",force=true)
install.packages("seqinr",force=true)
packageurl <- "https://cran.r-project.org/src/contrib/ggdendro_0.2.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
download.file(url = "https://cran.r-project.org/src/contrib/Archive/taxize/taxize_0.9.100.1.tar.gz", destfile = "taxize_0.9.100.1.tar.gz")
install.packages(pkgs="taxize_0.9.100.1.tar.gz",  type="source",  repos=NULL)
BiocManager::install("pwalign", force=TRUE)
lapply(c("vegan","ape","ggdendro","seqinr",
         "Biostrings", "tidyverse","dendextend","msa","ggmsa"),library, character.only = TRUE)
library(Biostrings)
library(ape)
library(dendextend)
library(ggmsa)
library(msa)
library(vegan)
library(ggplot2)
library(pwalign)
library(seqinr)
library(stringr)
library(tidyverse)
arch <- seqinr::read.fasta(file = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/Arch.fa")
#Note, using the function lapply to iterate over the list returned by 
unlist(base::lapply(arch, function(x) length(x)))
seqinr::dotPlot(arch[[1]], arch[[3]], col = c("black", "white"))
seqinr::dotPlot(arch[[2]], arch[[3]], col = c("black", "white"))

#1c
seqinr::dotPlot(arch[[1]],arch[[3]], wsize = 50, nmatch = 30)
seqinr::dotPlot(arch[[2]],arch[[3]], wsize = 50, nmatch = 30)

#2a 
selA <- seqinr::read.fasta(file = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/selA.fa")
unlist(base::lapply(selA, function(x) length(x)))
seqinr::dotPlot(selA[[1]],selA[[5]],col=c("black","white"))
#2b
seqinr::dotPlot(selA[[1]],selA[[5]],wsize=10,nmatch=5)
#2c
seqinr::dotPlot(selA[[1]],selA[[5]],wsize=60,nmatch=12)
#3a 
arch_16s <- unlist(lapply(arch[(length(arch)-2):length(arch)],
                          function(x) paste0(x,collapse = "")))

arch_strinS <- DNAStringSet(arch_16s)
selAset <- AAStringSet(unlist(lapply((lapply(selA, base::paste, collapse = "")), stringr::str_to_upper), '[[',1))
nucsubmat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1,type = "DNA")
globalseq13<-pwalign::pairwiseAlignment(arch_16s[[1]],arch_16s[[3]], type = "global") ##seq 1 and 3 nuc's 
globalseq23<-pwalign::pairwiseAlignment(arch_16s[[2]],arch_16s[[3]], type = "global") #seq 2 and 3 nucs
seqlist <- list(globalseq13, globalseq23)

aapw<-pwalign::pairwiseAlignment(selAset[[2]],selAset[[3]], type = "global")## default pairwise alignment for aa's 2-3
aapw62<-pwalign::pairwiseAlignment(selAset[[2]],selAset[[3]], type = "global",substitutionMatrix = "BLOSUM62")#blossum 62 
aapw80<-pwalign::pairwiseAlignment(selAset[[2]],selAset[[3]], type = "global",substitutionMatrix = "BLOSUM80")#blossum 80
aapw120<-pwalign::pairwiseAlignment(selAset[[2]],selAset[[3]], type = "global",substitutionMatrix = "PAM120")#pam 120
aapw40<-pwalign::pairwiseAlignment(selAset[[2]],selAset[[3]], type = "global",substitutionMatrix = "PAM40")#pam 40
aalist<-list(aapw,aapw62, aapw80, aapw120, aapw40)
#3b
globalseq13
globalseq23
#scores are -151.4236 and 3058.172
#I would not expect a negative score for the srna alignment score however i am wrong with my prediction this could be the case since we are dealing with archaea here
#3c
aapw
aapw62
aapw80
aapw120
aapw40
#yes there was a trend all of these global alignment scores are negative. higher level BLOSSUm gave a more negative score and lower level pams gave more negative score which makes sense
#4a
localseq13<-pwalign::pairwiseAlignment(arch_16s[[1]],arch_16s[[3]], type = "local") ##seq 1 and 3 nuc's 
localseq23<-pwalign::pairwiseAlignment(arch_16s[[2]],arch_16s[[3]], type = "local") #seq 2 and 3 nucs
localseqlist <- list(globalseq13, globalseq23)
localaapw<-pwalign::pairwiseAlignment(selAset[[1]],selAset[[5]], type = "local")## default pairwise alignment for aa
localaapw62<-pwalign::pairwiseAlignment(selAset[[1]],selAset[[5]], type = "local",substitutionMatrix = "BLOSUM62")#blossum 62 
localaapw80<-pwalign::pairwiseAlignment(selAset[[1]],selAset[[5]], type = "local",substitutionMatrix = "BLOSUM80")#blossum 80
localaapw120<-pwalign::pairwiseAlignment(selAset[[1]],selAset[[5]], type = "local",substitutionMatrix = "PAM120")#pam 120
localaapw40<-pwalign::pairwiseAlignment(selAset[[1]],selAset[[5]], type = "local",substitutionMatrix = "PAM40")#pam 40
localaalist<-list(aapw,aapw62,aapw80,aapw120,aapw40)
#4b
aligned_pattern <- pattern(localseq13)
aligned_subject <- subject(localseq13)
aligned_pattern_no_gaps <- gsub("-", "", aligned_pattern)
aligned_subject_no_gaps <- gsub("-", "", aligned_subject)

# Get the length of the aligned sequences without gaps
aligned_pattern_length <- nchar(aligned_pattern_no_gaps)
aligned_subject_length <- nchar(aligned_subject_no_gaps)
aligned_pattern_length#length of sequence 1 
aligned_subject_length #length of sequence 3
pattern_positions <- start(aligned_pattern):end(aligned_pattern)
subject_positions <- start(aligned_subject):end(aligned_subject)
pattern_positions#aligned positions of 1
subject_positions#aligned positions of 3

#4c
arch[1]
#yes it is an arachaea 
#Aeropyrum camini is the species in the top result 
#4d aa alignment lengths, they did not change by more than 15 but they did slightly change
nchar(aapw)
nchar(aapw62)
nchar(aapw80)
nchar(aapw120)
nchar(aapw40)
#5a
nucs <- readDNAStringSet(filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBA.fasta")
aas <- readAAStringSet(filepath = "https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/HBA_aa.fasta")
names(nucs) <- apply(str_split_fixed(names(nucs),  pattern = " ",4)[,2:3],1,base::paste, collapse = " ")
names(aas) <- str_remove_all(str_split_fixed(names(aas),pattern = "\\[", n=2)[,2],"\\]")
#5b
hbanucs <- msa(nucs)  # MSA for nucleotide sequences of HBA
hbaaas <- msa(aas)    # MSA for amino acid sequences of HBA

hbbnucs <- msa(nucs)  # MSA for nucleotide sequences of HBB
hbbaas <- msa(aas)    # MSA for amino acid sequences of HBB

hba_aa_msa_ms <- AAStringSet(hbaaas)  # Convert MSA result to AAStringSet
writeXStringSet(hba_aa_msa_ms, filepath = "HBA_aa_alignment.fasta")

hba_n_msa_ms <- DNAStringSet(hbanucs)  # Convert MSA result to DNAStringSet
writeXStringSet(hba_n_msa_ms, filepath = "HBA_nuc_alignment.fasta")

# HBB (similar approach for HBB if desired)
hbb_aa_msa_ms <- AAStringSet(hbbaas)
writeXStringSet(hbb_aa_msa_ms, filepath = "HBB_aa_alignment.fasta")

hbb_n_msa_ms <- DNAStringSet(hbbnucs)
writeXStringSet(hbb_n_msa_ms, filepath = "HBB_nuc_alignment.fasta")


#5c) 
hba_n_msa_consensus <- msaConsensusSequence(hbanucs)
hba_n_msa_ms[length(hba_n_msa_ms) + 1] <- DNAStringSet(
  gsub("\\?", "-", hba_n_msa_consensus))
names(hba_n_msa_ms)[length(hba_n_msa_ms)] <- "Consensus"
print(hba_n_msa_ms[length(hba_n_msa_ms)])

hba_aa_msa_consensus <- msaConsensusSequence(hbaaas)
hba_aa_msa_ms[length(hba_aa_msa_ms) + 1] <- AAStringSet(
  gsub("\\?", "*", hba_aa_msa_consensus))
names(hba_aa_msa_ms)[length(hba_aa_msa_ms)] <- "Consensus"
print(hba_aa_msa_ms[length(hba_aa_msa_ms)])

hbb_n_msa_consensus <- msaConsensusSequence(hbbnucs)
hbb_n_msa_ms[length(hbb_n_msa_ms) + 1] <- DNAStringSet(
  gsub("\\?", "-", hbb_n_msa_consensus))
names(hbb_n_msa_ms)[length(hbb_n_msa_ms)] <- "Consensus"
print(hbb_n_msa_ms[length(hbb_n_msa_ms)])

hbb_aa_msa_consensus <- msaConsensusSequence(hbbaas)
hbb_aa_msa_ms[length(hbb_aa_msa_ms) + 1] <- AAStringSet(
  gsub("\\?", "*", hbb_aa_msa_consensus))
names(hbb_aa_msa_ms)[length(hbb_aa_msa_ms)] <- "Consensus"
print(hbb_aa_msa_ms[length(hbb_aa_msa_ms)])

##aa shows way more detail than nucleotide and show more stuff other than gaps. 
#5d
# HBA
for(i in seq(1, 144, by = 50)) {
  p <- ggmsa(hba_aa_msa_ms, start = i, end = min(i + 49, width(hba_aa_msa_ms)),
             color = "Chemistry_AA", font = NULL) +
    ggtitle(paste("HBA Alignment Positions", i, "to", min(i + 49, width(hba_aa_msa_ms))))
  
  ggsave(filename = paste0("HBA_alignment_positions_", i, "to", min(i + 49, width(hba_aa_msa_ms)), ".pdf"),
         plot = p, device = "pdf", width = 8, height = 6)
}

# HBB
for(i in seq(1, 148, by = 50)) {
  p <- ggmsa(hbb_aa_msa_ms, start = i, end = min(i + 49, width(hbb_aa_msa_ms)),
             color = "Chemistry_AA", font = NULL) +
    ggtitle(paste("HBB Alignment Positions", i, "to", min(i + 49, width(hbb_aa_msa_ms))))
  
  ggsave(filename = paste0("HBB_alignment_positions_", i, "to", min(i + 49, width(hbb_aa_msa_ms)), ".pdf"),
         plot = p, device = "pdf", width = 8, height = 6)
}

#5e
HBAaaAlignSeqinr <- msaConvert(hbaaas, type = "seqinr::alignment")
HBAaaAlignSeqinr
HBAnAlignSeqinr <- msaConvert(hbanucs, type = "seqinr::alignment")
HBAnAlignSeqinr
HBBaaAlignSeqinr <- msaConvert(hbbaas, type = "seqinr::alignment")
HBBaaAlignSeqinr
HBBnAlignSeqinr <- msaConvert(hbbnucs, type = "seqinr::alignment")
HBBnAlignSeqinr

hba_aa_msa_m_dist <- seqinr::dist.alignment(HBAaaAlignSeqinr)
hba_aa_msa_m_dist
hba_n_msa_m_dist <- seqinr::dist.alignment(HBAnAlignSeqinr)
hba_n_msa_m_dist
hbb_aa_msa_m_dist <- seqinr::dist.alignment(HBBaaAlignSeqinr)
hbb_aa_msa_m_dist
hbb_n_msa_m_dist <- seqinr::dist.alignment(HBBnAlignSeqinr)
hbb_n_msa_m_dist

#5f 
hba_aa_msa_m_dend <- as.dendrogram(hclust(hba_aa_msa_m_dist))
plot(hba_aa_msa_m_dend)
hba_n_msa_m_dend <- as.dendrogram(hclust(hba_n_msa_m_dist))
plot(hba_n_msa_m_dend)
hbb_aa_msa_m_dend <- as.dendrogram(hclust(hbb_aa_msa_m_dist))
plot(hbb_aa_msa_m_dend)
hbb_n_msa_m_dend <- as.dendrogram(hclust(hbb_n_msa_m_dist))
plot(hbb_n_msa_m_dend)

tanglegram(hba_aa_msa_m_dend, hba_aa_msa_m_dend, sort = TRUE, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)


#6a 
tanglegram(hba_aa_msa_m_dend, hba_n_msa_m_dend, sort = TRUE, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)
#the top 4 lines, and the bottom two lines are orderly but i can't read what the names are for the most part but its very web like from what i am seeing which is a sign of horizontal gene transfer
#everything but top 4 and bottom two are discrepancies 
#6b
tanglegram(hbb_aa_msa_m_dend, hbb_n_msa_m_dend, sort = TRUE, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)
#its identical to HBA's tanglegram, a mess. 
#6c 
tanglegram(hbb_aa_msa_m_dend, hba_aa_msa_m_dend, sort = TRUE, 
           common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)
#i don't see any discrepancy here it looks perfect!


