#1. Write a function that reads fasta files and counts how many sequences are present in the files and how long is each sequence
library(tidyverse)
#Bioconductor require its own package manager
install.packages("BiocManager")
install.packages("seqinr")
library(seqinr)
#Once the manager is installed, can install bioconductor packaged
BiocManager::install("Biostrings")
library(Biostrings)
nuc_fa <- readLines(
  file("https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/practice_nuc.fa", open = "r"))
length(str_which(nuc_fa, ">"))
(length(hdr <- str_which(nuc_fa, ">")))## number of sequences from a file i am just using ido's 
str_which(nuc_fa, ">")[2]
sequences<-list()
for (i in seq_along(hdr)) {
  if(i<5){
  startSeq<-hdr[i]+ 1
  endSeq<-hdr[i+1]-1
  }
  else{
    startSeq<-41
    endSeq<-47
  }
 seqlength<- paste(nuc_fa[startSeq:endSeq])
sequences[[i]]<-seqlength} 
sequences |> lapply(FUN = nchar) |> lapply(FUN = sum)## native pipe rather than %>%
sequences
#2. Write a function that takes a nucleotide sequence as a string and reverse complements it. Make sure that the function differentiate between DNA input and RNA input. Make sure the output is a DNA string
Dseq <- paste(nuc_fa[2:9], collapse = "")  
Rseq<-str_replace_all(Dseq,"T","U")
revComp <- function(seq) {
  # Split the sequence into a character vector
  sq2c <- str_split(seq, "")[[1]]
  
  sq2c <- if_else(sq2c == "A", "T",
                  if_else(sq2c == "T", "A", 
                          if_else(sq2c == "G", "C", 
                                  if_else(sq2c == "C", "G", 
                                          if_else(sq2c == "U", "A", sq2c)))))

  answer <- paste(sq2c, collapse = "")
  
  return(rev(strsplit(answer, NULL)[[1]])) 
}
revComp(Rseq)
revComp(Dseq)

#3. Write a for loop that generates five 500 amino acid strings, make sure to only use the standard 20 AA alphabet. Make sure that the strings are NOT identical
#The loop should store each string as an element of a vector named my_fave_proteins (hint make the vector before the loop).
AA_ALPHABET
my_fave_proteins <- vector("character", 5)
for (i in seq_along(my_fave_proteins)){ 
  my_fave_proteins[i] =paste(sample(AA_ALPHABET,500, replace = T), collapse = "")
  
}
my_fave_proteins
#4. Assign the vector my_fave_proteins to an XStringSet object, name that object my_fave_proteins_set
sq1Set <- AAStringSet(my_fave_proteins)
sq1Set
#5. Change the first AA of each in my_fave_proteins_set protein to M and give each protein a creative name
names(sq1Set) <- c("Beelzebub", "Lucifer", "Belial","Behemoth","Asmodeus")
for (i in seq_along(sq1Set)){
  sq1Set[[i]][1]="M"
}
sq1Set


#6. Get the frequency of each AA in each protein stored in my_fave_proteins_set and plot a box plot of the frequency of each AA in all the proteins (plot using ggplot2)
# Export the plot as PDF
sq1Set.frq <- as.data.frame(alphabetFrequency(sq1Set, baseOnly=T, as.prob = T)*100)

#Let's orgenize it a bit, add the names and plot
sq1Set.frq <- sq1Set.frq %>%  mutate(Sqs = names(sq1Set)) %>% select(Sqs,everything(), -other)%>%
  gather(key = "Bases", value = "Frequency",2:5)

#Not the pretiest, but you can see that the one at the top left has a much lower GC contant than 
#the rest
ggplot(sq1Set.frq, aes(x=Bases,y=Frequency))+
  geom_bar(stat = "identity") +
  facet_wrap(~Sqs)+
  theme_bw()
#7. Export my_fave_proteins_set as a FASTA file, name the file yourfirstname_fave_prots.fasta
writeXStringSet(sq1Set, "ilay_fave_prots.fasta")

#8. Use seqinr to import the first 200 Archaea associated sequences from genbank, how many Archaeal sequences are there in total?
bank <- seqinr::choosebank(bank="genbank", infobank = T)
bank
archealiq <- seqinr::query(listname = "genbank", query = "sp=Archaea")
seqNum<-length(archaeliq$req)
seqNum


#  9. Assign the imported sequences to an XStringSet object named my_fave_arch_set and reverse complement the sequences
#Do not forget to include the sequences’ names
#10. What is the the length of the longest sequence out of the 200 you imported, show your work

#11. Export my_fave_bac_set as a FASTA file, name the file yourfirstname_fave_bac.fa