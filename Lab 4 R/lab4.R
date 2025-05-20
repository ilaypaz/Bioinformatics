
install.packages(c("vegan","ape","ggdendro","seqinr"))
lapply(c("vegan","ape","ggdendro","seqinr",
         "Biostrings", "tidyverse"),library, character.only = TRUE)

#1. Convert the distance matrices generated for the dune data-set (euclidean, and both Bray-Curtis) to an hclust object.
data("dune")
dune.dist.euc <- dist(dune, method = "euclidean")
dune.vd.bry <- vegdist(dune, method = "bray", binary = F)
dune.vd.bry.bin <- vegdist(dune, method = "bray", binary = T)
hc1 <- hclust(dune.dist.euc, method = "average")
hc2 <- hclust(dune.vd.bry, method = "average")
hc3<-hclust(dune.vd.bry.bin,method="average")
hc1
hc2
hc3
#2. Convert the distance matrices generated for the dune data-set to an R matrix (you should have 3 matrices) and save them as CSV files. Make sure the file names are descriptive.
set.seed(1)
dune.euc.nmds <- metaMDS(dune.dist.euc)
# non binary bray
set.seed(1)
dune.bry.nmds <- metaMDS(dune.vd.bry)
#binary bray
set.seed(1)
dune.bry.nmds1 <- metaMDS(dune.vd.bry.bin)

dune.euc.points <- as.data.frame(dune.euc.nmds$points) %>% 
  mutate(Samp = base::rownames(as.data.frame(dune.euc.nmds$points))) %>%
  select(Samp, everything())

dune.bry.points <- as.data.frame(dune.bry.nmds$points) %>% 
  mutate(Samp = base::rownames(as.data.frame(dune.bry.nmds$points))) %>%
  select(Samp, everything())

dune.Binbry.points <- as.data.frame(dune.bry.nmds$points) %>% 
  mutate(Samp = base::rownames(as.data.frame(dune.bry.nmds$points))) %>%
  select(Samp, everything())

write.csv(dune.euc.points, "euclidean.csv")
write.csv(dune.bry.points, "nonbinaryBray.csv")
write.csv(dune.Binbry.points, "binaryBray.csv")


#3. Convert the class hclast variables you have created to class dendrogram
hcd1<- dendro_data(hc1)
hcd2<- dendro_data(hc2)
hcd3<- dendro_data(hc3)
hcd1
hcd2
hcd3

#4. Plot both the class dendrogram and class hclass variables you created for the dune data-set using base R, give title using the argument main. Save plots as jpeg using the RStudio plots tab.
#Do you see any difference in clustering between each of the distance measures?
plot(hc1, xlab="",ylab="", main = "hclust euclidian")
plot(hc2, xlab="",ylab="", main = "hclust bray non binary")
plot(hc3, xlab="",ylab="", main = "hclust bray  binary")
par(mar=c(9.7,4,4,1)+0.1)
plot(as.dendrogram(hc1),xlab="",ylab="", main= "dendrogram euclidian")
plot(as.dendrogram(hc2), xlab="",ylab="", main = "dendrogram non binary bray")
plot(as.dendrogram(hc3), xlab="", ylab="", main = "dendrogram  binary bray")

 
# 5. Convert the dendrogram objects to dendro objects and plot using ggplot2
p.euc <- ggplot(dune.euc.points, aes(x=MDS1, y = MDS2))+
  geom_point(size=2)+
  theme_bw(base_size = 16)+
  ggtitle("NMDS Dune dataset\neuclidean distances")

p.euc


p.bry <- ggplot(dune.bry.points, aes(x=MDS1, y = MDS2))+
  geom_point(size=2)+
  theme_bw(base_size = 16)+
  ggtitle("NMDS Dune dataset\nBray-Curtis distances")

p.bry

p.Binbry <- ggplot(dune.Binbry.points, aes(x=MDS1, y = MDS2))+
  geom_point(size=2)+
  theme_bw(base_size = 16)+
  ggtitle("NMDS Dune dataset\n Binary Bray-Curtis distances")
p.Binbry

#6. Export plots as pdf using ggsave
ggsave("NMDS Dune datasetneuclidean distances.pdf", plot= p.euc)
ggsave("NMDS Dune datasetBray-Curtis distances.pdf", plot= p.bry)
ggsave("NMDS Dune datasetBinary Bray-Curtis distances.pdf",plot=p.Binbry)

#7. Add the column altitude to the label dataframe that’s part of the dednro class variables, have the values “low”, “high”, “sea level” added to each sample randomly.
altitudes <- c("low", "high", "sea level")
hcd1$labels$Altitude <- sample(altitudes, length(hcd1$labels$label), replace = TRUE)
hcd2$labels$Altitude <- sample(altitudes, length(hcd2$labels$label), replace = TRUE)
hcd3$labels$Altitude <- sample(altitudes, length(hcd3$labels$label), replace = TRUE)

#8. Add the column Groups to the segment dataframe of the dendro class variables you generated. Divide the two Bray-Curtis based dendrograms to 3 groups and the euclidean to 2 groups based on the order of the samples.
hcd1$segments <- hcd1$segments %>%mutate(Group = if_else(x < 10, "1", "2"))
hcd2$segments<-hcd2$segments%>%mutate(Group = if_else(x < 7, "1", "2")) %>%mutate(Group = if_else(x > 14, "3", Group))
hcd3$segments<-hcd3$segments%>%mutate(Group=if_else(x<7,"1","2")) %>% mutate(Group=if_else(x>14, "3", Group))

#9. Plot each of the dendro class objects, color the labels based on altitude and the branches based on the groups column. If this was real data would you see any trends in clustering of these samples based on altitude? Make sure to export these plots as pdfs as well
#yeah you probably would use a trend rather than randomizing it
PLOT1 <- ggplot(data = hcd1$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = Group), #create legend 
         size = 1.25, alpha=0.5) + 
geom_text(data = hcd1$labels, aes(x, y-0.01, label = label),
            show.legend = F,
            hjust = 1, angle = 90, size = 3, alpha=0.5)+
theme_bw(base_size = 16)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.93,0.89),
        legend.background = element_blank())+
  scale_color_manual(values = c("dark red", "black", "dark green"))+
  xlab(NULL)+
  ylab(NULL)
PLOT1

PLOT2 <- ggplot(data = hcd2$segments) +  
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = Group), 
               size = 1.25, alpha=0.5) + 
  geom_text(data = hcd2$labels, aes(x, y-0.01, label = label),
            show.legend = F,
            hjust = 1, angle = 90, size = 3, alpha=0.5)+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.93,0.89),
        legend.background = element_blank())+
  scale_color_manual(values = c("dark red", "black", "dark green"))+
  xlab(NULL)+
  ylab(NULL)
PLOT2
# PLOT3
PLOT3 <-ggplot(data = hcd3$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = Group), #create legend 
               size = 1.25, alpha=0.5) + 
  geom_text(data = hcd3$labels, aes(x, y-0.01, label = label),
            show.legend = F,
            hjust = 1, angle = 90, size = 3, alpha=0.5)+
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.93,0.89),
        legend.background = element_blank())+
  scale_color_manual(values = c("dark red", "black", "dark green"))+
  xlab(NULL)+
  ylab(NULL)

PLOT3
ggsave("number 9 plot1.pdf", plot= PLOT1)
ggsave("number 9 plot2.pdf", plot= PLOT2)
ggsave("number 9 plot3.pdf", plot= PLOT3)
#10. Use the ?write.tree and ?write.nexus to export the phylogenetic tree generated in this lab. Give a descriptive file name so that it is clear which is a nexus and which is a newick format.
#Which functions would you have used to import each of the files you exported? for each, if else, and many more
hcp <- as.phylo(hc1)  
write.tree(hcp, file = "prnthc.newick")
newickTree <- read.tree("prnthc.newick")
plot(newickTree)
write.tree(hcp,file="nxs.nexus")
nexusTree<-read.tree("nxs.nexus")
plot(nexusTree)
#11
csvreader <- function(folderPath){
  files <- list.files(path = folderPath, pattern = "\\csv$", full.names = TRUE)
  matrix <- vector("list", length(files))
  i=1
  for(i in seq_along(files)){
    data <- read.csv(files[i])
    log_transform <- log(data)
    e <- metaMDS(log_transform, distance = "bray", binary = T)
    matrix[[i]] <- e
  }
  return(matrix)
}
d=csvreader("/Users/ilaypaz/Downloads/BIOL 3315")
