library(ggplot2)
library(stringr)
library(dplyr)

setwd('final_files_Nov23/')
pyscan1_90 <- read.csv('../boutspermin/boutspermin_scanone_input_scanoneEmp.csv', header=T, sep="\t") # Output from scanone.py
scan1Cut90BinPerms <- read.csv('../boutspermin/boutspermin_combined_perms.csv', header=T, sep="\t", row.names=NULL) # Output from scanoneBinPerm.py

##from TR script 
p005 = -log10(0.05)
labeladjust = 0.2
## Creating a variable to keep the windows in order
markers <- str_split(pyscan1_90$Marker1,'_')
markers <- as.data.frame(do.call(rbind, markers))
names(markers) <- c("Chr","start","stop")
pyscan1_90 <- cbind(pyscan1_90,markers$Chr)
pyscan1_90 <- cbind(pyscan1_90,markers$start)
markers <- cbind(markers, seq(1,length(markers[,1])))
names(markers) <- c("Chr","start","stop","markerorder")
pyscan1_90 <- cbind(pyscan1_90,markers$markerorder)
names(pyscan1_90) <- c("Marker1", "LOD", "LODnum", "n","coef","r2","a","d","mZZ","m1","m2","v0","v1","v2","n0","n1","n2", "Chr", "Start", "markerorder")
pyscan1_90 <- cbind(pyscan1_90,(pyscan1_90$n0*2+pyscan1_90$n1)/(2*(pyscan1_90$n0+pyscan1_90$n1+pyscan1_90$n2)))
names(pyscan1_90) <- c("Marker1", "LOD", "LODnum", "n","coef","r2","a","d","mZZ","m1","m2","v0","v1","v2","n0","n1","n2", "Chr", "Start", "markerorder","zProp")
pyscan1_90$zProp[pyscan1_90$n0==-0.25] = -0.25
pyscan1_90$zProp[is.na(pyscan1_90$zProp)] <- -1.1 # replacing NA with
pyscan1_90$LOD[is.na(pyscan1_90$LOD)] <- -1.1 # replacing NA with
pyscan1_90$LODnum[is.na(pyscan1_90$LODnum)] <- -1.1 # replacing NA with
pyscan1_90$Chr <- factor(pyscan1_90$Chr, levels = c("X","2L","2R","3L","3R"))
pvCatVector <- c()
pvNumVector <- c()
for(i in 1:length(pyscan1_90$Chr)){ # whole data set
  if (pyscan1_90$zProp[i] <= 0.1 | pyscan1_90$zProp[i] >= 0.9) {
    totalPerms = 10
    aboveObs = 11
    aboveObsNum = 11
  } else if (pyscan1_90$zProp[i] <= 0.2 | pyscan1_90$zProp[i] >= 0.8) { # bi
    totalPerms = length(scan1Cut90BinPerms$LOD90) # total number of perms
    aboveObs = length(scan1Cut90BinPerms$LOD90[scan1Cut90BinPerms$LOD90 > pyscan1_90$LOD[i]]) # number of perms higher than observation
    aboveObsNum = length(scan1Cut90BinPerms$LODnum90[scan1Cut90BinPerms$LODnum90 > pyscan1_90$LODnum[i]]) # number of perms higher than observation
  } else if (pyscan1_90$zProp[i] <= 0.3 | pyscan1_90$zProp[i] >= 0.7) {
    totalPerms = length(scan1Cut90BinPerms$LOD80) # total number of perms
    aboveObs = length(scan1Cut90BinPerms$LOD80[scan1Cut90BinPerms$LOD80 > pyscan1_90$LOD[i]]) # number of perms higher than observation
    aboveObsNum = length(scan1Cut90BinPerms$LODnum80[scan1Cut90BinPerms$LODnum80 > pyscan1_90$LODnum[i]]) # number of perms higher than observation
  } else if (pyscan1_90$zProp[i] <= 0.4 | pyscan1_90$zProp[i] >= 0.6) {
    totalPerms = length(scan1Cut90BinPerms$LOD70) # total number of perms
    aboveObs = length(scan1Cut90BinPerms$LOD70[scan1Cut90BinPerms$LOD70 > pyscan1_90$LOD[i]]) # number of perms higher than observation
    aboveObsNum = length(scan1Cut90BinPerms$LODnum70[scan1Cut90BinPerms$LODnum70 > pyscan1_90$LODnum[i]]) # number of perms higher than observation
  } else {
    totalPerms = length(scan1Cut90BinPerms$LOD60) # total number of perms
    aboveObs = length(scan1Cut90BinPerms$LOD60[scan1Cut90BinPerms$LOD60 > pyscan1_90$LOD[i]]) # number of perms higher than observation
    aboveObsNum = length(scan1Cut90BinPerms$LODnum60[scan1Cut90BinPerms$LODnum60 > pyscan1_90$LODnum[i]]) # number of perms higher than observation
  } 
  
  pvCat = aboveObs/totalPerms
  pvNum = aboveObsNum/totalPerms
  pvCatVector <- c(pvCatVector, pvCat)
  pvNumVector <- c(pvNumVector, pvNum)
}
pyscan1_90 <- cbind(pyscan1_90,pvCatVector)
pyscan1_90 <- cbind(pyscan1_90,pvNumVector)
pyscan1_90$pvCatVector[pyscan1_90$n0==-0.25] = 1.1
pyscan1_90$pvNumVector[pyscan1_90$n0==-0.25] = 1.1

pyscan1_90 <- cbind(pyscan1_90,-log10(pvCatVector))
pyscan1_90 <- cbind(pyscan1_90,-log10(pvNumVector))
names(pyscan1_90) <- c("Marker1", "LOD", "LODnum", "n","coef","r2","a","d","mZZ","m1","m2","v0","v1","v2","n0","n1","n2", "Chr", "Start", "markerorder", "zProp", "pvalueLODCat", "pvalueLODNum", "minusLogPvCat", "minusLogPvNum")


#ML plotting
pyscan1_90$LODnum[pyscan1_90$LODnum==-0.25] <- as.numeric(0)
pyscan1_90$minusLogPvNum[is.infinite(pyscan1_90$minusLogPvNum)==TRUE] <- as.numeric(4.1)
pyscan1_90$minusLogPvNum[pyscan1_90$minusLogPvNum<0] <- as.numeric(0)
end_point <- tail(pyscan1_90$markerorder,1)
#breaks_x <- c(185,493,790,1179,1464,1785,2106,2378,2674,2947,3212,3519,3818,4124,4378,4630,4885,5189,5467)
#labels_x <- c("5","10", "15", "0", "5", "10", "15", "0", "10", "15", "0", "5", "10", "15", "0", "10", "15", "20", "25")
breaks_x <- c(1,493,1179,1785,2378,2781,3212,3818,4378,4885,5467)
labels_x <- c("X:0", "10", "2L:0", "10", "2R:0", "12", "3L:0", "10", "3R:0", "15", "25")

my_transform <- 4.1/2.5
lab_1 <- 4
lab_2 <- 3
lab_3 <- 2

## Numeric
my_plot <- ggplot(data=pyscan1_90, aes(y=minusLogPvNum, x=markerorder)) + 
  geom_hline(yintercept=p005, color="black", linetype='dashed') +
  geom_area(aes(x=markerorder, y=(LODnum/my_transform)),fill='#cccccc',  outline.type='lower') +
  geom_point(aes(fill=Chr),shape=21,size=1, color='black', stroke=0.2) +
  scale_fill_manual(values = c("#332288", "#DDCC77","#117733", "#88CCEE", "#882255")) +
  geom_segment(aes(x=end_point-1000, xend=end_point, y=(lab_1/my_transform), yend=(lab_1/my_transform)),linewidth=0.5) +
  geom_segment(aes(x=end_point-30, xend=end_point, y=(lab_2/my_transform), yend=(lab_2/my_transform)),linewidth=0.5) +
  geom_segment(aes(x=end_point-30, xend=end_point, y=(lab_3/my_transform), yend=(lab_3/my_transform)),linewidth=0.5) +
  scale_x_continuous(breaks=breaks_x, labels=labels_x, expand=c(0,0),name=NULL) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +

  geom_vline(xintercept=end_point, linewidth=0.8) +
  labs(y=NULL) +
  theme_classic() +
  theme(axis.text=element_text(size=10, family='Arial', color='black')) + #smaller
  #theme(axis.text=element_text(size=16, family='Arial', color='black')) +
  theme(legend.position="none")

my_plot
ggsave("boutspermin_100423.jpeg", path='/run/media/matt/WD_BLACK/for_SONG_PAPER', plot=my_plot + theme(), device="jpeg", scale=1, width=1400, height=700, units="px", dpi=300)
#ggsave("boutspermin_100423.jpeg", path='/run/media/matt/WD_BLACK/for_SONG_PAPER', plot=my_plot + theme(), device="jpeg", scale=1, width=2800, height=1200, units="px", dpi=300)

