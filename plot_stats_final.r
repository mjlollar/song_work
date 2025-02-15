library(ggplot2)
library(stringr)
library(dplyr)

song <- read.csv('medianpulseamp_scanone_input_scanoneEmp.csv', header=T, sep="\t") # Output from scanone.py
songPerm <- read.csv('medianpulseamp_perms.csv', header=T, sep="\t")  # Output scanoneBinPerm.py

process_pyscan_data <- function(pyscan1_90) {
  # Check if columns "a" and "d" exist, if yes, remove them
  if ("a" %in% colnames(pyscan1_90) && "d" %in% colnames(pyscan1_90)) {
    pyscan1_90 <- pyscan1_90[, !(colnames(pyscan1_90) %in% c("a", "d"))]
  }
  
  # Display the first row of the data
  cat("First row of pyscan1_90:\n")
  print(pyscan1_90[1,])
  
  # Creating a variable to keep the windows in order
  markers <- strsplit(pyscan1_90$Marker1, '_')
  markers <- as.data.frame(do.call(rbind, markers))
  names(markers) <- c("Chr", "start", "stop")
  pyscan1_90 <- cbind(pyscan1_90, markers$Chr)
  pyscan1_90 <- cbind(pyscan1_90, markers$start)
  markers <- cbind(markers, seq(1, nrow(markers)))
  names(markers) <- c("Chr", "start", "stop", "markerorder")
  pyscan1_90 <- cbind(pyscan1_90, markers$markerorder)
  
  names(pyscan1_90) <- c("Marker1", "LOD", "LODnum", "n", "coef", "r2", "m0", "m1", "m2ZZ", "v0", "v1", "v2", "n0", "n1", "n2", "Chr", "Start", "markerorder")
  pyscan1_90 <- cbind(pyscan1_90, (pyscan1_90$n0 * 2 + pyscan1_90$n1) / (2 * (pyscan1_90$n0 + pyscan1_90$n1 + pyscan1_90$n2)))
  names(pyscan1_90) <- c("Marker1", "LOD", "LODnum", "n", "coef", "r2", "m0", "m1", "m2ZZ", "v0", "v1", "v2", "n0", "n1", "n2", "Chr", "Start", "markerorder", "ancProp")
  
  pyscan1_90$ancProp[pyscan1_90$n0 == -0.25] <- -0.25
  pyscan1_90$ancProp[is.na(pyscan1_90$ancProp)] <- -1.1
  pyscan1_90$LOD[is.na(pyscan1_90$LOD)] <- -1.1
  pyscan1_90$LODnum[is.na(pyscan1_90$LODnum)] <- -1.1
  
  pyscan1_90$Chr <- factor(pyscan1_90$Chr, levels = c("X", "2L", "2R", "3L", "3R"))
  
  return(pyscan1_90)
}

song <- process_pyscan_data(song)

lodnumPvalue <- function(pyscan1_90, scan1Cut90BinPerms) {
  pvCatVector <- c()
  pvNumVector <- c()
  
  for (i in 1:length(pyscan1_90$Chr)) {
    if (pyscan1_90$ancProp[i] <= 0.1 || pyscan1_90$ancProp[i] >= 0.9) {
      totalPerms = 10
      aboveObs = 11
    } else {
      totalPerms = length(scan1Cut90BinPerms$LOD)
      aboveObs = length(scan1Cut90BinPerms$LOD[scan1Cut90BinPerms$LOD > pyscan1_90$LODnum[i]])
    } 
    
    pvCat = aboveObs / totalPerms
    pvCatVector <- c(pvCatVector, pvCat)
  }
  
  pyscan1_90 <- cbind(pyscan1_90, pvCatVector)
  
  pyscan1_90$pvCatVector[pyscan1_90$n0 == -0.25] <- 1.1
  
  pyscan1_90 <- cbind(pyscan1_90, -log(pvCatVector))
  
  names(pyscan1_90) <- c("Marker1", "LOD", "LODnum", "n", "coef", "r2", "m0", "m1", "m2ZZ", "v0", "v1", "v2", "n0", "n1", "n2", "Chr", "Start", "markerorder", "ancProp", "pvalueLOD",  "minusLogPvNum")
  return(pyscan1_90)
}

song <- lodnumPvalue(song, songPerm)

#ML plotting
song$LODnum[song$LODnum==-0.25] <- as.numeric(0)
song$minusLogPvNum[is.infinite(song$minusLogPvNum)==TRUE] <- as.numeric(9.2) #represent as max y
song$minusLogPvNum[song$minusLogPvNum<0] <- as.numeric(0)
end_point <- tail(song$markerorder,1)
breaks_x <- c(1,185,493,790,1179,1464,1785,2106,2378,2674,2947,3212,3519,3818,4124,4378,4630,4885,5189,5467)
labels_x <- c("0", "5", "10", "15", "0", "5", "10", "15", "0", "10", "15", "0", "5", "10", "15", "0", "10", "15", "20","25")
breaks_y <- c(0,2,4,6,8) #edit as needed
labels_y <- c('0','2','4','6','8') #edit as needed
## Numeric (Last edit:05/02/24)
my_plot <- ggplot(data=song, aes(y=LODnum, x=markerorder)) + 
  geom_area(aes(fill=Chr), color="black",linewidth=0.3) +
  scale_fill_manual(values = c("#117733", "#88CCEE", "#DDCC40", "#882255","#332288")) +
  coord_cartesian(ylim=c(0,8)) + #edit as needed
  scale_x_continuous(breaks=breaks_x, labels=labels_x, expand=c(0,0),name=NULL) +
  scale_y_continuous(breaks=breaks_y,name=NULL) +
  geom_hline(yintercept = 3.41287269132781, color = "red", linetype = "dashed") + #edit intercept from below options
  labs(y=NULL) +
  theme_classic() +
  theme(axis.text.x=element_text(size=12, family='Arial', color='black')) + #smaller(10) larger(12)
  theme(axis.text.y=element_text(size=16, family='Arial', color='black')) + #smaller
  #theme(axis.text=element_text(size=16, family='Arial', color='black')) +
  theme(legend.position="none")

my_plot
ggsave(".jpeg", plot=my_plot + theme(), device="jpeg", scale=1, width=2200, height=700, units="px", dpi=300)
#ggsave("slowtototalpulses_013125.jpeg", plot=my_plot + theme(), device="jpeg", scale=1, width=2300, height=1000, units="px", dpi=300) #slowtofasstlarger

#0.05% lod
#3.37468507002206 bouts
#3.49451914692038 nulltosine
#3.43155364585916 medianpulsetrain
#3.39580258534267 mediansinetrain
#3.39285262436019 modepulsemfft
#3.63039121405853 modesinemfft
#3.41287269132781 medianpulseamp
#3.4393777573736 mediansineamp
#3.44537914458102 slowtototal

#residuals
#3.39604532892948 sineamp-pulseamp
#3.4094659524772 sinetrain-pulseamp
#3.44168235068599 pulseamp-sinetrain
#3.40762662719376 pulsemfft-pulseamp
