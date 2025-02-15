library("ggplot2")
library("ggthemes")
library("EnvStats")
library("viridis")
library("tidyverse")
library("Hmisc") #for corr matrix pvalues
library("corrplot") #for heatmap

## To adjust size of circles and postition of significance marks, adjust the base code of "corrplot" package with trace() (see online for details).

df_inbred <- read.csv("inbred_R_work_orderedfiltered.csv", header=T, as.is=T,check.names=F) #filtered data
df_inbred$genotype <- as.factor(df_inbred$genotype)
df_inbred$population <- as.factor(df_inbred$population)
df_ril <- read.csv("RIL_R_work_orderedfiltered.csv", header=T, as.is=T,check.names=F) #filtered data
df_ril$genotype <- as.factor(df_ril$genotype)

###subset tables
#df_inbred_sub <- df_inbred[, c("genotype", "population", "BoutsPerMin", "NullToSine",
#                               "MedianSineTrainLength", "ModeSineCarrierFreq", "MedianSineAmplitudes",
#                               "MedianPulseTrainLength", "ModePulseCarrierFreq", "MedianPulseAmplitudes",
#                               "Slow2FastandSlowPulses")]
#df_ril_sub <- df_ril[, c("genotype", "BoutsPerMin", "NullToSine",
#                         "MedianSineTrainLength", "ModeSineCarrierFreq", "MedianSineAmplitudes",
#                         "MedianPulseTrainLength", "ModePulseCarrierFreq", "MedianPulseAmplitudes",
#                         "Slow2FastandSlowPulses")]
df_inbred_sub <- df_inbred[2:29]
df_inbred_fr <- filter(df_inbred_sub, population == "FR")
df_inbred_sub_fr <- df_inbred_fr[2:28]
df_inbred_zi <- filter(df_inbred_sub, population == "ZI")
df_inbred_sub_zi <- df_inbred_zi[2:28]

df_ril_sub <- df_ril[2:28]

###Pearsons
inbred_corr_fr <- cor(df_inbred_sub_fr, use="complete.obs")
inbred_corr_zi <- cor(df_inbred_sub_zi, use="complete.obs")
inbred_p_fr = cor.mtest(df_inbred_sub_fr, conf.level=0.95)
inbred_p_fr <- inbred_p_fr$p
inbred_p_zi = cor.mtest(df_inbred_sub_zi, conf.level=0.95)
inbred_p_zi <- inbred_p_zi$p

ril_corr <- cor(df_ril_sub, use="complete.obs")
ril_p = cor.mtest(df_ril_sub, conf.level=0.95)
ril_p <- ril_p$p

###combine ril/inbred data. ril top right (upper)
#inbred_corr[upper.tri(inbred_corr)] <- ril_corr[upper.tri(ril_corr)]
#inbred_p[upper.tri(inbred_p)] <- ril_p[upper.tri(ril_p)]
inbred_corr_fr[upper.tri(inbred_corr_fr)] <- inbred_corr_zi[upper.tri(inbred_corr_zi)]
inbred_p_fr[upper.tri(inbred_p_fr)] <- inbred_p_zi[upper.tri(inbred_p_zi)]

#Clean names
#colnames(inbred_corr_fr) <- c("Bouts per Min", "Null to Sine", "Sine Train Length", "Sine Carrier Freq.", "Sine Amplitude",
#                              "Pulse Train Length", "Pulse Carrier Freq.", "Pulse Amplitude", "Slow to Total Pulses")
#rownames(inbred_corr_fr) <- c("Bouts per Min", "Null to Sine", "Sine Train Length", "Sine Carrier Freq.", "Sine Amplitude",
#                              "Pulse Train Length", "Pulse Carrier Freq.", "Pulse Amplitude", "Slow to Total Pulses")
#colnames(inbred_p_fr) <- c("Bouts per Min", "Null to Sine", "Sine Train Length", "Sine Carrier Freq.", "Sine Amplitude",
#                           "Pulse Train Length", "Pulse Carrier Freq.", "Pulse Amplitude", "Slow to Total Pulses")
#rownames(inbred_p_fr) <- c("Bouts per Min", "Null to Sine", "Sine Train Length", "Sine Carrier Freq.", "Sine Amplitude",
#                           "Pulse Train Length", "Pulse Carrier Freq.", "Pulse Amplitude", "Slow to Total Pulses")

#colnames(ril_corr) <- c("Bouts per Min", "Null to Sine", "Sine Train Length", "Sine Carrier Freq.", "Sine Amplitude",
#                        "Pulse Train Length", "Pulse Carrier Freq.", "Pulse Amplitude", "Slow to Total Pulses")
#rownames(ril_corr) <- c("Bouts per Min", "Null to Sine", "Sine Train Length", "Sine Carrier Freq.", "Sine Amplitude",
#                        "Pulse Train Length", "Pulse Carrier Freq.", "Pulse Amplitude", "Slow to Total Pulses")
#colnames(ril_p) <- c("Bouts per Min", "Null to Sine", "Sine Train Length", "Sine Carrier Freq.", "Sine Amplitude",
#                     "Pulse Train Length", "Pulse Carrier Freq.", "Pulse Amplitude", "Slow to Total Pulses")
#rownames(ril_p) <- c("Bouts per Min", "Null to Sine", "Sine Train Length", "Sine Carrier Freq.", "Sine Amplitude",
#                     "Pulse Train Length", "Pulse Carrier Freq.", "Pulse Amplitude", "Slow to Total Pulses")

tiff(file = "/Users/matt/Desktop/10_14_24_new_fig_S1/full_fr_zi_split_inbred_111424.tiff",   # The directory you want to save the file in
     width = 2300, height = 1800, units='px', res=300)
corrplot(inbred_corr_fr, type='lower', insig="label_sig", sig.level=c(0.000000001,0.000001,0.001), p.mat=inbred_p_fr,
         tl.srt=45, tl.col="black", tl.cex=0.7, order='original', diag=FALSE,
         col=viridis::viridis(200), cl.ratio=0.1, cl.cex=0.8,
         cl.pos='n', pch.cex=0.5, pch.col='red',
         addgrid.col="black", family="Palatino")  
dev.off()

#single
tiff(file = "/Users/matt/Desktop/10_14_24_new_fig_S1/full_zi_inbred_111424.tiff",   # The directory you want to save the file in
     width = 2000, height = 1500, units='px', res=300)
corrplot(inbred_corr_zi, type='lower', insig="label_sig", sig.level=c(0.000000001,0.000001,0.001), p.mat=inbred_p_zi,
         tl.srt=45, tl.col="black", tl.cex=0.7, order='original', diag=FALSE,
         col=viridis::viridis(200), cl.ratio=0.1, cl.cex=0.8,
         cl.pos='n', pch.cex=0.5, pch.col='red',
         addgrid.col="black", family="Palatino")  
dev.off()

tiff(file = "/Users/matt/Desktop/10_14_24_new_fig_S1/full_ril_114224.tiff",   # The directory you want to save the file in
     width = 2300, height = 1800, units='px', res=300)
corrplot(ril_corr, type='lower', insig="label_sig", sig.level=c(0.000000001,0.000001,0.001), p.mat=ril_p,
         tl.srt=45, tl.col="black", tl.cex=0.8, order='original', diag=FALSE,
         col=viridis::viridis(200), cl.ratio=0.1, cl.cex=0.8,
         cl.pos='n', pch.cex=0.5, pch.col='red',
         addgrid.col="black", family="Palatino")  
dev.off()
