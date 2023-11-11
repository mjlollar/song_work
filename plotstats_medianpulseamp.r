library('tidyverse')
library('ggplot2')
library('ggthemes')
library('dplyr')
library('tidyr')
library('reshape2')
library('magrittr')
library('hrbrthemes') #for theme ipsum

setwd("medianpulseamp/")
df_fst <- read.csv('medianpulseamp_peak10_fst_input_0930.csv', sep='\t')
df_chi <- read.csv('medianpulseamp_peak10_chi_input_0930.csv', sep='\t')
df_lod <- read.csv('medianpulseamp_peak10_lod_input_0930.csv', sep='\t')
#df_narrow <- read.csv('narrow_peak_ids_medianpulseamp.txt', sep='\t')

df_lod <- df_lod %>%
  separate_wider_delim(Marker, delim="_", names= c('chrom', 'Start', 'Stop'))
df_lod$Start <- as.numeric(df_lod$Start)
df_lod$Stop <- as.numeric(df_lod$Stop)
df_lod$LODnum <- as.numeric(df_lod$LODnum)

scaleFUNC <- function(x){
  value <- x/1000000
  value <- format(round(value, 2), nsmall=2)
  return(value) #function to print x axis in MB
}

#get midpoint of base
midpoint = function(x, output){
  if (length(x) == 7){ #FST
    a = (as.numeric(x[1]) + as.numeric(x[2])) / 2
  } else if (length(x) == 10){ #chi
    a = (as.numeric(x[2]) + as.numeric(x[3])) / 2
  } else { #lod
    a = (as.numeric(x[2]) + as.numeric(x[3])) / 2
  }
  return(a)
}

#create midpoint columns for plotting
mid_fst <- apply(df_fst, 1, midpoint)
df_fst["midpoint"] = mid_fst
mid_chi <- apply(df_chi, 1, midpoint)
df_chi["midpoint"] = mid_chi
mid_lod <- apply(df_lod, 1, midpoint)
df_lod["midpoint"] = mid_lod
rm(mid_fst, mid_chi, mid_lod, midpoint)

#Get (-log) Quantile
win_fst_log = function(x, output){
  value = -log10(as.numeric(x[6]))
  return(value)
}
chi_log = function(x, output){
  value = -log10(as.numeric(x[10]))
  return(value)
}
snp_fst_log = function(x, output){
  value = -log10(as.numeric(x[7]))
  return(value)
}
final_quant_winfst <- apply(df_fst, 1, win_fst_log)
df_fst['logp_win'] = final_quant_winfst
final_quant_chi <- apply(df_chi, 1, chi_log)
df_chi['logp_chi'] = final_quant_chi
final_quant_snpfst <- apply(df_fst, 1, snp_fst_log)
df_fst['logp_snp'] = final_quant_snpfst
rm(final_quant_chi, final_quant_snpfst, final_quant_winfst, win_fst_log, chi_log, snp_fst_log)

##plotting variables
window_start <- as.numeric(df_fst$Start[1])
window_stop <- as.numeric(tail(df_fst$End,1))
list_windows <- c(seq(window_start, window_stop, by=1))

#set below 1 to one for plot
logp_filtered_fst <- vector(mode="numeric", length=(length(df_fst$Start)))
for (index in 1: length(logp_filtered_fst)){
  if (df_fst$logp_win[index] < 1){
    logp_filtered_fst[index] <- as.numeric(1.0)
  } else{
    logp_filtered_fst[index] <- df_fst$logp_win[index]
  }
}
df_fst['logp_win_cut'] <- logp_filtered_fst

logp_filtered_chi <- vector(mode="numeric", length=(length(df_chi$WinStart)))
for (index in 1: length(logp_filtered_chi)){
  if (df_chi$logp_chi[index] < 1){
    logp_filtered_chi[index] <- as.numeric(1.0)
  } else{
    logp_filtered_chi[index] <- df_chi$logp_chi[index]
  }
}
df_chi['logp_chi_cut'] <- logp_filtered_chi

logp_filtered_snp <- vector(mode="numeric", length=(length(df_fst$Start)))
for (index in 1: length(logp_filtered_snp)){
  if (df_fst$logp_snp[index] < 1){
    logp_filtered_snp[index] <- as.numeric(1.0)
  } else{
    logp_filtered_snp[index] <- df_fst$logp_snp[index]
  }
}
df_fst['logp_snp_cut'] <- logp_filtered_snp

#narrow_peaks <- data.frame(margins= sort(c(dplyr::pull(df_narrow, start), dplyr::pull(df_narrow, stop))))

#breaks_set <- c(5160000, 5500000, 5800000, 6100000, 6400000, 6700000, 7000000, 7400000) #slowtofasstpeak3

#extend lod
df_lod$midpoint[1] <- window_start
df_lod$midpoint[36] <- window_stop # 18, 139, 28, 25, 37, 52, 26, 43, 27, 36

#plot NARROW
ggplot()+
  ##LOD
  geom_line(data=df_lod, aes(x=midpoint, y=(LODnum/(4.1/3.1))),color='black') +
  geom_segment(aes(x=window_stop-1000, xend=window_stop, y=3.1, yend=3.1)) +
  geom_segment(aes(x=window_stop-1000, xend=window_stop, y=2.2683, yend=2.2683)) +
  geom_segment(aes(x=window_stop-1000, xend=window_stop, y=1.5122, yend=1.5122)) +
  
  ###narrow marks
  #geom_segment(data=narrow_peaks, aes(x=margins, xend=margins, y=0, yend=1), color='black', linetype='dotted') +
  
  #lod intercept
  geom_vline(xintercept=tail(df_lod$Stop,1), linewidth=0.8) +
  
  ##Broad Qtl
  geom_line(data=df_fst, aes(x=midpoint, y=logp_snp_cut), color='#EAC659', alpha=0.8, linewidth=.8) +
  geom_line(data=df_fst, aes(x=midpoint, y=logp_win_cut), color='#D81B60',alpha=0.8, linewidth=.8) +
  geom_line(data=df_chi, aes(x=midpoint, y=logp_chi_cut), color='#1E88E5', alpha=0.8, linewidth=.8) +
  
  ###Gene marks, peak 1
  #geom_segment(aes(x = 4803077, xend = 4805274 , y = 3.1, yend = 3.1),  color='black', linewidth=1) + #CG7509
  #geom_segment(aes(x = 4861269, xend = 4895749 , y = 3.1, yend = 3.1),  color='black', linewidth=1) + #Dnah3
  #geom_segment(aes(x = 4910368, xend = 4920198 , y = 3.1, yend = 3.1),  color='black', linewidth=1) + #Rh50
  
  #Graph design
  #scale_x_continuous(limits=c(window_start, window_stop), breaks=breaks_set, labels=scaleFUNC, name=NULL, expand=c(0,0)) +
  scale_x_continuous(limits=c(window_start, window_stop), labels=scaleFUNC, name=NULL, expand=c(0,0)) +
  coord_cartesian(ylim=c(1.0, 3.1)) +
  labs(y=NULL) +
  theme_classic() +
  theme(axis.text=element_text(size=20, family='Bakersville')) +
  theme(legend.position = "none")

