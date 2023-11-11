library('tidyverse')
library('ggplot2')
library('ggthemes')
library('dplyr')
library('tidyr')
library('reshape2')
library('magrittr')
library('hrbrthemes') #for theme ipsum

setwd("slow2fastandslow/")
df_fst <- read.csv('slow2fast/slow2fast_peak1_fst_input_0928.txt', sep='\t')
df_chi <- read.csv('slow2fast_peak1_chi_input_0928.txt', sep='\t')
df_lod <- read.csv('slow2fast_peak1_lod_input_0928.txt', sep='\t')
#df_narrow <- read.csv('../narrowpeak_ids.txt')

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

#breaks_set <- c(8650000, 8750000, 8850000, 8950000, 9050000, 9150000, 9250000) #slowtofasstpeak1
#breaks_set <- c(20090000, 20250000, 20500000, 20750000, 21000000, 21250000, 21500000, 21750000, 22000000, 22250000, 22500000) #slowtofasstpeak2
breaks_set <- c(24820000, 24900000, 25000000, 25100000, 25200000, 25300000, 25370000) #slowtofasstpeak3

#extend lod
df_lod$midpoint[1] <- window_start #slow2fast
#df_lod$midpoint[146] <- window_stop #slow2fast peak2
#df_lod$midpoint[38] <- window_stop #slow2fast peak1
df_lod$midpoint[34] <- window_stop #slow2fast peak3

my_transform <- 3.1/3.5
lab_1 <- 3
lab_2 <- 2
lab_3 <- 1

#plot
my_plot <- ggplot()+
  ##LOD
  geom_line(data=df_lod, aes(x=midpoint, y=(LODnum/my_transform)),color='black') +
  #lod ticks
  geom_segment(aes(x=window_stop-5000, xend=window_stop, y=lab_1/my_transform, yend=lab_1/my_transform)) +
  geom_segment(aes(x=window_stop-5000, xend=window_stop, y=lab_2/my_transform, yend=lab_2/my_transform)) +
  geom_segment(aes(x=window_stop-5000, xend=window_stop, y=lab_3/my_transform, yend=lab_3/my_transform)) +
  
  ###narrow marks, slow2fast
  #geom_segment(data=narrow_peaks, aes(x=margins, xend=margins, y=0, yend=1), color='black', linetype='dotted') +
  
  #lod intercept
  geom_vline(xintercept=tail(df_lod$Stop,1), linewidth=0.8) +
  
  ##Broad Qtl
  geom_line(data=df_fst, aes(x=midpoint, y=logp_snp_cut), color='#EAC659', alpha=0.8, linewidth=.8) +
  geom_line(data=df_fst, aes(x=midpoint, y=logp_win_cut), color='#D81B60',alpha=0.8, linewidth=.8) +
  geom_line(data=df_chi, aes(x=midpoint, y=logp_chi_cut), color='#1E88E5', alpha=0.8, linewidth=.8) +
  
  ###Gene marks, peak 2 slowtofast
  #geom_segment(aes(x = 20334950, xend = 20340828 , y = 4, yend = 4),  color='black', linewidth=1) + #nAChRβ2
  #geom_segment(aes(x = 20350516, xend = 20353255, y = 4, yend = 4),  color='black', linewidth=1) + #Jhbp12
  #geom_segment(aes(x = 20411375, xend = 20449621, y = 4, yend = 4),  color='black', linewidth=1) + #mld
  #geom_segment(aes(x = 21065902, xend = 21067390, y = 4, yend = 4),  color='black', linewidth=1) + #to
  #geom_segment(aes(x = 21173153, xend = 21298780, y = 4, yend = 4),  color='black', linewidth=1) + #fur1
  ##geom_segment(aes(x = 21417224, xend = 21434849, y = 3.98, yend = 3.98),  color='black', linewidth=1) + #Ymp
  ##geom_segment(aes(x = 21681876, xend = 21698392, y = 3.48, yend = 3.48),  color='black', linewidth=1) + #CG5890
  #geom_segment(aes(x = 21749221, xend = 21781992, y = 4, yend = 4),  color='black', linewidth=1) + #CcapR
  #geom_segment(aes(x = 21808735, xend = 21821459, y = 4, yend = 4),  color='black', linewidth=1) + #Nf1
  #geom_segment(aes(x = 21858633, xend = 21859769, y = 4, yend = 4),  color='black', linewidth=1) +#E(spl)m6-BFM
  #geom_segment(aes(x = 21867458, xend = 21875633, y = 4, yend = 4),  color='black', linewidth=1) + #gro
  ##geom_segment(aes(x = 21932590, xend = 21956396, y = 3.5, yend = 3.5),  color='black', linewidth=1) + #Cipc
  #geom_segment(aes(x = 22055970, xend = 22060422, y = 4, yend = 4),  color='black', linewidth=1) + #CG12290
  #geom_segment(aes(x = 22088384, xend = 22102772, y = 4, yend = 4),  color='black', linewidth=1) + #CG6154
  #geom_segment(aes(x = 22227631, xend = 22229378, y = 4, yend = 4),  color='black', linewidth=1) + #ppk15

  ###Gene marks, peak 1 slowtofast
  #geom_segment(aes(x = 8725958, xend = 8738410, y = 2, yend = 2),  color='black', linewidth=1) + # Prm
  #geom_segment(aes(x = 8820605, xend = 8884292, y = 2, yend = 2),  color='black', linewidth=1) + # dally
  #geom_segment(aes(x = 8930200, xend = 8947563, y = 2, yend = 2),  color='black', linewidth=1) + # orb2
  #geom_segment(aes(x = 8998296, xend = 9000349, y = 2, yend = 2),  color='black', linewidth=1) + # doc3
  #geom_segment(aes(x = 9005786, xend = 9012346, y = 2, yend = 2),  color='black', linewidth=1) + # doc2
  #geom_segment(aes(x = 9034475, xend = 9038142, y = 2, yend = 2),  color='black', linewidth=1) + # doc1
  #geom_segment(aes(x = 9138938, xend = 9175249, y = 2, yend = 2),  color='black', linewidth=1) + # Rdl
  #geom_segment(aes(x = 9041881, xend = 9059127, y = 2, yend = 2),  color='black', linewidth=1) + # Argk1
  
  ###genes peak3 slow2fast
  geom_segment(aes(x = 24852161, xend = 24852896, y = 3.5, yend = 3.5),  color='black', linewidth=1) + # eIF4E6
  geom_segment(aes(x = 24864878, xend = 24874870, y = 3.5, yend = 3.5),  color='black', linewidth=1) + # pkc98E
  geom_segment(aes(x = 24938793, xend = 24941919, y = 3.5, yend = 3.5),  color='black', linewidth=1) + # Pdhb
  geom_segment(aes(x = 24949092, xend = 24953969, y = 3.5, yend = 3.5),  color='black', linewidth=1) + # Ctl2
  geom_segment(aes(x = 24996913, xend = 25020716, y = 3.5, yend = 3.5),  color='black', linewidth=1) + # jus
  geom_segment(aes(x = 25047939, xend = 25049150, y = 3.5, yend = 3.5),  color='black', linewidth=1) + # Cisd2
  geom_segment(aes(x = 25140976, xend = 25145366, y = 3.5, yend = 3.5),  color='black', linewidth=1) + # Cnx99A
  geom_segment(aes(x = 25382108, xend = 25390966, y = 3.5, yend = 3.5),  color='black', linewidth=1) + # Dr
  
  #Graph design
  scale_x_continuous(limits=c(window_start, window_stop), breaks=breaks_set, labels=scaleFUNC, name=NULL, expand=c(0,0)) +
  #slow2fast peak1 (1,2) peak2 (1,4), peak3 (1,3.5)
  coord_cartesian(ylim=c(1.0, 3.5)) +
  labs(y=NULL) +
  theme_classic() +
  theme(axis.text=element_text(size=10, family='Arial', color='black')) +
  theme(legend.position = "none")

my_plot
#ggsave("slow2fast_peak1_100623.jpeg", path='/run/media/matt/WD_BLACK/for_SONG_PAPER', plot=my_plot + theme(), device="jpeg", scale=1, width=2400, height=1200, units="px", dpi=300)
ggsave("slow2fast_peak3_100623.jpeg", path='/run/media/matt/WD_BLACK/for_SONG_PAPER', plot=my_plot + theme(), device="jpeg", scale=1, width=1200, height=800, units="px", dpi=300)






######
#plot ZOOM
ggplot()+
  ##LOD
  #transformations(slow2fast: peak2 (/2.375)) // peak1 (/1.5)) // peak3 (*1.129)
  geom_line(data=df_lod, aes(x=midpoint, y=(LODnum/2.375)),color='black') +
  #peak2
  geom_segment(aes(x=22425391-10, xend=22425391, y=2.105263, yend=2.105263)) +
  geom_segment(aes(x=22425391-10, xend=22425391, y=2.947368, yend=2.947368)) +
  geom_segment(aes(x=22425391-10, xend=22425391, y=4, yend=4)) +
  ###narrow marks, slow2fast
  geom_segment(data=narrow_peaks, aes(x=margins, xend=margins, y=0, yend=1), color='black', linetype='dotted') +
  #lod intercept
  geom_vline(xintercept=22425391, linewidth=0.8) +
  ##Broad Qtl
  geom_line(data=df_fst[(df_fst['Start']>=21725712 & df_fst['End']<=22425391),], aes(x=midpoint, y=logp_snp_cut), color='#EAC659', alpha=0.8, linewidth=.8) +
  geom_line(data=df_fst[(df_fst['Start']>=21725712 & df_fst['End']<=22425391),], aes(x=midpoint, y=logp_win_cut), color='#D81B60',alpha=0.8, linewidth=.8) +
  geom_line(data=df_chi[(df_chi['WinStart']>=21725712 & df_chi['WinStop']<=22425391),], aes(x=midpoint, y=logp_chi_cut), color='#1E88E5', alpha=0.8, linewidth=.8) +
  ###Gene marks, peak 2 slowtofast
  geom_segment(aes(x = 21808735, xend = 21821459, y = 3.98, yend = 3.98),  color='black', linewidth=1) + #Nf1
  geom_segment(aes(x = 21858633, xend = 21859769, y = 4, yend = 4),  color='black', linewidth=1) +#E(spl)m6-BFM
  geom_segment(aes(x = 21867458, xend = 21875633, y = 3.98, yend = 3.98),  color='black', linewidth=1) + #gro
  geom_segment(aes(x = 22227631, xend = 22229378, y = 4, yend = 4),  color='black', linewidth=1) + #ppk1
  #Graph design
  scale_x_continuous(limits=c(21725712, 22425391), labels=scaleFUNC, name=NULL, expand=c(0,0)) +
  #slow2fast peak1 (1,2) peak2 (1,4), peak3 (1,3.5)
  coord_cartesian(ylim=c(1.0, 4)) +
  labs(y=NULL) +
  theme_classic() +
  theme(axis.text=element_text(size=20, family='Bakersville')) +
  theme(legend.position = "none")
