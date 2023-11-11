library('tidyverse')
library('ggplot2')
library('ggthemes')
library('dplyr')
library('tidyr')
library('reshape2')
library("EnvStats")
library('magrittr')
library('stringr')
library('hrbrthemes') #for theme ipsum

setwd("/home/matt/song_secondpass/R_work/workspace_finalized_behaviorset/10_24_23_work/")
df <- read.csv("inbred_R_work_orderedfiltered.csv", header=T, as.is=T)


#Convert 0 values to NA so they are excluded from calcs
#Convert NaN explicitly to NA to avoid possible wonky behavior
df <- df %>% mutate_all(~ifelse(is.nan(.), NA, .))
df[df==0] <- NA
#set as factor
df$population <- as.factor(df$population) #column not in RIL or parental
df$genotype <- as.factor(df$genotype)

### Get mean and sd for data
#Initialize dfs
df_summary_mean <- data.frame(genotype=unique(df$genotype))
df_summary_sd <- data.frame(genotype=unique(df$genotype))
df_summary_median <- data.frame(genotype=unique(df$genotype))

###Output sample sizes
df %>% group_by(genotype) %>% count()

#inbred index 3:29
for (x in (3:29)){
  #intermediate names
  col.name <- colnames(df)[x]
  mean.name <- paste(col.name, "mean", sep="_")
  sd.name <- paste(col.name, "sd", sep="_")
  median.name <- paste(col.name, "median", sep="_")
  #Get mean and SD by genotype
  list.mean <- data.frame(setNames(aggregate(df[x], list(df$genotype), FUN='mean', na.rm=T), c("genotype", col.name)))
  list.sd <- data.frame(setNames(aggregate(df[x], list(df$genotype), FUN='sd', na.rm=T), c("genotype", col.name)))
  list.median <- data.frame(setNames(aggregate(df[x], list(df$genotype), FUN='median', na.rm=T), c("genotype", col.name)))
  #Join summaries
  int.list.mean <- list(df_summary_mean, list.mean)
  int.list.sd <- list(df_summary_sd, list.sd)
  int.list.median <- list(df_summary_median, list.median)
  df_summary_mean <- data.frame(reduce(int.list.mean, full_join, by="genotype"))
  df_summary_sd <- data.frame(reduce(int.list.sd, full_join, by="genotype"))
  df_summary_median <- data.frame(reduce(int.list.median, full_join, by="genotype"))
  
  rm(col.name, sd.name, mean.name, median.name, list.mean, list.sd, list.median, int.list.mean, int.list.sd, int.list.median) #Sanity clear int variables
}
rm(x)

#add population category, for inbred
#This is index dependent, works with me but maybe not for thee
population = c("FR","FR","FR","FR","FR","FR","FR","FR","FR","FR","ZI","ZI","ZI","ZI","ZI","ZI","ZI")
df_summary_mean["population"] <- population
df_summary_median["population"] <- population
df_summary_sd["population"] <- population
rm(population) #sanity clear
#Output to csv
#write.csv(df_summary_mean,'parental_mean_102823.csv', row.names=FALSE)
#write.csv(df_summary_median,'parental_median_102823.csv', row.names=FALSE)
#write.csv(df_summary_sd,'parental_sd_102823.csv', row.names=FALSE)

### Plotting
my_plot <- ggplot(data=df, aes(x=genotype, y=MedianPulseAmplitudes))+
  geom_violin(aes(fill=population), width=1, alpha=0.8, na.rm=TRUE)+#, draw_quantiles=c(0.5)) +
  scale_fill_manual(values=c('#1E88E5','#D81B60')) +
  geom_boxplot(width=.2, color="black", alpha=0.2, outlier.shape=NA) +
  #coord_cartesian(ylim=c(0.5, 1.0)) + #slow2fast
  coord_cartesian(ylim=c(0.0, 0.25)) + #pulse amps
  labs(y=NULL, x=NULL) +
  theme_classic() +
  theme(axis.text=element_text(size=8, family='Arial', color='black')) +
  theme(legend.position = "none")

ggsave("medianpulseamp_inbredmeans_102623.jpeg", plot=my_plot + theme(), device="jpeg", scale=1, width=2400, height=1200, units="px", dpi=300)
rm(my_plot)


  
  
  
  
