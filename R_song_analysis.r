library("EnvStats")
library("tidyverse")

setwd("/home/matt/song_secondpass/R_work/workspace_finalized_behaviorset/")
df <- read.csv("inbred_R_work_orderedfiltered.csv", header=T, as.is=T)
df <- read.csv("parental_R_work_orderedfiltered.csv", header=T, as.is=T)

#set as factor
df$population <- as.factor(df$population) #column not in RIL or parental
df$genotype <- as.factor(df$genotype)

### Get mean and sd for data
#Initialize dfs
df_summary_mean <- data.frame(genotype=unique(df$genotype))
df_summary_sd <- data.frame(genotype=unique(df$genotype))
df_summary_median <- data.frame(genotype=unique(df$genotype))
#inbred index 3:29
#ril index 2:37
#parental index 3:28
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
population = c("FR","FR","FR","FR","FR","FR","FR","FR","FR","FR","ZI","ZI","ZI","ZI","ZI","ZI","ZI")
df_summary_mean["population"] <- population
df_summary_median["population"] <- population
df_summary_sd["population"] <- population
rm(population) #sanity clear
#counts
list.counts <- df %>% group_by(genotype) %>% count()
df_summary_mean["count"] <- list.counts$n
df_summary_median["count"] <- list.counts$n
df_summary_sd["count"] <- list.counts$n
rm(list.counts) #sanity clear
#Output to csv
write.csv(df_summary_mean,'parental_mean_062823.csv', row.names=FALSE)
write.csv(df_summary_median,'parental_median_062823.csv', row.names=FALSE)
write.csv(df_summary_sd,'parental_sd_062823.csv', row.names=FALSE)

#Get IQR and Mann-Whitney values for inbred population
df_summary_iqr <- data.frame(population=unique(df$population))
#df_summary_iqr <- data.frame(population=unique(df$genotype)) #inbred
mw_pvalue <- vector(mode='logical', length=27)
phenotype <- vector(mode='logical', length=27)
#inbred index: 2:28
for (x in (2:28)){
  #intermediate names
  col.name <- colnames(df_summary_mean)[x]
  iqr.name <- paste(col.name, "iqr", sep="_")
  #Get interquartile value and perform wilcox test (Mann Whitney)
  list.iqr <- data.frame(setNames(aggregate(df_summary_mean[x], list(df_summary_mean$population), FUN=IQR, na.rm=T), c("population", col.name)))
  #Get Mann-Whitney
  df_fr <- df_summary_mean %>% filter(population=='FR') %>% select(all_of(col.name))
  df_fr <- dplyr::pull(df_fr, col.name)
  df_zi <- df_summary_mean %>% filter(population=='ZI') %>% select(all_of(col.name))
  df_zi <- dplyr::pull(df_zi, col.name)
  test.mw <- wilcox.test(x=df_fr, y=df_zi, paired=FALSE, na.rm=T)
  #store IQR
  int.list.iqr <- list(df_summary_iqr, list.iqr)
  df_summary_iqr <- data.frame(reduce(int.list.iqr, full_join, by="population"))
  #store MW pvalue
  y = x - 1
  mw_pvalue[y] <- test.mw$p.value
  phenotype[y] <- col.name
  rm(col.name, iqr.name, list.iqr, df_fr, df_zi, test.mw, int.list.iqr, y) #Sanity clear int variables
}
rm(x) #sanity clear
df_mw <- data.frame(phenotype, mw_pvalue)
#Output to csv
write.csv(df_summary_iqr, 'inbred_iqr_062823.csv', row.names=FALSE)
write.csv(df_mw, 'inbred_mannwhitney_062823.csv', row.names=FALSE)

###MW multiple Test Corrections
pvalues <- vector("logical",10000) #Initiate p.value list
#Permute Mann-Whitney for each trait
for (x in (1:10000)){
  #Shuffle line population identity
  df.shuffle <- transform(df, population = sample(population))
  df.shuffle$population <- as.factor(df.shuffle$population)
  # Calculate M.W. for each test
  p.value.tmp <- vector("logical",27)
  i=1
  for (y in (3:29)){
    col.name <- colnames(df.shuffle)[y]
    df_fr <- df.shuffle %>% filter(population=='FR') %>% select(all_of(col.name))
    df_fr <- dplyr::pull(df_fr, col.name)
    df_zi <- df.shuffle %>% filter(population=='ZI') %>% select(all_of(col.name))
    df_zi <- dplyr::pull(df_zi, col.name)
    test.mw <- wilcox.test(x=df_fr, y=df_zi, paired=FALSE, na.rm=T)
    p.value.tmp[i] <- test.mw$p.value
    i = i + 1
    rm(col.name, df_fr, df_zi, test.mw)
  }
  #Get minimum p.value for replicate tests
  min.pvalue <- min(p.value.tmp)
  pvalues[x] <- min.pvalue
  rm(df.shuffle, i, p.value.tmp)
}
rm(x) #Sanity clear
#sort pvalues, get pvalue @ alpha=0.05, print to output
pvalues <- sort(pvalues)
print(pvalues[500]) # 95% conf p-value, output 06/28/23: 0.002903162
write.csv(pvalues, file="MW_multipletest_perm_pvalues_061623.csv", row.names=FALSE)



#Parental calculations
df_summary_iqr <- data.frame(genotype=unique(df$genotype))
mw_pvalue <- vector(mode='logical', length=27)
phenotype <- vector(mode='logical', length=27)
for (x in (2:28)){
  col.name <- colnames(df)[x]
  #IQR
  list.iqr <- data.frame(setNames(aggregate(df[x], list(df$genotype), FUN=IQR, na.rm=T), c("genotype", col.name)))
  #Mann-Whitney
  fr.group <- df %>% filter(genotype=='FR320N') %>% select(all_of(col.name))
  fr.group <- dplyr::pull(fr.group, col.name)
  zi.group <- df %>% filter(genotype=='ZI418N') %>% select(all_of(col.name))
  zi.group <- dplyr::pull(zi.group, col.name)
  test.mw <- wilcox.test(x=fr.group, y=zi.group, paired=FALSE, na.rm=T)
  y = x - 1
  mw_pvalue[y] <- test.mw$p.value
  phenotype[y] <- col.name
  int.list.iqr <- list(df_summary_iqr, list.iqr)
  df_summary_iqr <- data.frame(reduce(int.list.iqr, full_join, by="genotype"))
  rm(col.name, fr.group, zi.group, test.mw, y, list.iqr, int.list.iqr)
}
rm(x)
df_mw <- data.frame(phenotype, mw_pvalue)
write.csv(df_summary_iqr, 'parental_iqr_062823.csv', row.names=FALSE)
write.csv(df_mw, 'parental_mannwhitney_062823.csv', row.names=FALSE)

