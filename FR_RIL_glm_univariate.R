df <- read.csv('<input>', header=T, as.is=T) #removed low/high song (1<x>55), high ambient/other, aberrant sine/pulse per min measures (first_pass_data_removal.py)

#load some libs
library("ggplot2")
library("ggthemes")
library("EnvStats")
library("tidyverse")
library("dplyr")
library('arm')
library('viridis')
options(max.print=10000) #ensure summary() calls print full lists to stdout, for larger models

### Code mostly set up to work for any phenotype, use Find and Replace All to quickly alternate between phenos
### (#ADJUST comment marks where alterations are needed when switching between phenos

### Visualize data (Slow2FastandSlowPulses used in this script)
summary(df$Slow2FastandSlowPulses)
ggplot(df, aes(x=genotype, y=Slow2FastandSlowPulses)) +
  geom_boxplot(show.legend=FALSE, outlier.shape=NA)

#### Univariate models

# NOTE: Univariate models find the category who's mean trait value is closest to the mean,
# and sets this category as the reference category in each model. When changing phenotypes,
# you will need to take the output from which.min() [first instance in line 34] and adjuat the call to Group.1 in the 
# following two lines (after each which.min() call) to get the appropriate index to the category 
# closest to the mean, and re-level the factor.

### univariate genotype
df$genotype <- as.factor(df$genotype)
group_genotype <- aggregate(df$Slow2FastandSlowPulses,list(df$genotype),mean)
mean.genotype <- mean(group_genotype$x)
diff.genotype <- abs(group_genotype$x - mean.genotype) 
which.min(diff.genotype)  # example output: "index" (numeric)
group_genotype$Group.1['index'] # example output: "category" (str, category closest to mean)
df$genotype <- relevel(df$genotype,ref='category')
glm.uni.genotype <- glm(formula=Slow2FastandSlowPulses ~ genotype, data=df)
summary(glm.uni.genotype)

### Chamber model (1-32)
df$chamber <- as.factor(df$chamber)
group_chamber <- aggregate(df$Slow2FastandSlowPulses,list(df$chamber),mean)
mean.chamber <- mean(group_chamber$x)
diff.chamber <- abs(group_chamber$x - mean.chamber) 
which.min(diff.chamber) 
group_chamber$Group.1["index"] # ADJUST
df$chamber <- relevel(df$chamber,ref="category") #ADJUST
glm.uni.chamber <- glm(formula=Slow2FastandSlowPulses ~ chamber, data=df)
summary(glm.uni.chamber)
### Assay Group (1-3, assayed minutes from lights on: 0-30,31-60,61-90)
df$assay_group <- as.factor(df$assay_group)
group_assay.group <- aggregate(df$Slow2FastandSlowPulses,list(df$assay_group),mean)
mean.assay.group <- mean(group_assay.group$x)
diff.assay.group <- abs(group_assay.group$x - mean.assay.group) 
which.min(diff.assay.group) 
group_assay.group$Group.1['index'] # ADJUST
df$assay_group <- relevel(df$assay_group,ref='category') #ADJUST
glm.uni.assay_group <- glm(formula=Slow2FastandSlowPulses ~ assay_group, data=df)
summary(glm.uni.chamber)
### Incubator
df$incubator <- as.factor(df$incubator)
group_incubator <- aggregate(df$Slow2FastandSlowPulses,list(df$incubator),mean)
mean.incubator <- mean(group_incubator$x)
diff.incubator <- abs(group_incubator$x - mean.incubator) 
which.min(diff.incubator) 
group_incubator$Group.1['index'] #ADJUST
df$incubator <- relevel(df$incubator,ref='category') #ADJUST
glm.uni.incubator <- glm(formula=Slow2FastandSlowPulses ~ incubator, data=df)
summary(glm.uni.incubator)
### Date (Day of assay, 1-37 from past -> present June2022-July2022)
df$date <- as.factor(df$date)
group_date <- aggregate(df$Slow2FastandSlowPulses,list(df$date),mean)
mean.date <- mean(group_date$x)
diff.adate <- abs(group_date$x - mean.date) 
which.min(diff.date) 
group_date$Group.1['index'] # ADJUST
df$date <- relevel(df$date,ref='category') # ADJUST
glm.date <- glm(formula=Slow2FastandSlowPulses ~ date, data=df)
summary(glm.uni.chamber)
