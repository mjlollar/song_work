library("ggplot2")
library("ggthemes")
library("EnvStats")
library("viridis")
library("tidyverse")

setwd(<YOUR_PATH>)
df_data <- read.csv(<YOUR_SONG_DATA_CSV>', header=T, as.is=T)
df_data$population <- as.factor(df_data$population) #column listing population (FR or ZI)
df_data$genotype <- as.factor(df_data$genotype) #Column listing line number

### TIP: cnrl+F and 'replace all' for 'SinesPerMinNew' with column header of choice to quickly change between measures
#### RUN ESD test on each genotype in genotype colum, will test up to 11 outliers.
#### Prints genotype followed by a sequential tally of outlier number where test statistic > critical value
#### Highest number used as the number of outliers to remove for that genotype
#### Uncomment Line 30 to produce ggplots for each genotype
genos <- as.vector(unique(df_data$genotype))
alpha = 0.05
lam = c(1:11)
R = c(1:11)
for (x in genos) {
  y <- df_data$SinesPerMinNew[df_data$genotype==x]
  ######ONLINE CODE BLOCK###### (Lines 29-62 copied or minimally adapted from: https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm)
  ## Generate normal probability plot.
  #qqnorm(y)
  
  ## Create function to compute the test statistic.
  rval = function(y){
    ares = abs(y - mean(y))/sd(y)
    df = data.frame(y, ares)
    r = max(df$ares)
    list(r, df)}
  n = length(y)
  for (i in 1:11){
    
    if(i==1){
      rt = rval(y)
      R[i] = unlist(rt[1])
      df = data.frame(rt[2])
      newdf = df[df$ares!=max(df$ares),]}
    
    else if(i!=1){
      rt = rval(newdf$y)
      R[i] = unlist(rt[1])
      df = data.frame(rt[2])
      newdf = df[df$ares!=max(df$ares),]}
    
    ## Compute critical value.
    p = 1 - alpha/(2*(n-i+1))
    t = qt(p,(n-i-1))
    lam[i] = t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))
    
  }
  ## Print results.
  newdf = data.frame(c(1:11),R,lam)
  names(newdf)=c("No. Outliers","Test Stat.", "Critical Val.")
#################################################
  
  print(x)
  for (z in 1:11){
    j = newdf$`Test Stat.`[z]
    k = newdf$`Critical Val.`[z]
    if(j > k){
    print(newdf$`No. Outliers`[z])
    }
  }
}
