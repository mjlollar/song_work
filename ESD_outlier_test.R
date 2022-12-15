library("ggplot2")
library("ggthemes")
library("EnvStats")
library("viridis")
library("tidyverse")

setwd(<YOUR_PATH>)
df_data <- read.csv(<YOUR_SONG_DATA_CSV>, header=T, as.is=T)

#### RUN ESD test on each genotype in genotype colum, will test up to 10 outliers.
#### Prints behavior followed by genotype if a significant number of outliers (>=1) detected.
#### Prints greatest number of outliers where test statistic > crtitical value(0.05) below genotype.
#### Currently prints information to STDOUT

genos <- as.vector(unique(df_data$genotype))
alpha = 0.05
lam = c(1:10)
R = c(1:10)

for (calls in colnames(df_data[4:length(df_data)])){    ## Script assumes first three columns of input include [genotype,population,chamber]
  print(calls)
  for (x in genos) {
    y <- df_data[which(df_data$genotype==x), calls]
    ######ONLINE CODE BLOCK###### (Lines 25-57 copied from: https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm, R script 'eda35h3.r')
    ## Generate normal probability plot.
    #qqnorm(y) ## Uncomment to produce ggplots for each genotype
    
    ## Create function to compute the test statistic.
    rval = function(y){
      ares = abs(y - mean(y))/sd(y)
      df = data.frame(y, ares)
      r = max(df$ares)
      list(r, df)}
    n = length(y)
    for (i in 1:10){
      
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
    newdf = data.frame(c(1:10),R,lam)
    names(newdf)=c("No. Outliers","Test Stat.", "Critical Val.")
    #################################################
    out_vec <- c()
    for (z in 1:5){
      j = newdf$`Test Stat.`[z]
      k = newdf$`Critical Val.`[z]
      if(j > k){
        outlier_count <- newdf$`No. Outliers`[z]
        out_vec <- append(out_vec, outlier_count)
      }
    }
    if (length(out_vec) > 0){
      print(x)
      print(out_vec[length(out_vec)])
    }
  }
}
