#### RUN ESD test on each genotype in genotype colum, will test up to 10 outliers.
#### Removes greatest number of outliers where test statistic > crtitical value(0.05) below genotype.
#### Outputs csv with outliers replaced as NA values
#### Modify where <THIS_STYLE_IS_PRESENT>
### Last update MJL 12.16.22

library("ggplot2")
library("ggthemes")
library("EnvStats")
library("viridis")
library("tidyverse")

setwd(<YOUR_PATH>)
df_data <- read.csv(<YOUR_SONG_DATA_CSV>, header=T, as.is=T)

#For ESD test, adjust as desired
alpha = 0.05
lam = c(1:10)
R = c(1:10)

#requires as a genotype column and first three columns as non-data columns (ex. genotype,population,chamber)
genos <- as.vector(unique(df_data$genotype))
phenos <- colnames(df_data[4:length(df_data)])

#function to get distance from mean
get_diff <- function(x){
  z = abs(x - mean_set)
  return(z)}

#get max row possible in df
max_row = nrow(df_data)

#initialise output dataframe
data_genotypes <- as.vector(df_data$genotype)
df_out <- data.frame(genotype=data_genotypes)

#test for outliers within each genotype for each phenotype
for (calls in phenos){
  count_out <- vector(length=max_row)
  ticker = 1
  for (x in genos) {
    y <- df_data[which(df_data$genotype==x), calls]
    og_y_length <- length(y)
    ######ONLINE CODE BLOCK###### (Lines 44-71 copied from: https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm)
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
    newdf = data.frame(c(1:10),R,lam)
    names(newdf)=c("No. Outliers","Test Stat.", "Critical Val.")
    #################################################
    
    out_vec = 0 #Highest outlier number
    #Get maximum outlier count
    for (z in 1:10){
      j = newdf$`Test Stat.`[z]
      k = newdf$`Critical Val.`[z]
      if(j > k){
        out_vec = as.numeric(newdf$`No. Outliers`[z])
      }
    }
    
    #Replace outliers if detected with NA values in dataset
    #Find outlier index/indecies if present 
    if (out_vec > 1)
      mean_set <- mean(y) #get mean of set
      y_diff <- unlist(lapply(y, FUN=get_diff)) #calculate distances
      y_sort <- sort(y_diff, decreasing=TRUE, index.return=TRUE)$ix[1:out_vec]#Get index of highest diffs
      y <- y[- y_sort] #remove outliers by indexed position
      out_mean <- mean(y)
      print(out_mean,quote=FALSE)}
    else if (out_vec == 1){
      mean_set <- mean(y)
      y_diff <- unlist(lapply(y, FUN=get_diff)) #calculate distances
      y_sort <- sort(y_diff, decreasing=TRUE, index.return=TRUE)$ix[1] #Get index of highest diffs
      y <- y[-y_sort] #remove outliers by indexed position
      out_mean <- mean(y)
      print(out_mean,quote=FALSE)} 
    else {
      out_mean <- mean(y)
      print(out_mean,quote=FALSE)}
    
    # Add values to output vector
    y <- c(y, rep(NA, og_y_length - length(y)))
    for (d in y){
      count_out[ticker] <- d
      ticker = ticker + 1
    }
  }
  df_out[calls] <- count_out #add corrected phenotype column to output df
}

write.csv(df_out, row.names=False, file=<YOUR_OUTFILE>) #output df to cwd
