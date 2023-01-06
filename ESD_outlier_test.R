library("ggplot2")
library("ggthemes")
library("EnvStats")
library("viridis")
library("tidyverse")

df_data <- read.csv('<YOUR_DATA.csv>', header=T, as.is=T)

alpha = 0.05
lam = c(1:10)
R = c(1:10)

genos <- as.vector(unique(df_data$genotype))
phenos <- colnames(df_data[4:length(df_data)])

#function to get distance from mean
get_diff <- function(x){
  z = abs(x - mean_set)
  return(z)}

#get max row possible in df
max_row = nrow(df_data)

#Initialise output dataframes
data_genotypes <- as.vector(df_data$genotype)
data_genotypes_sub <- unique(df_data$genotype)
df_out <- data.frame(genotype=data_genotypes)
df_out_means <- data.frame(genotype=data_genotypes_sub, group_mean=rep(NA, length(data_genotypes_sub)))


for (calls in phenos){
  count_out <- vector(length=max_row)
  ticker = 1
  ticker_sub = 1
  df_out_means[calls] <- rep(NA,length(data_genotypes_sub)) #Make mean column for phenotype
  for (x in genos) {
    y <- df_data[which(df_data$genotype==x), calls]
    og_y_length <- length(y)
    ######ONLINE CODE BLOCK###### (Lines 40-70 copied from: https://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm)
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
    
    out_vec = 0 #Highest outlier number
    #Get maximum outlier count
    for (z in 1:10){
      j = newdf$`Test Stat.`[z]
      k = newdf$`Critical Val.`[z]
      if(j > k){
        out_vec = as.numeric(newdf$`No. Outliers`[z])
      }
    }
    
    #Remove outliers if detected
    if (out_vec > 1){
      #Find outliers 
      mean_set <- mean(y) #get mean of set
      y_diff <- unlist(lapply(y, FUN=get_diff)) #calculate distances
      y_sort <- sort(y_diff, decreasing=TRUE, index.return=TRUE)$ix[1:out_vec]#Get index of highest diffs
      y <- y[- y_sort] #remove outliers by indexed position
      out_mean <- mean(y)
      if (out_vec == 10){
        print("10 outliers detected. Behavior/Genotype:")
        print(calls, quote=FALSE)
        print(x, quote=FALSE)}}
    else if (out_vec == 1){
      mean_set <- mean(y)
      y_diff <- unlist(lapply(y, FUN=get_diff)) #calculate distances
      y_sort <- sort(y_diff, decreasing=TRUE, index.return=TRUE)$ix[1] #Get index of highest diffs
      y <- y[-y_sort] #remove outliers by indexed position
      out_mean <- mean(y)}
    else {
      out_mean <- mean(y)}
    
    # Add values to output vector
    y <- c(y, rep(NA, og_y_length - length(y)))
    for (d in y){
      count_out[ticker] <- d
      ticker = ticker + 1
    }
    df_out_means[ticker_sub,calls] <- out_mean #Output mean for each unique genotype
    ticker_sub = ticker_sub + 1
  }
  df_out[calls] <- count_out #add data to output dataframe with outliers replaced as NA values.
}

write.csv(df_out, row.names=FALSE, file='<YOUR_FILE>')
write.csv(df_out_means, row.names=FALSE, file='<YOUR_FILE>')
