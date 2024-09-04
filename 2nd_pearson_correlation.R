##UPDATET PEARSON CORRELATION##

library(dplyr)
library(GGally)
library(tidyverse)
library(lubridate)
library(ggplot2)
#install.packages("psych")
library(psych)
#install.packages("factoextra")
library(factoextra)
library(cluster)
install.packages("ggfortify")
library(ggfortify)

##--- Importing File ====

file <- read.csv("/path/training_reduced_ORIGINAL.csv", header = TRUE)

col <- ncol(file)
col

flag = c(1:col)
flag[1:col]=0

Rlist = c()

col <- col-1
col


#============== Create a list of columns that we want to remove ==================
for(i in 1:col)
{
  for(j in 1:col)
  {
    if(flag[i]==0)
    {
      col2= file[[i]]
      col1 = file[[j]]
      i
      j
      col2
      col1
      corre<-abs(cor(col1,col2, method = "pearson"))
      res <- c(corre, i, j)
      if(corre >= 0.90 && i!=j)
      {
        flag[j]=1
        Rlist<-c(Rlist,j)
        #Rlist[index]=j
        #index=index+1
      }
    }
  }
}


##--- Remove the duplicates in our list ====

Rlist <- unique(Rlist)
Rlist
Rlist<-sort(Rlist)
Rlist



##--- Removing the correlated list from the original ====

corrFST <- file
corrFST <- corrFST[,-Rlist]

length(corrFST)

write.csv(corrFST, file = "/path/new_validated_pearson_corr_reduced_FS.csv", row.names = FALSE)


file_read <-read.csv("/path/new_validated_pearson_corr_reduced_FS.csv", header = TRUE, sep = ",")
length(file_read)
file_read$Output


