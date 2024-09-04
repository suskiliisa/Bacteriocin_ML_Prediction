## Reduced feature set with mean decrease Gini (rf) ##

library(caret)
library(randomForest)


##--- Importing File ====

##header = TRUE --> uses the column headers
file <- read.csv("C:\\Users\\Dell PC\\Desktop\\CORRECT\\new_validated_pearson_corr_reduced_FS.csv", header = TRUE)

##--- Reading the features ====

col <- ncol(file) 
col

file_backup<-file
length(file_backup$Output)

c(col)

file<-file[-c(col)]
file

file<-as.data.frame(scale(file))
length(file)


#adding the output column back
file$Output<-file_backup$Output


##--- Importing Pearson corr feature set back for MGD ====

corrFSTnum <-file

colnames(corrFSTnum)
ncol(corrFSTnum)
corrFSTnum$Output

len_num_after<-length(corrFSTnum)-1

corrFSTnum$Output<-as.factor(corrFSTnum$Output)
str(corrFSTnum$Output)

## model fitting
set.seed(123)

fit_rf <- randomForest(corrFSTnum$Output~., data=corrFSTnum)

vi<-floor(importance(fit_rf))
#vi changed the floor to ceiling
#vi<-ceiling(importance(fit_rf))
length(vi)
vi
#varImp(fit_rf)
plot(vi)
#vi

Rlist<-c()

for(k in 1:len_num_after){
  if(vi[k]==0){
    Rlist<-c(Rlist,k)
    Rlist
  }
}

length(Rlist)

if(length(Rlist)>0){
  corrFSTnum <- corrFSTnum[,-Rlist]
  length(Rlist)
  vi<-vi[-Rlist]
}

ncol(corrFSTnum)

lngt<-ncol(corrFSTnum)

corrFSTnum_vi<-corrFSTnum[,-lngt]

vi_col<-colnames(corrFSTnum_vi)
vi[1]

length(vi)

length(corrFSTnum)

dtfr<-data.frame(vi_col,vi)
# ????? variables including the 'Output' var

write.csv(corrFSTnum, file = "C:\\Users\\Dell PC\\Desktop\\MDGI_featureset_bapres.csv" )

