#install.packages('rJava', type = 'source', INSTALL_opts='--merge-multiarch')
#install.package(c("stringr","rattle", "tibble", "bitops","RWeka", "ROSE", "dplyr"))
#install.packages("rattle")
#install.packages("tibble")
#install.packages("bitops")
#install.packages("RWeka")
#install.packages("ROSE")

library(dplyr)
library(stringr)
#library(GGally)
#library(finalfit)
library(tibble)
library(bitops)
library(rattle)
library(RWeka)
library(rpart)

library(ROSE)


##--- Setting up ====

WPM("refresh-cache")
WPM("install-package", "alternatingDecisionTrees")

WPM("load-package", "alternatingDecisionTrees")
cpath <- "weka/classifiers/trees/ADTree"
ADT <- make_Weka_classifier(cpath)

ADT

WOW(ADT)

##--- Training ====

df_trainX<-read.csv("/path/to/pearson/corr/featu/set/new_validated_pearson_corr_reduced_FS.csv", header = TRUE, sep = ",")

length(df_trainX)

##output balance
summary(as.factor(df_trainX$Output))

##subset of the training data 
subtrainX<-df_trainX
subtrainX

summary(df_trainX)

subtrainX$Output<-as.factor(subtrainX$Output)
subtrainX$Output
str(subtrainX$Output)
#str(df_trainX$Output)

##
## testing ##
##



set.seed(123)
#possiblecValue <- round(seq(from = 1, to = 100, length.out = 200),0)
possiblecValue <- round(seq(from = 5, to = 50, length.out = 55),0)
possiblecValue
numModels <- 50
cValue <- sample(possiblecValue, numModels)
cValue
cValue<-round(cValue)
cValue

pctCorrect <- MAE <- Kappa <- rep(0, numModels)
pctCorrect

for(i in 1:numModels){
  print(i)
  audit.adt <- ADT(Output ~ ., data=subtrainX, control = Weka_control(B =cValue[i]))
  e <- evaluate_Weka_classifier(audit.adt,
                                numFolds = 10, complexity = TRUE,
                                seed = 123, class = TRUE)
  evaluation <- e$details
  evaluation
  pctCorrect[i] <- evaluation["pctCorrect"]
  Kappa[i] <- evaluation["kappa"]
  #MAE[i] <- evaluation["meanAbsoluteError"]
  MAE[i] <- evaluation["rootMeanSquaredError"]
}

ind<-which.min(MAE)
ind
min(MAE)


t <- ADT(Output ~ ., data=subtrainX, control = Weka_control(B =cValue[ind]))
t

summary(t)


library("Rgraphviz")
ff <- tempfile()
write_to_dot(t, ff)
write_to_dot(t, "/path/feature_evaluation/220624_3_adt_tree.dot")

plot(agread(ff)) 

