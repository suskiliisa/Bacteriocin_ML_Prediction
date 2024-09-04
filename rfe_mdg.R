
library(doMC)
registerDoMC(cores = 10)
library(caret)
library(foreach)
library(randomForest)


## RECURSIVE FEATURE EXTRACTION 

##--- Data Import ====
# MDGI data
file_n <- "/nobackup/kbxz52/R_scripts/MDGI_featureset.csv"
xtdm <- read.csv(file_n, header = TRUE)

col_names <- names(xtdm)
col_names
target <- "Output"

xtdm$Output[xtdm$Output == 1] <- "yes"
xtdm$Output[xtdm$Output == 0] <- "no"

xtdm$Output <- as.factor(xtdm$Output)

feature_names = col_names[col_names != target]
feature_names

caretFuncs$summary <- twoClassSummary

set.seed(123)

ctrl <- rfeControl(functions = caretFuncs, 
                   method = "cv",
                   number = 10,
                   repeats = 5)

trainctrl <- trainControl(classProbs = TRUE,
                          summaryFunction = twoClassSummary)

# Get results
rfe_fit = rfe(xtdm[, feature_names], xtdm[, target],
              sizes = c(1:231), # t-test = c(1:140), RF = c(1:231)
              rfeControl = ctrl,
              method = "rf",
              ntree = 400,
              metric = "ROC",
              # Additional arguments to train method here
              trControl = trainctrl
)

# Summarize the results
print(rfe_fit)

rfe_fit$results

rfe_fit$bestSubset
# df <- as.data.frame(rfe_fit$variables, rfe_fit$)

rfe_fit$optVariables

file <- "/nobackup/kbxz52/R_outputs/summary_MDG_RF.csv"
write.csv(rfe_fit$results, file, row.names = FALSE)

file <- "/nobackup/kbxz52/R_outputs/optva_MDG_RF.csv"
write.csv(rfe_fit$optVariables, file, row.names = FALSE)

