#library(ROSE) 
library(e1071)
library(caret)
library(ROCR)
#library(pROC)
#library(klaR)
library(mltools)
library(pROC)
library(party)
library(partykit)
library(dplyr)

##--- selecting ADTree features ====

#ADTree reduced dataset - originally 590 variables, reduced to 39 after adtree
#selected and reduced features after ADTree

adt_features <- c("amphipseudo_15","amphipseudo_2", "dist_72", "aac_5", "dipep_203", "comp_16", "dipep_148", "dipep_302", "dipep_68", "dipep_290", "dipep_40", "comp_3", "trans_18", "comp_19", "dist_12", "dist_76", "dist_75", "dipep_198", "dipep_292", "aac_12", "dist_100", "amphipseudo_11", "dist_82", "dist_18","dipep_380", "dipep_297", "dist_10", "dist_27", "dist_23", "dist_58", "dipep_169", "dist_79", "dipep_398", "aac_3", "dipep_294", "dist_38",
"dipep_295", "dipep_274", "dist_91", "Output")

ADT_selected <- df_trainX %>% 
  select(all_of(adt_features))

write.csv(ADT_selected, file = "/path/ADT_selected_features.csv", row.names = FALSE)


# extracting the adtree features from the test data
test_data <-read.table("/path/training_testing_data/test_DATA.txt", header = TRUE, sep = "")

ADT_selected_test <- test_data %>% 
  select(all_of(adt_features))

write.csv(ADT_selected_test, file = "/Users/susannasutki/Library/CloudStorage/OneDrive-DurhamUniversity/DISSERTATION/R_scripts/feature_evaluation/ADT_train_test/testing_ADT_selected_features.csv", row.names = FALSE)

###
org_check <- read.table("/path/train_DATA.txt", header = TRUE, sep = "")
org_check_test <- read.table("/path/test_DATA.txt", header = TRUE, sep = "")
###

##--- Importing ADTree feature sets ====

## WITH ADTREE SELECTED FEATURES

# training
x<-"/path/training_ADT_selected_features.csv"

#testing
y<-"/path/testing_ADT_selected_features.csv"

##--- Preparation ====

# training
file_n <- x

#testing
file_other <- y

## TRAINING DATA
data1 <- read.csv(file_n, header = TRUE)
data1$Output

# column number, ie. var number
cols<-ncol(data1)

# backup training
data_backup<-data1

nrows_training <- nrow(data1)


## TESTING DATA
data_other <- read.csv(file_other, header = TRUE)
data_other$Output

# backup testing
data_other_backup <- data_other

data_other2<-data_other
cols2 <- ncol(data_other)
cols2

# identifies the columns in data_other that have no variability 
which(apply(data_other, 2, var) == 0)



##--- DATA NORMALIZATION ====


# Excluding the output column from the training and testing data for normalization
output_col <- "Output"
data1_features <- data1[, !colnames(data1) %in% output_col]
data_other_features <- data_other[, !colnames(data_other) %in% output_col]

# Normalizing the feature columns based on training data
data_norm <- preProcess(data1_features, method = c("center", "scale"))

# Applying the normalization to both training and testing data
data1_normalized <- predict(data_norm, data1_features)
data_other_normalized <- predict(data_norm, data_other_features)

# Adding the output column back to the normalized data
data1_normalized[output_col] <- data_backup[, output_col]
data_other_normalized[output_col] <- data_other_backup[, output_col]

# Converting the output column to a factor with valid R variable names
data1_normalized$Output <- factor(data1_normalized$Output, levels = c("0", "1"), labels = c("Class0", "Class1"))
data_other_normalized$Output <- factor(data_other_normalized$Output, levels = c("0", "1"), labels = c("Class0", "Class1"))


##--- Model training ====

train<-data1_normalized
test<-data_other_normalized

### TESTING DIFFERENT APPROACH

# Set up train control for cross-validation
fitControl <- trainControl(method = "cv", 
                           number = 5,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE)

# Define the tuning grid for hyperparameters
svmGrid <- expand.grid(sigma = seq(0, 1, 0.1), C = 2^(2:7))

# Set seed for reproducibility
set.seed(123)

# Train the model with hyperparameter tuning
new_tmodel2 <- train(Output ~ ., 
                 data = train, 
                 method = "svmRadial", 
                 trControl = fitControl, 
                 tuneGrid = svmGrid, 
                 metric = "ROC")
new_tmodel2
plot(new_tmodel2)


# final model after training
new_tmodel2$finalModel


variable_importance_svm <- varImp(new_tmodel2)
variable_importance_svm 
plot(varImp(new_tmodel2), main="Variable Importance with ADTree feautre set & SVM")

# predictions for TRAINING data
svm_preds_train <- predict(new_tmodel2, train, type = "raw") # 'type=raw' would return raw values

# prediction on unseen !TEST data
svm_preds_test <- predict(new_tmodel2, newdata = test, type = "raw")


##--- Evaluating the predictive performance ====

# Predict the probabilities for the testing data
probabilities_train <- predict(new_tmodel2, train, type = "prob")
probabilities_test <- predict(new_tmodel2, newdata = test, type = "prob")


# Write the testing probabilities to a CSV file
write.csv(probabilities_test, file = "/path/ADTfs_probability_test_SVM.csv", row.names = FALSE)


# Check column names of the predicted probabilities
colnames(probabilities_train)
colnames(probabilities_test)


# Confusion Matrix for Training Data
conf_matrix_train <- confusionMatrix(svm_preds_train, train$Output)
conf_matrix_train

# Confusion Matrix for Testing Data
conf_matrix_test <- confusionMatrix(svm_preds_test, test$Output)
conf_matrix_test

# ROC Curve and AUC for Training Data
roc_curve_train <- roc(train$Output, probabilities_train[,"Class1"])
plot(roc_curve_train, col = "blue", main = "ROC Curve - Training Data")
auc_train <- auc(roc_curve_train)
print(paste("AUC for Training Data:", auc_train))

# ROC Curve and AUC for Testing Data
roc_curve_test <- roc(test$Output, probabilities_test[,"Class1"])
plot(roc_curve_test, col = "red", main = "ROC Curve - Testing Data")
auc_test <- auc(roc_curve_test)
print(paste("AUC for Testing Data:", auc_test))

# Calculate additional performance metrics for training data
precision_train <- posPredValue(svm_preds_train, train$Output, positive = "Class1")
recall_train <- sensitivity(svm_preds_train, train$Output, positive = "Class1")
f1_score_train <- (2 * precision_train * recall_train) / (precision_train + recall_train)
mcc_train <- mcc(train$Output, svm_preds_train)

# Calculate additional performance metrics for testing data
precision_test <- posPredValue(svm_preds_test, test$Output, positive = "Class1")
recall_test <- sensitivity(svm_preds_test, test$Output, positive = "Class1")
f1_score_test <- (2 * precision_test * recall_test) / (precision_test + recall_test)
mcc_test <- mcc(test$Output, svm_preds_test)

# Print additional metrics
print(paste("Precision for Training Data:", precision_train))
print(paste("Recall for Training Data:", recall_train))
print(paste("F1 Score for Training Data:", f1_score_train))
print(paste("MCC for Training Data:", mcc_train))

print(paste("Precision for Testing Data:", precision_test))
print(paste("Recall for Testing Data:", recall_test))
print(paste("F1 Score for Testing Data:", f1_score_test))
print(paste("MCC for Testing Data:", mcc_test))


# Define a function to draw the confusion matrix
draw_confusion_matrix <- function(cm) {
  layout(matrix(1))
  par(mar=c(2,2,2,2))
  plot(c(123, 345), c(300, 452), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  
  # Create the matrix
  rect(150, 430, 240, 370, col='#287D8EFF')
  text(195, 435, "Class0", cex=1.2)
  rect(250, 430, 340, 370, col='#DCE319FF')
  text(295, 435, "Class1", cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col='#DCE319FF')
  rect(250, 305, 340, 365, col='#287D8EFF')
  text(140, 400, "Class0", cex=1.2, srt=90)
  text(140, 335, "Class1", cex=1.2, srt=90)
  
  # Add capital letter to the top left corner
  text(125, 445, "B", cex=2, font=2)
  
  # Add in the cm results
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='black')
  text(195, 335, res[2], cex=1.6, font=2, col='black')
  text(295, 400, res[3], cex=1.6, font=2, col='black')
  text(295, 335, res[4], cex=1.6, font=2, col='black')
}
draw_confusion_matrix(conf_matrix_test)


