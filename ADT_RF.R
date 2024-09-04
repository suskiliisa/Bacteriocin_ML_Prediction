#library(ROSE) 
library(e1071)
library(caret)
library(ROCR)
#library(pROC)
#library(klaR)
library(mltools)
library(pROC)
#install.packages("plyr")
library(plyr)
library(randomForest)
library(data.table)
library(ggplot2)

##--- Data Import ====

# training
x<-"/path/ADT_train_test/training_ADT_selected_features.csv"

#testing
y<-"/path/ADT_train_test/testing_ADT_selected_features.csv"

##--- Data Preparation ==== 

#training
file_n<-x

#testing
file_other<-y

##reading files

# training
data1 <- read.csv(file_n, header = TRUE)
####changing/testing for corresponfing var names!
# now matches with the predicted features
data1 <- data1 %>% rename(trans_18 = tran_18)
data1$Output
#testing
data_other <- read.csv(file_other, header = TRUE)
####changing/testing for corresponfing var names!
# now matches with the predicted features
data_other <- data_other %>% rename(trans_18 = tran_18)

cols<-ncol(data1)


##creating back-ups

#training
data_backup <- data1
#testing
data_other_backup <- data_other
data_other2 <- data_other


nrows_training <- nrow(data1)
nrows_testing<-nrow(data_other)

cols <- ncol(data1)
cols2<-ncol(data_other)


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


#data_norm <- readRDS("/Users/susannasutki/Desktop/to_hamilton/upd_ADT_data_norm.rds")
#data_norm_old <- readRDS("/Users/susannasutki/Desktop/to_hamilton/ADT_data_norm.rds")
#model_1 <- readRDS("/Users/susannasutki/Desktop/to_hamilton/FINAL_ADT_RF_model.rds")


##--- MODEL TRAINING ====

train<-data1_normalized
test<-data_other_normalized

counts_tain <- table(train$Output)
# Class0 Class1 
# 1588   1559 

counts_test <- table(test$Output)
# Class0 Class1 
# 379    408 

set.seed(123)

# Set up train control for cross-validation
fitControl <- trainControl(method = "cv", 
                           number = 5,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE)


# Set the number of trees (ntree) values to try
ntree_values <- c(400, 500)
results <- list()

for (ntree in ntree_values) {
  # Train the model with hyperparameter tuning for each ntree value
  rf_tmodel <- train(Output ~ ., 
                     data = train, 
                     method = "rf", 
                     trControl = fitControl, 
                     tuneGrid = expand.grid(mtry = c(6,5)), 
                     metric = "ROC",
                     ntree = ntree)
  results[[paste("ntree", ntree, sep = "_")]] <- rf_tmodel
}

# Check results for each ntree value
results
## best performance with mtry=5 & ntree 400!!!!!!!!!!!! checked 0608, this is chocen 

best_ntree <- 400

opt_rf_tmodel <- train(Output ~ ., 
                   data = train, 
                   method = "rf", 
                   trControl = fitControl, 
                   tuneGrid = expand.grid(mtry = c(6,5)), 
                   metric = "ROC",
                   ntree = best_ntree)

opt_rf_tmodel

### SAVING the MODEL
saveRDS(opt_rf_tmodel, "/path/upd_FINAL_ADT_RF_model.rds")

### SAVING THE NORMALISATION
#saveRDS(data_norm, "/path/ADT_data_norm.rds")

###testing
loaded_model <- readRDS("/path/FINAL_ADT_RF_model.rds")
rf_preds_train2 <- predict(loaded_model, train, type = "raw") 
###testing ends

# Loading the saved model
opt_rf_tmodel <- readRDS("/path/upd_FINAL_ADT_RF_model.rds")

variable_importance_rf <- varImp(opt_rf_tmodel)
variable_importance_rf
plot(varImp(opt_rf_tmodel), main="Variable importance with ADTree feature set & Random Forest")

# Convert to a data frame and arrange by importance
importance_df <- as.data.frame(variable_importance_rf$importance)
importance_df$Variable <- rownames(importance_df)
importance_df <- importance_df[order(-importance_df$Overall), ]



# Plot using ggplot2
ggplot(importance_df, aes(x = reorder(Variable, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "#287D8EFF") +
  coord_flip() + # Flip coordinates for better readability
  theme_minimal() +
  labs(title = "Variable Importance with ADTree feature set & Random Forest",
       x = "Variables",
       y = "Importance") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"))



# predictions for TRAINING data
rf_preds_train <- predict(opt_rf_tmodel, train, type = "raw") # 'type=raw' would return raw values

# prediction on unseen !TEST data
rf_preds_test <- predict(opt_rf_tmodel, newdata = test, type = "raw")

### testing
fwrite(rf_probabilities_test, "/path/TESTADTfs_probability_new_RF.csv")

rf_preds_test <- as.data.frame(rf_preds_test)

# If you want to print the predicted classes instead of probabilities
fwrite(rf_preds_test, "/path/TESTADTfs_predictions_new_RF.csv")



##--- Evaluating the predictive performance ====

# Predict the probabilities for the testing data
rf_probabilities_train <- predict(opt_rf_tmodel, train, type = "prob")
rf_probabilities_test <- predict(opt_rf_tmodel, newdata = test, type = "prob")

# Write the testing probabilities to a CSV file
write.csv(rf_probabilities_test, file = "/path/ADTfs_probability_test_RF.csv", row.names = FALSE)


# Check column names of the predicted probabilities
colnames(rf_probabilities_train)
colnames(rf_probabilities_test)


# Confusion Matrix for Training Data
conf_matrix_train <- confusionMatrix(rf_preds_train, train$Output)
conf_matrix_train

# Confusion Matrix for Testing Data
conf_matrix_test <- confusionMatrix(rf_preds_test, test$Output)
conf_matrix_test

# ROC Curve and AUC for Training Data
roc_curve_train <- roc(train$Output, rf_probabilities_train[,"Class1"])
plot(roc_curve_train, col = "blue", main = "ROC Curve - Training Data")
auc_train <- auc(roc_curve_train)
print(paste("AUC for Training Data:", auc_train))

# ROC Curve and AUC for Testing Data
roc_curve_test <- roc(test$Output, rf_probabilities_test[,"Class1"])
plot(roc_curve_test, col = "red", main = "ROC Curve - Testing Data")
auc_test <- auc(roc_curve_test)
print(paste("AUC for Testing Data:", auc_test))

# Calculate additional performance metrics for training data
precision_train <- posPredValue(rf_preds_train, train$Output, positive = "Class1")
recall_train <- sensitivity(rf_preds_train, train$Output, positive = "Class1")
f1_score_train <- (2 * precision_train * recall_train) / (precision_train + recall_train)
mcc_train <- mcc(train$Output, rf_preds_train)

# Calculate additional performance metrics for testing data
precision_test <- posPredValue(rf_preds_test, test$Output, positive = "Class1")
recall_test <- sensitivity(rf_preds_test, test$Output, positive = "Class1")
f1_score_test <- (2 * precision_test * recall_test) / (precision_test + recall_test)
mcc_test <- mcc(test$Output, rf_preds_test)

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
  text(125, 445, "A", cex=2, font=2)
  
  # Add in the cm results
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='black')
  text(195, 335, res[2], cex=1.6, font=2, col='black')
  text(295, 400, res[3], cex=1.6, font=2, col='black')
  text(295, 335, res[4], cex=1.6, font=2, col='black')
}
draw_confusion_matrix(conf_matrix_test)

