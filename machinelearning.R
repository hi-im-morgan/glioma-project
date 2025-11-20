install.packages("randomForest")
install.packages("xgboost")

library(caret)
library(randomForest)
library(xgboost)

# Prepare the final dataset
model_data <- training_matrix %>%
  select(expression, cn, mut_flag, essential_binary) %>%
  filter(!is.na(essential_binary)) %>%  # Remove NA targets
  mutate(essential_binary = as.factor(essential_binary))  # Convert to factor for classification

print(paste("Final dataset size:", nrow(model_data)))
print("Class distribution:")
print(table(model_data$essential_binary))

# Set seed for reproducibility
set.seed(123)

# Create train/test split (80/20)
train_index <- createDataPartition(model_data$essential_binary, p = 0.8, list = FALSE)
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

print(paste("Training set size:", nrow(train_data)))
print(paste("Test set size:", nrow(test_data)))


rf_model <- randomForest(
  essential_binary ~ expression + cn + mut_flag,
  data = train_data,
  ntree = 100,
  importance = TRUE
)

print("Random Forest Model Summary:")
print(rf_model)

# Feature importance
importance(rf_model)
varImpPlot(rf_model)

# Logistic Regression (baseline)
logit_model <- glm(essential_binary ~ expression + cn + mut_flag,
                   data = train_data, family = "binomial")

# XGBoost
xgb_model <- train(
  essential_binary ~ expression + cn + mut_flag,
  data = train_data,
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 5)
)

# Predict on test set
rf_predictions <- predict(rf_model, test_data)
logit_predictions <- predict(logit_model, test_data, type = "response")
xgb_predictions <- predict(xgb_model, test_data)

# Convert logistic regression probabilities to classes
logit_predictions_class <- ifelse(logit_predictions > 0.5, 1, 0)

# Create confusion matrices
rf_cm <- confusionMatrix(rf_predictions, test_data$essential_binary)
logit_cm <- confusionMatrix(as.factor(logit_predictions_class), test_data$essential_binary)
xgb_cm <- confusionMatrix(xgb_predictions, test_data$essential_binary)

print("Random Forest Performance:")
print(rf_cm)
print("Logistic Regression Performance:")
print(logit_cm)
print("XGBoost Performance:")
print(xgb_cm)

# Collect performance metrics
performance_comparison <- data.frame(
  Model = c("Random Forest", "Logistic Regression", "XGBoost"),
  Accuracy = c(rf_cm$overall["Accuracy"], 
               logit_cm$overall["Accuracy"],
               xgb_cm$overall["Accuracy"]),
  Sensitivity = c(rf_cm$byClass["Sensitivity"],
                  logit_cm$byClass["Sensitivity"], 
                  xgb_cm$byClass["Sensitivity"]),
  Specificity = c(rf_cm$byClass["Specificity"],
                  logit_cm$byClass["Specificity"],
                  xgb_cm$byClass["Specificity"])
)

print("Model Performance Comparison:")
print(performance_comparison)

# Detailed feature importance
feature_importance <- importance(rf_model)
print("Feature Importance (Random Forest):")
print(feature_importance)

# Example: Predict essentiality for new gene-cell line pairs
new_predictions <- predict(rf_model, test_data, type = "prob")

# Add predictions to test data
test_results <- test_data %>%
  mutate(
    predicted_prob = new_predictions[,2],  # Probability of being essential
    predicted_class = rf_predictions
  )

# View results
head(test_results)

#seeing genes with a predicted essentiality probability higher than the mean
threshold <- mean(test_results$predicted_prob)
high_prob_genes <- test_results %>%
  filter(predicted_prob > threshold)


ranked_predictions <- test_results %>%
  arrange(desc(predicted_prob))