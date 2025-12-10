# Install and load required packages
if (!require("randomForest")) install.packages("randomForest")
if (!require("caret")) install.packages("caret")
if (!require("dplyr")) install.packages("dplyr")

library(randomForest)
library(caret)
library(dplyr)

# Function to prepare model data
prepare_model_data <- function(data) {
  # Select and clean the data
  model_data <- data %>%
    select(expression, cn, mut_flag, essential_binary) %>%
    filter(!is.na(essential_binary)) %>%
    mutate(essential_binary = as.factor(essential_binary))
  
  # Print dataset information
  cat("Final dataset size:", nrow(model_data), "\n")
  cat("Class distribution:\n")
  print(table(model_data$essential_binary))
  
  return(model_data)
}

# Function to split data into train/test sets
create_train_test_split <- function(data, train_ratio = 0.8) {
  set.seed(123)  # For reproducibility
  
  train_index <- createDataPartition(
    data$essential_binary, 
    p = train_ratio, 
    list = FALSE
  )
  
  train_data <- data[train_index, ]
  test_data <- data[-train_index, ]
  
  cat("Training set size:", nrow(train_data), "\n")
  cat("Test set size:", nrow(test_data), "\n")
  
  return(list(train = train_data, test = test_data))
}

# Function to train random forest model
train_random_forest <- function(train_data, n_trees = 100) {
  formula <- essential_binary ~ expression + cn + mut_flag
  
  rf_model <- randomForest(
    formula = formula,
    data = train_data,
    ntree = n_trees,
    importance = TRUE,
    na.action = na.omit
  )
  
  cat("Random Forest Model Summary:\n")
  print(rf_model)
  
  return(rf_model)
}

# Function to evaluate model performance
evaluate_model <- function(model, test_data) {
  # Make predictions
  predictions <- predict(model, test_data)
  
  # Create confusion matrix
  cm <- confusionMatrix(predictions, test_data$essential_binary)
  
  # Extract performance metrics
  performance <- data.frame(
    Accuracy = cm$overall["Accuracy"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    Kappa = cm$overall["Kappa"]
  )
  
  return(list(
    predictions = predictions,
    confusion_matrix = cm,
    performance = performance
  ))
}

# Function to display feature importance
display_feature_importance <- function(model) {
  cat("\nFeature Importance (Random Forest):\n")
  imp <- importance(model)
  print(imp)
  
  # Create variable importance plot
  varImpPlot(
    model,
    main = "Random Forest - Feature Importance",
    col = "steelblue",
    pch = 19
  )
  
  return(imp)
}

# Main execution
main <- function() {
  # Check if training_matrix exists
  if (!exists("training_matrix")) {
    stop("Error: 'training_matrix' not found. Please load your data first.")
  }
  
  # Step 1: Prepare data
  cat("=== Step 1: Preparing Data ===\n")
  model_data <- prepare_model_data(training_matrix)
  
  # Step 2: Create train/test split
  cat("\n=== Step 2: Creating Train/Test Split ===\n")
  split_data <- create_train_test_split(model_data)
  train_data <- split_data$train
  test_data <- split_data$test
  
  # Step 3: Train model
  cat("\n=== Step 3: Training Random Forest Model ===\n")
  rf_model <- train_random_forest(train_data, n_trees = 100)
  
  # Step 4: Evaluate model
  cat("\n=== Step 4: Evaluating Model ===\n")
  results <- evaluate_model(rf_model, test_data)
  
  # Print performance metrics
  cat("\nModel Performance:\n")
  print(results$performance)
  
  # Print confusion matrix
  cat("\nConfusion Matrix:\n")
  print(results$confusion_matrix$table)
  
  # Step 5: Feature importance
  cat("\n=== Step 5: Feature Importance Analysis ===\n")
  feature_importance <- display_feature_importance(rf_model)
  
  # Return all results
  return(list(
    model = rf_model,
    train_data = train_data,
    test_data = test_data,
    predictions = results$predictions,
    performance = results$performance,
    feature_importance = feature_importance
  ))
}



plot_simple_confusion_matrix <- function(actual, predicted, title = "Confusion Matrix") {
  library(caret)
  
  cm <- confusionMatrix(predicted, actual)
  cm_table <- as.data.frame(cm$table)
  
  ggplot(cm_table, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Freq), color = "white") +
    geom_text(aes(label = Freq), color = "white", size = 6, fontface = "bold") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(title = title,
         subtitle = paste("Accuracy:", round(cm$overall["Accuracy"] * 100, 1), "%"),
         x = "Actual",
         y = "Predicted") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 10)) +
    coord_fixed()
}


plot_simple_confusion_matrix(
  actual = test_data$essential_binary,
  predicted = rf_predictions,
  title = "Random Forest Performance"
)

#major class imbalance!
