



# Command line arguments ----
#------------------------------------------------#
#                                                #
#          SET UP COMMAND LINE ARGUMENTS         # 
#                                                #
#------------------------------------------------#

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)             
cat(args, sep = "\n")


# Test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file)", call. = FALSE)
} 

# Declare transformation
transformation_id <- as.numeric(args[1])

# Declare model name
model_name <- args[2]

# Declare grid size
grid_size <- as.numeric(args[3])


# Load packages
library("tidymodels")
library("dplyr")
library("stringr")


# Specify transformation used 
transformations_all <- c("ILR", "CLR", "rCLR", "ALR1", "ALR2", "ALR3", "ALR4", "ALR5", "PA", "aSIN", "TSS", "logTSS")
transformation <- transformations_all[transformation_id]







# ML modelling ----
#------------------------------------------------#
#                                                #
#               TIDYMODELS FUNCTION              # 
#                                                #
#------------------------------------------------#

run_tidymodels <- function(p_inputData){
  
  
  # Initial data split ----
  #------------------------------------------------#
  
  # Split the data according to the data_split variable
  data_train <- p_inputData %>% 
    dplyr::filter(data_split == "train") %>% 
    dplyr::select(-one_of("data_split"))
  
  data_test <- p_inputData %>% 
    dplyr::filter(data_split == "test") %>% 
    dplyr::select(-one_of("data_split"))
  
  
  # Create CV sets
  set.seed(1)
  data_cv <- rsample::vfold_cv(data_train, v = 5, repeats = 4, strata = target_var) # NB - allows strata! 
  
  
  
  # Preprocessing steps for the data ----
  #------------------------------------------------#
  
  # Data Preprocessing
  data_preprocessing <- 
    recipes::recipe(target_var ~ ., data = data_train) %>%
    recipes::step_dummy(all_nominal(), -all_outcomes()) %>%
    recipes::step_normalize(all_predictors()) %>%
    recipes::step_zv(all_predictors()) %>%
    recipes::step_nzv(all_predictors())
  
  
  
  # Define the model used ----
  #------------------------------------------------#
  
  # Define the model - no training done yet
  if (model_name == "RF"){
    
    model <- parsnip::rand_forest() %>%
      parsnip::set_args(
        trees = 100,   
        mtry = tune::tune(),
        min_n = tune::tune()) %>%
      set_engine("ranger") %>% 
      set_mode("classification")
    
  }else if (model_name == "LASSO"){
    
    model <- parsnip::logistic_reg() %>%
      parsnip::set_args(
        penalty = tune::tune(),  
        mixture = tune::tune()) %>%
      set_engine("glmnet") %>% 
      set_mode("classification")
  }else if (model_name == "XGB"){
    
    model <- parsnip::boost_tree() %>%
      parsnip::set_args(
        trees = 500, 
        mtry = tune::tune(),
        tree_depth = tune::tune(), 
        learn_rate = tune(),
        min_n = tune::tune(),
        sample_size = tune::tune()) %>%
      set_engine("xgboost") %>% 
      set_mode("classification")
  }
  
  
  # Construct a workflow that combines your recipe and your model
  ml_workflow <- workflows::workflow() %>%
    workflows::add_recipe(data_preprocessing) %>%
    workflows::add_model(model)
  
  
  # Start tuning the model ----
  #------------------------------------------------#
  
  # Tune the models
  set.seed(1)
  tuned_model <- ml_workflow %>%
    tune::tune_grid(resamples = data_cv,
                    grid = grid_size, 
                    metrics = yardstick::metric_set(yardstick::roc_auc))
  
  
  # Return best parameters
  best_params <-
    tuned_model %>%
    tune::select_best(metric = "roc_auc")
  
  
  
  # Evaluate model performance ----
  #------------------------------------------------#
  
  # Describe the best models according to the best parameters gathered earlier
  best_wf <-
    ml_workflow %>%
    tune::finalize_workflow(best_params) 
  
  
  # Fit the models
  model_final <- best_wf %>% 
    parsnip::fit(data_train)
  
  
  # Collect CV metrics
  CV_performance <- tuned_model %>% tune::show_best(n = 1, metric = "roc_auc") 
  
  
  # Collect test metrics
  test_auc <- yardstick::roc_auc_vec(truth = data_test$target_var, 
                                     estimate = predict(model_final, data_test, type = "prob") %>% dplyr::pull(.pred_1), 
                                     event_level = "second")
  
  test_precision <- yardstick::precision_vec(truth = data_test$target_var, 
                                             estimate = predict(model_final, data_test, type = "class") %>% dplyr::pull(.pred_class), 
                                             event_level = "second")
  
  test_sensitivity <- yardstick::sens_vec(truth = data_test$target_var, 
                                          estimate = predict(model_final, data_test, type = "class") %>% dplyr::pull(.pred_class), 
                                          event_level = "second")
  
  test_specificity <- yardstick::spec_vec(truth = data_test$target_var, 
                                          estimate = predict(model_final, data_test, type = "class") %>% dplyr::pull(.pred_class), 
                                          event_level = "second")
  
  test_fscore <- yardstick::f_meas_vec(truth = data_test$target_var, 
                                       estimate = predict(model_final, data_test, type = "class") %>% dplyr::pull(.pred_class), 
                                       event_level = "second")
  
  
  # Output ----
  #------------------------------------------------#
  
  # Collect model performance data
  model_performance <- data.frame(target = target, 
                                  transformation = transformation,
                                  grid_size = grid_size,
                                  algoritm = model_name,
                                  set = c("CV", "test"), 
                                  AUC = c(CV_performance$mean, test_auc),
                                  fscore = c(NA, test_fscore), 
                                  precision = c(NA, test_precision), 
                                  sensitivity = c(NA, test_sensitivity), 
                                  test_specificity = c(NA, test_specificity))
  
  return(model_performance)  
}








# Apply the models to all LOSO setting ----
#------------------------------------------------#
#                                                #
#                  APPLY MODELS                  # 
#                                                #
#------------------------------------------------#

input_folder = "Datasets/Datasets_LOSO"
out_folder = "RData/Results_LOSO/"

for (target in c("BMI", "CRC")){
  
  # Prepare ML data
  metadata_file <- file.path(input_folder, paste("LOSO_", target, "_metadata.rds", sep = ""))
  countdata_file <- file.path(input_folder, paste("LOSO_", target, "_transformedData_", transformation, ".rds", sep = ""))
  
  metadata <- readRDS(metadata_file)
  countdata <- readRDS(countdata_file) %>% 
    tibble::rownames_to_column(var = "Sample_id")
  
  ml_data <- metadata %>% 
    dplyr::left_join(countdata, by = "Sample_id") %>% 
    dplyr::filter(complete.cases(.)) %>% 
    tibble::column_to_rownames(var = "Sample_id") %>% 
    dplyr::mutate(target_var = factor(target_var, levels = c(0, 1)))
  
  
  # Run over all LOSO settings
  datasets_ids <- unique(ml_data$LOSO_group)
  for (j in datasets_ids){
    
    ml_data_run = ml_data %>% 
      dplyr::mutate(data_split = ifelse(LOSO_group == j, "test", "train")) %>% 
      dplyr::select(-one_of("LOSO_group"))
    
    # train the model
    model <- run_tidymodels(p_inputData = ml_data_run) %>% 
      dplyr::mutate(LOSOset = j)
    
    # Save the output
    output_name <- paste(out_folder, 
                         "LOSO", 
                         "_target", target, 
                         "_LOSOset", j, 
                         "_trans", transformation,
                         "_algo", model_name,  
                         "_grid", grid_size, 
                         ".rds", sep = "")
    
    saveRDS(model, output_name)
  }
}



