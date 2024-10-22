

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

# Declare input dataset folder
file_id <- args[1]

# Declare input dataset folder
filelist <- args[2]

# Declare input dataset folder
input_folder <- args[3]

# Declare model name
model_name <- args[4]

# Declare grid size
grid_size <- as.numeric(args[5])

# Declare whether to calculate SHAP
calculate_SHAP <- args[6]

# Declare seed start
seed_start = as.numeric(args[7])

# Declare seed start
seed_end = as.numeric(args[8])

# Declare output folder
output_folder <- args[9]


# Load packages
library("tidymodels")
library("dplyr")
library("stringr")
library("readr")





# ML modelling ----
#------------------------------------------------#
#                                                #
#               TIDYMODELS FUNCTION              # 
#                                                #
#------------------------------------------------#

run_tidymodels <- function(p_inputData, p_seed, p_transformation, p_studyname, p_calculate_SHAP){
  
  # Initial data split ----
  #------------------------------------------------#
  
  # Random data split
  set.seed(p_seed)
  data_split <- rsample::initial_split(p_inputData, prop = .75, strata = target_var)
  
  # Split the data
  data_train <- rsample::training(data_split)
  data_test <- rsample::testing(data_split)
  
  
  # Create CV sets
  set.seed(p_seed)
  data_cv <- rsample::vfold_cv(data_train, v = 5, repeats = 1, strata = target_var)
  
  
  
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
        trees = 500,   
        mtry = tune::tune(),  
        min_n = tune::tune()) %>%
      set_engine("ranger") %>%
      set_mode("classification")
    
  } else if (model_name == "LASSO"){
    
    model <- parsnip::logistic_reg() %>%
      parsnip::set_args(
        penalty = tune::tune(),  
        mixture = tune::tune()) %>%
      set_engine("glmnet") %>% 
      set_mode("classification")
    
  } else if (model_name == "XGB"){
    
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
  set.seed(p_seed)
  tuned_model <- suppressMessages(ml_workflow %>%
                                    tune::tune_grid(resamples = data_cv,
                                                    grid = grid_size, 
                                                    metrics = yardstick::metric_set(yardstick::roc_auc)))
  
  
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
  model_performance <- data.frame(study = p_studyname, 
                                  transformation = p_transformation,
                                  grid_size = grid_size,
                                  seed = p_seed,
                                  algoritm = model_name,
                                  set = c("CV", "test"), 
                                  AUC = c(CV_performance$mean, test_auc),
                                  fscore = c(NA, test_fscore), 
                                  precision = c(NA, test_precision), 
                                  sensitivity = c(NA, test_sensitivity), 
                                  test_specificity = c(NA, test_specificity))
  
  
  # Explain the model ----
  #------------------------------------------------#
  
  if (p_calculate_SHAP == "TRUE"){
    
    library("fastshap")
    set.seed(p_seed)
    
    # Prediction wrapper
    pfun <- function(object, newdata) {  # needs to return a numeric vector
      predict(object, new_data = newdata, type = "prob")$.pred_1
    }
    
    shap <- fastshap::explain(object = model_final, 
                              X = data_train %>% dplyr::select(-one_of("target_var")), 
                              newdata = data_test %>% dplyr::select(-one_of("target_var")), 
                              nsim = 10, pred_wrapper = pfun)
    
    rownames(shap) <- rownames(data_test)
    
    
    return(list(shap, model_performance))
    
  } else{
    return(list(NULL, model_performance))
    
  }
}








# Apply the models to several seeds/LOSO setting ----
#------------------------------------------------#
#                                                #
#                  APPLY MODELS                  # 
#                                                #
#------------------------------------------------#

# Run over all files for this run
file_name <- read_table(paste(filelist, ".txt", sep = ""), col_names = FALSE) %>% 
  dplyr::filter(X1 == file_id) %>% 
  dplyr::pull(X2)


for (i in file_name){
  
  current_file_path <- file.path("Datasets", input_folder, i)
  
  # Metadata
  study_name <- str_split(i, "_")[[1]][1]
  study_year <- str_split(i, "_")[[1]][2]
  
  if (str_detect(i, "filtered")){
    transformation <- str_split(i, "_")[[1]][3]
  } else if(str_detect(i, "rarefied")){
    transformation <- str_split(i, "_")[[1]][3]
  }else{
    transformation <- substr(str_split(i, "_")[[1]][3], 1, nchar(str_split(i, "_")[[1]][3]) - 4)
  }
  
  if (str_detect(i, "filtered")){
    
    filtering_p <- substr(str_split(i, "_")[[1]][4], 1, nchar(str_split(i, "_")[[1]][4]) - 4)
    
    study_name_final <- paste(study_name, "_", study_year, "_", filtering_p, sep = "")
  } else{
    study_name_final <- paste(study_name, "_", study_year, sep = "")
  }
  
  print(study_name_final)
  
  # Read count matrix used as input (genes, species etc)
  ml_data <- read_csv(current_file_path) %>% 
    tibble::column_to_rownames(var = "sample_id") %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::mutate(target_var = factor(target_var, levels = c(0, 1)))
  
  
  # Run over all seeds
  for (j in seed_start:seed_end){
    
    model = run_tidymodels(p_studyname = study_name_final, 
                           p_inputData = ml_data, 
                           p_seed = j,
                           p_transformation = transformation, 
                           p_calculate_SHAP = calculate_SHAP)
    
    # Save output
    output_name = paste("WithinStudy_", 
                        study_name_final, 
                        "_trans", transformation,
                        "_algo", model_name,  
                        "_grid", grid_size, 
                        "_seed", j, 
                        "_SHAP", calculate_SHAP, 
                        ".rds", sep = "")
    
    saveRDS(model, file.path("RData", output_folder, output_name))
    
  } # End for seed
} # End for file




