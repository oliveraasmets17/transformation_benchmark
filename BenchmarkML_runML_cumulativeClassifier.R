

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
file_id <- as.numeric(args[1])



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

run_tidymodels <- function(p_inputData, p_model_name, p_seed, p_calculate_SHAP, p_datasplit){
  
  # Initial data split ----
  #------------------------------------------------#
  
  # Random data split
  set.seed(p_seed)
  data_split <- rsample::initial_split(p_inputData, prop = p_datasplit, strata = target_var) # NB - allows strata! 
  
  # Split the data
  data_train <- rsample::training(data_split)
  data_test <- rsample::testing(data_split)
  
  
  # Create CV sets
  set.seed(p_seed)
  data_cv <- rsample::vfold_cv(data_train, v = 5, repeats = 1, strata = target_var) # NB - allows strata! 
  
  
  
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
  if (p_model_name == "RF"){
    
    model <- parsnip::rand_forest() %>%
      parsnip::set_args(
        trees = 500,    
        mtry = tune::tune(), 
        min_n = tune::tune()) %>%
      set_engine("ranger") %>% 
      set_mode("classification")
    
  } else if (p_model_name == "LASSO"){
    
    model <- parsnip::logistic_reg() %>%
      parsnip::set_args(
        penalty = tune::tune(),  
        mixture = tune::tune()) %>%
      set_engine("glmnet") %>% 
      set_mode("classification")
    
  } else if (p_model_name == "XGB"){
    
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
                                                    grid = 10, 
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
  
  
  
  # Output ----
  #------------------------------------------------#
  
  # Collect model performance data
  model_performance <- data.frame(algoritm = p_model_name,
                                  set = c("CV", "test"), 
                                  AUC = c(CV_performance$mean, test_auc)) 
  
  
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
    
    
    return(list(shap, data_test, model_performance))
    
  } else{
    return(list(NULL, NULL, model_performance))
    
  }
}








# Apply the models to several seeds/LOSO setting ----
#------------------------------------------------#
#                                                #
#                  APPLY MODELS                  # 
#                                                #
#------------------------------------------------#


# Make filelist
filelist <- data.frame(filename = rep(c("EstBMI30metaphlan_2022_PA.csv", "EstAB90metaphlan_2022_PA.csv", 
                                        "EstBMI30metaphlan_2022_PA.csv", "EstAB90metaphlan_2022_PA.csv", 
                                        "EstBMI30metaphlan_2022_PA.csv", "EstAB90metaphlan_2022_PA.csv"), each = 10), 
                       algorithm = rep(c("XGB", "LASSO", "RF"), each = 20), 
                       seed = rep(1:10, 6)) %>% 
  tibble::rowid_to_column(var = "id")

filelist <- data.frame(filename = rep("EstDepressionmetaphlan_2022_PA.csv", 30), 
                       algorithm = rep(c("XGB", "LASSO", "RF"), each = 10), 
                       seed = rep(1:10, 3)) %>% 
  tibble::rowid_to_column(var = "id")


# Define setting used in current run
setup_to_run <- filelist %>% dplyr::filter(id == file_id)


# Run analysis
run_filename = setup_to_run %>% dplyr::pull(filename)
run_algorithm = setup_to_run %>% dplyr::pull(algorithm)
run_seed = setup_to_run %>% dplyr::pull(seed)


# Run stage 1
datasets_path <- "Datasets/Datasets_EstMB/"
stage1_data = read.csv(file.path(datasets_path, run_filename)) %>% 
  tibble::column_to_rownames(var = "sample_id") %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::mutate(target_var = factor(target_var, levels = c(0, 1)))

results_featureimp <- run_tidymodels(p_inputData = stage1_data, 
                                     p_model_name = run_algorithm, 
                                     p_seed = run_seed, 
                                     p_calculate_SHAP = T, 
                                     p_datasplit = 0.5)

# Run stage 2
shap_results <- results_featureimp[[1]]
shap_df <- shap_results %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "skood") %>% 
  tidyr::gather(f_name, value, -skood) %>% 
  dplyr::group_by(f_name) %>% 
  dplyr::summarize(shap_value = mean(abs(value))) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(shap_value)) %>% 
  tibble::rowid_to_column(var = "id")
  
for (subset in c(5, 10, 25, 50, 75, 100, 200)){
  
  top_taxa = shap_df %>% 
    dplyr::filter(id <= subset) %>% 
    dplyr::pull(f_name)
  
  ml_data = results_featureimp[[2]]
  ml_data_subset = ml_data %>% 
    dplyr::select(target_var, all_of(top_taxa)) %>% 
    dplyr::mutate(target_var = factor(target_var, levels = c(0, 1)))
  
  for (seed2 in 1:10){
    results_evaluation <- run_tidymodels(p_inputData = ml_data_subset, 
                                         p_model_name = run_algorithm, 
                                         p_seed = seed2, 
                                         p_calculate_SHAP = FALSE, 
                                         p_datasplit = 0.75)
    
    run_performance <- results_evaluation[[3]] %>% 
      dplyr::mutate(filename = run_filename, 
                    subset_size = subset, 
                    algorithm = run_algorithm,
                    fimp_seed = run_seed, 
                    eval_seed = seed2)
    
    outfile = paste(substr(run_filename, 1, nchar(run_filename)-4), 
                    "_subset", subset, 
                    "_algo", run_algorithm, 
                    "_fimp", run_seed, 
                    "_eval", seed2, 
                    ".rds", sep = "") 
    saveRDS(run_performance, file.path("RData/Results_cumulative", outfile))
  }
}
