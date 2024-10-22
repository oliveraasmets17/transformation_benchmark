

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("tidyverse")
library("stringr")







# Within-study setting ----
#------------------------------------------------#
#                                                #
#              WITHIN-STUDY SETTING              # 
#                                                #
#------------------------------------------------#

# Read and gather data
#-----------------------------------#
collect_WS_data <- FALSE
if (collect_WS_data == TRUE){
  WS_result_df <- data.frame()
  WS_fimp_df <- data.frame()
  
  filelist <- list.files("/Bioinf/BenchmarkML/Rev_results_WS/")
  for (i in filelist){
    i_split = str_split(i, pattern = "_")
    run_result = readRDS(file.path("/Bioinf/BenchmarkML/Rev_results_WS/", i))
    
    run_performance = run_result[[2]]
    run_fimp = run_result[[1]] %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "id") %>% 
      tidyr::gather(f_name, value, -id) %>% 
      dplyr::group_by(f_name) %>% 
      dplyr::summarize(shap_value = mean(abs(value))) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(study = paste(i_split[[1]][2], i_split[[1]][3], sep = "_"), 
                    transformation = substring(i_split[[1]][4], 6), 
                    algorithm = ifelse(str_detect(i, "algoLASSO") == TRUE, "ENET", ifelse(str_detect(i, "algoRF") == TRUE, "RF", "XGB")),
                    fold = substring(i_split[[1]][7], 5))
    
    WS_result_df = dplyr::bind_rows(WS_result_df, run_performance)
    WS_fimp_df = dplyr::bind_rows(WS_fimp_df, run_fimp)
  }
  saveRDS(WS_result_df, "WS_result_df.rds")
  saveRDS(WS_fimp_df, "WS_fimp_df.rds")
} else{
  WS_result_df <- readRDS("WS_result_df.rds")
  WS_fimp_df <- readRDS("WS_fimp_df.rds")
}




# Check availability of the models
#-----------------------------------#
length(unique(WS_result_df$study))
length(unique(WS_result_df$transformation))
length(unique(WS_result_df$seed))

table(WS_result_df$algoritm, WS_result_df$seed)
table(WS_result_df$algoritm)








# Within-study setting (rarefied) ----
#------------------------------------------------#
#                                                #
#              WITHIN-STUDY SETTING              # 
#                                                #
#------------------------------------------------#

# Read and gather data
#-----------------------------------#
collect_WS_R_data <- FALSE
if (collect_WS_R_data == TRUE){
  WS_R_result_df <- data.frame()
  WS_R_fimp_df <- data.frame()
  
  filelist <- list.files("/Bioinf/BenchmarkML/Rev_results_WS_R/")
  for (i in filelist){
    i_split = str_split(i, pattern = "_")
    run_result = readRDS(file.path("/Bioinf/BenchmarkML/Rev_results_WS_R/", i))
    
    run_performance = run_result[[2]]
    run_fimp = run_result[[1]] %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "id") %>% 
      tidyr::gather(f_name, value, -id) %>% 
      dplyr::group_by(f_name) %>% 
      dplyr::summarize(shap_value = mean(abs(value))) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(study = paste(i_split[[1]][2], i_split[[1]][3], sep = "_"), 
                    transformation = substring(i_split[[1]][4], 6), 
                    algorithm = ifelse(str_detect(i, "algoLASSO") == TRUE, "ENET", ifelse(str_detect(i, "algoRF") == TRUE, "RF", "XGB")),
                    fold = substring(i_split[[1]][7], 5))
    
    WS_R_result_df = dplyr::bind_rows(WS_R_result_df, run_performance)
    WS_R_fimp_df = dplyr::bind_rows(WS_R_fimp_df, run_fimp)
  }
  saveRDS(WS_R_result_df, "WS_R_result_df.rds")
  saveRDS(WS_R_fimp_df, "WS_R_fimp_df.rds")
} else{
  WS_R_result_df <- readRDS("WS_R_result_df.rds")
  WS_R_fimp_df <- readRDS("WS_R_fimp_df.rds")
}




# Check availability of the models
#-----------------------------------#
length(unique(WS_R_result_df$study))
length(unique(WS_R_result_df$transformation))
length(unique(WS_R_result_df$seed))

table(WS_R_result_df$algoritm, WS_R_result_df$seed)
table(WS_R_result_df$algoritm)







# Within-study setting (EstMB) ----
#------------------------------------------------#
#                                                #
#              WITHIN-STUDY SETTING              # 
#                                                #
#------------------------------------------------#

# Read and gather data
#-----------------------------------#
collect_EstMB_data <- FALSE
if (collect_EstMB_data == TRUE){
  EstMB_result_df <- data.frame()
  EstMB_fimp_df <- data.frame()
  
  filelist <- list.files("/Bioinf/BenchmarkML/Rev_results_EstMB/")
  for (i in filelist){
    i_split = str_split(i, pattern = "_")
    run_result = readRDS(file.path("/Bioinf/BenchmarkML/Rev_results_EstMB/", i))
    
    run_performance = run_result[[2]]
    run_fimp = run_result[[1]] %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "id") %>% 
      tidyr::gather(f_name, value, -id) %>% 
      dplyr::group_by(f_name) %>% 
      dplyr::summarize(shap_value = mean(abs(value))) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(study = paste(i_split[[1]][2], i_split[[1]][3], sep = "_"), 
                    transformation = substring(i_split[[1]][4], 6), 
                    algorithm = ifelse(str_detect(i, "algoLASSO") == TRUE, "ENET", ifelse(str_detect(i, "algoRF") == TRUE, "RF", "XGB")),
                    fold = substring(i_split[[1]][7], 5))
    
    EstMB_result_df = dplyr::bind_rows(EstMB_result_df, run_performance)
    EstMB_fimp_df = dplyr::bind_rows(EstMB_fimp_df, run_fimp)
  }
  saveRDS(EstMB_result_df, "EstMB_result_df.rds")
  saveRDS(EstMB_fimp_df, "EstMB_fimp_df.rds")
} else{
  EstMB_result_df <- readRDS("EstMB_result_df.rds")
  EstMB_fimp_df <- readRDS("EstMB_fimp_df.rds")
}




# Check availability of the models
#-----------------------------------#
length(unique(EstMB_result_df$study))
length(unique(EstMB_result_df$transformation))
length(unique(EstMB_result_df$seed))

table(EstMB_result_df$algoritm, EstMB_result_df$seed)
table(EstMB_result_df$algoritm)







# Within-study setting (EstMB rarefied) ----
#------------------------------------------------#
#                                                #
#              WITHIN-STUDY SETTING              # 
#                                                #
#------------------------------------------------#

# Read and gather data
#-----------------------------------#
collect_EstMB_R_data <- FALSE
if (collect_EstMB_R_data == TRUE){
  EstMB_R_result_df <- data.frame()
  EstMB_R_fimp_df <- data.frame()
  
  filelist <- list.files("/Bioinf/BenchmarkML/Rev_results_EstMB_R/")
  for (i in filelist){
    i_split = str_split(i, pattern = "_")
    run_result = readRDS(file.path("/Bioinf/BenchmarkML/Rev_results_EstMB_R/", i))
    
    run_performance = run_result[[2]]
    run_fimp = run_result[[1]] %>% 
      as.data.frame() %>% 
      tibble::rownames_to_column(var = "id") %>% 
      tidyr::gather(f_name, value, -id) %>% 
      dplyr::group_by(f_name) %>% 
      dplyr::summarize(shap_value = mean(abs(value))) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(study = paste(i_split[[1]][2], i_split[[1]][3], sep = "_"), 
                    transformation = substring(i_split[[1]][4], 6), 
                    algorithm = ifelse(str_detect(i, "algoLASSO") == TRUE, "ENET", ifelse(str_detect(i, "algoRF") == TRUE, "RF", "XGB")),
                    fold = substring(i_split[[1]][7], 5))
    
    EstMB_R_result_df = dplyr::bind_rows(EstMB_R_result_df, run_performance)
    EstMB_R_fimp_df = dplyr::bind_rows(EstMB_R_fimp_df, run_fimp)
  }
  saveRDS(EstMB_R_result_df, "EstMB_R_result_df.rds")
  saveRDS(EstMB_R_fimp_df, "EstMB_R_fimp_df.rds")
} else{
  EstMB_R_result_df <- readRDS("EstMB_R_result_df.rds")
  EstMB_R_fimp_df <- readRDS("EstMB_R_fimp_df.rds")
}




# Check availability of the models
#-----------------------------------#
length(unique(EstMB_R_result_df$study))
length(unique(EstMB_R_result_df$transformation))
length(unique(EstMB_R_result_df$seed))

table(EstMB_R_result_df$algoritm, EstMB_R_result_df$seed)
table(EstMB_R_result_df$algoritm)







# Within-study setting (filtered) ----
#------------------------------------------------#
#                                                #
#              WITHIN-STUDY SETTING              # 
#                                                #
#------------------------------------------------#

# Read and gather data
#-----------------------------------#
collect_WS_F_data <- FALSE
if (collect_WS_F_data == TRUE){
  WS_F_result_df <- data.frame()

  filelist <- list.files("/Bioinf/BenchmarkML/Results_WS_F")
  for (i in filelist){
    run_result = readRDS(file.path("/Bioinf/BenchmarkML/Results_WS_F", i))
    
    run_performance = run_result[[2]] %>% 
      dplyr::mutate(filtering = ifelse(length(str_split(study, "_")[[1]]) ==3, str_split(study, "_")[[1]][3], "full"), 
                    study = paste(str_split(study, "_")[[1]][1], "_", str_split(study, "_")[[1]][2], sep = ""))
    
    WS_F_result_df = dplyr::bind_rows(WS_F_result_df, run_performance)
  }
  saveRDS(WS_F_result_df, "WS_F_result_df.rds")
} else{
  WS_F_result_df <- readRDS("WS_F_result_df.rds")
}







# Within-study setting (rarefied and filtered) ----
#------------------------------------------------#
#                                                #
#              WITHIN-STUDY SETTING              # 
#                                                #
#------------------------------------------------#

# Read and gather data
#-----------------------------------#
collect_WS_FR_data <- FALSE
if (collect_WS_FR_data == TRUE){
  WS_FR_result_df <- data.frame()
  
  filelist <- list.files("/Bioinf/BenchmarkML/Results_WS_FR")
  for (i in filelist){
    run_result = readRDS(file.path("/Bioinf/BenchmarkML/Results_WS_FR", i))
    
    run_performance = run_result[[2]] %>% 
      dplyr::mutate(filtering = ifelse(length(str_split(study, "_")[[1]]) ==3, str_split(study, "_")[[1]][3], "full"), 
                    study = paste(str_split(study, "_")[[1]][1], "_", str_split(study, "_")[[1]][2], sep = ""))
    
    WS_FR_result_df = dplyr::bind_rows(WS_FR_result_df, run_performance)
  }
  saveRDS(WS_FR_result_df, "WS_FR_result_df.rds")
} else{
  WS_FR_result_df <- readRDS("WS_FR_result_df.rds")
}








# Within-study setting (grid analysis) ----
#------------------------------------------------#
#                                                #
#              WITHIN-STUDY SETTING              # 
#                                                #
#------------------------------------------------#

# Read and gather data
#-----------------------------------#
collect_grid_data <- FALSE
if (collect_grid_data == TRUE){
  grid_result_df <- data.frame()
  
  filelist <- list.files("/Bioinf/BenchmarkML/ML_grid/")
  for (i in filelist){
    run_result = readRDS(file.path("/Bioinf/BenchmarkML/ML_grid/", i))
    
    grid_result_df = dplyr::bind_rows(grid_result_df, run_result)
  }
  saveRDS(grid_result_df, "grid_result_df.rds")
} else{
  grid_result_df <- readRDS("grid_result_df.rds")
}








# LOSO setting ----
#------------------------------------------------#
#                                                #
#                   LOSO SETTING                 # 
#                                                #
#------------------------------------------------#

# Read and gather data
#-----------------------------------#
collect_LOSO_data <- FALSE
if (collect_LOSO_data == TRUE){
  LOSO_result_df <- data.frame()
  
  filelist <- list.files("/Bioinf/BenchmarkML/Rev_results_LOSO/")
  for (i in filelist){
    run_result = readRDS(file.path("/Bioinf/BenchmarkML/Rev_results_LOSO/", i))
    
    LOSO_result_df = dplyr::bind_rows(LOSO_result_df, run_result)
  }
  
  # Save
  saveRDS(LOSO_result_df, "LOSO_result_df.rds")
} else{
  LOSO_result_df <- readRDS("LOSO_result_df.rds")
}











# ALR sensitivity analysis ----
#------------------------------------------------#
#                                                #
#                ALR SENSITIVITY                 # 
#                                                #
#------------------------------------------------#

collect_ALR_data <- FALSE
if (collect_ALR_data == TRUE){
  
  ALRsensitivity_results <- data.frame()
  ALRsens_filelist <- list.files("/Bioinf/BenchmarkML/Rev_results_ALRsensitivity/")
  
  for (i in ALRsens_filelist){
    
    transformation = substring(str_split(i, "_")[[1]][3], 6)
    
    run_data = readRDS(file = file.path("/Bioinf/BenchmarkML/Rev_results_ALRsensitivity/", i)) %>% 
      dplyr::mutate(transformation)
    
    ALRsensitivity_results = dplyr::bind_rows(ALRsensitivity_results, run_data)
  }
  # Save
  saveRDS(ALRsensitivity_results, "ALRsensitivity_results.rds")
} else{
  ALRsensitivity_results <- readRDS("ALRsensitivity_results.rds")
}











# Cumulative classifier ----
#------------------------------------------------#
#                                                #
#            CUMULATIVE CLASSIFIER               # 
#                                                #
#------------------------------------------------#

run_cumulative <- FALSE
if (run_cumulative == TRUE){
  
  cumulative_filelist <- list.files("/Bioinf/BenchmarkML/Rev_results_cumulative/")
  cumulative_results = data.frame()
  
  for (i in cumulative_filelist){
    
    run_data <- readRDS(file.path("/Bioinf/BenchmarkML/Rev_results_cumulative/", i)) 
    
    cumulative_results = dplyr::bind_rows(cumulative_results, run_data)
  }
  # Save
  saveRDS(cumulative_results, "cumulative_results.rds")
} else{
  cumulative_results <- readRDS("cumulative_results.rds")
}

