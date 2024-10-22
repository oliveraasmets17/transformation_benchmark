

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------#

# Load packages
library("curatedMetagenomicData")
library("tidyverse")
library("mia")
library("SummarizedExperiment")
library("compositions")






# Define the function for making data transformations ----
#------------------------------------------------#
#                                                #
#     FUNCTION TO MAKE DATA TRANSFORMATIONS      # 
#                                                #
#------------------------------------------------#

transform_countData <- function(p_metadata, p_studyName, p_rarefy = FALSE, p_filterTop = FALSE, p_filterTopPercentage = 1){
  
  # Return the assays from the curatedMetagenomicData package - compliant with the metadata
  metagenomics_object <- suppressMessages(p_metadata %>% returnSamples("relative_abundance", counts = TRUE, rownames = "long"))
  
  # Filter taxa by prevalence where needed
  if (p_filterTop == TRUE){
    
    prevalence <- mia::getPrevalence(metagenomics_object, assay.type = "relative_abundance")
    prevalent_topTaxa <- prevalence[prevalence >= p_filterTopPercentage/100] %>% names()
    
    metagenomics_object <- metagenomics_object[prevalent_topTaxa, ]
  }
  
  # Rarefy
  if (p_rarefy == TRUE){
    
    set.seed(0)
    metagenomics_object_rarefied = mia::subsampleCounts(metagenomics_object,
                                                        min_size = 100000, 
                                                        assay_name = "relative_abundance")
    
    # Retrieve the count data
    p_countdata <- assays(metagenomics_object_rarefied)[[2]] %>% t() %>% as.data.frame()
  } else{
    # Retrieve the count data
    p_countdata <- assays(metagenomics_object)[[1]] %>% t() %>% as.data.frame()
  }
  

  # Rename the columns of the count table. Raw features incude full taxonomy - bad feature name
  colnames_new <- data.frame(taxa = colnames(p_countdata)) %>% 
    tidyr::separate(taxa, into = c("k", "p", "c", "o", "f", "g", "s"), sep = "\\|") %>% 
    dplyr::pull(s)
  
  colnames(p_countdata) <- colnames_new
  
  
  # Define pseudocount
  p_pseudocount = 0.5
  p_pseudocount_rel = min((p_countdata/rowSums(p_countdata))[p_countdata/rowSums(p_countdata) != 0])/2
  count_data_zeroimputed = p_countdata
  count_data_zeroimputed[count_data_zeroimputed == 0] = p_pseudocount
  
  
  # Make SummarizedExperiment object
  se_object = SummarizedExperiment(assays = SimpleList(counts = t(as.matrix(p_countdata))))
  
  
  # Carry out transformations within SE object - use the options offered by the "mia" package 
  se_object <- mia::transformAssay(se_object, method = "pa")                                   # Presence-absence transformation
  se_object <- mia::transformAssay(se_object, method = "relabundance")                         # Relative abundance transformation
  se_object <- mia::transformAssay(se_object, method = "clr", 
                                   assay.type = "relabundance", pseudocount = p_pseudocount_rel)    # CLR transformation
  se_object <- mia::transformAssay(se_object, method = "rclr", 
                                   assay.type = "relabundance")    # Robust CLR transformation
  se_object <- mia::transformAssay(se_object, assay.type = "relabundance", 
                                   method = "log10", pseudocount = p_pseudocount_rel)    # Log10 transformation
  
  
  # Extract abundance data
  count_data_PA <- assay(se_object, "pa") %>% t() %>% as.data.frame()               # Presence-absence transformation
  count_data_TSS <- assay(se_object, "relabundance") %>% t() %>% as.data.frame()    # Relative abundance transformation
  count_data_CLR <- assay(se_object, "clr") %>% t() %>% as.data.frame()             # CLR transformation
  count_data_rCLR <- assay(se_object, "rclr") %>% t() %>% as.data.frame()             # Robust CLR transformation
  count_data_logTSS <- assay(se_object, "log10") %>% t() %>% as.data.frame()             # Robust CLR transformation
  
  
  # ALR transformation
  set.seed(0)
  ALR_ref_elements = sample(1:ncol(count_data_zeroimputed), 5)
  for (i in 1:5){
    ALR_run <- compositions::alr(count_data_zeroimputed, ivar = ALR_ref_elements[i]) %>% as.data.frame()
    assign(paste("count_data_ALR", i, sep = ""), ALR_run)
  }
  
  
  # Arcsine transformation
  count_data_aSIN = asin(sqrt(count_data_TSS))
  
  
  # ILR transformation
  count_data_ILR = compositions::ilr(count_data_TSS) %>% 
    as.matrix() %>% 
    as.data.frame()
  
  
  # Merge with metadata and save the files
  for (dataset in c("count_data_PA", "count_data_ILR", "count_data_TSS", "count_data_CLR", "count_data_rCLR", "count_data_logTSS", 
                    "count_data_aSIN",  paste("count_data_ALR", 1:5, sep = ""))){
    
    # Define output name
    transformation_name = substring(dataset, 12)
    if (p_filterTop == TRUE & p_rarefy == TRUE){
      output_name = paste("Datasets_WS_FR/",p_studyName, "_", transformation_name, "_filteredRarefied", p_filterTopPercentage, ".csv", sep = "")
    } else if (p_filterTop == TRUE & p_rarefy == FALSE){
      output_name = paste("Datasets_WS_F/",p_studyName, "_", transformation_name, "_filtered", p_filterTopPercentage, ".csv", sep = "")
    } else if (p_filterTop == FALSE & p_rarefy == TRUE){
      output_name = paste("Datasets_WS_R/",p_studyName, "_", transformation_name, "_rarefied", ".csv", sep = "")
    } else if (p_filterTop == FALSE & p_rarefy == FALSE){
      output_name = paste("Datasets_WS/",p_studyName, "_", transformation_name, ".csv", sep = "")
    }
    
    # Merge abundance data with metadata
    dataset_id = get(dataset) %>% 
      tibble::rownames_to_column(var = "sample_id") 
    
    final_data <- p_metadata %>% 
      dplyr::select(sample_id, target_var) %>% 
      dplyr::left_join(dataset_id, by = "sample_id")
    
    # Save the output
    write.csv(x = final_data, 
              file = file.path("Datasets/", output_name),
              row.names = FALSE)
  }
}



# Helper
apply_transformations <- function(dataset, name){
  
  transform_countData(p_metadata = dataset, p_studyName = name)
  transform_countData(p_metadata = dataset, p_studyName = name, p_filterTop = TRUE, p_filterTopPercentage = 90)
  transform_countData(p_metadata = dataset, p_studyName = name, p_filterTop = TRUE, p_filterTopPercentage = 75)
  transform_countData(p_metadata = dataset, p_studyName = name, p_filterTop = TRUE, p_filterTopPercentage = 50)
  transform_countData(p_metadata = dataset, p_studyName = name, p_filterTop = TRUE, p_filterTopPercentage = 25)
  transform_countData(p_metadata = dataset, p_studyName = name, p_filterTop = TRUE, p_filterTopPercentage = 10)
  print("Unrarefied DONE")
  
  transform_countData(p_metadata = dataset, p_studyName = name, p_rarefy = TRUE)
  transform_countData(p_metadata = dataset, p_studyName = name, p_rarefy = TRUE, p_filterTop = TRUE, p_filterTopPercentage = 90)
  transform_countData(p_metadata = dataset, p_studyName = name, p_rarefy = TRUE, p_filterTop = TRUE, p_filterTopPercentage = 75)
  transform_countData(p_metadata = dataset, p_studyName = name, p_rarefy = TRUE, p_filterTop = TRUE, p_filterTopPercentage = 50)
  transform_countData(p_metadata = dataset, p_studyName = name, p_rarefy = TRUE, p_filterTop = TRUE, p_filterTopPercentage = 25)
  transform_countData(p_metadata = dataset, p_studyName = name, p_rarefy = TRUE, p_filterTop = TRUE, p_filterTopPercentage = 10)
  print("Rarefied DONE")
}







# Apply the function to many datasets ----
#------------------------------------------------#
#                                                #
#               APPLY THE FUNCTION               # 
#                                                #
#------------------------------------------------#

# ---------- QinJ_2012 ----------

# Get the phenotype data used for the analysis
metadata1 <- sampleMetadata %>% 
  dplyr::filter(study_name == "QinJ_2012") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "T2D"))) - 1)

table(metadata1$target_var)
table(metadata1$age_category)

# Apply the function to create transformed datasets
apply_transformations(metadata1, "QinJ_2012")





# ---------- YuJ_2015 ----------

# Get the phenotype data used for the analysis
metadata2 <- sampleMetadata %>% 
  dplyr::filter(study_name == "YuJ_2015") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "CRC"))) - 1)

table(metadata2$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata2, "YuJ_2015")




# ---------- YachidaS_2019 ----------

# Get the phenotype data used for the analysis
metadata3 <- sampleMetadata %>% 
  dplyr::filter(study_name == "YachidaS_2019") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "CRC"))) - 1)

table(metadata3$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata3, "YachidaS_2019")





# ---------- VogtmannE_2016 ----------

# Get the phenotype data used for the analysis
metadata4 <- sampleMetadata %>% 
  dplyr::filter(study_name == "VogtmannE_2016") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "CRC"))) - 1)

table(metadata4$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata4, "VogtmannE_2016")




# ---------- WirbelJ_2018 ----------

# Get the phenotype data used for the analysis
metadata5 <- sampleMetadata %>% 
  dplyr::filter(study_name == "WirbelJ_2018") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "CRC"))) - 1)

table(metadata5$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata5, "WirbelJ_2018")






# ---------- ZhuF_2020 ----------

# Get the phenotype data used for the analysis
metadata6 <- sampleMetadata %>% 
  dplyr::filter(study_name == "ZhuF_2020") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "schizophrenia"))) - 1)

table(metadata6$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata6, "ZhuF_2020")






# ---------- ZellerG_2014 ----------

# Get the phenotype data used for the analysis
metadata7 <- sampleMetadata %>% 
  dplyr::filter(study_name == "ZellerG_2014") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::filter(study_condition %in% c("control", "CRC")) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "CRC"))) - 1)

table(metadata7$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata7, "ZellerG_2014")





# ---------- RubelMA_2020 ----------

# Get the phenotype data used for the analysis
metadata8 <- sampleMetadata %>% 
  dplyr::filter(study_name == "RubelMA_2020") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "STH"))) - 1)

table(metadata8$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata8, "RubelMA_2020")




# ---------- QinN_2014 ----------

# Get the phenotype data used for the analysis
metadata9 <- sampleMetadata %>% 
  dplyr::filter(study_name == "QinN_2014") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "cirrhosis"))) - 1)

table(metadata9$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata9, "QinN_2014")






# ---------- NielsenHB_2014 ----------

# Get the phenotype data used for the analysis
metadata10 <- sampleMetadata %>% 
  dplyr::filter(study_name == "NielsenHB_2014") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(days_from_first_collection == 0) %>% 
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "IBD"))) - 1)

table(metadata10$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata10, "NielsenHB_2014")






# ---------- NagySzakalD_2017 ----------

# Get the phenotype data used for the analysis
metadata11 <- sampleMetadata %>% 
  dplyr::filter(study_name == "NagySzakalD_2017") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "ME/CFS"))) - 1)

table(metadata11$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata11, "NagySzakalD_2017")





# ---------- LiJ_2017 ----------

# Get the phenotype data used for the analysis
metadata12 <- sampleMetadata %>% 
  dplyr::filter(study_name == "LiJ_2017") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::filter(study_condition %in% c("pre-hypertension", "hypertension")) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("pre-hypertension", "hypertension"))) - 1)

table(metadata12$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata12, "LiJ_2017")









# ---------- JieZ_2017 ----------

# Get the phenotype data used for the analysis
metadata14 <- sampleMetadata %>% 
  dplyr::filter(study_name == "JieZ_2017") %>% 
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>% 
  dplyr::mutate(target_var = as.numeric(factor(study_condition, levels = c("control", "ACVD"))) - 1)

table(metadata14$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata14, "JieZ_2017")







#--------------------------------------------------------#
#                ADDITIONAL DATASETS FOR BMI             #
#--------------------------------------------------------#

# ---------- LifeLinesDeep_2016----------

# Get the phenotype data used for the analysis
metadata15 <- sampleMetadata %>%
  dplyr::filter(study_name == "LifeLinesDeep_2016") %>%
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>%
  dplyr::mutate(target_var = as.numeric(as.factor(BMI >= 30))-1)

table(metadata15$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata15, "LifeLinesDeep_2016")




# ---------- AsnicarF_2021----------

# Get the phenotype data used for the analysis
metadata16 <- sampleMetadata %>%
  dplyr::filter(study_name == "AsnicarF_2021") %>%
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>%
  dplyr::mutate(target_var = as.numeric(as.factor(BMI >= 30))-1)

table(metadata16$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata16, "AsnicarF_2021")




# ---------- LiJ_2014----------

# Get the phenotype data used for the analysis
metadata17 <- sampleMetadata %>%
  dplyr::filter(study_name == "LiJ_2014") %>%
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>%
  dplyr::mutate(target_var = as.numeric(as.factor(BMI >= 30))-1)

table(metadata17$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata17, "LiJ_2014")




# ---------- LeChatelierE_2013 ----------

# Get the phenotype data used for the analysis
metadata18 <- sampleMetadata %>%
  dplyr::filter(study_name == "LeChatelierE_2013") %>%
  dplyr::filter(body_site == "stool") %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>%
  dplyr::mutate(target_var = as.numeric(as.factor(BMI >= 30))-1)

table(metadata18$target_var)

# Apply the function to create transformed datasets
apply_transformations(metadata18, "LeChatelierE_2013")



