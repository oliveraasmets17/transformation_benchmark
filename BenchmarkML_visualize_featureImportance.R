

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------# 

# Read the packages
library("tidyverse")
library("ggthemes")
library("ggpubr")



# Feature metadata
EstMB_metadata <- readRDS("/Utilities/Data/MetaPhlan_transformedData_TSS.rds") %>% 
  tibble::rownames_to_column(var = "sample") %>% 
  tidyr::gather(taxa, value, -sample) %>% 
  dplyr::mutate(taxa = ifelse(substr(taxa, 1, 3) != "s__", paste("s__", taxa, sep = ""), taxa)) %>% 
  dplyr::group_by(taxa) %>% 
  dplyr::summarise(prevalence = sum(value > 0)/2509*100, 
                   relabundance = mean(value)*100) %>% 
  dplyr::ungroup()

EstMB_metadata_df <- dplyr::bind_rows(EstMB_metadata, EstMB_metadata, EstMB_metadata, EstMB_metadata) %>% 
  dplyr::mutate(study = c(rep("EstMB_Depression_2022", nrow(EstMB_metadata)),
                          rep("EstMB_AB90_2022", nrow(EstMB_metadata)),
                          rep("EstMB_BMI30_2022", nrow(EstMB_metadata)),
                          rep("EstMB_T2D_2022", nrow(EstMB_metadata))))

feature_metadata_df <- readRDS("/RData/Dataset_characteristics.rds") %>% 
  dplyr::mutate(taxa = ifelse(substr(taxa, 1, 3) != "s__", paste("s__", taxa, sep = ""), taxa))  %>% 
  dplyr::bind_rows(EstMB_metadata_df)







# Read raw data ----
#------------------------------------------------#
#                                                #
#                   READ THE DATA                # 
#                                                #
#------------------------------------------------# 

# Read data 
#---------------------------------------#
# New run
WS_fimp_df <- readRDS("WS_fimp_df.rds") %>% dplyr::mutate(rarefied = FALSE)
WS_R_fimp_df <- readRDS("WS_R_fimp_df.rds") %>% dplyr::mutate(rarefied = TRUE)
EstMB_fimp_df <- readRDS("EstMB_fimp_df.rds") %>% dplyr::mutate(rarefied = FALSE)
EstMB_R_fimp_df <- readRDS("EstMB_R_fimp_df.rds") %>% dplyr::mutate(rarefied = TRUE)


# Combine results 
FI_results_WS <- dplyr::bind_rows(WS_fimp_df, WS_R_fimp_df, EstMB_fimp_df, EstMB_R_fimp_df)



# First data transformations
#---------------------------------------#
FI_results_final <- FI_results_WS %>% 
  dplyr::filter(transformation != "ILR") %>% 
  dplyr::mutate(transformation2 = ifelse(substr(transformation, 1, 3) == "ALR", "ALR", as.character(transformation)),
                transformation2 = factor(transformation2, levels = c("PA", "TSS", "logTSS", "aSIN", "rCLR", "CLR", "ALR"))) %>%
  dplyr::left_join(feature_metadata_df, by = c("study", "f_name" = "taxa")) %>% 
  dplyr::mutate(target = case_when(study %in% c("LifeLinesDeep_2016", "AsnicarF_2021", "LiJ_2014", "LeChatelierE_2013", "EstBMI30metaphlan_2022") ~ "BMI", 
                                   study %in% c("YuJ_2015", "YachidaS_2019", "VogtmannE_2016", "WirbelJ_2018", 
                                                "ZellerG_2014", "ThomasAM_2019c", "ThomasAM_2018b", "ThomasAM_2018a", "HanniganGD_2017", 
                                                "GuptaA_2019", "FengQ_2015") ~ "CRC", 
                                   study %in% c("QinJ_2012", "EstT2Dmetaphlan_2022") ~ "T2D", 
                                   study == "ZhuF_2020" ~ "schizophrenia", 
                                   study == "RubelMA_2020" ~ "STH", 
                                   study == "NielsenHB_2014" ~ "IBD", 
                                   study == "NagySzakalD_2017" ~ "fatigue", 
                                   study == "LiJ_2017" ~ "hypertension", 
                                   study == "KeohaneDM_2020" ~ "smoking", 
                                   study == "JieZ_2017" ~ "ACD", 
                                   study == "QinN_2014" ~ "cirrhosis",
                                   study == "EstDepressionmetaphlan_2022" ~ "depression", 
                                   study == "EstAB90metaphlan_2022" ~ "antibiotics",
                                   TRUE ~ "")) %>% 
  dplyr::mutate(study = case_when(study == "EstAB90metaphlan_2022" ~ "EstMB_AB90_2022",
                                  study == "EstBMI30metaphlan_2022" ~ "EstMB_BMI30_2022",
                                  study == "EstDepressionmetaphlan_2022" ~ "EstMB_Depression_2022",
                                  study == "EstT2Dmetaphlan_2022" ~ "EstMB_T2D_2022",
                                  TRUE ~ study))








# Ordination figures ----
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#   

# Run PCoA on the same ordination - nonrarefied
#------------------------------------------#
FI_results_pcoa <- FI_results_final %>%  
  dplyr::filter(rarefied == FALSE) %>% 
  dplyr::group_by(study, target, algorithm, f_name, transformation2) %>% 
  dplyr::summarise(mean_shap = mean(shap_value)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(f_name, mean_shap, fill = 0) 

# Run PCA
dist_both <- vegan::vegdist(FI_results_pcoa[ ,5:ncol(FI_results_pcoa)], method = "euclidean")
pcoa_both <- stats::cmdscale(dist_both)

# Visualize
pca_res_data <- data.frame(pcoa_both) %>% 
  dplyr::bind_cols(FI_results_pcoa[ ,1:4])

p_pcoa_sameOrd <- ggplot(pca_res_data %>% dplyr::filter(transformation2 %in% c("CLR", "logTSS", "PA", "aSIN", "ALR", "rCLR", "TSS")), 
                         aes(x = X1, y = X2, color = target, shape = transformation2)) + 
  geom_point(size = 2, stroke = 1.5) + 
  scale_shape_manual(name = "Transformation", values = c(0, 1, 2, 3, 4, 5, 10, 12)) + 
  scale_color_manual(name = "Target", values = c("#a6cee3", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#1f78b4", "#b2df8a", "#33a02c", "#b15928")) + 
  facet_wrap(vars(factor(algorithm, levels = c("ENET", "RF", "XGB")))) + 
  
  theme_bw() + 
  xlab("PC 1") + 
  ylab("PC 2") + 
  theme(axis.title = element_text(size = 20), 
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14))
p_pcoa_sameOrd

ggsave(plot = p_pcoa_sameOrd, "/Bioinf/BenchmarkML/Figures/p_pcoa_sameOrdination_eucliden.png", width = 18, height = 6)
ggsave(plot = p_pcoa_sameOrd, "/Bioinf/BenchmarkML/Figures/p_pcoa_sameOrdination_eucliden.pdf", width = 18, height = 6)




# Run PCoA on the same ordination - rarefied
#------------------------------------------#
FI_results_pcoa_rarefied <- FI_results_final %>%  
  dplyr::filter(rarefied == TRUE) %>% 
  dplyr::group_by(study, target, algorithm, f_name, transformation2) %>% 
  dplyr::summarise(mean_shap = mean(shap_value)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(f_name, mean_shap, fill = 0) 

# Run PCA
dist_both_rarefied <- vegan::vegdist(FI_results_pcoa_rarefied[ ,5:ncol(FI_results_pcoa_rarefied)], method = "euclidean")
pcoa_both_rarefied <- stats::cmdscale(dist_both_rarefied)

# Visualize
pca_res_data_rarefied <- data.frame(pcoa_both_rarefied) %>% 
  dplyr::bind_cols(FI_results_pcoa_rarefied[ ,1:4])

p_pcoa_sameOrd_rarefied <- ggplot(pca_res_data_rarefied %>% dplyr::filter(transformation2 %in% c("CLR", "logTSS", "PA", "aSIN", "ALR", "rCLR", "TSS")), 
                         aes(x = X1, y = X2, color = target, shape = transformation2)) + 
  geom_point(size = 2, stroke = 1.5) + 
  scale_shape_manual(name = "Transformation", values = c(0, 1, 2, 3, 4, 5, 10, 12)) + 
  scale_color_manual(name = "Target", values = c("#a6cee3", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#1f78b4", "#b2df8a", "#33a02c", "#b15928")) + 
  facet_wrap(vars(factor(algorithm, levels = c("ENET", "RF", "XGB")))) + 
  
  theme_bw() + 
  xlab("PC 1") + 
  ylab("PC 2") + 
  theme(axis.title = element_text(size = 20), 
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14))
p_pcoa_sameOrd_rarefied

ggsave(plot = p_pcoa_sameOrd_rarefied, "/Bioinf/BenchmarkML/Figures/p_pcoa_sameOrdination_eucliden_rarefied.png", width = 20, height = 8)
ggsave(plot = p_pcoa_sameOrd_rarefied, "/Bioinf/BenchmarkML/Figures/p_pcoa_sameOrdination_eucliden_rarefied.pdf", width = 20, height = 8)












# Similarity figures ----
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#   

# Calculate distances between the combinations
#------------------------------------------#
run_distance_combinations = FALSE
if (run_distance_combinations == TRUE){
  
  FI_results_pcoa <- FI_results_final %>%  
    dplyr::filter(rarefied == FALSE) %>% 
    dplyr::group_by(study, target, algorithm, f_name, transformation2) %>% 
    dplyr::summarise(mean_shap = mean(shap_value)) %>% 
    dplyr::ungroup() %>% 
    tidyr::spread(f_name, mean_shap, fill = 0) %>% 
    dplyr::mutate(setting = paste(study, "_p_", target, "_p_", algorithm, "_p_", transformation2, sep = ""))
  
  dist_mat <- vegan::vegdist(FI_results_pcoa[ ,5:(ncol(FI_results_pcoa)-1)], method = "euclidean") %>% 
    as.matrix()
  rownames(dist_mat) <- FI_results_pcoa$setting
  colnames(dist_mat) <- FI_results_pcoa$setting
  
  output_df <- data.frame()
  
  settings = unique(FI_results_pcoa$setting)
  for (i in 1:length(settings)){
    for(j in 1:length(settings)){
      
      setting1 <- settings[i] %>% str_split(pattern = "_p_")
      setting2 <- settings[j] %>% str_split(pattern = "_p_")
      
      study1 <- setting1[[1]][1]
      study2 <- setting2[[1]][1]
      
      study = ifelse(study1 == study2, "same", "different")
      
      algo1 <- setting1[[1]][3]
      algo2 <- setting2[[1]][3]
      
      algo = ifelse(algo1 == algo2, "same", "different")
      
      trans1 <- setting1[[1]][4]
      trans2 <- setting2[[1]][4]
      
      trans = ifelse(trans1 == trans2, "same", "different")
      
      dist_submat <- dist_mat[settings[i], settings[j]]
      
      run_df <- data.frame(study1, study2, study, algo1, algo2, algo, trans1, trans2, trans, dist = dist_submat)
      output_df <- dplyr::bind_rows(output_df, run_df)
    }
  }
  saveRDS(output_df, "RData/dist_euclidean_df.rds")
} else{
  output_df <- readRDS("RData/dist_euclidean_df.rds")
}

help_df1 <- output_df %>% 
  dplyr::filter(dist != 0) %>% 
  dplyr::filter(algo == "same") %>% 
  dplyr::mutate(group = case_when(study == "same" & trans == "different" ~ "Same study\ndifferent transformation",
                                  study == "different" & trans == "different" ~ "Different study\ndifferent transformation",
                                  study == "different" & trans == "same" ~ "Different study\nsame transformation"))

my_comparisons = list(c("Same study\ndifferent transformation", "Different study\ndifferent transformation"),
                      c("Same study\ndifferent transformation", "Different study\nsame transformation"),
                      c("Different study\nsame transformation", "Different study\ndifferent transformation"))

p_distcomparison <- ggplot(help_df1, aes(x = factor(group, levels = c("Same study\ndifferent transformation", "Different study\ndifferent transformation", "Different study\nsame transformation")), y = dist)) + 
  geom_boxplot() + 
  facet_wrap(vars(algo1)) + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  xlab("") + 
  ylab("Euclidean distance") + 
  theme(strip.text = element_text(size = 20), 
        axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
        strip.background = element_rect(fill = "white"))
p_distcomparison

ggsave(plot = p_distcomparison, "/Bioinf/BenchmarkML/Figures/distcomparison.png", width = 18, height = 8)
ggsave(plot = p_distcomparison, "/Bioinf/BenchmarkML/Figures/distcomparison.pdf", width = 18, height = 8)



# Calculate CohensD
library(effsize)
x = help_df1 %>% filter(algo1 == "ENET" & group == "Same study\ndifferent transformation") %>% dplyr::pull(dist)
y = help_df1 %>% filter(algo1 == "ENET" & group == "Different study\ndifferent transformation") %>% dplyr::pull(dist)

cohen.d(x, y)

x = help_df1 %>% filter(algo1 == "RF" & group == "Same study\ndifferent transformation") %>% dplyr::pull(dist)
y = help_df1 %>% filter(algo1 == "RF" & group == "Different study\ndifferent transformation") %>% dplyr::pull(dist)

cohen.d(x, y)

x = help_df1 %>% filter(algo1 == "XGB" & group == "Same study\ndifferent transformation") %>% dplyr::pull(dist)
y = help_df1 %>% filter(algo1 == "XGB" & group == "Different study\ndifferent transformation") %>% dplyr::pull(dist)

cohen.d(x, y)






# Correlation between SHAP values and feature metadata
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#  
feature_metadata_imp <- FI_results_final %>% 
  dplyr::group_by(rarefied, f_name, study, transformation2, algorithm) %>% 
  dplyr::summarize(mean_shap = mean(shap_value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(feature_metadata_df, by = c("study", "f_name" = "taxa")) %>% 
  dplyr::group_by(rarefied, algorithm, transformation2) %>% 
  dplyr::summarise(cor_prevalence = cor(mean_shap, prevalence, use = "pairwise.complete.obs"), 
                   cor_abundance = cor(mean_shap, relabundance, use = "pairwise.complete.obs"))

p_featureMetadata <- ggplot(feature_metadata_imp, aes(x = cor_prevalence, y = cor_abundance, color = transformation2, shape = algorithm)) + 
  geom_hline(yintercept = 0, linetype = 3, size = 1) + 
  geom_vline(xintercept = 0, linetype = 3, size = 1) + 
  geom_point(size = 7, alpha = 0.8) +
  scale_shape_manual(values = c(17, 16, 18), name = "Algorithm") + 
  scale_color_calc(name = "Transformation") +
  xlab("Correlation with taxa prevalence") + 
  ylab("Correlation with taxa abundance") + 
  facet_wrap(vars(factor(rarefied, levels = c(FALSE, TRUE), labels = c("Non-rarefied", "Rarefied")))) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 12))

p_featureMetadata
ggsave(plot = p_featureMetadata, "/Bioinf/BenchmarkML/Figures/p_featureMetadata.png", width = 16, height = 5)
ggsave(plot = p_featureMetadata, "/Bioinf/BenchmarkML/Figures/p_featureMetadata.pdf", width = 16, height = 5)




# correlation between SHAP values and feature metadata - BY STUDY
#------------------------------------------#
feature_metadata_imp_study <- FI_results_final %>% 
  dplyr::group_by(f_name, study, transformation2, algorithm) %>% 
  dplyr::summarize(mean_shap = mean(shap_value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(feature_metadata_df, by = c("study", "f_name" = "taxa")) %>% 
  dplyr::group_by(algorithm, study, transformation2) %>% 
  dplyr::summarise(cor_prevalence = cor(mean_shap, prevalence), 
                   cor_abundance = cor(mean_shap, relabundance))

p_featureMetadata_byStudy <- ggplot(feature_metadata_imp_study, aes(x = cor_prevalence, y = cor_abundance, color = transformation2, shape = algorithm)) + 
  geom_hline(yintercept = 0, linetype = 3, size = 1) + 
  geom_vline(xintercept = 0, linetype = 3, size = 1) + 
  geom_point(size = 3, alpha = 0.7) +
  scale_shape_manual(values = c(16, 17, 18)) + 
  scale_color_calc(name = "Transformation") +
  facet_wrap(vars(target)) + 
  xlab("Correlation with taxa prevalence") + 
  ylab("Correlation with taxa abundance") + 
  theme_bw() + 
  facet_wrap(vars(study), ncol = 3) + 
  theme(axis.title = element_text(size = 16),
        strip.background = element_rect(fill = "white"), 
        axis.text = element_text(size = 12))

p_featureMetadata_byStudy
ggsave(plot = p_featureMetadata_byStudy, "/Bioinf/BenchmarkML/Figures/p_featureMetadata_byStudy.png", width = 12, height = 12)
ggsave(plot = p_featureMetadata_byStudy, "/Bioinf/BenchmarkML/Figures/p_featureMetadata_byStudy.pdf", width = 12, height = 12)






# Percentage of features selected per target x transformation
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#  

n_imp <- FI_results_final %>% 
  dplyr::group_by(study, fold, algorithm, transformation2) %>% 
  dplyr::summarise(n_chosen = sum(abs(shap_value) > 0),
                   n = n())  %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(proc = n_chosen/n,
                model = paste(algorithm, transformation2, sep = "_")) %>% 
  dplyr::group_by(study, algorithm, transformation2, model) %>% 
  dplyr::summarise(mean_proc = mean(proc))  %>% 
  dplyr::ungroup()


# Define order of the rows/cols 
n_imp_wide <- n_imp %>% 
  dplyr::filter(algorithm == "ENET") %>% 
  dplyr::select(study, transformation2, mean_proc) %>% 
  tidyr::spread(transformation2, mean_proc) %>% 
  tibble::column_to_rownames(var = "study") 

hclust_nimp_study <- hclust(dist(n_imp_wide %>% filter(complete.cases(.))))
hclust_nimp_trans <- hclust(dist(t(n_imp_wide)))

study_order <- hclust_nimp_study$labels[hclust_nimp_study$order]
trans_order <- hclust_nimp_trans$labels[hclust_nimp_trans$order]


# Plot
p_procFeaturesChosen <- ggplot(n_imp %>% dplyr::filter(study %in% study_order), 
                               aes(x = factor(transformation2, levels = trans_order), 
                                   y = factor(study, levels = study_order), 
                                   fill = 100*mean_proc, label = round(100*mean_proc, 1))) + 
  geom_tile(color = "white") + 
  geom_text(size = 6) + 
  scale_fill_gradient(name = "Percentage of features with abs(SHAP) > 0", low = "white", high = "deepskyblue4") +  
  theme_classic()  + 
  xlab("") + 
  facet_wrap(vars(algorithm)) + 
  ylab("") + 
  theme(axis.text = element_text(size = 16),
        legend.position = "bottom", 
        legend.key.width = unit(1, 'cm'),
        legend.text = element_text(size = 16), 
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20))
p_procFeaturesChosen
ggsave(plot = p_procFeaturesChosen, "/Bioinf/BenchmarkML/Figures/p_procFeaturesChosen.png", width = 20, height = 10)
ggsave(plot = p_procFeaturesChosen, "/Bioinf/BenchmarkML/Figures/p_procFeaturesChosen.pdf", width = 20, height = 10)






# Overlap of top 25 features
#---------------------------------------#
top25_features <- FI_results_final %>%
  dplyr::group_by(rarefied, study, transformation, f_name, algorithm) %>% 
  dplyr::summarise(mean_SHAP = mean(abs(shap_value))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(help = 1) %>% 
  dplyr::arrange(study, transformation, algorithm, desc(mean_SHAP)) %>% 
  dplyr::group_by(rarefied, study, transformation, algorithm) %>% 
  dplyr::mutate(counter = cumsum(help)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(counter <= 25) %>% 
  dplyr::mutate(transformation = as.character(transformation))

# Calculate the overlap
transformation_list <- c("PA", "aSIN", "TSS", "logTSS", "CLR", "rCLR", "ALR1", "ALR2", "ALR3", "ALR4", "ALR5")

overlap_df <- data.frame()
interaction_df <- data.frame()

for (i in transformation_list){
  for (j in transformation_list){
    
    if (i != j){
      
      within_algorithm_overlap = top25_features %>% 
        dplyr::filter(transformation %in% c(i, j)) %>% 
        dplyr::mutate(transformation = ifelse(transformation == i, "one", "two")) %>% 
        dplyr::select(study, transformation, algorithm, rarefied, f_name, help) %>% 
        tidyr::spread(transformation, help, fill = 0) %>% 
        dplyr::group_by(algorithm, rarefied, study) %>% 
        dplyr::summarise(n_overlap = sum(one == 1 & two == 1)) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(percentage = n_overlap/25*100) %>% 
        dplyr::mutate(tr1 = i, 
                      tr2 = j)
      
      between_algorithm_overlap = top25_features %>% 
        dplyr::filter(transformation %in% c(i, j)) %>% 
        dplyr::mutate(transformation = ifelse(transformation == i, "one", "two"),
                      setting = paste(algorithm, transformation, sep = "_")) %>% 
        dplyr::select(rarefied, study, setting, f_name, help) %>% 
        tidyr::spread(setting, help, fill = 0) %>% 
        dplyr::group_by(rarefied, study) %>% 
        dplyr::summarise(n_overlap_ENET_RF = sum(ENET_one == 1 & RF_two == 1),
                         n_overlap_ENET_XGB = sum(ENET_one == 1 & XGB_two == 1),
                         n_overlap_XGB_RF = sum(XGB_one == 1 & RF_two == 1)) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(percentage_ENET_RF = n_overlap_ENET_RF/25*100,
                      percentage_ENET_XGB = n_overlap_ENET_XGB/25*100,
                      percentage_XGB_RF = n_overlap_XGB_RF/25*100) %>% 
        dplyr::mutate(tr1 = i, 
                      tr2 = j, 
                      algorithm = "Interaction")
      
      overlap_df = dplyr::bind_rows(overlap_df, within_algorithm_overlap)
      interaction_df = dplyr::bind_rows(interaction_df, between_algorithm_overlap)
      
    } else{
      between_algorithm_overlap = top25_features %>%
        dplyr::select(-one_of("counter", "mean_SHAP")) %>%
        dplyr::filter(transformation == i) %>%
        dplyr::select(rarefied, study, algorithm, f_name, help) %>%
        tidyr::spread(algorithm, help, fill = 0) %>%
        dplyr::group_by(rarefied, study) %>%
        dplyr::summarise(n_overlap_ENET_RF = sum(ENET == 1 & RF == 1),
                         n_overlap_ENET_XGB = sum(ENET == 1 & XGB == 1),
                         n_overlap_XGB_RF = sum(XGB == 1 & RF == 1)) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(percentage_ENET_RF = n_overlap_ENET_RF/25*100,
                      percentage_ENET_XGB = n_overlap_ENET_XGB/25*100,
                      percentage_XGB_RF = n_overlap_XGB_RF/25*100) %>% 
        dplyr::ungroup() %>%
        dplyr::mutate(tr1 = i,
                      tr2 = j,
                      algorithm = "Interaction")
      interaction_df = dplyr::bind_rows(interaction_df, between_algorithm_overlap)
    }
  }
}

overlap_df_aggr = overlap_df %>% 
  dplyr::mutate(tr1_mod = ifelse(substr(tr1, 1, 3) == "ALR", "ALR", tr1),
                tr2_mod = ifelse(substr(tr2, 1, 3) == "ALR", "ALR", tr2)) %>% 
  dplyr::group_by(rarefied, algorithm, tr1_mod, tr2_mod) %>% 
  dplyr::summarise(percentage = mean(percentage)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(print = as.character(round(percentage, 3))) %>% 
  dplyr::bind_rows(data.frame(tr1_mod = c("PA", "aSIN", "TSS", "logTSS", "CLR", "rCLR"), 
                              tr2_mod = c("PA", "aSIN", "TSS", "logTSS", "CLR", "rCLR"), 
                              print = "-")) %>% 
  dplyr::filter(!(tr1_mod == "ALR" & tr2_mod == "ALR"))


# Define the ordering - ENET
overlap_df_aggr_wide <- overlap_df_aggr %>% 
  dplyr::filter(rarefied == FALSE) %>% 
  dplyr::filter(algorithm == "ENET") %>% 
  dplyr::select(-one_of("print", "algorithm")) %>% 
  tidyr::spread(tr2_mod, percentage) %>% 
  tibble::column_to_rownames(var = "tr1_mod")

order_hclust <- hclust(dist(overlap_df_aggr_wide))
order_labels_ENET <- order_hclust$labels[order_hclust$order]


# Plot
p_featureOverlapPercentage25 <- ggplot(overlap_df_aggr %>% 
                                         dplyr::filter(rarefied == FALSE) %>% 
                                         dplyr::filter(is.na(algorithm) == F) %>% 
                                         dplyr::filter(algorithm != "Interaction"), 
                                       aes(x = factor(tr1_mod, levels = order_labels_ENET), 
                                           y = factor(tr2_mod, levels = order_labels_ENET), 
                                           fill = signif(percentage, 2))) + 
  geom_tile(color = "white") + 
  geom_text(aes(label = substr(print,1, 4)), size = 6) + 
  scale_fill_gradient(name = "Percentage of overlapping features", low = "white", high = "maroon3") + 
  xlab("") + 
  facet_wrap(vars(factor(algorithm, levels = c("ENET", "RF", "XGB", "Interaction")))) +
  ylab("") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 16), 
        legend.position = "bottom", 
        legend.key.width = unit(1, 'cm'),
        panel.grid = element_blank(),
        legend.text = element_text(size = 16), 
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(size = 20))

p_featureOverlapPercentage25
ggsave(plot = p_featureOverlapPercentage25, "/Bioinf/BenchmarkML/Figures/p_featureOverlapPercentage25.png", width = 20, height = 8)
ggsave(plot = p_featureOverlapPercentage25, "/Bioinf/BenchmarkML/Figures/p_featureOverlapPercentage25.pdf", width = 20, height = 8)




# Interaction plot
interaction_df_aggr = interaction_df %>% 
  dplyr::mutate(tr1_mod = ifelse(substr(tr1, 1, 3) == "ALR", "ALR", tr1),
                tr2_mod = ifelse(substr(tr2, 1, 3) == "ALR", "ALR", tr2)) %>% 
  dplyr::group_by(rarefied, algorithm, tr1_mod, tr2_mod) %>% 
  dplyr::summarise(percentage_ENET_RF = mean(percentage_ENET_RF), 
                   percentage_ENET_XGB = mean(percentage_ENET_XGB), 
                   percentage_XGB_RF = mean(percentage_XGB_RF)) %>% 
  dplyr::ungroup() %>% 
  tidyr::gather(key, percentage, -c("rarefied", "algorithm", "tr1_mod", "tr2_mod")) %>% 
  dplyr::mutate(key_print = case_when(key == "percentage_ENET_RF" ~ "ENET vs RF", 
                                      key == "percentage_ENET_XGB" ~ "ENET vs XGB",
                                      key == "percentage_XGB_RF" ~ "XGB vs RF"),
                print = as.character(round(percentage, 3))) %>% 
  dplyr::bind_rows(data.frame(tr1_mod = c("PA", "aSIN", "TSS", "logTSS", "CLR", "rCLR"), 
                              tr2_mod = c("PA", "aSIN", "TSS", "logTSS", "CLR", "rCLR"), 
                              print = "-"))

f_featureOverlapPercentage25_interaction <- function(algo1, algo2, p_rarefied){
  
  plotData <- interaction_df_aggr %>% 
    dplyr::filter(rarefied == p_rarefied) %>%
    dplyr::filter(is.na(algorithm) == F) %>% 
    dplyr::filter(key_print == paste(algo1, " vs ", algo2, sep = ""))
  
  p <- ggplot(plotData, 
              aes(x = factor(tr1_mod, levels = order_labels_ENET), 
                  y = factor(tr2_mod, levels = order_labels_ENET), 
                  fill = signif(percentage, 2))) + 
  geom_tile(color = "white") + 
  geom_text(aes(label = substr(print,1, 4)), size = 6) + 
  scale_fill_gradient(name = "Percentage of overlapping features", low = "white", high = "maroon3") + 
  xlab(algo1) + 
  facet_wrap(vars(factor(key_print))) +
  ylab(algo2) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 20),
        legend.position = "bottom", 
        legend.key.width = unit(1, 'cm'),
        panel.grid = element_blank(),
        legend.text = element_text(size = 16), 
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill = "white"),
        legend.title = element_text(size = 20))
  return (p)
}

p1 <- f_featureOverlapPercentage25_interaction(algo1 = "ENET", algo2  = "RF", p_rarefied = FALSE)
p2 <- f_featureOverlapPercentage25_interaction(algo1 = "ENET", algo2  = "XGB", p_rarefied = FALSE)
p3 <- f_featureOverlapPercentage25_interaction(algo1 = "XGB", algo2  = "RF", p_rarefied = FALSE)

# Combine figures
p_featureOverlapPercentage25_interaction <- ggpubr::ggarrange(p1, p2, p3, nrow = 1, common.legend = T, legend = "bottom")
p_featureOverlapPercentage25_interaction
ggsave(plot = p_featureOverlapPercentage25_interaction, "/Bioinf/BenchmarkML/Figures/p_featureOverlapPercentage25_interaction.png", width = 20, height = 8)
ggsave(plot = p_featureOverlapPercentage25_interaction, "/Bioinf/BenchmarkML/Figures/p_featureOverlapPercentage25_interaction.pdf", width = 20, height = 8)







# SHAP drop plot
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#  

Average_SHAP <- FI_results_final %>% 
  dplyr::group_by(rarefied, study, algorithm, transformation2, f_name) %>% 
  dplyr::summarise(mean_shap = mean(shap_value)) %>% 
  dplyr::mutate(help = 1) %>% 
  dplyr::arrange(study, algorithm, transformation2, desc(mean_shap)) %>% 
  dplyr::group_by(study, algorithm, transformation2) %>% 
  dplyr::mutate(counter = cumsum(help)) %>% 
  dplyr::filter(str_detect(study, "Est") == F)

p_SHAP_drop <- ggplot(Average_SHAP %>% dplyr::filter(counter <= 100) %>% dplyr::filter(rarefied == FALSE), 
                      aes(x = counter, y = mean_shap, color = transformation2, group = study)) + 
  geom_path(color = "gray70") + 
  geom_point() + 
  scale_color_calc(name = "Transformation") + 
  theme_bw() + 
  ylab("SHAP value") +
  xlab("Feature rank by SHAP values") + 
  facet_grid(algorithm ~ transformation2, scales = "free") + 
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 18),
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 12))
p_SHAP_drop

ggsave(plot = p_SHAP_drop, "/Bioinf/BenchmarkML/Figures/p_SHAP_drop.png", width = 20, height = 10)
ggsave(plot = p_SHAP_drop, "/Bioinf/BenchmarkML/Figures/p_SHAP_drop.pdf", width = 20, height = 10)





# Overrepresentation of abundant taxa among the pairs
#------------------------------------------#
taxa_distribution <- FI_results_final %>%
  dplyr::filter(rarefied == FALSE) %>% 
  dplyr::mutate(transformation2 = ifelse(substr(transformation, 1, 3) == "ALR", "ALR", as.character(transformation))) %>% 
  dplyr::group_by(study, transformation2, f_name, algorithm) %>% 
  dplyr::summarise(mean_SHAP = mean(abs(shap_value))) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(target = case_when(study %in% c("LifeLinesDeep_2016", "AsnicarF_2021", "LiJ_2014", "LeChatelierE_2013", "EstBMI30metaphlan_2022") ~ "BMI", 
                                   study %in% c("YuJ_2015", "YachidaS_2019", "VogtmannE_2016", "WirbelJ_2018", 
                                                "ZellerG_2014", "ThomasAM_2019c", "ThomasAM_2018b", "ThomasAM_2018a", "HanniganGD_2017", 
                                                "GuptaA_2019", "FengQ_2015") ~ "CRC", 
                                   study %in% c("QinJ_2012", "EstT2Dmetaphlan_2022") ~ "T2D", 
                                   study == "ZhuF_2020" ~ "schizophrenia", 
                                   study == "RubelMA_2020" ~ "STH", 
                                   study == "NielsenHB_2014" ~ "IBD", 
                                   study == "NagySzakalD_2017" ~ "fatigue", 
                                   study == "LiJ_2017" ~ "hypertension", 
                                   study == "KeohaneDM_2020" ~ "smoking", 
                                   study == "JieZ_2017" ~ "ACD", 
                                   study == "QinN_2014" ~ "cirrhosis",
                                   study == "EstDepressionmetaphlan_2022" ~ "depression", 
                                   study == "EstAB90metaphlan_2022" ~ "antibiotics",
                                   TRUE ~ ""),
                help = 1) %>% 
  dplyr::left_join(feature_metadata_summary, by = c("f_name" = "taxa")) %>% 
  dplyr::arrange(study, transformation2, algorithm, target, desc(mean_SHAP)) %>% 
  dplyr::group_by(study, transformation2, algorithm, target) %>% 
  dplyr::mutate(counter = cumsum(help)) %>% 
  dplyr::ungroup()


CRC_FI_generalTop <- taxa_distribution %>%
  dplyr::filter(transformation2 %in% c("PA", "TSS")) %>% 
  dplyr::select(study, transformation2, f_name, algorithm, mean_SHAP, counter) %>% 
  tidyr::pivot_wider(names_from = transformation2, values_from = c(mean_SHAP, counter)) %>% 
  dplyr::group_by(f_name, algorithm) %>% 
  dplyr::summarise(mean_delta_SHAP = mean(mean_SHAP_PA - mean_SHAP_TSS),
                   mean_delta_rank = mean(counter_PA - counter_TSS)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(algorithm, desc(abs(mean_delta_SHAP))) %>% 
  dplyr::left_join(feature_metadata_summary, by = c("f_name" = "taxa")) 


# ENET
taxa_order <- CRC_FI_generalTop %>% 
  dplyr::filter(algorithm == "ENET") %>%
  dplyr::mutate(z_delta_SHAP = (mean_delta_SHAP - mean(mean_delta_SHAP))/sd(mean_delta_SHAP)) %>% 
  dplyr::top_n(100, wt = abs(mean_delta_SHAP)) %>% 
  dplyr::arrange(mean_abund) 

abs_max <- taxa_order %>% 
  dplyr::filter(abs(z_delta_SHAP) >= 3)


p0 <- ggplot(taxa_order, aes(x = factor(f_name, levels = taxa_order$f_name), y = mean_abund)) + 
  geom_bar(stat = "identity", fill = "darkred") + 
  theme_classic() +
  xlab("") + 
  ylab("Relative \nabundance (%)") + 
  theme(axis.text.y = element_blank(), 
        axis.ticks = element_blank()) + 
  coord_flip()

p1 <- ggplot(CRC_FI_generalTop %>% dplyr::filter(algorithm == "ENET") %>% dplyr::filter(f_name %in% taxa_order$f_name), aes(x = factor(f_name, levels = taxa_order$f_name), y = mean_delta_SHAP)) + 
  geom_bar(stat = "identity", fill = "gray70") +
  geom_text(data = abs_max, aes(x = factor(f_name, levels = taxa_order$f_name), y = mean_delta_SHAP, label = f_name), size = 3) + 
  geom_hline(yintercept = 0, linetype = 3, size = 1) + 
  theme_bw() + 
  ggtitle("ENET - differences in SHAP values between PA and TSS") + 
  scale_y_continuous(limits = c(-0.01, 0.008)) + 
  annotate("text", x = 80, y = 0.006, label = "More important in PA", color = "blue", size = 5) + 
  annotate("text", x = 20, y = -0.005, label = "More important in TSS", color = "blue", size = 5) + 
  xlab("") + 
  ylab("Average difference in SHAP values (PA-TSS)") + 
  theme(axis.text.y = element_blank(), 
        title = element_text(size = 18), 
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 14), 
        panel.grid.major.y = element_blank()) + 
  coord_flip()
p1  

p2 <- ggarrange(p1, p0, 
                widths = c(10, 2),
                align = "h", 
                ncol = 2)
p2
ggsave(p2, file = "/Bioinf/BenchmarkML/Figures/ENET_featureDifference.png", width = 10, height = 11)
ggsave(p2, file = "/Bioinf/BenchmarkML/Figures/ENET_featureDifference.pdf", width = 10, height = 11)




# RF
taxa_order <- CRC_FI_generalTop %>% 
  dplyr::filter(algorithm == "RF") %>%
  dplyr::mutate(z_delta_SHAP = (mean_delta_SHAP - mean(mean_delta_SHAP))/sd(mean_delta_SHAP)) %>% 
  dplyr::top_n(100, wt = abs(mean_delta_SHAP)) %>% 
  dplyr::arrange(mean_abund) 

abs_max <- taxa_order %>% 
  dplyr::filter(abs(z_delta_SHAP) >= 3)


p0 <- ggplot(taxa_order, aes(x = factor(f_name, levels = taxa_order$f_name), y = mean_abund)) + 
  geom_bar(stat = "identity", fill = "darkred") + 
  theme_classic() +
  xlab("") + 
  ylab("Relative \nabundance (%)") + 
  theme(axis.text.y = element_blank(), 
        axis.ticks = element_blank()) + 
  coord_flip()

p1 <- ggplot(CRC_FI_generalTop %>% dplyr::filter(algorithm == "RF") %>% dplyr::filter(f_name %in% taxa_order$f_name), aes(x = factor(f_name, levels = taxa_order$f_name), y = mean_delta_SHAP)) + 
  geom_bar(stat = "identity", fill = "gray70") +
  geom_text(data = abs_max, aes(x = factor(f_name, levels = taxa_order$f_name), y = mean_delta_SHAP, label = f_name), size = 3) + 
  geom_hline(yintercept = 0, linetype = 3, size = 1) + 
  theme_bw() + 
  ggtitle("RF - differences in SHAP values between PA and TSS") + 
  scale_y_continuous(limits = c(-0.0035, 0.0035)) + 
  annotate("text", x = 80, y = 0.002, label = "More important in PA", color = "blue", size = 5) + 
  annotate("text", x = 20, y = -0.0025, label = "More important in TSS", color = "blue", size = 5) + 
  xlab("") + 
  ylab("Average difference in SHAP values (PA-TSS)") + 
  theme(axis.text.y = element_blank(), 
        title = element_text(size = 18), 
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 14), 
        panel.grid.major.y = element_blank()) + 
  coord_flip()
p1  

p2 <- ggarrange(p1, p0, 
                widths = c(10, 2),
                align = "h", 
                ncol = 2)
p2
ggsave(p2, file = "/Bioinf/BenchmarkML/Figures/RF_featureDifference.png", width = 10, height = 12)
ggsave(p2, file = "/Bioinf/BenchmarkML/Figures/RF_featureDifference.pdf", width = 10, height = 12)







# XGB
taxa_order <- CRC_FI_generalTop %>% 
  dplyr::filter(algorithm == "XGB") %>%
  dplyr::mutate(z_delta_SHAP = (mean_delta_SHAP - mean(mean_delta_SHAP))/sd(mean_delta_SHAP)) %>% 
  dplyr::top_n(100, wt = abs(mean_delta_SHAP)) %>% 
  dplyr::arrange(mean_abund) 

abs_max <- taxa_order %>% 
  dplyr::filter(abs(z_delta_SHAP) >= 3)


p0 <- ggplot(taxa_order, aes(x = factor(f_name, levels = taxa_order$f_name), y = mean_abund)) + 
  geom_bar(stat = "identity", fill = "darkred") + 
  theme_classic() +
  xlab("") + 
  ylab("Relative \nabundance (%)") + 
  theme(axis.text.y = element_blank(), 
        axis.ticks = element_blank()) + 
  coord_flip()

p1 <- ggplot(CRC_FI_generalTop %>% dplyr::filter(algorithm == "XGB") %>% dplyr::filter(f_name %in% taxa_order$f_name), aes(x = factor(f_name, levels = taxa_order$f_name), y = mean_delta_SHAP)) + 
  geom_bar(stat = "identity", fill = "gray70") +
  geom_text(data = abs_max, aes(x = factor(f_name, levels = taxa_order$f_name), y = mean_delta_SHAP, label = f_name), size = 3) + 
  geom_hline(yintercept = 0, linetype = 3, size = 1) + 
  theme_bw() + 
  ggtitle("XGB - differences in SHAP values between PA and TSS") + 
  scale_y_continuous(limits = c(-0.0035, 0.0035)) + 
  annotate("text", x = 80, y = 0.002, label = "More important in PA", color = "blue", size = 5) + 
  annotate("text", x = 20, y = -0.0025, label = "More important in TSS", color = "blue", size = 5) + 
  xlab("") + 
  ylab("Average difference in SHAP values (PA-TSS)") + 
  theme(axis.text.y = element_blank(), 
        title = element_text(size = 18), 
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 14), 
        panel.grid.major.y = element_blank()) + 
  coord_flip()
p1  

p2 <- ggarrange(p1, p0, 
                widths = c(10, 2),
                align = "h", 
                ncol = 2)
p2
ggsave(p2, file = "/Bioinf/BenchmarkML/Figures/XGB_featureDifference.png", width = 10, height = 12)
ggsave(p2, file = "/Bioinf/BenchmarkML/Figures/XGB_featureDifference.pdf", width = 10, height = 12)
