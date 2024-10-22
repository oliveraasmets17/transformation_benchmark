

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------# 

# Read the packages
library("tidyverse")
library("ggthemes")
library("ggnewscale")







# Read raw data ----
#------------------------------------------------#
#                                                #
#                   READ THE DATA                # 
#                                                #
#------------------------------------------------# 

# Read data
#----------------------------------#
WS_result_df <- readRDS("WS_result_df.rds") %>% dplyr::mutate(rarefied = F)
WS_R_result_df <- readRDS("WS_R_result_df.rds") %>% dplyr::mutate(rarefied = T)
WS_F_result_df <- readRDS("WS_F_result_df.rds") %>% dplyr::mutate(rarefied = F)
WS_FR_result_df <- readRDS("WS_FR_result_df.rds") %>% dplyr::mutate(rarefied = T)
EstMB_result_df <- readRDS("EstMB_result_df.rds") %>% dplyr::mutate(rarefied = F)
EstMB_R_result_df <- readRDS("EstMB_R_result_df.rds") %>% dplyr::mutate(rarefied = T)

# Combine results
ML_results_WS <- dplyr::bind_rows(WS_result_df, WS_R_result_df, WS_F_result_df, WS_FR_result_df, EstMB_result_df, EstMB_R_result_df) %>% 
  dplyr::mutate(filtering = case_when(is.na(filtering) == TRUE ~ "full", 
                                      str_detect(filtering, "filteredRarefied") == TRUE ~ paste("0.", substr(filtering, 17, nchar(filtering)), sep = ""), 
                                      str_detect(filtering, "filtered") == TRUE ~ paste("0.", substr(filtering, 9, nchar(filtering)), sep = "")), 
                algoritm = ifelse(algoritm == "LASSO", "ENET", algoritm))



# First data transformations
#----------------------------------#
WS_results_final <- ML_results_WS %>% 
  dplyr::mutate(transformation2 = ifelse(substr(transformation, 1, 3) == "ALR", "ALR", as.character(transformation)),
                transformation2 = factor(transformation2, levels = c("PA", "TSS", "logTSS", "aSIN", "rCLR", "CLR", "ALR", "ILR")),
                target = case_when(study %in% c("LifeLinesDeep_2016", "AsnicarF_2021", "LiJ_2014", "LeChatelierE_2013", "EstBMI30metaphlan_2022") ~ "BMI", 
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
                algoritm = ifelse(algoritm == "LASSO", "ENET", algoritm)) %>% 
  dplyr::filter(!(study %in% c("HanniganGD_2017", "ThomasAM_2018a", "ThomasAM_2018b", 
                               "ThomasAM_2019c", "GuptaA_2019", "FengQ_2015", "KeohaneDM_2020"))) %>% 
  dplyr::filter(set == "test") %>% 
  dplyr::rename("fold" = "seed", 
                "recall" = "sensitivity", 
                "algorithm" = "algoritm")






# Analyses ----
#------------------------------------------------#
#                                                #
#     RUN ANALYSES TO COMPARE THE PERFORMANCE    # 
#                                                #
#------------------------------------------------#   

# Paired T-test for WS setting 
#----------------------------------#
transformation_list <- c("PA", "aSIN", "TSS", "logTSS", "CLR", "rCLR", "ALR", "ILR")

calculate_paired_ttests <- function(dataset){
  paired_ttest_df <- data.frame()
  for (i in 1:length(transformation_list)){
    for (j in 1:length(transformation_list)){
      
      if (i != j){
        
        run_input = dataset %>%
          dplyr::filter(transformation2 %in% c(transformation_list[i], transformation_list[j])) %>% 
          dplyr::mutate(trans_help = ifelse(transformation2 == transformation_list[j], "A", "B")) %>% 
          dplyr::select(-one_of("transformation2")) %>% 
          tidyr::spread(trans_help, AUC)
        
        t_test_raw = t.test(run_input$A[run_input$rarefied == FALSE], run_input$B[run_input$rarefied == FALSE], paired = TRUE)
        w_test_raw = wilcox.test(run_input$A[run_input$rarefied == FALSE], run_input$B[run_input$rarefied == FALSE], paired = TRUE)
        
        t_test_rare = t.test(run_input$A[run_input$rarefied == TRUE], run_input$B[run_input$rarefied == TRUE], paired = TRUE)
        w_test_rare = wilcox.test(run_input$A[run_input$rarefied == TRUE], run_input$B[run_input$rarefied == TRUE], paired = TRUE)
        
        run_df1 = data.frame(tr1 = transformation_list[i], 
                            tr2 = transformation_list[j], 
                            mean_tr2_minus_tr1 = t_test_raw$estimate,
                            pval_t = t_test_raw$p.value,
                            pval_w= w_test_raw$p.value,
                            rarefied = FALSE)
        
        run_df2 = data.frame(tr1 = transformation_list[i], 
                             tr2 = transformation_list[j], 
                             mean_tr2_minus_tr1 = t_test_rare$estimate,
                             pval_t = t_test_rare$p.value,
                             pval_w = w_test_rare$p.value,
                             rarefied = TRUE)
        
        paired_ttest_df = dplyr::bind_rows(paired_ttest_df, run_df1, run_df2)
      }
    }
  }
  return(paired_ttest_df)
}


# Different representations - averaging etc
#----------------------------------#

# target weighted datasets by algorithm
WS_results_targetWeighted_RF <- WS_results_final %>% 
  dplyr::filter(filtering == "full") %>% 
  dplyr::filter(algorithm == "RF") %>% 
  dplyr::group_by(rarefied, transformation2, target) %>% 
  dplyr::summarise(AUC = mean(AUC)) %>% 
  dplyr::ungroup()

WS_results_targetWeighted_ENET <- WS_results_final %>% 
  dplyr::filter(filtering == "full") %>% 
  dplyr::filter(algorithm == "ENET") %>% 
  dplyr::group_by(rarefied, transformation2, target) %>% 
  dplyr::summarise(AUC = mean(AUC)) %>% 
  dplyr::ungroup()

WS_results_targetWeighted_XGB <- WS_results_final %>% 
  dplyr::filter(filtering == "full") %>% 
  dplyr::filter(algorithm == "XGB") %>% 
  dplyr::group_by(rarefied, transformation2, target) %>% 
  dplyr::summarise(AUC = mean(AUC)) %>% 
  dplyr::ungroup()

# Run t-test
ttest_paired_targetWeighted_RF <- calculate_paired_ttests(WS_results_targetWeighted_RF)
ttest_paired_targetWeighted_ENET <- calculate_paired_ttests(WS_results_targetWeighted_ENET)
ttest_paired_targetWeighted_XGB <- calculate_paired_ttests(WS_results_targetWeighted_XGB)





# Paired T-test between algorithms
#----------------------------------#
ttest_algos_targetweighted <- WS_results_final %>% 
  dplyr::filter(filtering == "full") %>% 
  dplyr::group_by(rarefied, target, algorithm, transformation2) %>%
  dplyr::summarise(mean_AUC = mean(AUC)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(algorithm, mean_AUC)

algos <- c("RF", "ENET", "XGB")

algo_comp_df = data.frame()
for (i in transformation_list){
  
  for  (j in c(TRUE, FALSE)){
    
    data = ttest_algos_targetweighted %>% 
      dplyr::filter(transformation2 == i) %>% 
      dplyr::filter(rarefied == j)
    
    for (algo1 in 2:length(algos)){
      for (algo2 in 1:(algo1-1)){
        weighted_w = wilcox.test(x = data %>% dplyr::pull(algos[algo1]), y = data %>% dplyr::pull(algos[algo2]), paired = T)
        weighted_t = t.test(x = data %>% dplyr::pull(algos[algo1]), y = data %>% dplyr::pull(algos[algo2]), paired = T)
        
        run_df <- data.frame(algorithm1 = algos[algo1], 
                             algorithm2 = algos[algo2], 
                             transformation = i, 
                             rarefied = j, 
                             est_diff = -weighted_t$estimate,
                             t_pval = weighted_t$p.value,
                             w_pval = weighted_w$p.value)
        
        algo_comp_df = dplyr::bind_rows(algo_comp_df, run_df)
      }
    }
  }
}

algo_comp_final <- algo_comp_df %>% 
  dplyr::select(algorithm1, algorithm2, transformation, rarefied, est_diff, w_pval) %>%
  dplyr::group_by(rarefied) %>% 
  dplyr::mutate(w_FDR = p.adjust(w_pval, method = "BH")) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(rarefied, transformation, algorithm1, algorithm2) %>% 
  as.data.frame()

#write.csv(algo_comp_final, "/Bioinf/BenchmarkML/Table 2. Algorithm comparison.csv", row.names = F, col.names = T)




# Paired T-test between rarefied vs not-rarefied
#----------------------------------#
ttest_rarefied_targetweighted <- WS_results_final %>% 
  dplyr::filter(filtering == "full") %>% 
  dplyr::group_by(rarefied, target, algorithm, transformation2) %>%
  dplyr::summarise(mean_AUC = mean(AUC)) %>% 
  dplyr::ungroup() %>% 
  tidyr::spread(rarefied, mean_AUC)

rarefied_comp_df = data.frame()
for (i in transformation_list){
  
  data <- ttest_rarefied_targetweighted %>% 
    dplyr::filter(transformation2 == i)

  weighted_w = wilcox.test(x = data %>% dplyr::pull(`FALSE`), 
                           y = data %>% dplyr::pull(`TRUE`), paired = T)
  
  weighted_t = t.test(x = data %>% dplyr::pull(`FALSE`), 
                      y = data %>% dplyr::pull(`TRUE`), paired = T)
        
  run_df <- data.frame(transformation = i, 
                       est_diff = -weighted_t$estimate,
                       t_pval = weighted_t$p.value,
                       w_pval = weighted_w$p.value)
        
  rarefied_comp_df = dplyr::bind_rows(rarefied_comp_df, run_df)
}


rarefied_comp_final <- rarefied_comp_df %>% 
  dplyr::mutate(w_FDR = p.adjust(w_pval, method = "BH"))







# Figures ----
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#   

# Within-study results in general
#----------------------------------#
plotData1 <- WS_results_final %>% 
  dplyr::filter(filtering == "full") %>% 
  dplyr::select(study, rarefied, algorithm, transformation2, filtering, fold, target, AUC, precision) %>% 
  tidyr::pivot_longer(cols = c("AUC", "precision")) %>% 
  dplyr::mutate(study = case_when(study == "EstAB90metaphlan_2022" ~ "EstMB_AB90_2022",
                                  study == "EstBMI30metaphlan_2022" ~ "EstMB_BMI30_2022",
                                  study == "EstDepressionmetaphlan_2022" ~ "EstMB_Depression_2022",
                                  study == "EstT2Dmetaphlan_2022" ~ "EstMB_T2D_2022",
                                  TRUE ~ study)) 

plotData1_aggreg <- plotData1 %>% 
  dplyr::group_by(study, rarefied, algorithm, transformation2, filtering, target, name) %>% 
  dplyr::summarise(value = mean(value))

p_WS_main <- ggplot(plotData1_aggreg %>% dplyr::filter(rarefied == FALSE), 
                    aes(x = transformation2, y = value, fill = algorithm)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#2b8cbe", "orange", "azure4"), name = "") + 
  xlab("") + 
  ylab("AUROC") + 
  theme_classic() +
  facet_wrap(vars(factor(name, levels = c("AUC", "precision"), labels = c("AUROC", "Precision"))), scales = "free", ncol = 2) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "right",
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 16))

p_WS_main
ggsave(plot = p_WS_main, "/Bioinf/BenchmarkML/Figures/p_WS_main.png", width = 16, height = 5)
ggsave(plot = p_WS_main, "/Bioinf/BenchmarkML/Figures/p_WS_main.pdf", width = 16, height = 5)




p_WS_main_rarefied <- ggplot(plotData1_aggreg %>% dplyr::filter(rarefied == TRUE), 
                    aes(x = transformation2, y = value, fill = algorithm)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("#2b8cbe", "orange", "azure4"), name = "") + 
  xlab("") + 
  ylab("AUROC") + 
  theme_classic() +
  facet_wrap(vars(factor(name, levels = c("AUC", "precision"), labels = c("AUROC", "Precision"))), scales = "free", ncol = 2) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "right",
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 16))

p_WS_main_rarefied
ggsave(plot = p_WS_main_rarefied, "/Bioinf/BenchmarkML/Figures/p_WS_main_rarefied.png", width = 16, height = 5)
ggsave(plot = p_WS_main_rarefied, "/Bioinf/BenchmarkML/Figures/p_WS_main_rarefied.pdf", width = 16, height = 5)




# Within-study results by study
#----------------------------------#
p_WS_byStudy_AUC <- ggplot(plotData1 %>% filter(rarefied == FALSE), 
                           aes(x = transformation2, y = value, fill = algorithm)) + 
  geom_boxplot() + 
  ylab("AUROC") +
  scale_fill_manual(values = c("#2b8cbe", "orange", "azure4"), name = "") + 
  facet_wrap(vars(study), scales = "free", ncol = 3) + 
  xlab("") + 
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        strip.text = element_text(size = 14), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        panel.grid.major = element_line())

p_WS_byStudy_AUC
ggsave(plot = p_WS_byStudy_AUC, "/Bioinf/BenchmarkML/Figures/p_WS_byStudy_AUC.png", width = 16, height = 16)
ggsave(plot = p_WS_byStudy_AUC, "/Bioinf/BenchmarkML/Figures/p_WS_byStudy_AUC.pdf", width = 16, height = 16)

p_WS_byStudy_AUC_rarefied <- ggplot(plotData1 %>% filter(rarefied == TRUE), 
                           aes(x = transformation2, y = value, fill = algorithm)) + 
  geom_boxplot() + 
  ylab("AUROC") +
  scale_fill_manual(values = c("#2b8cbe", "orange", "azure4"), name = "") + 
  facet_wrap(vars(study), scales = "free", ncol = 3) + 
  xlab("") + 
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        strip.text = element_text(size = 14), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        panel.grid.major = element_line())

p_WS_byStudy_AUC_rarefied
ggsave(plot = p_WS_byStudy_AUC_rarefied, "/Bioinf/BenchmarkML/Figures/p_WS_byStudy_AUC_rarefied.png", width = 16, height = 16)
ggsave(plot = p_WS_byStudy_AUC_rarefied, "/Bioinf/BenchmarkML/Figures/p_WS_byStudy_AUC_rarefied.pdf", width = 16, height = 16)







# Heatmaps ----
#------------------------------------------------#
#                                                #
#                    HEATMAPS                    # 
#                                                #
#------------------------------------------------#  

# P-values heatmap for WS setting
#----------------------------------#

# Make a combined matrix of p-values and t-values
help_df1 <- data.frame(tr1 = transformation_list, tr2 = transformation_list, mean_tr2_minus_tr1 = 0)
help_df_RF <- help_df1 %>% dplyr::mutate(algorithm = "RF")
help_df_ENET <- help_df1 %>% dplyr::mutate(algorithm = "ENET")
help_df_XGB <- help_df1 %>% dplyr::mutate(algorithm = "XGB")


# Help
ttest_paired_targetWeighted_RF <- ttest_paired_targetWeighted_RF %>% 
  dplyr::mutate(algorithm = "RF")

ttest_paired_targetWeighted_ENET <- ttest_paired_targetWeighted_ENET %>% 
  dplyr::mutate(algorithm = "ENET")

ttest_paired_targetWeighted_XGB <- ttest_paired_targetWeighted_XGB %>% 
  dplyr::mutate(algorithm = "XGB")

# Calculate FDR
FDR_data <- dplyr::bind_rows(ttest_paired_targetWeighted_ENET, ttest_paired_targetWeighted_RF, ttest_paired_targetWeighted_XGB) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(pair = paste(sort(c(tr1, tr2)), collapse = "_")) %>% 
  dplyr::distinct(rarefied, pair, algorithm, pval_w ) %>% 
  dplyr::ungroup() %>% 
  dplyr::distinct(rarefied, pair, algorithm, pval_w ) %>% 
  dplyr::group_by(rarefied) %>% 
  dplyr::mutate(FDR = p.adjust(pval_w, method = "BH")) %>% 
  dplyr::ungroup()

# Plot data for the heatmap
plotData2 <- ttest_paired_targetWeighted_ENET %>% 
  dplyr::bind_rows(ttest_paired_targetWeighted_RF) %>% 
  dplyr::bind_rows(ttest_paired_targetWeighted_XGB) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(pair = paste(sort(c(tr1, tr2)), collapse = "_")) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(FDR_data[ ,c("pair", "algorithm", "FDR", "rarefied")], by = c("rarefied", "pair", "algorithm")) %>% 
  dplyr::bind_rows(help_df_RF, help_df_ENET, help_df_XGB) %>% 
  dplyr::mutate(signficance = case_when(FDR <= 0.05 ~ "**",
                                        pval_w <= 0.05 ~ "*",
                                        TRUE ~ ""),
                print = case_when(tr1 == tr2 ~ "-",
                                  TRUE ~ paste(round(mean_tr2_minus_tr1, 3), signficance, sep = "")))

p_ttest_heatmap <- ggplot(plotData2 %>% dplyr::filter(rarefied == FALSE), 
                          aes(x = factor(tr1, levels = c(c("PA", "TSS", "logTSS", "aSIN", "rCLR", "CLR", "ALR", "ILR"))), 
                              y= factor(tr2, levels = c(c("PA", "TSS", "logTSS", "aSIN", "rCLR", "CLR", "ALR", "ILR"))),
                              fill = mean_tr2_minus_tr1)) +
  geom_tile(color = "white") + 
  geom_text(aes(label = print), size = 4) +
  scale_fill_gradient2_tableau(palette = "Orange-Blue Diverging",
                               trans = "reverse",
                               name = "") + 
  xlab("Transformation 1") + 
  ylab("Transformation 2") + 
  facet_grid(. ~ algorithm, scales = "free", space = "free") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 14), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 28), 
        legend.key.size = unit(0.9, "cm"), 
        legend.text = element_text(size = 12), 
        legend.position = "bottom")
p_ttest_heatmap

ggsave(plot = p_ttest_heatmap, "/Bioinf/BenchmarkML/Figures/p_comparison_heatmap.png", width = 18, height = 7)
ggsave(plot = p_ttest_heatmap, "/Bioinf/BenchmarkML/Figures/p_comparison_heatmap.pdf", width = 18, height = 7)



# Rarefied
p_ttest_heatmap_rarefied <- ggplot(plotData2 %>% dplyr::filter(rarefied == TRUE), 
                          aes(x = factor(tr1, levels = c(c("PA", "TSS", "logTSS", "aSIN", "rCLR", "CLR", "ALR", "ILR"))), 
                              y= factor(tr2, levels = c(c("PA", "TSS", "logTSS", "aSIN", "rCLR", "CLR", "ALR", "ILR"))),
                              fill = mean_tr2_minus_tr1)) +
  geom_tile(color = "white") + 
  geom_text(aes(label = print), size = 4) +
  scale_fill_gradient2_tableau(palette = "Orange-Blue Diverging",
                               trans = "reverse",
                               name = "") + 
  xlab("Transformation 1") + 
  ylab("Transformation 2") + 
  facet_grid(. ~ algorithm, scales = "free", space = "free") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 16), 
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 28), 
        legend.key.size = unit(1, "cm"), 
        legend.text = element_text(size = 14), 
        legend.position = "bottom")
p_ttest_heatmap_rarefied

ggsave(plot = p_ttest_heatmap_rarefied, "/Bioinf/BenchmarkML/Figures/p_comparison_heatmap_rarefied.png", width = 18, height = 7)
ggsave(plot = p_ttest_heatmap_rarefied, "/Bioinf/BenchmarkML/Figures/p_comparison_heatmap_rarefied.pdf", width = 18, height = 7)







# Within-study results in by feature filtering
#----------------------------------#
WS_filtered_plotData <- WS_results_final %>% 
  dplyr::mutate(filtering = factor(filtering, 
                                   levels = rev(c("full", "0.10", "0.25", "0.50", "0.75", "0.90")),
                                   labels = rev(c("No filtering", "10% prevalence", "25% prevalence", "50% prevalence", 
                                                  "75% prevalence", "90% prevalence")))) %>% 
  dplyr::group_by(rarefied, transformation2, algorithm, filtering) %>%  
  dplyr::summarise(value = mean(AUC, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(key = "AUC")


WS_filtered_plotData_AUC <- WS_results_final %>% 
  dplyr::mutate(filtering = factor(filtering, 
                                   levels = rev(c("full", "0.10", "0.25", "0.50", "0.75", "0.90")),
                                   labels = rev(c("No filtering", "10% prevalence", "25% prevalence", "50% prevalence", 
                                                  "75% prevalence", "90% prevalence")))) %>% 
  dplyr::mutate(AUC = ifelse(is.na(AUC_old) == FALSE, AUC_old, AUC)) %>% 
  dplyr::group_by(rarefied, transformation2, algorithm, filtering) %>%  
  dplyr::summarise(value = mean(AUC, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(key = "AUC")
  
p_WS_filtering <- ggplot(WS_filtered_plotData_AUC %>% dplyr::filter(rarefied == F), 
                         aes(x = filtering, y = value, color = transformation2, group = transformation2)) + 
  geom_line(size = 2, alpha = 0.8) + 
  geom_point(size = 5, alpha = 0.8) + 
  scale_color_calc(name = "") + 
  facet_wrap(vars(algorithm), nrow = 1, scales = "free") + 
  xlab("") + 
  ylab("") + 
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 20),
        legend.position = "right",
        strip.background = element_rect(fill = "white"), 
        legend.text = element_text(size = 14))

p_WS_filtering
ggsave(plot = p_WS_filtering, "/Bioinf/BenchmarkML/Figures/p_WS_filtering.png", width = 18, height = 8)
ggsave(plot = p_WS_filtering, "/Bioinf/BenchmarkML/Figures/p_WS_filtering.pdf", width = 18, height = 8)

