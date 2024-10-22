

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------# 

# Read the packages
library("tidyverse")
library("ggsci")
library("pheatmap")
library("ggthemes")
library("readr")
library("viridis")
library("ggpubr")

# Read raw data 
cumulative_results <- readRDS("cumulative_results.rds")




# First data transformations
#---------------------------------------#
cumulative_results_final <- cumulative_results %>% 
  dplyr::filter(set == "test")%>% 
  dplyr::mutate(target = case_when(filename == "EstAB90metaphlan_2022_PA.csv" ~ "AB",
                                   filename == "EstBMI30metaphlan_2022_PA.csv" ~ "BMI",
                                   filename == "EstDepressionmetaphlan_2022_PA.csv" ~ "Depression"),
                algorithm = ifelse(algorithm == "LASSO", "ENET", algorithm))

cumulative_results_final %>% 
  dplyr::group_by(algorithm, subset_size) %>% 
  dplyr::summarise(mean_AUC = mean(AUC))






# Figures ----
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#   

# Baseline AUC
baseline_AUC_df <- EstMB_result_df %>% 
  dplyr::filter(set == "test") %>% 
  dplyr::filter(transformation == "PA") %>% 
  dplyr::filter(rarefied == FALSE) %>% 
  dplyr::filter(study %in% c("EstAB90metaphlan_2022", "EstBMI30metaphlan_2022", "EstDepressionmetaphlan_2022")) %>%
  dplyr::group_by(study, algoritm) %>% 
  dplyr::summarise(mean_AUC = mean(AUC)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(target = case_when(study == "EstAB90metaphlan_2022" ~ "AB",
                                   study == "EstBMI30metaphlan_2022" ~ "BMI",
                                   study == "EstDepressionmetaphlan_2022" ~ "Depression"))




# Cumulative classifier plot
#----------------------------------#
p_cumulativePA_AUC <- ggplot(cumulative_results_final, 
                             aes(x = factor(subset_size), 
                                 y = AUC, 
                                 fill = algorithm)) + 
  geom_hline(data = baseline_AUC_df %>%
               dplyr::filter(algoritm == "LASSO"),
             aes(yintercept = mean_AUC), linetype = 3, linewidth = 1, color = "#2b8cbe") +
  geom_hline(data = baseline_AUC_df %>%
               dplyr::filter(algoritm == "RF"),
             aes(yintercept = mean_AUC), linetype = 3, linewidth = 1, color = "orange") +
  geom_hline(data = baseline_AUC_df %>%
               dplyr::filter(algoritm == "XGB"),
             aes(yintercept = mean_AUC), linetype = 3, linewidth = 1, color = "azure4") +
  geom_boxplot() + 
  ylab("AUROC") +
  xlab("Number of selected features") + 
  theme_classic() +
  scale_fill_manual(name = "", values = c("#2b8cbe", "orange", "azure4")) + 
  facet_wrap(vars(factor(target, levels = c("AB", "BMI", "T2D", "Depression"), 
                         labels = c("Antibiotics usage last 90 days", "Body mass index", "Type 2 diabetes", "Depression"))), scales = "free") + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20), 
        legend.position = "bottom",
        strip.text = element_text(size = 20),
        panel.grid.major = element_line())

p_cumulativePA_AUC
ggsave(plot = p_cumulativePA_AUC, "Bioinf/BenchmarkML/Figures/p_cumulativePA_AUC.png", width = 20, height = 7)
ggsave(plot = p_cumulativePA_AUC, "Bioinf/BenchmarkML/Figures/p_cumulativePA_AUC.pdf", width = 20, height = 7)

