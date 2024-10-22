

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


# Read results
ALRsensitivity_results <- readRDS("ALRsensitivity_results.rds")







# Figures ----
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#   


# Variability of ALR results dependent on the number of reference elements chosen
#----------------------------------#
ALR_sensitivity_data <- ALRsensitivity_results %>% 
  dplyr::filter(set == "test")

ALR_sensitivity_res <- data.frame()
for (i in 1:100){
  for (j in 1:20){
    sample <- sample(x = 1:20, size = j, replace = F)
    
    ALR_perf <- ALR_sensitivity_data %>% 
      dplyr::filter(transformation %in% paste("ALR", sample, sep = "")) %>% 
      dplyr::group_by(algoritm) %>% 
      dplyr::summarize(mean_AUC = mean(AUC)) %>% 
      dplyr::ungroup()
    
    run_df <- ALR_perf %>% 
      dplyr::mutate(n_elements = j, 
                    seed = i)
    
    ALR_sensitivity_res = dplyr::bind_rows(ALR_sensitivity_res, run_df)
  }
}

mean_LASSO = ALR_sensitivity_res %>% dplyr::filter(algoritm == "LASSO") %>% dplyr::summarise(mean_AUC = mean(mean_AUC))
mean_RF = ALR_sensitivity_res %>% dplyr::filter(algoritm == "RF") %>% dplyr::summarise(mean_AUC = mean(mean_AUC))
mean_XGB = ALR_sensitivity_res %>% dplyr::filter(algoritm == "XGB") %>% dplyr::summarise(mean_AUC = mean(mean_AUC))

mean_df <- data.frame(algoritm = c("LASSO", "RF", "XGB"), 
                      mean = c(mean_LASSO$mean_AUC, mean_RF$mean_AUC, mean_XGB$mean_AUC))

p_ALRsensitivity <- ggplot(ALR_sensitivity_res, aes(x = factor(n_elements), y = mean_AUC, fill = algoritm)) + 
  geom_boxplot() + 
  scale_fill_manual(name = "", values = c("#2b8cbe", "orange", "azure4"), guide = F) + 
  xlab("Number of different reference elements for averaging") +
  ylab("Mean AUROC") + 
  geom_hline(data = mean_df, aes(yintercept = mean), color = "red", linetype = 2, size = 1) + 
  facet_wrap(vars(algoritm), ncol = 3, scales = "free") + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 16))
p_ALRsensitivity

ggsave(plot = p_ALRsensitivity, "/BenchmarkML/Figures/p_ALRsensitivity.png", width = 20, height = 6)
ggsave(plot = p_ALRsensitivity, "/BenchmarkML/Figures/p_ALRsensitivity.pdf", width = 20, height = 6)

