

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








# Read raw data ----
#------------------------------------------------#
#                                                #
#                   READ THE DATA                # 
#                                                #
#------------------------------------------------# 

# First data transformations
#---------------------------------------#
LOSO_result_df <- readRDS("LOSO_result_df.rds")

LOSO_result_final <- LOSO_result_df %>% 
  dplyr::filter(set == "test") %>% 
  dplyr::mutate(transformation2 = ifelse(substr(transformation, 1, 3) == "ALR", "ALR", as.character(transformation)),
                transformation2 = factor(transformation2, levels = c("PA", "TSS", "logTSS", "aSIN", "ANCOMBC", "rCLR", "CLR", "ALR", "ILR")),
                algoritm = ifelse(algoritm == "LASSO", "ENET", algoritm))







# Figures ----
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#   

# LOSO filtered
#----------------------------------#
LOSO_plotData <- LOSO_result_final %>% 
  dplyr::group_by(transformation2, LOSOset, algoritm, target) %>% 
  dplyr::summarize(mean_AUC = mean(AUC), 
                   mean_precision = mean(precision, na.rm = T)) 


p_LOSO <- ggplot(LOSO_plotData,
                 aes(x = transformation2, y = mean_AUC, fill = algoritm)) + 
  geom_boxplot() + 
  facet_wrap(vars(factor(target, levels = c("CRCfiltered", "BMIfiltered"), 
                         labels = c("Colorectal cancer (LOSO-CV)", "Body mass index (LOSO-CV)"))), scales = "free", ncol = 2) + 
  scale_fill_manual(values = c("#2b8cbe", "orange", "azure4"), name = "") + 
  xlab("") + 
  ylab("AUROC") + 
  theme_classic() + 
  theme(strip.text = element_text(size = 20),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))
p_LOSO

ggsave(plot = p_LOSO, "/Bioinf/BenchmarkML/Figures/p_LOSOfiltered.png", width = 16, height = 6)
ggsave(plot = p_LOSO, "/Bioinf/BenchmarkML/Figures/p_LOSOfiltered.pdf", width = 16, height = 6)


