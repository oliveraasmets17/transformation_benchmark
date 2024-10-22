

# Program setup phase ----
#------------------------------------------------#
#                                                #
#              PROGRAM SETUP PHASE               # 
#                                                #
#------------------------------------------------# 

# Read the packages
library("tidyverse")
library("ggsci")
library("ggthemes")
library("ggpubr")


# Read raw data 
grid_result_df <- readRDS("grid_result_df.rds")







# Figures ----
#------------------------------------------------#
#                                                #
#                     FIGURES                    # 
#                                                #
#------------------------------------------------#   

# Grid figure for AB usage and effect of sample size and feature number
#----------------------------------#
grid_plot_avg <- grid_result_df %>% 
  dplyr::filter(set == "test") %>% 
  dplyr::mutate(transformation2 = ifelse(substr(transformation, 1, 3) == "ALR", "ALR", as.character(transformation)),
                transformation2 = factor(transformation2, levels = c("PA", "TSS", "logTSS", "aSIN", "rCLR", "CLR", "ALR", "ILR")),
                algoritm = ifelse(algoritm == "LASSO", "ENET", algoritm)) %>% 
  dplyr::group_by(target, subsampling, bioinf, transformation2, algoritm) %>% 
  dplyr::summarise(mean_AUC = mean(AUC)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(subsampling = factor(subsampling, 
                                     levels = c("sub20", "sub40", "sub60", "sub80", "sub100"),
                                     labels = c("20%", "40%", "60%", "80%", "100%")),
                bioinf = factor(bioinf, 
                                levels = c("1p", "10p", "50p"), 
                                labels = c("1% prevalence", "10% prevalence", "50% prevalence"))) %>% 
  dplyr::filter(is.na(subsampling) == FALSE)

p_grid <- ggplot(data = grid_plot_avg) + 
  geom_line(data = grid_plot_avg %>% filter(bioinf == "1% prevalence"), 
            aes(x = subsampling, y = mean_AUC, color = algoritm, group = paste(bioinf, algoritm, sep = "_")), linetype = 3, size = 1, alpha = 0.8) + 
  geom_line(data = grid_plot_avg %>% filter(bioinf == "10% prevalence"), 
            aes(x = subsampling, y = mean_AUC, color = algoritm, group = paste(bioinf, algoritm, sep = "_")), linetype = 1, size = 1, alpha = 0.8) + 
  geom_line(data = grid_plot_avg %>% filter(bioinf == "50% prevalence"), 
            aes(x = subsampling, y = mean_AUC, color = algoritm, group = paste(bioinf, algoritm, sep = "_")), linetype = 6, size = 1, alpha = 0.8) + 
  facet_grid(target ~ transformation2, scales = "free", shrink = T) + 
  geom_point(data = grid_plot_avg, 
             aes(x = subsampling, y = mean_AUC, color = algoritm, group = paste(bioinf, algoritm, sep = "_"), shape = bioinf), size = 3) + 
  scale_color_calc(name = "Algorithm") + 
  scale_shape_manual(name = "Feature filtering", values = c(16, 17, 18)) + 
  ylab("AUROC") + 
  xlab("Sample subsetting") + 
  theme_bw() + 
  theme(legend.position = "right", 
        axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16),
        panel.grid.major.y = element_line(), 
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 16), 
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12))
p_grid

ggsave(plot = p_grid, "Bioinf/BenchmarkML/Figures/p_grid.png", width = 18, height = 8)
ggsave(plot = p_grid, "Bioinf/BenchmarkML/Figures/p_grid.pdf", width = 18, height = 8)

