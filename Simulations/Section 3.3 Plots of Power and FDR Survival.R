# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Read CSV files
df_simulations_results_1 <- read.csv("df_nondvine_results_Survival_nonGaussian_X.csv")
df_simulations_results_2 <- read.csv("df_dvine_results_Survival_nonGaussian_X.csv")
df_simulations_results_3 <- read.csv("df_noncvine_results_Survival_nonGaussian_X.csv")
df_simulations_results_4 <- read.csv("df_cvine_results_Survival_nonGaussian_X.csv")
df_simulations_results_5 <- read.csv("df_second_order_results_Survival_nonGaussian_X.csv")
df_simulations_results_6 <- read.csv("df_sequential_results_Survival_nonGaussian_X.csv")
df_simulations_results_7 <- read.csv("df_lasso_results_Survival_nonGaussian_X.csv")
df_simulations_results_8 <- read.csv("df_SCAD_results_Survival_nonGaussian_X.csv")

#Renaming
colnames(df_simulations_results_1) <- c("Censoring","Iteration","nonparDvine_power","nonparDvine_FDP" )
colnames(df_simulations_results_2) <- c("Censoring","Iteration","parDvine_power","parDvine_FDP" )
colnames(df_simulations_results_3) <- c("Censoring","Iteration","nonparCvine_power","nonparCvine_FDP" )
colnames(df_simulations_results_4) <- c("Censoring","Iteration","parCvine_power","parCvine_FDP" )
colnames(df_simulations_results_5) <- c("Censoring","Iteration","2nd_order_power","2nd_order_FDP" )
colnames(df_simulations_results_6) <- c("Censoring","Iteration","Sequential_power","Sequential_FDP")
colnames(df_simulations_results_7) <- c("Censoring","Iteration","lasso_power","lasso_FDP" )
colnames(df_simulations_results_8) <- c("Censoring","Iteration","SCAD_power","SCAD_FDP" )

# Concatenate data frames
df_simulations_results <- cbind(df_simulations_results_1, df_simulations_results_2[,3:4], 
                                    df_simulations_results_3[,3:4],  df_simulations_results_4[,3:4],
                                df_simulations_results_5[,3:4], df_simulations_results_6[,3:4],
                                df_simulations_results_7[,3:4], df_simulations_results_8[,3:4])

#View(df_simulations_results)

results_power <- df_simulations_results %>%
  group_by(Censoring) %>%
  summarise(across(ends_with("er"), mean, na.rm = TRUE))
#View(results_power)

results_FDR <- df_simulations_results %>%
  group_by(Censoring) %>%
  summarise(across(ends_with("FDP"), mean, na.rm = TRUE))
#View(results_FDR)

colnames(results_FDR) <- c("Censoring", "nonparDvine_FDR", "parDvine_FDR", "nonparCvine_FDR","parCvine_FDR", "2nd_order_FDR","Sequential_FDR", "lasso_FDR", "SCAD_FDR")  


# Plot Power and FDP side-by-side
library(gridExtra)
library(reshape2)
library(tidyr)

# Melt data for ggplot2 compatibility (power)
df_melted_results_power <- melt(results_power, id.vars = "Censoring", variable.name = "Method", value.name = "Power")
View(df_melted_results_power)

# Define custom colors and labels
custom_colors_power <- c("nonparDvine_power" = "chartreuse4", "parDvine_power" = "deepskyblue3", 
                   "2nd_order_power" = "orange", "Sequential_power" = "darkmagenta",
                   "lasso_power"="red3","SCAD_power"="burlywood4",
                   "nonparCvine_power" = "palegreen2", "parCvine_power" = "deepskyblue1")


custom_labels <- c("nonparDvine_power" = "Non-DTLDCKe", 
                   "parDvine_power" = "DTLDCKe", 
                   "2nd_order_power" = "2nd Order",
                   "Sequential_power" = "Sequential",
                   "lasso_power"="LASSO",
                   "SCAD_power"="SCAD",
                   "nonparCvine_power" = "Non-DLCCKe", 
                   "parCvine_power" = "DLCCKe")


# Plotting
power_plot<- ggplot(df_melted_results_power, aes(x = Censoring, y = Power, color=Method, shape = Method)) +
  geom_point(size = 2) +
  geom_line() +
  scale_color_manual(values = custom_colors_power,labels = custom_labels) +  # Custom colors
  scale_shape_manual(values = c(16, 17, 15, 18, 19, 20, 16, 17), labels = custom_labels) +  # Custom shapes and labels
  scale_y_continuous(limits = c(0, 100), breaks = seq(0,100, by=10))+
  scale_x_continuous(breaks = results_power$Censoring) +
  labs(x = "Censoring rate", y = "Average Power (%)") +
  theme_light() +
  theme(plot.title = element_text(size = 18, hjust = 0.5))

power_plot



# Melt data for ggplot2 compatibility (FDR)
df_melted_results_FDR <- melt(results_FDR, id.vars = "Censoring", variable.name = "Method", value.name = "FDR")
View(df_melted_results_FDR)

# Define custom colors and labels
custom_colors_FDR <- c("nonparDvine_FDR" = "chartreuse4", "parDvine_FDR" = "deepskyblue3", 
                   "2nd_order_FDR" = "orange","Sequential_FDR" = "darkmagenta", 
                   "lasso_FDR"="red3","SCAD_FDR"="burlywood4",
                   "nonparCvine_FDR" = "palegreen2", "parCvine_FDR" = "deepskyblue1")

custom_labels <- c("nonparDvine_FDR" = "Non-DTLDCKe", 
                   "parDvine_FDR" = "DTLDCKe", 
                   "2nd_order_FDR" = "2nd Order",
                   "Sequential_FDR" = "Sequential",
                   "lasso_FDR"="LASSO",
                   "SCAD_FDR"="SCAD",
                   "nonparCvine_FDR" = "Non-DLCCKe", 
                   "parCvine_FDR" = "DLCCKe")

# Plotting
fdr_plot<- ggplot(df_melted_results_FDR, aes(x = Censoring, y = FDR, color=Method, shape = Method)) +
  geom_point(size = 2) +
  geom_hline(yintercept=20, linetype="dashed", color = "darkslategrey")+
  geom_line() +
  scale_color_manual(values = custom_colors_FDR,labels = custom_labels) +  # Custom colors
  scale_shape_manual(values = c(16, 17, 15, 18, 19, 20, 16, 17), labels = custom_labels) +  # Custom shapes and labels
  scale_y_continuous(limits = c(0, 100), breaks = seq(10,100, by=10)) +
  scale_x_continuous(breaks = results_FDR$Censoring) +
  labs(x = "Censoring rate", y = "FDR (%)") +
  theme_light() +
  theme(plot.title = element_text(size = 18, hjust = 0.5))

fdr_plot2 <- fdr_plot + theme(legend.position = "none")
fdr_plot2

grid.arrange(power_plot, fdr_plot2,widths=c(1.3,1), ncol = 2)

# Save plots
ggsave("Survival.eps", width = 20, height = 8, plot = grid.arrange(power_plot, fdr_plot2,widths=c(1.45,1), ncol = 2), 
       device = "eps", dpi = 1200, units = "cm")

ggsave("Survival.jpg", width = 20, height = 8, plot = grid.arrange(power_plot, fdr_plot2,widths=c(1.45,1), ncol = 2), 
       device = "jpg", dpi = 300, units = "cm")

