library(tidyverse)   
library(zoo)         
library(patchwork)   
library(matrixStats)
library(ggplot2)
library(dplyr)
library(tidyr)

pSS_02_001_Raw_Data = read.csv('Data/Figure_2_Panel_E/pSS-02-001_Data.csv', header = TRUE)
keep_cols <- c("exp_time", "od_measured", "LED_395nm_setpoint", "FP1_emit1", "growth_rate")
pSS_02_001_Subset <- pSS_02_001_Raw_Data[, keep_cols]

#Setting up inducer concentration columns
pSS_02_001_Subset$`OC6 Concentration nM` <- ifelse(pSS_02_001_Subset$LED_395nm_setpoint == 0.111, 139, 0)
pSS_02_001_Subset$`Van Concentration uM` <- ifelse(pSS_02_001_Subset$LED_395nm_setpoint != 0.1, 5.5, 0)


#Convert time to hours 
pSS_02_001_Subset$exp_time = pSS_02_001_Subset$exp_time/60/60


# Subtracting background fluorescence
transition_row <- which(diff(pSS_02_001_Subset$`Van Concentration uM`) == 5.5)[1] + 1
last30 <- pSS_02_001_Subset[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_001_Subset$FP1_emit1 = pSS_02_001_Subset$FP1_emit1 - mean_FP1_last30

# Find average GFP at steady state: 
transition_row <- which(diff(pSS_02_001_Subset$`OC6 Concentration nM`) == 139)[1] + 1
last30 <- pSS_02_001_Subset[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_001_Subset$FP1_emit1 = pSS_02_001_Subset$FP1_emit1 / mean_FP1_last30

# Now calculate moving averages for relevant variables 
pSS_02_001_Subset$od_measured_MA21    <- rollapply(pSS_02_001_Subset$od_measured, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_001_Subset$FP1_emit1_MA21      <- rollapply(pSS_02_001_Subset$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_001_Subset$growth_rate_MA21    <- rollapply(pSS_02_001_Subset$growth_rate, width = 21, FUN = mean, align = "center", fill = NA)

# Setting negative fluorescence values to zero 
pSS_02_001_Subset$FP1_emit1_MA21[pSS_02_001_Subset$FP1_emit1_MA21 < 0] <- 0

# Plotting 
transition_row <- which(diff(pSS_02_001_Subset$`OC6 Concentration nM`) == 139)[1] + 1
transition_time <- pSS_02_001_Subset$exp_time[transition_row]

p1 <- ggplot(pSS_02_001_Subset, aes(x = exp_time, y = FP1_emit1_MA21)) +
  geom_line(linewidth = 2, color = "green") +
  geom_vline(xintercept = transition_time_1, linetype = "dotted", color = "black", linewidth = 2) + 
  geom_vline(xintercept = transition_time, linetype = "dashed", color = "darkgreen", linewidth = 2) +
  annotate("text", x = 20 - 8.5, y = 0.75, label = "[Van]\n5.5 µM", 
           color = "black", size = 7, hjust = 0) +
  annotate("text", x = transition_time + 1.5, y = 0.75, label = "[OC6]\n139 nM", 
           color = "darkgreen", size = 7, hjust = 0) + 
  theme_minimal() + 
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20, hjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(colour = "grey90", linewidth = 0.5, linetype = "dashed")
  ) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 50)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2)) + 
  labs(
    x = "Time (Hours)",
    y = "Relative GFPmut3 Fluorescence (AU)",
    title = "pSS-02-001 GFPmut3 Chi.Bio Fluorescence Timecourse"
  ) 

print(p1)
svg("pSS-02-001_GFPmut3.svg", width = 9.49, height = 6.53, bg = "white")
print(p1)
dev.off()