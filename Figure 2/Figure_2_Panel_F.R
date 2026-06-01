library(tidyverse)   
library(zoo)         
library(patchwork)   
library(matrixStats)
library(ggplot2)
library(dplyr)
library(tidyr)

Data = read.csv('Data/Figure_2_Panel_F/pSS-02-24.csv', header = TRUE)

# Subtracting background fluorescence
transition_row <- which(diff(Data$Van) == 5.5)[1] + 1
last30 <- Data[(transition_row - 30):(transition_row - 1), ]
mean_FP2_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
Data$FP1_emit1 = Data$FP1_emit1 - mean_FP2_last30


# Normalise intitial output to 1 
transition_row1 <- which(diff(Data$OC6) == 139)[1] + 1
last30 <- Data[(transition_row1 - 30):(transition_row1 - 1), ]
mean_FP2_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
Data$FP1_emit1 = Data$FP1_emit1 / mean_FP2_last30


# Calculate moving averages: 
Data$FP2_emit1_MA21      <- rollapply(Data$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)

# Setting negative fluorescence values to zero 
Data$FP2_emit1_MA21[Data$FP2_emit1_MA21 < 0] <- 0

transition_time_1 <- Data$Exp_time[transition_row]
transition_time  = Data$Exp_time[transition_row1]

p1 <- ggplot(Data, aes(x = Exp_time, y = FP2_emit1_MA21)) +
  geom_line(linewidth = 3, color = "darkgreen") +
  geom_vline(xintercept = transition_time_1, linetype = "dotted", color = "blue", linewidth = 2) + 
  geom_vline(xintercept = transition_time, linetype = "dashed", color = "black", linewidth = 2) +
  annotate("text", x = transition_time_1 - 1, y = 1.5, label = "[Van]\n5.5 µM", 
           color = "blue", size = 7, hjust = 0.5) +
  annotate("text", x = transition_time - 1, y = 1.5, label = "[OC6]\n139 nM", 
           color = "black", size = 7, hjust = 0.5) +
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
  scale_x_continuous(expand = c(0, 0), limits = c(10, 25.2), breaks = seq(10, 25, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2), breaks = seq(0, 2, 0.5)) + 
  labs(
    x = "Time (Hours)",
    y = "Relative GFPmut3 Fluorescence (AU)",
    title = "pSS-02-024 GFPmut3 Fluorescence Timecourse in Chi.Bio"
  ) 



print(p1)
svg("pSS-02-024_GFPmut3.pdf", width = 7.5*1.5, height = 4.55*1.5, bg = "white")
print(p1)
dev.off()
