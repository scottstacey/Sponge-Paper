library(tidyverse)  
library(zoo)         
library(patchwork)   
library(matrixStats)
library(ggplot2)
library(dplyr)
library(tidyr)


pSS_02_092_Raw_Data = read.csv('Data/Figure_2_Panel_G/pSS_02_092_ChiBio.csv', header = TRUE)
keep_cols <- c("exp_time", "od_measured", "OC6..nM.", "IPTG..uM.", "FP2_emit1", "growth_rate")

Data = pSS_02_092_Raw_Data[, keep_cols]


# Subtracting background fluorescence
transition_row <- which(diff(Data$IPTG..uM.) == 30)[1] + 1
last30 <- Data[(transition_row - 30):(transition_row - 1), ]
mean_FP2_last30 <- mean(last30$FP2_emit1, na.rm = TRUE)


Data$FP2_emit1 = Data$FP2_emit1 - mean_FP2_last30

# Normalise intitial output to 1 
transition_row <- which(diff(Data$OC6..nM.) == 60)[1] + 1
last30 <- Data[(transition_row - 30):(transition_row - 1), ]
mean_FP2_last30 <- mean(last30$FP2_emit1, na.rm = TRUE)
Data$FP2_emit1 = Data$FP2_emit1 / mean_FP2_last30


# Calculate moving averages: 
Data$FP2_emit1_MA21      <- rollapply(Data$FP2_emit1, width = 21, FUN = mean, align = "center", fill = NA)

# Setting negative fluorescence values to zero 
Data$FP2_emit1_MA21[Data$FP2_emit1_MA21 < 0] <- 0

transition_row_1 <- which(diff(Data$IPTG..uM.) == 30)[1] + 1
transition_time_1 <- Data$exp_time[transition_row_1]

transition_time  = Data$exp_time[transition_row]

p1 <- ggplot(Data, aes(x = exp_time, y = FP2_emit1_MA21)) +
  geom_line(linewidth = 2, color = "darkred") +
  geom_vline(xintercept = transition_time_1, linetype = "dotted", color = "blue", linewidth = 2) + 
  geom_vline(xintercept = transition_time, linetype = "dashed", color = "black", linewidth = 2) +
  annotate("text", x = transition_time_1 - 4.5, y = 0.75, label = "[IPTG]\n30 µM", 
           color = "blue", size = 7, hjust = 0) +
  annotate("text", x = transition_time + 1.5, y = 0.75, label = "[OC6]\n60 nM", 
           color = "black", size = 7, hjust = 0) +
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
  scale_x_continuous(expand = c(0, 0), limits = c(0, 40 + 0.5), breaks = seq(0, 40, 10)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2)) + 
  labs(
    x = "Time (Hours)",
    y = "Relative mScarlet-I Fluorescence (AU)",
    title = "pSS-02-092 mScarlet-I Chi.Bio Fluorescence Timecourse"
  ) 


print(p1)
svg("pSS-02-092_mScarlet.svg", width = 7.5*1.5, height = 4.55*1.5, bg = "white")
print(p1)
dev.off()

