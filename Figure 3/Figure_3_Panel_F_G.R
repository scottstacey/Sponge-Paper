library(tidyverse)   
library(zoo)        
library(patchwork)   
library(matrixStats)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

pSS_02_044 = read.csv('Data/Figure_3_Panel_F_G/pSS-02-044.csv', header = TRUE)
pSS_02_045 = read.csv('Data/Figure_3_Panel_F_G/pSS-02-045.csv', header = TRUE)

pSS_02_044$exp_time = pSS_02_044$exp_time/60/60
pSS_02_045$exp_time = pSS_02_045$exp_time/60/60

# Remove background fluorescence: 
transition_row  = which(diff(pSS_02_044$Van) == 5.5)[1] + 1
last30          = pSS_02_044[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 = mean(last30$FP1_emit1, na.rm = TRUE)
mean_FP2_last30 = mean(last30$FP2_emit1, na.rm = TRUE)
pSS_02_044$FP1_emit1  = pSS_02_044$FP1_emit1 - mean_FP1_last30
pSS_02_044$FP2_emit1  = pSS_02_044$FP2_emit1 - mean_FP2_last30

transition_row  = which(diff(pSS_02_045$Van) == 5.5)[1] + 1
last30          = pSS_02_045[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 = mean(last30$FP1_emit1, na.rm = TRUE)
mean_FP2_last30 = mean(last30$FP2_emit1, na.rm = TRUE)
pSS_02_045$FP1_emit1  = pSS_02_045$FP1_emit1 - mean_FP1_last30
pSS_02_045$FP2_emit1  = pSS_02_045$FP2_emit1 - mean_FP2_last30

# Make Ratios
pSS_02_044$Ratio = pSS_02_044$FP2_emit1/pSS_02_044$FP1_emit1
pSS_02_045$Ratio = pSS_02_045$FP2_emit1/pSS_02_045$FP1_emit1


# Normalising GFPmut3 and Correcting mScarlet-I 
transition_row1 <- which(diff(pSS_02_044$OC6) == 50)[1] + 1
pSS_02_044$growth_rate_ma      <- rollapply(pSS_02_044$growth_rate, width = 31, FUN = mean, align = "center", fill = NA)
pSS_02_044$GrowthNormalisedGFP = pSS_02_044$FP1_emit1 * pSS_02_044$growth_rate_ma
last30 <- pSS_02_044[(transition_row1 - 30):(transition_row1 - 1), ]
mean_ratio = mean(last30$Ratio, na.rm = TRUE)
pSS_02_044$FP2_emit1 = pSS_02_044$FP2_emit1 - mean_ratio*pSS_02_044$FP1_emit1
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
mean_FP1_GrowthNornalised = mean(last30$GrowthNormalisedGFP, na.rm = TRUE)
pSS_02_044$FP1_emit1 = pSS_02_044$FP1_emit1 / mean_FP1_last30
pSS_02_044$GrowthNormalisedGFP = pSS_02_044$GrowthNormalisedGFP / mean_FP1_GrowthNornalised

transition_row1 <- which(diff(pSS_02_045$OC6) == 50)[1] + 1
pSS_02_045$growth_rate_ma      <- rollapply(pSS_02_045$growth_rate, width = 31, FUN = mean, align = "center", fill = NA)
pSS_02_045$GrowthNormalisedGFP = pSS_02_045$FP1_emit1 * pSS_02_045$growth_rate_ma
last30 <- pSS_02_045[(transition_row1 - 30):(transition_row1 - 1), ]
mean_ratio = mean(last30$Ratio, na.rm = TRUE)
pSS_02_045$FP2_emit1 = pSS_02_045$FP2_emit1 - mean_ratio*pSS_02_045$FP1_emit1
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
mean_FP1_GrowthNornalised = mean(last30$GrowthNormalisedGFP, na.rm = TRUE)
pSS_02_045$FP1_emit1 = pSS_02_045$FP1_emit1 / mean_FP1_last30
pSS_02_045$GrowthNormalisedGFP = pSS_02_045$GrowthNormalisedGFP / mean_FP1_GrowthNornalised


# Moving Averages 
pSS_02_044$FP1_emit1_MA21             <- rollapply(pSS_02_044$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_044$FP2_emit1_MA21             <- rollapply(pSS_02_044$FP2_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_044$GrowthNormalisedGFP_MA     <- rollapply(pSS_02_044$GrowthNormalisedGFP, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_044$FP1_emit1_MA21[pSS_02_044$FP1_emit1_MA21 < 0] <- 0
pSS_02_044$FP2_emit1_MA21[pSS_02_044$FP2_emit1_MA21 < 0] <- 0
pSS_02_044$GrowthNormalisedGFP_MA[pSS_02_044$GrowthNormalisedGFP_MA < 0] <- 0

pSS_02_045$FP1_emit1_MA21             <- rollapply(pSS_02_045$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_045$FP2_emit1_MA21             <- rollapply(pSS_02_045$FP2_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_045$GrowthNormalisedGFP_MA     <- rollapply(pSS_02_045$GrowthNormalisedGFP, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_045$FP1_emit1_MA21[pSS_02_045$FP1_emit1_MA21 < 0] <- 0
pSS_02_045$FP2_emit1_MA21[pSS_02_045$FP2_emit1_MA21 < 0] <- 0
pSS_02_045$GrowthNormalisedGFP_MA[pSS_02_045$GrowthNormalisedGFP_MA < 0] <- 0


## Plotting pSS-02-044
First_transition_row  = which(diff(pSS_02_044$Van) == 5.5)[1] + 1
First_transition_time = pSS_02_044$exp_time[First_transition_row]
Second_transition_row  = which(diff(pSS_02_044$OC6) == 50)[1] + 1
Second_transition_time = pSS_02_044$exp_time[Second_transition_row]
Third_transition_row  = which(diff(pSS_02_044$IPTG) == 1000)[1] + 1
Third_transition_time = pSS_02_044$exp_time[Third_transition_row]

factor <- 1.4 / 0.12

p1 <- ggplot(pSS_02_044, aes(x = exp_time)) +
  geom_line(aes(y = FP1_emit1_MA21, color = "GFPmut3"), linewidth = 2) +
  geom_line(aes(y = FP2_emit1_MA21 * factor, color = "mScarlet-I"), linewidth = 2) +
  geom_vline(xintercept = First_transition_time, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = Second_transition_time, linetype = "dotted", color = "black", linewidth = 1) +
  geom_vline(xintercept = Third_transition_time, linetype = "twodash", color = "black", linewidth = 1) +
  annotate("text", x = First_transition_time + 0.5, y = 1.3, label = "[Van]\n5.5 µM", 
           color = "black", size = 7, hjust = 0) +
  annotate("text", x = Second_transition_time + 0.5, y = 1.3, label = "[OC6]\n50 nM", 
           color = "black", size = 7, hjust = 0) +
  annotate("text", x = Third_transition_time + 0.5, y = 1.3, label = "[IPTG]\n1000 µM",
           color = "black", size = 7, hjust = 0) + 
  
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20, hjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 16, hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.position = "inside",                        
    legend.position.inside = c(0.1, 0.175),             
    legend.background = element_rect(fill = "white", color = "black", linewidth = 1),  
    legend.key = element_rect(fill = "white", color = "white", linewidth = 1),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(colour = "grey90", linewidth = 0.5, linetype = "dashed")         
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(5, 33), breaks = seq(5, 33, 4)) +
  scale_y_continuous(
    name = "Relative GFPmut3 Fluorescence (AU)",
    expand = c(0, 0),
    limits = c(0, 1.4),
    breaks = seq(0, 1.4, 0.2),
    sec.axis = sec_axis(~ . / factor, name = "mScarlet-I Fluorescence (AU)")
  ) +
  scale_color_manual(
    name = "Signal",
    values = c("mScarlet-I" = "darkred", "GFPmut3" = "green")
  ) +
  labs(
    x = "Time (Hours)",
    title = "pSS-02-044 Fluorescence Timecourse in Chi.Bio"
  )

print(p1)
svg("pSS-02-044.svg", width = 11.2559, height = 7.7362, bg = "white")
print(p1) 
dev.off()

## Plotting pSS-02-045
First_transition_row  = which(diff(pSS_02_045$Van) == 5.5)[1] + 1
First_transition_time = pSS_02_045$exp_time[First_transition_row]
Second_transition_row  = which(diff(pSS_02_045$OC6) == 50)[1] + 1
Second_transition_time = pSS_02_045$exp_time[Second_transition_row]
Third_transition_row  = which(diff(pSS_02_045$IPTG) == 1000)[1] + 1
Third_transition_time = pSS_02_045$exp_time[Third_transition_row]


factor <- 1.4 / 0.12  

p2 <- ggplot(pSS_02_045, aes(x = exp_time)) +
  geom_line(aes(y = FP1_emit1_MA21, color = "GFPmut3"), linewidth = 2) +
  geom_line(aes(y = FP2_emit1_MA21 * factor, color = "mScarlet-I"), linewidth = 2) +
  geom_vline(xintercept = First_transition_time, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = Second_transition_time, linetype = "dotted", color = "black", linewidth = 1) +
  geom_vline(xintercept = Third_transition_time, linetype = "twodash", color = "black", linewidth = 1) +
  annotate("text", x = First_transition_time + 0.5, y = 1.3, label = "[Van]\n5.5 µM", 
           color = "black", size = 7, hjust = 0) +
  annotate("text", x = Second_transition_time + 0.5, y = 1.3, label = "[OC6]\n50 nM", 
           color = "black", size = 7, hjust = 0) +
  annotate("text", x = Third_transition_time + 0.5, y = 1.3, label = "[IPTG]\n1000 µM",
           color = "black", size = 7, hjust = 0) + 
  
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20, hjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 16, hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.position = "inside",                        
    legend.position.inside = c(0.1, 0.175),              
    legend.background = element_rect(fill = "white", color = "black", linewidth = 1),  
    legend.key = element_rect(fill = "white", color = "white", linewidth = 1),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(colour = "grey90", linewidth = 0.5, linetype = "dashed")        
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(5, 33), breaks = seq(5, 33, 4)) +
  scale_y_continuous(
    name = "Relative GFPmut3 Fluorescence (AU)",
    expand = c(0, 0),
    limits = c(0, 1.4),
    breaks = seq(0, 1.4, 0.2),
    sec.axis = sec_axis(~ . / factor, name = "mScarlet-I Fluorescence (AU)")
  ) +
  scale_color_manual(
    name = "Signal",
    values = c("mScarlet-I" = "darkred", "GFPmut3" = "green")
  ) +
  labs(
    x = "Time (Hours)",
    title = "pSS-02-045 Fluorescence Timecourse in Chi.Bio"
  )

print(p2)
svg("pSS-02-045.svg", width = 11.2559, height = 7.7362, bg = "white")
print(p2)  
dev.off()