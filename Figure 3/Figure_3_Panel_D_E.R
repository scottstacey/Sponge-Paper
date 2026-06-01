library(tidyverse)   # for data manipulation and ggplot2
library(zoo)         # for rollmean and rollapply functions
library(patchwork)   # for combining multiple plots
library(matrixStats)
library(ggplot2)
library(dplyr)
library(tidyr)

# Reading in data
pSS_02_005 = read.csv('Data/Figure_3_Panel_D_E/pSS-02-005.csv', header = TRUE)
pSS_02_006 = read.csv('Data/Figure_3_Panel_D_E/pSS-02-006.csv', header = TRUE)

# Converting time to hours
pSS_02_005$exp_time = pSS_02_005$exp_time/60/60
pSS_02_006$exp_time = pSS_02_006$exp_time/60/60

# Subtracting background fluorescence
## pSS-02-005
transition_row  = which(diff(pSS_02_005$Van) == 5.5)[1] + 1
last30          = pSS_02_005[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 = mean(last30$FP1_emit1, na.rm = TRUE)
mean_FP2_last30 = mean(last30$FP2_emit1, na.rm = TRUE)
pSS_02_005$FP1_emit1  = pSS_02_005$FP1_emit1 - mean_FP1_last30
pSS_02_005$FP2_emit1  = pSS_02_005$FP2_emit1 - mean_FP2_last30

## pSS-02-006
transition_row  = which(diff(pSS_02_006$Van) == 5.5)[1] + 1
last30          = pSS_02_006[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 = mean(last30$FP1_emit1, na.rm = TRUE)
mean_FP2_last30 = mean(last30$FP2_emit1, na.rm = TRUE)
pSS_02_006$FP1_emit1  = pSS_02_006$FP1_emit1 - mean_FP1_last30
pSS_02_006$FP2_emit1  = pSS_02_006$FP2_emit1 - mean_FP2_last30

# Ratio columns 
pSS_02_005$Ratio = pSS_02_005$FP2_emit1/pSS_02_005$FP1_emit1
pSS_02_006$Ratio = pSS_02_006$FP2_emit1/pSS_02_006$FP1_emit1

# Correcting mScarlet-I and Normalising GFP   
## pSS-02-005
transition_row1 <- which(diff(pSS_02_005$OC6) == 139)[1] + 1
last30 <- pSS_02_005[(transition_row1 - 30):(transition_row1 - 1), ]
mean_ratio = mean(last30$Ratio, na.rm = TRUE)
pSS_02_005$FP2_emit1 = pSS_02_005$FP2_emit1 - mean_ratio*pSS_02_005$FP1_emit1
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_005$FP1_emit1 = pSS_02_005$FP1_emit1 / mean_FP1_last30

## pSS-02-006
transition_row1 <- which(diff(pSS_02_006$OC6) == 139)[1] + 1
last30 <- pSS_02_006[(transition_row1 - 30):(transition_row1 - 1), ]
mean_ratio = mean(last30$Ratio, na.rm = TRUE)
pSS_02_006$FP2_emit1 = pSS_02_006$FP2_emit1 - mean_ratio*pSS_02_006$FP1_emit1
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_006$FP1_emit1 = pSS_02_006$FP1_emit1 / mean_FP1_last30



# Moving Averages 
## pSS-02-005
pSS_02_005$FP1_emit1_MA21      <- rollapply(pSS_02_005$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_005$FP2_emit1_MA21      <- rollapply(pSS_02_005$FP2_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_005$FP1_emit1_MA21[pSS_02_005$FP1_emit1_MA21 < 0] <- 0
pSS_02_005$FP2_emit1_MA21[pSS_02_005$FP2_emit1_MA21 < 0] <- 0

## pSS-02-006
pSS_02_006$FP1_emit1_MA21      <- rollapply(pSS_02_006$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_006$FP2_emit1_MA21      <- rollapply(pSS_02_006$FP2_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_006$FP1_emit1_MA21[pSS_02_006$FP1_emit1_MA21 < 0] <- 0
pSS_02_006$FP2_emit1_MA21[pSS_02_006$FP2_emit1_MA21 < 0] <- 0


## Plotting pSS-02-005

First_transition_row  = which(diff(pSS_02_005$Van) == 5.5)[1] + 1
First_transition_time = pSS_02_005$exp_time[First_transition_row]
Second_transition_row  = which(diff(pSS_02_005$OC6) == 139)[1] + 1
Second_transition_time = pSS_02_005$exp_time[Second_transition_row]
Third_transition_row  = which(diff(pSS_02_005$IPTG) == 1000)[1] + 1
Third_transition_time = pSS_02_005$exp_time[Third_transition_row]

factor <- 1.2 / 0.1 

p1 <- ggplot(pSS_02_005, aes(x = exp_time)) +
  geom_line(aes(y = FP1_emit1_MA21, color = "Green Fluorescence"), linewidth = 2) +
  geom_line(aes(y = FP2_emit1_MA21 * factor, color = "Red Fluorescence"), linewidth = 2) +
  
  geom_vline(xintercept = First_transition_time, linetype = "dotted", color = "blue", linewidth = 1) +
  geom_vline(xintercept = Second_transition_time, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = Third_transition_time, linetype = "dotted", color = "grey", linewidth = 1) +
  annotate("text", x = First_transition_time + 2.3 , y = 1.1, label = "[Van]\n5.5 µM", 
           color = "blue", size = 7, hjust = 0.5) +
  annotate("text", x = Second_transition_time + 2.35, y = 1.1, label = "[OC6]\n139 nM", 
           color = "black", size = 7, hjust = 0.5) +
  annotate("text", x = Third_transition_time + 2.75, y = 1.1, label = "[IPTG]\n1000 µM",
           color = "grey", size = 7, hjust = 0.5) + 
  
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20, hjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 16, hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.position = "inside",                        
    legend.position.inside = c(0.875, 0.875),              
    legend.background = element_rect(fill = "white", color = "black", linewidth = 1),  
    legend.key = element_rect(fill = "white", color = "white", linewidth = 1),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(colour = "grey90", linewidth = 0.5, linetype = "dashed")         
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(5, 50), breaks = seq(5, 50, 5)) +
  scale_y_continuous(
    name = "Relative GFPmut3 Fluorescence (AU)",
    expand = c(0, 0),
    limits = c(0, 1.2),
    breaks = seq(0, 1.2, 0.2),
    sec.axis = sec_axis(~ . / factor, name = "Red Fluorescence (AU)")
  ) +
  scale_color_manual(
    name = "Signal",
    values = c("Red Fluorescence" = "darkred", "Green Fluorescence" = "green")
  ) +
  labs(
    x = "Time (Hours)",
    title = "pSS-02-005 GFPmut3 Fluorescence Timecourse in Chi.Bio"
  )

print(p1)
svg("pSS-02-005.svg", width = 11.2559, height = 7.7362, bg = "white")
print(p1) 
dev.off()



## Plotting pSS-02-006
First_transition_row  = which(diff(pSS_02_006$Van) == 5.5)[1] + 1
First_transition_time = pSS_02_006$exp_time[First_transition_row]
Second_transition_row  = which(diff(pSS_02_006$OC6) == 139)[1] + 1
Second_transition_time = pSS_02_006$exp_time[Second_transition_row]
Third_transition_row  = which(diff(pSS_02_006$IPTG) == 1000)[1] + 1
Third_transition_time = pSS_02_006$exp_time[Third_transition_row]


factor <- 1.2 / 0.1  

p2 <- ggplot(pSS_02_006, aes(x = exp_time)) +
  geom_line(aes(y = FP1_emit1_MA21, color = "Green Fluorescence"), linewidth = 2) +
  geom_line(aes(y = FP2_emit1_MA21 * factor, color = "Red Fluorescence"), linewidth = 2) +
  geom_vline(xintercept = First_transition_time, linetype = "dotted", color = "blue", linewidth = 1) +
  geom_vline(xintercept = Second_transition_time, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = Third_transition_time, linetype = "dotted", color = "grey", linewidth = 1) +
  annotate("text", x = First_transition_time + 2.3 , y = 1.1, label = "[Van]\n5.5 µM", 
           color = "blue", size = 7, hjust = 0.5) +
  annotate("text", x = Second_transition_time + 2.35, y = 1.1, label = "[OC6]\n139 nM", 
           color = "black", size = 7, hjust = 0.5) +
  annotate("text", x = Third_transition_time + 2.75, y = 1.1, label = "[IPTG]\n1000 µM",
           color = "grey", size = 7, hjust = 0.5) + 
  
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20, hjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 16, hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.position = "inside",                        
    legend.position.inside = c(0.875, 0.875),              
    legend.background = element_rect(fill = "white", color = "black", linewidth = 1),  
    legend.key = element_rect(fill = "white", color = "white", linewidth = 1),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(colour = "grey90", linewidth = 0.5, linetype = "dashed")         
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(5, 50), breaks = seq(5, 50, 5)) +
  scale_y_continuous(
    name = "Relative GFPmut3 Fluorescence (AU)",
    expand = c(0, 0),
    limits = c(0, 1.2),
    breaks = seq(0, 1.2, 0.2),
    sec.axis = sec_axis(~ . / factor, name = "mScarlet-I Fluorescence (AU)")
  ) +
  scale_color_manual(
    name = "Signal",
    values = c("Red Fluorescence" = "darkred", "Green Fluorescence" = "green")
  ) +
  labs(
    x = "Time (Hours)",
    title = "pSS-02-006 GFPmut3 Fluorescence Timecourse in Chi.Bio"
  )

print(p2)
svg("pSS-02-006.svg", width = 11.2559, height = 7.7362, bg = "white")
print(p2) 
dev.off()


