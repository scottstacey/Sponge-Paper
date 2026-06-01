library(tidyverse)   
library(zoo)         
library(patchwork)   
library(matrixStats)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)


pSS_02_005 = read.csv('Data/Figure_5_Panel_G/pSS-02-005.csv', header = TRUE)
pSS_02_103 = read.csv('Data/Figure_5_Panel_G/pSS-02-103.csv', header = TRUE)
pSS_02_104 = read.csv('Data/Figure_5_Panel_G/pSS-02-104.csv', header = TRUE)
pSS_02_105 = read.csv('Data/Figure_5_Panel_G/pSS-02-105.csv', header = TRUE)
pSS_02_106 = read.csv('Data/Figure_5_Panel_G/pSS-02-106.csv', header = TRUE)
pSS_02_005$exp_time = pSS_02_005$exp_time/60/60
pSS_02_103$exp_time = pSS_02_103$exp_time/60/60
pSS_02_104$exp_time = pSS_02_104$exp_time/60/60
pSS_02_105$exp_time = pSS_02_105$exp_time/60/60
pSS_02_106$exp_time = pSS_02_106$exp_time/60/60
pSS_02_005$Circuit = "pSS-02-005"
pSS_02_103$Circuit = "pSS-02-103"
pSS_02_104$Circuit = "pSS-02-104"
pSS_02_105$Circuit = "pSS-02-105"
pSS_02_106$Circuit = "pSS-02-106"


# Remove background fluorescence: 
transition_row  = which(diff(pSS_02_005$Van) == 5.5)[1] + 1
last30          = pSS_02_005[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 = mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_005$FP1_emit1  = pSS_02_005$FP1_emit1 - mean_FP1_last30

transition_row  = which(diff(pSS_02_103$Van) == 5.5)[1] + 1
last30          = pSS_02_103[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 = mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_103$FP1_emit1  = pSS_02_103$FP1_emit1 - mean_FP1_last30

transition_row  = which(diff(pSS_02_104$Van) == 5.5)[1] + 1
last30          = pSS_02_104[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 = mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_104$FP1_emit1  = pSS_02_104$FP1_emit1 - mean_FP1_last30

transition_row  = which(diff(pSS_02_105$Van) == 5.5)[1] + 1
last30          = pSS_02_105[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 = mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_105$FP1_emit1  = pSS_02_105$FP1_emit1 - mean_FP1_last30

transition_row  = which(diff(pSS_02_106$Van) == 5.5)[1] + 1
last30          = pSS_02_106[(transition_row - 30):(transition_row - 1), ]
mean_FP1_last30 = mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_106$FP1_emit1  = pSS_02_106$FP1_emit1 - mean_FP1_last30

# Normalising GFPmut3
transition_row1 <- which(diff(pSS_02_005$OC6) == 139)[1] + 1
last30 <- pSS_02_005[(transition_row1 - 30):(transition_row1 - 1), ]
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_005$FP1_emit1 = pSS_02_005$FP1_emit1 / mean_FP1_last30

transition_row1 <- which(diff(pSS_02_103$OC6) == 139)[1] + 1
last30 <- pSS_02_103[(transition_row1 - 30):(transition_row1 - 1), ]
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_103$FP1_emit1 = pSS_02_103$FP1_emit1 / mean_FP1_last30

transition_row1 <- which(diff(pSS_02_104$OC6) == 139)[1] + 1
last30 <- pSS_02_104[(transition_row1 - 30):(transition_row1 - 1), ]
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_104$FP1_emit1 = pSS_02_104$FP1_emit1 / mean_FP1_last30

transition_row1 <- which(diff(pSS_02_105$OC6) == 139)[1] + 1
last30 <- pSS_02_105[(transition_row1 - 30):(transition_row1 - 1), ]
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_105$FP1_emit1 = pSS_02_105$FP1_emit1 / mean_FP1_last30

transition_row1 <- which(diff(pSS_02_106$OC6) == 139)[1] + 1
last30 <- pSS_02_106[(transition_row1 - 30):(transition_row1 - 1), ]
mean_FP1_last30 <- mean(last30$FP1_emit1, na.rm = TRUE)
pSS_02_106$FP1_emit1 = pSS_02_106$FP1_emit1 / mean_FP1_last30

# Moving Averages 
pSS_02_005$FP1_emit1_MA21      <- rollapply(pSS_02_005$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_005$growth_rate      <- rollapply(pSS_02_005$growth_rate, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_005$FP1_emit1_MA21[pSS_02_005$FP1_emit1_MA21 < 0] <- 0

pSS_02_103$FP1_emit1_MA21      <- rollapply(pSS_02_103$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_103$growth_rate      <- rollapply(pSS_02_103$growth_rate, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_103$FP1_emit1_MA21[pSS_02_103$FP1_emit1_MA21 < 0] <- 0

pSS_02_104$FP1_emit1_MA21      <- rollapply(pSS_02_104$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_104$growth_rate      <- rollapply(pSS_02_104$growth_rate, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_104$FP1_emit1_MA21[pSS_02_104$FP1_emit1_MA21 < 0] <- 0

pSS_02_105$FP1_emit1_MA21      <- rollapply(pSS_02_105$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_105$growth_rate      <- rollapply(pSS_02_105$growth_rate, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_105$FP1_emit1_MA21[pSS_02_105$FP1_emit1_MA21 < 0] <- 0

pSS_02_106$FP1_emit1_MA21      <- rollapply(pSS_02_106$FP1_emit1, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_106$growth_rate      <- rollapply(pSS_02_106$growth_rate, width = 21, FUN = mean, align = "center", fill = NA)
pSS_02_106$FP1_emit1_MA21[pSS_02_106$FP1_emit1_MA21 < 0] <- 0


TotalData = rbind(pSS_02_005, pSS_02_103, pSS_02_104, pSS_02_105, pSS_02_106)

library(RColorBrewer)
color_palette <- brewer.pal(min(9, length(unique(TotalData$Circuit))), "Dark2") 

First_transition_row  = which(diff(pSS_02_005$Van) == 5.5)[1] + 1
First_transition_time = pSS_02_005$exp_time[First_transition_row]
Second_transition_row  = which(diff(pSS_02_005$OC6) == 139)[1] + 1
Second_transition_time = pSS_02_005$exp_time[Second_transition_row]
Fourth_transition_row  = which(diff(pSS_02_005$IPTG) == 1000)[1] + 1
Fourth_transition_time = pSS_02_005$exp_time[Fourth_transition_row]



p1 <- ggplot(TotalData, aes(x = exp_time, y = FP1_emit1_MA21, color = Circuit, group = Circuit)) +
  geom_line(linewidth =2) +
  scale_color_manual(values = color_palette) +
  geom_vline(xintercept = First_transition_time, linetype = "dashed", color = "black", size = 1) + 
  geom_vline(xintercept = Second_transition_time, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = Fourth_transition_time, linetype = "twodash", color = "black", size = 1) +
  annotate("text", x = First_transition_time + 2, y = 1.1, label = "[Van]\n5.5 µM", 
           color = "black", size = 7, hjust = 0.5) +
  annotate("text", x = Second_transition_time + 2, y = 1.1, label = "[OC6]\n139 nM", 
           color = "black", size = 7, hjust = 0.5) +
  annotate("text", x = Fourth_transition_time + 2, y = 1.1, label = "[IPTG]\n1000 µM", 
           color = "black", size = 7, hjust = 0.5) +
  labs(
    x = "Time (Hr)",
    y = "Relative GFPmut3 Fluorescence (AU)",
    title = "GFPmut3 Fluorescence Timecourses in Chi.Bio Engineered Sponges"
  ) + 
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
    panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed")
  ) + 
  scale_x_continuous(expand = c(0, 0), limits = c(5, 35), breaks = seq(5, 35, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2)) + 
  labs(
    x = "Time (Hr)",
    y = "Relative GFPmut3 Fluorescence (AU)",
    title = "Sponge Sequence Tuning"
  ) + 
  # scale_color_manual(
  #   values = scales::gradient_n_pal(c("lightgreen", "green", "darkgreen"))(seq(0, 1, length.out = length(unique(TotalData$Circuit))))
  # ) +
  guides(color = guide_legend(title = "Circuit", labels = as.character(1:length(unique(TotalData$Circuit)))))

print(p1)
pdf("Figure_3.49.pdf", width = 11.2559*1.5, height = 7.7362, bg = "white")
print(p1)  # Replace with your ggplot variable
dev.off()
