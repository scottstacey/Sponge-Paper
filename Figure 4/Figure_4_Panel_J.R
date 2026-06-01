library(tidyverse)   
library(zoo)         
library(patchwork)   
library(matrixStats)
library(ggplot2)
library(dplyr)
library(tidyr)


# Reading in data
Data = read.csv('Data/Figure_4_Panel_J/Data.csv', header = TRUE)


# Remove background fluorescence 
transition_row  = which(diff(Data$Van) == 5.5)[1] + 1
last30          = Data[(transition_row - 30):(transition_row - 1), ]
Data$Sponge.Circuit...mutagenised.mScarlet.I = Data$Sponge.Circuit...mutagenised.mScarlet.I - mean(last30$Sponge.Circuit...mutagenised.mScarlet.I)
Data$Sponge.Circuit...mScarlet.I = Data$Sponge.Circuit...mScarlet.I - mean(last30$Sponge.Circuit...mScarlet.I)
Data$Negative.Control.Sponge.Circuit...mScarlet = Data$Negative.Control.Sponge.Circuit...mScarlet - mean(last30$Negative.Control.Sponge.Circuit...mScarlet)
Data$Sponge.Circuit = Data$Sponge.Circuit - mean(last30$Sponge.Circuit)


# Normalising GFP 
transition_row1 <- which(diff(Data$OC6) == 50)[1] + 1
last30 <- Data[(transition_row1 - 30):(transition_row1 - 1), ]
Data$Sponge.Circuit...mutagenised.mScarlet.I = Data$Sponge.Circuit...mutagenised.mScarlet.I/mean(last30$Sponge.Circuit...mutagenised.mScarlet.I)
Data$Sponge.Circuit...mScarlet.I = Data$Sponge.Circuit...mScarlet.I/mean(last30$Sponge.Circuit...mScarlet.I)
Data$Sponge.Circuit = Data$Sponge.Circuit/mean(last30$Sponge.Circuit)
Data$Negative.Control.Sponge.Circuit...mScarlet = Data$Negative.Control.Sponge.Circuit...mScarlet/mean(last30$Negative.Control.Sponge.Circuit...mScarlet)

# Moving Averages
Data$Sponge.Circuit...mutagenised.mScarlet.I_MA      <- rollapply(Data$Sponge.Circuit...mutagenised.mScarlet.I, width = 21, FUN = mean, align = "center", fill = NA)
Data$Sponge.Circuit...mScarlet.I_MA      <- rollapply(Data$Sponge.Circuit...mScarlet.I, width = 21, FUN = mean, align = "center", fill = NA)
Data$Sponge.Circuit_MA      <- rollapply(Data$Sponge.Circuit, width = 21, FUN = mean, align = "center", fill = NA)
Data$Negative.Control.Sponge.Circuit...mScarlet_MA      <- rollapply(Data$Negative.Control.Sponge.Circuit...mScarlet, width = 21, FUN = mean, align = "center", fill = NA)

Data$Negative.Control.Sponge.Circuit...mScarlet_MA[Data$Negative.Control.Sponge.Circuit...mScarlet_MA < 0] <- 0
Data$Sponge.Circuit_MA[Data$Sponge.Circuit_MA < 0] <- 0
Data$Sponge.Circuit...mScarlet.I_MA[Data$Sponge.Circuit...mScarlet.I_MA < 0] <- 0
Data$Sponge.Circuit...mutagenised.mScarlet.I_MA[Data$Sponge.Circuit...mutagenised.mScarlet.I_MA < 0] <- 0


First_transition_time = Data$Time[transition_row]
Second_transition_time = Data$Time[transition_row1]
transition_row2 = which(diff(Data$IPTG) == 1000)[1] + 1
Third_transition_time = Data$Time[transition_row2]

greens <- c(
  "pSS-02-027" = "#005a32",
  "pSS-02-009" = "#1b9e77",
  "pSS-02-006" = "#66a61e"
)

p1 = ggplot(Data, aes(x = Time)) + 
  geom_line(aes(y = Sponge.Circuit...mutagenised.mScarlet.I_MA, 
                color = "pSS-02-009"), linewidth = 2) +
  geom_line(aes(y = Sponge.Circuit...mScarlet.I_MA, 
                color = "pSS-02-006"), linewidth = 2) + 
  geom_line(aes(y = Negative.Control.Sponge.Circuit...mScarlet_MA, 
                color = "pSS-02-027"), linewidth = 2) +
  theme_minimal()+
  
  # Add vertical lines and annotations as before
  geom_vline(xintercept = First_transition_time, linetype = "dotted", color = "black", linewidth = 1) +
  geom_vline(xintercept = Second_transition_time, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = Third_transition_time, linetype = "dotted", color = "black", linewidth = 1) +
  annotate("text", x = First_transition_time +1.3 , y = 1.5, label = "[Van]\n5.5 µM", 
           color = "black", size = 6, hjust = 0.5) +
  annotate("text", x = Second_transition_time + 1.3, y = 1.5, label = "[OC6]\n50 nM", 
           color = "black", size = 6, hjust = 0.5) +
  annotate("text", x = Third_transition_time + 1.3, y = 1.5, label = "[IPTG]\n1000 µM",
           color = "black", size = 6, hjust = 0.5) +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20, hjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 16, hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.position = "inside",
    legend.position.inside = c(0.82, 0.82),
    legend.background = element_rect(fill = "transparent", color = "NA", linewidth = 0),
    legend.key = element_rect(fill = "transparent", color = "NA", linewidth = 0),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(colour = "grey90", linewidth = 0.5, linetype = "dashed")
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(12, 36.5), breaks = seq(12, 36, 6)) +
  scale_y_continuous(
    name = "Relative GFPmut3 Fluorescence (AU)",
    expand = c(0, 0),
    limits = c(0, 1.62),
    breaks = seq(0, 1.6, 0.2)
  ) +
  scale_color_manual(values = greens, name = "") +
  labs(x = "Time (Hours)", title = "GFPmut3 Chi.Bio timecourse")
  

print(p1)

