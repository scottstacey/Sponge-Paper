library(ggplot2)
library(dplyr)
library(purrr)
library(zoo)

process_data <- function(file) {
  df <- read.csv(file)
  df$exp_time <- df$exp_time / 60 / 60
  
  van_row <- which(diff(df$Van) == 5.5)[1] + 1
  bg <- mean(df$FP1_emit1[(van_row - 30):(van_row - 1)], na.rm = TRUE)
  df$FP1_emit1 <- df$FP1_emit1 - bg
  
  oc6_row <- which(diff(df$OC6) == 50)[1] + 1
  norm_factor <- mean(df$FP1_emit1[(oc6_row - 30):(oc6_row - 1)], na.rm = TRUE)
  df$FP1_emit1 <- df$FP1_emit1 / norm_factor
  
  df$FP1_emit1_MA21 <- rollapply(
    df$FP1_emit1,
    width = 21,
    FUN = mean,
    align = "center",
    fill = NA
  )
  
  df$FP1_emit1_MA21[df$FP1_emit1_MA21 < 0] <- 0
  
  df
}

datasets <- list(
  CC8_M2 = process_data("Data/Figure_5_Panel_E/CC8_M2.csv"),
  CC8_M3 = process_data("Data/Figure_5_Panel_E/CC8_M3.csv"),
  CC9_M0 = process_data("Data/Figure_5_Panel_E/CC9_M0.csv"),
  CC9_M1 = process_data("Data/Figure_5_Panel_E/CC9_M1.csv"),
  CC9_M2 = process_data("Data/Figure_5_Panel_E/CC9_M2.csv"),
  CC9_M3 = process_data("Data/Figure_5_Panel_E/CC9_M3.csv")
)

combined_df <- map2_dfr(datasets, names(datasets), function(df, name) {
  df$condition <- name
  df$max_IPTG <- max(df$IPTG, na.rm = TRUE)
  df
})

ref_df <- datasets$CC8_M2

van_time  <- ref_df$exp_time[which(diff(ref_df$Van) == 5.5)[1] + 1]
oc6_time  <- ref_df$exp_time[which(diff(ref_df$OC6) == 50)[1] + 1]
iptg_time <- ref_df$exp_time[which(diff(ref_df$IPTG) == 160)[1] + 1]

p_all <- ggplot(combined_df, aes(x = exp_time, y = FP1_emit1_MA21, group = condition)) +
  geom_line(aes(color = factor(max_IPTG)), linewidth = 1.5) +
  
  geom_vline(xintercept = van_time, linetype = "dotted", color = "black", linewidth = 1) +
  geom_vline(xintercept = oc6_time, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = iptg_time, linetype = "dotdash", color = "black", linewidth = 1) +
  
  annotate("text", x = van_time + 3, y = 1.2, label = "[Van]\n5.5 µM",
           color = "black", size = 8, hjust = 0.5) +
  annotate("text", x = oc6_time + 2.5, y = 1.2, label = "[OC6]\n50 nM",
           color = "black", size = 8, hjust = 0.5) +
  annotate("text", x = iptg_time + 3, y = 1.2, label = "Sponge\nInduced",
           color = "black", size = 8, hjust = 0.5) +
  
  scale_color_manual(
    name = "Max IPTG (µM)",
    values = c(
      "0" = "black",
      "25" = "#006400",
      "100" = "#228B22",
      "160" = "#32CD32",
      "200" = "#7CFC00",
      "1000" = "#ADFF2F"
    )
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(5, 41),
    breaks = seq(5, 40, 5)
  ) +
  scale_y_continuous(
    name = "Relative GFPmut3 Fluorescence (AU)",
    expand = c(0, 0),
    limits = c(0, 1.4),
    breaks = seq(0, 1.4, 0.2)
  ) +
  labs(
    x = "Time (Hours)",
    title = "Chi.Bio characterisation of pSS-02-005 \nwith varying sponge induction"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20, hjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.text = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p_all)