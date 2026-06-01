library(ggplot2)
library(dplyr)
library(zoo)

process_GFP <- function(file, circuit_name) {
  
  df <- read.csv(file, header = TRUE)
  df$exp_time <- df$exp_time / 60 / 60
  
  van_row <- which(diff(df$Van) == 5.5)[1] + 1
  oc6_row <- which(diff(df$OC6) == 139)[1] + 1
  iptg_row <- which(diff(df$IPTG) == 1000)[1] + 1
  
  bg_rows <- (van_row - 30):(van_row - 1)
  df$FP1_emit1 <- df$FP1_emit1 - mean(df$FP1_emit1[bg_rows], na.rm = TRUE)
  
  gfp_norm_rows <- (oc6_row - 30):(oc6_row - 1)
  gfp_norm <- mean(df$FP1_emit1[gfp_norm_rows], na.rm = TRUE)
  df$FP1_emit1 <- df$FP1_emit1 / gfp_norm
  
  df$GFP_MA21 <- rollapply(
    df$FP1_emit1,
    width = 21,
    FUN = mean,
    align = "center",
    fill = NA
  )
  
  df$GFP_MA21[df$GFP_MA21 < 0] <- 0
  
  df$circuit <- circuit_name
  df$van_row <- van_row
  df$oc6_row <- oc6_row
  df$iptg_row <- iptg_row
  
  df
}

get_sRNA_pre_spRNA_mean <- function(df) {
  
  iptg_row <- unique(df$iptg_row)
  rows <- (iptg_row - 15):(iptg_row - 1)
  
  mean(df$GFP_MA21[rows], na.rm = TRUE)
}

pSS_02_001 <- process_GFP("Data/Figure_4_Panel_I/pSS-02-001.csv", "pSS-02-001")
pSS_02_005 <- process_GFP("Data/Figure_4_Panel_I/pSS-02-005.csv", "pSS-02-005")
pSS_02_006 <- process_GFP("Data/Figure_4_Panel_I/pSS-02-006.csv", "pSS-02-006")

sRNA_levels <- c(
  "pSS-02-001" = get_sRNA_pre_spRNA_mean(pSS_02_001),
  "pSS-02-005" = get_sRNA_pre_spRNA_mean(pSS_02_005),
  "pSS-02-006" = get_sRNA_pre_spRNA_mean(pSS_02_006)
)

target_level <- max(sRNA_levels, na.rm = TRUE)
scale_factors <- target_level / sRNA_levels

pSS_02_001$GFP_MA21 <- pSS_02_001$GFP_MA21 * scale_factors["pSS-02-001"]
pSS_02_005$GFP_MA21 <- pSS_02_005$GFP_MA21 * scale_factors["pSS-02-005"]
pSS_02_006$GFP_MA21 <- pSS_02_006$GFP_MA21 * scale_factors["pSS-02-006"]

plot_data <- bind_rows(
  pSS_02_001,
  pSS_02_005,
  pSS_02_006
)

van_time <- pSS_02_001$exp_time[unique(pSS_02_001$van_row)]
oc6_time <- pSS_02_001$exp_time[unique(pSS_02_001$oc6_row)]
iptg_time <- pSS_02_001$exp_time[unique(pSS_02_001$iptg_row)]

p <- ggplot(plot_data, aes(x = exp_time, y = GFP_MA21, color = circuit)) +
  geom_line(linewidth = 1.5) +
  
  geom_vline(xintercept = van_time, linetype = "dotted", color = "black", linewidth = 1) +
  geom_vline(xintercept = oc6_time, linetype = "dashed", color = "black", linewidth = 1) +
  geom_vline(xintercept = iptg_time, linetype = "dotdash", color = "black", linewidth = 1) +
  
  annotate(
    "text",
    x = van_time + 2.3,
    y = 1.48,
    label = "[Van]\n5.5 µM",
    color = "black",
    size = 7,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = oc6_time + 2.3,
    y = 1.48,
    label = "[OC6]\n139 nM",
    color = "black",
    size = 7,
    hjust = 0.5
  ) +
  annotate(
    "text",
    x = iptg_time + 2.7,
    y = 1.48,
    label = "[IPTG]\n1000 µM",
    color = "black",
    size = 7,
    hjust = 0.5
  ) +
  
  scale_color_manual(
    name = "Circuit",
    values = c(
      "pSS-02-001" = "darkgreen",
      "pSS-02-005" = "limegreen",
      "pSS-02-006" = "chartreuse"
    )
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(5, 45),
    breaks = seq(5, 45, 5)
  ) +
  scale_y_continuous(
    name = "Relative GFPmut3 Fluorescence (AU)",
    expand = c(0, 0),
    limits = c(0, 1.6),
    breaks = seq(0, 2, 0.2)
  ) +
  labs(
    x = "Time (Hours)",
    title = "GFPmut3 Fluorescence Timecourse in Chi.Bio"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 20, hjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 16, hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.position = "inside",
    legend.position.inside = c(0.825, 0.825),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 1),
    legend.key = element_rect(fill = "white", color = "white", linewidth = 1),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.text = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)