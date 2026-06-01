library(matrixStats)
library(ggplot2)
library(dplyr)
library(tidyr)

sRNA_validation_raw_data = read.csv('Data/Figure_2_Panel_D/sRNA_PlateReader.csv', header = TRUE)
sRNA_validation_raw_data = sRNA_validation_raw_data[c(1:9216),]
sRNA_validation_raw_data = sRNA_validation_raw_data[,c(1:6)]

sRNA_Data_18 = subset(sRNA_validation_raw_data, Van_Concentration_uM  == 18)
sRNA_Data_18$OC6_Concentration_nM = sRNA_Data_18$OC6_Concentration_nM

M9               = subset(sRNA_Data_18, Genotype == "M9")
Backbone         = subset(sRNA_Data_18, Genotype == "Marionette-Clo pTU2 Backbone")
GFP_sRNA_18_D   = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_18_D")
GFP_sRNA_18_U   = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_18_U")
GFP_sRNA_20_D   = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_20_D")
GFP_sRNA_20_U   = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_20_U")
GFP_sRNA_22_RD  = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_22_RD")
GFP_sRNA_22_U   = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_22_U")
GFP_sRNA_24_RD  = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_24_RD")
GFP_sRNA_24_RU  = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_24_RU")
GFP_sRNA_24_U   = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_24_U")
GFP_sRNA_25_D   = subset(sRNA_Data_18, Genotype == "Marionette-Clo PVanCC-GFP PLuxB-sRNA_25_D")


GFP_sRNA_18_D$Normalised_GFP = (GFP_sRNA_18_D$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_18_D$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 
GFP_sRNA_18_U$Normalised_GFP = (GFP_sRNA_18_U$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_18_U$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 
GFP_sRNA_20_D$Normalised_GFP = (GFP_sRNA_20_D$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_20_D$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 
GFP_sRNA_20_U$Normalised_GFP = (GFP_sRNA_20_U$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_20_U$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 
GFP_sRNA_22_RD$Normalised_GFP = (GFP_sRNA_22_RD$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_22_RD$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 
GFP_sRNA_22_U$Normalised_GFP = (GFP_sRNA_22_U$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_22_U$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 
GFP_sRNA_24_RD$Normalised_GFP = (GFP_sRNA_24_RD$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_24_RD$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 
GFP_sRNA_24_RU$Normalised_GFP = (GFP_sRNA_24_RU$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_24_RU$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 
GFP_sRNA_24_U$Normalised_GFP = (GFP_sRNA_24_U$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_24_U$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 
GFP_sRNA_25_D$Normalised_GFP = (GFP_sRNA_25_D$GFP_Fluorescence - M9$GFP_Fluorescence) / (GFP_sRNA_25_D$OD600 - M9$OD600) - (Backbone$GFP_Fluorescence - M9$GFP_Fluorescence) / (Backbone$OD600 - M9$OD600) 


# A = ggplot(GFP_sRNA_18_D, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
#   geom_line(size = 2) + 
#   labs(
#     x = "Time (Hr)", 
#     y = "Normalised GFP Fluorescence Intensity (AU)", 
#     color = expression("[OC6] (" * "n" * "M)"),  
#     fill = expression("[OC6] (" * "n" * "M)")   
#   )  + 
#   ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-18D-GFP") +  
#   theme_minimal() + 
#   theme(
#     text = element_text(size = 20),          
#     axis.title = element_text(size = 20, hjust = 0.5),  
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
#     legend.title = element_text(size = 16),  
#     legend.text = element_text(size = 14),   
#     axis.line = element_line(size = 1, color = "black"),  
#     axis.ticks = element_line(size = 1, color = "black"),  
#     panel.grid.major = element_blank(),      
#     panel.grid.minor = element_blank(),      
#     panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
#   ) + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
#   scale_fill_gradientn(colors = c("green", "darkgreen", "black")) 
# 
# print(A)
# 
# 
# 
# svg("Syn-ChiX-18D-GFP.svg", width = 11.25, height = 8, bg = "white")
# print(A) 
# dev.off()
# 
# 
# A = ggplot(GFP_sRNA_18_U, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
#   geom_line(size = 2) + 
#   labs(
#     x = "Time (Hr)", 
#     y = "Normalised GFP Fluorescence Intensity (AU)", 
#     color = expression("[OC6] (" * "n" * "M)"),  
#     fill = expression("[OC6] (" * "n" * "M)")   
#   )  + 
#   ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-18U-GFP") +  
#   theme_minimal() + 
#   theme(
#     text = element_text(size = 20),          
#     axis.title = element_text(size = 20, hjust = 0.5),  
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
#     legend.title = element_text(size = 16),  
#     legend.text = element_text(size = 14),   
#     axis.line = element_line(size = 1, color = "black"),  
#     axis.ticks = element_line(size = 1, color = "black"),  
#     panel.grid.major = element_blank(),      
#     panel.grid.minor = element_blank(),      
#     panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
#   ) + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
#   scale_fill_gradientn(colors = c("green", "darkgreen", "black"))    
# 
# A
# 
# 
# 
# svg("Syn-ChiX-18U-GFP.svg", width = 11.25, height = 8, bg = "white")
# print(A) 
# dev.off()


A = ggplot(GFP_sRNA_20_D, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
  geom_line(size = 2) + 
  labs(
    x = "Time (Hr)", 
    y = "Normalised GFP Fluorescence Intensity (AU)", 
    color = expression("[OC6] (" * "n" * "M)"),  
    fill = expression("[OC6] (" * "n" * "M)")   
  )  + 
  ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-20D-GFP") +  
  theme_minimal() + 
  theme(
    text = element_text(size = 20),          
    axis.title = element_text(size = 20, hjust = 0.5),  
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14),   
    axis.line = element_line(size = 1, color = "black"),  
    axis.ticks = element_line(size = 1, color = "black"),  
    panel.grid.major = element_blank(),      
    panel.grid.minor = element_blank(),      
    panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
  ) + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
  scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
  scale_fill_gradientn(colors = c("green", "darkgreen", "black"))    

A



svg("Syn-ChiX-20D-GFP.svg", width = 11.25, height = 8, bg = "white")
print(A) 
dev.off()



# A = ggplot(GFP_sRNA_20_U, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
#   geom_line(size = 2) + 
#   labs(
#     x = "Time (Hr)", 
#     y = "Normalised GFP Fluorescence Intensity (AU)", 
#     color = expression("[OC6] (" * "n" * "M)"),  
#     fill = expression("[OC6] (" * "n" * "M)")   
#   )  + 
#   ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-20U-GFP") +  
#   theme_minimal() + 
#   theme(
#     text = element_text(size = 20),          
#     axis.title = element_text(size = 20, hjust = 0.5),  
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
#     legend.title = element_text(size = 16),  
#     legend.text = element_text(size = 14),   
#     axis.line = element_line(size = 1, color = "black"),  
#     axis.ticks = element_line(size = 1, color = "black"),  
#     panel.grid.major = element_blank(),      
#     panel.grid.minor = element_blank(),      
#     panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
#   ) + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
#   scale_fill_gradientn(colors = c("green", "darkgreen", "black"))    
# 
# A
# 
# 
# 
# svg("Syn-ChiX-20U-GFP.svg", width = 11.25, height = 8, bg = "white")
# print(A) 
# dev.off()
# 
# A = ggplot(GFP_sRNA_22_RD, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
#   geom_line(size = 2) + 
#   labs(
#     x = "Time (Hr)", 
#     y = "Normalised GFP Fluorescence Intensity (AU)", 
#     color = expression("[OC6] (" * "n" * "M)"),  
#     fill = expression("[OC6] (" * "n" * "M)")   
#   )  + 
#   ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-22RD-GFP") +  
#   theme_minimal() + 
#   theme(
#     text = element_text(size = 20),          
#     axis.title = element_text(size = 20, hjust = 0.5),  
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
#     legend.title = element_text(size = 16),  
#     legend.text = element_text(size = 14),   
#     axis.line = element_line(size = 1, color = "black"),  
#     axis.ticks = element_line(size = 1, color = "black"),  
#     panel.grid.major = element_blank(),      
#     panel.grid.minor = element_blank(),      
#     panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
#   ) + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
#   scale_fill_gradientn(colors = c("green", "darkgreen", "black"))    
# 
# A
# 
# 
# 
# svg("Syn-ChiX-22RD-GFP.svg", width = 11.25, height = 8, bg = "white")
# print(A) 
# dev.off()
# 
# A = ggplot(GFP_sRNA_22_U, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
#   geom_line(size = 2) + 
#   labs(
#     x = "Time (Hr)", 
#     y = "Normalised GFP Fluorescence Intensity (AU)", 
#     color = expression("[OC6] (" * "n" * "M)"),  
#     fill = expression("[OC6] (" * "n" * "M)")   
#   )  + 
#   ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-22U-GFP") +  
#   theme_minimal() + 
#   theme(
#     text = element_text(size = 20),          
#     axis.title = element_text(size = 20, hjust = 0.5),  
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
#     legend.title = element_text(size = 16),  
#     legend.text = element_text(size = 14),   
#     axis.line = element_line(size = 1, color = "black"),  
#     axis.ticks = element_line(size = 1, color = "black"),  
#     panel.grid.major = element_blank(),      
#     panel.grid.minor = element_blank(),      
#     panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
#   ) + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
#   scale_fill_gradientn(colors = c("green", "darkgreen", "black"))    
# 
# A
# 
# 
# 
# svg("Syn-ChiX-22U-GFP.svg", width = 11.25, height = 8, bg = "white")
# print(A) 
# dev.off()
# 
# A = ggplot(GFP_sRNA_24_RD, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
#   geom_line(size = 2) + 
#   labs(
#     x = "Time (Hr)", 
#     y = "Normalised GFP Fluorescence Intensity (AU)", 
#     color = expression("[OC6] (" * "n" * "M)"),  
#     fill = expression("[OC6] (" * "n" * "M)")   
#   )  + 
#   ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-24RD-GFP") +  
#   theme_minimal() + 
#   theme(
#     text = element_text(size = 20),          
#     axis.title = element_text(size = 20, hjust = 0.5),  
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
#     legend.title = element_text(size = 16),  
#     legend.text = element_text(size = 14),   
#     axis.line = element_line(size = 1, color = "black"),  
#     axis.ticks = element_line(size = 1, color = "black"),  
#     panel.grid.major = element_blank(),      
#     panel.grid.minor = element_blank(),      
#     panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
#   ) + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
#   scale_fill_gradientn(colors = c("green", "darkgreen", "black"))    
# 
# A
# 
# 
# 
# svg("Syn-ChiX-24RD-GFP.svg", width = 11.25, height = 8, bg = "white")
# print(A) 
# dev.off()
# 
# A = ggplot(GFP_sRNA_24_RU, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
#   geom_line(size = 2) + 
#   labs(
#     x = "Time (Hr)", 
#     y = "Normalised GFP Fluorescence Intensity (AU)", 
#     color = expression("[OC6] (" * "n" * "M)"),  
#     fill = expression("[OC6] (" * "n" * "M)")   
#   )  + 
#   ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-24RU-GFP") +  
#   theme_minimal() + 
#   theme(
#     text = element_text(size = 20),          
#     axis.title = element_text(size = 20, hjust = 0.5),  
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
#     legend.title = element_text(size = 16),  
#     legend.text = element_text(size = 14),   
#     axis.line = element_line(size = 1, color = "black"),  
#     axis.ticks = element_line(size = 1, color = "black"),  
#     panel.grid.major = element_blank(),      
#     panel.grid.minor = element_blank(),      
#     panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
#   ) + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
#   scale_fill_gradientn(colors = c("green", "darkgreen", "black"))    
# 
# A
# 
# 
# 
# svg("Syn-ChiX-24RU-GFP.svg", width = 11.25, height = 8, bg = "white")
# print(A) 
# dev.off()
# 
# A = ggplot(GFP_sRNA_24_U, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
#   geom_line(size = 2) + 
#   labs(
#     x = "Time (Hr)", 
#     y = "Normalised GFP Fluorescence Intensity (AU)", 
#     color = expression("[OC6] (" * "n" * "M)"),  
#     fill = expression("[OC6] (" * "n" * "M)")   
#   )  + 
#   ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-24U-GFP") +  
#   theme_minimal() + 
#   theme(
#     text = element_text(size = 20),          
#     axis.title = element_text(size = 20, hjust = 0.5),  
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
#     legend.title = element_text(size = 16),  
#     legend.text = element_text(size = 14),   
#     axis.line = element_line(size = 1, color = "black"),  
#     axis.ticks = element_line(size = 1, color = "black"),  
#     panel.grid.major = element_blank(),      
#     panel.grid.minor = element_blank(),      
#     panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
#   ) + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
#   scale_fill_gradientn(colors = c("green", "darkgreen", "black"))    
# 
# A
# 
# 
# 
# svg("Syn-ChiX-24U-GFP.svg", width = 11.25, height = 8, bg = "white")
# print(A) 
# dev.off()
# 
# A = ggplot(GFP_sRNA_25_D, aes(x = Time_Hr, y = Normalised_GFP, group = OC6_Concentration_nM, color = OC6_Concentration_nM)) +
#   geom_line(size = 2) + 
#   labs(
#     x = "Time (Hr)", 
#     y = "Normalised GFP Fluorescence Intensity (AU)", 
#     color = expression("[OC6] (" * "n" * "M)"),  
#     fill = expression("[OC6] (" * "n" * "M)")   
#   )  + 
#   ggtitle("PVanCC_GFPmut3 PLuxB_Syn-ChiX-25D-GFP") +  
#   theme_minimal() + 
#   theme(
#     text = element_text(size = 20),          
#     axis.title = element_text(size = 20, hjust = 0.5),  
#     plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),  
#     legend.title = element_text(size = 16),  
#     legend.text = element_text(size = 14),   
#     axis.line = element_line(size = 1, color = "black"),  
#     axis.ticks = element_line(size = 1, color = "black"),  
#     panel.grid.major = element_blank(),      
#     panel.grid.minor = element_blank(),      
#     panel.grid = element_line(colour = "grey90", size = 0.5, linetype = "dashed") 
#   ) + 
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +  
#   scale_color_gradientn(colors = c("green", "darkgreen", "black")) +  
#   scale_fill_gradientn(colors = c("green", "darkgreen", "black"))    
# 
# A
# 
# 
# 
# svg("Syn-ChiX-25D-GFP.svg", width = 11.25, height = 8, bg = "white")
# print(A) 
# dev.off()
# 
