getwd()
workdir <- "C:/Users/hebam/OneDrive - University of Calgary/Mitacs/"
workdir <- "C:/Users/Owner/OneDrive - University of Calgary/Mitacs/"
setwd(workdir)

library(sf)
library(ggplot2)
library(stars)
library(dplyr)
library(rgdal)
library(ggthemes)
library(stringr)
library(ncdf4)
library(ncdf4.helpers)
library(PCICt)

library(Lmoments)
library(hydroGOF)
library(imager)

# out_WG <- readRDS("./eval/eval_metrics_WGANGonly_400dayNZ.rds")
# out_G <- readRDS("./eval/eval_metrics_Gonly_400dayNZ.rds")



outSS_WG <- readRDS("./eval/eval_metrics_WGAN.rds")
outSS_G <- readRDS("./eval/eval_metrics_Generator.rds")
outSS_G$model <- "UNET"
outSS_WGAN_NZ <- readRDS("./eval/eval_metrics_WGAN_NZ.rds")
outSS_WGAN_NZ_DT <- readRDS("./eval/eval_metrics_WGAN_NZ_DT.rds")

outSS_WGAN_NZ_DT$model <- "WGAN_NZ_DT"


outSS_G <- na.omit(outSS_G)
outSS_WG <- na.omit(outSS_WG)
outSS_WGAN_NZ <- na.omit(outSS_WGAN_NZ)
outSS_WGAN_NZ_DT <- na.omit(outSS_WGAN_NZ_DT)

sum(outSS_G$ssim < 0.8)/nrow(outSS_G) * 100
sum(outSS_WG$ssim < 0.8)/nrow(outSS_WG) * 100
sum(outSS_WGAN_NZ$ssim < 0.8)/nrow(outSS_WGAN_NZ) * 100
sum(outSS_WGAN_NZ_DT$ssim < 0.8)/nrow(outSS_WGAN_NZ_DT) * 100

meanI_g <- mean(outSS_G$ssim)
meanI_wgan <- mean(outSS_WG$ssim)
meanI_wganNZ <- mean(outSS_WGAN_NZ$ssim)
meanI_wganDT <- mean(outSS_WGAN_NZ_DT$ssim)

SSIM_all <- rbind(outSS_WG, outSS_G, outSS_WGAN_NZ, outSS_WGAN_NZ_DT)
SSIM_all$mod <- factor(SSIM_all$mod, levels = c("UNET", "WGAN", "WGAN_NZ","WGAN_NZ_DT"))

ggplot() +
  geom_violin(data = SSIM_all, aes(x = model, y = ssim, fill = model), alpha = 0.7) + 
  geom_boxplot(data = SSIM_all, aes(x = model, y = ssim, fill = model), width=0.05,
               color="grey50", alpha=0.2, outlier.shape = NA) + 
  scale_fill_manual("", values =  c("UNET" = "#60a862",
                                    "WGAN" = "#727cce",
                                    "WGAN_NZ" = "#b4943e",
                                    "WGAN_NZ_DT" ="#cb5a4c"),
                    breaks = c("UNET" ,
                               "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  xlab("Model") + 
  ylab("SSIM") + 
  # scale_y_continuous(limits = c(0.9, 1)) + 
  theme(#aspect.ratio = 0.3,
    legend.text = element_text(size = 7),
    axis.title.x = element_text(size = 9, colour = gray(0.25)),
    axis.title.y = element_text(size = 9, colour = gray(0.25)),
    legend.title = element_text(size = 7, colour = gray(0.25)),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = gray(0.25), size = 0.3),
    axis.text.x = element_text(size = 8, colour = gray(0.25)),
    axis.text.y = element_text(size = 8, colour = gray(0.25)),
    plot.title = element_text(hjust = 0, size = 6),
    plot.margin=unit(c(3,3,0.5,0.5), "mm"),
    axis.ticks.length.x = unit(-0.5, "mm"),
    axis.ticks.length.y = unit(-0.5, "mm"),
    axis.ticks = element_line(color = gray(0.25), size = 0.3),
    legend.position="none",
    legend.spacing.x = unit(3.5, "mm"),
    legend.spacing.y = unit(0.1, "mm"),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    # legend.background = element_rect(color = "transparent", fill = gray(0.95)),
    legend.key.size = unit(3,"mm"),
    legend.key.width = unit(3,"mm"),
    # legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
    text = element_text(family = "Calibri", color = gray(0.25)),
    legend.box.margin=margin(0.5,0.5,0.5,0.5)) 

ggsave("./eval/SSIM_ALL_plot_final_namescor.jpg",  height = 12, width = 15, units = "cm", dpi = 400)

