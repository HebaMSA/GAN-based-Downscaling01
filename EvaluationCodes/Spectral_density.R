getwd()
workdir <- "C:/Users/hebam/OneDrive - University of Calgary/Mitacs/"
setwd(workdir)

library(sf)
library(ggplot2)
library(stars)
library(dplyr)
library(ggthemes)
library(stringr)
library(ncdf4)
library(ncdf4.helpers)
library(PCICt)

library(Lmoments)
library(hydroGOF)
library(imager)

library(seewave)
library(tseries)


nctodf <- function(nc){
  fn <- nc
  YourBrick <- raster::brick(fn)
  
  points <- raster::rasterToPoints(YourBrick)
  pts <- as.data.frame(points)
  
  #Removing NA points
  pts1 <- subset(pts,pts$x != "NA",)
  pts <- subset(pts1, pts$y != "NA",)
  
  return(pts)
}
##############################################

files_origALL <- as.list(list.files(pattern = "\\.nc$",
                                    path = "./data/Orig_all/"))

# files_orig <- files_origALL[2801:3200]
files_orig <- files_origALL[13612:20416]



out_df <- data.frame()

for (j in 1:length(files_orig)){
  
  # test_pred <- nctodf(paste0("./codes/nc_files_wganG30_hope3/sample_",j,".nc"))
  test_pred <- nctodf(paste0("./codes/nc_files_wganG30_ALL/sample_",j,".nc"))
  test_orig <- nctodf(paste0("./data/Orig_all/",files_orig[[j]]))
  
  stats_pred <- test_pred %>% 
    select(-c(x,y)) %>% 
    summarise_all(~list(min = as.numeric(quantile(as.numeric(unlist(.)), 0.025, na.rm = T)),
                        max = as.numeric(quantile(as.numeric(unlist(.)), 0.975, na.rm = T)),
                        mean = as.numeric(mean(as.numeric(unlist(.)), na.rm = T))))
  
  
  stats_orig <- test_orig %>% 
    select(-c(x,y)) %>% 
    summarise_all(~list(min = as.numeric(quantile(as.numeric(unlist(.)), 0.025, na.rm = T)) * 3600,
                        max = as.numeric(quantile(as.numeric(unlist(.)), 0.975, na.rm = T)) * 3600,
                        mean = as.numeric(mean(as.numeric(unlist(.)), na.rm = T)) * 3600))
  
  time_seriesP <- as.numeric(unlist(stats_pred[3,c(1:ncol(stats_pred))]))
  time_seriesO <- as.numeric(unlist(stats_orig[3,c(1:ncol(stats_orig))]))
  
  # time_seriesP <- as.numeric(unlist(test_pred[1,c(3:ncol(test_pred))]))
  # time_seriesO <- as.numeric(unlist(test_orig[1,c(3:ncol(test_orig))])) * 3600
  # Compute and plot spectral density
  spectrum_resultP <- stats::spectrum(time_seriesP, main = "Spectral Density",
                                      log = "no", # Use log scale for frequency axis if needed
                                      span = 3)   # Span controls smoothing
  
  
  
  spectrum_resultO <- stats::spectrum(time_seriesO, main = "Spectral Density",
                                      log = "no", # Use log scale for frequency axis if needed
                                      span = 3)   # Span controls smoothing
  

  
  
  # Use the 'coh' function from the 'seewave' package
  coherence <- coh(time_seriesO, time_seriesP, f = 1, plot = TRUE)
  
  ksRes <- ks.test(spectrum_resultO$spec, spectrum_resultP$spec)
  
  
  acfResO <- acf(time_seriesO, lag = 5)
  acfResP <- acf(time_seriesP, lag = 5)
  
  bias_acf <- mean(acfResP[["acf"]] - acfResO[["acf"]])
  rmse_acf <- rmse(as.numeric(unlist(acfResP[["acf"]])), as.numeric(unlist(acfResO[["acf"]])))
 
  out_df[j,1] <- mean(coherence[1:12,2])
  out_df[j,2] <- ksRes[["statistic"]][["D"]]
  out_df[j,3] <- bias_acf
  out_df[j,4] <- rmse_acf
}

colnames(out_df) <- c("coh_avg","ks_Pvalue","bias_acf","rmse_acf")
out_df$mod <- "WGAN"

out_df_WGAN <- out_df
saveRDS(out_df_WGAN, "./eval/spec_WGAN_stat.rds")





out_df <- data.frame()

for (j in 1:length(files_orig)){
  
  test_pred <- nctodf(paste0("./codes/nc_files_wganG30_ALLNZ/sample_",j,".nc"))
  # test_pred <- nctodf(paste0("./codes/nc_filesALL/sample_",j,".nc"))
  test_orig <- nctodf(paste0("./data/Orig_all/",files_orig[[j]]))
  
  stats_pred <- test_pred %>% 
    select(-c(x,y)) %>% 
    summarise_all(~list(min = as.numeric(quantile(as.numeric(unlist(.)), 0.025, na.rm = T)),
                        max = as.numeric(quantile(as.numeric(unlist(.)), 0.975, na.rm = T)),
                        mean = as.numeric(mean(as.numeric(unlist(.)), na.rm = T))))
  
  
  stats_orig <- test_orig %>% 
    select(-c(x,y)) %>% 
    summarise_all(~list(min = as.numeric(quantile(as.numeric(unlist(.)), 0.025, na.rm = T)) * 3600,
                        max = as.numeric(quantile(as.numeric(unlist(.)), 0.975, na.rm = T)) * 3600,
                        mean = as.numeric(mean(as.numeric(unlist(.)), na.rm = T)) * 3600))
  
  time_seriesP <- as.numeric(unlist(stats_pred[3,c(1:ncol(stats_pred))]))
  time_seriesO <- as.numeric(unlist(stats_orig[3,c(1:ncol(stats_orig))]))
  
  # time_seriesP <- as.numeric(unlist(test_pred[1,c(3:ncol(test_pred))]))
  # time_seriesO <- as.numeric(unlist(test_orig[1,c(3:ncol(test_orig))])) * 3600
  # Compute and plot spectral density
  spectrum_resultP <- stats::spectrum(time_seriesP, main = "Spectral Density",
                                      log = "no", # Use log scale for frequency axis if needed
                                      span = 3)   # Span controls smoothing
  
  
  
  spectrum_resultO <- stats::spectrum(time_seriesO, main = "Spectral Density",
                                      log = "no", # Use log scale for frequency axis if needed
                                      span = 3)   # Span controls smoothing
  
  
  
  
  # Use the 'coh' function from the 'seewave' package
  coherence <- coh(time_seriesO, time_seriesP, f = 1, plot = FALSE)
  
  ksRes <- ks.test(spectrum_resultO$spec, spectrum_resultP$spec)
  
  
  acfResO <- acf(time_seriesO, lag = 5, plot = FALSE)
  acfResP <- acf(time_seriesP, lag = 5, plot = FALSE)
  
  bias_acf <- mean(acfResP[["acf"]] - acfResO[["acf"]])
  rmse_acf <- rmse(as.numeric(unlist(acfResP[["acf"]])), as.numeric(unlist(acfResO[["acf"]])))
  
  out_df[j,1] <- mean(coherence[1:12,2])
  out_df[j,2] <- ksRes[["statistic"]][["D"]] #if above 0.05, therefore the two spectral densities highly likely comes from the same distribution at 95% confidence interval
  out_df[j,3] <- bias_acf
  out_df[j,4] <- rmse_acf
}


colnames(out_df) <- c("coh_avg","ks_Pvalue","bias_acf","rmse_acf")

out_df$mod <- "WGAN_NZ"

out_df_G <- out_df
saveRDS(out_df_G, "./eval/spec_WGAN_NZ_stat.rds")




out_df <- data.frame()

for (j in 1:length(files_orig)){
  
  # test_pred <- nctodf(paste0("./codes/nc_files_wganG30_hope3/sample_",j,".nc"))
  test_pred <- nctodf(paste0("./codes/nc_filesALL/sample_",j,".nc"))
  test_orig <- nctodf(paste0("./data/Orig_all/",files_orig[[j]]))
  
  stats_pred <- test_pred %>% 
    select(-c(x,y)) %>% 
    summarise_all(~list(min = as.numeric(quantile(as.numeric(unlist(.)), 0.025, na.rm = T)),
                        max = as.numeric(quantile(as.numeric(unlist(.)), 0.975, na.rm = T)),
                        mean = as.numeric(mean(as.numeric(unlist(.)), na.rm = T))))
  
  
  stats_orig <- test_orig %>% 
    select(-c(x,y)) %>% 
    summarise_all(~list(min = as.numeric(quantile(as.numeric(unlist(.)), 0.025, na.rm = T)) * 3600,
                        max = as.numeric(quantile(as.numeric(unlist(.)), 0.975, na.rm = T)) * 3600,
                        mean = as.numeric(mean(as.numeric(unlist(.)), na.rm = T)) * 3600))
  
  time_seriesP <- as.numeric(unlist(stats_pred[3,c(1:ncol(stats_pred))]))
  time_seriesO <- as.numeric(unlist(stats_orig[3,c(1:ncol(stats_orig))]))
  
  # time_seriesP <- as.numeric(unlist(test_pred[1,c(3:ncol(test_pred))]))
  # time_seriesO <- as.numeric(unlist(test_orig[1,c(3:ncol(test_orig))])) * 3600
  # Compute and plot spectral density
  spectrum_resultP <- stats::spectrum(time_seriesP, main = "Spectral Density",
                                      log = "no", # Use log scale for frequency axis if needed
                                      span = 3)   # Span controls smoothing
  
  
  
  spectrum_resultO <- stats::spectrum(time_seriesO, main = "Spectral Density",
                                      log = "no", # Use log scale for frequency axis if needed
                                      span = 3)   # Span controls smoothing
  
  
  
  
  # Use the 'coh' function from the 'seewave' package
  coherence <- coh(time_seriesO, time_seriesP, f = 1, plot = TRUE)
  
  ksRes <- ks.test(spectrum_resultO$spec, spectrum_resultP$spec)
  
  
  acfResO <- acf(time_seriesO, lag = 5, plot = FALSE)
  acfResP <- acf(time_seriesP, lag = 5, plot = FALSE)
  
  bias_acf <- mean(acfResP[["acf"]] - acfResO[["acf"]])
  rmse_acf <- rmse(as.numeric(unlist(acfResP[["acf"]])), as.numeric(unlist(acfResO[["acf"]])))
  
  out_df[j,1] <- mean(coherence[1:12,2])
  out_df[j,2] <- ksRes[["statistic"]][["D"]] #if above 0.05, therefore the two spectral densities highly likely comes from the same distribution at 95% confidence interval
  out_df[j,3] <- bias_acf
  out_df[j,4] <- rmse_acf
}

colnames(out_df) <- c("coh_avg","ks_Pvalue","bias_acf","rmse_acf")
out_df$mod <- "Generator_MSE"

out_df$mod <- "WGAN_h3"

out_df_G <- out_df
saveRDS(out_df_G, "./eval/spec_G_stat.rds")



###PLOT

out_Gen <- readRDS( "./eval/spec_G_stat.rds")
out_Gen$mod <- "UNET"
out_WGAN <- readRDS( "./eval/spec_WGAN_stat.rds")
out_WGAN_NZ <- readRDS( "./eval/spec_WGAN_NZ_stat.rds")
out_WGAN_NZ_DT <- readRDS( "./eval/spec_WGAN_h3_stat.rds")
out_WGAN_NZ_DT$mod <- "WGAN_NZ_DT"

out_all <- rbind(out_Gen, out_WGAN, out_WGAN_NZ, out_WGAN_NZ_DT)
out_all$mod <- factor(out_all$mod, levels = c("UNET", "WGAN", "WGAN_NZ","WGAN_NZ_DT"))


g1 <- ggplot() + 
  geom_boxplot(
    data = out_all, 
    aes(
      x = as.factor(mod), 
      y = coh_avg, 
      group = as.factor(mod),  # Group by both lags and mod
      fill = mod  # Color by mod
    ), 
    linewidth = 0.1,
    alpha = 0.7,
    outlier.shape = NA
  ) + 
  scale_fill_manual("", values =  c("UNET" = "#60a862",
                                    "WGAN" = "#727cce",
                                    "WGAN_NZ" = "#b4943e",
                                    "WGAN_NZ_DT" ="#cb5a4c"),
                    breaks = c("UNET" ,
                               "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  xlab("") + 
  ylab("Coherence coefficient") +
  theme(#aspect.ratio = 0.3,
    legend.text = element_text(size = 6),
    axis.title.x = element_text(size = 8, colour = gray(0.25)),
    axis.title.y = element_text(size = 8, colour = gray(0.25)),
    legend.title = element_text(size = 7, colour = gray(0.25)),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = gray(0.25), size = 0.3),
    axis.text.x = element_text(size = 6, colour = gray(0.25)),
    axis.text.y = element_text(size = 6, colour = gray(0.25)),
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

g1




g2 <- ggplot() + 
  geom_boxplot(
    data = out_all, 
    aes(
      x = as.factor(mod), 
      y = ks_Pvalue, 
      group = as.factor(mod),  # Group by both lags and mod
      fill = mod  # Color by mod
    ), 
    linewidth = 0.1,
    alpha = 0.7,
    outlier.shape = NA
  ) + 
  scale_fill_manual("", values =  c("UNET" = "#60a862",
                                    "WGAN" = "#727cce",
                                    "WGAN_NZ" = "#b4943e",
                                    "WGAN_NZ_DT" ="#cb5a4c"),
                    breaks = c("UNET" ,
                               "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  xlab("") + 
  ylab("KS statistic") +
  theme(#aspect.ratio = 0.3,
    legend.text = element_text(size = 6),
    axis.title.x = element_text(size = 8, colour = gray(0.25)),
    axis.title.y = element_text(size = 8, colour = gray(0.25)),
    legend.title = element_text(size = 7, colour = gray(0.25)),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = gray(0.25), size = 0.3),
    axis.text.x = element_text(size = 6, colour = gray(0.25)),
    axis.text.y = element_text(size = 6, colour = gray(0.25)),
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

g2






g3 <- ggplot() + 
  geom_boxplot(
    data = out_all, 
    aes(
      x = as.factor(mod), 
      y = bias_acf, 
      group = as.factor(mod),  # Group by both lags and mod
      fill = mod  # Color by mod
    ), 
    linewidth = 0.1,
    alpha = 0.7,
    outlier.shape = NA
  ) + 
  scale_fill_manual("", values =  c("UNET" = "#60a862",
                                    "WGAN" = "#727cce",
                                    "WGAN_NZ" = "#b4943e",
                                    "WGAN_NZ_DT" ="#cb5a4c"),
                    breaks = c("UNET" ,
                               "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  scale_y_continuous(limits = c(-0.2, 0.2))  + 
  xlab("") + 
  ylab("ACF bias") +
  theme(#aspect.ratio = 0.3,
    legend.text = element_text(size = 6),
    axis.title.x = element_text(size = 8, colour = gray(0.25)),
    axis.title.y = element_text(size = 8, colour = gray(0.25)),
    legend.title = element_text(size = 7, colour = gray(0.25)),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = gray(0.25), size = 0.3),
    axis.text.x = element_text(size = 6, colour = gray(0.25)),
    axis.text.y = element_text(size = 6, colour = gray(0.25)),
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

g3





g4 <- ggplot() + 
  geom_boxplot(
    data = out_all, 
    aes(
      x = as.factor(mod), 
      y = rmse_acf, 
      group = as.factor(mod),  # Group by both lags and mod
      fill = mod  # Color by mod
    ), 
    linewidth = 0.1,
    alpha = 0.7,
    outlier.shape = NA
  ) + 
  scale_fill_manual("", values =  c("UNET" = "#60a862",
                                    "WGAN" = "#727cce",
                                    "WGAN_NZ" = "#b4943e",
                                    "WGAN_NZ_DT" ="#cb5a4c"),
                    breaks = c("UNET" ,
                               "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  scale_y_continuous(limits = c(-0.05, 0.3))  +
  xlab("") + 
  ylab("ACF RMSE") +
  theme(#aspect.ratio = 0.3,
    legend.text = element_text(size = 6),
    axis.title.x = element_text(size = 8, colour = gray(0.25)),
    axis.title.y = element_text(size = 8, colour = gray(0.25)),
    legend.title = element_text(size = 7, colour = gray(0.25)),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = gray(0.25), size = 0.3),
    axis.text.x = element_text(size = 6, colour = gray(0.25)),
    axis.text.y = element_text(size = 6, colour = gray(0.25)),
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

g4



library(ggpubr)
ggarrange(g1, g2, g3,g4,
          nrow = 2, 
          ncol = 2,
          labels = c("(a)",
                     "(b)",
                     "(c)",
                     "(d)"),
          align = "hv",
          vjust = c(1.5,1.5,1.5,1.5),
          hjust = c(-3.2,-3.2,-3.2,-3.2),
          # widths = c(1,0.7, 0.7),
          font.label = list(size = 8,
                            color = gray(0.25), 
                            face = "plain",
                            base_family = 'Calibri'))

ggsave("./eval/ACF_Alternatives_stats_namescor.jpg",  height = 12, width = 13.5, units = "cm", dpi = 400)
