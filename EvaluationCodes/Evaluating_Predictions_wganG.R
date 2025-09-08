getwd()
workdir <- "C:/Users/hebam/OneDrive - University of Calgary/Mitacs/"
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


# Define an SSIM function
calculate_ssim_WGAN <- function(img1, img2) {
  mean1 <- mean(img1)
  mean2 <- mean(img2)
  var1 <- var(img1)
  var2 <- var(img2)
  cov12 <- cov(img1, img2)
  
  # Constants for stability
  #L <- 22  # Dynamic range of pixel values, change based on image type (For 400 days)
  L <- 39  # Dynamic range of pixel values, change based on image type (For 6805 days)
  C1 <- (0.01 * L)^2
  C2 <- (0.03 * L)^2
  
  # SSIM calculation
  ssim <- ((2 * mean1 * mean2 + C1) * (2 * cov12 + C2)) /
    ((mean1^2 + mean2^2 + C1) * (var1 + var2 + C2))
  return(ssim)
}

#visual test

out_df <- readRDS("./eval/eval_metrics_WGANGonly_6805dayNZ.rds")
#lowest ssim = 19880716.nc_10
#second lowes ssim = 19900810.nc_9

which(files_orig == "day_monthly_P19900810.nc")

test_pred <- nctodf("./codes/nc_files_wganG30_ALL/sample_1204.nc")
test_orig <- nctodf("./data/Orig_all/day_monthly_P19900810.nc")

# min_val <- min(test_pred$X0, test_orig$X1987.11.15.8 * 3600)
# max_val <- max(test_pred$X0, test_orig$X1987.11.15.8 * 3600)

min_val <- min(test_pred$X8, test_orig$X1990.08.10.9 * 3600)
max_val <- max(test_pred$X8, test_orig$X1990.08.10.9 * 3600)

p1 <- ggplot(test_pred, aes(x = x, y = y, fill = `X8`)) +
  geom_tile() +
  # scale_fill_viridis_c() +
  xlab("longitude") + 
  ylab("latitude") + 
  scale_fill_gradientn("Hourly Precipitation (mm)",
                       colors = c(
                         "#f7f6f5",
                         "#ece8e5",
                         "#d0b5ae",
                         "#b78b81",
                         "#a56f67",
                         "#965c57",
                         "#864d4b",
                         "#773f42",
                         "#693139",
                         "#592631",
                         "#4a1c2c",
                         "#391a3c",
                         "#32215b",
                         "#4943a0",
                         "#5054ae",
                         "#485486",
                         "#4c6083",
                         "#54708d",
                         "#6a8da3",
                         "#85aab7",
                         "#c3b799",
                         "#edb376",
                         "#e98a4a",
                         "#e3662e",
                         "#d83e1d",
                         "#c5281a",
                         "#ac1b1e",
                         "#961a24",
                         "#831b2b"                       ),
                       # alpha = 0.7,
                       trans = scales::pseudo_log_trans(),
                       limits = c(min_val,max_val),
                       # breaks = c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06),
                       guide = "colourbar") +
  labs(title = "Predicted Data") +
  coord_fixed() +
  theme_minimal() + theme(legend.text = element_text(size = 5),
                          legend.title = element_text(size = 6),
                          axis.title = element_text(size = 6),
                          axis.text = element_text(size = 5),
                          legend.text.align = 0.5,
                          legend.title.align=0.5,
                          legend.box.just = "center",
                          legend.justification = "center",
                          legend.position='bottom',
                          legend.key.size = unit(0.2, "cm"), legend.key.width = unit(1.2,"cm"),
                          plot.margin=unit(c(-0.5,-0.5,5,-0.5), "mm"),
                          legend.box.margin=unit(c(-3,-1.5,-1.5,-1.5), "mm"),
                          plot.title = element_text(size = 6, vjust = -2),
                          # legend.key.size = unit(7,"mm"),
                          # legend.background = element_rect(color = "transparent",
                          #                                  fill = gray(0.95)),
                          # legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
                          text = element_text(color = gray(0.25)))+ 
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

p1
# Plot the upscaled data (lower resolution)
p2 <- ggplot(test_orig, aes(x = x, y = y, fill = `X1990.08.10.9` * 3600)) +
  geom_tile() +
  xlab("longitude") + 
  ylab("latitude") + 
  scale_fill_gradientn("Hourly Precipitation (mm)",
                       colors = c(
                         "#f7f6f5",
                         "#ece8e5",
                         "#d0b5ae",
                         "#b78b81",
                         "#a56f67",
                         "#965c57",
                         "#864d4b",
                         "#773f42",
                         "#693139",
                         "#592631",
                         "#4a1c2c",
                         "#391a3c",
                         "#32215b",
                         "#4943a0",
                         "#5054ae",
                         "#485486",
                         "#4c6083",
                         "#54708d",
                         "#6a8da3",
                         "#85aab7",
                         "#c3b799",
                         "#edb376",
                         "#e98a4a",
                         "#e3662e",
                         "#d83e1d",
                         "#c5281a",
                         "#ac1b1e",
                         "#961a24",
                         "#831b2b"),
                       # alpha = 0.7,
                       trans = scales::pseudo_log_trans(),
                       limits = c(min_val,max_val),
                       # breaks = c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06),
                       guide = "colourbar") +
  labs(title = "Original Data") +
  coord_fixed() +
  theme_minimal() + theme(legend.text = element_text(size = 5),
                          legend.title = element_text(size = 6),
                          axis.title = element_text(size = 6),
                          axis.text = element_text(size = 5),
                          legend.text.align = 0.5,
                          legend.title.align=0.5,
                          legend.box.just = "center",
                          legend.justification = "center",
                          legend.position='bottom',
                          legend.key.size = unit(0.2, "cm"), legend.key.width = unit(1.2,"cm"),
                          plot.margin=unit(c(-0.5,-0.5,5,-0.5), "mm"),
                          legend.box.margin=unit(c(-3,-1.5,-1.5,-1.5), "mm"),
                          plot.title = element_text(size = 6, vjust = -2),
                          # legend.key.size = unit(7,"mm"),
                          # legend.background = element_rect(color = "transparent",
                          #                                  fill = gray(0.95)),
                          # legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
                          text = element_text(color = gray(0.25))) + 
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

p2

library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

ggsave("./eval/LowSSIM_eg2.jpg",  height = 4, width = 8, units = "cm", dpi = 400)












#Comprehensive Evaluation


# files_pred <- as.list(list.files(pattern = "\\.nc$",
#                                  path = "./codes/nc_files_wganG30_ALL/"))

files_origALL <- as.list(list.files(pattern = "\\.nc$",
                                    path = "./data/Orig_all/"))

files_orig <- files_origALL[13612:20416]
#files_orig <- files_origALL[2801:3200]

out_list <- list()

for (j in 1:length(files_orig)){
  
  # test_pred <- nctodf(paste0("./codes/nc_filesALL/sample_",j,".nc"))
  # test_pred <- nctodf(paste0("./codes/nc_files_wganG30_ALL/sample_",j,".nc"))
  # test_pred <- nctodf(paste0("./codes/nc_files_wganG30_ALLNZ/sample_",j,".nc"))
  test_pred <- nctodf(paste0("./codes/nc_files_wganG30_hope3/sample_",j,".nc"))
  
  test_orig <- nctodf(paste0("./data/Orig_all/",files_orig[[j]]))
  
  x_coords <- unique(test_pred$x)
  y_coords <- unique(test_pred$y)
  # Determine spatial dimensions
  nx <- length(x_coords)
  ny <- length(y_coords)
  
  # Reshape into a 3D array (nx x ny x time)
  time_data <- as.matrix(test_pred[, -(1:2)])  # Drop x and y columns
  arrayP <- array(time_data, dim = c(nx, ny, ncol(time_data)))
  
  time_data1 <- as.matrix(test_orig[, -(1:2)]) * 3600  # Drop x and y columns
  arrayO <- array(time_data1, dim = c(nx, ny, ncol(time_data1)))
  
  
  metrics_df <- data.frame()
  
  for (i in 1:24){
    
    one_tsP <- arrayP[,,i]
    one_tsO <- arrayO[,,i]
    
    p0_P <- sum(one_tsP <= 0.1)/length(one_tsP)
    p0_O <- sum(one_tsO <= 0.1)/length(one_tsO)
    
    NZone_tsP <-  one_tsP[one_tsP > 0.1]
    NZone_tsO <-  one_tsO[one_tsO > 0.1]
    
    if (length(NZone_tsO) < 30){
      lmo_O <- rep(NA,4)
    } else {
    lmo_O <- Lcoefs(NZone_tsO)
    }
    
    if (length(NZone_tsP) < 30){
      lmo_P <- rep(NA,4)
    } else {
      lmo_P <- Lcoefs(NZone_tsP)
    }
    
    # lmo_P <- Lcoefs(NZone_tsP)
    
    img1 <- as.cimg(one_tsO)  # Convert time step i to cimg format
    img2 <- as.cimg(one_tsP)  # Convert time step i to cimg format
    
    # Calculate SSIM
    ssim_value <- calculate_ssim_WGAN(img1, img2)
    
    metrics_df[i, 1] <- p0_P
    metrics_df[i, 2] <- p0_O
    metrics_df[i, 3] <- lmo_P[1]
    metrics_df[i, 4] <- lmo_O[1]
    
    metrics_df[i, 5] <- lmo_P[2]
    metrics_df[i, 6] <- lmo_O[2]
    
    metrics_df[i, 7] <- lmo_P[3]
    metrics_df[i, 8] <- lmo_O[3]
    
    metrics_df[i, 9] <- lmo_P[4]
    metrics_df[i, 10] <- lmo_O[4]
    
    metrics_df[i, 11] <- ssim_value
    
    metrics_df[i, 12] <- paste0(files_orig[[j]],"_",i)
    
  }
  
  colnames(metrics_df) <- c("P0_p","P0_O","mean_P","mean_O","l2_P","l2_O","t3_P","t3_O","t4_P","t4_O","ssim","name")
  
  out_list[[j]] <- metrics_df
  
}

out_df <- bind_rows(out_list)
# out_df$model <- "Generator_MSE"
out_df$model <- "WGAN_NZ_DT"

saveRDS(out_df, "./eval/eval_metrics_WGAN_NZ_DT.rds")
# saveRDS(out_df, "./eval/eval_metrics_WGANGonly_6805dayNZ.rds")
# 
#out_df <- readRDS("./eval/eval_metrics_WGANGonly_6805dayNZ.rds")

sum(is.na(out_df$l2_O))


out_df <- readRDS("./eval/eval_metrics_WGAN_NZ_DT.rds")
out_df <- na.omit(out_df)

bias_p0 <- mean(out_df$P0_p - out_df$P0_O)
bias_mean <- mean(out_df$mean_P - out_df$mean_O)
bias_l2 <- mean(out_df$l2_P - out_df$l2_O)
bias_t3 <- mean(out_df$t3_P - out_df$t3_O)
bias_t4 <- mean(out_df$t4_P - out_df$t4_O)

library(hydroGOF)
rmse_p0 <- rmse(out_df$P0_p,out_df$P0_O)
rmse_mean <- rmse(out_df$mean_P,out_df$mean_O)
rmse_l2 <- rmse(out_df$l2_P,out_df$l2_O)
rmse_t3 <- rmse(out_df$t3_P,out_df$t3_O)
rmse_t4 <- rmse(out_df$t4_P,out_df$t4_O)

g1 <- ggplot() + 
  geom_bin2d(data = out_df, aes(x = P0_O, y = P0_p), alpha = 0.8,
             binwidth = c(0.01, 0.01)) +
  geom_abline(slope = 1, linewidth = 0.1, linetype = "dashed") +
  geom_text(aes(label = paste0("RMSE = ", round(rmse_p0, 2))), x = 0, y = 1, size = 1.5, hjust = -0.1) + 
  geom_text(aes(label = paste0("BIAS = ", round(bias_p0, 4))), x = 0, y = 1, size = 1.5, hjust = -0.1, vjust = -0.9) + 
  xlab("Original") +
  ylab("Simulated") + 
  ggtitle("Probability of zero") + 
  scale_x_continuous(limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.25,0.5,0.75,1)) +
  scale_fill_gradientn("",
                       colors = c(
                         "#081d58",
                         "#253494",
                         "#225ea8",
                         "#1d91c0",
                         "#41b6c4",
                         "#7fcdbb",
                         "#c7e9b4",
                         "#edf8b1",
                         "#ffffd9",
                         "#ffeda0",
                         "#fed976",
                         "#feb24c",
                         "#fd8d3c",
                         "#fc4e2a",
                         "#e31a1c",
                         "#bd0026",
                         "#800026")) +
  # alpha = 0.7,
  # trans = scales::pseudo_log_trans(),
  # limits = c(-0.06,0.06),
  # breaks = c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06)) +
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 7, colour = gray(0.25)),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = gray(0.25), size = 0.2),
        axis.text.x = element_text(size = 5.5, color = gray(0.25)),
        axis.text.y = element_text(size = 5.5, color = gray(0.25)),
        plot.title = element_text(hjust = 0, size = 6.5),
        plot.margin=unit(c(5,1,0.5,1), "mm"),
        axis.ticks.length.x = unit(-0.5, "mm"),
        axis.ticks.length.y = unit(-0.5, "mm"),
        axis.ticks = element_line(color = gray(0.25)),
        legend.position='none',
        legend.spacing.x = unit(2, "mm"),
        legend.spacing.y = unit(2, "mm"),
        legend.background = element_rect(color = "transparent", fill = gray(0.95)),
        legend.key.size = unit(3,"mm"),
        legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
        text = element_text(family = "Calibri", color = gray(0.25))) 

g1


g2 <- ggplot() + 
  geom_bin2d(data = out_df, aes(x = mean_O, y = mean_P), alpha = 0.8,
             binwidth = c(0.01, 0.01)) +
  geom_abline(slope = 1, linewidth = 0.1, linetype = "dashed") +
  geom_text(aes(label = paste0("RMSE = ", round(rmse_mean, 2))), x = 0, y = 3.2, size = 1.5, hjust = -0.1) + 
  geom_text(aes(label = paste0("BIAS = ", round(bias_mean, 4))), x = 0, y = 3.2, size = 1.5, hjust = -0.1, vjust = -0.9) + 
  xlab("Original") +
  ylab("Simulated") + 
  ggtitle("Mean") +
  scale_x_continuous(limits = c(0,3.5)) +
  scale_y_continuous(limits = c(0,3.5)) +
  scale_fill_gradientn("",
                       colors = c(
                         "#081d58",
                         "#253494",
                         "#225ea8",
                         "#1d91c0",
                         "#41b6c4",
                         "#7fcdbb",
                         "#c7e9b4",
                         "#edf8b1",
                         "#ffffd9",
                         "#ffeda0",
                         "#fed976",
                         "#feb24c",
                         "#fd8d3c",
                         "#fc4e2a",
                         "#e31a1c",
                         "#bd0026",
                         "#800026")) +
  # alpha = 0.7,
  # trans = scales::pseudo_log_trans(),
  # limits = c(-0.06,0.06),
  # breaks = c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06)) +
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 7, colour = gray(0.25)),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = gray(0.25), size = 0.2),
        axis.text.x = element_text(size = 5.5, color = gray(0.25)),
        axis.text.y = element_text(size = 5.5, color = gray(0.25)),
        plot.title = element_text(hjust = 0, size = 6.5),
        plot.margin=unit(c(5,1,0.5,1), "mm"),
        axis.ticks.length.x = unit(-0.5, "mm"),
        axis.ticks.length.y = unit(-0.5, "mm"),
        axis.ticks = element_line(color = gray(0.25)),
        legend.position='none',
        legend.spacing.x = unit(2, "mm"),
        legend.spacing.y = unit(2, "mm"),
        legend.background = element_rect(color = "transparent", fill = gray(0.95)),
        legend.key.size = unit(3,"mm"),
        legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
        text = element_text(family = "Calibri", color = gray(0.25))) 
g2


g3 <- ggplot() + 
  geom_bin2d(data = out_df, aes(x = l2_O, y = l2_P), alpha = 0.8,
             binwidth = c(0.01, 0.01)) +
  geom_abline(slope = 1, linewidth = 0.1, linetype = "dashed") +
  geom_text(aes(label = paste0("RMSE = ", round(rmse_l2, 2))), x = 0, y = 2.7, size = 1.5, hjust = -0.1) + 
  geom_text(aes(label = paste0("BIAS = ", round(bias_l2, 4))), x = 0, y = 2.7, size = 1.5, hjust = -0.1, vjust = -0.9) + 
  xlab("Original") +
  ylab("Simulated") + 
  ggtitle("Second L-moment") +
  scale_x_continuous(limits = c(0,3)) +
  scale_y_continuous(limits = c(0,3)) +
  scale_fill_gradientn("",
                       colors = c(
                         "#081d58",
                         "#253494",
                         "#225ea8",
                         "#1d91c0",
                         "#41b6c4",
                         "#7fcdbb",
                         "#c7e9b4",
                         "#edf8b1",
                         "#ffffd9",
                         "#ffeda0",
                         "#fed976",
                         "#feb24c",
                         "#fd8d3c",
                         "#fc4e2a",
                         "#e31a1c",
                         "#bd0026",
                         "#800026")) +
  # alpha = 0.7,
  # trans = scales::pseudo_log_trans(),
  # limits = c(-0.06,0.06),
  # breaks = c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06)) +
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 7, colour = gray(0.25)),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = gray(0.25), size = 0.2),
        axis.text.x = element_text(size = 5.5, color = gray(0.25)),
        axis.text.y = element_text(size = 5.5, color = gray(0.25)),
        plot.title = element_text(hjust = 0, size = 6.5),
        plot.margin=unit(c(5,1,0.5,1), "mm"),
        axis.ticks.length.x = unit(-0.5, "mm"),
        axis.ticks.length.y = unit(-0.5, "mm"),
        axis.ticks = element_line(color = gray(0.25)),
        legend.position='none',
        legend.spacing.x = unit(2, "mm"),
        legend.spacing.y = unit(2, "mm"),
        legend.background = element_rect(color = "transparent", fill = gray(0.95)),
        legend.key.size = unit(3,"mm"),
        legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
        text = element_text(family = "Calibri", color = gray(0.25))) 

g3



g4 <- ggplot() + 
  geom_bin2d(data = out_df, aes(x = t3_O, y = t3_P), alpha = 0.8,
             binwidth = c(0.01, 0.01)) +
  geom_abline(slope = 1, linewidth = 0.1, linetype = "dashed") +
  geom_text(aes(label = paste0("RMSE = ", round(rmse_t3, 2))), x = -0.2, y = 0.9, size = 1.5, hjust = -0.1) + 
  geom_text(aes(label = paste0("BIAS = ", round(bias_t3, 4))), x = -0.2, y = 0.9, size = 1.5, hjust = -0.1, vjust = -0.9) + 
  xlab("Original") +
  ylab("Simulated") + 
  ggtitle("L-skewness") +
  scale_x_continuous(limits = c(-0.2,1), breaks = c(-0.2, 0, 0.5,1)) +
  scale_y_continuous(limits = c(-0.2,1), breaks = c(-0.2, 0, 0.5,1)) +
  scale_fill_gradientn("",
                       colors = c(
                         "#081d58",
                         "#253494",
                         "#225ea8",
                         "#1d91c0",
                         "#41b6c4",
                         "#7fcdbb",
                         "#c7e9b4",
                         "#edf8b1",
                         "#ffffd9",
                         "#ffeda0",
                         "#fed976",
                         "#feb24c",
                         "#fd8d3c",
                         "#fc4e2a",
                         "#e31a1c",
                         "#bd0026",
                         "#800026")) +
  # alpha = 0.7,
  # trans = scales::pseudo_log_trans(),
  # limits = c(-0.06,0.06),
  # breaks = c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06)) +
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 7, colour = gray(0.25)),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = gray(0.25), size = 0.2),
        axis.text.x = element_text(size = 5.5, color = gray(0.25)),
        axis.text.y = element_text(size = 5.5, color = gray(0.25)),
        plot.title = element_text(hjust = 0, size = 6.5),
        plot.margin=unit(c(5,1,0.5,1), "mm"),
        axis.ticks.length.x = unit(-0.5, "mm"),
        axis.ticks.length.y = unit(-0.5, "mm"),
        axis.ticks = element_line(color = gray(0.25)),
        legend.position='none',
        legend.spacing.x = unit(2, "mm"),
        legend.spacing.y = unit(2, "mm"),
        legend.background = element_rect(color = "transparent", fill = gray(0.95)),
        legend.key.size = unit(3,"mm"),
        legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
        text = element_text(family = "Calibri", color = gray(0.25))) 
g4






g5 <- ggplot() + 
  geom_bin2d(data = out_df, aes(x = t4_O, y = t4_P), alpha = 0.8,
             binwidth = c(0.01, 0.01)) +
  geom_abline(slope = 1, linewidth = 0.1, linetype = "dashed") +
  geom_text(aes(label = paste0("RMSE = ", round(rmse_t4, 2))), x = -0.2, y = 0.75, size = 1.5, hjust = -0.1) + 
  geom_text(aes(label = paste0("BIAS = ", round(bias_t4, 4))), x = -0.2, y = 0.75, size = 1.5, hjust = -0.1, vjust = -0.9) + 
  xlab("Original") +
  ylab("Simulated") + 
  ggtitle("L-kurtosis") +
  scale_x_continuous(limits = c(-0.2,0.85)) +
  scale_y_continuous(limits = c(-0.2,0.85)) +
  scale_fill_gradientn("",
                       colors = c(
                         "#081d58",
                         "#253494",
                         "#225ea8",
                         "#1d91c0",
                         "#41b6c4",
                         "#7fcdbb",
                         "#c7e9b4",
                         "#edf8b1",
                         "#ffffd9",
                         "#ffeda0",
                         "#fed976",
                         "#feb24c",
                         "#fd8d3c",
                         "#fc4e2a",
                         "#e31a1c",
                         "#bd0026",
                         "#800026")) +
  # alpha = 0.7,
  # trans = scales::pseudo_log_trans(),
  # limits = c(-0.06,0.06),
  # breaks = c(-0.06,-0.04,-0.02,0,0.02,0.04,0.06)) +
  theme(aspect.ratio = 1,
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 7, colour = gray(0.25)),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = gray(0.25), size = 0.2),
        axis.text.x = element_text(size = 5.5, color = gray(0.25)),
        axis.text.y = element_text(size = 5.5, color = gray(0.25)),
        plot.title = element_text(hjust = 0, size = 6.5),
        plot.margin=unit(c(5,1,0.5,1), "mm"),
        axis.ticks.length.x = unit(-0.5, "mm"),
        axis.ticks.length.y = unit(-0.5, "mm"),
        axis.ticks = element_line(color = gray(0.25)),
        legend.position='none',
        legend.spacing.x = unit(2, "mm"),
        legend.spacing.y = unit(2, "mm"),
        legend.background = element_rect(color = "transparent", fill = gray(0.95)),
        legend.key.size = unit(3,"mm"),
        legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
        text = element_text(family = "Calibri", color = gray(0.25))) 

g5

library(ggpubr)
gg1 <- ggarrange(g1,g2,g3,g4,g5, nrow = 1, ncol = 5,
                labels = c("(d) WGAN_NZ_DT","","","",""),
                vjust = 1.2,
                hjust = -0.03,
                font.label = list(size = 7,
                                  color = gray(0.25), 
                                  face = "plain",
                                  base_family = 'Calibri'))
gg1

# ggsave(plot = gg1, "./eval/evalNZ_WGAN_G_400d.jpg",  height = 3.3, width = 16.5, units = "cm", dpi = 400)
ggsave(plot = gg1, "./eval/evalNZ_WGAN_NZ_DT_final.jpg",  height = 3.3, width = 16.5, units = "cm", dpi = 400)

