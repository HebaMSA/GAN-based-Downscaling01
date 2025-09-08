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
##############################################

#Spatial Correlation
# Example 3D array (replace with your actual data)
test_pred <- nctodf("./codes/nc_files/sample_1.nc")
test_orig <- nctodf("./data/Orig/day_monthly_P19570904.nc")

# Extract x, y, and time data
x_coords <- unique(test_pred$x)
y_coords <- unique(test_pred$y)

# Determine spatial dimensions
nx <- length(x_coords)
ny <- length(y_coords)

# Reshape into a 3D array (nx x ny x time)
time_data <- as.matrix(test_orig[, -(1:2)])  # Drop x and y columns
arrayO <- array(time_data, dim = c(nx, ny, ncol(time_data)))

# Initialize lists to store results
horizontal_correlations <- vector("list", dim(arrayO)[3])  # One for each time step
vertical_correlations <- vector("list", dim(arrayO)[3])    # One for each time step

# Function to calculate lagged correlation for a single time slice
calculate_spatial_correlation <- function(slice, max_lagH = 1, max_lagV = 1) {
  nrows <- nrow(slice)
  ncols <- ncol(slice)
  
  # Ensure lags are within valid bounds
  if (max_lagH >= ncols || max_lagV >= nrows) {
    stop("Maximum lags must be smaller than the respective grid dimensions")
  }
  
  # Initialize vectors to store correlations for all lags
  horiz_corrs <- numeric(max_lagH)
  vert_corrs <- numeric(max_lagV)
  
  # Calculate horizontal correlations for each lag
  for (lagH in 1:max_lagH) {
    horizontal <- slice[, 1:(ncols - lagH)]  # Exclude last `lagH` columns
    horizontal_shifted <- slice[, (1 + lagH):ncols]  # Exclude first `lagH` columns
    horiz_corrs[lagH] <- cor(as.vector(horizontal), as.vector(horizontal_shifted), use = "complete.obs")
  }
  
  # Calculate vertical correlations for each lag
  for (lagV in 1:max_lagV) {
    vertical <- slice[1:(nrows - lagV), ]  # Exclude last `lagV` rows
    vertical_shifted <- slice[(1 + lagV):nrows, ]  # Exclude first `lagV` rows
    vert_corrs[lagV] <- cor(as.vector(vertical), as.vector(vertical_shifted), use = "complete.obs")
  }
  
  # Return correlations for all lags
  list(horizontal = horiz_corrs, vertical = vert_corrs)
}

# Loop through each time slice and calculate correlations

###For original files 

files_origALL <- as.list(list.files(pattern = "\\.nc$",
                                    path = "./data/Orig_all/"))

# files_orig <- files_origALL[2801:3200]
files_orig <- files_origALL[13612:20416]

out_listH <- list()
out_listV <- list()

for (j in 1:length(files_orig)){
  
  test_orig <- nctodf(paste0("./data/Orig_all/",files_orig[[j]]))
  x_coords <- unique(test_orig$x)
  y_coords <- unique(test_orig$y)
  
  # Determine spatial dimensions
  nx <- length(x_coords)
  ny <- length(y_coords)
  
  # Reshape into a 3D array (nx x ny x time)
  time_data <- as.matrix(test_orig[, -(1:2)])  # Drop x and y columns
  arrayO <- array(time_data, dim = c(nx, ny, ncol(time_data)))
  
  
  time_listH <- list()
  time_listV <- list()
  for (t in 1:dim(arrayO)[3]) {
    slice <- arrayO[, , t]  # Extract the time slice
    correlations <- calculate_spatial_correlation(slice, max_lagH = 10, max_lagV = 5)
    horizontal_correlations <- data.frame(lags  = c(1:10), cor = correlations$horizontal)
    vertical_correlations <- data.frame(lags  = c(1:5), cor = correlations$vertical)
    time_listH[[t]] <- horizontal_correlations
    time_listV[[t]] <- vertical_correlations
  }
  
  time_dfH <- bind_rows(time_listH)
  time_dfV <- bind_rows(time_listV)
  
  out_listH[[j]] <- time_dfH
  out_listV[[j]] <- time_dfV
}

out_listH_orig <- bind_rows(out_listH)
out_listV_orig <- bind_rows(out_listV)

out_listH_orig$mod <- "Original"
out_listV_orig$mod <- "Original"

saveRDS(out_listH_orig,"./eval/out_listH_orig.rds")
saveRDS(out_listV_orig,"./eval/out_listV_orig.rds")
###For predictions from Generator_MSE files 

# files_pred <- as.list(list.files(pattern = "\\.nc$",
#                                  path = "./codes/nc_filesALL/"))



out_listH <- list()
out_listV <- list()

for (j in 1:length(files_orig)){
  
  test_pred <- nctodf(paste0("./codes/nc_filesALL/sample_",j,".nc"))
  x_coords <- unique(test_pred$x)
  y_coords <- unique(test_pred$y)
  
  # Determine spatial dimensions
  nx <- length(x_coords)
  ny <- length(y_coords)
  
  # Reshape into a 3D array (nx x ny x time)
  time_data <- as.matrix(test_pred[, -(1:2)])  # Drop x and y columns
  arrayO <- array(time_data, dim = c(nx, ny, ncol(time_data)))
  
  
  time_listH <- list()
  time_listV <- list()
  for (t in 1:dim(arrayO)[3]) {
    slice <- arrayO[, , t]  # Extract the time slice
    correlations <- calculate_spatial_correlation(slice, max_lagH = 10, max_lagV = 5)
    horizontal_correlations <- data.frame(lags  = c(1:10), cor = correlations$horizontal)
    vertical_correlations <- data.frame(lags  = c(1:5), cor = correlations$vertical)
    time_listH[[t]] <- horizontal_correlations
    time_listV[[t]] <- vertical_correlations
  }
  
  time_dfH <- bind_rows(time_listH)
  time_dfV <- bind_rows(time_listV)
  
  out_listH[[j]] <- time_dfH
  out_listV[[j]] <- time_dfV
}

out_listH_predG <- bind_rows(out_listH)
out_listV_predG <- bind_rows(out_listV)

out_listH_predG$mod <- "Generator_MSE"
out_listV_predG$mod <- "Generator_MSE"

saveRDS(out_listH_predG,"./eval/out_listH_predG.rds")
saveRDS(out_listV_predG,"./eval/out_listV_predG.rds")

###For predictions from WGAB files 

# files_pred <- as.list(list.files(pattern = "\\.nc$",
#                                  path = "./codes/nc_files_wganG30_ALL/"))



out_listH <- list()
out_listV <- list()

for (j in 1:length(files_orig)){
  
  test_pred <- nctodf(paste0("./codes/nc_files_wganG30_ALL/sample_",j,".nc"))
  
  x_coords <- unique(test_pred$x)
  y_coords <- unique(test_pred$y)
  
  # Determine spatial dimensions
  nx <- length(x_coords)
  ny <- length(y_coords)
  
  # Reshape into a 3D array (nx x ny x time)
  time_data <- as.matrix(test_pred[, -(1:2)])  # Drop x and y columns
  arrayO <- array(time_data, dim = c(nx, ny, ncol(time_data)))
  
  
  time_listH <- list()
  time_listV <- list()
  for (t in 1:dim(arrayO)[3]) {
    slice <- arrayO[, , t]  # Extract the time slice
    correlations <- calculate_spatial_correlation(slice, max_lagH = 10, max_lagV = 5)
    horizontal_correlations <- data.frame(lags  = c(1:10), cor = correlations$horizontal)
    vertical_correlations <- data.frame(lags  = c(1:5), cor = correlations$vertical)
    time_listH[[t]] <- horizontal_correlations
    time_listV[[t]] <- vertical_correlations
  }
  
  time_dfH <- bind_rows(time_listH)
  time_dfV <- bind_rows(time_listV)
  
  out_listH[[j]] <- time_dfH
  out_listV[[j]] <- time_dfV
}

out_listH_predWG <- bind_rows(out_listH)
out_listV_predWG <- bind_rows(out_listV)

out_listH_predWG$mod <- "WGAN"
out_listV_predWG$mod <- "WGAN"

saveRDS(out_listH_predWG,"./eval/out_listH_predWG.rds")
saveRDS(out_listV_predWG,"./eval/out_listV_predWG.rds")
##



out_listH <- list()
out_listV <- list()

for (j in 1:length(files_orig)){
  
  test_pred <- nctodf(paste0("./codes/nc_files_wganG30_ALLNZ/sample_",j,".nc"))
  
  x_coords <- unique(test_pred$x)
  y_coords <- unique(test_pred$y)
  
  # Determine spatial dimensions
  nx <- length(x_coords)
  ny <- length(y_coords)
  
  # Reshape into a 3D array (nx x ny x time)
  time_data <- as.matrix(test_pred[, -(1:2)])  # Drop x and y columns
  arrayO <- array(time_data, dim = c(nx, ny, ncol(time_data)))
  
  
  time_listH <- list()
  time_listV <- list()
  for (t in 1:dim(arrayO)[3]) {
    slice <- arrayO[, , t]  # Extract the time slice
    correlations <- calculate_spatial_correlation(slice, max_lagH = 10, max_lagV = 5)
    horizontal_correlations <- data.frame(lags  = c(1:10), cor = correlations$horizontal)
    vertical_correlations <- data.frame(lags  = c(1:5), cor = correlations$vertical)
    time_listH[[t]] <- horizontal_correlations
    time_listV[[t]] <- vertical_correlations
  }
  
  time_dfH <- bind_rows(time_listH)
  time_dfV <- bind_rows(time_listV)
  
  out_listH[[j]] <- time_dfH
  out_listV[[j]] <- time_dfV
}

out_listH_predWG_NZ <- bind_rows(out_listH)
out_listV_predWG_NZ <- bind_rows(out_listV)

out_listH_predWG_NZ$mod <- "WGAN_NZ"
out_listV_predWG_NZ$mod <- "WGAN_NZ"

saveRDS(out_listH_predWG_NZ,"./eval/out_listH_predWG_NZ.rds")
saveRDS(out_listV_predWG_NZ,"./eval/out_listV_predWG_NZ.rds")

##
out_listH <- list()
out_listV <- list()

for (j in 1:length(files_orig)){
  
  # test_pred <- nctodf(paste0("./codes/nc_files_wganG30_ALL/",files_pred[[j]]))
  # test_pred <- nctodf(paste0("./codes/nc_files_wganG30_hope3/",files_pred[[j]]))
  test_pred <- nctodf(paste0("./codes/nc_files_wganG30_hope3/sample_",j,".nc"))
  
  x_coords <- unique(test_pred$x)
  y_coords <- unique(test_pred$y)
  
  # Determine spatial dimensions
  nx <- length(x_coords)
  ny <- length(y_coords)
  
  # Reshape into a 3D array (nx x ny x time)
  time_data <- as.matrix(test_pred[, -(1:2)])  # Drop x and y columns
  arrayO <- array(time_data, dim = c(nx, ny, ncol(time_data)))
  
  
  time_listH <- list()
  time_listV <- list()
  for (t in 1:dim(arrayO)[3]) {
    slice <- arrayO[, , t]  # Extract the time slice
    correlations <- calculate_spatial_correlation(slice, max_lagH = 10, max_lagV = 5)
    horizontal_correlations <- data.frame(lags  = c(1:10), cor = correlations$horizontal)
    vertical_correlations <- data.frame(lags  = c(1:5), cor = correlations$vertical)
    time_listH[[t]] <- horizontal_correlations
    time_listV[[t]] <- vertical_correlations
  }
  
  time_dfH <- bind_rows(time_listH)
  time_dfV <- bind_rows(time_listV)
  
  out_listH[[j]] <- time_dfH
  out_listV[[j]] <- time_dfV
}

out_listH_predWG_NZ_DT <- bind_rows(out_listH)
out_listV_predWG_NZ_DT <- bind_rows(out_listV)

out_listH_predWG_NZ_DT$mod <- "WGAN_NZ_DT"
out_listV_predWG_NZ_DT$mod <- "WGAN_NZ_DT"

saveRDS(out_listH_predWG_NZ_DT,"./eval/out_listH_predWG_NZ_DT.rds")
saveRDS(out_listV_predWG_NZ_DT,"./eval/out_listV_predWG_NZ_DT.rds")


out_H_all <- rbind(out_listH_orig, out_listH_predG, out_listH_predWG, out_listH_predWG_NZ, out_listH_predWG_NZ_DT)
out_V_all <- rbind(out_listV_orig, out_listV_predG, out_listV_predWG, out_listV_predWG_NZ, out_listV_predWG_NZ_DT)

out_H_all$mod <- factor(out_H_all$mod, levels = c("Original", "Generator_MSE", "WGAN","WGAN_NZ","WGAN_NZ_DT"))
out_V_all$mod <- factor(out_V_all$mod, levels = c("Original", "Generator_MSE", "WGAN","WGAN_NZ","WGAN_NZ_DT"))

saveRDS(out_H_all,"./eval/SCor_H_6805_ALL.rds")
saveRDS(out_V_all,"./eval/SCor_V_6805_ALL.rds")


#####
out_H_all <- readRDS("./eval/SCor_H_6805_ALL.rds")
out_V_all <- readRDS("./eval/SCor_V_6805_ALL.rds")


out_statsH <- out_H_all %>% 
  group_by(lags, mod) %>% 
  summarise(min = quantile(cor, 0.05, na.rm = T),
            max = quantile(cor, 0.95, na.rm = T),
            med = median(cor, 0.5, na.rm = T),
            mean = mean(cor, na.rm = T))

out_statsH$mod <- as.character(out_statsH$mod)
out_statsH$mod[out_statsH$mod == "Generator_MSE"] <- "UNET"

out_statsH$mod <- factor(out_statsH$mod, levels = c("Original", "UNET", "WGAN", "WGAN_NZ","WGAN_NZ_DT"))


out_statsV <- out_V_all %>% 
  group_by(lags, mod) %>% 
  summarise(min = quantile(cor, 0.05, na.rm = T),
            max = quantile(cor, 0.95, na.rm = T),
            med = median(cor, 0.5, na.rm = T),
            mean = mean(cor, na.rm = T))


out_statsV$mod <- as.character(out_statsV$mod)
out_statsV$mod[out_statsV$mod == "Generator_MSE"] <- "UNET"
out_statsV$mod <- factor(out_statsV$mod, levels = c("Original", "UNET", "WGAN", "WGAN_NZ","WGAN_NZ_DT"))

g1 <- ggplot() + 
  geom_point(
    data = out_statsH, 
    aes(
      x = as.factor(lags), 
      y = mean, 
      group = interaction(as.factor(lags), mod),  # Group by both lags and mod
      color = mod  # Color by mod
    ), position = position_dodge2(0.5), size = 0.6) + 
  geom_linerange(data = out_statsH, 
                 aes(
                   x = as.factor(lags), 
                   ymin = min,
                   ymax = max,
                   group = interaction(as.factor(lags), mod),  # Group by both lags and mod
                   color = mod  # Color by mod
                 ), position = position_dodge2(0.5), linewidth = 0.3) +
  scale_color_manual("", values =  c("Original" =  "#c15ca5",
                                     "UNET" = "#60a862",
                                     "WGAN" = "#727cce",
                                     "WGAN_NZ" = "#b4943e",
                                     "WGAN_NZ_DT" ="#cb5a4c"),
                     breaks = c("Original" ,
                                "UNET" ,
                                "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  scale_x_discrete(labels = c("50","100","150","200","250","300","350","400","450","500")) + 
  xlab("Horizontal distance lag (km)") + 
  ylab("Correlation") +
  theme(#aspect.ratio = 0.3,
    legend.text = element_text(size = 5),
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
    legend.position=c(0.6,0.98),
    legend.spacing.x = unit(0.01, "mm"),
    legend.spacing.y = unit(0.01, "mm"),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    # legend.background = element_rect(color = "transparent", fill = gray(0.95)),
    legend.key.size = unit(2,"mm"),
    legend.key.width = unit(2,"mm"),
    # legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
    text = element_text(family = "Calibri", color = gray(0.25)),
    legend.box.margin=margin(0.5,0.5,0.5,0.5)) +
  guides(color = guide_legend(nrow = 1, ncol = 5))

g1
ggsave("./eval/SC_H_cor_hope3.jpg",  height = 12, width = 15, units = "cm", dpi = 400)

g11 <- ggplot() + 
  geom_point(
    data = out_statsV, 
    aes(
      x = as.factor(lags), 
      y = mean, 
      group = interaction(as.factor(lags), mod),  # Group by both lags and mod
      color = mod  # Color by mod
    ), position = position_dodge2(0.5), size = 0.5) + 
  geom_linerange(data = out_statsV, 
                 aes(
                   x = as.factor(lags), 
                   ymin = min,
                   ymax = max,
                   group = interaction(as.factor(lags), mod),  # Group by both lags and mod
                   color = mod  # Color by mod
                 ), position = position_dodge2(0.5), linewidth = 0.3) +
  scale_color_manual("", values =  c("Original" =  "#c15ca5",
                                     "UNET" = "#60a862",
                                     "WGAN" = "#727cce",
                                     "WGAN_NZ" = "#b4943e",
                                     "WGAN_NZ_DT" ="#cb5a4c"),
                     breaks = c("Original" ,
                                "UNET" ,
                                "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  scale_x_discrete(labels = c("50","100","150","200","250","300","350","400","450","500")) + 
  scale_y_continuous(limits = c(-0.2,1)) + 
  xlab("Vertical distance lag (km)") + 
  ylab("Correlation") +
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

g11
ggsave("./eval/SC_V_cor_hope3.jpg",  height = 12, width = 15, units = "cm", dpi = 400)

out_listH_orig <- readRDS("./eval/out_listH_orig.rds")
out_listH_predG <- readRDS("./eval/out_listH_predG.rds")
out_listH_predG$mod <- "UNET"
out_listH_predWG <- readRDS("./eval/out_listH_predWG.rds")
out_listH_predWG_NZ <- readRDS("./eval/out_listH_predWG_NZ.rds")
out_listH_predWG_NZ_DT <- readRDS("./eval/out_listH_predWG_NZ_DT.rds")

H_wide <- cbind(out_listH_orig, out_listH_predG,out_listH_predWG,out_listH_predWG_NZ, out_listH_predWG_NZ_DT)
colnames(H_wide) <- c("lagsO","corO","modO",
                      "lagsG","corG","modG",
                      "lagsWG","corWG","modWG",
                      "lagsWG_NZ","corWG_NZ","modWG_NZ",
                      "lagsWG_NZ_DT","corWG_NZ_DT","modWG_NZ_DT")

H_wide_nNA <- na.omit(H_wide)

df_evalG <- data.frame()

for (i in c(1:10)){
  lags1 <- H_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corG - lags1$corO)
  rmse_1 <- rmse(lags1$corG, lags1$corO)
  
  df_evalG[i,1] <- i
  df_evalG[i,2] <- bias_1
  df_evalG[i,3] <- rmse_1
}

colnames(df_evalG) <- c("lags","bias","rmse")
df_evalG$mod <- "UNET"




df_evalWG <- data.frame()

for (i in c(1:10)){
  lags1 <- H_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corWG - lags1$corO)
  rmse_1 <- rmse(lags1$corWG, lags1$corO)
  
  df_evalWG[i,1] <- i
  df_evalWG[i,2] <- bias_1
  df_evalWG[i,3] <- rmse_1
}

colnames(df_evalWG) <- c("lags","bias","rmse")
df_evalWG$mod <- "WGAN"


df_evalWG_NZ <- data.frame()

for (i in c(1:10)){
  lags1 <- H_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corWG_NZ - lags1$corO)
  rmse_1 <- rmse(lags1$corWG_NZ, lags1$corO)
  
  df_evalWG_NZ[i,1] <- i
  df_evalWG_NZ[i,2] <- bias_1
  df_evalWG_NZ[i,3] <- rmse_1
}

colnames(df_evalWG_NZ) <- c("lags","bias","rmse")
df_evalWG_NZ$mod <- "WGAN_NZ"


df_evalWG_NZ_DT <- data.frame()

for (i in c(1:10)){
  lags1 <- H_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corWG_NZ_DT - lags1$corO)
  rmse_1 <- rmse(lags1$corWG_NZ_DT, lags1$corO)
  
  df_evalWG_NZ_DT[i,1] <- i
  df_evalWG_NZ_DT[i,2] <- bias_1
  df_evalWG_NZ_DT[i,3] <- rmse_1
}

colnames(df_evalWG_NZ_DT) <- c("lags","bias","rmse")
df_evalWG_NZ_DT$mod <- "WGAN_NZ_DT"

df_eval_H <- rbind(df_evalG, df_evalWG, df_evalWG_NZ, df_evalWG_NZ_DT)

saveRDS(df_eval_H,"df_eval_H_all.rds")


g2 <- ggplot() + 
  geom_line(data = df_eval_H, aes(as.factor(lags), bias, group = mod, color = mod), linewidth = 0.3)+
  geom_point(data = df_eval_H, aes(as.factor(lags), bias, group = mod, color = mod), size = 0.5)+
  scale_color_manual("", values =  c("Generator_MSE" = "#60a862",
                                     "WGAN" = "#727cce",
                                     "WGAN_NZ" = "#b4943e",
                                     "WGAN_NZ_DT" ="#cb5a4c"),
                     breaks = c("Generator_MSE" ,
                                "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  scale_x_discrete(labels = c("50","100","150","200","250","300","350","400","450","500")) + 
  xlab("Horizontal distance lag (km)") + 
  ylab("Bias") +
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
  geom_line(data = df_eval_H, aes(as.factor(lags), rmse, group = mod, color = mod), linewidth  = 0.3)+
  geom_point(data = df_eval_H, aes(as.factor(lags), rmse, group = mod, color = mod), size = 0.6)+
  scale_color_manual("", values =  c("Generator_MSE" = "#60a862",
                                     "WGAN" = "#727cce",
                                     "WGAN_NZ" = "#b4943e",
                                     "WGAN_NZ_DT" ="#cb5a4c"),
                     breaks = c("Generator_MSE" ,
                                "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  scale_x_discrete(labels = c("50","100","150","200","250","300","350","400","450","500")) + 
  xlab("Horizontal distance lag (km)") + 
  ylab("RMSE") +
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




out_listV_orig <- readRDS("./eval/out_listV_orig.rds")
out_listV_predG <- readRDS("./eval/out_listV_predG.rds")
out_listV_predG$mod <- "UNET"
out_listV_predWG <- readRDS("./eval/out_listV_predWG.rds")
out_listV_predWG_NZ <- readRDS("./eval/out_listV_predWG_NZ.rds")
out_listV_predWG_NZ_DT <- readRDS("./eval/out_listV_predWG_NZ_DT.rds")


V_wide <- cbind(out_listV_orig, out_listV_predG,out_listV_predWG,out_listV_predWG_NZ,out_listV_predWG_NZ_DT)
colnames(V_wide) <- c("lagsO","corO","modO","lagsG","corG","modG","lagsWG","corWG","modWG",
                      "lagsWG_NZ","corWG_NZ","modWG_NZ",
                      "lagsWG_NZ_DT","corWG_NZ_DT","modWG_NZ_DT")

V_wide_nNA <- na.omit(V_wide)

df_evalG <- data.frame()

for (i in c(1:5)){
  lags1 <- V_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corG - lags1$corO)
  rmse_1 <- rmse(lags1$corG, lags1$corO)
  
  df_evalG[i,1] <- i
  df_evalG[i,2] <- bias_1
  df_evalG[i,3] <- rmse_1
}

colnames(df_evalG) <- c("lags","bias","rmse")
df_evalG$mod <- "UNET"




df_evalWG <- data.frame()

for (i in c(1:5)){
  lags1 <- V_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corWG - lags1$corO)
  rmse_1 <- rmse(lags1$corWG, lags1$corO)
  
  df_evalWG[i,1] <- i
  df_evalWG[i,2] <- bias_1
  df_evalWG[i,3] <- rmse_1
}

colnames(df_evalWG) <- c("lags","bias","rmse")
df_evalWG$mod <- "WGAN"


df_evalWG_NZ <- data.frame()

for (i in c(1:5)){
  lags1 <- V_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corWG_NZ - lags1$corO)
  rmse_1 <- rmse(lags1$corWG_NZ, lags1$corO)
  
  df_evalWG_NZ[i,1] <- i
  df_evalWG_NZ[i,2] <- bias_1
  df_evalWG_NZ[i,3] <- rmse_1
}

colnames(df_evalWG_NZ) <- c("lags","bias","rmse")
df_evalWG_NZ$mod <- "WGAN_NZ"


df_evalWG_NZ_DT <- data.frame()

for (i in c(1:5)){
  lags1 <- V_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corWG_NZ_DT - lags1$corO)
  rmse_1 <- rmse(lags1$corWG_NZ_DT, lags1$corO)
  
  df_evalWG_NZ_DT[i,1] <- i
  df_evalWG_NZ_DT[i,2] <- bias_1
  df_evalWG_NZ_DT[i,3] <- rmse_1
}

colnames(df_evalWG_NZ_DT) <- c("lags","bias","rmse")
df_evalWG_NZ_DT$mod <- "WGAN_NZ_DT"


df_eval_V <- rbind(df_evalG, df_evalWG, df_evalWG_NZ, df_evalWG_NZ_DT)

saveRDS(df_eval_V,"df_eval_V_all.rds")


g22 <- ggplot() + 
  geom_line(data = df_eval_V, aes(as.factor(lags), bias, group = mod, color = mod), linewidth = 0.3)+
  geom_point(data = df_eval_V, aes(as.factor(lags), bias, group = mod, color = mod), size = 0.6)+
  scale_color_manual("", values =  c("Generator_MSE" = "#60a862",
                                     "WGAN" = "#727cce",
                                     "WGAN_NZ" = "#b4943e",
                                     "WGAN_NZ_DT" ="#cb5a4c"),
                     breaks = c("Generator_MSE" ,
                                "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  scale_x_discrete(labels = c("50","100","150","200","250","300","350","400","450","500")) + 
  
  xlab("Vertical distance lag (km)") + 
  ylab("Bias") +
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

g22




g33 <- ggplot() + 
  geom_line(data = df_eval_V, aes(as.factor(lags), rmse, group = mod, color = mod), linewidth = 0.3)+
  geom_point(data = df_eval_V, aes(as.factor(lags), rmse, group = mod, color = mod), size = 0.6)+
  scale_color_manual("", values =  c("Generator_MSE" = "#60a862",
                                     "WGAN" = "#727cce",
                                     "WGAN_NZ" = "#b4943e",
                                     "WGAN_NZ_DT" ="#cb5a4c"),
                     breaks = c("Generator_MSE" ,
                                "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  scale_x_discrete(labels = c("50","100","150","200","250","300","350","400","450","500")) + 
  
  xlab("Vertical distance lag (km)") + 
  ylab("RMSE") +
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

g33


library(ggpubr)
ggarrange(g1, g2, g3,
          g11, g22, g33,
          nrow = 2, 
          ncol = 3,
          labels = c("(a)",
                     "(b)",
                     "(c)",
                     "(d)",
                     "(e)",
                     "(f)"),
          common.legend = TRUE,
          legend = "bottom",
          align = "hv",
          vjust = c(1.5,1.5,1.5,1.5,1.5,1.5),
          hjust = c(-3.2,-3.2,-3.2,-3.2,-3.2,-4),
          widths = c(1,0.7, 0.7),
          font.label = list(size = 8,
                            color = gray(0.25), 
                            face = "plain",
                            base_family = 'Calibri'))

ggsave("./eval/Cor_spatial_6805all_ALL_namescor.jpg",  height = 8, width = 16.5, units = "cm", dpi = 400)
