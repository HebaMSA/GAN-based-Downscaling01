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




calculate_temporal_correlation_nonzero <- function(data_array, max_lag = 1) {
  # Get dimensions of the 3D array
  nrows <- dim(data_array)[1]
  ncols <- dim(data_array)[2]
  ntime <- dim(data_array)[3]
  
  # Ensure max_lag is within valid bounds
  if (max_lag >= ntime) {
    stop("max_lag must be smaller than the number of time steps")
  }
  
  # Initialize a 4D array to store temporal correlations for each grid point
  # Dimensions: (nrows x ncols x max_lag)
  temporal_corrs <- array(NA, dim = c(nrows, ncols, max_lag))
  
  # Loop through each lag
  for (lag in 1:max_lag) {
    # Loop through each grid point
    for (i in 1:nrows) {
      for (j in 1:ncols) {
        # Extract the time series for the grid point
        ts1 <- data_array[i, j, 1:(ntime - lag)]        # Original time series
        ts2 <- data_array[i, j, (1 + lag):ntime]        # Lagged time series
        
        # Filter out zero values from both time series
        nonzero_indices <- which(ts1 > 0 & ts2 > 0)
        ts1_nonzero <- ts1[nonzero_indices]
        ts2_nonzero <- ts2[nonzero_indices]
        
        # Calculate correlation only if there are enough nonzero values
        if (length(ts1_nonzero) > 1 && length(ts2_nonzero) > 1) {
          temporal_corrs[i, j, lag] <- cor(ts1_nonzero, ts2_nonzero, use = "complete.obs")
        } else {
          temporal_corrs[i, j, lag] <- NA  # Not enough nonzero values for correlation
        }
      }
    }
  }
  
  return(temporal_corrs)
}


calculate_temporal_correlation <- function(data_array, max_lag = 1) {
  # Get dimensions of the 3D array
  nrows <- dim(data_array)[1]
  ncols <- dim(data_array)[2]
  ntime <- dim(data_array)[3]
  
  # Ensure max_lag is within valid bounds
  if (max_lag >= ntime) {
    stop("max_lag must be smaller than the number of time steps")
  }
  
  # Initialize a 4D array to store temporal correlations for each grid point
  # Dimensions: (nrows x ncols x max_lag)
  temporal_corrs <- array(NA, dim = c(nrows, ncols, max_lag))
  
  # Loop through each lag
  for (lag in 1:max_lag) {
    # Loop through each grid point
    for (i in 1:nrows) {
      for (j in 1:ncols) {
        # Extract the time series for the grid point
        ts1 <- data_array[i, j, 1:(ntime - lag)]        # Original time series
        ts2 <- data_array[i, j, (1 + lag):ntime]        # Lagged time series
        
        # Calculate correlation and store it
        temporal_corrs[i, j, lag] <- cor(ts1, ts2, use = "complete.obs")
      }
    }
  }
  
  return(temporal_corrs)
}


library(reshape2)

# Assuming temporal_corrs is the 4D array with dimensions (nrows, ncols, max_lag)
# Convert the array to a long-format dataframe
convert_to_dataframe <- function(temporal_corrs) {
  nrows <- dim(temporal_corrs)[1]
  ncols <- dim(temporal_corrs)[2]
  max_lag <- dim(temporal_corrs)[3]
  
  # Melt the array into a long-format dataframe
  df <- melt(temporal_corrs, varnames = c("row", "col", "lag"), value.name = "cor")
  
  # Add a unique grid ID for each grid cell (row, col combination)
  df$grid_id <- (df$row - 1) * ncols + df$col  # Unique ID for each grid
  
  # Convert lag to numeric
  df$lag <- as.numeric(df$lag)
  
  # Reorder columns
  df <- df[, c("grid_id", "lag", "cor")]
  
  return(df)
}


testC <- calculate_temporal_correlation_nonzero(arrayO, max_lag = 5)

testC_df <- convert_to_dataframe(testC)






###For original files 

files_origALL <- as.list(list.files(pattern = "\\.nc$",
                                    path = "./data/Orig_all/"))

# files_orig <- files_origALL[2801:3200]
files_orig <- files_origALL[13612:20416]

out_list <- list()

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
  
  Cor_arr <- calculate_temporal_correlation_nonzero(arrayO, max_lag = 5)
  Cor_df <- convert_to_dataframe(Cor_arr)
  
  
  out_list[[j]] <- Cor_df
  
}

out_list_orig <- bind_rows(out_list)

out_list_orig$mod <- "Original"

saveRDS(out_list_orig,"./eval/Tcor_orig_NZ.rds")

out_list_orig <- readRDS("./eval/Tcor_orig_NZ.rds")

remove(out_list)

###For predictions from Generator_MSE files 

out_list <- list()

for (j in 1:length(files_orig)){
  
  # test_orig <- nctodf(paste0("./codes/nc_filesALL/sample_",j,".nc"))
  test_orig <- nctodf(paste0("./codes/nc_files_hope3/sample_",j,".nc"))
  x_coords <- unique(test_orig$x)
  y_coords <- unique(test_orig$y)
  
  # Determine spatial dimensions
  nx <- length(x_coords)
  ny <- length(y_coords)
  
  # Reshape into a 3D array (nx x ny x time)
  time_data <- as.matrix(test_orig[, -(1:2)])  # Drop x and y columns
  arrayO <- array(time_data, dim = c(nx, ny, ncol(time_data)))
  
  Cor_arr <- calculate_temporal_correlation_nonzero(arrayO, max_lag = 5)
  Cor_df <- convert_to_dataframe(Cor_arr)
  
  
  out_list[[j]] <- Cor_df
  
}

out_list_predG <- bind_rows(out_list)

out_list_predG$mod <- "Generator_MSE"


saveRDS(out_list_predG,"./eval/Tcor_predG_NZ_h3.rds")
out_list_predG <- readRDS("./eval/Tcor_predG_NZ_h3.rds")


###For predictions from Generator_MSE files 
# 
# files_pred <- as.list(list.files(pattern = "\\.nc$",
#                                  path = "./codes/nc_files_wganG/"))

out_list <- list()

for (j in 1:length(files_orig)){
  
  
  test_orig <- nctodf(paste0("./codes/nc_files_wganG30_ALLNZ/sample_",j,".nc"))
  # test_orig <- nctodf(paste0("./codes/nc_files_wganG30_hope3/sample_",j,".nc"))
  x_coords <- unique(test_orig$x)
  y_coords <- unique(test_orig$y)
  
  # Determine spatial dimensions
  nx <- length(x_coords)
  ny <- length(y_coords)
  
  # Reshape into a 3D array (nx x ny x time)
  time_data <- as.matrix(test_orig[, -(1:2)])  # Drop x and y columns
  arrayO <- array(time_data, dim = c(nx, ny, ncol(time_data)))
  
  Cor_arr <- calculate_temporal_correlation_nonzero(arrayO, max_lag = 5)
  Cor_df <- convert_to_dataframe(Cor_arr)
  
  
  out_list[[j]] <- Cor_df
  
}

out_list_predWG <- bind_rows(out_list)

out_list_predWG$mod <- "WGAN"
saveRDS(out_list_predWG,"./eval/Tcor_predWG_NZ_h2.rds")


##

out_list_orig <- readRDS("./eval/Tcor_orig_NZ.rds")
out_list_predG <- readRDS("./eval/Tcor_predG_NZ_h3.rds")
out_list_predG$mod <- "UNET"

out_list_predWG <- readRDS("./eval/Tcor_predWG_NZ.rds")
out_list_predWG_NZ <- readRDS("./eval/Tcor_predWG_NZ_h2.rds")
out_list_predWG_NZ$mod <- "WGAN_NZ"
out_list_predWG_NZ_DT <- readRDS("./eval/Tcor_predWG_NZ_h3.rds")
out_list_predWG_NZ_DT$mod <- "WGAN_NZ_DT"

out_all <- rbind(out_list_orig, out_list_predG, out_list_predWG, out_list_predWG_NZ, out_list_predWG_NZ_DT)

out_stats <- out_all %>% 
  group_by(lag, mod) %>% 
  summarise(min = quantile(cor, 0.05, na.rm = T),
            max = quantile(cor, 0.95, na.rm = T),
            med = median(cor, 0.5, na.rm = T),
            mean = mean(cor, na.rm = T))


out_stats$mod <- factor(out_stats$mod, levels = c("Original", "UNET", "WGAN", "WGAN_NZ","WGAN_NZ_DT"))


saveRDS(out_all,"./eval/Tcor_all_namescor.rds")

g4 <- ggplot() + 
  geom_point(
    data = out_stats, 
    aes(
      x = as.factor(lag), 
      y = mean, 
      group = interaction(as.factor(lag), mod),  # Group by both lags and mod
      color = mod  # Color by mod
    ), position = position_dodge2(0.5), size = 0.6) + 
  geom_linerange(data = out_stats, 
                aes(
                  x = as.factor(lag), 
                  ymin = min,
                  ymax = max,
                  group = interaction(as.factor(lag), mod),  # Group by both lags and mod
                  color = mod  # Color by mod
                ), position = position_dodge2(0.5), linewidth = 0.3) +
  scale_color_manual("", values =  c("Original" =  "#c15ca5",
                                    "Generator_MSE" = "#60a862",
                                    "WGAN" = "#727cce",
                                    "WGAN_NZ" = "#b4943e",
                                    "WGAN_NZ_DT" ="#cb5a4c"),
                    breaks = c("Original" ,
                               "Generator_MSE" ,
                               "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  xlab("Temporal lag (hr)") + 
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
    legend.position=c(0.6,1),
    legend.spacing.x = unit(0.01, "mm"),
    legend.spacing.y = unit(0.01, "mm"),
    legend.box.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    # legend.background = element_rect(color = "transparent", fill = gray(0.95)),
    legend.key.size = unit(3,"mm"),
    legend.key.width = unit(3,"mm"),
    # legend.key = element_rect(fill = gray(0.95), colour = gray(0.95)),
    text = element_text(family = "Calibri", color = gray(0.25)),
    legend.box.margin=margin(0.5,0.5,0.5,0.5)) +
  guides(color = guide_legend(nrow=1, ncol = 5))

g4



all_wide <- cbind(out_list_orig, out_list_predG,out_list_predWG, out_list_predWG_NZ, out_list_predWG_NZ_DT)
colnames(all_wide) <- c("idO","lagsO","corO","modO",
                        "idG", "lagsG","corG","modG",
                        "idWG", "lagsWG","corWG","modWG",
                        "idWG_NZ", "lagsWG_NZ","corWG_NZ","modWG_NZ",
                        "idWG_NZDT", "lagsWG_NZDT","corWG_NZDT","modWG_NZDT")
remove(out_list_orig, out_list_predG,out_list_predWG, out_list_predWG_NZ, out_list_predWG_NZ_DT)

all_wide_nNA <- na.omit(all_wide)

df_evalG <- data.frame()

for (i in c(1:5)){
  lags1 <- all_wide_nNA %>% 
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
  lags1 <- all_wide_nNA %>% 
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
  lags1 <- all_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corWG_NZ - lags1$corO)
  rmse_1 <- rmse(lags1$corWG_NZ, lags1$corO)
  
  df_evalWG_NZ[i,1] <- i
  df_evalWG_NZ[i,2] <- bias_1
  df_evalWG_NZ[i,3] <- rmse_1
}

colnames(df_evalWG_NZ) <- c("lags","bias","rmse")
df_evalWG_NZ$mod <- "WGAN_NZ"


df_evalWG_NZDT <- data.frame()

for (i in c(1:5)){
  lags1 <- all_wide_nNA %>% 
    filter(lagsO == i) 
  
  bias_1 <- mean(lags1$corWG_NZDT - lags1$corO)
  rmse_1 <- rmse(lags1$corWG_NZDT, lags1$corO)
  
  df_evalWG_NZDT[i,1] <- i
  df_evalWG_NZDT[i,2] <- bias_1
  df_evalWG_NZDT[i,3] <- rmse_1
}

colnames(df_evalWG_NZDT) <- c("lags","bias","rmse")
df_evalWG_NZDT$mod <- "WGAN_NZ_DT"

df_eval <- rbind(df_evalG, df_evalWG, df_evalWG_NZ, df_evalWG_NZDT)



remove(df_evalG, df_evalWG, df_evalWG_NZ, df_evalWG_NZDT)


g5 <- ggplot() + 
  geom_line(data = df_eval, aes(as.factor(lags), bias, group = mod, color = mod), linewidth = 0.3)+
  geom_point(data = df_eval, aes(as.factor(lags), bias, group = mod, color = mod), size  = 0.6)+
  # scale_color_manual("", values =  c("Generator_MSE" = "#50b085",
  #                                    "WGAN" = "#74893d"),
  #                    breaks = c("Generator_MSE" ,
  #                               "WGAN" )) +
  scale_color_manual("", values =  c("Generator_MSE" = "#60a862",
                                    "WGAN" = "#727cce",
                                    "WGAN_NZ" = "#b4943e",
                                    "WGAN_NZ_DT" ="#cb5a4c"),
                    breaks = c("Generator_MSE" ,
                               "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  xlab("Temporal lag (hr)") + 
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

g5





g6 <- ggplot() + 
  geom_line(data = df_eval, aes(as.factor(lags), rmse, group = mod, color = mod), linewidth = 0.3)+
  geom_point(data = df_eval, aes(as.factor(lags), rmse, group = mod, color = mod), size = 0.6)+
  # scale_color_manual("", values =  c("Generator_MSE" = "#50b085",
  #                                    "WGAN" = "#74893d"),
  #                    breaks = c("Generator_MSE" ,
  #                               "WGAN" )) +
  scale_color_manual("", values =  c("Generator_MSE" = "#60a862",
                                    "WGAN" = "#727cce",
                                    "WGAN_NZ" = "#b4943e",
                                    "WGAN_NZ_DT" ="#cb5a4c"),
                    breaks = c("Generator_MSE" ,
                               "WGAN",  "WGAN_NZ", "WGAN_NZ_DT")) +
  xlab("Temporal lag (hr)") + 
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

g6

library(ggpubr)
gg <- ggarrange(g4, g5, g6,
          nrow = 1, 
          ncol = 3,
          labels = c("(g)",
                     "(h)",
                     "(i)"),
          align = "hv",
          common.legend = TRUE, 
          legend = "bottom",
          vjust = c(1.5,1.5,1.5),
          hjust = c(-3.2,-3.2,-4.1),
          widths = c(1,0.7, 0.7),
          font.label = list(size = 8,
                            color = gray(0.25), 
                            face = "plain",
                            base_family = 'Calibri'))
gg
ggsave(plot = gg, "./eval/Cor_temporal_NZ_ALL.jpg",  height = 5, width = 16.5, units = "cm", dpi = 400)






library(ggpubr)
ggarrange(g1, g2, g3,
          g11, g22, g33,
          g4, g5, g6,
          nrow = 3, 
          ncol = 3,
          labels = c("(a)",
                     "(b)",
                     "(c)",
                     "(d)",
                     "(e)",
                     "(f)",
                     "(g)",
                     "(h)",
                     "(i)"),
          align = "hv",
          vjust = c(1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5),
          hjust = c(-3.2,-3.2,-3.2,-3.2,-3.2,-4,-3.2,-3.2,-4),
          widths = c(1,0.7, 0.7),
          font.label = list(size = 8,
                            color = gray(0.25), 
                            face = "plain",
                            base_family = 'Calibri'))

# ggsave("Cor_all.jpg",  height = 12, width = 16.5, units = "cm", dpi = 400)


