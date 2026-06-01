#Author: Nayan Talmale
#This code helps in computing the Gof and BIC of the copula model

library(parallel)
library(foreach)
library(doParallel)

library(VineCopula)
library(ncdf4)
library(ggplot2)
library(dplyr)
library(maps)

set.seed(271)

# Load data
Tmax_file <- "/g/data/gb02/nt4273/Copula_Data/Detrend_ERA5_Tmax_Jul.nc"
Precip_file <- "/g/data/gb02/nt4273/Copula_Data/Detrend_ERA5_P_Jul.nc"

Tmax_nc <- nc_open(Tmax_file)
Precip_nc <- nc_open(Precip_file)

lon <- ncvar_get(Tmax_nc, "longitude")
lat <- ncvar_get(Tmax_nc, "latitude")
T <- ncvar_get(Tmax_nc, "Tmax_ERA5_Jul")
P <- ncvar_get(Precip_nc, "P_ERA5_Jul")

# Initialize output matrices
best_copula <- matrix(NA, nrow = length(lon), ncol = length(lat))
LMF_map <- matrix(NA, nrow = length(lon), ncol = length(lat))

# Define copula families and names
copula_families <- c(4, 3, 5, 6, 24, 23, 26)
copula_names <- c("Gumbel", "Clayton", "Frank", "Joe",
                  "Gumbel90", "Clayton90", "Joe90")


#Parallel
ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)

# -------------------------------
# Flatten grid indices (each task = one grid point)
# -------------------------------
grid_idx <- expand.grid(i = 1:length(lon), j = 1:length(lat))

# -------------------------------
# Run parallel loop
# -------------------------------
results <- foreach(idx = 1:nrow(grid_idx), .combine = rbind,
                   .packages = "VineCopula") %dopar% {
                     i <- grid_idx$i[idx]
                     j <- grid_idx$j[idx]
                     
                     u <- pobs(T[i, j, ])
                     v <- pobs(P[i, j, ])
                     
                     # Goodness-of-fit tests
                     gof_list <- lapply(copula_families, function(fam) {
                       BiCopGofTest(u, v, family = fam, method = "kendall")
                     })
                     p_values <- sapply(gof_list, function(g) g$p.value.CvM)
                     valid_idx <- which(p_values > 0.05)
                     
                     if (length(valid_idx) > 0) {
                       # Fit valid copulas
                       fits_list <- lapply(copula_families[valid_idx], function(fam) {
                         BiCopEst(u, v, family = fam, method = "mle")
                       })
                       BIC_values <- sapply(fits_list, function(fit) fit$BIC)
                       best_idx <- which.min(BIC_values)
                       best_fit <- fits_list[[best_idx]]
                       best_name <- copula_names[valid_idx][best_idx]
                       
                       # Compute LMF
                       u_thresh <- 0.75
                       v_thresh <- 0.25
                       Pr <- v_thresh - BiCopCDF(u_thresh, v_thresh,
                                                 family = best_fit$family,
                                                 par    = best_fit$par)
                       LMF <- Pr / (0.25 * 0.25)
                     } else {
                       best_name <- NA
                       LMF <- NA
                     }
                     
                     c(i, j, best_name, LMF)
                   }

# -------------------------------
# Stop parallel cluster
# -------------------------------
stopCluster(cl)

# -------------------------------
# Convert results into matrices
# -------------------------------
best_copula <- matrix(NA, nrow = length(lon), ncol = length(lat))
LMF_map     <- matrix(NA, nrow = length(lon), ncol = length(lat))

for (r in 1:nrow(results)) {
  i <- as.numeric(results[r, 1])
  j <- as.numeric(results[r, 2])
  best_copula[i, j] <- results[r, 3]
  LMF_map[i, j]     <- as.numeric(results[r, 4])
}

# Convert to dataframe for plotting the best copula
df_best_copula <- expand.grid(lon = lon, lat = lat)
df_best_copula$best_copula <- as.vector(best_copula)
df_best_copula <- na.omit(df_best_copula)  # Remove NA rows

# Convert to dataframe for plotting LMF
df_LMF <- expand.grid(lon = lon, lat = lat)
df_LMF$LMF <- as.vector(LMF_map)
df_LMF <- na.omit(df_LMF)  # Remove NA rows

# Define colors for the copulas
copula_colors <- c(
  "Gumbel"    = "#ADD8E6",
  "Clayton"   = "#90EE90",
  "Frank"     = "#FFFF99",
  "Joe"       = "#4682B4",
  "Gumbel90"  = "#FFA07A",  # Light Salmon
  "Clayton90" = "#DDA0DD",  # Plum
  "Joe90"     = "#FFD700"   # Gold
)

# Plot the best copula
p_best_copula <- ggplot(df_best_copula, aes(x = lon, y = lat, fill = best_copula)) +
  geom_raster() +
  scale_fill_manual(values = copula_colors, name = "Best Copula") +
  borders("world", colour = "black", size = 0.3) +
  coord_cartesian(xlim = range(lon), ylim = range(lat)) +  # <- Restrict plot to input data extent
  theme_minimal() +
  labs(title = "Best-Fitting Copula (Jan All Data)",
       x = "Longitude", y = "Latitude") +
  theme(legend.position = "right")


# Plot the LMF
p_LMF <- ggplot(df_LMF, aes(x = lon, y = lat, fill = LMF)) +
  geom_raster() +
  scale_fill_viridis_c(name = "LMF (Likelihood Multiplier Factor)") +
  borders("world", colour = "black", size = 0.3) +
  coord_cartesian(xlim = range(lon), ylim = range(lat)) +  # <- Restrict plot to input data extent
  theme_minimal() +
  labs(title = "Likelihood Multiplier Factor (LMF) for Hot-Dry Events (Jul)",
       x = "Longitude", y = "Latitude") +
  theme(legend.position = "right")

# Print the plot
print(p_LMF)

print(p_best_copula)

# Save Best Copula as CSV
write.csv(df_best_copula, "/g/data/gb02/nt4273/Copula_Data/Global_Best_Copula_Jul(HD).csv",
          row.names = FALSE)

# Save LMF as CSV
write.csv(df_LMF, "/g/data/gb02/nt4273/Copula_Data/Global_LMF_Jul(HD).csv",
          row.names = FALSE)




