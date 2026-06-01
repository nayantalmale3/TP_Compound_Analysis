# Author: Nayan Talmale
# This code selects the best copula model by RMSE (2x2 grid) for Sydney

# Load required packages
library(copula)
library(VineCopula)
library(ggplot2)
library(dplyr)
library(gridExtra)

set.seed(271)

# Load the data over Sydney from 1940-2020 for March
data <- read.table("/Users/NTALMALE/Documents/Topic_I_like/My_PhD_work/My_code/Copulas/Before_All_Code/Point_Sydney_Mar.csv",
                   header = TRUE, sep = ",")

# Step 1: Separate even and odd years data (Training and Testing)
Even_data <- data[seq(1, nrow(data), by = 2), ]
Odd_data  <- data[seq(2, nrow(data), by = 2), ]

# Step 2: Extract temperature and precipitation for both datasets
Temp_test  <- Even_data[, 2]
Rain_test  <- Even_data[, 3]
Temp_train <- Odd_data[, 2]
Rain_train <- Odd_data[, 3]

# Step 3: Transform the data using probability integral transform
Temp_odd_u  <- pobs(Temp_train)
Precip_odd_v <- pobs(Rain_train)
Temp_even_u <- pobs(Temp_test)
Precip_even_v <- pobs(Rain_test)

# Number of observations
n <- length(Temp_test)

# ------------- Copula Fitting -------------
fit.f <- BiCopEst(Temp_odd_u, Precip_odd_v, family = 5, method="mle")  # Frank
fit.g <- BiCopEst(Temp_odd_u, Precip_odd_v, family = 4, method="mle")  # Gumbel
fit.c <- BiCopEst(Temp_odd_u, Precip_odd_v, family = 3, method="mle")  # Clayton
fit.j <- BiCopEst(Temp_odd_u, Precip_odd_v, family = 6, method="mle")  # Joe

# ------------- Simulate data from fitted copulas -------------
simulate_copula <- function(fit, n) {
  sim_data <- BiCopSim(n, family = fit$family, par = fit$par)
  sim_temp <- quantile(Temp_train, probs = sim_data[,1])
  sim_precip <- quantile(Rain_train, probs = sim_data[,2])
  df <- data.frame(U = pobs(sim_temp), V = pobs(sim_precip))
  return(df)
}

df_F <- simulate_copula(fit.f, n)
df_G <- simulate_copula(fit.g, n)
df_C <- simulate_copula(fit.c, n)
df_j <- simulate_copula(fit.j, n)
df_emp <- data.frame(U = Temp_even_u, V = Precip_even_v)  # Testing (empirical)

# ------------- Assign 2x2 grid indices -------------
assign_grid_index <- function(U, V, U_median, V_median) {
  cut_U <- cut(U, breaks = c(-Inf, U_median, Inf), labels = 1:2, include.lowest = TRUE)
  cut_V <- cut(V, breaks = c(-Inf, V_median, Inf), labels = 1:2, include.lowest = TRUE)
  grid_index <- as.numeric(cut_U) + 2 * (as.numeric(cut_V) - 1)
  return(grid_index)
}

# Compute medians
U_q <- median(Temp_even_u)
V_q <- median(Precip_even_v)
U_q_F <- median(df_F$U); V_q_F <- median(df_F$V)
U_q_G <- median(df_G$U); V_q_G <- median(df_G$V)
U_q_C <- median(df_C$U); V_q_C <- median(df_C$V)
U_q_j <- median(df_j$U); V_q_j <- median(df_j$V)

# Assign grid indices
df_emp$Grid_Index <- assign_grid_index(df_emp$U, df_emp$V, U_q, V_q)
df_F$Grid_Index <- assign_grid_index(df_F$U, df_F$V, U_q_F, V_q_F)
df_G$Grid_Index <- assign_grid_index(df_G$U, df_G$V, U_q_G, V_q_G)
df_C$Grid_Index <- assign_grid_index(df_C$U, df_C$V, U_q_C, V_q_C)
df_j$Grid_Index <- assign_grid_index(df_j$U, df_j$V, U_q_j, V_q_j)

# ------------- Compute RMSE for each 2x2 box -------------
compute_rmse <- function(empirical, theoretical) {
  boxes <- 1:4
  rmse <- numeric(4)
  emp_counts <- table(empirical)
  theo_counts <- table(theoretical)
  for (b in boxes) {
    O <- ifelse(is.na(emp_counts[as.character(b)]), 0, emp_counts[as.character(b)])
    E <- ifelse(is.na(theo_counts[as.character(b)]), 0, theo_counts[as.character(b)])
    rmse[b] <- sqrt((O - E)^2)
  }
  return(matrix(rmse, nrow = 2, ncol = 2, byrow = TRUE))
}

rmse_F <- compute_rmse(df_emp$Grid_Index, df_F$Grid_Index)
rmse_G <- compute_rmse(df_emp$Grid_Index, df_G$Grid_Index)
rmse_C <- compute_rmse(df_emp$Grid_Index, df_C$Grid_Index)
rmse_j <- compute_rmse(df_emp$Grid_Index, df_j$Grid_Index)

# Convert RMSE to data frames for plotting
rmse_df <- function(rmse_matrix) {
  df <- as.data.frame(as.table(rmse_matrix))
  colnames(df) <- c("Row", "Column", "RMSE")
  return(df)
}

rmse_df_F <- rmse_df(rmse_F)
rmse_df_G <- rmse_df(rmse_G)
rmse_df_C <- rmse_df(rmse_C)
rmse_df_j <- rmse_df(rmse_j)

# ------------- Plot function for copula (no colorbar) -------------
plot_copula <- function(name, df_syn, df_rmse) {
  # Synthetic scatter
  p1 <- ggplot(df_syn, aes(x = U, y = V)) +
    geom_point(color = "blue") +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    labs(title = paste(name, "Synthetic_Odd (Train)"), y = "Precip (V)", x = ifelse(name=="Joe","Temp (U)","")) +
    theme_minimal()
  
  # Empirical scatter
  p2 <- ggplot(df_emp, aes(x = U, y = V)) +
    geom_point(color = "red") +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    labs(title = paste("Empirical_Even (Test)"), x = ifelse(name=="Joe","Temp (U)","")) +
    theme_minimal()
  
  # RMSE heatmap without colorbar
  df_rmse$RMSE_disc <- cut(df_rmse$RMSE, breaks = c(-Inf,0.5,1.5,2.5,3.5,Inf), labels = 0:4)
  df_rmse$RMSE_disc <- as.numeric(as.character(df_rmse$RMSE_disc))
  
  p3 <- ggplot(df_rmse, aes(x = Column, y = Row, fill = RMSE_disc)) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(RMSE,2)), color = "black") +
    scale_fill_gradientn(colors = c("lightblue","yellow","lavender","yellow","darkorange"), 
                         limits = c(0,4), guide = "none") +  # remove colorbar
    labs(title = paste(name, "RMSE")) +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  
  return(list(p1,p2,p3))
}

# ------------- Generate plots for standard copulas -------------
copulas_list <- list(
  Frank=list(synthetic=df_F, rmse=rmse_df_F),
  Gumbel=list(synthetic=df_G, rmse=rmse_df_G),
  Clayton=list(synthetic=df_C, rmse=rmse_df_C),
  Joe=list(synthetic=df_j, rmse=rmse_df_j)
)

all_plots <- lapply(names(copulas_list), function(n) plot_copula(n, copulas_list[[n]]$synthetic, copulas_list[[n]]$rmse))
plots_flat <- unlist(all_plots, recursive = FALSE)

# Arrange plots in a 3-column grid
grid.arrange(grobs = plots_flat, ncol = 3)





