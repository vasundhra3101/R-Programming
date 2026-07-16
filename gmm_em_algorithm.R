# ===================================================================
# VECTORIZED EXPECTATION-MAXIMIZATION (EM) ALGORITHM
# FOR GAUSSIAN MIXTURE MODELS FROM SCRATCH (IMPROVED)
# ===================================================================

# Clear environment
rm(list = ls())
cat("\014")

# Create Output directory
if (!dir.exists("Output")) {
  dir.create("Output")
}

# ===================================================================
# PART 1: Generate Dataset
# ===================================================================

cat("\n========================================\n")
cat("PART 1: Generating Synthetic Dataset\n")
cat("========================================\n")

set.seed(123)
n <- 200

# Cluster 1
mu1_true <- c(2, 3)
sigma1_true <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)

# Cluster 2
mu2_true <- c(7, 8)
sigma2_true <- matrix(c(1, -0.3, -0.3, 1), nrow = 2)

# Generate data manually
generate_data <- function(n, mean, sigma) {
  L <- chol(sigma)
  Z <- matrix(rnorm(n * 2), ncol = 2)
  X <- Z %*% L + matrix(mean, n, 2, byrow = TRUE)
  return(X)
}

data1 <- generate_data(n/2, mu1_true, sigma1_true)
data2 <- generate_data(n/2, mu2_true, sigma2_true)
X <- rbind(data1, data2)

write.csv(X, "synthetic_data.csv", row.names = FALSE)
cat("Dataset saved to synthetic_data.csv\n")
cat("First 6 rows:\n")
print(head(X))

# ===================================================================
# PART 2: Initialize Parameters (Improved)
# ===================================================================

cat("\n========================================\n")
cat("PART 2: Initializing Parameters (Improved)\n")
cat("========================================\n")

K <- 2
n <- nrow(X)
d <- ncol(X)

# Better initialization: Use k-means-like initialization
set.seed(456)

# Initialize with points far apart
# Find two points that are far from each other
dist_matrix <- as.matrix(dist(X))
max_dist_idx <- which(dist_matrix == max(dist_matrix), arr.ind = TRUE)[1, ]
mu_k <- X[max_dist_idx, ]

# If we got the same point twice, use random initialization
if (nrow(mu_k) < K) {
  mu_k <- X[sample(1:n, K), ]
}

pi_k <- rep(1/K, K)

# Initialize covariance matrices with slight variation
sigma_k <- list()
for (k in 1:K) {
  # Use variance of data for each dimension
  sigma_k[[k]] <- diag(apply(X, 2, var)) * 0.5 + diag(1e-6, d)
}

cat("Initial pi:", pi_k, "\n")
cat("Initial mu:\n")
print(mu_k)

# ===================================================================
# PART 3: Gaussian Density Function
# ===================================================================

cat("\n========================================\n")
cat("PART 3: Gaussian Density Function\n")
cat("========================================\n")

gaussian_density <- function(X, mu, sigma) {
  n <- nrow(X)
  d <- ncol(X)
  density <- numeric(n)
  
  # Check if matrix is singular
  if (rcond(sigma) < 1e-15) {
    sigma <- sigma + diag(1e-6, d)
  }
  
  sigma_inv <- solve(sigma)
  sigma_det <- det(sigma)
  
  # Numerical stability
  if (sigma_det < 1e-300) {
    sigma_det <- 1e-300
  }
  
  for (i in 1:n) {
    x_mu <- X[i, ] - mu
    quad_form <- t(x_mu) %*% sigma_inv %*% x_mu
    # Prevent overflow
    if (quad_form > 700) quad_form <- 700
    density[i] <- (1 / ((2*pi)^(d/2) * sqrt(sigma_det))) *
      exp(-0.5 * quad_form)
  }
  return(density)
}

# ===================================================================
# PART 4: Expectation Step (E-Step)
# ===================================================================

cat("\n========================================\n")
cat("PART 4: Expectation Step (E-Step)\n")
cat("========================================\n")

E_step <- function(X, pi_k, mu_k, sigma_k) {
  n <- nrow(X)
  K <- length(pi_k)
  gamma <- matrix(0, n, K)
  
  for (k in 1:K) {
    gamma[, k] <- pi_k[k] * gaussian_density(X, mu_k[k, ], sigma_k[[k]])
  }
  
  # Normalize with numerical stability
  row_sum <- rowSums(gamma)
  # Avoid division by zero
  row_sum[row_sum < 1e-300] <- 1e-300
  gamma <- gamma / row_sum
  
  return(gamma)
}

# ===================================================================
# PART 5: Maximization Step (M-Step)
# ===================================================================

cat("\n========================================\n")
cat("PART 5: Maximization Step (M-Step)\n")
cat("========================================\n")

M_step <- function(X, gamma) {
  n <- nrow(X)
  d <- ncol(X)
  K <- ncol(gamma)
  
  Nk <- colSums(gamma)
  
  # Prevent division by zero
  Nk[Nk < 1e-10] <- 1e-10
  
  pi_k <- Nk / n
  mu_k <- matrix(0, K, d)
  sigma_k <- list()
  
  for (k in 1:K) {
    # Update mean
    mu_k[k, ] <- colSums(gamma[, k] * X) / Nk[k]
    
    # Update covariance
    sigma <- matrix(0, d, d)
    for (i in 1:n) {
      x_mu <- matrix(X[i, ] - mu_k[k, ], ncol = 1)
      sigma <- sigma + gamma[i, k] * (x_mu %*% t(x_mu))
    }
    sigma_k[[k]] <- sigma / Nk[k] + diag(1e-6, d)
  }
  
  return(list(pi_k = pi_k, mu_k = mu_k, sigma_k = sigma_k))
}

# ===================================================================
# PART 6: Log Likelihood
# ===================================================================

cat("\n========================================\n")
cat("PART 6: Log Likelihood Function\n")
cat("========================================\n")

log_likelihood <- function(X, pi_k, mu_k, sigma_k) {
  n <- nrow(X)
  K <- length(pi_k)
  logL <- 0
  
  for (i in 1:n) {
    temp <- 0
    for (k in 1:K) {
      temp <- temp + pi_k[k] * gaussian_density(X[i, , drop=FALSE], mu_k[k, ], sigma_k[[k]])
    }
    # Prevent log of zero
    if (temp < 1e-300) temp <- 1e-300
    logL <- logL + log(temp)
  }
  return(logL)
}

# ===================================================================
# PART 7: EM Algorithm (Improved)
# ===================================================================

cat("\n========================================\n")
cat("PART 7: Running EM Algorithm\n")
cat("========================================\n")

max_iter <- 50  # Increased max iterations
logL_values <- numeric(max_iter)
tol <- 1e-6     # Tighter tolerance

cat("\nIteration\tLog-Likelihood\t\tChange\n")
cat("----------------------------------------\n")

for (iter in 1:max_iter) {
  # E-Step
  gamma <- E_step(X, pi_k, mu_k, sigma_k)
  
  # M-Step
  params <- M_step(X, gamma)
  pi_k <- params$pi_k
  mu_k <- params$mu_k
  sigma_k <- params$sigma_k
  
  # Log-Likelihood
  logL_values[iter] <- log_likelihood(X, pi_k, mu_k, sigma_k)
  
  # Print progress
  if (iter == 1) {
    cat(iter, "\t\t", logL_values[iter], "\t\t--\n")
  } else {
    change <- logL_values[iter] - logL_values[iter-1]
    cat(iter, "\t\t", logL_values[iter], "\t\t", change, "\n")
  }
  
  # Convergence check
  if (iter > 1 && abs(logL_values[iter] - logL_values[iter-1]) < tol) {
    cat("\nConvergence achieved at iteration", iter, "\n")
    break
  }
}

logL_values <- logL_values[1:iter]

# ===================================================================
# PART 8: Final Cluster Assignment
# ===================================================================

cat("\n========================================\n")
cat("PART 8: Final Cluster Assignment\n")
cat("========================================\n")

clusters <- apply(gamma, 1, which.max)
cat("Cluster distribution:\n")
print(table(clusters))

cat("\nFinal parameters:\n")
cat("Mixing coefficients (pi):", pi_k, "\n")
cat("Mean vectors (mu):\n")
print(mu_k)

# ===================================================================
# PART 9: Visualization
# ===================================================================

cat("\n========================================\n")
cat("PART 9: Visualization\n")
cat("========================================\n")

# Dataset Plot
png("Output/dataset.png", width = 800, height = 600)
plot(X, col = "blue", pch = 16, main = "Generated Dataset", xlab = "X1", ylab = "X2")
dev.off()
cat("Dataset plot saved to Output/dataset.png\n")

# Final Clusters
png("Output/final_clusters.png", width = 800, height = 600)
plot(X, col = clusters, pch = 16, main = "Final Clusters", xlab = "X1", ylab = "X2")
points(mu_k, col = 1:K, pch = 8, cex = 2, lwd = 2)
legend("topright", legend = c(paste("Cluster", 1:K), "Centers"), 
       col = c(1:K, "black"), pch = c(16, 16, 8))
dev.off()
cat("Final clusters plot saved to Output/final_clusters.png\n")

# Log Likelihood
png("Output/log_likelihood.png", width = 800, height = 600)
plot(logL_values, type = "l", main = "Log-Likelihood Convergence", 
     xlab = "Iteration", ylab = "Log-Likelihood", col = "blue", lwd = 2)
grid()
points(logL_values, col = "red", pch = 16)
dev.off()
cat("Log-likelihood plot saved to Output/log_likelihood.png\n")

# ===================================================================
# PART 10: True vs Estimated Comparison
# ===================================================================

cat("\n========================================\n")
cat("PART 10: True vs Estimated Comparison\n")
cat("========================================\n")

true_clusters <- c(rep(1, n/2), rep(2, n/2))

# Fix label switching if needed
if (sum(clusters == 1 & true_clusters == 1) < sum(clusters == 1 & true_clusters == 2)) {
  # Swap cluster labels
  clusters_swapped <- clusters
  clusters_swapped[clusters == 1] <- 2
  clusters_swapped[clusters == 2] <- 1
  clusters <- clusters_swapped
}

png("Output/true_vs_estimated.png", width = 1000, height = 500)
par(mfrow = c(1, 2))

# True clusters
plot(X, col = true_clusters, pch = 16,
     main = "True Clusters",
     xlab = "X1", ylab = "X2")
points(rbind(mu1_true, mu2_true), col = 1:2, pch = 8, cex = 2, lwd = 2)
legend("topright", legend = c("Cluster 1", "Cluster 2", "Centers"), 
       col = c(1, 2, "black"), pch = c(16, 16, 8))

# Estimated clusters
plot(X, col = clusters, pch = 16,
     main = "Estimated Clusters (EM)",
     xlab = "X1", ylab = "X2")
points(mu_k, col = 1:K, pch = 8, cex = 2, lwd = 2)
legend("topright", legend = c("Cluster 1", "Cluster 2", "Centers"), 
       col = c(1, 2, "black"), pch = c(16, 16, 8))

dev.off()
cat("True vs Estimated comparison plot saved to Output/true_vs_estimated.png\n")

# ===================================================================
# RESULTS SUMMARY
# ===================================================================

cat("\n========================================\n")
# ===================================================================
# RESULTS SUMMARY
# ===================================================================

cat("\n========================================\n")
cat("RESULTS SUMMARY\n")
cat("========================================\n")

cat("\nEstimated Parameters:\n")
cat("Mixing Coefficients:\n")
print(pi_k)

cat("\nEstimated Means:\n")
print(mu_k)

cat("\nEstimated Covariance Matrices:\n")
for(k in 1:K){
  cat("\nCluster", k, "Covariance Matrix:\n")
  print(sigma_k[[k]])
}

# ===================================================================
# FINAL OUTPUT
# ===================================================================

cat("\n========================================\n")
cat("EM ALGORITHM COMPLETED SUCCESSFULLY\n")
cat("========================================\n")

cat("\nOutput Files Generated:\n")
cat("1. Output/dataset.png\n")
cat("2. Output/final_clusters.png\n")
cat("3. Output/log_likelihood.png\n")

cat("\nSummary:\n")
cat("Total Data Points :", n, "\n")
cat("Number of Clusters :", K, "\n")
cat("Iterations :", length(logL_values), "\n")
cat("Final Log-Likelihood :", tail(logL_values, 1), "\n")

cat("\n========================================\n")