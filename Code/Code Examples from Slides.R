# -----------------------------
# Optimal Transport
# Code Examples
# 
# Author: Anthony Christidis
# CCB, DBMI, HMS
# 
# -----------------------------

# _____________________
# Monge Example
# _____________________

# Use assignment problem solver from clue package
library(clue)

# Define cost matrix for assignment
cost_matrix <- matrix(
    c(8, 2, 5, 9, 7,
      6, 4, 3, 7, 5,
      7, 8, 6, 5, 6,
      4, 6, 7, 3, 9,
      9, 7, 5, 2, 8), 
    nrow = 5, byrow = TRUE
)

# Solve the optimal assignment problem
assignment <- solve_LSAP(cost_matrix, maximum = FALSE)

# Print the optimal assignment
print(assignment)

# Calculate total cost of optimal assignment
total_cost <- sum(cost_matrix[cbind(seq_along(assignment), assignment)])
print(paste("Total Cost:", total_cost))


# _____________________
# Karantovich Example
# _____________________

# Define the source and target distributions as probability vectors
a <- c(0.2, 0.2, 0.2, 0.2, 0.2)
b <- c(0.15, 0.25, 0.25, 0.15, 0.2)

# Define the cost matrix (example values for demonstration)
cost_matrix <- matrix(c(
    1, 2, 3, 4, 5,
    2, 1, 3, 4, 5,
    3, 3, 1, 2, 5,
    4, 4, 2, 1, 3,
    5, 5, 5, 3, 1
), nrow = 5, byrow = TRUE)

# Solve the Kantorovich problem using the transport function
solution <- transport(a, b, cost_matrix, method = "network")

# Print the optimal transport plan
print("Optimal Transport Plan:")
round(solution, 2)

# Compute the total transport cost
total_cost <- sum(solution * cost_matrix)
print(paste("Total Transport Cost:", total_cost))

# _____________________
# Continuous Example
# _____________________

# Library to simulate data
library(MASS) 

# Set seed for reproducibility
set.seed(1)

# Number of samples
n_samples <- 2000

# Source distribution: Gaussian centered at (0, 0)
mu_source <- c(0, 0)
sigma_source <- diag(2)
source_samples <- mvrnorm(n_samples, mu = mu_source, Sigma = sigma_source)

# Target distribution: Gaussian centered at (1, 1)
mu_target <- c(1, 1)
sigma_target <- diag(2)
target_samples <- mvrnorm(n_samples, mu = mu_target, Sigma = sigma_target)

# Discretize by creating histograms or empirical measures
source_weights <- rep(1 / n_samples, n_samples)
target_weights <- rep(1 / n_samples, n_samples)

# Compute cost matrix (Euclidean distance squared)
cost_matrix <- as.matrix(dist(rbind(source_samples, target_samples))^2)[1:n_samples, (n_samples+1):(2*n_samples)]

# Solve the Kantorovich OT problem
ot_plan <- transport(source_weights, target_weights, cost_matrix, method = "network")

# Display the resulting transport plan
print("Optimal Transport Plan:")
print(head(ot_plan, 10))  

# Compute the total transport cost
total_cost <- sum(ot_plan * cost_matrix)
cat("Total Transport Cost:", total_cost, "\n")

# Check specific transformations
print("Expected Shift Verification:")

# Compare the first few mappings
for (i in 1:10) {  
    source_index <- ot_plan$from[i]
    target_index <- ot_plan$to[i]
    
    # Compare the original source and target coordinates
    original_source <- source_samples[source_index,]
    intended_target <- target_samples[target_index,]
    
    # Print results
    cat(sprintf("Source (%.2f, %.2f) mapped to Target (%.2f, %.2f)\n",
                original_source[1], original_source[2],
                intended_target[1], intended_target[2]))
}

# Create a data frame for plotting
transport_map <- data.frame(
    SourceX = source_samples[ot_plan$from, 1],
    SourceY = source_samples[ot_plan$from, 2],
    TargetX = target_samples[ot_plan$to, 1],
    TargetY = target_samples[ot_plan$to, 2]
)

# Plot with ggplot2
library(ggplot2)
ggplot(transport_map) + 
    geom_point(aes(x = SourceX, y = SourceY), color = "blue", alpha = 0.5) +
    geom_point(aes(x = TargetX, y = TargetY), color = "red", alpha = 0.5) +
    geom_segment(aes(x = SourceX, y = SourceY, xend = TargetX, yend = TargetY), 
                 arrow = arrow(length = unit(0.2, "cm")), alpha = 0.3) +
    labs(title = "Optimal Transport Mapping", x = "X", y = "Y") +
    theme_minimal()

# __________________________
# scDiagnostics Application
# __________________________

# Load data
data("reference_data")
data("query_data")

# Extract CD4 cells
ref_data_subset <- reference_data[, which(reference_data$expert_annotation == "CD4")]
query_data_subset <- query_data[, which(query_data$expert_annotation == "CD4")]

# Selecting highly variable genes (can be customized by the user)
ref_top_genes <- scran::getTopHVGs(ref_data_subset, n = 500)
query_top_genes <- scran::getTopHVGs(query_data_subset, n = 500)

# Intersect the gene symbols to obtain common genes
common_genes <- intersect(ref_top_genes, query_top_genes)
ref_data_subset <- ref_data_subset[common_genes,]
query_data_subset <- query_data_subset[common_genes,]

# Run PCA on reference data
ref_data_subset <- scater::runPCA(ref_data_subset)

# Compute Wasserstein distances and compare using quantile-based permutation test
wasserstein_data <- calculateWassersteinDistance(query_data = query_data_subset,
                                                 reference_data = ref_data_subset,
                                                 query_cell_type_col = "expert_annotation",
                                                 ref_cell_type_col = "expert_annotation",
                                                 pc_subset = 1:5,
                                                 n_resamples = 300)
plot(wasserstein_data)
