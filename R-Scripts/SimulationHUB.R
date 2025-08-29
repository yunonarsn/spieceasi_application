# ===============================
# 1. Load libraries
# ===============================
library(igraph)
library(MASS)
library(pheatmap)
library(huge)
library(pulsar)
library(MLmetrics)
library(brew)
library(SpiecEasi)
library(tidyverse)
library(ggplot2)

set.seed(123)
# ======================================
# 2. Helper functions to simulate topology and replace non zero values in precision matrix
# ======================================
make_two_hub_adj <- function(n_per_hub = 20, connect_hubs = FALSE, w = 1) {
  p  <- 2 * n_per_hub
  A  <- matrix(0, p, p)
  c1 <- 1L
  c2 <- n_per_hub + 1L

  A[c1, 2:n_per_hub] <- w
  A[2:n_per_hub, c1] <- w

  A[c2, (n_per_hub + 2):p] <- w
  A[(n_per_hub + 2):p, c2] <- w

  diag(A) <- 0
  A
}
replace_non_zero_values_from_the_selected_topology <- function(topology, p = 5) {
  topology[topology != 0] <- 1/p
  return(topology)
}

# ======================================
# 3. Simulate two-hub topology and create plot
# ======================================
hub_topology = make_two_hub_adj(n_per_hub = 20)
g <- graph_from_adjacency_matrix(hub_topology, mode = "undirected", diag = FALSE)
deg <- degree(g)
max_deg <- max(deg)
hub_nodes <- which(deg == max_deg)

vertex_colors <- rep("#377EB8", vcount(g))
vertex_colors[hub_nodes] <- "#E41A1C"

two_hub <- plot(g, vertex.size = 5, vertex.label = NA, vertex.color = vertex_colors, layout = layout_with_fr)
ggsave(filename = "Plots/two_hub_network.pdf", plot = two_hub, height = 9, width = 13)

# ======================================
# 4. Create precision and covariance matrix to simulate data
# ======================================
True_network <- hub_topology

precision_matrix <- replace_non_zero_values_from_the_selected_topology(True_network, p= 5)
diag(precision_matrix) <- 1
precision <- precision_matrix
cov_mat <- prec2cov(precision)

sim_data <- function(n, cov_mat) {
  MASS::mvrnorm(n = n, mu = rep(0, nrow(cov_mat)), Sigma = cov_mat)
}

# ======================================
# 5. Helper functions for hamming distance and jaccard index
# ======================================
hamming <- function(est, true) {
  sum(est != true)
}

jaccard <- function(est, true) {
  intersection <- sum(est & true)
  union <- sum((est + true) > 0)
  intersection / union
}

# ======================================
# 6.Actual Simulation
# ======================================
# --- Parameters ---
n_reps <- 10
methods <- c("MB", "glasso")
thresholds <- seq(0.01, 1, by = 0.01)

# Initialize where to store all results
results <- data.frame(
  method = rep(methods, each = n_reps),
  hamming = NA,
  jaccard = NA,
  f1 = NA,
  sim = rep(1:n_reps, times = length(methods))
)

all_lambda_df <- data.frame()
recall_df <- data.frame()
edge_prob_df <- data.frame()
opt_lambdas <- data.frame()

# --- Simulation loop ---
for (i in 1:n_reps) {
  n <- 200
  df <- sim_data(n, cov_mat)
  lmax <- getMaxCov(df, cov = FALSE)
  lams <- getLamPath(lmax, lmax*0.01, len=30)

  for (m in methods) {
    hugeargs <- list(lambda = lams, verbose = FALSE, method = tolower(m))

    se <- pulsar(
      as.matrix(df),
      fun = huge,
      fargs = hugeargs,
      rep.num = 30,
      criterion = "stars",
      lb.stars = TRUE,
      ub.stars = TRUE,
      seed = 10010,
      thresh = 0.05
    )

    # --- Refit network ---
    fit <- refit(se)
    est <- as.matrix(fit$refit$stars)

    # --- Metrics ---
    idx <- which(results$method == m & results$sim == i)
    results$hamming[idx] <- hamming(est, hub_topology)
    results$jaccard[idx] <- jaccard(est, hub_topology)
    results$f1[idx] <- F1_Score(est, hub_topology)

    # --- Sparsity & instability per lambda ---
    lambda_df <- data.frame(
      lambda = lams,
      sparsity = se$est$sparsity,
      instability = se$stars$summary,
      method = m,
      sim = i
    )
    all_lambda_df <- rbind(all_lambda_df, lambda_df)

    # --- Optimal lambda ---
    opt_idx <- se$stars$opt.index

    opt_lambdas <- rbind(
      opt_lambdas,
      data.frame(
        method = m,
        sim = i,
        opt_lambda = se$est$lambda[opt_idx]
      )
    )
    # --- Edge probabilities ---
    edge_probs <- as.matrix(se$stars$merge[[se$stars$opt.index]])

    edge_prob_df <- rbind(
      edge_prob_df,
      data.frame(
        edge_prob = as.vector(edge_probs),
        method = m,
        sim = i
      )
    )

    # --- Recall vs threshold ---
    true_edges <- hub_topology != 0
    recall_vec <- sapply(thresholds, function(th) {
      est_edges <- edge_probs >= th
      sum(est_edges & true_edges) / max(sum(est_edges), 1)  # avoid div by 0
    })

    recall_df <- rbind(
      recall_df,
      data.frame(
        threshold = thresholds,
        recall = recall_vec,
        method = m,
        sim = i
      )
    )
    print("done with method")
  }
}

# --- Results ready for plotting ---
# Lambda vs sparsity / instability: all_lambda_df
# Recall vs threshold: recall_df
# Edge probability distributions: edge_prob_df
# Hamming/Jaccard/F1 per method: results


# ======================================
# 7. Lambda vs. Sparisty plot
# ======================================
avg_lambda <- all_lambda_df %>%
  group_by(method, lambda) %>%
  summarise(sparsity = mean(sparsity), .groups = "drop")

p5.plot.sparsity <- ggplot(avg_lambda, aes(x=lambda, y=sparsity, color=method)) +
  geom_line(size=1.3) +
  labs(title="Lambda vs. Sparsity", x="Lambda", y="Sparsity", subtitle = "p = 40, n = 200") +
  theme_bw() +
  scale_color_brewer(palette = "Set1", name = "Method")+
  theme(
    plot.title = element_text(
      size = 28, face = "bold", hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 24, hjust = 0.5
    ),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size=26),
    legend.text = element_text(size=24)
  ) + ylim(0, 1)

# ======================================
# 8. Lambda vs. Edge Instability plot
# ======================================
avg_instability <- all_lambda_df %>%
  group_by(method, lambda) %>%
  summarise(instability = mean(instability), .groups = "drop")

p5.plot.instability <- ggplot(avg_instability, aes(x=lambda, y=instability, color=method)) +
  geom_line(size=1.3) +
  labs(title="Lambda vs Edge Instability (STARS)", x="Lambda", y="Edge Instability", subtitle = "p = 40, n = 200") +
  theme_bw() +
  scale_color_brewer(palette = "Set1", name = "Method")+
  theme(
    plot.title = element_text(
      size = 28, face = "bold", hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 24, hjust = 0.5
    ),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size=26),
    legend.text = element_text(size=24)
  ) +
  ylim(0, 0.4) +
  geom_hline(yintercept = 0.05, linetype="dashed", color = "black", size=1)

# ======================================
# 9. Recall vs. Edge Threshold plot
# ======================================
avg_recall <- recall_df %>%
  group_by(method, threshold) %>%
  summarise(recall = mean(recall), .groups = "drop")

p5.plot.recall <- ggplot(avg_recall, aes(x=threshold, y=recall, color=method)) +
  geom_line(size=1.3) +
  labs(title="Recall vs Edge Probability Threshold",
       x="Edge Probability Threshold", y="Recall", subtitle = "p = 40, n = 200") +
  theme_bw()+
  scale_color_brewer(palette = "Set1", name = "Method")+
  theme(
    plot.title = element_text(
      size = 28, face = "bold", hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 24, hjust = 0.5
    ),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size=26),
    legend.text = element_text(size=24)
  ) + ylim(0, 1)


# ======================================
# 10. Edge probability distribution plot
# ======================================
p5.plot.edgedist <- ggplot(edge_prob_df, aes(x=edge_prob, fill=method, color=method)) +
  geom_density(alpha=0.7) +
  labs(title="Edge Probability Distributions",
       x="Edge Probability", y="Density", subtitle = "p = 40, n = 200") +
  facet_wrap(~method) +
  theme_bw()+
  scale_fill_brewer(palette = "Set1")+
  scale_color_brewer(palette = "Set1")+
  theme(
    plot.title = element_text(
      size = 28, face = "bold", hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 24, hjust = 0.5
    ),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# ======================================
# 11. PLot with Hamming distance, jaccard index and F1 Score
# ======================================
results_long <- results %>%
  tidyr::pivot_longer(cols = c(hamming, jaccard, f1),
                      names_to = "metric", values_to = "value")
results_long$metric <- factor(
  results_long$metric,
  levels = c("hamming", "jaccard", "f1"),
  labels = c("Hamming Distance", "Jaccard Index", "F1 Score")
)

p5.plot.distance <- ggplot(results_long, aes(x=method, y=value, fill=method)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~metric, scales="free_y") +
  labs(
    title = "Hamming Distance, Jaccard Index and F1 Score across Simulations",
    y = NULL, subtitle = "p = 40, n = 200"
  ) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1", name = "Method")  +
  theme(
    plot.title = element_text(
      size = 28, face = "bold", hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 24, hjust = 0.5
    ),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    strip.text = element_text(size = 22),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# --- Average optimal lambda for each method ---
opt_lambdas %>%
  group_by(method) %>%
  summarise(mean_opt_lambda = mean(opt_lambda))

p5.plot.sparsity
p5.plot.instability
p5.plot.recall
p5.plot.edgedist
p5.plot.distance
# ======================================
# 12. Save all plots as pdf
# ======================================
ggsave(filename = "Plots/p40n200sparisty.pdf", plot = p5.plot.sparsity, height = 9, width = 10)
ggsave(filename = "Plots/p40n200instability.pdf", plot = p5.plot.instability, height = 9, width = 10)
ggsave(filename = "Plots/p40n200recall.pdf", plot = p5.plot.recall, height = 9, width = 10)
ggsave(filename = "Plots/p40n200edgedist.pdf", plot = p5.plot.edgedist, height = 8, width = 13)
ggsave(filename = "Plots/p40n200distance.pdf", plot = p5.plot.distance, height = 9, width = 13)

# --- Heatmap of network ---
pheatmap(est,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#E41A1C"))(50),
         main = "Estimated Network Heatmap",
         fontsize = 12)

# --- median values for hamming distance and the other metrics---
results_long %>%
  group_by(method, metric) %>%
  summarise(median = median(value))
