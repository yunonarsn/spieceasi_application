# ===============================
# 1. Load libraries
# ===============================
library(igraph)
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
make_ten_hub_adj <- function(n_per_hub = 20, connect_hubs = FALSE, w = 1) {
  n_hubs <- 10
  p <- n_hubs * n_per_hub
  A <- matrix(0, p, p)

  # Loop over each hub
  for (hub in 1:n_hubs) {
    hub_center <- (hub - 1) * n_per_hub + 1
    other_nodes <- ((hub - 1) * n_per_hub + 2):(hub * n_per_hub)

    if (length(other_nodes) > 0) {
      A[hub_center, other_nodes] <- w
      A[other_nodes, hub_center] <- w
    }
  }

  # Optionally connect hubs
  if (connect_hubs) {
    hub_centers <- seq(1, p, by = n_per_hub)
    A[hub_centers, hub_centers] <- w
    diag(A) <- 0
  }

  diag(A) <- 0
  return(A)
}
replace_non_zero_values_from_the_selected_topology <- function(topology, p = 5) {
  topology[topology != 0] <- 1/p
  return(topology)
}

# ======================================
# 3. Simulate ten-hub topology and create plot
# ======================================
hub_topology <- make_ten_hub_adj(n_per_hub = 20, connect_hubs = FALSE, w = 1)
g <- graph_from_adjacency_matrix(hub_topology, mode = "undirected", diag = FALSE)
deg <- degree(g)
max_deg <- max(deg)
hub_nodes <- which(deg == max_deg)

vertex_colors <- rep("#377EB8", vcount(g))
vertex_colors[hub_nodes] <- "#E41A1C"

ten_hub <- plot(g, vertex.size = 5, vertex.label = NA, vertex.color = vertex_colors, layout = layout_with_fr)
ggsave(filename = "Plots/ten_hub_network.pdf", plot = ten_hub, height = 9, width = 10)

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
results_10 <- data.frame(
  method = rep(methods, each = n_reps),
  hamming = NA,
  jaccard = NA,
  f1 = NA,
  sim = rep(1:n_reps, times = length(methods))
)

all_lambda_df_10 <- data.frame()
recall_df_10 <- data.frame()
edge_prob_df_10 <- data.frame()
opt_lambdas_10 <- data.frame()

# --- Simulation loop ---
for (i in 1:n_reps) {
  n <- 200
  df <- sim_data(n, cov_mat)
  lmax <- getMaxCov(df, cov = FALSE)
  lams <- getLamPath(lmax, lmax*0.001, len=60)

  for (m in methods) {
    hugeargs <- list(lambda = lams, verbose = FALSE, method = tolower(m))

    se <- pulsar(
      as.matrix(df),
      fun = huge,
      fargs = hugeargs,
      rep.num = 60,
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
    idx <- which(results_10$method == m & results_10$sim == i)
    results_10$hamming[idx] <- hamming(est, hub_topology)
    results_10$jaccard[idx] <- jaccard(est, hub_topology)
    results_10$f1[idx] <- F1_Score(est, hub_topology)

    # --- Sparsity & instability per lambda ---
    lambda_df_10 <- data.frame(
      lambda = lams,
      sparsity = se$est$sparsity,
      instability = se$stars$summary,
      method = m,
      sim = i
    )
    all_lambda_df_10 <- rbind(all_lambda_df_10, lambda_df_10)

    # --- Optimal lambda ---
    opt_idx_10 <- se$stars$opt.index

    opt_lambdas_10 <- rbind(
      opt_lambdas_10,
      data.frame(
        method = m,
        sim = i,
        opt_lambda = se$est$lambda[opt_idx_10]
      )
    )
    # --- Edge probabilities ---
    edge_probs <- as.matrix(se$stars$merge[[se$stars$opt.index]])

    edge_prob_df_10 <- rbind(
      edge_prob_df_10,
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

    recall_df_10 <- rbind(
      recall_df_10,
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
# Lambda vs sparsity / instability: all_lambda_df_10
# Recall vs threshold: recall_df_10
# Edge probability distributions: edge_prob_df_10
# Hamming/Jaccard/F1 per method: results_10


# ======================================
# 7. Lambda vs. Sparisty plot
# ======================================
avg_lambda_10 <- all_lambda_df_10 %>%
  group_by(method, lambda) %>%
  summarise(sparsity = mean(sparsity), .groups = "drop")

p200.plot.sparsity <- ggplot(avg_lambda_10, aes(x=lambda, y=sparsity, color=method)) +
  geom_line(size=1.3) +
  labs(title="Lambda vs. Sparsity", x="Lambda", y="Sparsity", subtitle = "p = 200, n = 200") +
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
  )

# ======================================
# 8. Lambda vs. Edge Instability plot
# ======================================
avg_instability_10 <- all_lambda_df_10 %>%
  group_by(method, lambda) %>%
  summarise(instability = mean(instability), .groups = "drop")

p200.plot.instability <- ggplot(avg_instability_10, aes(x=lambda, y=instability, color=method)) +
  geom_line(size=1.3) +
  labs(title="Lambda vs Edge Instability (STARS)", x="Lambda", y="Edge Instability", subtitle = "p = 200, n = 200") +
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
  )+
  geom_hline(yintercept = 0.05, linetype="dashed", color = "black", size=1)


# ======================================
# 9. Recall vs. Edge Threshold plot
# ======================================
avg_recall_10 <- recall_df_10 %>%
  group_by(method, threshold) %>%
  summarise(recall = mean(recall), .groups = "drop")

p200.plot.recall <- ggplot(avg_recall_10, aes(x=threshold, y=recall, color=method)) +
  geom_line(size=1.3) +
  labs(title="Recall vs Edge Probability Threshold",
       x="Edge Probability Threshold", y="Recall", subtitle = "p = 200, n = 200") +
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
  )

# ======================================
# 10. Edge probability distribution plot
# ======================================
p200.plot.edgedist <- ggplot(edge_prob_df_10, aes(x=edge_prob, fill=method, color=method)) +
  geom_density(alpha=0.7) +
  labs(title="Edge Probability Distributions",
       x="Edge Probability", y="Density", subtitle = "p = 200, n = 200") +
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
results_10_long <- results_10 %>%
  tidyr::pivot_longer(cols = c(hamming, jaccard, f1),
                      names_to = "metric", values_to = "value")
results_10_long$metric <- factor(
  results_10_long$metric,
  levels = c("hamming", "jaccard", "f1"),
  labels = c("Hamming Distance", "Jaccard Index", "F1 Score")
)

p200.plot.distance <- ggplot(results_10_long, aes(x=method, y=value, fill=method)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~metric, scales="free_y") +
  labs(
    title = "Hamming Distance, Jaccard Index and F1 Score across Simulations",
    y = NULL, subtitle = "p = 200, n = 200"
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
opt_lambdas_10 %>%
  group_by(method) %>%
  summarise(mean_opt_lambda = mean(opt_lambda))

p200.plot.sparsity
p200.plot.instability
p200.plot.recall
p200.plot.edgedist
p200.plot.distance
# ======================================
# 12. Save all plots as pdf
# ======================================
ggsave(filename = "Plots/p200n200sparisty.pdf", plot = p200.plot.sparsity, height = 9, width = 10)
ggsave(filename = "Plots/p200n200instability.pdf", plot = p200.plot.instability, height = 9, width = 10)
ggsave(filename = "Plots/p200n200recall.pdf", plot = p200.plot.recall, height = 9, width = 10)
ggsave(filename = "Plots/p200n200edgedist.pdf", plot = p200.plot.edgedist, height = 8, width = 13)
ggsave(filename = "Plots/p200n200distance.pdf", plot = p200.plot.distance, height = 9, width = 13)

# --- Heatmap of example network (last network that ran in simulation)---
heatmap <-pheatmap(est,
         cluster_rows = FALSE,  # keep nodes in original order
         cluster_cols = FALSE,
         color = colorRampPalette(c("#FAFAFA", "#E41A1C"))(50),
         main = "Estimated Network Heatmap",
         fontsize = 20,
         border_color = "grey")

ggsave(filename = "Plots/heatmap.pdf", plot = heatmap, height = 10, width = 12)

# --- median values for hamming distance and the other metrics---
results_10_long %>%
  group_by(method, metric) %>%
  summarise(median = median(value))
