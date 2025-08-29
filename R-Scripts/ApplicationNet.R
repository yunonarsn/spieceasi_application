# ===============================
# 1. Load libraries
# ===============================
library(phyloseq)
library(microbiome)
library(SpiecEasi)
library(igraph)
library(NetCoMi)
library(limma)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
# ===============================
# 2. Collapse to species level
# ===============================
agp_species <- tax_glom(agp_obeseunderweight_processed, taxrank = "Rank7", NArm = FALSE)

# ===============================
# 3. Prevalence filtering
# ===============================
prev_filter <- function(x, cutoff = 0.1) {
  keep <- apply(otu_table(x), 1, function(row) {
    sum(row > 0) / nsamples(x) >= cutoff
  })
  prune_taxa(keep, x)
}

agp_species_filt <- prev_filter(agp_species, cutoff = 0.1)


# ===============================
# 4. Define groups
# ===============================
sample_data(agp_species_filt)$group <- sample_data(agp_species_filt)$bmi

agp_obese <- subset_samples(agp_species_filt, group == 1)   # obese
agp_under <- subset_samples(agp_species_filt, group == 0)   # underweight

# ===============================
# 5. NetCoMi analysis
# ===============================

net_group <- netConstruct(data = agp_obese, data2 = agp_under,
                          measure = "spieceasi",
                          filtTax = "none",
                          sparsMethod = "threshold",   # or "topology"
                          thresh = 0.3,
                          verbose = 3)

net_groups_props <- netAnalyze(net_group)

network <- plot(net_groups_props,
     sameLayout = TRUE,
     title1 = "Obese",
     title2 = "Underweight",
     type = "net",
     layout = "spring",
     vertex.size = 5,
     edge.color = "grey70",
     vertex.color = tax_table(agp_species_filt)[, "Rank4"],
     nodeSize = "degree",
     nodeColorLegendTitle = "Taxa",
     fontsize = 20)


net_grous_comp1 <- netCompare(net_groups_props, permTest = FALSE, seed = 10010)
summary(net_grous_comp1)

# ===============================
# 6. GCD plot
# ===============================
adja1 <- net_group$adjaMat1
adja2 <- net_group$adjaMat2
gcd <- calcGCD(adja1, adja2)
gcmtest <- testGCM(gcd)
plotHeat(gcmtest$gcm1,
         pmat = gcmtest$pAdjust1,
         type = "mixed")
title(main = "GCM1 (obese)", cex.main = 1.5)

plotHeat(gcmtest$gcm2,
         pmat = gcmtest$pAdjust2,
         type = "mixed")
title(main = "GCM2 (underweight)", cex.main = 1.5)

plotHeat(gcmtest$diff,
         pmat = gcmtest$pAdjustDiff,
         type = "mixed",
         textLow = "code")
title(main = "GCM Difference Heatmap (GCM1 - GCM2)", cex.main = 1.5)
# ===============================
# 7. Abundance Heatmap
# ===============================

abund_mat <- as.matrix(otu_table(agp_species_filt))

abund_mat_log <- log1p(abund_mat)

# Simple heatmap
abundances <- pheatmap(abund_mat_log,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("white", "red"))(50),
         main = "Abundance heatmap (log transformed)",
         fontsize = 20,
         name = "Abundance")
ggsave(abundances, filename = "Plots/Abundance_heatmap_log_scaled.pdf", width = 10, height = 8)


# Gruppierung
groups <- sample_data(agp_species_filt)$group
groups <- factor(groups, levels = c("0", "1"))
names(groups) <- colnames(abund_mat_log)
col_order <- names(sort(groups))
# Reorder the abundance matrix
abund_mat_log_ordered <- abund_mat_log[, col_order]

# Reorder the annotation factor
groups_ordered <- groups[col_order]# Farben für die Annotation
ann_colors <- list(Group = c("0" = "blue", "1" = "red"))

# Annotation-Objekt für die Spalten
ha <- HeatmapAnnotation(
  Group = groups_ordered,
  col = ann_colors,
  annotation_legend_param = list(Group = list(title = "Group"))
)

# Farbpalette für Abundanzen
col_fun <- colorRamp2(
  c(min(abund_mat_log_ordered), 0, max(abund_mat_log_ordered)),
  c("white", "pink", "red")
)

# Heatmap zeichnen
Heatmap(abund_mat_log_ordered,
        name = "Abundance",
        col = col_fun,
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_title = "Abundance heatmap by group (log transformed)",
        )


