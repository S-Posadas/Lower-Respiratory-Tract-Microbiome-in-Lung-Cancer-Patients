# PICRUSt2 pathway analysis #
# Predicted MetaCyc pathway abundances (PICRUSt2, scripts/picrust2) in the
# tumour-bearing lung of NSCLC patients, compared with LinDA (ggpicrust2):
# - M stage, M0 vs M1:      heatmap (Figure 4C)
# - N stage, N0 vs N1-3:    combined bar plot (Supplementary Figure 10E)
# Pathways present in >= 10% of the samples with mean relative abundance
# >= 1e-4 are tested; significant: BH-adjusted p < 0.1.

#### Package setup ####

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggpicrust2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(vegan)
})

#### End ####

#### Load data ####

res.dir <- "Results_direct"
picrust.dir <- file.path(res.dir, "picrust2")

# Sample data of the filtered phyloseq object
rdata <- new.env()
load(file.path(res.dir, "RData", "5.phyloseq.filtered.RData"), envir = rdata)
sd_all <- data.frame(sample_data(rdata$physeq_filtered), check.names = F)
rm(rdata)

# Unstratified pathway abundances (PICRUSt2)
pathway_abundance <- read.delim(file.path(picrust.dir, "picrust2_out", "pathways_out", "path_abun_unstrat.tsv"),
                                check.names = F) %>%
  column_to_rownames("pathway")

#### End ####

#### Settings ####

prev_threshold <- 0.1    # present in >= 10% of the samples
mean_threshold <- 1e-4   # mean relative abundance
padj_threshold <- 0.1    # BH-adjusted p-value
n_top <- 20              # pathways in the plots (lowest p-values)

group_colors <- c("skyblue", "firebrick")      # bar plots
heatmap_colors <- c("lightblue", "firebrick")  # heatmap annotation

#### End ####

#### Functions ####

# Samples: tumour-bearing lung, NSCLC, stage known
stage_metadata <- function(stage) {
  md <- sd_all %>%
    filter(Sample_or_Control == "True sample", Isolation == "Direct_Isolation", LibrarySizeDecontam > 0,
           Lung == "Diseased", Diagnosis == "NSCLC", !is.na(.data[[stage]])) %>%
    rownames_to_column("SampleID")
  # Samples in the PICRUSt2 output, in its order
  md[match(intersect(colnames(pathway_abundance), md$SampleID), md$SampleID), ]
}

# Prevalence and mean relative abundance filter
filter_pathways <- function(abund) {
  rel <- sweep(abund, 2, colSums(abund), "/")
  keep <- rowSums(abund > 0) >= prev_threshold * ncol(abund) & rowMeans(rel) >= mean_threshold
  cat("Pathways:", nrow(abund), "->", sum(keep), "after filtering\n")
  abund[keep, ]
}

# LinDA (reference: first level) and MetaCyc descriptions. log2_ratio is
# log2(mean relative abundance group 2 / group 1)
run_linda <- function(abund, md, group, reference) {
  res <- pathway_daa(abundance = abund, metadata = md, group = group, daa_method = "LinDA",
                     reference = reference, p_adjust_method = "BH", include_abundance_stats = TRUE)
  res <- pathway_annotation(pathway = "MetaCyc", daa_results_df = res)
  res$description[is.na(res$description)] <- res$feature[is.na(res$description)]
  res %>%
    mutate(log2_ratio = log2(mean_rel_abundance_group2 / mean_rel_abundance_group1)) %>%
    arrange(p_values)
}

# Combined bar plot: mean relative abundance (+- SE), log2 fold change and
# -log10(adjusted p) of the top pathways (lowest p at the top)
plot_combined_bar <- function(top, md, group, labels, file) {
  n <- table(md[[group]])[c(1, 2)]
  feature_order <- rev(top$description)

  abund_long <- top %>%
    transmute(description,
              g1_mean = mean_rel_abundance_group1, g1_se = sd_rel_abundance_group1 / sqrt(n[1]),
              g2_mean = mean_rel_abundance_group2, g2_se = sd_rel_abundance_group2 / sqrt(n[2])) %>%
    pivot_longer(-description, names_to = c("group", ".value"), names_sep = "_") %>%
    mutate(group = factor(labels[group], levels = labels),
           description = factor(description, levels = feature_order))

  p_abund <- ggplot(abund_long, aes(x = description, y = mean, fill = group)) +
    geom_col(position = position_dodge(0.9), color = "black", width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.9), width = 0.2) +
    coord_flip() +
    labs(x = "", y = "Mean relative abundance", fill = sub("_grouped", "", group)) +
    scale_fill_manual(values = setNames(group_colors, labels)) +
    theme_minimal() +
    theme(legend.position = "left", axis.text.y = element_text(size = 12, face = "bold"))

  up <- paste("Up in", labels[2]); down <- paste("Down in", labels[2])
  logfc <- top %>%
    mutate(description = factor(description, levels = feature_order),
           direction = case_when(p_adjust < padj_threshold & log2_ratio > 0 ~ up,
                                 p_adjust < padj_threshold & log2_ratio < 0 ~ down,
                                 TRUE ~ "Not significant"))

  p_logfc <- ggplot(logfc, aes(x = description, y = log2_ratio, fill = direction)) +
    geom_col(color = "black", width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    coord_flip() +
    labs(x = "", y = "Log2 fold change") +
    scale_fill_manual(values = setNames(c(group_colors[2], group_colors[1], "gray"), c(up, down, "Not significant"))) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.y = element_blank())

  p_pval <- ggplot(logfc, aes(x = description, y = -log10(p_adjust), group = 1)) +
    geom_line(color = "darkgreen", linewidth = 0.8) +
    geom_point(color = "darkgreen", size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    coord_flip() +
    labs(x = "", y = "-log10(adjusted p-value)") +
    theme_minimal() +
    theme(axis.text.y = element_blank())

  p <- p_abund + p_logfc + p_pval + plot_layout(ncol = 3, guides = "collect")
  ggsave(file, p, width = 16, height = max(6, nrow(top) * 0.4))
}

# Heatmap of log2(abundance + 1); samples and pathways clustered by
# Bray-Curtis distance (average linkage)
plot_heatmap <- function(top, abund, md, group, colors, file) {
  features <- intersect(rownames(abund), top$feature)  # order of the abundance table
  mat <- as.matrix(abund[features, md$SampleID, drop = F])
  rownames(mat) <- top$description[match(features, top$feature)]

  # Pseudocount (half the smallest non-zero value) for the distances
  mat_pseudo <- mat + min(mat[mat > 0]) / 2
  dist_cols <- as.matrix(vegdist(t(mat_pseudo), method = "bray"))
  dist_cols[is.nan(dist_cols)] <- 0
  dist_rows <- vegdist(mat_pseudo, method = "bray")

  annotation_col <- data.frame(md[[group]], row.names = md$SampleID)
  colnames(annotation_col) <- group
  annotation_colors <- list(setNames(colors, levels(md[[group]])))
  names(annotation_colors) <- group

  # Larger legends (set before the heatmap is built)
  ComplexHeatmap::ht_opt(legend_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
                         legend_labels_gp = grid::gpar(fontsize = 13),
                         legend_grid_width = grid::unit(7, "mm"),
                         legend_grid_height = grid::unit(7, "mm"))

  hmp <- ComplexHeatmap::pheatmap(log2(mat + 1),
                                  annotation_col = annotation_col,
                                  annotation_colors = annotation_colors,
                                  show_rownames = T, show_colnames = F, fontsize_row = 14,
                                  clustering_distance_cols = as.dist(dist_cols),
                                  clustering_distance_rows = dist_rows,
                                  clustering_method = "average",
                                  name = "log2\n(abundance + 1)",
                                  heatmap_legend_param = list(legend_height = grid::unit(5, "cm")),
                                  # full pathway names (default maximum: 6 cm)
                                  row_names_max_width = ComplexHeatmap::max_text_width(rownames(mat), gp = grid::gpar(fontsize = 14)))

  # Larger name of the annotation row (e.g. "M")
  hmp@top_annotation@anno_list[[group]]@name_param$gp <- grid::gpar(fontsize = 16, fontface = "bold")

  svg(file, width = 16, height = max(8, nrow(mat) * 0.3))
  ComplexHeatmap::draw(hmp, heatmap_legend_side = "right", annotation_legend_side = "right",
                       legend_gap = grid::unit(10, "mm"))
  dev.off()
  ComplexHeatmap::ht_opt(RESET = TRUE)
}

# Whole comparison: filter, LinDA, tables and plots
pathway_comparison <- function(md, group, labels, out_dir) {
  dir.create(out_dir, recursive = T, showWarnings = F)
  md[[group]] <- factor(md[[group]])
  cat("\n", group, ": ", paste(names(table(md[[group]])), table(md[[group]]), sep = " n=", collapse = ", "), "\n", sep = "")

  abund <- filter_pathways(pathway_abundance[, md$SampleID])
  res <- run_linda(abund, md, group, reference = levels(md[[group]])[1])
  write.csv(res, file.path(out_dir, "pathways_LinDA_all.csv"), row.names = F)

  sig <- filter(res, p_adjust < padj_threshold)
  cat("Significant pathways (adjusted p <", padj_threshold, "):", nrow(sig), "\n")
  write.csv(sig, file.path(out_dir, "pathways_LinDA_significant.csv"), row.names = F)
  if (nrow(sig) == 0) return(invisible(res))

  top <- head(sig, n_top)
  plot_combined_bar(top, md, group, labels, file.path(out_dir, "pathways_barplot_combined.svg"))
  plot_heatmap(top, abund, md, group, heatmap_colors, file.path(out_dir, "pathways_heatmap.svg"))
  invisible(res)
}

#### End ####

#### M0 vs M1 (Figure 4C) ####

md_M <- stage_metadata("M")
md_M$M <- as.character(md_M$M)
res_M <- pathway_comparison(md_M, "M", labels = c(g1 = "M0", g2 = "M1"),
                            out_dir = file.path(picrust.dir, "M"))

#### End ####

#### N0 vs N1-3 (Supplementary Figure 10E) ####

md_N <- stage_metadata("N")
md_N$N_grouped <- ifelse(md_N$N == 0, "0", "1-3")
res_N <- pathway_comparison(md_N, "N_grouped", labels = c(g1 = "N0", g2 = "N1-3"),
                            out_dir = file.path(picrust.dir, "N"))

#### End ####
