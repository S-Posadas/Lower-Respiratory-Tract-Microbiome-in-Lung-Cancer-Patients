# HAllA plots: diagnosis and smoking vs ASVs #
# Run through scripts/halla, after HAllA.
#
# - Copies of the HAllA output with the ASV names labelled with their taxonomy,
#   for hallagram (run by scripts/halla):
#   labelled/            all significant clusters (Supplementary Figure 6A)
#   labelled_diagnosis/  clusters associated with diagnosis (Figure 2E)
# - Violin plots of the relative abundance of the ASVs in the clusters
#   associated with diagnosis (Supplementary Figure 6B)

#### Package setup ####

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

theme_set(theme_bw())

#### End ####

#### Load data ####

res.dir <- "Results_direct"
halla.res.dir <- file.path(res.dir, "halla", "diagnosis", "ASVs_nmi_uniform")

rdata <- new.env()
load(file.path(res.dir, "RData", "6.phyloseq.corrected.RData"), envir = rdata)
physeq_halla <- rdata$physeq_re_bc
diagnosis_col <- rdata$diagnosis_col
rm(rdata)

# Same samples as in halla_input.R
sd <- data.frame(sample_data(physeq_halla))
keep <- sd$Sample_or_Control == "True sample" &
  sd$Isolation == "Direct_Isolation" &
  sd$Lung == "Diseased" &
  sd$Diagnosis %in% c("Benign", "NSCLC", "SCLC") &
  !is.na(sd$History.of.smoking.y.n)
physeq_halla <- prune_samples(keep, physeq_halla)
physeq_halla <- prune_taxa(taxa_sums(physeq_halla) > 0, physeq_halla)

sig_clusters <- read.delim(file.path(halla.res.dir, "sig_clusters.txt"), stringsAsFactors = FALSE)

# One row per ASV and significant cluster
cluster_info <- sig_clusters %>%
  separate_rows(cluster_X, sep = ";") %>%
  select(ASV = cluster_X, cluster_rank, cluster_condition = cluster_Y, pvalue = best_adjusted_pvalue) %>%
  distinct()

#### End ####

#### Labelled copies of the HAllA output for hallagram ####

# Label: "ASV412 (Micrococcus spp.)" - species ("<Genus> spp." if several
# species are possible), or the lowest assigned rank
tax <- as.data.frame(tax_table(physeq_halla))
lowest_rank <- apply(tax[, c("Species", "Genus", "Family", "Order", "Class", "Phylum")], 1, function(r) {
  if (!is.na(r[["Species"]]) && grepl("/", r[["Species"]]) && !is.na(r[["Genus"]])) return(paste(r[["Genus"]], "spp."))
  r[!is.na(r) & r != ""][1]
})
asv_label <- setNames(paste0(rownames(tax), " (", lowest_rank, ")"), rownames(tax))

relabel <- function(x) unname(ifelse(x %in% names(asv_label), asv_label[x], x))

write_labelled_copy <- function(out_name, clusters) {
  out_dir <- file.path(halla.res.dir, out_name)
  dir.create(out_dir, showWarnings = FALSE)
  # Data tables with the ASVs (first column) relabelled
  for (f in c("X.tsv", "X_original.tsv")) {
    d <- read.delim(file.path(halla.res.dir, f), check.names = FALSE, stringsAsFactors = FALSE)
    d[[1]] <- relabel(d[[1]])
    write.table(d, file.path(out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  for (f in c("Y.tsv", "Y_original.tsv", "performance.txt")) file.copy(file.path(halla.res.dir, f), out_dir, overwrite = TRUE)
  assoc <- read.delim(file.path(halla.res.dir, "all_associations.txt"), stringsAsFactors = FALSE)
  assoc$X_features <- relabel(assoc$X_features)
  write.table(assoc, file.path(out_dir, "all_associations.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  clusters$cluster_X <- sapply(strsplit(clusters$cluster_X, ";"), function(a) paste(relabel(a), collapse = ";"))
  write.table(clusters, file.path(out_dir, "sig_clusters.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

write_labelled_copy("labelled", sig_clusters)
write_labelled_copy("labelled_diagnosis", sig_clusters[sig_clusters$cluster_Y == "Diagnosis", ])

#### End ####

#### Supplementary Figure 6B: ASVs of the clusters associated with diagnosis ####

# Facet label: cluster rank, ASV and species (genus if the species is not
# resolved, i.e. "spp." or several species)
facet_rank <- apply(tax[, c("Species", "Genus", "Family", "Order", "Class", "Phylum")], 1, function(r) {
  if (!is.na(r[["Species"]]) && !grepl("spp|/", r[["Species"]])) return(r[["Species"]])
  r[-1][!is.na(r[-1])][1]
})

diag_clusters <- cluster_info %>% filter(cluster_condition == "Diagnosis")
ps_plot <- prune_taxa(taxa_names(physeq_halla) %in% diag_clusters$ASV, physeq_halla)

violin_df <- psmelt(ps_plot) %>%
  mutate(Abundance = Abundance * 100) %>%
  filter(Abundance > 0) %>%
  left_join(diag_clusters %>% select(ASV, cluster_rank), by = c("OTU" = "ASV")) %>%
  mutate(Cluster = paste("Cluster Nr. ", cluster_rank),
         ASV_tax = paste0(OTU, "\n", facet_rank[OTU]))
violin_df$Cluster <- factor(violin_df$Cluster, levels = unique(violin_df$Cluster[order(violin_df$cluster_rank)]))

p_violin <- ggplot(violin_df, aes(x = Diagnosis, y = Abundance, color = Diagnosis, fill = Diagnosis)) +
  geom_violin(alpha = 0.7) +
  geom_jitter(width = 0.2) +
  facet_wrap(~ Cluster + ASV_tax, nrow = 3) +
  scale_fill_manual(values = diagnosis_col) + scale_color_manual(values = diagnosis_col) +
  scale_x_discrete(drop = FALSE) +
  scale_y_log10(labels = scales::label_number(decimal.mark = ",", big.mark = ".")) +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme(text = element_text(size = 30),
        strip.text = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")
ggsave(file.path(halla.res.dir, "diagnosis_violin.svg"), p_violin, dpi = 300, width = 14, height = 12)

write.table(cluster_info %>% mutate(label = asv_label[ASV]), file.path(halla.res.dir, "sig_clusters_per_ASV.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#### End ####
