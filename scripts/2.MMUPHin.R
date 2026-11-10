# Batch correction with MMUPHin, for the directly isolated or the culture-enriched samples #
# Set `dataset` below ("direct" or "cult"); run the script once for each.
# Input: 5.phyloseq.filtered.RData (1.DADA2_direct.R or 1.DADA2_culture.R)
# Output: 6.phyloseq.corrected.RData and physeq_re_bc.rds (batch-corrected relative
# abundances of the true samples), used by beta.R, the taxonomy scripts and merge_phyloseq.R.
# Supplementary Document 1, Figure 5: variance explained by sequencing batch and
# other variables before and after correction (Bray-Curtis and weighted UniFrac).
# The culture-enriched object has no tree, so only Bray-Curtis is used there.

#### Package setup ####

Sys.setenv(language = "EN")

library(dplyr)
library(ggplot2)
library(phyloseq)
library(MMUPHin)
library(vegan)
library(openxlsx)

theme_set(theme_bw())

#### End ####

#### Plot theme settings ####

text_size = 20

thm.r <- theme(text = element_text(size = text_size),
               axis.ticks.x = element_blank(),
               legend.position = "right")

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

dataset <- "direct"   # "direct" or "cult"

mmuphin.res.dir <- switch(dataset,
                          direct = "Results_direct",
                          cult   = "Results_cult")
load(file.path(mmuphin.res.dir, "RData", "5.phyloseq.filtered.RData"))

# Set the folders after load(): the checkpoint contains its own res.dir/qc.dir/R.dir
# (relative paths in the older culture checkpoints)
res.dir <- mmuphin.res.dir
qc.dir <- file.path(res.dir, "1.QC")
R.dir <- file.path(res.dir, "RData")

#### End ####

#### Batch correction ####

# Remove controls and samples with 0 counts from relative abundance phyloseq object
physeq_re_samples = subset_samples(physeq_re, sample_names(physeq_re) %in% c(rownames(sample_info[!sample_info$Sample_or_Control %in% c("Negative control sample", "Positive control sample"),])))
physeq_re_samples = subset_samples(physeq_re_samples, !is.na(colSums(otu_table(physeq_re_samples))))

# Batch correction
sample_data(physeq_re_samples)$Batch <- as.factor(sample_data(physeq_re_samples)$Batch)
fit_adjust_batch <- adjust_batch(feature_abd = otu_table(physeq_re_samples),
                                 batch = "Batch",
                                 data = data.frame(sample_data(physeq_re_samples)),
                                 control = list(verbose = FALSE))

# Create corrected phyloseq object (with the tree only if there is one: none for the culture-enriched samples)
physeq_re_bc <- phyloseq(otu_table(fit_adjust_batch$feature_abd_adj, taxa_are_rows = T),
                         sample_data(physeq_re_samples),
                         tax_table(physeq_re_samples),
                         refseq(physeq_re_samples))
if (!is.null(phy_tree(physeq_re_samples, errorIfNULL = FALSE))) {
  physeq_re_bc <- merge_phyloseq(physeq_re_bc, phy_tree(physeq_re_samples))
}

#### End ####

#### Variance explained before and after correction ####
# PERMANOVA with batch, lung, diagnosis and smoking status (marginal tests).
# Samples with missing values are excluded (na.exclude).

run_adonis <- function(d) {
  set.seed(123)
  adonis2(d ~ Batch + Lung + Diagnosis + History.of.smoking.y.n,
          na.action = na.exclude,
          data = data.frame(sample_data(physeq_re_samples)),
          permutations = 999,
          by = "margin")
}

# R2 (in %) and p-value of the four variables, plus the residual variance
prepare_df <- function(fit, correction) {
  r2 <- fit[1:4, "R2"]
  data.frame(Variable = c(rownames(fit)[1:4], "Residual"),
             Variance = 100 * c(r2, 1 - sum(r2)),
             Pval = c(fit[1:4, "Pr(>F)"], NA),
             Correction = correction)
}

batch_variance <- list()

vmethods <- if (is.null(phy_tree(physeq_re_samples, errorIfNULL = FALSE))) "bray" else c("bray", "wunifrac")

for (vmethod in vmethods) {
  D_before <- phyloseq::distance(physeq_re_samples, method = vmethod)
  D_after <- phyloseq::distance(physeq_re_bc, method = vmethod)

  df <- rbind(prepare_df(run_adonis(D_before), "Before correction"),
              prepare_df(run_adonis(D_after), "After correction")) %>%
    mutate(Variable = recode(Variable, "History.of.smoking.y.n" = "History of smoking"),
           Variable = factor(Variable, levels = c("Residual", "Batch", "Lung", "Diagnosis", "History of smoking")),
           Correction = factor(Correction, levels = c("Before correction", "After correction")))

  batch_variance[[vmethod]] <- df
  write.xlsx(df, file.path(qc.dir, paste0("MMUPHin_correction_", vmethod, ".xlsx")))

  # Plot
  ggplot(df, aes(x = Correction, y = Variance, fill = Variable)) +
    geom_col(position = "stack") +
    labs(x = "", y = "Variance Explained (%)", fill = NULL) +
    thm.r

  ggsave(file.path(qc.dir, paste0("MMUPHin_correction_", vmethod, ".svg")), width = 8, height = 8, dpi = 300)

  # Without the residual variance, zoomed-in y-axis
  ggplot(subset(df, Variable != "Residual"),
         aes(x = Correction, y = Variance, fill = Variable)) +
    geom_col(position = "stack") +
    coord_cartesian(ylim = c(0, 15)) +
    labs(x = "", y = "Variance Explained (%)", fill = NULL) +
    thm.r

  ggsave(file.path(qc.dir, paste0("MMUPHin_correction_no_residuals_", vmethod, ".svg")), width = 8, height = 4, dpi = 300)
}

batch_variance

#### End ####

#### Supplementary Document 1, Figure 5: both distances in one plot ####
# Bray-Curtis and weighted UniFrac side by side, without the residual variance
# (only when there is a tree, i.e. for the directly isolated samples)

if (length(batch_variance) == 2) {
  fig5_df <- bind_rows(lapply(names(batch_variance), function(v)
    mutate(batch_variance[[v]], Distance = v))) %>%
    filter(Variable != "Residual") %>%
    mutate(Distance = factor(recode(Distance, bray = "Bray-Curtis", wunifrac = "weighted UniFrac"),
                             levels = c("Bray-Curtis", "weighted UniFrac")))

  ggplot(fig5_df, aes(x = Correction, y = Variance, fill = Variable)) +
    geom_col(position = "stack") +
    facet_wrap(~ Distance) +
    coord_cartesian(ylim = c(0, 15)) +
    labs(x = "", y = "Variance Explained (%)", fill = NULL) +
    thm.r

  ggsave(file.path(qc.dir, "MMUPHin_correction_bray_wunifrac.svg"), width = 12, height = 4, dpi = 300)
}

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image(file.path(R.dir, "6.phyloseq.corrected.RData"))

# Save the batch-corrected phyloseq object (the filtered objects are saved by 1.DADA2_direct.R / 1.DADA2_culture.R)
saveRDS(physeq_re_bc, file.path(R.dir, "physeq_re_bc.rds"))

#### End ####
