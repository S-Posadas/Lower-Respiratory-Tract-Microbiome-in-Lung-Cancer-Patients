#### Package setup ####

Sys.setenv(language = "EN")

library(magrittr)
library(dplyr)
library(ggplot2)
library(MMUPHin)
library(vegan) 

theme_set(theme_bw())

#### End ####

#### Plot theme settings ####

text_size = 20

thm.r <- theme(text = element_text(size = text_size),
             #  axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               legend.position = "right")

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

res.dir <- "Results"
qc.dir <- file.path(res.dir, "1.QC")
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir,"5.phyloseq.filtered.RData"))

# Remove controls and samples with 0 counts from relative abundance phyloseq object
physeq_re_samples = subset_samples(physeq_re, sample_names(physeq_re) %in% c(rownames(sample_info[!sample_info$Sample_or_Control %in% c("Negative control sample", "Positive control sample"),])))
physeq_re_samples = subset_samples(physeq_re_samples, !is.na(colSums(otu_table(physeq_re_samples)))) 

# Batch correction
sample_data(physeq_re_samples)$Batch <- as.factor(sample_data(physeq_re_samples)$Batch)
fit_adjust_batch <- adjust_batch(feature_abd = otu_table(physeq_re_samples),
                                 batch = "Batch",
                                 data = data.frame(sample_data(physeq_re_samples)),
                                 control = list(verbose = FALSE))

# Create corrected phyloseq object
physeq_re_bc <- phyloseq(otu_table(fit_adjust_batch$feature_abd_adj, taxa_are_rows = T),
                         sample_data(physeq_re_samples), 
                         tax_table(physeq_re_samples),
                         phy_tree(physeq_re_samples),
                         refseq(physeq_re_samples))

# Variance analysis 
# Set method (e.g.: "bray", "wunifrac"...)
vmethod = "unifrac"
D_before <- phyloseq::distance(physeq_re_samples, method = vmethod)
D_after <- phyloseq::distance(physeq_re_bc, method = vmethod)

run_adonis <- function(d) {
  adonis2(d ~ Batch + Lung + Diagnosis + History.of.smoking.y.n, 
          na.action = na.exclude, 
          data = data.frame(sample_data(physeq_re_samples)), 
          permutations = 999, 
          by = "margin")
}

fit_adonis_before <- run_adonis(D_before)
fit_adonis_after <- run_adonis(D_after)

# Prepare plotting data
prepare_df <- function(fit, var = "R2") {
  vars <- fit[,"R2"][1:4]
  names(vars) <- rownames(fit)[1:4]
  c(vars, Residual = 1 - sum(vars))
}

df <- data.frame(
  Variable = c(names(prepare_df(fit_adonis_before)), names(prepare_df(fit_adonis_after))),
  Variance = c(prepare_df(fit_adonis_before) * 100, prepare_df(fit_adonis_after) * 100),
  Pval = c(prepare_df(fit_adonis_before, "Pr(>F)"), prepare_df(fit_adonis_after,"Pr(>F)")),
  Correction = rep(c("Before correction", "After correction"), each = 5)
  ) %>%
    mutate(Variable = recode(Variable, "History.of.smoking.y.n" = "History of smoking"),
           Variable = factor(Variable, levels = c("Residual", "Batch", "Lung", "Diagnosis", "History of smoking")),
           Correction = factor(Correction, levels = c("Before correction", "After correction")))

write.xlsx(df, file.path(qc.dir, paste0("MMUPHin_correction_", vmethod, ".xlsx")))

# Plot
ggplot(df, aes(x = Correction, y = Variance, fill = Variable)) +
    geom_col(position = "stack") +
    labs(x = "", y = "Variance Explained (%)", fill = NULL) +
    thm.r

ggsave(file.path(qc.dir, paste0("MMUPHin_correction_", vmethod, ".svg")), width = 8, height = 8, dpi = 300)  

# Filter out Residual and plot with zoomed-in y-axis
ggplot(subset(df, Variable != "Residual"), 
       aes(x = Correction, y = Variance, fill = Variable)) +
  geom_col(position = "stack") +
  coord_cartesian(ylim = c(0, 15)) +  # Adjust upper limit as needed
  labs(x = "", y = "Variance Explained (%)", fill = NULL) +
  thm.r

ggsave(file.path(qc.dir, paste0("MMUPHin_correction_no_residuals_", vmethod, ".svg")), width = 8, height = 4, dpi = 300)

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image(file.path(R.dir, "6.phyloseq.corrected.RData"))

#### End ####