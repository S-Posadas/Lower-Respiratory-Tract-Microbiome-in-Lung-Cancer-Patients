# MaAsLin2 differential abundance #
# #
# 1. Full-sample MaAsLin2 at ASV, phylum, family and genus level, with 
#    covariates (e.g. Figure 2F, 3C, Supp. Figures 4C, 10C).
# 2. Leave-one-out (LOO) sensitivity analyses (ASV level):
#    - Q3.1 diagnosis: one SCLC or one benign sample left out at a time, each
#      reference level (Supplementary Table 5)
#    - Q5.2 smoking in NSCLC: one never-smoker left out at a time
#      (Supplementary Table 7, Figure 5C)
# 3. Figure 3D (ASV24 by N stage), Figure 5D (never-smoker ASVs in NSCLC vs benign),
#    Supp. Figure 4D (Venn diagram of the ASVs per diagnosis in smokers) and
#    Supplementary Table 6 (ASVs found in only one diagnosis).

#### Package setup ####

Sys.setenv(language = "EN")

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(Maaslin2)
  library(ggvenn)
  library(openxlsx)
})

theme_set(theme_bw())

#### End ####

#### Load batch-corrected phyloseq object (2.MMUPHin.R) ####

res.dir <- "Results_direct"

rdata <- new.env()
load(file.path(res.dir, "RData", "6.phyloseq.corrected.RData"), envir = rdata)
physeq_maas0 <- rdata$physeq_re_bc
lung_col <- rdata$lung_col
diagnosis_col <- rdata$diagnosis_col
histology_col <- rdata$histology_col
smoker_col <- rdata$smoker_col
rm(rdata)

#### End ####

#### Settings and directories ####

maaslin.dir <- file.path(res.dir, "maaslin2")
dir.create(maaslin.dir, recursive = T, showWarnings = F)

min_abund <- 0.0005
min_prev  <- 0.02    
max_sig   <- 0.05    
tax_levels <- c("ASV", "Phylum", "Family", "Genus")

#### End ####

#### Plot theme settings ####

text_size = 30
thm <- theme(text = element_text(size = text_size),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "bottom")

# Axis labels with decimal comma and no scientific notation (as in the paper)
marks_no_sci <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)

#### End ####

#### Functions ####

# Taxonomy names without spaces (MaAsLin2 feature names); "spp." species as NA;
# "Incertae Sedis" families named after their order
clean_tax_names_maaslin <- function(physeq) {
  tax <- as.data.frame(tax_table(physeq))
  tax[] <- lapply(tax, gsub, pattern = " ", replacement = "_")
  tax$Species[grep("_spp.", tax$Species)] <- NA
  tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis", "Family"] <- paste(
    tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis", "Order"], "Incertae_Sedis", sep = "_")
  tax_table(physeq) <- tax_table(as.matrix(tax))
  physeq
}

# Agglomerate to a taxonomic level, taxa named after that level
agglomerate_taxa <- function(physeq, level) {
  phy_agg <- tax_glom(physeq, level)
  taxa_names(phy_agg) <- tax_table(phy_agg)[, level]
  phy_agg
}

# MaAsLin2 with the settings of the paper; returns all results
run_maaslin <- function(physeq, out_dir, fixed_effects, random_effects = NULL,
                        reference = NULL, plots = TRUE) {
  metadata <- data.frame(sample_data(physeq))
  metadata$Age <- as.numeric(as.character(metadata$Age))
  # All metadata as factors, except the numeric Age and Pack.years
  metadata[] <- lapply(colnames(metadata), function(col) {
    if (col %in% c("Age", "Pack.years")) metadata[[col]] else as.factor(metadata[[col]])
  })

  data <- data.frame(otu_table(physeq), check.names = FALSE)
  stopifnot(all(colnames(data) == rownames(metadata)))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)  # Maaslin2 creates only the last folder level
  invisible(capture.output(Maaslin2(
    input_data = data,
    input_metadata = metadata,
    output = out_dir,
    min_abundance = min_abund,
    min_prevalence = min_prev,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = max_sig,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    reference = reference,
    correction = "BH",
    standardize = FALSE,
    plot_heatmap = plots,
    plot_scatter = plots
  )))

  write.table(metadata, file.path(out_dir, "metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  read.delim(file.path(out_dir, "all_results.tsv"), stringsAsFactors = FALSE)
}

# Boxplots of the relative abundance (%) of the given features per group
# (samples with the feature present; log scale)
plot_features <- function(physeq, features, tax_level, group_var, palette, legend_title = group_var) {
  if (tax_level == "ASV") {
    df <- psmelt(prune_taxa(taxa_names(physeq) %in% features, physeq)) %>% mutate(feature = OTU)
  } else {
    df <- psmelt(prune_taxa(tax_table(physeq)[, tax_level] %in% features, physeq)) %>%
      group_by(Sample, feature = .data[[tax_level]], group = .data[[group_var]]) %>%
      summarise(Abundance = sum(Abundance), .groups = "drop") %>%
      rename(!!group_var := group)
  }
  df <- df %>% mutate(Abundance = Abundance * 100) %>% filter(Abundance > 0)
  n <- length(unique(df$feature))
  ggplot(df, aes(x = .data[[group_var]], y = Abundance, color = .data[[group_var]])) +
    geom_boxplot(fill = NA, outlier.shape = NA) +
    geom_point(size = 1.5, alpha = 0.6, position = position_jitter(width = 0.2, height = 0, seed = 1)) +
    facet_wrap(~ feature, ncol = min(5, n)) +
    scale_color_manual(values = palette, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_log10(labels = marks_no_sci) +
    labs(x = NULL, y = "Relative abundance (%)", color = legend_title) +
    thm
}

# Prevalence (share of samples with the feature) per group
prevalence_by_group <- function(physeq, group_var) {
  otu <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) otu <- t(otu)
  grp <- as.character(sample_data(physeq)[[group_var]])
  prev <- sapply(sort(unique(grp)), function(g) rowMeans(otu[, grp == g, drop = FALSE] > 0))
  prev <- as.data.frame(prev)
  colnames(prev) <- paste0("prevalence_", colnames(prev))
  prev$feature <- rownames(prev)
  prev
}

#### End ####

#### Questions to answer ####

#' Questions
#' Q1: Is there a difference between Diseased vs the parallel Normal Lung?
#' Q2: Is there a difference between Diseased vs the parallel Normal Lung without samples with synchronous tumor?
#' Q3.1: What are the differences between NSCLC vs SCLC vs Benign Tumor?
#' Q3.2/Q3.3: The same in nonsmokers / smokers
#' Q3.5: NSCLC vs Benign Tumor (Figure 5D)
#' Q4.1: What are the differences between NSCLC main histologies (Adenocarcinoma vs Squamous)?
#' Q4.2/Q4.3: The same in nonsmokers / smokers
#' Q5: Does history of smoking have an impact on the lung microbiome?
#' Q6.1: Does the T stage have an impact on the lung microbiome?
#' Q6.2: Does the N stage have an impact on the lung microbiome?
#' Q6.3: Does the M stage have an impact on the lung microbiome?

#' Covariates:
#' - Q1, Q2 (paired lungs): no covariates; patient (Study_Nr) as random effect.
#' - Q3.2, Q3.3, Q4.2, Q4.3, all Q5, Q6.3.1: Age + Sex (smoking status is the
#'   stratifying variable or the variable of interest).
#' - Everywhere else: History.of.smoking.y.n + Age + Sex.
#' Reference levels: one run per level in `refs` (NULL = first level).

smoking_age_sex <- c("History.of.smoking.y.n", "Age", "Sex")
age_sex <- c("Age", "Sex")
dx_refs <- c("Benign", "NSCLC", "SCLC")

params <- list(
  Q1.1 = list(name = "Q1.1.all_directly_isolated_lung", fixed = "Lung", col = lung_col, rand = "Study_Nr", refs = "Normal"),
  Q1.2 = list(name = "Q1.2.main-D_directly_isolated_lung", fixed = "Lung", col = lung_col, rand = "Study_Nr", refs = "Normal"),
  Q2.1 = list(name = "Q2.1.all_directly_isolated_non_synchronous_lung", fixed = "Lung", col = lung_col, rand = "Study_Nr", refs = "Normal"),
  Q2.2 = list(name = "Q2.2.main-D_directly_isolated_non_synchronous_lung", fixed = "Lung", col = lung_col, rand = "Study_Nr", refs = "Normal"),
  Q3.1 = list(name = "Q3.1.diagnosis", fixed = "Diagnosis", col = diagnosis_col, covariates = smoking_age_sex, refs = dx_refs),
  Q3.2 = list(name = "Q3.2.diagnosis_nonsmokers", fixed = "Diagnosis", col = diagnosis_col, covariates = age_sex, refs = "Benign"),
  Q3.3 = list(name = "Q3.3.diagnosis_smokers", fixed = "Diagnosis", col = diagnosis_col, covariates = age_sex, refs = dx_refs),
  Q3.5 = list(name = "Q3.5.diagnosis_NSCLCvsbenign", fixed = "Diagnosis", col = diagnosis_col, covariates = smoking_age_sex, refs = "Benign"),
  Q4.1 = list(name = "Q4.1.histology", fixed = "Histology.NSCLC", col = histology_col, lab_plot = "Histology", covariates = smoking_age_sex),
  Q4.2 = list(name = "Q4.2.histology_nonsmokers", fixed = "Histology.NSCLC", col = histology_col, lab_plot = "Histology", covariates = age_sex),
  Q4.3 = list(name = "Q4.3.histology_smokers", fixed = "Histology.NSCLC", col = histology_col, lab_plot = "Histology", covariates = age_sex),
  Q5.1 = list(name = "Q5.1.main-D_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking", covariates = age_sex, refs = "No"),
  Q5.2 = list(name = "Q5.2.NSCLC_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking", covariates = age_sex, refs = "No"),
  Q5.3 = list(name = "Q5.3.Benign_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking", covariates = age_sex, refs = "No"),
  Q5.4 = list(name = "Q5.4.Adenocarcinoma_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking", covariates = age_sex, refs = "No"),
  Q6.1 = list(name = "Q6.1.NSCLC_T", fixed = "T", col = scales::hue_pal()(5)[c(4,2,5,1)], covariates = smoking_age_sex),
  Q6.2 = list(name = "Q6.2.NSCLC_N", fixed = "N", col = scales::hue_pal()(5)[c(3:4,2,5,1)], covariates = smoking_age_sex),
  Q6.2.1 = list(name = "Q6.2.NSCLC_N_twolevels", fixed = "Ntwo", col = scales::hue_pal()(5)[c(3:4,2,5,1)], lab_plot = "N", covariates = smoking_age_sex),
  Q6.3 = list(name = "Q6.3.NSCLC_M", fixed = "M", col = c("lightblue", "firebrick"), covariates = smoking_age_sex),
  Q6.3.1 = list(name = "Q6.3.1.NSCLC_M_smokers", fixed = "M", col = c("lightblue", "firebrick"), covariates = age_sex)
)

# N stage in two levels: N0 vs N1-3
sample_data(physeq_maas0)$Ntwo <- factor(ifelse(sample_data(physeq_maas0)$N == 0, "N0", "N1-3"), levels = c("N0", "N1-3"))

tax_original <- as.data.frame(tax_table(physeq_maas0))  # taxonomy as in the phyloseq object (Supp. Table 6)
physeq_maas0 <- clean_tax_names_maaslin(physeq_maas0)

# Select the samples of one question (as in beta.R, plus Q3.5)
subset_physeq <- function(physeq, q){
  physeq_0 <- subset_samples(physeq, colSums(otu_table(physeq)) != 0 & !is.na(colSums(otu_table(physeq))))
  physeq_0 <- subset_samples(physeq_0, Sample_or_Control == "True sample" &
                               Isolation == "Direct_Isolation" &
                               LibrarySizeDecontam > 0)

  if(grepl("Q1|Q2", q)){
    # Patients with both lungs sampled
    study_nrs <- sample_data(physeq_0)$Study_Nr
    dup_studies <- study_nrs[duplicated(study_nrs) | duplicated(study_nrs, fromLast = TRUE)]
    physeq_sub <- prune_samples(sample_data(physeq_0)$Study_Nr %in% dup_studies, physeq_0)
    if (grepl("Q1.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    } else if (grepl("Q2", q)) {
      physeq_sub <- subset_samples(physeq_sub, Synchronous.tumor == "No")
      if (grepl("Q2.2", q)) physeq_sub <- subset_samples(physeq_sub, Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    }
  } else if (grepl("Q3", q)) {
    physeq_sub <- subset_samples(physeq_0, Lung == "Diseased" & Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    if (grepl("Q3.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    } else if (grepl("Q3.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    } else if (grepl("Q3.5", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis %in% c("Benign", "NSCLC"))
    }
  } else if (grepl("Q4", q)) {
    physeq_sub <- subset_samples(physeq_0, Lung == "Diseased" & Diagnosis == "NSCLC" &
                                   Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"))
    if (grepl("Q4.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    } else if (grepl("Q4.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    }
  } else if (grepl("Q5", q)) {
    physeq_sub <- subset_samples(physeq_0, Lung == "Diseased" & Diagnosis %in% c("Benign", "NSCLC", "SCLC") &
                                   !is.na(History.of.smoking.y.n))
    if (grepl("Q5.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis == "NSCLC")
    } else if (grepl("Q5.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis == "Benign")
    } else if (grepl("Q5.4", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis == "NSCLC" & Histology.NSCLC == "Adenocarcinoma")
    }
  } else if (grepl("Q6", q)) {
    physeq_sub <- subset_samples(physeq_0, Lung == "Diseased" & Diagnosis == "NSCLC")
    if (grepl("Q6.1", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(T))
    } else if (grepl("Q6.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(N))
    } else if (grepl("Q6.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(M))
      if (grepl("Q6.3.1", q)) {
        physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
      }
    }
  } else {
    physeq_sub <- physeq_0
  }

  prune_taxa(taxa_sums(physeq_sub) != 0, physeq_sub)
}

# Output folder of one run (one subfolder per reference level if there are several)
run_dir <- function(q, ref, level) {
  p <- params[[q]]
  if (length(p$refs) > 1) file.path(maaslin.dir, paste0("reference_", ref), paste0(p$name, "_", level))
  else file.path(maaslin.dir, paste0(p$name, "_", level))
}

#### End ####

#### 1. Full-sample MaAsLin2 per question ####

for (q in names(params)) {
  p <- params[[q]]
  physeq_sub <- subset_physeq(physeq_maas0, q)
  refs <- if (is.null(p$refs)) list(NULL) else as.list(p$refs)

  for (ref in refs) {
    for (level in tax_levels) {
      phy_input <- if (level == "ASV") physeq_sub else agglomerate_taxa(physeq_sub, level)
      out_dir <- run_dir(q, ifelse(is.null(ref), "", ref), level)

      res <- run_maaslin(phy_input, out_dir,
                         fixed_effects = c(p$fixed, p$covariates),
                         random_effects = p$rand,
                         reference = if (is.null(ref)) NULL else paste(p$fixed, ref, sep = ","))

      # Significant results with taxonomy
      sig <- res %>% filter(!is.na(qval), qval < max_sig)
      tax_df <- as.data.frame(tax_table(phy_input))
      tax_df$feature <- rownames(tax_df)
      write.xlsx(left_join(sig, tax_df, by = "feature"), file.path(out_dir, "significant_results.xlsx"))

      # Boxplots of the features significant for the variable of interest
      sig_fixed <- unique(sig$feature[sig$metadata == p$fixed])
      if (length(sig_fixed) > 0) {
        n_rows <- ceiling(length(sig_fixed) / 5)
        P <- plot_features(physeq_sub, sig_fixed, level, p$fixed, p$col,
                           legend_title = if (is.null(p$lab_plot)) p$fixed else p$lab_plot)
        ggsave(file.path(out_dir, "significant_results_boxplot.svg"), P,
               width = ifelse(length(sig_fixed) >= 5, 16, 4 + 3.2 * (length(sig_fixed) - 1)),
               height = 6 * n_rows, dpi = 300, limitsize = FALSE)
      }
    }
  }
}

#### End ####

#### 2a. Leave-one-out: Q3.1 diagnosis (Supplementary Table 5) ####

# One sample of the smaller groups (SCLC, Benign) left out at a time, each
# reference level. A feature x comparison is robust if it is significant in
# all iterations of a batch with the same direction, and present in > 5 % of
# the samples of at least one diagnosis

q <- "Q3.1"
p <- params[[q]]
loo.dir <- file.path(maaslin.dir, "Q3.1_LOO")
physeq_full <- subset_physeq(physeq_maas0, q)
loo_groups <- c("SCLC", "Benign")

prev_df <- prevalence_by_group(physeq_full, p$fixed)
tax_df <- as.data.frame(tax_table(physeq_full)); tax_df$feature <- rownames(tax_df)

q31_loo <- list()
for (ref in p$refs) {
  for (g in loo_groups) {
    loo_samples <- sample_names(physeq_full)[as.character(sample_data(physeq_full)[[p$fixed]]) == g]
    for (i in seq_along(loo_samples)) {
      s <- loo_samples[i]
      message("Q3.1 LOO | reference ", ref, " | batch ", g, " | ", i, "/", length(loo_samples), " | excluding ", s)
      physeq_iter <- prune_samples(sample_names(physeq_full) != s, physeq_full)
      physeq_iter <- prune_taxa(taxa_sums(physeq_iter) != 0, physeq_iter)
      res <- run_maaslin(physeq_iter, file.path(loo.dir, paste0("reference_", ref), g, paste0("iter", i, "_excl_", s)),
                         fixed_effects = c(p$fixed, p$covariates), reference = paste(p$fixed, ref, sep = ","),
                         plots = FALSE)
      q31_loo[[length(q31_loo) + 1]] <- res %>%
        filter(metadata == p$fixed) %>%
        mutate(reference = ref, loo_batch = g, n_iterations = length(loo_samples),
               iteration = i, excluded_sample = s)
    }
  }
}
q31_loo <- bind_rows(q31_loo)

# Full-sample results of the same comparisons
q31_full <- bind_rows(lapply(p$refs, function(ref) {
  read.delim(file.path(run_dir(q, ref, "ASV"), "all_results.tsv")) %>%
    filter(metadata == p$fixed) %>%
    transmute(reference = ref, feature, value, full_coef = coef, full_stderr = stderr, full_qval = qval)
}))

q31_summary <- q31_loo %>%
  filter(qval < max_sig) %>%
  group_by(reference, loo_batch, feature, value) %>%
  summarise(n_significant_iterations = n(),
            n_iterations = first(n_iterations),
            mean_qval = mean(qval), max_qval = max(qval), mean_coef = mean(coef),
            direction_consistent = length(unique(sign(coef))) == 1,
            .groups = "drop") %>%
  mutate(direction = ifelse(!direction_consistent, "inconsistent",
                            ifelse(mean_coef > 0, paste0("higher_in_", value), paste0("lower_in_", value))),
         pct_iterations_significant = n_significant_iterations / n_iterations) %>%
  left_join(q31_full, by = c("reference", "feature", "value")) %>%
  left_join(tax_df, by = "feature") %>%
  left_join(prev_df, by = "feature")

prev_cols <- grep("^prevalence_", colnames(q31_summary), value = TRUE)
q31_robust <- q31_summary %>%
  filter(n_significant_iterations == n_iterations, direction_consistent,
         if_any(all_of(prev_cols), ~ . > 0.05))

write.xlsx(list(robust_features = q31_robust,
                summary_all_features = q31_summary,
                all_iterations = q31_loo),
           file.path(loo.dir, "Q3.1_LOO_summary.xlsx"))

#### End ####

#### 2b. Leave-one-out: Q5.2 smoking in NSCLC (Supplementary Table 7, Figure 5C) ####

# One never-smoker left out at a time (8 iterations). An ASV is LOO-robust if
# it is significant in at least 6 of the 8 iterations.

q <- "Q5.2"
p <- params[[q]]
loo.dir <- file.path(maaslin.dir, "Q5.2_LOO")
physeq_full <- subset_physeq(physeq_maas0, q)
loo_samples <- sample_names(physeq_full)[as.character(sample_data(physeq_full)[[p$fixed]]) == "No"]
n_iterations <- length(loo_samples)
loo_min_iterations <- 6

prev_df <- prevalence_by_group(physeq_full, p$fixed)
tax_df <- as.data.frame(tax_table(physeq_full)); tax_df$feature <- rownames(tax_df)
tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

q52_loo <- list()
for (i in seq_along(loo_samples)) {
  s <- loo_samples[i]
  message("Q5.2 LOO | ", i, "/", n_iterations, " | excluding ", s)
  physeq_iter <- prune_samples(sample_names(physeq_full) != s, physeq_full)
  physeq_iter <- prune_taxa(taxa_sums(physeq_iter) != 0, physeq_iter)
  res <- run_maaslin(physeq_iter, file.path(loo.dir, paste0("iter", i, "_excl_", s)),
                     fixed_effects = c(p$fixed, p$covariates), reference = paste(p$fixed, "No", sep = ","),
                     plots = FALSE)
  q52_loo[[i]] <- res %>% filter(metadata == p$fixed) %>% mutate(iteration = i, excluded_sample = s)
}
q52_loo <- bind_rows(q52_loo)

q52_summary <- q52_loo %>%
  filter(qval < max_sig) %>%
  group_by(feature) %>%
  summarise(n_significant_iterations = n(),
            mean_qval = mean(qval), max_qval = max(qval), mean_coef = mean(coef),
            direction_consistent = length(unique(sign(coef))) == 1,
            .groups = "drop") %>%
  mutate(direction = ifelse(!direction_consistent, "inconsistent",
                            ifelse(mean_coef > 0, "higher_in_smoker", "higher_in_never_smoker")))

# Full-sample results of every ASV, with taxonomy and the LOO summary
q52_combined <- read.delim(file.path(run_dir(q, "No", "ASV"), "all_results.tsv")) %>%
  filter(metadata == p$fixed) %>%
  transmute(feature, full_coef = coef, full_stderr = stderr, full_pval = pval, full_qval = qval,
            full_significant = !is.na(qval) & qval < max_sig) %>%
  left_join(tax_df[, c("feature", tax_cols)], by = "feature") %>%
  left_join(q52_summary, by = "feature") %>%
  left_join(prev_df, by = "feature") %>%
  mutate(n_significant_iterations = coalesce(n_significant_iterations, 0L),
         robust_LOO = n_significant_iterations >= loo_min_iterations & coalesce(direction_consistent, FALSE),
         # Plot label: genus, or the next assigned rank up
         label = coalesce(Genus, Family, Order, Class, Phylum, Kingdom, "Unclassified")) %>%
  arrange(full_qval)

write.xlsx(list(significant_and_LOO_robust = q52_combined %>% filter(full_significant, robust_LOO),
                all_features = q52_combined,
                all_iterations = q52_loo),
           file.path(loo.dir, "Q5.2_LOO_summary.xlsx"))

## Figure 5C: volcano plot of the full-sample results
q52_combined <- q52_combined %>%
  mutate(point_group = factor(case_when(
    full_significant & robust_LOO ~ paste0("LOO-robust (≥", loo_min_iterations, "/", n_iterations, ")"),
    full_significant ~ "Significant (all samples)",
    TRUE ~ "Not significant"),
    levels = c("Not significant", "Significant (all samples)",
               paste0("LOO-robust (≥", loo_min_iterations, "/", n_iterations, ")"))))

p_volcano <- ggplot(q52_combined, aes(x = full_coef, y = -log10(full_qval), color = point_group)) +
  geom_point(alpha = 0.75, size = 3) +
  scale_color_manual(values = setNames(c("grey70", "steelblue", "firebrick"), levels(q52_combined$point_group)),
                     name = NULL) +
  geom_hline(yintercept = -log10(max_sig), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_text_repel(data = filter(q52_combined, full_significant), aes(label = label),
                  size = 6, color = "black", max.overlaps = Inf, seed = 123) +
  labs(x = "MaAsLin2 coefficient (ever vs. never smokers)", y = expression(-log[10]("q-value"))) +
  theme_bw(base_size = 30) +
  theme(legend.position = "bottom")
ggsave(file.path(loo.dir, "Fig5C_volcano_full_vs_LOO.svg"), p_volcano, width = 10, height = 10, dpi = 300)

#### End ####

#### 3. Figures built on the MaAsLin2 results ####

fig.dir <- file.path(maaslin.dir, "figures")
dir.create(fig.dir, recursive = T, showWarnings = F)

# Forest plot of MaAsLin2 coefficients with 95 % intervals (coef +- 1.96 SE)
forest_plot <- function(df, title) {
  ggplot(df, aes(x = coef, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = coef - 1.96 * stderr, xmax = coef + 1.96 * stderr), height = 0.2) +
    labs(title = title, x = "MaAsLin2 coefficient with 95% CI", y = NULL) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), text = element_text(size = 15),
          plot.title = element_text(size = 15))
}

## Figure 3D: ASV24 (Streptococcus mitis) by N stage, N0 as reference
asv24 <- read.delim(file.path(run_dir("Q6.2", "", "ASV"), "all_results.tsv")) %>%
  filter(feature == "ASV24", metadata == "N") %>%
  transmute(label = paste0("N", value), coef, stderr, qval) %>%
  bind_rows(data.frame(label = "N0 (ref)", coef = 0, stderr = 0, qval = NA), .) %>%
  mutate(label = factor(label, levels = c("N0 (ref)", "N1", "N2", "N3")))
write.xlsx(asv24, file.path(fig.dir, "Fig3D_ASV24_N.xlsx"))
asv24_forest <- forest_plot(asv24, "Association between N and ASV24 abundance")
ggsave(file.path(fig.dir, "Fig3D_ASV24_N_forest.svg"),
       asv24_forest, width = 6, height = 4, dpi = 300)

## Figure 5D: never-smoker ASVs (Figure 5C, significant and LOO-robust, higher
## in never-smokers) in the NSCLC vs benign comparison (Q3.5)
never_smoker_asvs <- q52_combined %>% filter(full_significant, robust_LOO, full_coef < 0)
fig5d <- read.delim(file.path(run_dir("Q3.5", "Benign", "ASV"), "all_results.tsv")) %>%
  filter(feature %in% never_smoker_asvs$feature, metadata == "Diagnosis") %>%
  left_join(never_smoker_asvs %>% select(feature, label), by = "feature") %>%
  mutate(label = gsub("_", " ", paste0(label, " (", feature, ")")))
fig5d$label <- factor(fig5d$label, levels = fig5d$label[order(fig5d$coef)])
write.xlsx(fig5d, file.path(fig.dir, "Fig5D_never_smoker_ASVs_NSCLC_vs_benign.xlsx"))
fig5d_plot <- forest_plot(fig5d, "ASV association with Diagnosis (NSCLC vs Benign)") +
  theme(text = element_text(size = 20), axis.text.y = element_text(size = 20))
ggsave(file.path(fig.dir, "Fig5D_never_smoker_ASVs_NSCLC_vs_benign_forest.svg"),
       fig5d_plot, width = 15, height = 10, dpi = 300)

#### End ####

#### 4. ASVs per diagnosis in smokers (Supplementary Figure 4D, Supplementary Table 6) ####

# Relative abundances below 0.05 % are set to zero;
# ASVs are kept if present in more than 2 % of the samples of at least one
# diagnosis. An ASV belongs to a diagnosis if present in at least one of its
# samples.

venn_abundance_cutoff <- 0.0005
venn_prevalence_cutoff <- 0.02

physeq_venn <- subset_physeq(physeq_maas0, "Q3.3")
otu_venn <- as(otu_table(physeq_venn), "matrix")
if (!taxa_are_rows(physeq_venn)) otu_venn <- t(otu_venn)
otu_venn[otu_venn < venn_abundance_cutoff] <- 0

dx_venn <- as.character(sample_data(physeq_venn)$Diagnosis)
dx_groups <- sort(unique(dx_venn))
prev_venn <- sapply(dx_groups, function(g) rowMeans(otu_venn[, dx_venn == g, drop = FALSE] > 0))
otu_venn <- otu_venn[apply(prev_venn, 1, max) > venn_prevalence_cutoff, ]
prev_venn <- prev_venn[rownames(otu_venn), ]

asv_sets <- lapply(setNames(dx_groups, dx_groups), function(g) rownames(otu_venn)[rowSums(otu_venn[, dx_venn == g, drop = FALSE]) > 0])

# Set names with the number of samples
venn_sets <- setNames(asv_sets, paste0(dx_groups, "\nn = ", table(dx_venn)[dx_groups]))
p_venn <- ggvenn(venn_sets, fill_color = unname(diagnosis_col[dx_groups]), stroke_size = 1,
                 digits = 0, set_name_size = 6, text_size = 6)
ggsave(file.path(fig.dir, "SuppFig4D_ASVs_smokers_venn.svg"), p_venn, width = 8, height = 8, dpi = 300)

# Supplementary Table 6: ASVs found in only one diagnosis, with taxonomy,
# prevalence and mean relative abundance (%) in that diagnosis
unique_asvs <- bind_rows(lapply(dx_groups, function(g) {
  asvs <- setdiff(asv_sets[[g]], unlist(asv_sets[setdiff(dx_groups, g)]))
  if (length(asvs) == 0) return(NULL)
  data.frame(ASV = asvs, Diagnosis = g, tax_original[asvs, ], check.names = FALSE,
             Prevalence = prev_venn[asvs, g],
             MeanAbundance = 100 * rowMeans(otu_venn[asvs, dx_venn == g, drop = FALSE]))
}))
write.xlsx(unique_asvs, file.path(fig.dir, "SuppTable6_unique_ASVs_smokers.xlsx"))

#### End ####
