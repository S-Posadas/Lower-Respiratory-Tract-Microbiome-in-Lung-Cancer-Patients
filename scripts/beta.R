# Beta diversity #
# Bray-Curtis, weighted UniFrac and unweighted UniFrac distances on the
# batch-corrected relative abundances: NMDS and PCoA plots, PERMANOVA (adonis2)
# and a betadisper check per question. The paper reports weighted UniFrac.

#### Package setup ####

Sys.setenv(language = "EN")

library(dplyr)
library(ggplot2)
library(vegan)
library(phyloseq)
library(gridExtra)  # arrangeGrob()
library(cowplot)    # get_plot_component()
library(openxlsx)

theme_set(theme_bw())

#### End ####

#### Load batch-corrected phyloseq object (2.MMUPHin.R) ####

res.dir <- "Results_direct"

# Load into a separate environment and keep only what is needed
rdata <- new.env()
load(file.path(res.dir, "RData", "6.phyloseq.corrected.RData"), envir = rdata)
physeq_beta0 <- rdata$physeq_re_bc  # batch-corrected relative abundances, true samples
lung_col <- rdata$lung_col
diagnosis_col <- rdata$diagnosis_col
histology_col <- rdata$histology_col
smoker_col <- rdata$smoker_col
rm(rdata)

#### End ####

#### Settings and directories ####

b.stats <- file.path(res.dir, "beta_stats")
b.plots <- file.path(res.dir, "beta_plots")

dir.create(b.stats, recursive = T, showWarnings = F)
dir.create(b.plots, recursive = T, showWarnings = F)

n_perm <- 99999      # PERMANOVA permutations
n_perm_disp <- 9999  # betadisper permutations
n_cores <- 10        # permutations run in parallel; results do not depend on it

distances <- c(bray = "bray_dist", wUF = "wunifrac_dist", uwUF = "uunifrac_dist")

#### End ####

#### Plot theme settings ####

text_size = 30
thm.beta <- theme(text = element_text(size = text_size),
                  legend.position = "bottom")

#### End ####

#### Functions ####

# Bray-Curtis, weighted and unweighted UniFrac distances, NMDS and PCoA
beta_analysis <- function(physeq){
  # set.seed: with an unrooted tree phyloseq roots it on a random tip
  set.seed(123)
  bray_dist <- phyloseq::distance(physeq, method = "bray")
  set.seed(123)
  wunifrac_dist <- phyloseq::distance(physeq, method = "wunifrac")
  set.seed(123)
  uunifrac_dist <- phyloseq::distance(physeq, method = "uunifrac")
  
  set.seed(123)
  nmds_bray <- vegan::metaMDS(bray_dist, k = 2, trymax = 100, trace = FALSE)
  set.seed(123)
  nmds_wunifrac <- vegan::metaMDS(wunifrac_dist, k = 2, trymax = 100, trace = FALSE)
  set.seed(123)
  nmds_uunifrac <- vegan::metaMDS(uunifrac_dist, k = 2, trymax = 100, trace = FALSE)
  
  # PCoA of the same distance matrices
  pcoa_bray <- phyloseq::ordinate(physeq, method = "PCoA", distance = bray_dist)
  pcoa_wunifrac <- phyloseq::ordinate(physeq, method = "PCoA", distance = wunifrac_dist)
  pcoa_uunifrac <- phyloseq::ordinate(physeq, method = "PCoA", distance = uunifrac_dist)
  
  list(bray_dist = bray_dist, wunifrac_dist = wunifrac_dist, uunifrac_dist = uunifrac_dist,
       nmds_bray = nmds_bray, nmds_wunifrac = nmds_wunifrac, nmds_uunifrac = nmds_uunifrac,
       pcoa_bray = pcoa_bray, pcoa_wunifrac = pcoa_wunifrac, pcoa_uunifrac = pcoa_uunifrac)
}

# NMDS and PCoA side by side, with a shared legend below
plot_nmds_pcoa <- function(phy, beta, dist, var, palette = NULL, legend_title = NULL){
  # dist: "bray_dist", "wunifrac_dist" or "uunifrac_dist"
  NM <- beta[[paste0("nmds_", sub("_dist$", "", dist))]]
  PC <- beta[[paste0("pcoa_", sub("_dist$", "", dist))]]

  NMDS = plot_ordination(phy, NM, color = var) + geom_point(size = 8, alpha = 0.5) +
    stat_ellipse(type = "t") + thm.beta +
    annotate("text", x = max(scores(NM)[,"NMDS1"]), y = max(scores(NM)[,"NMDS2"]),
             label = paste("Stress:", round(NM$stress, 3)), hjust = 1, vjust = 1, size = 7)

  PCOA = plot_ordination(phy, PC, color = var) + geom_point(size = 8, alpha = 0.5) +
    stat_ellipse(type = "t") + thm.beta

  if (!is.null(palette)) {
    NMDS <- NMDS + scale_color_manual(values = palette)
    PCOA <- PCOA + scale_color_manual(values = palette)
  }

  if (!is.null(legend_title)) {
    NMDS <- NMDS + labs(color = legend_title)
    PCOA <- PCOA + labs(color = legend_title)
  }

  NP <- arrangeGrob(NMDS + ggtitle("NMDS") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)),
                    PCOA + ggtitle("PCoA") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)),
                    nrow = 1)

  arrangeGrob(NP, get_plot_component(PCOA, pattern = "guide-box-bottom"), nrow = 2, heights = c(10,1))
}

# PERMANOVA (adonis2, marginal tests) for one question. Samples with a missing
# value in any model variable are left out (complete cases).
run_permanova <- function(physeq, dist_obj, dist_name, fixed_var, covariates = NULL, strata_var = NULL) {
  formula_rhs <- paste(c(fixed_var, covariates), collapse = " + ")

  meta <- data.frame(sample_data(physeq))
  vars <- c(fixed_var, covariates, strata_var)
  meta_sub <- meta[complete.cases(meta[, vars, drop = FALSE]), , drop = FALSE]
  dist_sub <- as.dist(as.matrix(dist_obj)[rownames(meta_sub), rownames(meta_sub)])
  strata_vals <- if (!is.null(strata_var)) meta_sub[[strata_var]] else NULL

  set.seed(123)
  res <- adonis2(as.formula(paste("dist_sub ~", formula_rhs)),
                 data = meta_sub, permutations = n_perm, by = "margin",
                 strata = strata_vals, parallel = n_cores)

  res_df <- data.frame(Distance = dist_name,
                       Term = rownames(res),
                       Formula = formula_rhs,
                       Strata = ifelse(is.null(strata_var), "none", strata_var),
                       N_samples = nrow(meta_sub),
                       N_excluded_NA = nrow(meta) - nrow(meta_sub),
                       as.data.frame(res), check.names = FALSE)
  rownames(res_df) <- NULL

  list(result = res_df, data = meta_sub, dist = dist_sub)
}

# betadisper + permutest on one grouping variable, on the same samples as the
# PERMANOVA. A significant result means the groups differ in spread, which
# PERMANOVA can mistake for a difference in location.
check_betadisper <- function(dist_obj, data, factor_var, dist_name) {
  grp <- droplevels(factor(data[[factor_var]]))
  if (nlevels(grp) < 2) {
    return(data.frame(Distance = dist_name, Betadisp_Factor = factor_var,
                      Betadisp_F = NA, Betadisp_p = NA,
                      Dispersion_Warning = "N/A - only one group level"))
  }

  bd <- betadisper(dist_obj, grp)
  set.seed(123)
  pt <- permutest(bd, permutations = n_perm_disp, parallel = n_cores)

  p_val <- pt$tab$`Pr(>F)`[1]
  warn_text <- if (!is.na(p_val) && p_val < 0.05) {
    "CAUTION: significant dispersion heterogeneity (p < 0.05) - PERMANOVA result may partly reflect group spread rather than centroid differences."
  } else {
    "OK - no significant dispersion heterogeneity detected"
  }

  data.frame(Distance = dist_name, Betadisp_Factor = factor_var,
             Betadisp_F = round(pt$tab$F[1], 3), Betadisp_p = round(p_val, 5),
             Dispersion_Warning = warn_text)
}

#### End ####

#### Questions to answer ####

#' Questions
#' Q1: Is there a difference between Diseased vs the parallel Normal Lung?
#' Q2: Is there a difference between Diseased vs the parallel Normal Lung without samples with synchronous tumor?
#' Q3.1: What are the differences between NSCLC vs SCLC vs Benign Tumor?
#' Q3.2/Q3.3: The same in nonsmokers / smokers
#' Q4.1: What are the differences between NSCLC main histologies (Adenocarcinoma vs Squamous)?
#' Q4.2/Q4.3: The same in nonsmokers / smokers
#' Q5: Does history of smoking have an impact on the lung microbiome?
#' Q6.1: Does the T stage have an impact on the lung microbiome?
#' Q6.2: Does the N stage have an impact on the lung microbiome?
#' Q6.3: Does the M stage have an impact on the lung microbiome?

#' Covariate policy (same as alpha.R):
#' - Q1, Q2 (paired Diseased vs Normal lung, same patient): no covariates;
#'   permutations are restricted within patient (strata = Study_Nr).
#' - Q3.2, Q3.3, Q4.2, Q4.3, all Q5, Q6.3.1: adjust for Age + Sex only (smoking
#'   status is either the stratifying variable or the variable of interest).
#' - Everywhere else (Q3.1, Q4.1, Q6.1, Q6.2, Q6.3): adjust for
#'   History.of.smoking.y.n + Age + Sex.

smoking_age_sex <- c("History.of.smoking.y.n", "Age", "Sex")
age_sex <- c("Age", "Sex")

params <- list(
  Q1.1 = list(name = "Q1.1.all_directly_isolated_lung", fixed = "Lung", col = lung_col, strata = "Study_Nr"),
  Q1.2 = list(name = "Q1.2.main-D_directly_isolated_lung", fixed = "Lung", col = lung_col, strata = "Study_Nr"),
  Q2.1 = list(name = "Q2.1.all_directly_isolated_non_synchronous_lung", fixed = "Lung", col = lung_col, strata = "Study_Nr"),
  Q2.2 = list(name = "Q2.2.main-D_directly_isolated_non_synchronous_lung", fixed = "Lung", col = lung_col, strata = "Study_Nr"),
  Q3.1 = list(name = "Q3.1.diagnosis", fixed = "Diagnosis", col = diagnosis_col, covariates = smoking_age_sex),
  Q3.2 = list(name = "Q3.2.diagnosis_nonsmokers", fixed = "Diagnosis", col = diagnosis_col, covariates = age_sex),
  Q3.3 = list(name = "Q3.3.diagnosis_smokers", fixed = "Diagnosis", col = diagnosis_col, covariates = age_sex),
  Q4.1 = list(name = "Q4.1.histology", fixed = "Histology.NSCLC", col = histology_col, lab_plot = "Histology", covariates = smoking_age_sex),
  Q4.2 = list(name = "Q4.2.histology_nonsmokers", fixed = "Histology.NSCLC", col = histology_col, lab_plot = "Histology", covariates = age_sex),
  Q4.3 = list(name = "Q4.3.histology_smokers", fixed = "Histology.NSCLC", col = histology_col, lab_plot = "Histology", covariates = age_sex),
  Q5.1 = list(name = "Q5.1.main-D_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking", covariates = age_sex),
  Q5.2 = list(name = "Q5.2.NSCLC_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking", covariates = age_sex),
  Q5.3 = list(name = "Q5.3.Benign_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking", covariates = age_sex),
  Q5.4 = list(name = "Q5.4.Adenocarcinoma_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking", covariates = age_sex),
  Q6.1 = list(name = "Q6.1.NSCLC_T", fixed = "T", col = scales::hue_pal()(5)[c(4,2,5,1)], lab_plot = "T", covariates = smoking_age_sex),
  Q6.2 = list(name = "Q6.2.NSCLC_N", fixed = "N", col = scales::hue_pal()(5)[c(3:4,2,5,1)], covariates = smoking_age_sex),
  Q6.2.1 = list(name = "Q6.2.NSCLC_N_twolevels", fixed = "Ntwo", col = scales::hue_pal()(5)[c(3:4,2,5,1)], lab_plot = "N", covariates = smoking_age_sex),
  Q6.3 = list(name = "Q6.3.NSCLC_M", fixed = "M", col = c("lightblue", "firebrick"), covariates = smoking_age_sex),
  Q6.3.1 = list(name = "Q6.3.1.NSCLC_M_smokers", fixed = "M", col = c("lightblue", "firebrick"), covariates = age_sex)
)

# N stage in two levels: N0 vs N1-3
sample_data(physeq_beta0)$Ntwo <- factor(ifelse(sample_data(physeq_beta0)$N == 0, "N0", "N1-3"), levels = c("N0", "N1-3"))

# Select the samples of one question
subset_physeq <- function(physeq, q){
  # Remove samples with no reads
  physeq_0 <- subset_samples(physeq, colSums(otu_table(physeq)) != 0 & !is.na(colSums(otu_table(physeq))))

  # Directly isolated true samples only
  physeq_0 <- physeq_0 %>%
    subset_samples(Sample_or_Control == "True sample" &
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
      if (grepl("Q2.2", q)) {
        physeq_sub <- subset_samples(physeq_sub, Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
      }
    }

  } else if (grepl("Q3", q)) {
    physeq_sub <- subset_samples(physeq_0,
                                 Lung == "Diseased" &
                                   Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    if (grepl("Q3.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    } else if (grepl("Q3.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    }

  } else if (grepl("Q4", q)) {
    physeq_sub <- subset_samples(physeq_0,
                                 Lung == "Diseased" &
                                   Diagnosis == "NSCLC" &
                                   Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"))
    if (grepl("Q4.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    } else if (grepl("Q4.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    }

  } else if (grepl("Q5", q)) {
    physeq_sub <- subset_samples(physeq_0,
                                 Lung == "Diseased" &
                                   Diagnosis %in% c("Benign", "NSCLC", "SCLC") &
                                   !is.na(History.of.smoking.y.n))
    if (grepl("Q5.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis == "NSCLC")
    } else if (grepl("Q5.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis == "Benign")
    } else if (grepl("Q5.4", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis == "NSCLC" & Histology.NSCLC == "Adenocarcinoma")
    }

  } else if (grepl("Q6", q)) {
    physeq_sub <- subset_samples(physeq_0,
                                 Lung == "Diseased" &
                                   Diagnosis == "NSCLC")
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

#### End ####

#### Beta diversity: ordination plots and PERMANOVA ####

for (q in names(params)) {

  print(q)
  p <- params[[q]]

  physeq_beta <- subset_physeq(physeq_beta0, q)
  beta_list   <- beta_analysis(physeq_beta)

  pm_all <- list()
  disp_all <- list()

  for (d in names(distances)) {
    dist_name <- distances[[d]]

    # NMDS and PCoA plot
    NP <- plot_nmds_pcoa(physeq_beta, beta_list, dist = dist_name, var = p$fixed,
                         palette = p$col, legend_title = p$lab_plot)
    ggsave(plot = NP, file.path(b.plots, paste0(p$name, "_", d, ".svg")), width = 15, height = 8, dpi = 300)

    # PERMANOVA
    pm <- run_permanova(physeq_beta, beta_list[[dist_name]], dist_name, p$fixed,
                        covariates = p$covariates, strata_var = p$strata)
    pm_all[[d]] <- pm$result

    # betadisper on the variable of interest and, if it is a covariate, on
    # smoking status (the covariate with missing values)
    disp_all[[paste0(d, "_fixed")]] <- check_betadisper(pm$dist, pm$data, p$fixed, dist_name)
    if ("History.of.smoking.y.n" %in% p$covariates) {
      disp_all[[paste0(d, "_smoking")]] <- check_betadisper(pm$dist, pm$data, "History.of.smoking.y.n", dist_name)
    }
  }

  pm_all <- do.call(rbind, pm_all)
  disp <- do.call(rbind, disp_all)
  rownames(pm_all) <- NULL
  rownames(disp) <- NULL
  print(pm_all[pm_all$Term == p$fixed, ])

  # Excel file: PERMANOVA, betadisper (warnings highlighted in red), samples tested
  # (the same samples for all three distances)
  wb <- createWorkbook()
  for (sheet in c("PERMANOVA", "Betadispersion", "Data")) addWorksheet(wb, sheet)
  writeData(wb, "PERMANOVA", pm_all)
  writeData(wb, "Betadispersion", disp)
  writeData(wb, "Data", pm$data, rowNames = TRUE)

  warn_rows <- which(grepl("^CAUTION", disp$Dispersion_Warning)) + 1  # +1 for the header row
  if (length(warn_rows) > 0) {
    addStyle(wb, "Betadispersion", style = createStyle(fgFill = "#FFC7CE", fontColour = "#9C0006"),
             rows = warn_rows, cols = 1:ncol(disp), gridExpand = TRUE)
  }
  setColWidths(wb, "PERMANOVA", cols = 1:ncol(pm_all), widths = "auto")
  setColWidths(wb, "Betadispersion", cols = 1:ncol(disp), widths = "auto")

  saveWorkbook(wb, file.path(b.stats, paste0(p$name, ".xlsx")), overwrite = TRUE)
}

#### End ####
