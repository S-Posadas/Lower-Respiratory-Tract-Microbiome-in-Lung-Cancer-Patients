# LEfSe (lefser) #
# Directly isolated samples, batch-corrected relative abundances.
# For every question with two groups:
# Kruskal-Wallis test and LDA effect size at phylum, class, order, family and
# genus level (lefser, modified to also return the p-values:
# lefser_functions.R); thresholds p < 0.05 and |LDA score| >= 2.
# Supp. Figures 8E, 10D, 11B, 13B.

#### Package setup ####

Sys.setenv(language = "EN")

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(SummarizedExperiment)
  library(lefser)
  library(openxlsx)
})

source("/home/nanoteam/projects/LML/scripts/lefser_functions.R")

theme_set(theme_bw())

#### End ####

#### Load batch-corrected phyloseq object (2.MMUPHin.R) ####

res.dir <- "Results_direct"

rdata <- new.env()
load(file.path(res.dir, "RData", "6.phyloseq.corrected.RData"), envir = rdata)
physeq_lefse0 <- rdata$physeq_re_bc
lung_col <- rdata$lung_col
histology_col <- rdata$histology_col
smoker_col <- rdata$smoker_col
rm(rdata)

#### End ####

#### Settings and directories ####

lefse.dir <- file.path(res.dir, "lefse")
dir.create(lefse.dir, recursive = T, showWarnings = F)

#### End ####

#### Functions ####

# Taxonomy names as lefser needs them: no spaces or dots, "Incertae Sedis"
# families named after their order, no Kingdom
clean_tax_names_lefser <- function(physeq) {
  tax <- as.data.frame(tax_table(physeq))
  tax[] <- lapply(tax, gsub, pattern = " ", replacement = "_")
  tax[] <- lapply(tax, gsub, pattern = "\\.", replacement = "_")
  tax$Species[grep("_spp.", tax$Species)] <- NA
  tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis", "Family"] <- paste(
    tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis", "Order"], "Incertae_Sedis", sep = "_")
  tax$Kingdom <- NULL
  tax_table(physeq) <- tax_table(as.matrix(tax))
  physeq
}

# LEfSe for one question: marker table and plots (top 20 / 50 markers by LDA
# score and by p-value; LDA bar plot and dot plot)
run_lefser_analysis <- function(physeq, group_var, output_dir, col_palette) {
  groups <- sample_data(physeq)[[group_var]]
  # Colours of the two groups compared (lefserPlot needs exactly two)
  if (!is.null(names(col_palette))) col_palette <- col_palette[levels(droplevels(factor(groups)))]
  if (any(table(groups) < 2)) {
    writeLines("All groups must have at least 2 samples", file.path(output_dir, "my_analysis_log.txt"))
    return(NULL)
  }

  se <- SummarizedExperiment(assays = list(counts = unclass(otu_table(physeq))),
                             colData = as(sample_data(physeq), "data.frame"),
                             rowData = tax_table(physeq)[, c("Phylum", "Class", "Order", "Family", "Genus")])
  se <- lefser::relativeAb(se)

  marker <- tryCatch(my_lefserClades(se, classCol = group_var),
                     error = function(e) { message("LEfSe error: ", conditionMessage(e)); NULL })
  if (is.null(nrow(marker)) || nrow(marker) == 0) {
    writeLines(paste("No significant markers for group:", group_var), file.path(output_dir, "my_analysis_log.txt"))
    return(NULL)
  }

  write.xlsx(marker, file.path(output_dir, paste("marker_table", nrow(marker), "m.xlsx", sep = "_")), rowNames = TRUE)

  create_plots <- function(res_sub, prefix) {
    try({
      h <- nrow(res_sub) * 0.5
      p_val <- my_lefserPlot(res_sub, colors = col_palette, label.font.size = 6, other.font.size = 20)
      ggsave(file.path(output_dir, paste0(prefix, "_dot.svg")), p_val, width = 16, height = h, dpi = 300)
      p_lda <- lefser::lefserPlot(res_sub, colors = unname(col_palette), label.font.size = 6)
      ggsave(file.path(output_dir, paste0(prefix, "_lda_bar.svg")), p_lda, width = 16, height = h, dpi = 300)
    })
  }

  for (n in c(20, 50)) {
    by_lda <- marker[order(abs(marker$scores), decreasing = TRUE), ]
    if (nrow(by_lda) > n) create_plots(by_lda[1:n, ], paste0("top", n, "lda")) else create_plots(by_lda, "all_markers")
    by_p <- marker[order(marker$kw_p.value), ]
    if (nrow(by_p) > n) create_plots(by_p[1:n, ], paste0("top", n, "pval")) else create_plots(by_p, "all_markers")
  }
  marker
}

#### End ####

#### Questions with two groups ####

params <- list(
  Q1.1 = list(name = "Q1.1.all_directly_isolated_lung", fixed = "Lung", col = lung_col),
  Q1.2 = list(name = "Q1.2.main-D_directly_isolated_lung", fixed = "Lung", col = lung_col),
  Q2.1 = list(name = "Q2.1.all_directly_isolated_non_synchronous_lung", fixed = "Lung", col = lung_col),
  Q2.2 = list(name = "Q2.2.main-D_directly_isolated_non_synchronous_lung", fixed = "Lung", col = lung_col),
  Q4.1 = list(name = "Q4.1.histology", fixed = "Histology.NSCLC", col = histology_col),
  Q4.2 = list(name = "Q4.2.histology_nonsmokers", fixed = "Histology.NSCLC", col = histology_col),
  Q4.3 = list(name = "Q4.3.histology_smokers", fixed = "Histology.NSCLC", col = histology_col),
  Q5.1 = list(name = "Q5.1.main-D_smoking", fixed = "History.of.smoking.y.n", col = smoker_col),
  Q5.2 = list(name = "Q5.2.NSCLC_smoking", fixed = "History.of.smoking.y.n", col = smoker_col),
  Q5.3 = list(name = "Q5.3.Benign_smoking", fixed = "History.of.smoking.y.n", col = smoker_col),
  Q5.4 = list(name = "Q5.4.Adenocarcinoma_smoking", fixed = "History.of.smoking.y.n", col = smoker_col),
  Q6.2.1 = list(name = "Q6.2.NSCLC_N_twolevels", fixed = "Ntwo", col = c(N0 = "lightblue", `N1-3` = "firebrick")),
  Q6.3 = list(name = "Q6.3.NSCLC_M", fixed = "M", col = c(`0` = "lightblue", `1` = "firebrick")),
  Q6.3.1 = list(name = "Q6.3.1.NSCLC_M_smokers", fixed = "M", col = c(`0` = "lightblue", `1` = "firebrick"))
)

# N stage in two levels: N0 vs N1-3
sample_data(physeq_lefse0)$Ntwo <- factor(ifelse(sample_data(physeq_lefse0)$N == 0, "N0", "N1-3"), levels = c("N0", "N1-3"))

physeq_lefse0 <- clean_tax_names_lefser(physeq_lefse0)

# Select the samples of one question (as in beta.R)
subset_physeq <- function(physeq, q){
  cs <- colSums(otu_table(physeq))
  physeq_0 <- prune_samples(!is.na(cs) & cs != 0, physeq)
  sd0 <- data.frame(sample_data(physeq_0))
  physeq_0 <- prune_samples(sd0$Sample_or_Control == "True sample" & sd0$Isolation == "Direct_Isolation" &
                              sd0$LibrarySizeDecontam > 0, physeq_0)
  sd <- data.frame(sample_data(physeq_0))
  smoking <- sd$History.of.smoking.y.n
  diseased_nsclc <- sd$Lung == "Diseased" & sd$Diagnosis == "NSCLC"

  keep <- if (grepl("Q1|Q2", q)) {
    # Patients with both lungs sampled
    dup <- sd$Study_Nr[duplicated(sd$Study_Nr)]
    k <- sd$Study_Nr %in% dup
    if (grepl("Q1.2", q)) k <- k & sd$Diagnosis %in% c("Benign", "NSCLC", "SCLC")
    if (grepl("Q2", q)) k <- k & sd$Synchronous.tumor %in% "No"
    if (grepl("Q2.2", q)) k <- k & sd$Diagnosis %in% c("Benign", "NSCLC", "SCLC")
    k
  } else if (grepl("Q4", q)) {
    k <- diseased_nsclc & sd$Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma")
    if (grepl("Q4.2", q)) k <- k & smoking %in% "No"
    if (grepl("Q4.3", q)) k <- k & smoking %in% "Yes"
    k
  } else if (grepl("Q5", q)) {
    k <- sd$Lung == "Diseased" & sd$Diagnosis %in% c("Benign", "NSCLC", "SCLC") & !is.na(smoking)
    if (grepl("Q5.2", q)) k <- k & sd$Diagnosis == "NSCLC"
    if (grepl("Q5.3", q)) k <- k & sd$Diagnosis == "Benign"
    if (grepl("Q5.4", q)) k <- k & sd$Diagnosis == "NSCLC" & sd$Histology.NSCLC %in% "Adenocarcinoma"
    k
  } else if (grepl("Q6.2", q)) {
    diseased_nsclc & !is.na(sd$N)
  } else if (grepl("Q6.3", q)) {
    k <- diseased_nsclc & !is.na(sd$M)
    if (grepl("Q6.3.1", q)) k <- k & smoking %in% "Yes"
    k
  }

  physeq_sub <- prune_samples(keep, physeq_0)
  prune_taxa(taxa_sums(physeq_sub) != 0, physeq_sub)
}

#### End ####

#### LEfSe per question ####

for (q in names(params)) {
  p <- params[[q]]
  message(q)
  physeq_q <- subset_physeq(physeq_lefse0, q)
  out_dir <- file.path(lefse.dir, p$name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_lefser_analysis(physeq_q, group_var = p$fixed, output_dir = out_dir, col_palette = p$col)
}

#### End ####
