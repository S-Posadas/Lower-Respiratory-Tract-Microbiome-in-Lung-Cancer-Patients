# Culture enrichment shift and clinical variables #
# Supplementary Table 4: does the change in abundance after culture enrichment
# depend on N stage, M stage or history of smoking?
#
# For every patient with a directly isolated and a culture-enriched sample of
# the tumour-bearing lung: genus-level relative abundances, CLR transformation
# (pseudocount 1e-6) per sample, and delta = CLR(culture) - CLR(direct) per
# genus (positive: higher after culture). Genera present in >= 10% of the
# patients (in either sample type) are tested. Three sheets:
# 1. Wilcoxon_rank_sum-Kruskal_delta: delta per genus vs N (Kruskal-Wallis),
#    M and smoking (Wilcoxon rank-sum), BH correction per variable
# 2. MaAsLin2_delta: linear model per genus, delta ~ N + M + smoking (no
#    normalization or transformation; N and M as numbers, as published)
# 3. PERMANOVA_delta: Euclidean distance of the delta profiles,
#    ~ N + M + smoking (marginal effects), patients with all three variables

#### Package setup ####

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Maaslin2)
  library(openxlsx)
})

#### End ####

#### Load merged filtered phyloseq object (merge_phyloseq.R) ####

merged.dir <- "Results_merged"
physeq_delta <- readRDS(file.path(merged.dir, "RData", "physeq_filtered.rds"))

delta.dir <- file.path(merged.dir, "delta")
dir.create(delta.dir, recursive = T, showWarnings = F)

#### End ####

#### Settings ####

meta_vars <- c("N", "M", "History.of.smoking.y.n")
pseudocount <- 1e-6           # added to relative abundances before the CLR
prevalence_threshold <- 0.10  # genus present in >= 10% of the patients

#### End ####

#### Paired samples ####

# Tumour-bearing lung, true samples with reads
lib_size <- colSums(otu_table(physeq_delta))
sd <- data.frame(sample_data(physeq_delta))
physeq_delta <- prune_samples(sd$Sample_or_Control == "True sample" & sd$Lung == "Diseased" &
                                !is.na(lib_size) & lib_size != 0, physeq_delta)

# Genus level, genera as taxa names
physeq_delta <- tax_glom(physeq_delta, "Genus")
taxa_names(physeq_delta) <- as.character(tax_table(physeq_delta)[, "Genus"])

counts <- as(otu_table(physeq_delta), "matrix")
if (!taxa_are_rows(physeq_delta)) counts <- t(counts)
sd <- data.frame(sample_data(physeq_delta))

# One direct and one culture-enriched sample per patient
direct  <- sd[sd$Isolation == "Direct_Isolation", ]
culture <- sd[sd$Isolation == "Culture_Enriched", ]
direct$Study_Nr <- as.character(direct$Study_Nr)
culture$Study_Nr <- as.character(culture$Study_Nr)
patients <- sort(intersect(direct$Study_Nr, culture$Study_Nr))
stopifnot(!anyDuplicated(direct$Study_Nr), !anyDuplicated(culture$Study_Nr))

counts_direct  <- counts[, rownames(direct)[match(patients, direct$Study_Nr)]]
counts_culture <- counts[, rownames(culture)[match(patients, culture$Study_Nr)]]
colnames(counts_direct) <- colnames(counts_culture) <- patients

# Patients with no genus-level reads in one of the two samples cannot be compared
no_reads <- colSums(counts_direct) == 0 | colSums(counts_culture) == 0
if (any(no_reads)) cat("Patients without genus-level reads in one sample type:", patients[no_reads], "\n")
patients <- patients[!no_reads]
counts_direct  <- counts_direct[, patients]
counts_culture <- counts_culture[, patients]

# All genera of the object are kept for the CLR (its mean is taken over all genera)
cat("Patients with both sample types:", length(patients), "| genera:", nrow(counts_direct), "\n")

#### End ####

#### Delta of the CLR-transformed abundances ####

clr <- function(counts) {
  log_rel <- log(sweep(counts, 2, colSums(counts), "/") + pseudocount)
  sweep(log_rel, 2, colMeans(log_rel), "-")
}
delta <- clr(counts_culture) - clr(counts_direct)  # genera x patients

# Prevalence filter: present in >= 10% of the patients, in either sample type
prevalence <- rowMeans(counts_direct > 0 | counts_culture > 0)
delta <- delta[prevalence >= prevalence_threshold, ]
cat("Genera tested (present in >= ", 100 * prevalence_threshold, "% of patients): ", nrow(delta), "\n", sep = "")

# Clinical variables of the patients (from the directly isolated samples)
meta <- direct[match(patients, direct$Study_Nr), meta_vars]
rownames(meta) <- patients
meta[] <- lapply(meta, function(x) as.character(x))
meta$N <- ifelse(is.na(meta$N), NA, paste0("N", meta$N))
meta$M <- ifelse(is.na(meta$M), NA, paste0("M", meta$M))
cat("Missing values:", paste(meta_vars, colSums(is.na(meta)), sep = " = ", collapse = ", "), "\n")

#### End ####

#### 1. Non-parametric tests per genus ####

var_labels <- c(N = "N", M = "M", History.of.smoking.y.n = "Smoking")

test_one <- function(x, group) {
  ok <- !is.na(group)
  x <- x[ok]; group <- group[ok]
  if (length(unique(group)) == 2) {
    t <- suppressWarnings(wilcox.test(x ~ group))
    c(p_value = t$p.value, statistic = unname(t$statistic), test = "wilcoxon")
  } else {
    t <- kruskal.test(x, factor(group))
    c(p_value = t$p.value, statistic = unname(t$statistic), test = "kruskal")
  }
}

tests <- bind_rows(lapply(meta_vars, function(v) {
  res <- t(apply(delta, 1, test_one, group = meta[[v]]))
  tibble(species = rownames(delta),
         p_value = as.numeric(res[, "p_value"]),
         statistic = as.numeric(res[, "statistic"]),
         test = res[, "test"],
         variable = var_labels[[v]]) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    arrange(fdr)
}))

#### End ####

#### 2. MaAsLin2: delta ~ N + M + smoking ####

delta_input <- as.data.frame(t(delta))
delta_meta <- meta
delta_meta$N <- as.numeric(sub("N", "", delta_meta$N))
delta_meta$M <- as.numeric(sub("M", "", delta_meta$M))

maaslin_dir <- file.path(delta.dir, "maaslin2_delta")
invisible(capture.output(
  Maaslin2(input_data = delta_input, input_metadata = delta_meta, output = maaslin_dir,
           fixed_effects = meta_vars, normalization = "NONE", transform = "NONE",
           analysis_method = "LM", standardize = FALSE,
           min_prevalence = 0, min_abundance = -Inf,  # delta can be negative
           plot_heatmap = FALSE, plot_scatter = FALSE)
))

maaslin <- read.delim(file.path(maaslin_dir, "all_results.tsv")) %>%
  mutate(feature = rownames(delta)[match(feature, make.names(rownames(delta)))],  # genus names back
         metadata = unname(var_labels[metadata])) %>%
  rename(`N not 0` = N.not.0)

#### End ####

#### 3. PERMANOVA on the delta profiles ####

complete <- complete.cases(meta)
dist_delta <- vegdist(t(delta[, complete]), method = "euclidean")
meta_complete <- meta[complete, ]
colnames(meta_complete) <- unname(var_labels[colnames(meta_complete)])

set.seed(1)
permanova <- adonis2(dist_delta ~ N + M + Smoking, data = meta_complete, permutations = 999, by = "margin")
permanova <- rownames_to_column(as.data.frame(permanova), " ")
cat("PERMANOVA: ", sum(complete), " of ", length(patients), " patients with N, M and smoking\n", sep = "")
print(permanova)

#### End ####

#### Supplementary Table 4 ####

write.xlsx(list(`Wilcoxon_rank_sum-Kruskal_delta` = tests,
                MaAsLin2_delta = maaslin,
                PERMANOVA_delta = permanova),
           file.path(delta.dir, "SuppTable4_delta_culture_vs_direct.xlsx"))

#### End ####
