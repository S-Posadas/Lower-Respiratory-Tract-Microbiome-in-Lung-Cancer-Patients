# Alpha diversity #

#### Package setup ####

Sys.setenv(language = "EN")

# Tidyverse (load dplyr AFTER plyr)
library(dplyr)
library(tidyr)    
library(ggplot2)

# Other packages
library(phyloseq)
library(ggpubr)
library(ggstatsplot)
library(ggpattern)
library(grid)
library(gridExtra)
library(openxlsx)
library(cowplot)
library(car)      # Type III F-test in the permutation ANCOVA
library(emmeans)  # Pairwise post-hoc contrasts in the permutation ANCOVA

theme_set(theme_bw())

#### End ####

#### Load output from DADA2 pipeline with original phyloseq object ####

res.dir <- "Results_direct"
load(file.path(res.dir, "RData", "4.physeq.decontam.RData"))

#### End ####

#### Create directories for results ####

a.stats <- file.path(res.dir, "alpha_stats")
a.plots <- file.path(res.dir, "alpha_plots")

dir.create(a.stats, recursive = T)
dir.create(a.plots, recursive = T)

#### End ####

#### Plot theme settings ####
text_size = 30
thm <- theme(text = element_text(size = text_size),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "bottom")

thm.x <- theme(text = element_text(size = text_size),
               axis.text.x = element_text(angle = -30, hjust = 0),
               legend.position = "bottom")
#### End ####

#### Select samples ####

physeq_alpha = physeq
physeq_alpha <- subset_samples(physeq_alpha, 
                               Sample_or_Control == "True sample" &
                                 Isolation == "Direct_Isolation"# &
                               # Diagnosis %in% c("Benign", "SCLC", "NSCLC")
)
physeq_alpha = subset_samples(physeq_alpha, !is.na(colSums(otu_table(physeq_alpha))) & colSums(otu_table(physeq_alpha)) != 0) 

#### End ####

#### Convert sample data to data frame ####

rich <- data.frame(sample_data(physeq_alpha))
rich$sample <- rownames(rich)
dim(rich)

# Non-smokers have 0 pack-years by definition; impute where the value is
# missing rather than letting complete.cases() drop them from adjusted models
rich$Pack.years[rich$History.of.smoking.y.n == "No" & is.na(rich$Pack.years)] <- 0

# Convert to long data frame for the facets
rich_long <- rich %>% pivot_longer(
  cols = c(Chao1, Shannon, InvSimpson),
  names_to = "alpha_measure",
  values_to = "alpha_value"
) %>% mutate(alpha_measure = factor(alpha_measure, levels = c("Chao1", "Shannon", "InvSimpson")))

#### End ####

#### Explore correlation of sequencing depth and alpha diversity ####

for(index in c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")){
  
  ggscatterstats(
    data = rich,
    x = !!sym(index), 
    y = LibrarySizeDecontam, 
    xlab = index,
    ylab = "Library Size",
    point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
    marginal = F,
    type ="n",
    title = "Directly isolated samples",
  ) + theme(text = element_text(size = 32),
            plot.title = element_text(size = 28),
            plot.subtitle = element_text(size = 18))
  
  ggsave(file.path(qc.dir, paste0("LibSize_", index, "_direct.svg")), width = 11, height = 10, dpi = 300)
  
}

#### End ####

#### Analysis Functions ####
run_wilcox_test <- function(data, div_measure, group, paired = FALSE, exact = TRUE) {
  results <- list()
  for (x in div_measure) {
    if (paired) {
      # Paired test requires two separate vectors
      group_levels <- unique(data[[group]])
      vec1 <- data[data[[group]] == group_levels[1], x]
      vec2 <- data[data[[group]] == group_levels[2], x]
      
      test <- wilcox.test(vec1, vec2, paired = TRUE, exact = exact)
    } else {
      # Independent test can use formula
      test <- wilcox.test(as.formula(paste(x, "~", group)), 
                          data = data, exact = exact)
    }
    
    z <- abs(qnorm(test$p.value/2))
    r <- z/sqrt(nrow(data))  # Effect size calculation
    results[[x]] <- data.frame(
      Test = test$method,
      Variable1 = x,
      Variable2 = group,
      Statistic = test$statistic,
      p.value = test$p.value,
      Effect.size = r
    )
  }
  results_df <- as.data.frame(do.call(rbind, results)) 
  results_df$p.adj.BH <- p.adjust(results_df$p.value, method = "BH")
  results_df$p.adj.bonferroni <- p.adjust(results_df$p.value, method = "bonferroni")
  
  return(results_df)
}

#' Permutation ANCOVA (Freedman-Lane procedure): tests the group effect on the
#' original (unranked) scale, so covariate relationships are kept as fitted,
#' without assuming normal residuals (the null distribution of the F-statistic
#' is built by permutation). Computes the omnibus F-test and all pairwise
#' post-hoc contrasts from one set of permuted model refits.
run_permutation_ancova_full <- function(data, div_measure, group, covariates = NULL, n_perm = 999, seed = 123, p.adjust.method = "BH") {
  omnibus_list <- list()
  posthoc_list <- list()
  
  for (x in div_measure) {
    keep_cols <- c(x, group, covariates)
    model_data <- data[, keep_cols, drop = FALSE]
    model_data <- model_data[complete.cases(model_data), ]
    model_data[[group]] <- droplevels(as.factor(model_data[[group]]))
    n_obs <- nrow(model_data)
    
    covar_terms <- if (!is.null(covariates) && length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
    full_formula <- if (!is.null(covar_terms)) paste(x, "~", group, "+", covar_terms) else paste(x, "~", group)
    reduced_formula <- if (!is.null(covar_terms)) paste(x, "~", covar_terms) else paste(x, "~ 1")
    
    full_model <- lm(as.formula(full_formula), data = model_data)
    reduced_model <- lm(as.formula(reduced_formula), data = model_data)
    
    aov_table <- car::Anova(full_model, type = 3)
    row_idx <- which(rownames(aov_table) == group)
    resid_idx <- which(rownames(aov_table) == "Residuals")
    f_obs <- aov_table$`F value`[row_idx]
    df1 <- aov_table$Df[row_idx]
    df2 <- aov_table$Df[resid_idx]
    
    emm_obs <- emmeans::emmeans(full_model, specs = as.formula(paste("~", group)))
    pw_obs <- as.data.frame(pairs(emm_obs, adjust = "none"))
    obs_t <- setNames(pw_obs$t.ratio, pw_obs$contrast)
    
    resid_reduced <- residuals(reduced_model)
    fitted_reduced <- fitted(reduced_model)
    
    set.seed(seed)
    f_perm <- numeric(n_perm)
    perm_t <- matrix(NA, nrow = n_perm, ncol = length(obs_t), dimnames = list(NULL, names(obs_t)))
    
    for (i in seq_len(n_perm)) {
      y_perm <- fitted_reduced + sample(resid_reduced)
      perm_data <- model_data
      perm_data[[x]] <- y_perm
      perm_model <- lm(as.formula(full_formula), data = perm_data)
      
      # Both statistics come from the SAME permuted refit
      perm_aov <- car::Anova(perm_model, type = 3)
      f_perm[i] <- perm_aov$`F value`[which(rownames(perm_aov) == group)]
      
      emm_perm <- emmeans::emmeans(perm_model, specs = as.formula(paste("~", group)))
      pw_perm <- as.data.frame(pairs(emm_perm, adjust = "none"))
      perm_t[i, ] <- pw_perm$t.ratio[match(names(obs_t), pw_perm$contrast)]
    }
    
    p_omnibus <- (1 + sum(f_perm >= f_obs)) / (n_perm + 1)
    partial_eta2 <- (f_obs * df1) / (f_obs * df1 + df2)
    
    omnibus_list[[x]] <- data.frame(
      Test = paste0("Permutation ANCOVA (Freedman-Lane, ", n_perm, " permutations)"),
      Variable1 = x,
      Variable2 = group,
      Covariates = if (is.null(covariates)) "None" else paste(covariates, collapse = ", "),
      N = n_obs,
      F.value = f_obs,
      df1 = df1,
      df2 = df2,
      p.value = p_omnibus,
      Partial.eta.sq = partial_eta2
    )
    
    p_perm_pw <- sapply(names(obs_t), function(cn) {
      (1 + sum(abs(perm_t[, cn]) >= abs(obs_t[cn]))) / (n_perm + 1)
    })
    pw_obs$p.perm <- p_perm_pw
    pw_obs$Diversity_measure <- x
    pw_obs$Grouping_variable <- group
    pw_obs$Covariates <- if (is.null(covariates)) "None" else paste(covariates, collapse = ", ")
    pw_obs$Test <- paste0("Permutation ANCOVA post-hoc (pooled model, shared permutations, ", n_perm, ")")
    
    posthoc_list[[x]] <- pw_obs
  }
  
  omnibus_df <- as.data.frame(do.call(rbind, omnibus_list))
  omnibus_df$p.adj.BH <- p.adjust(omnibus_df$p.value, method = "BH")
  omnibus_df$p.adj.bonferroni <- p.adjust(omnibus_df$p.value, method = "bonferroni")
  
  posthoc_df <- do.call(rbind, posthoc_list)
  posthoc_df$p.perm.adj.BH <- NA
  for (measure in div_measure) {
    idx <- posthoc_df$Diversity_measure == measure
    posthoc_df$p.perm.adj.BH[idx] <- p.adjust(posthoc_df$p.perm[idx], method = p.adjust.method)
  }
  
  return(list(omnibus = omnibus_df, posthoc = posthoc_df))
}

save_results <- function(results, file_name) {
  write.xlsx(list(results, rich_sub), file.path(a.stats, file_name))
}

#' Save the permutation ANCOVA (omnibus + post-hoc) and the sample data of one
#' question as separate sheets in one xlsx file.
save_permutation <- function(permutation, file_name) {
  sheets <- list(Permutation_omnibus = permutation$omnibus,
                 Permutation_posthoc = permutation$posthoc,
                 Data = rich_sub)
  write.xlsx(sheets, file.path(a.stats, file_name))
}

alpha_plot <- function(data, x_var, col_palette = NULL, y_var = "alpha_value", fill_var = x_var, color_var = x_var, facet_var = "alpha_measure", scales = "free_y", legend_lab = x_var) {
  if(is.null(col_palette)){
    col_palette <- scales::hue_pal()(nlevels(rich[,x_var]))
  }
  
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], 
                   fill = .data[[fill_var]], color = .data[[color_var]])) +
    
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
    geom_boxplot(alpha = 0.5) +
    facet_wrap(as.formula(paste("~", facet_var)), scales = scales) +
    scale_fill_manual(values = col_palette) + scale_color_manual(values = col_palette) +
    labs(y = NULL, x = NULL, color = legend_lab, fill = legend_lab) + thm 
  
}

#### Questions to answer ####

#' Questions
#' Q1: Is there a difference between Diseased vs the parallel Normal Lung?
#' Q2: Is there a difference between Diseased vs the parallel Normal Lung without samples with synchronous tumor?
#' Q3.1: What is the difference between NSCLC vs SCLC vs Bening Tumor?
#' Q3.2: What is the difference between NSCLC vs Bening Tumor in nonsmokers?
#' Q4.1: What is the difference between NSCLC main histologies (Adenocarcinoma vs Squamous)?
#' Q4.2: What is the difference between NSCLC main histologies (Adenocarcinoma vs Squamous) in nonsmokers?
#' Q4:.3 What is the difference between NSCLC main histologies (Adenocarcinoma vs Squamous) in smokers?
#' Q5: Does history of smoking have an impact on the lung microbiome?
#' Q6.1: Does the T stage have an impact on the lung microbiome?
#' Q6.2: Does the N stage have an impact on the lung microbiome?
#' Q6.3: Does the M stage have an impact on the lung microbiome?

#' Covariate policy:
#' - Q1, Q2 (paired Diseased vs Normal lung, same patient): no covariates.
#'   Age/sex/smoking status don't vary within a patient, so pairing already
#'   controls for them.
#' - Q3.2, Q3.3, Q4.2, Q4.3, all Q5: adjust for Age + Sex only (smoking status
#'   is either the stratifying variable or the variable of interest here, so
#'   it is not also entered as a covariate).
#' - Everywhere else (Q3.1, Q4.1, Q6.x): adjust for History.of.smoking.y.n + Age + Sex.
#'
#' Every covariate-adjusted question below is tested with a permutation
#' ANCOVA (Freedman-Lane, 999 permutations; omnibus F-test and pairwise
#' post-hoc contrasts). Samples with a missing covariate are left out.
#' The parametric and rank ANCOVA and their diagnostics are kept in a
#' separate working copy (not part of this repository).

# Maaslin2 analysis parameters
params <- list(
  Q1.1 = list(name = "Q1.1.all_directly_isolated_lung", fixed = "Lung", rand = NULL, ref = "Lung,Normal", col = lung_col, covariates = NULL),
  Q1.2 = list(name = "Q1.2.main-D_directly_isolated_lung", fixed = "Lung", rand = NULL, ref = "Lung,Normal", col = lung_col, covariates = NULL),
  Q1.3 = list(name = "Q1.3.main-D-individually_directly_isolated_lung", covariates = NULL),
  Q2.1 = list(name = "Q2.1.all_directly_isolated_non_synchronous_lung", fixed = "Lung", rand = NULL, ref = "Lung,Normal", col = lung_col, covariates = NULL),
  Q2.2 = list(name = "Q2.2.main-D_directly_isolated_non_synchronous_lung", fixed = "Lung", rand = NULL, ref = "Lung,Normal", col = lung_col, covariates = NULL),
  Q2.3 = list(name = "Q2.3.main-D-individually_directly_isolated_non_synchronous_lung", covariates = NULL),
  Q3.1 = list(name = "Q3.1.diagnosis", fixed = "Diagnosis", rand = NULL, ref = "Diagnosis,Benign", col = diagnosis_col, covariates = c("History.of.smoking.y.n", "Age", "Sex")),
  Q3.2 = list(name = "Q3.2.diagnosis_nonsmokers", fixed = "Diagnosis", rand = NULL, ref = "Diagnosis,Benign", col = diagnosis_col, covariates = c("Age", "Sex")),
  Q3.3 = list(name = "Q3.3.diagnosis_smokers", fixed = "Diagnosis", rand = NULL, ref = "Diagnosis,Benign", col = diagnosis_col, covariates = c("Age", "Sex")),
  Q4.1 = list(name = "Q4.1.histology", fixed = "Histology.NSCLC", rand = NULL, ref = NULL, col = histology_col, covariates = c("History.of.smoking.y.n", "Age", "Sex")),
  Q4.2 = list(name = "Q4.2.histology_nonsmokers", fixed = "Histology.NSCLC", rand = NULL, ref = NULL, col = histology_col, covariates = c("Age", "Sex")),
  Q4.3 = list(name = "Q4.3.histology_smokers", fixed = "Histology.NSCLC", rand = NULL, ref = NULL, col = histology_col, covariates = c("Age", "Sex")),
  Q5.0 = list(name = "Q5.0.all_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, covariates = c("Age", "Sex")),
  Q5.1 = list(name = "Q5.1.main-D_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, covariates = c("Age", "Sex")),
  Q5.2 = list(name = "Q5.2.NSCLC_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, covariates = c("Age", "Sex")),
  Q5.3 = list(name = "Q5.3.Benign_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, covariates = c("Age", "Sex")),
  Q5.4 = list(name = "Q5.4.Adenocarcinoma_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, lab_plot = "History of smoking", covariates = c("Age", "Sex")),
  Q6.1 = list(name = "Q6.1.NSCLC_T", fixed = "T", rand = NULL, ref = NULL, covariates = c("History.of.smoking.y.n", "Age", "Sex")),
  Q6.2 = list(name = "Q6.2.NSCLC_N", fixed = "N", rand = NULL, ref = NULL, covariates = c("History.of.smoking.y.n", "Age", "Sex")),
  Q6.2.1 = list(name = "Q6.2.NSCLC_N_twolevels", fixed = "N", rand = NULL, ref = NULL, covariates = c("History.of.smoking.y.n", "Age", "Sex")),
  Q6.3 = list(name = "Q6.3.NSCLC_M", fixed = "M", rand = NULL, ref = NULL, covariates = c("History.of.smoking.y.n", "Age", "Sex")),
  Q6.3.1 = list(name = "Q6.3.1.NSCLC_M_smokers", fixed = "M", rand = NULL, ref = NULL, covariates = c("Age", "Sex"))
)

# Function to subset phyloseq object
subset_data <- function(data, q){
  
  # Custom subsetting for specific questions
  if (grepl("Q1.2", q)) {
    data_sub <- dplyr::filter(data, Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
  } else if (grepl("Q2", q)) {
    data_sub <- dplyr::filter(data, Synchronous.tumor == "No")
    if (grepl("Q2.2", q)) {
      data_sub <- dplyr::filter(data_sub, Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    }
  } else if (grepl("Q3", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" &
                                Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    if (grepl("Q3.2", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    }else if (grepl("Q3.3", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    }
  } else if (grepl("Q4", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" & 
                                Diagnosis == "NSCLC" &
                                Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma")
    )
    if (grepl("Q4.2", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    }else if (grepl("Q4.3", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    }
  } else if (grepl("Q5.0", q)) {
    # Whole cohort: all diagnoses (incl. carcinoids and metastases)
    data_sub <- dplyr::filter(data, Lung == "Diseased" & !is.na(History.of.smoking.y.n))
  } else if (grepl("Q5", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" & 
                                Diagnosis %in% c("Benign", "NSCLC", "SCLC") &
                                !is.na(History.of.smoking.y.n)
    )
    if (grepl("Q5.2", q)) {
      data_sub <- dplyr::filter(data_sub, Diagnosis == "NSCLC")
    }else if (grepl("Q5.3", q)) {
      data_sub <- dplyr::filter(data_sub, Diagnosis == "Benign")
    }else if (grepl("Q5.4", q)) {
      data_sub <- dplyr::filter(data_sub, Diagnosis == "NSCLC" & Histology.NSCLC == "Adenocarcinoma")
    }
  } else if (grepl("Q6", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" & 
                                Diagnosis == "NSCLC"
    )
    if (grepl("Q6.1", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(T))
    }else if (grepl("Q6.2", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(N))
    }else if (grepl("Q6.3", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(M))
      if(grepl("Q6.3.1", q)){
        data_sub <- dplyr::filter(data_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
      }
    }
  } else {
    data_sub <- data
  }
  
  return(data_sub)
}

#### Q1: Is there a difference between Diseased vs the parallel Normal Lung? ####
#### (paired within patient -> no covariates needed, see policy note above) ####

### Q1.1.all_directly_isolated_lung (paired test):
q = "Q1.1"

rich_sub <- subset_data(rich, q)
rich_sub <- rich_sub[rich_sub$Study_Nr %in% rich_sub$Study_Nr[duplicated(rich_sub$Study_Nr)], ]
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

# First ensure paired data is properly ordered
paired_data <- rich_sub %>%
  filter(Lung %in% c("Diseased", "Normal")) %>%
  arrange(Study_Nr, Lung)  # Critical for paired tests

wt_results <- run_wilcox_test(paired_data, c("Shannon", "InvSimpson"), "Lung", paired = TRUE)

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Lung", col_palette = params[[q]]$col)

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q1.2.main-D_directly_isolated_lung
q = "Q1.2"

rich_sub <- subset_data(rich, q)
rich_sub <- rich_sub[rich_sub$Study_Nr %in% rich_sub$Study_Nr[duplicated(rich_sub$Study_Nr)], ]
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

# First ensure paired data is properly ordered
paired_data <- rich_sub %>%
  filter(Lung %in% c("Diseased", "Normal")) %>%
  arrange(Study_Nr, Lung)  # Critical for paired tests

wt_results <- run_wilcox_test(paired_data, c("Shannon", "InvSimpson"), "Lung", paired = TRUE)
wt_results
save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Lung", col_palette = params[[q]]$col)

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q1.3.main-D-individually_directly_isolated_lung
# In each of the three main diagnoses
q = "Q1.3"

wt <- list()
for (d in c("NSCLC", "SCLC", "Benign")) {  
  rich_sub <- subset_data(rich, q)  %>% filter(Diagnosis == d)
  rich_sub <- rich_sub[rich_sub$Study_Nr %in% rich_sub$Study_Nr[duplicated(rich_sub$Study_Nr)], ]
  paired_data <- rich_sub %>%
    filter(Lung %in% c("Diseased", "Normal")) %>%
    arrange(Study_Nr, Lung)  # Critical for paired tests
  
  wt_results <- run_wilcox_test(paired_data, c("Shannon", "InvSimpson"), "Lung", paired = TRUE)
  
  wt[[d]] <- wt_results
}
wt <- as.data.frame(do.call(rbind, wt))

save_results(wt, paste0(params[[q]]$name, ".xlsx"))

Plot <- ggplot(rich_sub_long, aes(x = Lung, y = alpha_value, color = Diagnosis, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = diagnosis_col) + scale_color_manual(values = diagnosis_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "Lung", color = "Diagnosis", fill = "Diagnosis") + thm.x

Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Plot all diagnoses + individual main diagoses 
rich_sub <- subset_data(rich, "Q1.1")  # Q1.1.all_directly_isolated_lung
rich_sub <- rich_sub[rich_sub$Study_Nr %in% rich_sub$Study_Nr[duplicated(rich_sub$Study_Nr)], ]
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]
rich_all <- rich_sub_long
rich_all$Diagnosis <- "All diagnoses"
rich_all <- rbind(rich_all, rich_sub_long[rich_sub_long$Diagnosis %in% c("NSCLC", "SCLC", "Benign"),])

Plot_extra <- ggplot(rich_all, aes(x = Diagnosis, y = alpha_value, color = Lung, fill = Lung)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = lung_col) + scale_color_manual(values = lung_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "Lung", fill = "Lung") + thm.x

Plot_extra

ggsave(plot = Plot_extra, file.path(a.plots, "Q1.all.svg"), width = 10, height = 8, dpi = 300)

#### Q2: Is there a difference between Diseased vs the parallel Normal Lung without samples with synchronous tumor? ####
#### (paired within patient -> no covariates needed, see policy note above) ####

### Q2.1.all_directly_isolated_non_synchronous_lung (paired test):
q = "Q2.1"

rich_sub <- subset_data(rich, q)
rich_sub <- rich_sub[rich_sub$Study_Nr %in% rich_sub$Study_Nr[duplicated(rich_sub$Study_Nr)], ]
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

# First ensure paired data is properly ordered
paired_data <- rich_sub %>%
  filter(Lung %in% c("Diseased", "Normal")) %>%
  arrange(Study_Nr, Lung)  # Critical for paired tests

wt_results <- run_wilcox_test(paired_data, c("Shannon", "InvSimpson"), "Lung", paired = TRUE)
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Lung", col_palette = params[[q]]$col)

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q2.2.main-D_directly_isolated_non_synchronous_lung
q = "Q2.2"

rich_sub <- subset_data(rich, q)
rich_sub <- rich_sub[rich_sub$Study_Nr %in% rich_sub$Study_Nr[duplicated(rich_sub$Study_Nr)], ]
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

# First ensure paired data is properly ordered
paired_data <- rich_sub %>%
  filter(Lung %in% c("Diseased", "Normal")) %>%
  arrange(Study_Nr, Lung)  # Critical for paired tests

wt_results <- run_wilcox_test(paired_data, c("Shannon", "InvSimpson"), "Lung", paired = TRUE)
wt_results
save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Lung", col_palette = params[[q]]$col)

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q2.3.main-D-individually_directly_isolated_non_synchronous_lung
# In each of the three main diagnoses
q = "Q2.3"

wt <- list()
for (d in c("NSCLC", "SCLC")) {  
  rich_sub <- subset_data(rich, "Q2.1")  %>% filter(Diagnosis == d)
  rich_sub <- rich_sub[rich_sub$Study_Nr %in% rich_sub$Study_Nr[duplicated(rich_sub$Study_Nr)], ]
  paired_data <- rich_sub %>%
    filter(Lung %in% c("Diseased", "Normal")) %>%
    arrange(Study_Nr, Lung)  # Critical for paired tests
  
  wt_results <- run_wilcox_test(paired_data, c("Shannon", "InvSimpson"), "Lung", paired = TRUE)
  
  wt[[d]] <- wt_results
}
wt <- as.data.frame(do.call(rbind, wt))
wt

save_results(wt, paste0(params[[q]]$name, ".xlsx"))

Plot <- ggplot(rich_sub_long, aes(x = Lung, y = alpha_value, color = Diagnosis, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = diagnosis_col) + scale_color_manual(values = diagnosis_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "Lung", color = "Diagnosis", fill = "Diagnosis") + thm.x

Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Plot all diagnoses + individual main diagoses 
rich_sub <- subset_data(rich, "Q2.1")  # Q2.1.all_directly_isolated_non_synchronous_lung
rich_sub <- rich_sub[rich_sub$Study_Nr %in% rich_sub$Study_Nr[duplicated(rich_sub$Study_Nr)], ]
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]
rich_all <- rich_sub_long
rich_all$Diagnosis <- "All diagnoses"
rich_all <- rbind(rich_all, rich_sub_long[rich_sub_long$Diagnosis %in% c("NSCLC", "SCLC", "Benign"),])

Plot_extra <- ggplot(rich_all, aes(x = Diagnosis, y = alpha_value, color = Lung, fill = Lung)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = lung_col) + scale_color_manual(values = lung_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = params[[q]]$fixed, fill = params[[q]]$fixed) + thm.x

Plot_extra

ggsave(plot = Plot_extra, file.path(a.plots, "Q2.all.svg"), width = 10, height = 8, dpi = 300)

#### Conclusion ####
#' There is no significant difference between diseased and normal lung, thus from now on
#' we will only consider the diseased lung

#### Q3: What is the difference between NSCLC vs SCLC vs Bening Tumor? ####
#### Q3.1 adjusted for History.of.smoking.y.n + Age + Sex; Q3.2/Q3.3 adjusted for Age + Sex only ####
### Permutation ANCOVA (Freedman-Lane) ###

### Q3.1.diagnosis
q = "Q3.1"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Chao1", "Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col)
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Plot divided by history of smoking
rich_sub_long <- dplyr::filter(rich_sub_long, !is.na(History.of.smoking.y.n))

ggplot(rich_sub_long, aes(x = Diagnosis, y = alpha_value, 
                          fill = History.of.smoking.y.n, color = History.of.smoking.y.n)) +
  geom_boxplot(alpha = 0.5) +
  # geom_point() +
  facet_wrap(as.formula(paste("~", "alpha_measure")), scales = "free_y") +
  scale_fill_manual(values = smoker_col) + scale_color_manual(values = smoker_col) +
  labs(y = NULL, x = NULL, color = "History of smoking", fill = "History of smoking") + thm.x

ggsave(file.path(a.plots, paste0(params[[q]]$name, "_smoking_noNA.svg")), width = 10, height = 8, dpi = 300)


### Q3.2.diagnosis_nonsmokers
q = "Q3.2"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Chao1", "Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col)
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q3.3.diagnosis_smokers
q = "Q3.3"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Chao1", "Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col)
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

#### Q4: What are the difference between NSCLC main histologies (Adeno vs Squamous)? ####
#### Q4.1 adjusted for History.of.smoking.y.n + Age + Sex; Q4.2/Q4.3 adjusted for Age + Sex only ####
### Permutation ANCOVA (Freedman-Lane) ###

### Q4.1.histology
q = "Q4.1"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Shannon", "InvSimpson"), "Histology.NSCLC", covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Histology.NSCLC", col_palette = params[[q]]$col, legend_lab = "Histology")
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Plot divided by history of smoking
rich_sub_long <- dplyr::filter(rich_sub_long, !is.na(History.of.smoking.y.n))

ggplot(rich_sub_long, aes(x = Histology.NSCLC, y = alpha_value, 
                          fill = History.of.smoking.y.n, color = History.of.smoking.y.n)) +
  geom_boxplot(alpha = 0.5) +
  # geom_point() +
  facet_wrap(as.formula(paste("~", "alpha_measure")), scales = "free_y") +
  scale_fill_manual(values = smoker_col) + scale_color_manual(values = smoker_col) +
  labs(y = NULL, x = NULL, color = "History of smoking", fill = "History of smoking") + thm.x

ggsave(file.path(a.plots, paste0(params[[q]]$name, "_smoking_noNA.svg")), width = 10, height = 8, dpi = 300)


### Q4.2.histology_nonsmokers
q = "Q4.2"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Shannon", "InvSimpson"), "Histology.NSCLC", covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Histology.NSCLC", col_palette = params[[q]]$col, legend_lab = "Histology")
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q4.3.histology_smokers
q = "Q4.3"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Shannon", "InvSimpson"), "Histology.NSCLC", covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Histology.NSCLC", col_palette = params[[q]]$col, legend_lab = "Histology")
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

#### Q5: Does history of smoking have an impact on the lung microbiome in NSCLC and Adenocarcinoma? ####
#### All Q5 sub-questions adjusted for Age + Sex only (smoking is the variable of interest) ####
#### Permutation ANCOVA (Freedman-Lane) ####

### Q5.0.all_smoking (Supplementary Figure 12A: whole study cohort, all diagnoses)
q = "Q5.0"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q5.1.main-D_smoking
q = "Q5.1"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

dim(rich_sub_long)

rich_sub_long$Diagnosis <- droplevels(factor(rich_sub_long$Diagnosis))
Plot <- alpha_plot(rich_sub_long, x_var = "Diagnosis", fill_var = params[[q]]$fixed, color_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot$layers[[1]]$position <- position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)  # points
Plot$layers[[2]]$position <- position_dodge2(preserve = "single")                          # boxes
Plot + thm.x
ggsave(plot = Plot + thm.x, file.path(a.plots, paste0(params[[q]]$name, "_diagnosis.svg")), width = 10, height = 8, dpi = 300)

### Q5.2.NSCLC_smoking
q = "Q5.2"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q5.3.Benign_smoking
q = "Q5.3"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q5.4.Adenocarcinoma_smoking
q = "Q5.4"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Correlation between packyears and alpha diversity (raw, unadjusted association;
# not part of the covariate-controlled testing above)
for (q in c("Q5.1", "Q5.2", "Q5.3", "Q5.4")) {
  rich_sub <- subset_data(rich, q)
  pack.years <- list()
  for(index in c("Shannon", "InvSimpson")){
    
    pack.years[[index]] <- ggscatterstats(
      data = rich_sub,
      y = !!sym(index), 
      x = Pack.years, 
      ylab = index,
      xlab = "Pack years",
      point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
      marginal = F,
      type ="n",
      title = NULL
    ) + theme(text = element_text(size = 32),
              plot.subtitle = element_text(size = 16))
  }
  
  plot_grid(top = textGrob("Relationship between alpha diversity and pack years",gp=gpar(fontsize=30, fontface="bold")),
            nrow= 2, rel_heights = c(1,8),
            plot_grid(pack.years[[1]], pack.years[[2]]))
  
  ggsave(file.path(a.plots, paste0(params[[q]]$name, "_pack.years.svg")), width = 18, height = 10, dpi = 300)
}



#### Q6: Does the N and M stage of NSCLC have an impact on the lung microbiome ####
#### All Q6 sub-questions adjusted for History.of.smoking.y.n + Age + Sex (Q6.3.1: Age + Sex only) ####
### Permutation ANCOVA (Freedman-Lane) ###

# Q6.1.NSCLC_T
q = "Q6.1"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Chao1", "Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = scales::hue_pal()(5)[c(4,2,5,1)])
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Q6.2.NSCLC_N
q = "Q6.2"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Chao1", "Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = scales::hue_pal()(5)[c(3:4,2,5,1)])
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Q6.2.NSCLC_N_twolevels
q = "Q6.2.1"

rich_sub <- subset_data(rich, q)
rich_sub$Ntwo <- as.numeric(as.character(rich_sub$N))
rich_sub <- rich_sub %>%
  mutate(Ntwo = case_when(
    Ntwo == 0 ~ "0",
    Ntwo > 0 ~ "1-3"))
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]
rich_sub_long$Ntwo <- as.numeric(as.character(rich_sub_long$N))
rich_sub_long <- rich_sub_long %>%
  mutate(Ntwo = case_when(
    Ntwo == 0 ~ "0",
    Ntwo > 0 ~ "1-3"))

perm_full <- run_permutation_ancova_full(rich_sub, c("Chao1", "Shannon", "InvSimpson"), "Ntwo", covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Ntwo", col_palette = scales::hue_pal()(5)[c(3:4,2,5,1)], legend_lab = "N")
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Plot <- alpha_plot(rich_sub_long, x_var = "Ntwo", col_palette = colors_N[c(2,4)], legend_lab = "N")
# Plot
# 
# ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, "gradient.svg")), width = 10, height = 8, dpi = 300)

# Q6.3.NSCLC_M
q = "Q6.3"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Chao1", "Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = c("lightblue", "firebrick"))
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, "2.svg")), width = 10, height = 8, dpi = 300)

# Q6.3.NSCLC_M_smokers
q = "Q6.3.1"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

perm_full <- run_permutation_ancova_full(rich_sub, c("Chao1", "Shannon", "InvSimpson"), params[[q]]$fixed, covariates = params[[q]]$covariates)
perm_full

save_permutation(perm_full, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = c("lightblue", "firebrick"))
Plot

ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, "2.svg")), width = 10, height = 8, dpi = 300)

#### End ####