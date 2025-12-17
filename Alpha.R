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

theme_set(theme_bw())

#### End ####

#### Load output from DADA2 pipeline with original phyloseq object ####

res.dir <- "Results_trimmed_cutadapt3_20241216/2025"
load(file.path(res.dir, "RData", "4.physeq.decontam.RData"))

#### End ####

#### Create directories for results ####

a.stats <- file.path(res.dir, "2.Alpha_stats_excluded")
a.plots <- file.path(res.dir, "3.Alpha_plots_excluded")

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

# Convert to long data frame for the facets
rich_long <- rich %>% pivot_longer(
  cols = c(Shannon, InvSimpson),
  names_to = "alpha_measure",
  values_to = "alpha_value"
) %>% mutate(alpha_measure = factor(alpha_measure, levels = c("Shannon", "InvSimpson")))

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
  
  ggsave(file.path(qc.dir, paste0("LibSize_", index, "_direct.tiff")), width = 11, height = 10, dpi = 300)
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

run_kruskal_test <- function(data, div_measure, group) {
  results <- list()
  for (x in div_measure) {
    test <- kruskal.test(as.formula(paste(x, "~", group)), data = data)
    eta_squared <- (test$statistic - 3 + 1)/(nrow(data) - 3)
    f <- sqrt(eta_squared/(1-eta_squared))
    results[[x]] <- data.frame(
      Test = test$method,
      Variable1 = x,
      Variable2 = group,
      Chi2 = test$statistic,
      df = test$parameter,
      p.value = test$p.value,
      Effect.size = f
    )
  }
  results_df <- as.data.frame(do.call(rbind, results)) 
  results_df$p.adj.BH <- p.adjust(results_df$p.value, method = "BH")
  results_df$p.adj.bonferroni <- p.adjust(results_df$p.value, method = "bonferroni")
  
  return(results_df)

}

run_pairwise_mannwhitney <- function(data, div_measure, group, p.adjust.method = "BH") {
  results_list <- list()
  group_vector <- data[[group]]
  
  for (x in div_measure) {
    # Extract diversity measure
    div_vector <- data[[x]]
    valid_indices <- !is.na(div_vector) & !is.na(group_vector)
    div_vector <- div_vector[valid_indices]
    group_vector_clean <- factor(group_vector[valid_indices])
    groups <- levels(group_vector_clean)
    
    # Generate UNIQUE pairwise combinations (without duplicates)
    group_pairs <- combn(groups, 2, simplify = FALSE)
    
    # Run pairwise Wilcoxon test
    pairwise_result <- pairwise.wilcox.test(
      x = div_vector,
      g = group_vector_clean,
      p.adjust.method = "none",
      exact = FALSE
    )
    p_matrix <- pairwise_result$p.value
    
    for (pair in group_pairs) {
      g1 <- pair[1]
      g2 <- pair[2]
      
      # Get p-value (handles matrix orientation)
      p_val <- if(g1 %in% rownames(p_matrix) && g2 %in% colnames(p_matrix)) {
        p_matrix[g1, g2]
      } else if(g2 %in% rownames(p_matrix) && g1 %in% colnames(p_matrix)) {
        p_matrix[g2, g1]
      } else NA
      
      if (!is.na(p_val)) {
        vals1 <- div_vector[group_vector_clean == g1]
        vals2 <- div_vector[group_vector_clean == g2]
        
        # Calculate effect size
        n1 <- length(vals1)
        n2 <- length(vals2)
        ranks <- rank(c(vals1, vals2))
        R1 <- sum(ranks[1:n1])
        U1 <- R1 - n1 * (n1 + 1) / 2
        r_effect <- 2 * (U1 / (n1 * n2)) - 1
        
        results_list[[length(results_list) + 1]] <- data.frame(
          Diversity_measure = x,
          Grouping_variable = group,
          Group1 = g1,
          Group2 = g2,
          p.unadj = p_val,
          Effect.size = r_effect
        )
      }
    }
  }
  
  results_df <- do.call(rbind, results_list)
  
  # Adjust p-values within each diversity measure
  results_df$p.adj <- NA
  for (measure in div_measure) {
    idx <- results_df$Diversity_measure == measure
    results_df$p.adj[idx] <- p.adjust(results_df$p.unadj[idx], method = p.adjust.method)
  }
  
  return(results_df)
}

save_results <- function(results, file_name) {
  write.xlsx(list(results, rich_sub), file.path(a.stats, file_name))
}

alpha_plot <- function(data, x_var, col_palette = NULL, y_var = "alpha_value", fill_var = x_var, color_var = x_var, facet_var = "alpha_measure", scales = "free_y", legend_lab = x_var) {
  if(is.null(col_palette)){
    col_palette <- scales::hue_pal()(nlevels(rich[,x_var]))
  }
  
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], 
                   fill = .data[[fill_var]], color = .data[[color_var]])) +
    geom_boxplot(alpha = 0.5) +
   # geom_point() +
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
#' Q7: Does the microbiome differ between patients who develop post operative pneumonia?

# Maaslin2 analysis parameters
params <- list(
   Q1.1a = list(name = "Q1.1a.all_directly_isolated_lung", fixed = "Lung", rand = NULL, ref = "Lung,Normal", col = lung_col),
   Q1.2a = list(name = "Q1.2a.main-D_directly_isolated_lung", fixed = "Lung", rand = NULL, ref = "Lung,Normal", col = lung_col),
   Q1.3 = list(name = "Q1.3a.main-D-individually_directly_isolated_lung"),
   Q2.1a = list(name = "Q2.1a.all_directly_isolated_non_synchronous_lung", fixed = "Lung", rand = NULL, ref = "Lung,Normal", col = lung_col),
   Q2.2a = list(name = "Q2.2a.main-D_directly_isolated_non_synchronous_lung", fixed = "Lung", rand = NULL, ref = "Lung,Normal", col = lung_col),
   Q2.3 = list(name = "Q2.3a.main-D-individually_directly_isolated_non_synchronous_lung"),
   Q3.1a = list(name = "Q3.1a.diagnosis", fixed = "Diagnosis", rand = NULL, ref = "Diagnosis,Benign", col = diagnosis_col),
   Q3.2a = list(name = "Q3.2a.diagnosis_nonsmokers", fixed = "Diagnosis", rand = NULL, ref = "Diagnosis,Benign", col = diagnosis_col),
   Q3.3a = list(name = "Q3.3a.diagnosis_smokers", fixed = "Diagnosis", rand = NULL, ref = "Diagnosis,Benign", col = diagnosis_col),
   Q4.1a = list(name = "Q4.1a.histology", fixed = "Histology.NSCLC", rand = NULL, ref = NULL, col = histology_col),
   Q4.2a = list(name = "Q4.2a.histology_nonsmokers", fixed = "Histology.NSCLC", rand = NULL, ref = NULL, col = histology_col),
   Q4.3a = list(name = "Q4.3a.histology_smokers", fixed = "Histology.NSCLC", rand = NULL, ref = NULL, col = histology_col),
   Q5.1a = list(name = "Q5.1a.main-D_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col),
   Q5.2a = list(name = "Q5.2a.NSCLC_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col),
   Q5.3a = list(name = "Q5.3a.Benign_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col),
   Q5.4a = list(name = "Q5.4a.Adenocarcinoma_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, lab_plot = "History of smoking"),
   Q6.1a = list(name = "Q6.1a.NSCLC_T", fixed = "T", rand = NULL, ref = NULL),
   Q6.2a = list(name = "Q6.2a.NSCLC_N", fixed = "N", rand = NULL, ref = NULL),
   Q6.2.1a = list(name = "Q6.2a.NSCLC_N_twolevels", fixed = "N", rand = NULL, ref = NULL),
   Q6.3a = list(name = "Q6.3a.NSCLC_M", fixed = "M", rand = NULL, ref = NULL),
   Q7a = list(name = "Q7a.NSCLC_PostopPneumonia", fixed = "Postop..pneumonia", rand = NULL, ref = NULL),
   Q7b = list(name = "Q7b.NSCLC_PostopPneumonia_batch", fixed = "Postop..pneumonia", rand = "Batch", ref = NULL)
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
    }
  } else if (grepl("Q7", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" & 
                                Diagnosis == "NSCLC" &
                                !is.na(Postop..pneumonia)
    )
  } else {
    data_sub <- data
  }
  
  return(data_sub)
}

#### Q1: Is there a difference between Diseased vs the parallel Normal Lung? ####

### Q1.1a.all_directly_isolated_lung (paired test):
q = "Q1.1a"

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

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q1.2a.main-D_directly_isolated_lung
q = "Q1.2a"

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

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q1.3a.main-D-individually_directly_isolated_lung
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

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Plot all diagnoses + individual main diagoses 
rich_sub <- subset_data(rich, "Q1.1a")  # Q1.1a.all_directly_isolated_lung
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

#ggsave(plot = Plot_extra, file.path(a.plots, "Q1.all.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot_extra, file.path(a.plots, "Q1.all.svg"), width = 10, height = 8, dpi = 300)

#### Q2: Is there a difference between Diseased vs the parallel Normal Lung without samples with synchronous tumor? ####

### Q2.1a.all_directly_isolated_non_synchronous_lung (paired test):
q = "Q2.1a"

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

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q2.2a.main-D_directly_isolated_non_synchronous_lung
q = "Q2.2a"

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

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q2.3a.main-D-individually_directly_isolated_non_synchronous_lung
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

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Plot all diagnoses + individual main diagoses 
rich_sub <- subset_data(rich, "Q2.1a")  # Q2.1a.all_directly_isolated_non_synchronous_lung
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

#ggsave(plot = Plot_extra, file.path(a.plots, "Q2.all.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot_extra, file.path(a.plots, "Q2.all.svg"), width = 10, height = 8, dpi = 300)

#### Conclusion ####
#' There is no significant difference between diseased and normal lung, thus from now on
#' we will only consider the diseased lung

#### Q3: What is the difference between NSCLC vs SCLC vs Bening Tumor? ####

### Q3.1a.diagnosis
q = "Q3.1a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

kw_results <- run_kruskal_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
kw_results

save_results(kw_results, paste0(params[[q]]$name, ".xlsx"))

kwph_results <- run_pairwise_mannwhitney(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
kwph_results

save_results(kwph_results, paste0(params[[q]]$name, "posthoc.xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col)
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
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


### Q3.2a.diagnosis_nonsmokers
q = "Q3.2a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col)
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q3.3a.diagnosis_smokers
q = "Q3.3a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

kw_results <- run_kruskal_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
kw_results

save_results(kw_results, paste0(params[[q]]$name, ".xlsx"))

kwph_results <- run_pairwise_mannwhitney(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
kwph_results

save_results(kwph_results, paste0(params[[q]]$name, "posthoc.xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col)
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

#### Q4: What are the difference between NSCLC main histologies (Adeno vs Squamous)? ####

### Q4.1a.histology
q = "Q4.1a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), "Histology.NSCLC")
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Histology.NSCLC", col_palette = params[[q]]$col, legend_lab = "Histology")
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
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
  

### Q4.2a.histology_nonsmokers
q = "Q4.2a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), "Histology.NSCLC")
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Histology.NSCLC", col_palette = params[[q]]$col, legend_lab = "Histology")
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q4.3a.histology_smokers
q = "Q4.3a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), "Histology.NSCLC")
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Histology.NSCLC", col_palette = params[[q]]$col, legend_lab = "Histology")
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

#### Q5: Does history of smoking have an impact on the lung microbiome in NSCLC and Adenocarcinoma? ####

### Q5.1a.main-D_smoking
q = "Q5.1a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

dim(rich_sub_long)

supplement <- rich_sub_long[1:2,]
supplement$Diagnosis <- "SCLC"
supplement$History.of.smoking.y.n <- "No"
rich_sub_long_sup <- rbind(rich_sub_long, supplement)
Plot <- alpha_plot(rich_sub_long_sup, x_var = "Diagnosis", fill_var = params[[q]]$fixed, color_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot + thm.x
ggsave(plot = Plot + thm.x, file.path(a.plots, paste0(params[[q]]$name, "_diagnosis.svg")), width = 10, height = 8, dpi = 300)

### Q5.2a.NSCLC_smoking
q = "Q5.2a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q5.3a.Benign_smoking
q = "Q5.3a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

### Q5.4a.Adenocarcinoma_smoking
q = "Q5.4a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = params[[q]]$col, legend_lab = "History of smoking")
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Correlation between packyears and alpha diversity
for (q in c("Q5.1a", "Q5.2a", "Q5.3a", "Q5.4a")) {
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
  
  #ggsave(file.path(a.plots, paste0(params[[q]]$name, "_pack.years.tiff")), width = 18, height = 10, dpi = 300)
  ggsave(file.path(a.plots, paste0(params[[q]]$name, "_pack.years.svg")), width = 18, height = 10, dpi = 300)
}



#### Q6: Does the N and M stage of NSCLC have an impact on the lung microbiome ####

# Q6.1a.NSCLC_T
q = "Q6.1a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

kw_results <- run_kruskal_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
kw_results

save_results(kw_results, paste0(params[[q]]$name, ".xlsx"))

kwph_results <- run_pairwise_mannwhitney(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
kwph_results

save_results(kwph_results, paste0(params[[q]]$name, "posthoc.xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = scales::hue_pal()(5)[c(4,2,5,1)])
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Q6.2a.NSCLC_N
q = "Q6.2a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

kw_results <- run_kruskal_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
kw_results

save_results(kw_results, paste0(params[[q]]$name, ".xlsx"))

kwph_results <- run_pairwise_mannwhitney(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
kwph_results

save_results(kwph_results, paste0(params[[q]]$name, "posthoc.xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = scales::hue_pal()(5)[c(3:4,2,5,1)])
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Q6.2a.NSCLC_N_twolevels
q = "Q6.2.1a"

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

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), "Ntwo")
wt_results

save_results(wt_results, paste0(params[[q]]$name, "posthoc.xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = "Ntwo", col_palette = scales::hue_pal()(5)[c(3:4,2,5,1)], legend_lab = "N")
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

# Plot <- alpha_plot(rich_sub_long, x_var = "Ntwo", col_palette = colors_N[c(2,4)], legend_lab = "N")
# Plot
# 
# #ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
# ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, "gradient.svg")), width = 10, height = 8, dpi = 300)

# Q6.3a.NSCLC_M
q = "Q6.3a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = scales::hue_pal()(5)[3:4])
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

#### Q7: Can the microbiome predict the development of postoperative pneumonia? ####

# Q7a.NSCLC_PostopPneumonia
q = "Q7a"

rich_sub <- subset_data(rich, q)
rich_sub_long <- rich_long[rich_long$sample %in% rich_sub$sample,]

rich_sub_long <- rich_sub_long %>%
  mutate(Postop..pneumonia = case_when(
    Postop..pneumonia == 0 ~ "No",
    Postop..pneumonia == 1 ~ "Yes"))
rich_sub_long$Postop..pneumonia <- factor(rich_sub_long$Postop..pneumonia)

wt_results <- run_wilcox_test(rich_sub, c("Shannon", "InvSimpson"), params[[q]]$fixed)
wt_results

save_results(wt_results, paste0(params[[q]]$name, ".xlsx"))

Plot <- alpha_plot(rich_sub_long, x_var = params[[q]]$fixed, col_palette = scales::hue_pal()(2)[2:1], legend_lab = "Pneumonia")
Plot

#ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".tiff")), width = 10, height = 8, dpi = 300)
ggsave(plot = Plot, file.path(a.plots, paste0(params[[q]]$name, ".svg")), width = 10, height = 8, dpi = 300)

#### End ####
