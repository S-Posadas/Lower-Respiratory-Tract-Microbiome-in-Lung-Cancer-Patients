# Most abundant taxa per question #
# For every question, at phylum and genus level:
# - boxplots of the relative abundance of the most abundant taxa per group
#   (e.g. Figure 2C, Supp. Figure 5B)
# - bar plots of the mean relative abundance per group (e.g. Supp. Figure 10A/B)
# - stacked bar plots of the mean relative abundance per group
# - ternary plots for questions with three groups (Figure 2D, Supp. Figure 5A)
# - mirrored composition plots for questions with two groups
#   (Supp. Figures 8C/D, 11A)

#### Package setup ####

Sys.setenv(language = "EN")

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(pals)    # glasbey palette
  library(openxlsx)
})

theme_set(theme_bw())

#### End ####

#### Load batch-corrected phyloseq object ####

res.dir <- "Results_direct"

rdata <- new.env()
load(file.path(res.dir, "RData", "6.phyloseq.corrected.RData"), envir = rdata)
physeq_tax0 <- rdata$physeq_re_bc
lung_col <- rdata$lung_col
diagnosis_col <- rdata$diagnosis_col
histology_col <- rdata$histology_col
smoker_col <- rdata$smoker_col
rm(rdata)

# Relative abundance in %
otu_table(physeq_tax0) <- otu_table(physeq_tax0) * 100

#### End ####

#### Settings and directories ####

tax.dir <- file.path(res.dir, "taxa_plots")
dir.create(tax.dir, recursive = T, showWarnings = F)

n_top    <- c(Phylum = 12, Genus = 20)  # boxplots, ternary and mirrored plots
n_top_bar <- 20                         # bar plots (both levels)

#### End ####

#### Plot theme settings ####

text_size = 30
thm <- theme(text = element_text(size = text_size),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "bottom")

# Axis labels with decimal comma and no scientific notation
marks_no_sci <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)

# Taxa colours: glasbey in alphabetical order of the taxa
taxa_colors <- function(taxa) setNames(as.character(pals::glasbey(length(taxa))), sort(unique(taxa)))

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

params <- list(
  Q1.1 = list(name = "Q1.1.all_directly_isolated_lung", fixed = "Lung", col = lung_col),
  Q1.2 = list(name = "Q1.2.main-D_directly_isolated_lung", fixed = "Lung", col = lung_col),
  Q2.1 = list(name = "Q2.1.all_directly_isolated_non_synchronous_lung", fixed = "Lung", col = lung_col),
  Q2.2 = list(name = "Q2.2.main-D_directly_isolated_non_synchronous_lung", fixed = "Lung", col = lung_col),
  Q3.1 = list(name = "Q3.1.diagnosis", fixed = "Diagnosis", col = diagnosis_col),
  Q3.2 = list(name = "Q3.2.diagnosis_nonsmokers", fixed = "Diagnosis", col = diagnosis_col),
  Q3.3 = list(name = "Q3.3.diagnosis_smokers", fixed = "Diagnosis", col = diagnosis_col),
  Q4.1 = list(name = "Q4.1.histology", fixed = "Histology.NSCLC", col = histology_col, lab_plot = "Histology"),
  Q4.2 = list(name = "Q4.2.histology_nonsmokers", fixed = "Histology.NSCLC", col = histology_col, lab_plot = "Histology"),
  Q4.3 = list(name = "Q4.3.histology_smokers", fixed = "Histology.NSCLC", col = histology_col, lab_plot = "Histology"),
  Q5.1 = list(name = "Q5.1.main-D_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking"),
  Q5.2 = list(name = "Q5.2.NSCLC_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking"),
  Q5.3 = list(name = "Q5.3.Benign_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking"),
  Q5.4 = list(name = "Q5.4.Adenocarcinoma_smoking", fixed = "History.of.smoking.y.n", col = smoker_col, lab_plot = "History of smoking"),
  Q6.1 = list(name = "Q6.1.NSCLC_T", fixed = "T", col = scales::hue_pal()(5)[c(4,2,5,1)], lab_plot = "T"),
  Q6.2 = list(name = "Q6.2.NSCLC_N", fixed = "N", col = scales::hue_pal()(5)[c(3:4,2,5,1)]),
  Q6.2.1 = list(name = "Q6.2.NSCLC_N_twolevels", fixed = "Ntwo", col = scales::hue_pal()(5)[c(3:4,2,5,1)], lab_plot = "N"),
  Q6.3 = list(name = "Q6.3.NSCLC_M", fixed = "M", col = c("lightblue", "firebrick")),
  Q6.3.1 = list(name = "Q6.3.1.NSCLC_M_smokers", fixed = "M", col = c("lightblue", "firebrick"))
)

# N stage in two levels: N0 vs N1-3
sample_data(physeq_tax0)$Ntwo <- factor(ifelse(sample_data(physeq_tax0)$N == 0, "N0", "N1-3"), levels = c("N0", "N1-3"))

# Select the samples of one question (as in beta.R)
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

#### End ####

#### Most abundant taxa per question ####

for (q in names(params)) {

  print(q)
  p <- params[[q]]
  legend_title <- if (is.null(p$lab_plot)) p$fixed else p$lab_plot

  out.dir <- file.path(tax.dir, p$name)
  dir.create(out.dir, recursive = T, showWarnings = F)

  physeq_q <- subset_physeq(physeq_tax0, q)
  grp <- droplevels(factor(sample_data(physeq_q)[[p$fixed]]))
  # T, N and M stages labelled as N0, N1, ...
  if (p$fixed %in% c("T", "N", "M")) levels(grp) <- paste0(p$fixed, levels(grp))
  sample_data(physeq_q)[[p$fixed]] <- grp
  groups <- levels(sample_data(physeq_q)[[p$fixed]])

  source_data <- list()
  mirror_plots <- list()

  for (tax_level in c("Phylum", "Genus")) {

    # Relative abundance (%) per sample and taxon
    melted <- psmelt(tax_glom(physeq_q, tax_level, NArm = TRUE)) %>%
      mutate(taxon = as.character(.data[[tax_level]]), group = .data[[p$fixed]])

    # Taxa ranked by their total relative abundance over all samples of the question
    ranking <- melted %>% group_by(taxon) %>% summarise(total = sum(Abundance)) %>% arrange(desc(total))
    top_taxa <- head(ranking$taxon, n_top[[tax_level]])
    top_taxa_bar <- head(ranking$taxon, n_top_bar)

    # Mean relative abundance per group (%), scaled to 100 % per group
    group_means <- melted %>%
      group_by(group, taxon) %>%
      summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
      group_by(group) %>%
      mutate(mean_abund = 100 * mean_abund / sum(mean_abund)) %>%
      ungroup()
    source_data[[tax_level]] <- group_means %>%
      pivot_wider(names_from = group, values_from = mean_abund) %>%
      arrange(match(taxon, ranking$taxon))

    ## Boxplots (log scale)
    box_df <- melted %>% filter(taxon %in% top_taxa, Abundance > 0)
    n_col <- ifelse(tax_level == "Phylum", 6, 5)
    n_row <- ceiling(length(top_taxa) / n_col)
    p_box <- ggplot(box_df, aes(x = group, y = Abundance, color = group)) +
      geom_boxplot(fill = NA, outlier.shape = NA) +  # all samples shown as points
      geom_point(size = 1.5, alpha = 0.6, position = position_jitter(width = 0.2, height = 0, seed = 1)) +
      facet_wrap(~ taxon, ncol = n_col) +
      scale_color_manual(values = p$col, drop = FALSE) +
      scale_x_discrete(drop = FALSE) +
      scale_y_log10(labels = marks_no_sci) +
      labs(x = NULL, y = "Relative abundance (%)", color = legend_title) +
      thm
    # Phylum: wide and low, the size of Figure 2C (1907 x 741 px in Inkscape = 19.87 x 7.72 in
    # for two rows of panels)
    box_width  <- ifelse(tax_level == "Phylum", 19.87, 18)
    box_height <- ifelse(tax_level == "Phylum", 7.72 / 2 * n_row, 3 + 6.5 * n_row)
    ggsave(file.path(out.dir, paste0(tax_level, "_boxplot.svg")), p_box, width = box_width, height = box_height,
           dpi = 300, limitsize = FALSE)

    ## Bar plots: mean relative abundance per group (pseudo-log scale)
    bar_df <- group_means %>% filter(taxon %in% top_taxa_bar)
    p_bar <- ggplot(bar_df, aes(x = taxon, y = mean_abund, fill = taxon)) +
      geom_col(color = "black", linewidth = 0.2) +
      facet_wrap(~ group, nrow = 1) +
      scale_fill_manual(values = taxa_colors(bar_df$taxon), name = tax_level) +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 1, 5, 10, 20, 30, 100), labels = marks_no_sci) +
      labs(x = NULL, y = "Relative abundance %") +
      theme(text = element_text(size = 22),
            axis.text = element_text(size = 18),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.text = element_text(size = 22),
            legend.title = element_text(size = 24))
    ggsave(file.path(out.dir, paste0(tax_level, "_barplot.svg")), p_bar, width = 13, height = 8, dpi = 300)
    
    ## Stacked bar plots: mean relative abundance per group
    stack_df <- group_means %>%
      mutate(taxon = ifelse(taxon %in% top_taxa_bar, taxon, "Other")) %>%
      group_by(group, taxon) %>%
      summarise(mean_abund = sum(mean_abund), .groups = "drop") %>%
      mutate(taxon = factor(taxon, levels = c(sort(setdiff(unique(taxon), "Other")), "Other")))
    p_stack <- ggplot(stack_df, aes(x = group, y = mean_abund, fill = taxon)) +
      geom_col(color = "black", linewidth = 0.2, width = 0.8) +
      scale_fill_manual(values = c(taxa_colors(setdiff(levels(stack_df$taxon), "Other")), Other = "grey80"),
                        name = tax_level) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.02)), labels = marks_no_sci) +
      labs(x = NULL, y = "Relative abundance (%)") +
      theme(text = element_text(size = 22),
            axis.text = element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 22)) +
      guides(fill = guide_legend(ncol = 1))
    ggsave(file.path(out.dir, paste0(tax_level, "_stacked_barplot.svg")), p_stack, width = 9, height = 9, dpi = 300)

    ## Mirrored composition (questions with two groups): first group left
    if (length(groups) == 2) {
      mirror_df <- group_means %>%
        mutate(taxon = ifelse(taxon %in% top_taxa, taxon, "Other")) %>%
        group_by(group, taxon) %>%
        summarise(mean_abund = sum(mean_abund) / 100, .groups = "drop") %>%
        mutate(value = ifelse(group == groups[1], -mean_abund, mean_abund))
      tax_order <- mirror_df %>% group_by(taxon) %>% summarise(t = sum(mean_abund)) %>% arrange(t) %>% pull(taxon)
      mirror_df$taxon <- factor(mirror_df$taxon, levels = tax_order)
      x_max <- ceiling(max(mirror_df$mean_abund) * 10) / 10

      mirror_plots[[tax_level]] <- ggplot(mirror_df, aes(x = value, y = taxon, fill = group)) +
        geom_col(width = 0.82, color = "white", linewidth = 0.2) +
        geom_vline(xintercept = 0, linewidth = 0.6, color = "black") +
        geom_vline(xintercept = seq(-x_max, x_max, by = 0.1), color = "grey88", linewidth = 0.3) +
        scale_fill_manual(values = p$col, name = NULL) +
        scale_x_continuous(limits = c(-x_max, x_max), breaks = seq(-x_max, x_max, by = 0.1),
                           labels = function(x) paste0(abs(round(x * 100)), "%")) +
        labs(title = paste0("Community composition (Top ", n_top[[tax_level]], " ",
                            ifelse(tax_level == "Phylum", "phyla", "genera"), " by abundance)"),
             x = "Mean relative abundance (%)", y = NULL) +
        ggplot2::theme_classic(base_size = 15) +
        theme(legend.position = "top", plot.title = element_text(size = 15, face = "bold"))
    }
  }

  if (length(mirror_plots) == 2) {
    p_mirror <- (mirror_plots$Phylum | mirror_plots$Genus) + plot_layout(guides = "collect") & theme(legend.position = "top")
    ggsave(file.path(out.dir, "mirrored_phylum_genus.svg"), p_mirror, width = 21, height = 6.5, dpi = 300)
  }

  # Source data: mean relative abundance (%) per group of all taxa
  write.xlsx(source_data, file.path(out.dir, "mean_relative_abundance.xlsx"))
}

#### End ####

#### Ternary plots (questions with three groups) ####

# Mean relative abundance (%) per group of the most abundant taxa. The taxa are
# ranked by the sum of the three group means (each group weighs the same, so
# small groups count as much as large ones).
# ggtern is loaded only here: once loaded, it changes how every plot is drawn.

library(ggtern)

for (q in names(params)) {
  p <- params[[q]]
  physeq_q <- subset_physeq(physeq_tax0, q)
  groups <- levels(droplevels(factor(sample_data(physeq_q)[[p$fixed]])))
  if (length(groups) != 3) next
  print(q)
  
  for (tax_level in c("Phylum", "Genus")) {
    tern_df <- psmelt(tax_glom(physeq_q, tax_level, NArm = TRUE)) %>%
      mutate(taxon = as.character(.data[[tax_level]]), group = as.character(.data[[p$fixed]])) %>%
      group_by(group, taxon) %>%
      summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
      pivot_wider(names_from = group, values_from = mean_abund, values_fill = 0) %>%
      mutate(total = rowSums(across(all_of(groups)))) %>%
      arrange(desc(total)) %>%
      slice_head(n = n_top[[tax_level]])
    
    p_tern <- ggtern(tern_df, aes(x = .data[[groups[2]]], y = .data[[groups[3]]], z = .data[[groups[1]]])) +
      geom_point(aes(size = total, color = taxon), alpha = 0.8) +
      scale_size_continuous(range = c(3, 8)) +
      scale_color_manual(values = taxa_colors(tern_df$taxon)) +
      labs(x = groups[2], y = groups[3], z = groups[1], size = "Relative abundance (%)", color = tax_level) +
      theme_bw() +
      theme(text = element_text(size = 15, face = "bold"),
            legend.text = element_text(face = "italic")) +
      guides(color = guide_legend(override.aes = list(size = 6), ncol = if (tax_level == "Genus") 2 else 1))
    ggsave(file.path(tax.dir, p$name, paste0(tax_level, "_ternary.svg")), p_tern, width = 9.5, height = 8, dpi = 300)
  }
}

#### End ####
