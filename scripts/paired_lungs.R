# Tumour-bearing vs contralateral lung (paired samples of the same patient) #
# Figure 1 A-C, Supplementary Figure 2A, Supplementary Figure 3A and
# Supplementary Table 3. 
# - Fig 1A: within-patient Bray-Curtis distance (tumour-bearing vs contralateral lung) vs a between-patient
#   null distribution (B samples re-paired across patients, 5,000 permutations)
# - Fig 1B / Supp Table 3: per patient, percentile rank of the within-patient
#   distance among the distances from the patient's A sample to all other
#   patients' A samples; one-sample Wilcoxon signed-rank test against 50
# - Supp Fig 2A: within-patient distance vs the closest and a random other
#   patient's A sample; paired Wilcoxon signed-rank tests
# - Fig 1C / Supp Fig 3A: mean relative abundance of the top 12 phyla and top
#   20 genera in A and B (all paired patients / unilateral tumours only)

#### Package setup ####

Sys.setenv(language = "EN")

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(openxlsx)
})

#### End ####

#### Load batch-corrected phyloseq object (2.MMUPHin.R) ####

res.dir <- "Results_direct"

rdata <- new.env()
load(file.path(res.dir, "RData", "6.phyloseq.corrected.RData"), envir = rdata)
physeq_pl <- rdata$physeq_re_bc
lung_col <- rdata$lung_col  # colours of the lungs: Normal and Diseased
rm(rdata)

#### End ####

#### Settings and directories ####

out.dir <- file.path(res.dir, "paired_lungs")
dir.create(out.dir, recursive = T, showWarnings = F)

patient_col <- "Study_Nr"
side_col    <- "Lung.code"  # column with the lung side (A/B)
lungA_value <- "A"  # tumour-bearing (diseased) lung
lungB_value <- "B"  # contralateral (normal) lung
dist_method <- "bray"
n_perm      <- 5000

topN_phylum <- 12
topN_genus  <- 20

#### End ####

#### Plot style ####

# Group colours of the distance plots (not lung sides; those use lung_col)
COL_BLUE <- "#0072B2"  # Okabe-Ito blue
COL_ORNG <- "#D55E00"  # Okabe-Ito vermillion
COL_GRN  <- "#009E73"  # Okabe-Ito green
COL_LINE <- "grey30"   # reference lines and their labels

theme_nature <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title      = element_text(face = "bold", size = base_size + 1),
      plot.subtitle   = element_text(size = base_size),
      axis.title      = element_text(face = "bold"),
      axis.text       = element_text(color = "black"),
      axis.line       = element_line(linewidth = 0.6),
      axis.ticks      = element_line(linewidth = 0.6),
      plot.margin     = margin(10, 18, 10, 10),
      legend.position = "top",
      legend.title    = element_blank(),
      text            = element_text(family = "Helvetica")
    )
}

save_fig <- function(p, base, w, h) {
  ggsave(file.path(out.dir, paste0(base, ".svg")), p, width = w, height = h, dpi = 600, limitsize = FALSE)
}

#### End ####

#### Functions ####

# Directly isolated true samples of the patients with both lungs sampled
select_paired <- function(ps) {
  ps <- subset_samples(ps, Sample_or_Control == "True sample" & Isolation == "Direct_Isolation")
  ps <- prune_samples(sample_sums(ps) > 0, ps)

  sd <- data.frame(sample_data(ps))
  sd <- sd[!is.na(sd[[patient_col]]) & sd[[side_col]] %in% c(lungA_value, lungB_value), ]
  patients <- intersect(sd[[patient_col]][sd[[side_col]] == lungA_value],
                        sd[[patient_col]][sd[[side_col]] == lungB_value])

  ps <- prune_samples(rownames(sd)[sd[[patient_col]] %in% patients], ps)
  prune_taxa(taxa_sums(ps) > 0, ps)
}

# Mean relative abundance per lung side of the top N taxa of one rank
# (top taxa = highest mean over both sides; all others pooled as "Other")
mirror_comp_table <- function(ps, rank_use, topN, other_label = "Other") {
  ps_r <- tax_glom(ps, taxrank = rank_use, NArm = TRUE)
  ps_r <- prune_taxa(taxa_sums(ps_r) > 0, ps_r)

  df <- psmelt(ps_r) %>%
    mutate(rank = as.character(.data[[rank_use]])) %>%
    group_by(side = .data[[side_col]], rank) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop")

  top_taxa <- df %>%
    group_by(rank) %>%
    summarise(global_mean = mean(mean_abund), .groups = "drop") %>%
    arrange(desc(global_mean)) %>%
    slice_head(n = topN) %>%
    pull(rank)

  df2 <- df %>%
    mutate(taxon = ifelse(rank %in% top_taxa, rank, other_label)) %>%
    group_by(side, taxon) %>%
    summarise(mean_abund = sum(mean_abund), .groups = "drop") %>%
    group_by(side) %>%
    mutate(mean_abund = mean_abund / sum(mean_abund)) %>%
    ungroup() %>%
    mutate(value = ifelse(side == lungA_value, -mean_abund, mean_abund))  # A to the left

  tax_order <- df2 %>%
    group_by(taxon) %>%
    summarise(total = sum(mean_abund), .groups = "drop") %>%
    arrange(total) %>%
    pull(taxon)

  df2$taxon <- factor(df2$taxon, levels = tax_order)
  df2
}

plot_mirror <- function(df_mirror, title, x_max) {
  ggplot(df_mirror, aes(x = value, y = taxon, fill = side)) +
    geom_col(width = 0.82, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.6, color = "black") +
    geom_vline(xintercept = seq(-x_max, x_max, by = 0.1), color = "grey88", linewidth = 0.3) +
    scale_fill_manual(values = c("A" = lung_col[["Diseased"]], "B" = lung_col[["Normal"]]),
                      labels = c("A" = "Diseased lung", "B" = "Normal lung")) +
    scale_x_continuous(limits = c(-x_max, x_max),
                       breaks = seq(-x_max, x_max, by = 0.1),
                       labels = function(x) paste0(abs(round(x * 100)), "%")) +
    labs(title = title, x = "Mean relative abundance (%)", y = NULL) +
    theme_nature() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 15),
          plot.title = element_text(size = 15),
          text = element_text(size = 15))
}

#### End ####

#### Paired samples ####

ps_pairs <- select_paired(physeq_pl)
sd <- data.frame(sample_data(ps_pairs)) %>% rownames_to_column("SampleID")

cat("Paired dataset:", nsamples(ps_pairs), "samples;", ntaxa(ps_pairs), "ASVs;",
    length(unique(sd[[patient_col]])), "patients\n")
stopifnot(all(table(sd[[patient_col]]) == 2))

# Bray-Curtis distances between all paired samples
Dm <- as.matrix(phyloseq::distance(ps_pairs, method = dist_method))

# One row per patient: A and B sample and their distance (within-patient)
within_df <- sd %>%
  select(SampleID, patient = all_of(patient_col), lung = all_of(side_col)) %>%
  group_by(patient) %>%
  summarise(A = SampleID[lung == lungA_value][1],
            B = SampleID[lung == lungB_value][1],
            .groups = "drop") %>%
  mutate(dist = Dm[cbind(A, B)])

#### End ####

#### Figure 1A: within-patient vs between-patient (permutation test) ####

# Null distribution: B samples re-paired at random across patients; a
# permutation that pairs a patient with itself is rotated (up to 10 times)
set.seed(123)
between_mean <- replicate(n_perm, {
  perm_idx <- sample(seq_len(nrow(within_df)))
  tries <- 0
  while (any(within_df$patient == within_df$patient[perm_idx]) && tries < 10) {
    perm_idx <- c(perm_idx[-1], perm_idx[1])
    tries <- tries + 1
  }
  mean(Dm[cbind(within_df$A, within_df$B[perm_idx])])
})

obs_mean_within <- mean(within_df$dist)
n_le <- sum(between_mean <= obs_mean_within)
p_perm <- (n_le + 1) / (n_perm + 1)  # one-sided: within < between
d_effect <- (mean(between_mean) - obs_mean_within) / sd(between_mean)

fig1a_summary <- data.frame(
  metric = c("n_patients", "mean_within", "mean_between_null", "n_perm",
             "n_null_means_le_within", "p_perm", "effect_size_null_SD", "dist_method"),
  value = c(nrow(within_df), obs_mean_within, mean(between_mean), n_perm,
            n_le, p_perm, d_effect, dist_method))
print(fig1a_summary)

# Group labels (reviewer: no "A vs B" for the lungs, it is confused with the panel letters)
within_label  <- "Within patient"
between_label <- "Between patients\n(random pairing)"

plot_df <- bind_rows(
  tibble(dist = within_df$dist, group = within_label),
  tibble(dist = between_mean,   group = between_label)
)
plot_df$group <- factor(plot_df$group, levels = c(within_label, between_label))

pA <- ggplot(plot_df, aes(x = group, y = dist, fill = group)) +
  geom_violin(trim = TRUE, linewidth = 0.25, alpha = 0.9) +
  coord_flip() +
  scale_fill_manual(values = setNames(c(COL_BLUE, COL_ORNG), c(within_label, between_label))) +
  labs(subtitle = paste0("Mean within = ", signif(obs_mean_within, 3),
                         " | Mean between = ", signif(mean(between_mean), 3),
                         " | Permutation p ", if (p_perm < 0.001) "< 0.001" else paste0("= ", signif(p_perm, 3))),
       x = NULL, y = "Bray-Curtis distance") +
  theme_nature() +
  theme(legend.position = "none",
        plot.title.position = "plot",  # subtitle starts at the left edge of the figure
        text = element_text(size = 13, face = "bold"),
        axis.text.y = element_text(size = 15))

save_fig(pA, "Fig1A_within_vs_between_distance", w = 6.5, h = 4.8)

#### End ####

#### Supplementary Figure 2A: within-patient vs closest and random other patient ####

lungA_ids <- within_df %>% select(patient, A)

# Closest A sample of another patient
lungA_ids$dist_between_closest <- sapply(seq_len(nrow(lungA_ids)), function(i) {
  min(Dm[lungA_ids$A[i], lungA_ids$A[-i]])
})

# A random A sample of another patient (one draw per patient)
set.seed(123)
lungA_ids$dist_between_random <- sapply(seq_len(nrow(lungA_ids)), function(i) {
  Dm[lungA_ids$A[i], sample(lungA_ids$A[-i], size = 1)]
})

supp2a_df <- within_df %>%
  select(patient, A, B, dist_within = dist) %>%
  left_join(lungA_ids %>% select(patient, dist_between_closest, dist_between_random), by = "patient")

wt_closest <- wilcox.test(supp2a_df$dist_within, supp2a_df$dist_between_closest, paired = TRUE)
wt_random  <- wilcox.test(supp2a_df$dist_within, supp2a_df$dist_between_random,  paired = TRUE)

supp2a_summary <- data.frame(
  metric = c("n_patients", "median_within", "median_closest", "median_random",
             "wilcoxon_V_closest", "wilcoxon_p_closest", "wilcoxon_V_random", "wilcoxon_p_random"),
  value = c(nrow(supp2a_df), median(supp2a_df$dist_within), median(supp2a_df$dist_between_closest),
            median(supp2a_df$dist_between_random),
            unname(wt_closest$statistic), wt_closest$p.value,
            unname(wt_random$statistic), wt_random$p.value))
print(supp2a_summary)

grp_levels <- c("Closest diseased lung\nbetween-patient", "Within-patient", "Random diseased lung\nbetween-patient")
supp2a_long <- supp2a_df %>%
  pivot_longer(c(dist_within, dist_between_closest, dist_between_random), names_to = "group", values_to = "dist") %>%
  mutate(group = factor(recode(group,
                               dist_within = grp_levels[2],
                               dist_between_closest = grp_levels[1],
                               dist_between_random = grp_levels[3]),
                        levels = grp_levels))

pS2A <- ggplot(supp2a_long, aes(x = group, y = dist, group = patient)) +
  geom_line(alpha = 0.25, color = "grey55") +
  geom_point(aes(color = group), size = 2.6) +
  scale_color_manual(values = setNames(c(COL_ORNG, COL_BLUE, COL_GRN), grp_levels)) +
  labs(x = NULL, y = "Bray-Curtis distance") +
  theme_nature() +
  theme(legend.position = "none",
        text = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 16))

# Test results as text to the right of the plot, centred in the free space
stats_text <- paste0("Within-patient vs closest diseased lung between-patient: p = ", signif(wt_closest$p.value, 3),
                     "\n\nWithin-patient vs random diseased lung between-patient: p = ", signif(wt_random$p.value, 3),
                     "\n\nn = ", nrow(supp2a_df), " patients")
text_panel <- ggdraw() +
  draw_label(stats_text, x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
             fontface = "bold", size = 13, fontfamily = "Helvetica", lineheight = 1.1) +
  theme(plot.background = element_rect(fill = "white", colour = NA))  # ggdraw() is transparent by default

pS2A_text <- plot_grid(pS2A, text_panel, nrow = 1, rel_widths = c(1.6, 1))

save_fig(pS2A_text, "SuppFig2A_within_vs_closest_random_distance", w = 16, h = 4.8)

#### End ####

#### Figure 1B and Supplementary Table 3: percentile rank per patient ####

# Percentage of the distances from the patient's A sample to all other patients'
# A samples that are larger than the within-patient distance (100 = the
# patient's own B sample is closer than any other patient)
percentile_df <- supp2a_df %>%
  select(patient, A, B, dist_within) %>%
  mutate(percentile_rank = sapply(seq_len(n()), function(i) mean(Dm[A[i], A[-i]] > dist_within[i]) * 100),
         n_other_patients = n() - 1)

wt_percentile <- wilcox.test(percentile_df$percentile_rank, mu = 50)
n_above_50 <- sum(percentile_df$percentile_rank > 50)

fig1b_summary <- data.frame(
  metric = c("n_patients", "median_percentile_rank", "mean_percentile_rank",
             "n_patients_above_50", "wilcoxon_V_vs_50", "wilcoxon_p_vs_50"),
  value = c(nrow(percentile_df), median(percentile_df$percentile_rank), mean(percentile_df$percentile_rank),
            n_above_50, unname(wt_percentile$statistic), wt_percentile$p.value))
print(fig1b_summary)

# ECDF of the percentile ranks vs the uniform null (no patient effect)
ecdf_df <- percentile_df %>% arrange(percentile_rank) %>% mutate(ecdf_y = row_number() / n())

p_ecdf <- ggplot(ecdf_df, aes(x = percentile_rank, y = ecdf_y)) +
  geom_line(data = tibble(x = c(0, 100), y = c(0, 1)), aes(x = x, y = y),
            linetype = "dashed", linewidth = 0.7, color = COL_LINE) +
  geom_step(color = COL_BLUE, linewidth = 1.0) +
  annotate("text", x = 15, y = 0.92, label = "Expected under\nno patient effect",
           size = 3.6, color = COL_LINE, hjust = 0, fontface = "italic") +
  labs(subtitle = paste0("Median rank = ", round(median(percentile_df$percentile_rank), 1),
                         " | Wilcoxon vs 50: p = ", signif(wt_percentile$p.value, 3),
                         " | n = ", nrow(percentile_df), " patients"),
       x = "Percentile rank of within-patient distance\n(among distances to all other diseased lung samples)",
       y = "Cumulative proportion of patients") +
  theme_nature() +
  theme(text = element_text(size = 13, face = "bold"))

# Percentile rank per patient, sorted
rank_sorted_df <- percentile_df %>%
  arrange(percentile_rank) %>%
  mutate(patient_order = row_number(), above_50 = percentile_rank > 50)

p_rank_sorted <- ggplot(rank_sorted_df, aes(x = patient_order, y = percentile_rank)) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.6, color = COL_LINE) +
  geom_segment(aes(xend = patient_order, y = 50, yend = percentile_rank, color = above_50),
               linewidth = 0.5, alpha = 0.6) +
  geom_point(aes(color = above_50), size = 1.9) +
  scale_color_manual(values = c(`TRUE` = COL_BLUE, `FALSE` = COL_ORNG), guide = "none") +
  labs(subtitle = paste0(n_above_50, " of ", nrow(percentile_df), " patients above the midpoint"),
       x = "Patients (sorted by percentile rank)",
       y = "Percentile rank of within-patient distance") +
  theme_nature() +
  theme(text = element_text(size = 13, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

save_fig(plot_grid(p_ecdf, p_rank_sorted, ncol = 2), "Fig1B_percentile_rank", w = 10.0, h = 4.8)

#### End ####

#### Figure 1C and Supplementary Figure 3A: composition by lung side ####

title_phylum <- paste0("Community composition by lung side (Top ", topN_phylum, " phyla by abundance)")
title_genus  <- paste0("Community composition by lung side (Top ", topN_genus, " genera by abundance)")

# Figure 1C: all patients with both lungs sampled
fig1c_phylum <- mirror_comp_table(ps_pairs, "Phylum", topN_phylum)
fig1c_genus  <- mirror_comp_table(ps_pairs, "Genus", topN_genus)

p_fig1c <- (plot_mirror(fig1c_phylum, title_phylum, x_max = 0.60) | plot_mirror(fig1c_genus, title_genus, x_max = 0.35)) +
  plot_layout(guides = "collect") & theme(legend.position = "top")
save_fig(p_fig1c, "Fig1C_composition_by_lung_side", w = 21, h = 6.5)

# Supplementary Figure 3A: unilateral tumours only (no synchronous tumour)
ps_unilateral <- select_paired(subset_samples(physeq_pl, Synchronous.tumor == "No"))

supp3a_phylum <- mirror_comp_table(ps_unilateral, "Phylum", topN_phylum)
supp3a_genus  <- mirror_comp_table(ps_unilateral, "Genus", topN_genus)

p_supp3a <- (plot_mirror(supp3a_phylum, title_phylum, x_max = 0.60) | plot_mirror(supp3a_genus, title_genus, x_max = 0.35)) +
  plot_layout(guides = "collect") & theme(legend.position = "top")
save_fig(p_supp3a, "SuppFig3A_composition_by_lung_side_unilateral", w = 21, h = 6.5)

# Source data: mean relative abundance per taxon and lung side (A = diseased, B = normal)
write.xlsx(list(Fig1C_phylum = fig1c_phylum, Fig1C_genus = fig1c_genus,
                SuppFig3A_phylum = supp3a_phylum, SuppFig3A_genus = supp3a_genus),
           file.path(out.dir, "composition_by_lung_side.xlsx"))

#### End ####

#### Save statistics and source data ####

write.xlsx(list(Fig1A_permutation = fig1a_summary,
                Fig1B_percentile_rank = fig1b_summary,
                SuppFig2A_wilcoxon = supp2a_summary),
           file.path(out.dir, "paired_lungs_statistics.xlsx"))

# Supplementary Table 3 and the per-patient distances of Supp Fig 2A
write.xlsx(list(SuppTable3_percentile_rank = percentile_df,
                SuppFig2A_distances = supp2a_df),
           file.path(out.dir, "paired_lungs_per_patient.xlsx"))

#### End ####
