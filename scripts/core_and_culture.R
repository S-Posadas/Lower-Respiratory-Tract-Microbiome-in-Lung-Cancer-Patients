# Core airway microbiome and culture enrichment #
# Supplementary Figure 2 B-E and Figure 3E.
## - Supp Fig 2B: core genera of the tumour-bearing lung (directly isolated
#   samples of the patients with both lungs sampled): prevalence vs mean
#   relative abundance; core = prevalence >= 20 % and mean >= 0.5 %
# - Supp Fig 2C: Bray-Curtis distance (genus level) between the directly
#   isolated and the culture-enriched sample of the same lung, vs a
#   between-patient null (culture samples re-paired across patients within
#   lung side, 5,000 permutations)
# - Supp Fig 2D: prevalence of the 30 most prevalent genera of the direct
#   samples in direct vs culture-enriched samples (tumour-bearing lung)
# - Supp Fig 2E: mean relative abundance of the core genera (2B) in direct vs
#   culture-enriched samples
# - Figure 3E: prevalence of Streptococcus mitis by N stage in direct vs
#   culture-enriched samples

#### Package setup ####

Sys.setenv(language = "EN")

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(scales)
  library(openxlsx)
})

#### End ####

#### Load data ####

res.dir <- "Results_merged"

# Direct + culture-enriched samples.
physeq_merged <- readRDS(file.path(res.dir, "RData", "physeq_re_bc.rds"))

#### End ####

#### Settings and directories ####

out.dir <- file.path(res.dir, "core_and_culture")
dir.create(out.dir, recursive = T, showWarnings = F)

core_prevalence <- 0.20   # core genus: present in >= 20 % of the patients
core_abundance  <- 0.005  # and mean relative abundance >= 0.5 %
n_perm          <- 5000
top_recovery    <- 30     # genera shown in the recovery plot

#### End ####

#### Plot style ####

COL_BLUE    <- "#0072B2"
COL_ORANGE  <- "#E69F00"
COL_RED     <- "#D55E00"
COL_GREY    <- "grey20"
COL_DIRECT  <- "#1b9e77"  # direct sequencing (Supp Fig 2E)
COL_CULTURE <- "#d95f02"  # culture enriched (Supp Fig 2E)

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

#### Supplementary Figure 2B: core genera of the tumour-bearing lung ####

# Directly isolated true samples of the patients with both lungs sampled;
# tumour-bearing lung (A) only, so each patient counts once
ps <- subset_samples(physeq_merged, Sample_or_Control == "True sample" & Isolation == "Direct_Isolation")
ps <- prune_samples(sample_sums(ps) > 0, ps)
sd <- data.frame(sample_data(ps))
patients_paired <- intersect(sd$Study_Nr[sd$Lung.code == "A"], sd$Study_Nr[sd$Lung.code == "B"])

ps_A <- prune_samples(sd$Lung.code == "A" & sd$Study_Nr %in% patients_paired, ps)
ps_A <- prune_taxa(taxa_sums(ps_A) > 0, ps_A)
cat("Supp Fig 2B:", nsamples(ps_A), "tumour-bearing lung samples\n")

# Mean relative abundance and prevalence per genus
ps_A_genus <- tax_glom(ps_A, taxrank = "Genus", NArm = TRUE)
ps_A_genus <- prune_taxa(taxa_sums(ps_A_genus) > 0, ps_A_genus)
otu_A <- as(otu_table(ps_A_genus), "matrix")
if (!taxa_are_rows(ps_A_genus)) otu_A <- t(otu_A)

genus_tbl <- tibble(genus = as.character(tax_table(ps_A_genus)[, "Genus"]),
                    mean_rel_abundance = rowMeans(otu_A),
                    prevalence = rowMeans(otu_A > 0)) %>%
  arrange(desc(mean_rel_abundance))

core_tbl <- genus_tbl %>%
  filter(prevalence >= core_prevalence, mean_rel_abundance >= core_abundance)
cat("Core genera:", nrow(core_tbl), "\n")

p_core <- ggplot(genus_tbl, aes(x = prevalence, y = mean_rel_abundance)) +
  geom_point(alpha = 0.65, size = 2.0, color = COL_BLUE) +
  geom_hline(yintercept = core_abundance, linetype = "dashed", linewidth = 0.5, color = "grey30") +
  geom_vline(xintercept = core_prevalence, linetype = "dashed", linewidth = 0.5, color = "grey30") +
  geom_text_repel(data = core_tbl, aes(label = genus), size = 4, box.padding = 0.5,
                  point.padding = 0.3, force = 2, max.overlaps = 40, segment.color = "grey50") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Prevalence across patients", y = "Mean relative abundance") +
  theme_nature() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"))

save_fig(p_core, "SuppFig2B_core_genera", w = 7.8, h = 5.2)

#### End ####

#### Patients with all four samples (direct and culture-enriched, both lungs) ####

sd_m <- data.frame(sample_data(physeq_merged)) %>% rownames_to_column("sample_id")
sd_m <- sd_m %>% filter(is.na(Sample_or_Control) | Sample_or_Control == "True sample")

patients_2x2 <- sd_m %>%
  filter(Lung.code %in% c("A", "B"), Isolation %in% c("Direct_Isolation", "Culture_Enriched")) %>%
  distinct(Study_Nr, Lung.code, Isolation) %>%
  count(Study_Nr) %>%
  filter(n == 4) %>%  # A and B, each direct and culture-enriched
  pull(Study_Nr)
cat("Patients with all four samples:", length(patients_2x2), "\n")

ps_22 <- prune_samples(sd_m$sample_id[sd_m$Study_Nr %in% patients_2x2], physeq_merged)
ps_22 <- prune_taxa(taxa_sums(ps_22) > 0, ps_22)
ps_22 <- prune_samples(sample_sums(ps_22) > 0, ps_22)

# Genus-level relative abundance matrix (samples x genera) and Bray-Curtis
ps_22_genus <- tax_glom(ps_22, taxrank = "Genus", NArm = TRUE)
taxa_names(ps_22_genus) <- tax_table(ps_22_genus)[, "Genus"]
Xg <- as(otu_table(ps_22_genus), "matrix")
if (taxa_are_rows(ps_22_genus)) Xg <- t(Xg)

meta_g <- data.frame(sample_data(ps_22_genus)) %>% rownames_to_column("sample_id")
Xg <- Xg[meta_g$sample_id, , drop = FALSE]

# Bray-Curtis for Supp Fig 2C.
ps_22_genus_all <- tax_glom(ps_22, taxrank = "Genus", NArm = FALSE)
X_all <- as(otu_table(ps_22_genus_all), "matrix")
if (taxa_are_rows(ps_22_genus_all)) X_all <- t(X_all)
Dmat <- as.matrix(vegdist(X_all[meta_g$sample_id, , drop = FALSE], method = "bray"))

#### End ####

#### Supplementary Figure 2C: direct vs culture-enriched, same lung (permutation test) ####

# One row per patient and lung: the direct and the culture-enriched sample
pairs_same_lung <- meta_g %>%
  select(sample_id, Study_Nr, Lung.code, Isolation) %>%
  pivot_wider(id_cols = c(Study_Nr, Lung.code), names_from = Isolation, values_from = sample_id) %>%
  filter(!is.na(Direct_Isolation), !is.na(Culture_Enriched)) %>%
  mutate(within_dist = Dmat[cbind(Direct_Isolation, Culture_Enriched)])

# Null distribution per lung side: culture samples re-paired across patients
set.seed(1)
null_df <- pairs_same_lung %>%
  group_by(Lung.code) %>%
  group_modify(~ tibble(null_mean = replicate(n_perm, {
    mean(Dmat[cbind(.x$Direct_Isolation, sample(.x$Culture_Enriched))])
  }))) %>%
  ungroup()

obs_mean <- mean(pairs_same_lung$within_dist)
n_le <- sum(null_df$null_mean <= obs_mean)
p_perm_str <- if (n_le == 0) paste0("< ", 1 / n_perm) else paste0("= ", signif(n_le / nrow(null_df), 3))

supp2c_summary <- data.frame(
  metric = c("n_patients", "n_pairs", "mean_within", "mean_between_null", "n_perm_per_lung",
             "n_null_means_le_within", "p_perm"),
  value = c(length(unique(pairs_same_lung$Study_Nr)), nrow(pairs_same_lung), obs_mean,
            mean(null_df$null_mean), n_perm, n_le, n_le / nrow(null_df)))
print(supp2c_summary)

plot_df <- bind_rows(
  tibble(group = "Within-patient\n(Direct vs Culture;\nsame lung)", value = pairs_same_lung$within_dist),
  tibble(group = "Between-patient null\n(mean over pairs)", value = null_df$null_mean)
) %>% mutate(group = factor(group, levels = unique(group)))

p_perm <- ggplot(plot_df, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = TRUE, linewidth = 0.4, color = COL_GREY) +
  geom_point(data = tibble(group = levels(plot_df$group)[1], value = obs_mean),
             aes(x = group, y = value), inherit.aes = FALSE, size = 3.2, color = COL_RED) +
  scale_fill_manual(values = c(COL_ORANGE, COL_BLUE)) +
  labs(title = "Culture enrichment preserves patient identity",
       subtitle = paste0("Permutation test: p ", p_perm_str, " \n Observed mean = ", round(obs_mean, 3)),
       x = NULL, y = "Bray–Curtis distance (genus-level)") +
  theme_nature(15) +
  theme(legend.position = "none", panel.grid = element_blank())

save_fig(p_perm, "SuppFig2C_direct_vs_culture_permutation", w = 7.4, h = 4.9)

#### End ####

#### Supplementary Figure 2D: recovery of common genera after culture (tumour-bearing lung) ####

lungA <- meta_g %>% filter(Lung.code == "A")
XA <- Xg[lungA$sample_id, , drop = FALSE]
is_direct <- lungA$Isolation == "Direct_Isolation"

# The most prevalent genera in the direct samples
prev_direct_all <- colMeans(XA[is_direct, , drop = FALSE] > 0)
top_genera <- names(sort(prev_direct_all, decreasing = TRUE))[1:top_recovery]

recovery_tbl <- tibble(genus = top_genera,
                       prev_direct = colMeans(XA[is_direct, top_genera, drop = FALSE] > 0),
                       prev_culture = colMeans(XA[!is_direct, top_genera, drop = FALSE] > 0)) %>%
  mutate(delta_prev = prev_culture - prev_direct) %>%
  arrange(desc(prev_culture))

# Labels: the 7 genera most often recovered and the 3 with the largest loss
label_genera <- unique(c(recovery_tbl %>% slice(1:7) %>% pull(genus),
                         recovery_tbl %>% arrange(delta_prev) %>% slice(1:3) %>% pull(genus)))

p_rec_scatter <- ggplot(recovery_tbl, aes(x = prev_direct, y = prev_culture)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, color = "grey80") +
  geom_point(size = 2.2, alpha = 0.85, color = COL_ORANGE) +
  geom_text_repel(data = recovery_tbl %>% filter(genus %in% label_genera), aes(label = genus),
                  size = 4, box.padding = 0.4, point.padding = 0.25, max.overlaps = 200,
                  segment.color = "grey60") +
  labs(x = "Prevalence in Direct (fraction > 0)",
       y = "Prevalence in Culture-enriched (fraction > 0)") +
  theme_nature(15) +
  theme(panel.grid = element_blank())

save_fig(p_rec_scatter, "SuppFig2D_recovery_scatter", w = 6.9, h = 5.4)

#### End ####

#### Supplementary Figure 2E: core genera in direct vs culture-enriched samples ####

# Mean relative abundance of each core genus per isolation method (share of
# the whole community; the core genera together are about 70 % of it)
ps_core <- prune_taxa(taxa_sums(ps_22) > 0, subset_samples(ps_22, LibrarySizeDecontam > 0))
ps_core <- prune_taxa(tax_table(ps_core)[, "Genus"] %in% core_tbl$genus, ps_core)
ps_core_genus <- tax_glom(ps_core, taxrank = "Genus", NArm = TRUE)

core_comp <- psmelt(ps_core_genus) %>%
  group_by(Isolation, taxon = as.character(Genus)) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(value = ifelse(Isolation == "Direct_Isolation", -mean_abund, mean_abund))  # direct to the left

tax_order <- core_comp %>% group_by(taxon) %>% summarise(total = sum(mean_abund)) %>% arrange(total) %>% pull(taxon)
core_comp$taxon <- factor(core_comp$taxon, levels = tax_order)

x_max <- 0.5
p_core_comp <- ggplot(core_comp, aes(x = value, y = taxon, fill = Isolation)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.6, color = "black") +
  geom_vline(xintercept = seq(-x_max, x_max, by = 0.1), color = "grey88", linewidth = 0.3) +
  scale_fill_manual(values = c("Direct_Isolation" = COL_DIRECT, "Culture_Enriched" = COL_CULTURE),
                    labels = c("Direct_Isolation" = "Direct sequencing", "Culture_Enriched" = "Culture enriched")) +
  scale_x_continuous(limits = c(-x_max, x_max), breaks = seq(-x_max, x_max, by = 0.1),
                     labels = function(x) paste0(abs(round(x * 100)), "%")) +
  labs(x = "Mean relative abundance (%)", y = NULL) +
  theme_nature() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15),
        text = element_text(size = 15))

save_fig(p_core_comp, "SuppFig2E_core_genera_direct_vs_culture", w = 10.5, h = 6)

#### End ####

#### Save statistics and source data ####

write.xlsx(list(SuppFig2B_genera = genus_tbl,
                SuppFig2B_core_genera = core_tbl,
                SuppFig2C_permutation = supp2c_summary,
                SuppFig2C_pairs = pairs_same_lung,
                SuppFig2D_recovery = recovery_tbl,
                SuppFig2E_core_composition = core_comp),
           file.path(out.dir, "core_and_culture.xlsx"))

#### End ####

#### Figure 3E: Streptococcus mitis by N stage, direct vs culture-enriched ####

# Figure 3E: Streptococcus mitis in direct vs culture-enriched samples by N
# stage (NSCLC, diseased lung). A sample counts as positive if
# any ASV with "mitis" in its species name is present. Only lungs with both a
# direct and a culture-enriched sample are counted (LML-093 and LML-110 have no
# reads in the direct sample and 1 read in the culture sample, so they are out).
ps_nsclc <- subset_samples(physeq_merged, Sample_or_Control == "True sample" & Lung == "Diseased" & Diagnosis == "NSCLC")
tax_m <- as.data.frame(tax_table(ps_nsclc))
mitis_asvs <- rownames(tax_m)[grepl("Streptococcus", tax_m$Genus) & grepl("mitis", tax_m$Species)]
otu_m <- as(otu_table(ps_nsclc), "matrix")
if (!taxa_are_rows(ps_nsclc)) otu_m <- t(otu_m)

mitis_df <- data.frame(sample_data(ps_nsclc)) %>%
  mutate(presence = colSums(otu_m[mitis_asvs, , drop = FALSE]) > 0,
         Isolation = as.character(Isolation), N = as.character(N)) %>%
  group_by(Study_Nr) %>%
  filter(n_distinct(Isolation) == 2) %>%   # only lungs with both a direct and a culture sample
  ungroup() %>%
  select(Isolation, N, presence) %>%
  filter(!is.na(N)) %>%
  mutate(Isolation = factor(recode(Isolation, Direct_Isolation = "Direct isolation", Culture_Enriched = "Culture enriched"),
                            levels = c("Direct isolation", "Culture enriched")),
         N = paste0("N", N))

fig3e <- mitis_df %>%
  group_by(N, Isolation) %>%
  summarise(n_samples = n(), n_positive = sum(presence), prevalence = 100 * n_positive / n_samples, .groups = "drop")
write.xlsx(list(prevalence = fig3e, mitis_ASVs = data.frame(ASV = mitis_asvs, Species = tax_m[mitis_asvs, "Species"])),
           file.path(out.dir, "Fig3E_S_mitis_prevalence.xlsx"))

p_fig3e <- ggplot(fig3e, aes(x = N, y = prevalence, fill = Isolation)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_text(aes(label = paste0(n_positive, "/", n_samples)), position = position_dodge(0.9), vjust = -0.3, size = 6) +
  scale_fill_manual(values = c("Direct isolation" = "#1b9e77", "Culture enriched" = "#d95f02")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Streptococcus mitis", x = NULL, y = "Prevalence (%)", fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(size = 20),
        plot.title = element_text(face = "italic", size = 20))
ggsave(file.path(out.dir, "Fig3E_S_mitis_prevalence_by_N.svg"), p_fig3e, width = 6, height = 6, dpi = 300)

#### End ####
