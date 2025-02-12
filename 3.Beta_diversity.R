
# Beta diversity #

#### Package setup ####

Sys.setenv(language = "EN")
library(phyloseq)
library(ggplot2)
library(ggordiplots)
library(plyr)
library(gridExtra)
library(openxlsx)
library(cowplot)
theme_set(theme_bw())

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

res.dir <- "Results_trimmed_cutadapt3_20241216"
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir,"5.phyloseq.filtered.RData"))

# Create directories for results
b.stats <- file.path(res.dir,"4.Beta_stats")
b.plots <- file.path(res.dir,"5.Beta_plots")

dir.create(b.stats, recursive = T)
dir.create(b.plots, recursive = T)

#### End ####

#### Functions ####

plot <- function(phy,dist,var, palette = NULL, legend_title = NULL, nrow = NULL){
  if(dist == "bray"){NM = nmds_bray ; PC = pcoa_bray}
  else if(dist == "wUF"){NM = nmds_wunifrac ; PC = pcoa_wunifrac}
  else if(dist == "uwUF"){NM = nmds_uunifrac ; PC = pcoa_uunifrac}
  
  NMDS = plot_ordination(phy, NM, color=var) + geom_point(size=8, alpha=0.5) +
    stat_ellipse(type = "t") + thm + #scale_color_manual(values = palette) +
    annotate("text", x = max(scores(NM)[,"NMDS1"]), y = max(scores(NM)[,"NMDS2"]),
             label = paste("Stress:", round(NM$stress, 3)), hjust = 1, vjust = 1, size = 7) 
  # Add color palette only if provided
  if (!is.null(palette)) {
    NMDS <- NMDS + scale_color_manual(values = palette)
  }

  if (!is.null(legend_title)) {
    NMDS <- NMDS + labs(color = legend_title)
  } 
  
  if (!is.null(nrow)) {
    NMDS <- NMDS + guides(color = guide_legend(nrow = nrow))
  } 

  PCOA = plot_ordination(phy, PC, color=var) + geom_point(size=8, alpha=0.5) +
    stat_ellipse(type = "t") + thm #+ scale_color_manual(values = palette)
    # Add color palette only if provided
  if (!is.null(palette)) {
    PCOA <- PCOA + scale_color_manual(values = palette)
  }
  
  if (!is.null(legend_title)) {
    PCOA <- PCOA + labs(color = legend_title)
  } 
  
  if (!is.null(nrow)) {
    PCOA <- PCOA + guides(color = guide_legend(nrow = nrow))
  } 
  
  NP <- grid.arrange(NMDS + ggtitle("NMDS") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)),
                     PCOA + ggtitle("PCoA") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)),
                     nrow=1)
  
  NP <- grid.arrange(NP, get_plot_component(PCOA, pattern = "guide-box-bottom"), nrow= 2, heights = c(10,1))
  
}

#### End ####

#### Set theme ####
text_size = 30
thm <- theme(text = element_text(size = text_size),
             #  axis.text.x = element_blank(),
             # axis.ticks.x = element_blank(),
             legend.position = "bottom")

thm.x <- theme(text = element_text(size = text_size),
               axis.text.x = element_text(angle = -30, hjust = 0),
               legend.position = "bottom")

#### All directly isolated samples ####

physeq_beta_0 = physeq_re

# Add LibrarySizeGroup
sample_data(physeq_beta_0)$LibrarySizeGroup <- ifelse(sample_data(physeq_beta_0)$LibrarySizeDecontam >= 1000, ">=1000", "<1000")

# Remove controls
physeq_beta_0 = subset_samples(physeq_beta_0, sample_names(physeq_beta_0) %in% c(rownames(sample_info[!sample_info$Sample_or_Control %in% c("Negative control sample", "Positive control sample"),]))) 

# Remove samples with no reads
physeq_beta_0 = subset_samples(physeq_beta_0, !is.na(colSums(otu_table(physeq_beta_0)))) 

# Keep only directly isolated samples
physeq_beta_0 = subset_samples(physeq_beta_0, sample_data(physeq_beta_0)$Isolation == "Direct_Isolation") 

# Remove outliers
#physeq_beta_0 = subset_samples(physeq_beta_0, !sample_names(physeq_beta_0) %in% c("LML_088_A1", "LML_110_B2")) 

# Select main diagnosis
physeq_beta_0 = subset_samples(physeq_beta_0, sample_data(physeq_beta_0)$Diagnosis %in% c("NSCLC", "SCLC", "Benign"))

# Remove not present taxa
physeq_beta_0 <- prune_taxa(taxa_sums(physeq_beta_0) != 0, physeq_beta_0)


#### All directly isolated samples - Lung ####

physeq_beta = physeq_beta_0

# Remove outliers
#physeq_beta = subset_samples(physeq_beta, !sample_names(physeq_beta) %in% c("LML_088_A1", "LML_110_B2")) 

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
set.seed(123)
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
# An outlier is detected in NMDS dimension 1.
# The outlier in NMDS dimension 1 is sample: LML_088_A1 
# An outlier is detected in NMDS dimension 2.
# The outlier in NMDS dimension 2 is sample: LML_110_B2 
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
# An outlier is detected in NMDS dimension 1.
# The outlier in NMDS dimension 1 is sample: LML_088_A1 
# An outlier is detected in NMDS dimension 2.
# The outlier in NMDS dimension 2 is sample: LML_073_B2 
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)
# No outliers detected in NMDS dimension 1.
# An outlier is detected in NMDS dimension 2.
# The outlier in NMDS dimension 2 is sample: LML_108_B2 

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin") 
  test.adonis2 <- adonis2(dists[[i]] ~ History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis3 <- adonis2(dists[[i]] ~ History.of.smoking.y.n + Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis4 <- adonis2(dists[[i]] ~ Lung + History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis5 <- adonis2(dists[[i]] ~ Diagnosis, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis6 <- adonis2(dists[[i]] ~ Diagnosis + Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis7 <- adonis2(dists[[i]] ~ Lung + Diagnosis, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis8 <- adonis2(dists[[i]] ~ History.of.smoking.y.n + Diagnosis + Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis9 <- adonis2(dists[[i]] ~ Batch, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 

  ad[[i]] <- rbind(test.adonis1, test.adonis2, test.adonis3, test.adonis4,
                   test.adonis5, test.adonis6, test.adonis7, test.adonis8, test.adonis9)
}
ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "bothlungs.xlsx"), rowNames = T)

# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "Lung", lung_col)

ggsave(plot = NP, file.path(b.plots, "Lung_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Lung_bray.svg"), width = 15, height = 8, dpi = 300)

NP <- plot(physeq_beta, dist = "bray", var = "Batch", palette)

ggsave(plot = NP, file.path(b.plots, "Batch_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Batch_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "Lung", lung_col)

ggsave(plot = NP, file.path(b.plots, "Lung_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Lung_wUF.svg"), width = 15, height = 8, dpi = 300)

NP <- plot(physeq_beta, dist = "wUF", var = "Batch", palette)

ggsave(plot = NP, file.path(b.plots, "Batch_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Batch_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "Lung", lung_col)

ggsave(plot = NP, file.path(b.plots, "Lung_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Lung_uwUF.svg"), width = 15, height = 8, dpi = 300)

NP <- plot(physeq_beta, dist = "uwUF", var = "Batch", palette)

ggsave(plot = NP, file.path(b.plots, "Batch_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Batch_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### All directly isolated samples without non synchronous tumor - Lung ####

physeq_beta = physeq_beta_0

# Remove non synchronous tumor
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Synchronous.tumor == "No")

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
set.seed(123)
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin") 
  test.adonis2 <- adonis2(dists[[i]] ~ History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis3 <- adonis2(dists[[i]] ~ History.of.smoking.y.n + Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis4 <- adonis2(dists[[i]] ~ Lung + History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis5 <- adonis2(dists[[i]] ~ Diagnosis, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis6 <- adonis2(dists[[i]] ~ Diagnosis + Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis7 <- adonis2(dists[[i]] ~ Lung + Diagnosis, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis8 <- adonis2(dists[[i]] ~ History.of.smoking.y.n + Diagnosis + Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)

  ad[[i]] <- rbind(test.adonis1, test.adonis2, test.adonis3, test.adonis4,
                   test.adonis5, test.adonis6, test.adonis7, test.adonis8)
}

ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "bothlungs_non_synchronous.xlsx"), rowNames = T)

# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "Lung", lung_col)

ggsave(plot = NP, file.path(b.plots, "Lung_non_synchronous_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Lung_non_synchronous_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "Lung", lung_col)

ggsave(plot = NP, file.path(b.plots, "Lung_non_synchronous_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Lung_non_synchronous_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "Lung", lung_col)

ggsave(plot = NP, file.path(b.plots, "Lung_non_synchronous_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Lung_non_synchronous_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### Directly isolated samples NSCLC ####

physeq_beta = physeq_beta_0

# Select NSCLC
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Diagnosis == "NSCLC")

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
set.seed(123)
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
  test.adonis2 <- adonis2(dists[[i]] ~ History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis3 <- adonis2(dists[[i]] ~ History.of.smoking.y.n + Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis4 <- adonis2(dists[[i]] ~ Lung + History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis5 <- adonis2(dists[[i]] ~ Histology.NSCLC, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis6 <- adonis2(dists[[i]] ~ Histology.NSCLC + Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis7 <- adonis2(dists[[i]] ~ Lung + Histology.NSCLC, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis8 <- adonis2(dists[[i]] ~ History.of.smoking.y.n + Batch + Lung, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)

  ad[[i]] <- rbind(test.adonis1, test.adonis2, test.adonis3, test.adonis4,
                   test.adonis5, test.adonis6, test.adonis7, test.adonis8)
}

ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "NSCLC_bothlungs.xlsx"), rowNames = T)

# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "Lung", lung_col)

ggsave(plot = NP, file.path(b.plots, "Lung_NSCLC_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Lung_NSCLC_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "Lung", lung_col)

ggsave(plot = NP, file.path(b.plots, "Lung_NSCLC_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Lung_NSCLC_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "Lung", lung_col)

ggsave(plot = NP, file.path(b.plots, "Lung_NSCLC_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Lung_NSCLC_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### Keep only diseased lung ####
physeq_beta_diseased = subset_samples(physeq_beta_0, sample_data(physeq_beta_0)$Lung == "Diseased") 

#### Three main diagnosis - diseased lung - directly isolated samples ####

physeq_beta = physeq_beta_diseased

# Select main diagnosis
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Diagnosis %in% c("NSCLC", "SCLC", "Benign"))

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
set.seed(123)
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ Diagnosis, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
  test.adonis2 <- adonis2(dists[[i]] ~ History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude) 
  test.adonis3 <- adonis2(dists[[i]] ~ History.of.smoking.y.n + Diagnosis, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis4 <- adonis2(dists[[i]] ~ Diagnosis + History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)

  ad[[i]] <- rbind(test.adonis1, test.adonis2, test.adonis3, test.adonis4)
}

ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "Three_main_Diagnosis.xlsx"), rowNames = T)

# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "Diagnosis", diagnosis_col)

ggsave(plot = NP, file.path(b.plots, "Diagnosis_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Diagnosis_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "Diagnosis", diagnosis_col)

ggsave(plot = NP, file.path(b.plots, "Diagnosis_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Diagnosis_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "Diagnosis", diagnosis_col)

ggsave(plot = NP, file.path(b.plots, "Diagnosis_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Diagnosis_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### NSCLC - diseased lung - directly isolated samples ####

physeq_beta = physeq_beta_diseased

# Select NSCLC
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Diagnosis == "NSCLC")

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
set.seed(1234)
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ Histology.NSCLC, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
  test.adonis2 <- adonis2(dists[[i]] ~ History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis3 <- adonis2(dists[[i]] ~ History.of.smoking.y.n + Histology.NSCLC, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis4 <- adonis2(dists[[i]] ~ Histology.NSCLC + History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis5 <- adonis2(dists[[i]] ~ T, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
  test.adonis6 <- adonis2(dists[[i]] ~ N, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
  test.adonis7 <- adonis2(dists[[i]] ~ M, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis8 <- adonis2(dists[[i]] ~ LibrarySizeGroup, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
  ad[[i]] <- rbind(test.adonis1, test.adonis2, test.adonis3, test.adonis4, test.adonis5, test.adonis6, test.adonis7, test.adonis8)
}

ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "NSCLC.xlsx"), rowNames = T)

## Histology
# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "Histology.NSCLC", histology_col, legend_title = "Histology")

ggsave(plot = NP, file.path(b.plots, "Histology_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "Histology.NSCLC", histology_col, legend_title = "Histology")

ggsave(plot = NP, file.path(b.plots, "Histology_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "Histology.NSCLC", histology_col, legend_title = "Histology")

ggsave(plot = NP, file.path(b.plots, "Histology_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_uwUF.svg"), width = 15, height = 8, dpi = 300)

## History.of.smoking
# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "History.of.smoking.y.n", smoker_col)

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "History.of.smoking.y.n", smoker_col)

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "History.of.smoking.y.n", smoker_col)

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_uwUF.svg"), width = 15, height = 8, dpi = 300)

## T
sample_data(physeq_beta)$T <- factor(sample_data(physeq_beta)$T)

# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "T", NULL)

ggsave(plot = NP, file.path(b.plots, "T.NSCLC_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "T.NSCLC_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "T", NULL)

ggsave(plot = NP, file.path(b.plots, "T.NSCLC_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "T.NSCLC_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "T", NULL)

ggsave(plot = NP, file.path(b.plots, "T.NSCLC_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "T.NSCLC_uwUF.svg"), width = 15, height = 8, dpi = 300)

## N
sample_data(physeq_beta)$N <- factor(sample_data(physeq_beta)$N)

# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "N", NULL)

ggsave(plot = NP, file.path(b.plots, "N.NSCLC_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "N.NSCLC_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "N", NULL)

ggsave(plot = NP, file.path(b.plots, "N.NSCLC_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "N.NSCLC_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "N", NULL)

ggsave(plot = NP, file.path(b.plots, "N.NSCLC_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "N.NSCLC_uwUF.svg"), width = 15, height = 8, dpi = 300)

# M
# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
physeq_beta.m = subset_samples(physeq_beta, !is.na(sample_data(physeq_beta)$M)) 
sample_data(physeq_beta.m)$M <- factor(sample_data(physeq_beta.m)$M)
set.seed(123)
bray_dist <- phyloseq::distance(physeq_beta.m, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta.m, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta.m, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta.m, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta.m, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta.m, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis7 <- adonis2(dists[[i]] ~ M, data = data.frame(sample_data(physeq_beta.m)), permutations = 99999, by = "margin")
  ad[[i]] <- rbind(test.adonis7)
}
ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "NSCLC.M.xlsx"), rowNames = T)

# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "M", NULL)

ggsave(plot = NP, file.path(b.plots, "M.NSCLC_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "M.NSCLC_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "M", NULL)

ggsave(plot = NP, file.path(b.plots, "M.NSCLC_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "M.NSCLC_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "M", NULL)

ggsave(plot = NP, file.path(b.plots, "M.NSCLC_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "M.NSCLC_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### NSCLC smokers - diseased lung - directly isolated samples ####

physeq_beta = physeq_beta_diseased

# Select NSCLC
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Diagnosis == "NSCLC")

# Select smokers
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$History.of.smoking.y.n == "Yes")

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
set.seed(123)
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ Histology.NSCLC, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin") %>% as.data.frame(.)
  rownames(test.adonis1)[1] <- "Histology.NSCLC"
  ad[[i]] <- rbind(test.adonis1)
}
ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "NSCLC.smokers.xlsx"), rowNames = T)

## Histology
# Plot Bray

# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "Histology.NSCLC", histology_col, legend_title = "Histology", nrow = 2)

ggsave(plot = NP, file.path(b.plots, "Histology_smokers_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_smokers_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "Histology.NSCLC", histology_col, legend_title = "Histology", nrow = 2)

ggsave(plot = NP, file.path(b.plots, "Histology_smokers_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_smokers_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "Histology.NSCLC", histology_col, legend_title = "Histology", nrow = 2)

ggsave(plot = NP, file.path(b.plots, "Histology_smokers_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_smokers_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### NSCLC - diseased lung - directly isolated samples - main histologies ####

physeq_beta = physeq_beta_diseased

# Select NSCLC
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Diagnosis == "NSCLC")

# Select histology
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"))

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ Histology.NSCLC, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
  test.adonis2 <- adonis2(dists[[i]] ~ History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis3 <- adonis2(dists[[i]] ~ History.of.smoking.y.n + Histology.NSCLC, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  test.adonis4 <- adonis2(dists[[i]] ~ Histology.NSCLC + History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  ad[[i]] <- rbind(test.adonis1, test.adonis2, test.adonis3, test.adonis4)
}
ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "NSCLC_main_h.xlsx"), rowNames = T)

## Histology
# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "Histology.NSCLC", histology_col, legend_title = "Histology")

ggsave(plot = NP, file.path(b.plots, "Histology_main_h_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_main_h_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "Histology.NSCLC", histology_col, legend_title = "Histology")

ggsave(plot = NP, file.path(b.plots, "Histology_main_h_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_main_h_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF
NP <- plot(physeq_beta, dist = "uwUF", var = "Histology.NSCLC", histology_col, legend_title = "Histology")

ggsave(plot = NP, file.path(b.plots, "Histology_main_h_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_main_h_uwUF.svg"), width = 15, height = 8, dpi = 300)

## History.of.smoking
# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "History.of.smoking.y.n", smoker_col, legend_title = "History of smoking")

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_main_h_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_main_h_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "History.of.smoking.y.n", smoker_col, legend_title = "History of smoking")

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_main_h_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_main_h_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "History.of.smoking.y.n", smoker_col, legend_title = "History of smoking")

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_main_h_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.NSCLC_main_h_uwUF.svg"), width = 15, height = 8, dpi = 300)

## T
# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "T")

ggsave(plot = NP, file.path(b.plots, "T_NSCLC_main_h_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "T_NSCLC_main_h_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "T")

ggsave(plot = NP, file.path(b.plots, "T_NSCLC_main_h_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "T_NSCLC_main_h_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF
NP <- plot(physeq_beta, dist = "uwUF", var = "T")

ggsave(plot = NP, file.path(b.plots, "T_NSCLC_main_h_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "T_NSCLC_main_h_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### NSCLC adenocarcinoma - diseased lung - directly isolated samples ####

physeq_beta = physeq_beta_diseased

# Select NSCLC
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Diagnosis == "NSCLC")

# Select histology
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Histology.NSCLC == "Adenocarcinoma")

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin", na.action = na.exclude)
  ad[[i]] <- rbind(test.adonis1)
}
ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "Adenocarcinoma.xlsx"), rowNames = T)

## History.of.smoking

# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "History.of.smoking.y.n", smoker_col, legend_title = "History of smoking")

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.adenocarcinoma_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.adenocarcinoma_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "History.of.smoking.y.n", smoker_col, legend_title = "History of smoking")

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.adenocarcinoma_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.adenocarcinoma_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "wUF", var = "History.of.smoking.y.n", smoker_col, legend_title = "History of smoking")

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.adenocarcinoma_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.adenocarcinoma_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### NSCLC smokers - diseased lung - directly isolated samples - main histologies ####

physeq_beta = physeq_beta_diseased

# Select NSCLC
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Diagnosis == "NSCLC")

# Select histology
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"))

# Select smokers
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$History.of.smoking.y.n == "Yes")

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ Histology.NSCLC, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
  ad[[i]] <- rbind(test.adonis1)
}
ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "NSCLC.smokers.main.h.xlsx"), rowNames = T)

## Histology
# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "Histology.NSCLC", histology_col, legend_title = "Histology")

ggsave(plot = NP, file.path(b.plots, "Histology_smokers_main_h_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_smokers_main_h_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "Histology.NSCLC", histology_col, legend_title = "Histology")

ggsave(plot = NP, file.path(b.plots, "Histology_smokers_main_h_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_smokers_main_h_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "Histology.NSCLC", histology_col, legend_title = "Histology")

ggsave(plot = NP, file.path(b.plots, "Histology_smokers_main_h_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "Histology_smokers_main_h_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### Benign - diseased lung - directly isolated samples ####

physeq_beta = physeq_beta_diseased

# Select NSCLC
physeq_beta = subset_samples(physeq_beta, sample_data(physeq_beta)$Diagnosis == "Benign")

# Remove not present taxa
physeq_beta = prune_taxa(taxa_sums(physeq_beta) != 0, physeq_beta)

# Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
set.seed(123)
bray_dist <- phyloseq::distance(physeq_beta, method = "bray")
wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
uunifrac_dist <- phyloseq::distance(physeq_beta, method = "uunifrac")

# Perform NMDS ordination
set.seed(123)
nmds_bray <- metaMDS(bray_dist, k = 2, trymax = 100)
nmds_wunifrac <- metaMDS(wunifrac_dist, k = 2, trymax = 100)
nmds_uunifrac <- metaMDS(uunifrac_dist, k = 2, trymax = 100)

# Perform PCoA ordination
pcoa_bray <- ordinate(physeq_beta, method = "PCoA", distance = "bray")
pcoa_wunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
pcoa_uunifrac <- ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")

# PERMANOVA
dists <- list(bray_dist=bray_dist, wunifrac_dist=wunifrac_dist, uunifrac_dist=uunifrac_dist)
ad <- list()
for(i in names(dists)){
  set.seed(123)
  test.adonis1 <- adonis2(dists[[i]] ~ History.of.smoking.y.n, data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
  ad[[i]] <- rbind(test.adonis1)
}
ad <- as.data.frame(do.call(rbind, ad))
ad
write.xlsx(ad, file.path(b.stats, "History.of.smoking.benign.xlsx"), rowNames = T)

## History.of.smoking
# Plot Bray

NP <- plot(physeq_beta, dist = "bray", var = "History.of.smoking.y.n", smoker_col, legend_title = "History of smoking")

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.benign_bray.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.benign_bray.svg"), width = 15, height = 8, dpi = 300)

# Plot wUF

NP <- plot(physeq_beta, dist = "wUF", var = "History.of.smoking.y.n", smoker_col, legend_title = "History of smoking")

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.benign_wUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.benign_wUF.svg"), width = 15, height = 8, dpi = 300)

# Plot uwUF

NP <- plot(physeq_beta, dist = "uwUF", var = "History.of.smoking.y.n", smoker_col, legend_title = "History of smoking")

ggsave(plot = NP, file.path(b.plots, "History.of.smoking.benign_uwUF.tiff"), width = 15, height = 8, dpi = 300)
ggsave(plot = NP, file.path(b.plots, "History.of.smoking.benign_uwUF.svg"), width = 15, height = 8, dpi = 300)

#### End ####

