
# 6.Taxonomy analysis #

# Bar and box plots 

#### Package setup ####

Sys.setenv(language = "EN")
library(phyloseq)
library(ggplot2)
library(gtable)
library(ggplotify)
library(openxlsx)
library(dplyr)
library(pgirmess)

theme_set(theme_bw())

#### End ####

#### Load RData with filtered and normalized phyloseq object ####

res.dir <- "Results_trimmed_cutadapt3_20241216"
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir,"5.phyloseq.filtered.RData"))

# Create directories for results
bar.dir <- file.path(res.dir, "6.Taxonomy/Bar_plots")
stats.dir <- file.path(res.dir, "6.Taxonomy/Stats")
dir.create(bar.dir, recursive = T)
dir.create(stats.dir, recursive = T)

#### End ####

#### Change Species names for visualization ####

# Absolute
tax <- as.data.frame(tax_table(physeq))   
tax[] <- lapply(tax, gsub, pattern = " ", replacement = "_")
tax$Species <- gsub("_spp.", "_spp", tax$Species)
tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis","Family"] <- paste(
  tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis","Order"], "Incertae_Sedis", sep = "_")

tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(physeq) <- tax; rm(tax)

# Relative
tax <- as.data.frame(tax_table(physeq_re))   
tax[] <- lapply(tax, gsub, pattern = " ", replacement = "_")
tax$Species <- gsub("_spp.", "_spp", tax$Species)
tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis","Family"] <- paste(
  tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis","Order"], "Incertae_Sedis", sep = "_")

tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(physeq_re) <- tax; rm(tax)

#### End ####

#### Functions ####
taxasum = function(physeq_object, taxa){
  tapply(taxa_sums(physeq_object), tax_table(physeq_object)[, taxa], sum, na.rm=TRUE) %>%
    sort(. , TRUE)
}

# Make phyloseq according to condition

phy_cond <- function(phy, condition, abundance = "relative"){
  
  # Transform all variables to factors just in case...
  df <- as.data.frame(lapply(sample_data(phy),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
  row.names(df) <- sample_names(phy)
  sample_data(phy) <- sample_data(df)
  
  # Merging samples according to condition
  ps.condition <- merge_samples(phy, condition)  
  sample_data(ps.condition)[,condition] <- rownames(sample_data(ps.condition))
  ps.condition_per <- transform_sample_counts(ps.condition, function(OTU) OTU/sum(OTU)*100) #2nd transformation to make it again in percentage
 # sample_data(ps.condition_per)[,condition] <- rownames(sample_data(ps.condition_per))
  
  if(abundance == "relative"){
  return(ps.condition_per)}else if(abundance == "absolute"){
    return(ps.condition)
  }
  
}

facets_condition <- function(phy, var, top, taxa, nrow = 1){
  sum = tapply(taxa_sums(phy), tax_table(phy)[, taxa], sum, na.rm=TRUE)
  top = names(sort(sum, TRUE))[1:top]
  top = prune_taxa((tax_table(phy)[, taxa] %in% top), phy)
  
  plot_bar(top, taxa, fill=taxa) + geom_bar(stat = "Identity", position = "stack") +
    facet_wrap(var, nrow = nrow) +
    ylab("Relative abundance") +
    theme(text = element_text(size=18),
          axis.text = element_text(size = 18),
          #axis.text.x.bottom = element_text(angle = -30),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 24)) 
  
}

plot_abundance = function(physeq, ylabn = "Relative abundance (%)",
                          Facet = "Species",
                          Color = "Diagnosis",
                          palette,
                          n = NULL,
                          legend_title = NULL){
  mphyseq = tax_glom(physeq, Facet)
  mphyseq = psmelt(mphyseq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  marks_no_sci <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = F)
  P <- ggplot(data = mphyseq,
         mapping = aes_string(x = Color, y = "Abundance",
                              color = Color, fill = Color)) + thm +
    # geom_point(size = 1, alpha = 0.1,
    #            position = position_jitter(width = 0.3)) +
    geom_boxplot(fill = NA) +
    facet_wrap(facets = Facet, nrow = n) + ylab(ylabn) +
   # stat_compare_means(method = "wilcox") +
    scale_fill_manual(values = palette) + scale_color_manual(values = palette) +
  #  scale_y_log10()
    scale_y_log10(labels = marks_no_sci)

  if (!is.null(legend_title)) {
    P <- P + labs(color = legend_title, fill = legend_title)
  } 
  print(P)
}

plot_abundance_ASV = function(physeq, ylabn = "Relative abundance (%)",
                          Facet = "Species",
                          Color = "Diagnosis",
                          palette,
                          n = NULL){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  marks_no_sci <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = F)
  ggplot(data = mphyseq,
         mapping = aes_string(x = Color, y = "Abundance",
                              color = Color, fill = Color)) + thm +
    # geom_point(size = 1, alpha = 0.1,
    #            position = position_jitter(width = 0.3)) +
    geom_boxplot(fill = NA) +
    facet_wrap(facets = Facet, nrow = n) + ylab(ylabn) +
    # stat_compare_means(method = "wilcox") +
    scale_fill_manual(values = palette) + scale_color_manual(values = palette) +
    #  scale_y_log10()
    scale_y_log10(labels = marks_no_sci)
}
# ps.phylum <- tax_glom(physeq_tax_re, "Phylum")
# taxa_names(ps.phylum) <- tax_table(ps.phylum)[,"Phylum"]
# ps.phylum_diseased <- subset_samples(ps.phylum, sample_data(ps.phylum)$Lung == "Diseased")
# median(as.numeric(otu_table(ps.phylum_diseased)["Actinomycetota",]))

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

thm.x.2 <- theme(text = element_text(size = text_size),
                   axis.text.x = element_text(angle = 0, hjust = 0),
                   legend.position = "bottom")
#### End ####

#### Select phyloseq object #### 

physeq_tax <- physeq
physeq_tax_re <- physeq_re

# Remove controls
physeq_tax <- subset_samples(physeq_tax, sample_names(physeq_tax) %in% c(rownames(sample_info[!sample_info$Sample_or_Control %in% c("Negative control sample", "Positive control sample"),]))) 

# Remove samples with no reads
physeq_tax <- subset_samples(physeq_tax, colSums(otu_table(physeq_tax)) > 0) 

# Select main diagnosis
physeq_tax = subset_samples(physeq_tax, sample_data(physeq_tax)$Diagnosis %in% c("NSCLC", "SCLC", "Benign"))

# Keep only directly isolated samples
physeq_tax <- subset_samples(physeq_tax, sample_data(physeq_tax)$Isolation == "Direct_Isolation") 

# Keep only samples read number > n
# n = 1000
# length(sample_sums(physeq_tax)[sample_sums(physeq_tax) >=n])
# less_n <- names(sample_sums(physeq_tax)[sample_sums(physeq_tax) <n])
# physeq_tax = subset_samples(physeq_tax, !sample_names(physeq_tax) %in% less_n) 
# x = length(sample_names(physeq_tax))

# Remove not present taxa
physeq_tax <- prune_taxa(taxa_sums(physeq_tax) != 0, physeq_tax)

# Keep same samples and taxa for relative abundance
physeq_tax_re <- prune_samples(sample_names(physeq_tax_re) %in% sample_names(physeq_tax), physeq_tax_re)
physeq_tax_re <- prune_taxa(taxa_names(physeq_tax_re) %in% taxa_names(physeq_tax), physeq_tax_re)
otu_table(physeq_tax_re) <- otu_table(physeq_tax_re)*100

physeq_tax
physeq_tax_re

#### End ####

#### Save matrix for correlations ####

# Agglomerate taxa and merge with sample data
matrix_absolute <- list()
for (taxon in c("Phylum", "Class", "Order", "Family", "Genus")) {
  physeq_tax_taxon <- tax_glom(physeq_tax, taxon)
  taxa_names(physeq_tax_taxon) <- tax_table(physeq_tax_taxon)[,taxon]
  matrix_absolute[[taxon]] <- merge(sample_data(physeq_tax), t(otu_table(physeq_tax_taxon)), by = 0)
}
matrix_absolute[["ASV"]] <- merge(sample_data(physeq_tax), t(otu_table(physeq_tax)), by = 0)

matrix_relative <- list()
for (taxon in c("Phylum", "Class", "Order", "Family", "Genus")) {
  physeq_tax_re_taxon <- tax_glom(physeq_tax_re, taxon)
  taxa_names(physeq_tax_re_taxon) <- tax_table(physeq_tax_re_taxon)[,taxon]
  matrix_relative[[taxon]] <- merge(sample_data(physeq_tax_re), t(otu_table(physeq_tax_re_taxon)), by = 0)
  }
matrix_relative[["ASV"]] <- merge(sample_data(physeq_tax_re), t(otu_table(physeq_tax_re)), by = 0)

# Save matrices
write.xlsx(matrix_absolute, file.path(stats.dir, "correlation_matrix_absolute.xlsx"))
write.xlsx(matrix_relative, file.path(stats.dir, "correlation_matrix_relative.xlsx"))

#### End ####

#### Q1a: Is there a difference between Diseased vs the parallel Normal Lung? ####

# Create phyloseq object for diseased and normal lung
physeq_A <- subset_samples(physeq_tax, sample_data(physeq_tax)$Lung == "Diseased")
physeq_B <- subset_samples(physeq_tax, sample_data(physeq_tax)$Lung == "Normal")
physeq_lung <- phy_cond(physeq_tax, "Lung", "absolute")
physeq_lung_re <- phy_cond(physeq_tax_re, "Lung")

# Define top n ASVs
n = 20
top <- names(sort(taxa_sums(physeq_tax), decreasing=TRUE))[1:n]  # adjust number to wished top ones
top_diseased <- names(sort(taxa_sums(physeq_A), decreasing=TRUE))[1:n] 
top_normal <- names(sort(taxa_sums(physeq_B), decreasing=TRUE))[1:n] 
topasv <- prune_taxa(top, physeq_lung)

# Absolute abundance
plot_bar(topasv, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Phylum"), fill="Phylum")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

# Top n Phyla
n = 10

top = names(taxasum(physeq_tax, "Phylum"))[1:n]
top_diseased = names(taxasum(physeq_A, "Phylum"))[1:n]
top_normal = names(taxasum(physeq_B, "Phylum"))[1:n]

topphy = prune_taxa((tax_table(physeq_lung)[, "Phylum"] %in% top), physeq_lung)
topphyrel = prune_taxa((tax_table(physeq_lung_re)[, "Phylum"] %in% top), physeq_lung_re)

# Absolute abundance
plot_bar(topphy, y =  "Abundance", title = paste("Absolute abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q1a_absolute_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q1a_absolute_phyla.svg"), width = 10, height = 8, dpi = 300)

#Relative Abundance
plot_bar(topphyrel, y =  "Abundance", title = paste("Relative abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q1a_relative_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q1a_relative_phyla.svg"), width = 10, height = 8, dpi = 300)

facets_condition(physeq_lung_re, "Lung", n, "Phylum")
ggsave(file.path(bar.dir, "P.Q1a_relative2_phyla.tiff"), width = 10, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "P.Q1a_relative2_phyla.svg"), width = 10, height = 7, dpi = 300)

toprel = prune_taxa((tax_table(physeq_tax_re)[, "Phylum"] %in% top), physeq_tax_re)
plot_abundance(toprel, Facet = "Phylum", Color = "Lung", palette = lung_col, n = 2) 

ggsave(file.path(bar.dir, "P.Q1a_relative_box_phyla.tiff"), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, "P.Q1a_relative_box_phyla.svg"), width = 16, height = 10, dpi = 300)

# Top n Genera
n = 10

top = names(taxasum(physeq_tax, "Genus"))[1:n]
top_diseased = names(taxasum(physeq_A, "Genus"))[1:n]
top_normal = names(taxasum(physeq_B, "Genus"))[1:n]

topgen = prune_taxa((tax_table(physeq_lung)[, "Genus"] %in% top), physeq_lung)
topgenrel = prune_taxa((tax_table(physeq_lung_re)[, "Genus"] %in% top), physeq_lung_re)

# Absolute abundance
plot_bar(topgen, y =  "Abundance", title = paste("Absolute abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q1a_absolute_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_absolute_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

#Relative Abundance
plot_bar(topgenrel, y =  "Abundance", title = paste("Relative abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q1a_relative_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_relative_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

facets_condition(physeq_lung_re, "Lung", n, "Genus")
ggsave(file.path(bar.dir, paste("P.Q1a_relative2_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_relative2_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

toprel = prune_taxa((tax_table(physeq_tax_re)[, "Genus"] %in% top), physeq_tax_re)
plot_abundance(toprel, Facet = "Genus", Color = "Lung", palette = lung_col, n = 2) 

ggsave(file.path(bar.dir, paste("P.Q1a_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)

#Significant taxa
mat <- read.delim(file.path(stats.dir,"Rhea_paired_two/no_median_cutoff/Main_diagnosis",
                            "Rhea_Lung_Phylum_correlation_matrix_relative_2025-01-15/Lung_Phylum_correlation_matrix_relative-sign_pairs.tab"),
                  stringsAsFactors=TRUE)

mat = mat[mat$corrected < 0.1,]
mat

taxon = "Family"
mat <- read.delim(file.path(stats.dir,"Rhea_paired_two/no_median_cutoff/Main_diagnosis", paste0("Rhea_Lung_", taxon,
                                                                                               "_correlation_matrix_relative_2025-01-15/Lung_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_tax_re)[, taxon] %in% sig), physeq_tax_re)
plot_abundance(topsig, Facet = taxon, Color = "Lung", palette = lung_col, n = 1) 
ggsave(file.path(bar.dir, paste("P.Q1a_relative_sig_", taxon, ".tiff")), width = 12, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_relative_sig_", taxon, ".tiff")), width = 12, height = 8, dpi = 300)

taxon = "Genus"
mat <- read.delim(file.path(stats.dir,"Rhea_paired_two/no_median_cutoff/Main_diagnosis", paste0("Rhea_Lung_", taxon,
                          "_correlation_matrix_relative_2025-01-15/Lung_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

sig = mat[mat$corrected < 0.1,"measure"]
topsig = prune_taxa((tax_table(physeq_tax_re)[, taxon] %in% sig), physeq_tax_re)
plot_abundance(topsig, Facet = taxon, Color = "Lung", palette = lung_col, n = 1) 
ggsave(file.path(bar.dir, paste("P.Q1a_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)

taxon = "ASV"
mat <- read.delim(file.path(stats.dir,"Rhea_paired_two/no_median_cutoff/Main_diagnosis", paste0("Rhea_Lung_", taxon,
                                                                                               "_correlation_matrix_relative_2025-01-15/Lung_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

sig = mat[mat$corrected < 0.1,"measure"]
topsig = prune_taxa((rownames(tax_table(physeq_tax_re)) %in% sig), physeq_tax_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "Lung", palette = lung_col, n = 1) 
ggsave(file.path(bar.dir, paste("P.Q1a_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)

#### End ####

#### Q1b: Is there a difference between Diseased vs the parallel Normal Lung in NSCLC? ####

physeq_q <- physeq_tax
physeq_q_re <- physeq_tax_re

# Select NSCLC
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis %in% c("NSCLC"))

# Remove not present taxa
physeq_q <- prune_taxa(taxa_sums(physeq_q) != 0, physeq_q)

# Keep same samples and taxa for relative abundance
physeq_q_re <- prune_samples(sample_names(physeq_q_re) %in% sample_names(physeq_q), physeq_q_re)
physeq_q_re <- prune_taxa(taxa_names(physeq_q_re) %in% taxa_names(physeq_q), physeq_q_re)

# Create phyloseq object for diseased and normal lung
physeq_A <- subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Diseased")
physeq_B <- subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Normal")
physeq_lung <- phy_cond(physeq_q, "Lung", "absolute")
physeq_lung_re <- phy_cond(physeq_q_re, "Lung")

# Define top n ASVs
n = 20
top <- names(sort(taxa_sums(physeq_q), decreasing=TRUE))[1:n]  # adjust number to wished top ones
top_diseased <- names(sort(taxa_sums(physeq_A), decreasing=TRUE))[1:n] 
top_normal <- names(sort(taxa_sums(physeq_B), decreasing=TRUE))[1:n] 
topasv <- prune_taxa(top, physeq_lung)

# Absolute abundance
plot_bar(topasv, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Phylum"), fill="Phylum")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

# Top n Phyla
n = 10

top = names(taxasum(physeq_q, "Phylum"))[1:n]
top_diseased = names(taxasum(physeq_A, "Phylum"))[1:n]
top_normal = names(taxasum(physeq_B, "Phylum"))[1:n]

topphy = prune_taxa((tax_table(physeq_lung)[, "Phylum"] %in% top), physeq_lung)
topphyrel = prune_taxa((tax_table(physeq_lung_re)[, "Phylum"] %in% top), physeq_lung_re)

# Absolute abundance
plot_bar(topphy, y =  "Abundance", title = paste("Absolute abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q1a_absolute_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q1a_absolute_phyla.svg"), width = 10, height = 8, dpi = 300)

#Relative Abundance
plot_bar(topphyrel, y =  "Abundance", title = paste("Relative abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q1a_relative_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q1a_relative_phyla.svg"), width = 10, height = 8, dpi = 300)

facets_condition(physeq_lung_re, "Lung", n, "Phylum")
ggsave(file.path(bar.dir, "P.Q1a_relative2_phyla.tiff"), width = 10, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "P.Q1a_relative2_phyla.svg"), width = 10, height = 7, dpi = 300)

# Top n Genera
n = 20

top = names(taxasum(physeq_q, "Genus"))[1:n]
top_diseased = names(taxasum(physeq_A, "Genus"))[1:n]
top_normal = names(taxasum(physeq_B, "Genus"))[1:n]

topgen = prune_taxa((tax_table(physeq_lung)[, "Genus"] %in% top), physeq_lung)
topgenrel = prune_taxa((tax_table(physeq_lung_re)[, "Genus"] %in% top), physeq_lung_re)

# Absolute abundance
plot_bar(topgen, y =  "Abundance", title = paste("Absolute abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q1a_absolute_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_absolute_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

#Relative Abundance
plot_bar(topgenrel, y =  "Abundance", title = paste("Relative abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q1a_relative_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_relative_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

facets_condition(physeq_lung_re, "Lung", n, "Genus")
ggsave(file.path(bar.dir, paste("P.Q1a_relative2_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1a_relative2_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

#Significant taxa
mat <- read.delim(file.path(stats.dir,"Rhea_paired_two/median_cutoff_0/NSCLC",
                            "Rhea_Lung_Phylum_correlation_matrix_relative_2025-01-15/Lung_Phylum_correlation_matrix_relative-sign_pairs.tab"),
                  stringsAsFactors=TRUE)

View(mat)
mat = mat[mat$corrected < 0.1,]

taxon = "Family"
mat <- read.delim(file.path(stats.dir,"Rhea_paired_two/median_cutoff_0/NSCLC", paste0("Rhea_Lung_", taxon,
                                                                                               "_correlation_matrix_relative_2025-01-15/Lung_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_tax_re)
plot_abundance(topsig, Facet = taxon, Color = "Lung", palette = lung_col, n = 1) 
ggsave(file.path(bar.dir, paste("P.Q1b_relative_sig_", taxon, ".tiff")), width = 12, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1b_relative_sig_", taxon, ".tiff")), width = 12, height = 8, dpi = 300)

taxon = "Genus"
mat <- read.delim(file.path(stats.dir,"Rhea_paired_two/median_cutoff_0/NSCLC", paste0("Rhea_Lung_", taxon,
                                                                                               "_correlation_matrix_relative_2025-01-15/Lung_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_tax_re)
plot_abundance(topsig, Facet = taxon, Color = "Lung", palette = lung_col, n = 1) 
ggsave(file.path(bar.dir, paste("P.Q1b_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1b_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)

taxon = "ASV"
mat <- read.delim(file.path(stats.dir,"Rhea_paired_two/median_cutoff_0/NSCLC", paste0("Rhea_Lung_", taxon,
                                                                                               "_correlation_matrix_relative_2025-01-15/Lung_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((rownames(tax_table(physeq_q_re)) %in% sig), physeq_q_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "Lung", palette = lung_col, n = 1) 
ggsave(file.path(bar.dir, paste("P.Q1b_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q1b_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)

#### End ####

#### Q2a: Is there a difference between Diseased vs the parallel Normal Lung without samples with synchronous tumor (after decontamination, only directly isolated samples) ####

physeq_q <- physeq_tax
physeq_q_re <- physeq_tax_re

# Remove non synchronous tumor
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Synchronous.tumor == "No")
physeq_q_re = subset_samples(physeq_q_re, sample_data(physeq_q_re)$Synchronous.tumor == "No")

# Remove not present taxa
physeq_q <- prune_taxa(taxa_sums(physeq_q) != 0, physeq_q)
physeq_q_re <- prune_taxa(taxa_names(physeq_q_re) %in% taxa_names(physeq_q), physeq_q_re)

# Create phyloseq object for diseased and normal lung
physeq_A <- subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Diseased")
physeq_B <- subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Normal")
physeq_lung <- phy_cond(physeq_q, "Lung", "absolute")
physeq_lung_re <- phy_cond(physeq_q_re, "Lung")

# Define top n ASVs
n = 100
top <- names(sort(taxa_sums(physeq_q), decreasing=TRUE))[1:n]  # adjust number to wished top ones
top_diseased <- names(sort(taxa_sums(physeq_A), decreasing=TRUE))[1:n] 
top_normal <- names(sort(taxa_sums(physeq_B), decreasing=TRUE))[1:n] 
topasv <- prune_taxa(top, physeq_lung)

# Absolute abundance
plot_bar(topasv, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Phylum"), fill="Phylum")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

# Top n Phyla
n = 10

top = names(taxasum(physeq_q, "Phylum"))[1:n]
top_diseased = names(taxasum(physeq_A, "Phylum"))[1:n]
top_normal = names(taxasum(physeq_B, "Phylum"))[1:n]

topphy = prune_taxa((tax_table(physeq_lung)[, "Phylum"] %in% top), physeq_lung)
topphyrel = prune_taxa((tax_table(physeq_lung_re)[, "Phylum"] %in% top), physeq_lung_re)

# Absolute abundance
plot_bar(topphy, y =  "Abundance", title = paste("Absolute abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q2a_absolute_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q2a_absolute_phyla.svg"), width = 10, height = 8, dpi = 300)

#Relative Abundance
plot_bar(topphyrel, y =  "Abundance", title = paste("Relative abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q2a_relative_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q2a_relative_phyla.svg"), width = 10, height = 8, dpi = 300)

facets_condition(physeq_lung_re, "Lung", n, "Phylum")
ggsave(file.path(bar.dir, "P.Q2a_relative2_phyla.tiff"), width = 10, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "P.Q2a_relative2_phyla.svg"), width = 10, height = 7, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Phylum"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Phylum", Color = "Lung", palette = lung_col, n = 2) 

ggsave(file.path(bar.dir, "P.Q2a_relative_box_phyla.tiff"), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, "P.Q2a_relative_box_phyla.svg"), width = 16, height = 10, dpi = 300)

# Top n Genera
n = 10

top = names(taxasum(physeq_q, "Genus"))[1:n]
top_diseased = names(taxasum(physeq_A, "Genus"))[1:n]
top_normal = names(taxasum(physeq_B, "Genus"))[1:n]

topgen = prune_taxa((tax_table(physeq_lung)[, "Genus"] %in% top), physeq_lung)
topgenrel = prune_taxa((tax_table(physeq_lung_re)[, "Genus"] %in% top), physeq_lung_re)

# Absolute abundance
plot_bar(topgen, y =  "Abundance", title = paste("Absolute abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q2a_absolute_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q2a_absolute_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

#Relative Abundance
plot_bar(topgenrel, y =  "Abundance", title = paste("Relative abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Lung, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q2a_relative_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q2a_relative_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

facets_condition(physeq_lung_re, "Lung", n, "Genus")
ggsave(file.path(bar.dir, paste("P.Q2a_relative2_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q2a_relative2_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Genus"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Genus", Color = "Lung", palette = lung_col, n = 2) 

ggsave(file.path(bar.dir, paste("P.Q2a_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q2a_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)

#### End ####

#### Q3: What are the difference between NSCLC vs SCLC vs Bening Tumor? ####
physeq_q <- physeq_tax
physeq_q_re <- physeq_tax_re

# Keep only diseased lung
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Diseased") 
#physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Diseased") 

# Select main diagnoses
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis %in% c("NSCLC", "SCLC", "Benign"))
#physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis %in% c("NSCLC", "SCLC", "Benign"))

# Remove not present taxa
physeq_q <- prune_taxa(taxa_sums(physeq_q) != 0, physeq_q)

# Keep same samples and taxa for relative abundance
physeq_q_re <- prune_samples(sample_names(physeq_q_re) %in% sample_names(physeq_q), physeq_q_re)
physeq_q_re <- prune_taxa(taxa_names(physeq_q_re) %in% taxa_names(physeq_q), physeq_q_re)

# Create phyloseq object for diseased and normal lung
physeq_NSCLC <- subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis == "NSCLC")
physeq_SCLC <- subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis == "SCLC")
physeq_benign <- subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis == "Benign")
physeq_diagnosis <- phy_cond(physeq_q, "Diagnosis", "absolute")
physeq_diagnosis_re <- phy_cond(physeq_q_re, "Diagnosis")

# Define top n ASVs
n = 100
top <- names(sort(taxa_sums(physeq_q), decreasing=TRUE))[1:n]  # adjust number to wished top ones
top_diseased <- names(sort(taxa_sums(physeq_A), decreasing=TRUE))[1:n] 
top_normal <- names(sort(taxa_sums(physeq_B), decreasing=TRUE))[1:n] 
topasv <- prune_taxa(top, physeq_diagnosis)

# Absolute abundance
plot_bar(topasv, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Phylum"), fill="Phylum")+ facet_wrap(~Diagnosis, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

# Top n Phyla
n = 10

top = names(taxasum(physeq_q, "Phylum"))[1:n]
top_NSCLC = names(taxasum(physeq_NSCLC, "Phylum"))[1:n]
top_SCLC = names(taxasum(physeq_SCLC, "Phylum"))[1:n]
top_benign = names(taxasum(physeq_benign, "Phylum"))[1:n]
print(cbind(top, top_NSCLC, top_SCLC, top_benign))

topphy = prune_taxa((tax_table(physeq_diagnosis)[, "Phylum"] %in% top), physeq_diagnosis)
topphyrel = prune_taxa((tax_table(physeq_diagnosis_re)[, "Phylum"] %in% top), physeq_diagnosis_re)

# Absolute abundance
plot_bar(topphy, y =  "Abundance", title = paste("Absolute abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Diagnosis, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q3_absolute_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q3_absolute_phyla.svg"), width = 10, height = 8, dpi = 300)

#Relative Abundance
plot_bar(topphyrel, y =  "Abundance", title = paste("Relative abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Diagnosis, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q3_relative_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q3_relative_phyla.svg"), width = 10, height = 8, dpi = 300)

facets_condition(physeq_diagnosis_re, "Diagnosis", n, "Phylum")
ggsave(file.path(bar.dir, "P.Q3_relative2_phyla.tiff"), width = 10, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "P.Q3_relative2_phyla.svg"), width = 10, height = 7, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Phylum"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Phylum", Color = "Diagnosis", palette = diagnosis_col, n = 2) 

ggsave(file.path(bar.dir, "P.Q3_relative_box_phyla.tiff"), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, "P.Q3_relative_box_phyla.svg"), width = 16, height = 10, dpi = 300)

# Top n Genera
n = 10

top = names(taxasum(physeq_q, "Genus"))[1:n]
top_NSCLC = names(taxasum(physeq_NSCLC, "Genus"))[1:n]
top_SCLC = names(taxasum(physeq_SCLC, "Genus"))[1:n]
top_benign = names(taxasum(physeq_benign, "Genus"))[1:n]
print(cbind(top, top_NSCLC, top_SCLC, top_benign))

topgen = prune_taxa((tax_table(physeq_diagnosis)[, "Genus"] %in% top), physeq_diagnosis)
topgenrel = prune_taxa((tax_table(physeq_diagnosis_re)[, "Genus"] %in% top), physeq_diagnosis_re)

# Absolute abundance
plot_bar(topgen, y =  "Abundance", title = paste("Absolute abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Diagnosis, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q3_absolute_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_absolute_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

#Relative Abundance
plot_bar(topgenrel, y =  "Abundance", title = paste("Relative abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Diagnosis, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q3_relative_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_relative_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

facets_condition(physeq_diagnosis_re, "Diagnosis", n, "Genus")
ggsave(file.path(bar.dir, paste("P.Q3_relative2_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_relative2_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Genus"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Genus", Color = "Diagnosis", palette = diagnosis_col, n = 2) 

ggsave(file.path(bar.dir, paste("P.Q3_relative_box_genus_", n, ".svg")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)

#Significant taxa
diagnosis_col["SCLC"] <- "#344CB7"
taxon = "Phylum"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Diagnosis_", taxon,
                                                                "_correlation_matrix_relative_2025-01-29/Diagnosis_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat = mat[mat$corrected < 0.1,]
#mat = mat[mat$corrected < 0.1 & mat$Group1 == " -",]
sig = unique(mat[mat$corrected < 0.1,"measure"])
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "Diagnosis", palette = diagnosis_col, n = 1) 
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".svg")), width = 12, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".tiff")), width = 12, height = 8, dpi = 300)

taxon = "Family"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Diagnosis_", taxon,
                                                                "_correlation_matrix_relative_2025-01-29/Diagnosis_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

#mat = mat[mat$corrected < 0.1 & mat$Group1 == " -",]
sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "Diagnosis", palette = diagnosis_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".svg")), width = 17, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".tiff")), width = 17, height = 10, dpi = 300)

taxon = "Genus"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Diagnosis_", taxon,
                                                                "_correlation_matrix_relative_2025-01-29/Diagnosis_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

#mat = mat[mat$corrected < 0.1 & mat$Group1 == " -",]
sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "Diagnosis", palette = diagnosis_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".svg")), width = 16, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".tiff")), width = 16, height = 8, dpi = 300)

taxon = "ASV"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Diagnosis_", taxon,
                                                                "_correlation_matrix_relative_2025-01-29/Diagnosis_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

sig = mat[mat$corrected < 0.05 & !is.na(mat$corrected) & mat$Group1 == " -","measure"]
sig
topsig = prune_taxa((rownames(tax_table(physeq_q_re)) %in% sig), physeq_q_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "Diagnosis", palette = diagnosis_col, n = 3) 
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".svg")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".tiff")), width = 16, height = 10, dpi = 300)

#### End ####

#### Q4: What are the difference between NSCLC main histologies (Adeno vs Squamous)? ####
physeq_q <- physeq_tax
physeq_q_re <- physeq_tax_re

# Keep only diseased lung
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Diseased") 

# Select NSCLC
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis == "NSCLC")

# Select histology
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"))

# Remove not present taxa
physeq_q <- prune_taxa(taxa_sums(physeq_q) != 0, physeq_q)

# Keep same samples and taxa for relative abundance
physeq_q_re <- prune_samples(sample_names(physeq_q_re) %in% sample_names(physeq_q), physeq_q_re)
physeq_q_re <- prune_taxa(taxa_names(physeq_q_re) %in% taxa_names(physeq_q), physeq_q_re)

# Create phyloseq object for diseased and normal lung
physeq_adeno <- subset_samples(physeq_q, sample_data(physeq_q)$Histology.NSCLC == "Adenocarcinoma")
physeq_squamous <- subset_samples(physeq_q, sample_data(physeq_q)$Histology.NSCLC == "Squamous cell carcinoma")
physeq_histology <- phy_cond(physeq_q, "Histology.NSCLC", "absolute")
physeq_histology_re <- phy_cond(physeq_q_re, "Histology.NSCLC")

# Define top n ASVs
n = 100
top <- names(sort(taxa_sums(physeq_q), decreasing=TRUE))[1:n]  # adjust number to wished top ones
top_diseased <- names(sort(taxa_sums(physeq_A), decreasing=TRUE))[1:n] 
top_normal <- names(sort(taxa_sums(physeq_B), decreasing=TRUE))[1:n] 
topasv <- prune_taxa(top, physeq_histology)

# Absolute abundance
plot_bar(topasv, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Phylum"), fill="Phylum")+ facet_wrap(~Histology.NSCLC, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

# Top n Phyla
n = 10

top = names(taxasum(physeq_q, "Phylum"))[1:n]
top_adeno = names(taxasum(physeq_adeno, "Phylum"))[1:n]
top_squamous = names(taxasum(physeq_squamous, "Phylum"))[1:n]
print(cbind(top, top_adeno, top_squamous))

topphy = prune_taxa((tax_table(physeq_histology)[, "Phylum"] %in% top), physeq_histology)
topphyrel = prune_taxa((tax_table(physeq_histology_re)[, "Phylum"] %in% top), physeq_histology_re)

# Absolute abundance
plot_bar(topphy, y =  "Abundance", title = paste("Absolute abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Histology.NSCLC, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q4_absolute_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q4_absolute_phyla.svg"), width = 10, height = 8, dpi = 300)

#Relative Abundance
plot_bar(topphyrel, y =  "Abundance", title = paste("Relative abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~Histology.NSCLC, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q4_relative_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q4_relative_phyla.svg"), width = 10, height = 8, dpi = 300)

facets_condition(physeq_histology_re, "Histology.NSCLC", n, "Phylum")
ggsave(file.path(bar.dir, "P.Q4_relative2_phyla.tiff"), width = 10, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "P.Q4_relative2_phyla.svg"), width = 10, height = 7, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Phylum"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Phylum", Color = "Histology.NSCLC", palette = histology_col, n = 2) 

ggsave(file.path(bar.dir, "P.Q4_relative_box_phyla.tiff"), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, "P.Q4_relative_box_phyla.svg"), width = 16, height = 10, dpi = 300)

# Top n Genera
n = 10

top = names(taxasum(physeq_q, "Genus"))[1:n]
top_adeno = names(taxasum(physeq_adeno, "Genus"))[1:n]
top_squamous = names(taxasum(physeq_squamous, "Genus"))[1:n]
print(cbind(top, top_adeno, top_squamous))

topgen = prune_taxa((tax_table(physeq_histology)[, "Genus"] %in% top), physeq_histology)
topgenrel = prune_taxa((tax_table(physeq_histology_re)[, "Genus"] %in% top), physeq_histology_re)

# Absolute abundance
plot_bar(topgen, y =  "Abundance", title = paste("Absolute abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Histology.NSCLC, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q4_absolute_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q4_absolute_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

#Relative Abundance
plot_bar(topgenrel, y =  "Abundance", title = paste("Relative abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~Histology.NSCLC, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q4_relative_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q4_relative_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

facets_condition(physeq_histology_re, "Histology.NSCLC", n, "Genus")
ggsave(file.path(bar.dir, paste("P.Q4_relative2_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q4_relative2_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Genus"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Genus", Color = "Histology.NSCLC", palette = histology_col, n = 2) 

ggsave(file.path(bar.dir, paste("P.Q4_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q4_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)

#Significant taxa
taxon = "Phylum"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Histology.NSCLC_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/Histology.NSCLC_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = unique(mat[mat$corrected < 0.1,"measure"])
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "Histology", palette = diagnosis_col, n = 1) 
ggsave(file.path(bar.dir, paste("P.Q4_relative_sig_", taxon, ".svg")), width = 12, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q4_relative_sig_", taxon, ".tiff")), width = 12, height = 8, dpi = 300)

taxon = "Family"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Histology.NSCLC_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/Histology.NSCLC_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "Diagnosis", palette = diagnosis_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".svg")), width = 17, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".tiff")), width = 17, height = 10, dpi = 300)

taxon = "Genus"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Histology.NSCLC_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/Histology.NSCLC_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "Diagnosis", palette = diagnosis_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".svg")), width = 16, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q3_relative_sig_", taxon, ".tiff")), width = 16, height = 8, dpi = 300)

taxon = "ASV"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Histology.NSCLC_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/Histology.NSCLC_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = mat[mat$corrected < 0.05 & !is.na(mat$corrected),"measure"]
sig
topsig = prune_taxa((rownames(tax_table(physeq_q_re)) %in% sig), physeq_q_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "Diagnosis", palette = diagnosis_col, n = 3) 
ggsave(file.path(bar.dir, paste("P.Q4_relative_sig_", taxon, ".svg")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q4_relative_sig_", taxon, ".tiff")), width = 16, height = 10, dpi = 300)

#### End ####

#### Q5b: Does history of smoking have an impact on the lung microbiome? ####
physeq_q <- physeq_tax
physeq_q_re <- physeq_tax_re

# Keep only diseased lung
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Diseased") 

# Remove NA
physeq_q = subset_samples(physeq_q, !is.na(sample_data(physeq_q)$History.of.smoking.y.n))

# Select NSCLC
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis == "NSCLC")

# Remove not present taxa
physeq_q <- prune_taxa(taxa_sums(physeq_q) != 0, physeq_q)

# Keep same samples and taxa for relative abundance
physeq_q_re <- prune_samples(sample_names(physeq_q_re) %in% sample_names(physeq_q), physeq_q_re)
physeq_q_re <- prune_taxa(taxa_names(physeq_q_re) %in% taxa_names(physeq_q), physeq_q_re)

# Create phyloseq object for diseased and normal lung
physeq_smoker <- subset_samples(physeq_q, sample_data(physeq_q)$History.of.smoking.y.n == "Yes")
physeq_nonsmoker <- subset_samples(physeq_q, sample_data(physeq_q)$History.of.smoking.y.n == "No")
physeq_smoking <- phy_cond(physeq_q, "History.of.smoking.y.n", "absolute")
physeq_smoking_re <- phy_cond(physeq_q_re, "History.of.smoking.y.n")

# Define top n ASVs
n = 100
top <- names(sort(taxa_sums(physeq_q), decreasing=TRUE))[1:n]  # adjust number to wished top ones
top_diseased <- names(sort(taxa_sums(physeq_A), decreasing=TRUE))[1:n] 
top_normal <- names(sort(taxa_sums(physeq_B), decreasing=TRUE))[1:n] 
topasv <- prune_taxa(top, physeq_smoking)

# Absolute abundance
plot_bar(topasv, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Phylum"), fill="Phylum")+ facet_wrap(~History.of.smoking.y.n, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

# Top n Phyla
n = 10

top = names(taxasum(physeq_q, "Phylum"))[1:n]
top_smoker = names(taxasum(physeq_smoker, "Phylum"))[1:n]
top_nonsmoker = names(taxasum(physeq_nonsmoker, "Phylum"))[1:n]
print(cbind(top, top_smoker, top_nonsmoker))

topphy = prune_taxa((tax_table(physeq_smoking)[, "Phylum"] %in% top), physeq_smoking)
topphyrel = prune_taxa((tax_table(physeq_smoking_re)[, "Phylum"] %in% top), physeq_smoking_re)

# Absolute abundance
plot_bar(topphy, y =  "Abundance", title = paste("Absolute abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~History.of.smoking.y.n, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q5b_absolute_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q5b_absolute_phyla.svg"), width = 10, height = 8, dpi = 300)

#Relative Abundance
plot_bar(topphyrel, y =  "Abundance", title = paste("Relative abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~History.of.smoking.y.n, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q5b_relative_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q5b_relative_phyla.svg"), width = 10, height = 8, dpi = 300)

facets_condition(physeq_smoking_re, "History.of.smoking.y.n", n, "Phylum")
ggsave(file.path(bar.dir, "P.Q5b_relative2_phyla.tiff"), width = 10, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "P.Q5b_relative2_phyla.svg"), width = 10, height = 7, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Phylum"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Phylum", Color = "History.of.smoking.y.n", palette = smoker_col, n = 2) 

ggsave(file.path(bar.dir, "P.Q5b_relative_box_phyla.tiff"), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, "P.Q5b_relative_box_phyla.svg"), width = 16, height = 10, dpi = 300)

# Top n Genera
n = 10

top = names(taxasum(physeq_q, "Genus"))[1:n]
top_smoker = names(taxasum(physeq_smoker, "Genus"))[1:n]
top_nonsmoker = names(taxasum(physeq_nonsmoker, "Genus"))[1:n]
print(cbind(top, top_smoker, top_nonsmoker))

topgen = prune_taxa((tax_table(physeq_smoking)[, "Genus"] %in% top), physeq_smoking)
topgenrel = prune_taxa((tax_table(physeq_smoking_re)[, "Genus"] %in% top), physeq_smoking_re)

# Absolute abundance
plot_bar(topgen, y =  "Abundance", title = paste("Absolute abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~History.of.smoking.y.n, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q5b_absolute_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5b_absolute_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

#Relative Abundance
plot_bar(topgenrel, y =  "Abundance", title = paste("Relative abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~History.of.smoking.y.n, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q5b_relative_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5b_relative_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

facets_condition(physeq_smoking_re, "History.of.smoking.y.n", n, "Genus")
ggsave(file.path(bar.dir, paste("P.Q5b_relative2_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5b_relative2_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Genus"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Genus", Color = "History.of.smoking.y.n", palette = smoker_col, n = 2) 

ggsave(file.path(bar.dir, paste("P.Q5b_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5b_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)

#Significant taxa
taxon = "Phylum"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_History.of.smoking.y.n_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/History.of.smoking.y.n_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

mat <- read.delim("C:/Users/posadas/Desktop/GitHub/LML/RESULT~2/63808~1.TAX/Stats/RH8B1F~1/NO_MED~1/RHEA_H~1.N_P/HISTOR~2.TAB", stringsAsFactors=TRUE)

mat = mat[mat$corrected < 0.1,]
#mat = mat[mat$corrected < 0.1 & mat$Group1 == " -",]
sig = unique(mat[mat$corrected < 0.1,"measure"])
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "History.of.smoking.y.n", palette = smoker_col, n = 1) 
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_", taxon, ".svg")), width = 8, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)

taxon = "Family"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_History.of.smoking.y.n_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/History.of.smoking.y.n_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat <- read.delim("C:/Users/posadas/Desktop/GitHub/LML/RESULT~2/63808~1.TAX/Stats/RH8B1F~1/NO_MED~1/RHEA_H~1.N_F/HISTOR~2.TAB", stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = mat[mat$corrected < 0.05,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "History.of.smoking.y.n", palette = smoker_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_0.05", taxon, ".svg")), width = 18, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_0.05", taxon, ".tiff")), width = 18, height = 10, dpi = 300)
sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "History.of.smoking.y.n", palette = smoker_col, n = 3) 
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_0.1", taxon, ".svg")), width = 18, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_0.1", taxon, ".tiff")), width = 18, height = 12, dpi = 300)

taxon = "Genus"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Histology.NSCLC_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/Histology.NSCLC_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat <- read.delim("C:/Users/posadas/Desktop/GitHub/LML/RESULT~2/63808~1.TAX/Stats/RH8B1F~1/NO_MED~1/RHEA_H~1.N_G/HISTOR~2.TAB", stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = mat[mat$corrected < 0.05,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "History.of.smoking.y.n", palette = smoker_col, n = 3) 
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_0.05_", taxon, ".svg")), width = 22, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_0.05_", taxon, ".tiff")), width = 22, height = 10, dpi = 300)

taxon = "ASV"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_Histology.NSCLC_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/Histology.NSCLC_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat <- read.delim("C:/Users/posadas/Desktop/GitHub/LML/RESULT~2/63808~1.TAX/Stats/RH8B1F~1/NO_MED~1/RHEA_H~1.N_A/HISTOR~2.TAB", stringsAsFactors=TRUE)
mat = mat[mat$corrected < 0.05,]
mat = mat[order(mat$corrected)[1:30],]
sig = mat[mat$corrected < 0.05 & !is.na(mat$corrected),"measure"]
sig
topsig = prune_taxa((rownames(tax_table(physeq_q_re)) %in% sig), physeq_q_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "History.of.smoking.y.n", palette = smoker_col, n = 3) 
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_0.05_", taxon, ".svg")), width = 18, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5_relative_sig_0.05_", taxon, ".tiff")), width = 18, height = 12, dpi = 300)


topsig = prune_taxa((rownames(tax_table(physeq_q_re)) %in% c("ASV10", "ASV34", "ASV52", "ASV53", "ASV56","ASV98",
                                                             "ASV120", "ASV147", "ASV188", "ASV230")), physeq_q_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "History.of.smoking.y.n", palette = smoker_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q5_relative_specific_", taxon, ".svg")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5_relative_specific_", taxon, ".tiff")), width = 16, height = 10, dpi = 300)
#### End ####

#### Q5aa: Does history of smoking have an impact on the lung microbiome in adenocarcinoma? ####
physeq_q <- physeq_tax
physeq_q_re <- physeq_tax_re

# Keep only diseased lung
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Diseased") 

# Remove NA
physeq_q = subset_samples(physeq_q, !is.na(sample_data(physeq_q)$History.of.smoking.y.n))

# Select NSCLC
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis == "NSCLC")

# Select histology
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Histology.NSCLC %in% c("Adenocarcinoma"))

# Remove not present taxa
physeq_q <- prune_taxa(taxa_sums(physeq_q) != 0, physeq_q)

# Keep same samples and taxa for relative abundance
physeq_q_re <- prune_samples(sample_names(physeq_q_re) %in% sample_names(physeq_q), physeq_q_re)
physeq_q_re <- prune_taxa(taxa_names(physeq_q_re) %in% taxa_names(physeq_q), physeq_q_re)

#Significant taxa
taxon = "Phylum"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_History.of.smoking.y.n_adeno_", taxon,
                                                                                 "/History.of.smoking.y.n_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

mat = mat[mat$corrected < 0.1,]
#mat = mat[mat$corrected < 0.1 & mat$Group1 == " -",]
sig = unique(mat[mat$corrected < 0.1,"measure"])
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "History.of.smoking.y.n", palette = smoker_col, n = 1, legend_title = "History of smoking") 
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_", taxon, ".svg")), width = 8, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_", taxon, ".tiff")), width = 8, height = 8, dpi = 300)

taxon = "Family"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_History.of.smoking.y.n_adeno_", taxon,
                                                                                 "/History.of.smoking.y.n_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = mat[mat$corrected < 0.05,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "History.of.smoking.y.n", palette = smoker_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_0.05", taxon, ".svg")), width = 18, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_0.05", taxon, ".tiff")), width = 18, height = 10, dpi = 300)

sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "History.of.smoking.y.n", palette = smoker_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_0.1", taxon, ".svg")), width = 18, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_0.1", taxon, ".tiff")), width = 18, height = 12, dpi = 300)

taxon = "Genus"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_History.of.smoking.y.n_adeno_", taxon,
                                                                                 "/History.of.smoking.y.n_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = mat[mat$corrected < 0.05,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "History.of.smoking.y.n", palette = smoker_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_0.05_", taxon, ".svg")), width = 22, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_0.05_", taxon, ".tiff")), width = 22, height = 10, dpi = 300)

taxon = "ASV"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_History.of.smoking.y.n_adeno_", taxon,
                                                                                 "/History.of.smoking.y.n_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat = mat[mat$corrected < 0.05,]
mat = mat[order(mat$corrected)[1:30],]
sig = mat[mat$corrected < 0.05 & !is.na(mat$corrected),"measure"]
sig
topsig = prune_taxa((rownames(tax_table(physeq_q_re)) %in% sig), physeq_q_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "History.of.smoking.y.n", palette = smoker_col, n = 3) 
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_0.05_", taxon, ".svg")), width = 18, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_sig_0.05_", taxon, ".tiff")), width = 18, height = 12, dpi = 300)


topsig = prune_taxa((rownames(tax_table(physeq_q_re)) %in% c("ASV10", "ASV34", "ASV52", "ASV53", "ASV56","ASV98",
                                                             "ASV120", "ASV147", "ASV188", "ASV230")), physeq_q_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "History.of.smoking.y.n", palette = smoker_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_specific_", taxon, ".svg")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q5aa_relative_specific_", taxon, ".tiff")), width = 16, height = 10, dpi = 300)
#### End ####

#### Q6b: Does the N stage of NSCLC have an impact on the lung microbiome ####
physeq_q <- physeq_tax
physeq_q_re <- physeq_tax_re

# Keep only diseased lung
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Diseased") 

# Remove NA
physeq_q = subset_samples(physeq_q, !is.na(sample_data(physeq_q)$N))

# Select NSCLC
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis == "NSCLC")

# Remove not present taxa
physeq_q <- prune_taxa(taxa_sums(physeq_q) != 0, physeq_q)

# Keep same samples and taxa for relative abundance
physeq_q_re <- prune_samples(sample_names(physeq_q_re) %in% sample_names(physeq_q), physeq_q_re)
physeq_q_re <- prune_taxa(taxa_names(physeq_q_re) %in% taxa_names(physeq_q), physeq_q_re)

# Create phyloseq object for diseased and normal lung
physeq_N0 <- subset_samples(physeq_q, sample_data(physeq_q)$N == 0)
physeq_N1 <- subset_samples(physeq_q, sample_data(physeq_q)$N == 1)
physeq_N2 <- subset_samples(physeq_q, sample_data(physeq_q)$N == 2)
physeq_N3 <- subset_samples(physeq_q, sample_data(physeq_q)$N == 3)
physeq_N <- phy_cond(physeq_q, "N", "absolute")
physeq_N_re <- phy_cond(physeq_q_re, "N")

# Define top n ASVs
n = 100
top <- names(sort(taxa_sums(physeq_q), decreasing=TRUE))[1:n]  # adjust number to wished top ones
top_diseased <- names(sort(taxa_sums(physeq_A), decreasing=TRUE))[1:n] 
top_normal <- names(sort(taxa_sums(physeq_B), decreasing=TRUE))[1:n] 
topasv <- prune_taxa(top, physeq_smoking)

# Absolute abundance
plot_bar(topasv, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Phylum"), fill="Phylum")+ facet_wrap(~M, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

# Top n Phyla
n = 10

top = names(taxasum(physeq_q, "Phylum"))[1:n]
top_N0 = names(taxasum(physeq_N0, "Phylum"))[1:n]
top_N1 = names(taxasum(physeq_N1, "Phylum"))[1:n]
top_N2 = names(taxasum(physeq_N0, "Phylum"))[1:n]
top_N3 = names(taxasum(physeq_N1, "Phylum"))[1:n]
print(cbind(top, top_N0, top_N1, top_N2, top_N3))

topphy = prune_taxa((tax_table(physeq_N)[, "Phylum"] %in% top), physeq_N)
topphyrel = prune_taxa((tax_table(physeq_N_re)[, "Phylum"] %in% top), physeq_N_re)

# Absolute abundance
plot_bar(topphy, y =  "Abundance", title = paste("Absolute abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~N, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q6b_absolute_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q6b_absolute_phyla.svg"), width = 10, height = 8, dpi = 300)

#Relative Abundance
plot_bar(topphyrel, y =  "Abundance", title = paste("Relative abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~N, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q6b_relative_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q6b_relative_phyla.svg"), width = 10, height = 8, dpi = 300)

facets_condition(physeq_N_re, "N", n, "Phylum")
ggsave(file.path(bar.dir, "P.Q6b_relative2_phyla.tiff"), width = 10, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "P.Q6b_relative2_phyla.svg"), width = 10, height = 7, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Phylum"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Phylum", Color = "N", palette = smoker_col, n = 2) + thm.x.2

ggsave(file.path(bar.dir, "P.Q6b_relative_box_phyla.tiff"), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, "P.Q6b_relative_box_phyla.svg"), width = 16, height = 10, dpi = 300)

# Top n Genera
n = 10

top = names(taxasum(physeq_q, "Genus"))[1:n]
top_N0 = names(taxasum(physeq_N0, "Genus"))[1:n]
top_N1 = names(taxasum(physeq_N1, "Genus"))[1:n]
top_N2 = names(taxasum(physeq_N0, "Genus"))[1:n]
top_N3 = names(taxasum(physeq_N1, "Genus"))[1:n]
print(cbind(top, top_N0, top_N1, top_N2, top_N3))

topgen = prune_taxa((tax_table(physeq_N)[, "Genus"] %in% top), physeq_N)
topgenrel = prune_taxa((tax_table(physeq_N_re)[, "Genus"] %in% top), physeq_N_re)

# Absolute abundance
plot_bar(topgen, y =  "Abundance", title = paste("Absolute abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~N, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q6b_absolute_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6b_absolute_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

#Relative Abundance
plot_bar(topgenrel, y =  "Abundance", title = paste("Relative abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~N, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q6b_relative_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6b_relative_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

facets_condition(physeq_N_re, "N", n, "Genus")
ggsave(file.path(bar.dir, paste("P.Q6b_relative2_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6b_relative2_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Genus"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Genus", Color = "N", palette = smoker_col, n = 2) + thm.x.2

ggsave(file.path(bar.dir, paste("P.Q6b_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6b_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)

#### End ####

#### Q6c: Does the M stage of NSCLC have an impact on the lung microbiome ####
physeq_q <- physeq_tax
physeq_q_re <- physeq_tax_re

# Keep only diseased lung
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Lung == "Diseased") 

# Remove NA
physeq_q = subset_samples(physeq_q, !is.na(sample_data(physeq_q)$M))

# Select NSCLC
physeq_q = subset_samples(physeq_q, sample_data(physeq_q)$Diagnosis == "NSCLC")

# Remove not present taxa
physeq_q <- prune_taxa(taxa_sums(physeq_q) != 0, physeq_q)

# Keep same samples and taxa for relative abundance
physeq_q_re <- prune_samples(sample_names(physeq_q_re) %in% sample_names(physeq_q), physeq_q_re)
physeq_q_re <- prune_taxa(taxa_names(physeq_q_re) %in% taxa_names(physeq_q), physeq_q_re)

# Create phyloseq object for diseased and normal lung
physeq_M0 <- subset_samples(physeq_q, sample_data(physeq_q)$M == 0)
physeq_M1 <- subset_samples(physeq_q, sample_data(physeq_q)$M == 1)
physeq_M <- phy_cond(physeq_q, "M", "absolute")
physeq_M_re <- phy_cond(physeq_q_re, "M")

# Define top n ASVs
n = 100
top <- names(sort(taxa_sums(physeq_q), decreasing=TRUE))[1:n]  # adjust number to wished top ones
topasv <- prune_taxa(top, physeq_M)

# Absolute abundance
plot_bar(topasv, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Phylum"), fill="Phylum")+ facet_wrap(~M, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

# Top n Phyla
n = 10

top = names(taxasum(physeq_q, "Phylum"))[1:n]
top_M0 = names(taxasum(physeq_M0, "Phylum"))[1:n]
top_M1 = names(taxasum(physeq_M1, "Phylum"))[1:n]
print(cbind(top, top_M0, top_M1))

topphy = prune_taxa((tax_table(physeq_M)[, "Phylum"] %in% top), physeq_M)
topphyrel = prune_taxa((tax_table(physeq_M_re)[, "Phylum"] %in% top), physeq_M_re)

# Absolute abundance
plot_bar(topphy, y =  "Abundance", title = paste("Absolute abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~M, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q6c_absolute_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q6c_absolute_phyla.svg"), width = 10, height = 8, dpi = 300)

#Relative Abundance
plot_bar(topphyrel, y =  "Abundance", title = paste("Relative abundance from top", n, "Phyla"), fill="Phylum")+ facet_wrap(~M, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, "P.Q6c_relative_phyla.tiff"), width = 10, height = 8, dpi = 300)
ggsave(file.path(bar.dir, "P.Q6c_relative_phyla.svg"), width = 10, height = 8, dpi = 300)

facets_condition(physeq_M_re, "M", n, "Phylum")
ggsave(file.path(bar.dir, "P.Q6c_relative2_phyla.tiff"), width = 10, height = 7, dpi = 300)
ggsave(file.path(bar.dir, "P.Q6c_relative2_phyla.svg"), width = 10, height = 7, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Phylum"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Phylum", Color = "M", palette = smoker_col, n = 2) + thm.x.2

ggsave(file.path(bar.dir, "P.Q6c_relative_box_phyla.tiff"), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, "P.Q6c_relative_box_phyla.svg"), width = 16, height = 10, dpi = 300)

# Top n Genera
n = 10

top = names(taxasum(physeq_q, "Genus"))[1:n]
top_M0 = names(taxasum(physeq_M0, "Genus"))[1:n]
top_M1 = names(taxasum(physeq_M1, "Genus"))[1:n]
print(cbind(top, top_M0, top_M1))

topgen = prune_taxa((tax_table(physeq_M)[, "Genus"] %in% top), physeq_M)
topgenrel = prune_taxa((tax_table(physeq_M_re)[, "Genus"] %in% top), physeq_M_re)

# Absolute abundance
plot_bar(topgen, y =  "Abundance", title = paste("Absolute abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~M, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q6c_absolute_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6c_absolute_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

#Relative Abundance
plot_bar(topgenrel, y =  "Abundance", title = paste("Relative abundance from\ntop", n, "Genera"), fill="Genus")+ facet_wrap(~M, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(bar.dir, paste("P.Q6c_relative_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6c_relative_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

facets_condition(physeq_M_re, "M", n, "Genus")
ggsave(file.path(bar.dir, paste("P.Q6c_relative2_genus_", n, ".tiff")), width = 10, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6c_relative2_genus_", n, ".svg")), width = 10, height = 12, dpi = 300)

toprel = prune_taxa((tax_table(physeq_q_re)[, "Genus"] %in% top), physeq_q_re)
plot_abundance(toprel, Facet = "Genus", Color = "M", palette = smoker_col, n = 2) + thm.x.2

ggsave(file.path(bar.dir, paste("P.Q6c_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6c_relative_box_genus_", n, ".tiff")), width = 16, height = 10, dpi = 300)

#Significant taxa
taxon = "Phylum"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_M_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/M_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)

mat[mat$corrected < 0.1,]
sig = unique(mat[mat$corrected < 0.1,"measure"])
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "M", palette = smoker_col, n = 1) + thm.x.2 
ggsave(file.path(bar.dir, paste("P.Q6c_relative_sig_", taxon, ".svg")), width = 12, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6c_relative_sig_", taxon, ".tiff")), width = 12, height = 8, dpi = 300)

taxon = "Family"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_M_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/M_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "M", palette = smoker_col, n = 2) + thm.x.2 
ggsave(file.path(bar.dir, paste("P.Q6c_relative_sig_0.1", taxon, ".svg")), width = 18, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6c_relative_sig_0.1", taxon, ".tiff")), width = 18, height = 12, dpi = 300)

taxon = "Genus"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_M_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/M_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat[mat$corrected < 0.1,]
sig = mat[mat$corrected < 0.1,"measure"]
sig
topsig = prune_taxa((tax_table(physeq_q_re)[, taxon] %in% sig), physeq_q_re)
plot_abundance(topsig, Facet = taxon, Color = "M", palette = smoker_col, n = 1) + thm.x.2 
ggsave(file.path(bar.dir, paste("P.Q6c_relative_sig_0.1_", taxon, ".svg")), width = 12, height = 8, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6c_relative_sig_0.1_", taxon, ".tiff")), width = 12, height = 8, dpi = 300)

taxon = "ASV"
mat <- read.delim(file.path(stats.dir,"Rhea_not_paired/no_median_cutoff", paste0("Rhea_M_", taxon,
                                                                                 "_correlation_matrix_relative_2025-01-29/M_", taxon, "_correlation_matrix_relative-sign_pairs.tab")),
                  stringsAsFactors=TRUE)
mat[mat$corrected < 0.05,]
mat = mat[order(mat$corrected)[1:30],]
sig = mat[mat$corrected < 0.1 & !is.na(mat$corrected),"measure"]
sig
topsig = prune_taxa((rownames(tax_table(physeq_q_re)) %in% sig), physeq_q_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "M", palette = smoker_col, n = 3) 
ggsave(file.path(bar.dir, paste("P.Q6c_relative_sig_0.05_", taxon, ".svg")), width = 18, height = 12, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6c_relative_sig_0.05_", taxon, ".tiff")), width = 18, height = 12, dpi = 300)


topsig = prune_taxa((rownames(tax_table(physeq_q_re)) %in% c("ASV10", "ASV34", "ASV52", "ASV53", "ASV56","ASV98",
                                                             "ASV120", "ASV147", "ASV188", "ASV230")), physeq_q_re)
plot_abundance_ASV(topsig, Facet = "OTU", Color = "M", palette = smoker_col, n = 2) 
ggsave(file.path(bar.dir, paste("P.Q6c_relative_specific_", taxon, ".svg")), width = 16, height = 10, dpi = 300)
ggsave(file.path(bar.dir, paste("P.Q6c_relative_specific_", taxon, ".tiff")), width = 16, height = 10, dpi = 300)

#### End ####

#### End ####

#### Plot abundance per condition ####

# Bar plot group Phylum
ps.group_plot = plot_bar(ps.group, y =  "Abundance", fill="Phylum") +
  facet_wrap(~group, scales="free_x", nrow = 2)

ps.group_plot + geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Relative abundance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.06, "npc")) 

# Bar plot group Species

ps.group_plot = plot_bar(ps.group, y =  "Abundance", fill="Species") +
  facet_wrap(~group, scales="free_x", nrow = 4)

ps.group_plot + geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") +
  ylab("Relative abundance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.06, "npc")) 

# Bar plot implant
ps.implant_plot = plot_bar(ps.implant, y =  "Abundance", fill="Phylum") +
  facet_wrap(~implant, scales="free_x")

ps.implant_plot + geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Relative abundance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.06, "npc")) 

# Bar plot tooth
ps.tooth_plot = plot_bar(ps.tooth, y =  "Abundance", fill="Phylum") +
  facet_wrap(~tooth, scales="free_x")

ps.tooth_plot + geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Relative abundance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.06, "npc")) 

# Bar plot time
ps.time_plot = plot_bar(ps.time, y =  "Abundance", fill="Phylum") +
  facet_wrap(~time, scales="free_x", nrow = 1)

ps.time_plot + geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack") +
  ylab("Relative abundance") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 20),
        legend.key.height = unit(0.06, "npc")) 



#### End ####

