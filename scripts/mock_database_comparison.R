#### Introduction ####

# Positive control: comparison of three 16S reference databases with a mock community
# (Supplementary Document 1, Figure 3).
# The ZymoBIOMICS Microbial Community Standard (Zymo Research) was sequenced in four runs
# (7 samples). For each database (SILVA, eHOMD, GTDB), the 40 most abundant ASVs of the
# mock samples are shown at genus level next to the expected composition.
# The taxonomy is the DADA2 assignment only (assignTaxonomy + exact matching with addSpecies).

#### End ####

#### Package setup ####

Sys.setenv(language = "EN")

library(phyloseq)
library(ggplot2)
library(dplyr)
library(openxlsx)

theme_set(theme_bw())

#### End ####

#### Set Directory ####

res.dir <- "Results_direct" # Same as res.dir in 1.DADA2_direct.R (the mock samples are in both runs; the figure uses this one)
R.dir <- file.path(res.dir, "RData")
pc.dir <- file.path(res.dir, "1.QC", "PositiveControl")
dir.create(pc.dir, recursive = T, showWarnings = F)

# Load only the objects needed from the checkpoint
dada2_objects <- new.env()
load(file.path(R.dir, "2.dada2_taxonomy.RData"), envir = dada2_objects)

#### End ####

#### Parameters ####

n_top <- 40   # Number of most abundant ASVs shown per database
databases <- c("SILVA", "eHOMD", "GTDB")   # Order of the panels

# Expected composition of the mock community (ZymoBIOMICS D6305, theoretical 16S rRNA gene copies, %),
# read from the shared phyloseq object (created by /mnt/data/databases/mocks/Zymo_Mock_Community_save_rds.R)
expected_physeq <- readRDS("/mnt/data/databases/mocks/physeq_expected_zymo_mock_D6305_16S.RDS")
expected <- data.frame(
  Genus     = as.vector(tax_table(expected_physeq)[, "Genus"]),
  Species   = as.vector(tax_table(expected_physeq)[, "Species"]),
  Abundance = as.vector(otu_table(expected_physeq)[, "Expected"])
)

text_size <- 30

#### End ####

#### Count table of the mock samples ####

count_tab <- dada2_objects$count_tab
sample_info <- dada2_objects$sample_info

mock_samples <- rownames(sample_info)[grepl("Mock", sample_info$Study_Nr)]
mock_samples <- intersect(mock_samples, colnames(count_tab))
length(mock_samples)

# Relative abundance (%) in each mock sample, computed on all ASVs of the sample.
# (Counts are first scaled to the median sequencing depth of all samples and rounded,
# as in the original analysis; this only affects the result through rounding.)
total <- median(colSums(count_tab))
count_norm <- apply(count_tab, 2, function(x) round(total * (x / sum(x))))
rel_mock <- sweep(count_norm[, mock_samples], 2, colSums(count_norm[, mock_samples]), "/") * 100

# The 40 most abundant ASVs (total counts over the mock samples); the same ASVs for every database
top_asvs <- names(sort(rowSums(count_tab[, mock_samples]), decreasing = TRUE))[1:n_top]

#### End ####

#### Taxonomy of each database ####

# Taxonomy tables of assignTaxonomy + addSpecies, with ASV names as rownames
asv_names <- setNames(sub(">", "", dada2_objects$asv_headers), dada2_objects$asv_seqs)
get_taxonomy <- function(tab) {
  tab <- as.data.frame(tab)
  if (all(rownames(tab) %in% names(asv_names))) rownames(tab) <- asv_names[rownames(tab)]
  tab
}
tax_tables <- list(SILVA = dada2_objects$taxa_species_SILVA,
                   GTDB  = dada2_objects$taxa_species_GTDB,
                   eHOMD = dada2_objects$taxa_species_eHOMD)
tax_tables <- lapply(tax_tables, get_taxonomy)
stopifnot(all(sapply(tax_tables, function(t) all(top_asvs %in% rownames(t)))))

#### End ####

#### Data for the figure ####

sample_label <- setNames(sample_info[mock_samples, "Study_Nr"], mock_samples)

observed <- bind_rows(lapply(databases, function(db) {
  data.frame(Database = paste(db, "database"),
             Sample = rep(paste(sample_label[mock_samples], db, sep = "_"), each = n_top),
             ASV = rep(top_asvs, times = length(mock_samples)),
             Genus = tax_tables[[db]][rep(top_asvs, times = length(mock_samples)), "Genus"],
             Abundance = as.vector(rel_mock[top_asvs, ]))
}))

plot_data <- bind_rows(
  data.frame(Database = "Expected", Sample = "Expected", ASV = expected$Species,
             Genus = expected$Genus, Abundance = expected$Abundance),
  observed)
plot_data$Database <- factor(plot_data$Database, levels = c("Expected", paste(databases, "database")))

write.xlsx(plot_data, file.path(pc.dir, "Mock_database_comparison_rel_genus.xlsx"))

#### End ####

#### Figure ####

# Each segment is one ASV, colored by its genus (NA = genus not assigned)
ggplot(plot_data, aes(x = Sample, y = Abundance, fill = Genus, color = Genus)) +
  geom_col(position = "stack") +
  facet_wrap(~Database, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Relative abundance (%)") +
  guides(fill = guide_legend(ncol = 4, title = NULL),
         color = guide_legend(ncol = 4, title = NULL)) +
  theme(text = element_text(size = text_size),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

ggsave(file.path(pc.dir, "Mock_database_comparison_rel_genus.svg"), width = 16, height = 8)

#### End ####
