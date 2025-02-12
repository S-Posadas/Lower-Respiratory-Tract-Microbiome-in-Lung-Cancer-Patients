
# Save phyloseq object as input file for lefse

#### Package setup ####

library(phyloseq)
library(phyloseqCompanion)
library(dplyr)
library(microViz)
#### End ####

#### Load RData with phyloseq object ####

res.dir <- "Results_trimmed_cutadapt3_20241216"
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir,"5.phyloseq.filtered.RData"))

# Change species names for visualozation in LEfSe

tax <- as.data.frame(tax_table(physeq_re))   
tax[] <- lapply(tax, gsub, pattern = " ", replacement = "_")
tax$Species[grep("_spp.", tax$Species)] <- NA

tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(physeq_re) <- tax; rm(tax)

#### End ####

#### Function to save lefse input without NA ####

lefse_noNA <- function(name){
  physeq_re_lfse <- read.delim(file.path(lefse.dir, paste0(name, ".txt")))
  physeq_re_lfse_noNA <- physeq_re_lfse[-grep("//|NA", physeq_re_lfse$Sample),]
  write.table(physeq_re_lfse_noNA, file.path(lefse.dir, paste0(name, "_noNA.txt")), sep = "\t", quote = F, row.names = F)
  }

#### End ####

#### Save as input for LEfSe #### 

#### Select phyloseq object ####
physeq2lefse0 <- physeq_re

# Add LibrarySizeGroup
sample_data(physeq2lefse0)$LibrarySizeGroup <- ifelse(sample_data(physeq2lefse0)$LibrarySizeDecontam >= 1000, ">=1000", "<1000")

# Remove controls
physeq2lefse0 <- subset_samples(physeq2lefse0, sample_names(physeq2lefse0) %in% c(rownames(sample_info[!sample_info$Sample_or_Control %in% c("Negative control sample", "Positive control sample"),]))) 

# Remove samples with no reads
physeq2lefse0 <- subset_samples(physeq2lefse0, !is.na(colSums(otu_table(physeq2lefse0)))) 

# Keep only directly isolated samples
physeq2lefse0 <- subset_samples(physeq2lefse0, sample_data(physeq2lefse0)$Isolation == "Direct_Isolation") 

# Select main diagnosis
physeq2lefse0 = subset_samples(physeq2lefse0, sample_data(physeq2lefse0)$Diagnosis %in% c("NSCLC", "SCLC", "Benign"))

# Select filter for minimum ASV count

n = 0
length(sample_sums(physeq)[sample_sums(physeq) >=n])
less_n <- names(sample_sums(physeq)[sample_sums(physeq) <n])
physeq2lefse0 = subset_samples(physeq2lefse0, !sample_names(physeq2lefse0) %in% less_n) 
x = length(sample_names(physeq2lefse0))

lefse.dir=file.path(res.dir, paste0("lefse_less", n, "_counts_", x, "_samples"))
dir.create(lefse.dir)

#### End ####

#### All directly isolated samples ####

phy2lefse <- physeq2lefse0

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q1a_lung_p-s"
phyloseq2lefse(phy2lefse, c("Lung"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q1a_lung_p-g"
phyloseq2lefse(phy2lefse, c("Lung"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q0_batch_p-s"
phyloseq2lefse(phy2lefse, c("Batch"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q0_batch_p-g"
phyloseq2lefse(phy2lefse, c("Batch"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q0.1a_libsize_p-s"
phyloseq2lefse(phy2lefse, c("LibrarySizeGroup"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q0.1a_libsize_p-g"
phyloseq2lefse(phy2lefse, c("LibrarySizeGroup"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

#### All directly isolated samples without non synchronous tumor ####

phy2lefse <- physeq2lefse0

# Remove non synchronous tumor
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Synchronous.tumor == "No")

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q2a_lung_p-s"
phyloseq2lefse(phy2lefse, c("Lung"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q2a_lung_p-g"
phyloseq2lefse(phy2lefse, c("Lung"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

#### Directly isolated samples NSCLC ####

phy2lefse <- physeq2lefse0

# Select NSCLC
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Diagnosis == "NSCLC")

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q1b_lung_NSCLC_p-s"

phyloseq2lefse(phy2lefse, c("Lung"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q1b_lung_NSCLC_p-g"

phyloseq2lefse(phy2lefse, c("Lung"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

# Remove non synchronous tumor
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Synchronous.tumor == "No")

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q2b_lung_NSCLC_p-s"

phyloseq2lefse(phy2lefse, c("Lung"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q2b_lung_NSCLC_p-g"

phyloseq2lefse(phy2lefse, c("Lung"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

#### Three main diagnosis - diseased lung - directly isolated samples ####

phy2lefse <- physeq2lefse0

# Keep only diseased lung
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Lung == "Diseased") 

# Select NSCLC
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Diagnosis %in% c("NSCLC", "SCLC", "Benign"))

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q3_diagnosis_p-s"
phyloseq2lefse(phy2lefse, c("Diagnosis"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q3_diagnosis_p-g"
phyloseq2lefse(phy2lefse, c("Diagnosis"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

#### NSCLC - diseased lung - directly isolated samples - history of smoking ####

phy2lefse <- physeq2lefse0

# Keep only diseased lung
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Lung == "Diseased") 

# Select NSCLC
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Diagnosis == "NSCLC")

# Remove NA
phy2lefse = subset_samples(phy2lefse, !is.na(sample_data(phy2lefse)$History.of.smoking.y.n))

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q5a_History.of.smoking.y.n_p-s"

phyloseq2lefse(phy2lefse, c("History.of.smoking.y.n"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q5a_History.of.smoking.y.n_p-g"

phyloseq2lefse(phy2lefse, c("History.of.smoking.y.n"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

#### NSCLC - diseased lung - directly isolated samples - main histologies ####

phy2lefse <- physeq2lefse0

# Select NSCLC
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Diagnosis == "NSCLC")

# Select histology
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"))
levels(sample_data(phy2lefse)$Histology.NSCLC) <- gsub(" ", "_", levels(sample_data(phy2lefse)$Histology.NSCLC))

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q4_histology_p-s"
phyloseq2lefse(phy2lefse, c("Histology.NSCLC"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q4_histology_p-g"
phyloseq2lefse(phy2lefse, c("Histology.NSCLC"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

### Adenocarcinoma - diseased lung - directly isolated samples - history of smoking ####

phy2lefse <- physeq2lefse0

# Keep only diseased lung
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Lung == "Diseased") 

# Select NSCLC
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Diagnosis == "NSCLC")

# Select histology
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Histology.NSCLC == "Adenocarcinoma")

# Remove NA
phy2lefse = subset_samples(phy2lefse, !is.na(sample_data(phy2lefse)$History.of.smoking.y.n))

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q5aa_History.of.smoking.y.n_p-s"

phyloseq2lefse(phy2lefse, c("History.of.smoking.y.n"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q5aa_History.of.smoking.y.n_p-g"

phyloseq2lefse(phy2lefse, c("History.of.smoking.y.n"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

#### NSCLC - diseased lung - directly isolated samples - M ####

phy2lefse <- physeq2lefse0

# Select NSCLC
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Diagnosis == "NSCLC")

# Remove NA
phy2lefse = subset_samples(phy2lefse, !is.na(sample_data(phy2lefse)$M))

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q6c_M_p-s"
phyloseq2lefse(phy2lefse, c("M"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q6c_M_batch_p-s"
phyloseq2lefse(phy2lefse, c("M", "Batch"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q6c_M_p-g"
phyloseq2lefse(phy2lefse, c("M"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q6c_M_batch_p-g"
phyloseq2lefse(phy2lefse, c("M", "Batch"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

#### NSCLC - diseased lung - directly isolated samples - N ####

phy2lefse <- physeq2lefse0

# Select NSCLC
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Diagnosis == "NSCLC")

# Remove NA
phy2lefse = subset_samples(phy2lefse, !is.na(sample_data(phy2lefse)$N))

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q6b_N_p-s"
phyloseq2lefse(phy2lefse, c("N"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q6b_N_batch_p-s"
phyloseq2lefse(phy2lefse, c("N", "Batch"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q6b_N_p-g"
phyloseq2lefse(phy2lefse, c("N"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q6b_N_batch_p-g"
phyloseq2lefse(phy2lefse, c("N", "Batch"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

#### NSCLC - diseased lung - directly isolated samples - T ####

phy2lefse <- physeq2lefse0

# Select NSCLC
phy2lefse = subset_samples(phy2lefse, sample_data(phy2lefse)$Diagnosis == "NSCLC")

# Remove NA
phy2lefse = subset_samples(phy2lefse, !is.na(sample_data(phy2lefse)$T))

# Remove not present taxa
phy2lefse <- prune_taxa(taxa_sums(phy2lefse) != 0, phy2lefse)

name <- "physeq_re_lfse_q6a_T_p-s"
phyloseq2lefse(phy2lefse, c("T"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q6a_T_batch_p-s"
phyloseq2lefse(phy2lefse, c("T", "Batch"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q6a_T_p-g"
phyloseq2lefse(phy2lefse, c("T"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

name <- "physeq_re_lfse_q6a_T_batch_p-g"
phyloseq2lefse(phy2lefse, c("T", "Batch"), file.name = file.path(lefse.dir, paste0(name, ".txt")),
               taxa.levels = c("Phylum", "Class", "Order", "Family", "Genus"),
               transpose.otus = T)
lefse_noNA(name)

#### End ####