####Introduction####

# This script infers Amplicon Sequence Variants (ASVs) using the DADA2 pipeline,
# for the directly isolated samples and the controls (culture-enriched samples: 1.DADA2_culture.R).
# Adapter sequences have been pre-trimmed using Cutadapt.
# ASVs are inferred per sequencing batch in two passes: a naive pass, and a second pass
# using the denoised forward and reverse sequences of all batches as priors.
# Taxonomy analysis combines DADA2 and Blast approaches.
# The inferred ASVs are then converted into phyloseq objects.
# The environment is saved at different stages of the workflow in RData files
# for subsequent analyses (e.g. unfiltered phyloseq objects for alpha diversity analysis;
# filtered phyloseq objects for taxonomy analyses).

# Checkpoints: each "Save RData" section loads the RData file if it already exists,
# replacing all objects in memory, and only saves it if it does not exist.
# To re-run a step with changed code or data, delete that RData file and all later ones first.

#### End ####

#### Package setup ####

Sys.setenv(language = "EN")

# Core DADA2 and sequence analysis
library(dada2); packageVersion("dada2")
library(Biostrings)

# Phylogenetics and alignment
library(ape)
library(msa)

# Microbiome analysis
library(vegan)
library(phyloseq)
library(decontam)
library(microViz)  

# Data manipulation
library(dplyr)
library(tibble)   
library(purrr)    

# Visualization
library(ggplot2)
library(gridExtra)
library(cowplot)
library(scales)   
library(gplots) 
library(ggstatsplot) 
theme_set(theme_bw())

# Excel
library(openxlsx)

#### End ####

#### Set Directory ####

path <- "Data/Fastq/trimmed_cutadapt"   # Reads after primer removal (cutadapt), one folder per sequencing run

# Create directories to save results of the analysis 

res.dir <- "Results_direct"
qc.dir <- file.path(res.dir, "1.QC")
tax.dir <- file.path(res.dir, "2.blast")
R.dir <- file.path(res.dir, "RData")

dir.create(qc.dir, recursive = T)
dir.create(tax.dir, recursive = T)
dir.create(R.dir, recursive = T)

#### End ####

#### Quality Control and Filtering####

# Helper function

# Aggregated quality profile (all files in one plot per read direction), saved as a pdf.
# Files with fewer than min_reads reads are skipped (plotQualityProfile fails on empty files).
# Missing files count as 0 reads (filterAndTrim does not write samples that lost all reads).
n_reads <- function(f, max = 100) {          # counts at most 'max' reads, so it is fast
  if (!file.exists(f)) return(0)
  fs <- ShortRead::FastqStreamer(f, n = max); on.exit(close(fs))
  length(ShortRead::yield(fs))
}
save_quality_profiles <- function(Fs = fnFs, Rs = fnRs, plot_name = "quality_profiles.pdf", dir = qc.dir, n = 5000, min_reads = 10){
  ok <- sapply(Fs, n_reads) >= min_reads & sapply(Rs, n_reads) >= min_reads
  if (any(!ok)) message("Skipped (< ", min_reads, " reads): ", paste(basename(Fs[!ok]), collapse = ", "))
  pdf(file.path(dir, plot_name), width = 10, height = 6)
  print(plotQualityProfile(Fs[ok], n = n, aggregate = TRUE) + ggtitle(paste("Forward reads -", sum(ok), "files")))
  print(plotQualityProfile(Rs[ok], n = n, aggregate = TRUE) + ggtitle(paste("Reverse reads -", sum(ok), "files")))
  dev.off()
}

# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq.gz and SAMPLENAME_R2.fastq.gz

# Directly isolated samples and controls; culture-enriched samples (*_Cul) are in 1.DADA2_culture.R
keep_sample <- function(x) !grepl("-Cul", x)

fnFs <- sort(list.files(path, pattern="_R1.fastq.gz$", full.names = TRUE, recursive = T))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz$", full.names = TRUE, recursive = T))
fnFs <- fnFs[keep_sample(sub("_R1.fastq.gz$", "", basename(fnFs)))]
fnRs <- fnRs[keep_sample(sub("_R2.fastq.gz$", "", basename(fnRs)))]

# Extract sample names

sample.names <- sub("_R1.fastq.gz$", "", basename(fnFs))
stopifnot(identical(sample.names, sub("_R2.fastq.gz$", "", basename(fnRs))))

# Inspect read quality profiles before primer removal (raw reads) and save as a pdf

fnFs.raw <- sort(list.files("Data/Fastq", pattern="_R1_001.fastq.gz$", full.names = TRUE, recursive = T))
fnFs.raw <- fnFs.raw[keep_sample(fnFs.raw)]
fnRs.raw <- sub("_R1_001.fastq.gz$", "_R2_001.fastq.gz", fnFs.raw)

save_quality_profiles(Fs = fnFs.raw, Rs = fnRs.raw, plot_name = "quality_profiles_before_trimming.pdf")

# Inspect read quality profiles after primer removal (cutadapt) and save as a pdf

save_quality_profiles(Fs = fnFs, Rs = fnRs, plot_name = "quality_profiles_after_cutadapt.pdf")

# Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim 

# Define truncation lengths 
truncLenF <- 230  # Forward trunc length
truncLenR <- 230  # Reverse trunc length

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncLenF,truncLenR),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Check quality after filtering and save as a pdf

save_quality_profiles(Fs = filtFs, Rs = filtRs, plot_name = "quality_profiles_after_dada2_filter.pdf")

#### End ####

#### Remove empty samples from list ####

# Since we are dealing with low biomass samples, some samples lost all reads after filtering.
filtFs <- sort(list.files(file.path(path, "filtered"), pattern="_F_filt.fastq$", full.names = TRUE))
filtRs <- sort(list.files(file.path(path, "filtered"), pattern="_R_filt.fastq$", full.names = TRUE))

# The filtered/ folder can also contain the samples of the other DADA2 script: keep the samples of this script
filtFs <- filtFs[keep_sample(filtFs)]
filtRs <- filtRs[keep_sample(filtRs)]

# Sample names as in metadata
sample.names <- sub("_.*", "", basename(filtFs))
stopifnot(identical(sample.names, sub("_.*", "", basename(filtRs))))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#### End ####

#### Learn the Error Rates ####
# In the BigData WF (https://benjjneb.github.io/dada2/bigdata.html),
# it is recommended to learn the error rates and do the sample inference for each run individually.
# Thus metadata is already added at this step to know the batch of each sample

sample_info <- read.xlsx("Data/000-LML_Metadata_clinical_Lab_final_20240114.xlsx", rowNames = T, sheet = "Combined_metadata")

# Create data frame containing filtFs, filtRs and batch

filt <- merge(as.data.frame(filtFs), as.data.frame(filtRs), by = 0) %>% merge(., sample_info[rownames(sample_info) %in% names(filtFs), c("nameF", "Batch")], by.x = "Row.names", by.y = 0)

# Learn the error rates
errF <- list()
errR <- list()
for(batch in unique(filt$Batch)){
  errF[[batch]] <- dada2::learnErrors(filt[filt$Batch == batch,"filtFs"], multithread=FALSE)
  errR[[batch]] <- dada2::learnErrors(filt[filt$Batch == batch,"filtRs"], multithread=FALSE)
}

# Plot errors and save as pdf

batchs <- unique(filt$Batch) # Select batch names

lapply(batchs, function(batch) {
  
  pdf(paste0(qc.dir, "/forward_error_rate_", batch, ".pdf"))
  print(plotErrors(errF[[batch]], nominalQ = TRUE))
  dev.off()
  
  pdf(paste0(qc.dir, "/reverse_error_rate_", batch, ".pdf"))
  print(plotErrors(errR[[batch]], nominalQ = TRUE))
  dev.off()
  
})

#### End ####

#### ASV inference, merge and first sequence table ####

# Apply the core sample inference algorithm to the the filtered and trimmed sequence data.
dadaFs <- list()
dadaRs <- list()
for(batch in unique(filt$Batch)){
  dadaFs[[batch]] <- dada(filt[filt$Batch == batch,"filtFs"], err=errF[[batch]], pool = F, multithread=FALSE)
  dadaRs[[batch]] <- dada(filt[filt$Batch == batch,"filtRs"], err=errR[[batch]], pool = F, multithread=FALSE)
}

# Merge paired reads
mergers.naive <- list()

for(batch in unique(filt$Batch)){
  mergers.naive[[batch]] <- mergePairs(dadaFs[[batch]], filt[filt$Batch == batch,"filtFs"],
                                       dadaRs[[batch]], filt[filt$Batch == batch,"filtRs"], verbose=TRUE) # min overlap is 12 as default, but can be adjusted
}

# Obtain sequence tables
seqtab.perbatch.naive <- list()
for(batch in unique(filt$Batch)){
  seqtab.perbatch.naive[[batch]] <- makeSequenceTable(mergers.naive[[batch]])
  message("Batch ", batch, ": ", dim(seqtab.perbatch.naive[[batch]])[1], " samples and ", dim(seqtab.perbatch.naive[[batch]])[2], " sequences")
}

# Merged all sequence tables into one
seqtab.naive <- seqtab.perbatch.naive[[1]]
for (i in 2:length(seqtab.perbatch.naive)) {
  seqtab.naive <- mergeSequenceTables(seqtab.naive, seqtab.perbatch.naive[[i]])
}

# Remove chimeras from the merged table of all batches
seqtab.nochim.naive <- removeBimeraDenovo(seqtab.naive, method="consensus", multithread=FALSE, verbose=TRUE)

dim(seqtab.naive)
dim(seqtab.nochim.naive)

# Percentage of not chimeric reads

sum(seqtab.nochim.naive)/sum(seqtab.naive)

# Sequence tables of the denoised forward and reverse reads (used as priors)
seqtabFs.perbatch <- list()
seqtabRs.perbatch <- list()
for(batch in unique(filt$Batch)){
  seqtabFs.perbatch[[batch]] <- makeSequenceTable(dadaFs[[batch]])
  seqtabRs.perbatch[[batch]] <- makeSequenceTable(dadaRs[[batch]])
  message("Batch ", batch, ": ", dim(seqtabFs.perbatch[[batch]])[2], " forward and ", dim(seqtabRs.perbatch[[batch]])[2], " reverse sequences")
}

# Merged all sequence tables into one
seqtabFs <- seqtabFs.perbatch[[1]]
seqtabRs <- seqtabRs.perbatch[[1]]
for (i in 2:length(seqtabFs.perbatch)) {
  seqtabFs <- mergeSequenceTables(seqtabFs, seqtabFs.perbatch[[i]])
  seqtabRs <- mergeSequenceTables(seqtabRs, seqtabRs.perbatch[[i]])
}

#### End ####

#### Second ASV inference with list of prior ASVs, merge and final sequence table ####

# Prepare priors for second pass: denoised forward and reverse sequences of all batches.
# Reverse priors come from the reverse reads, so they are already in the orientation of the reverse reads.
priorF <- colnames(seqtabFs)
priorR <- colnames(seqtabRs)

# Rerun dada with priors
dadaFs.prior <- list()
dadaRs.prior <- list()
for(batch in unique(filt$Batch)){
  dadaFs.prior[[batch]] <- dada(filt[filt$Batch == batch,"filtFs"], err=errF[[batch]], priors = priorF, multithread=FALSE)
  dadaRs.prior[[batch]] <- dada(filt[filt$Batch == batch,"filtRs"], err=errR[[batch]], priors = priorR, multithread=FALSE)
}

#Merge paired reads
mergers <- list()

for(batch in unique(filt$Batch)){
  mergers[[batch]] <- mergePairs(dadaFs.prior[[batch]], filt[filt$Batch == batch,"filtFs"],
                                 dadaRs.prior[[batch]], filt[filt$Batch == batch,"filtRs"], verbose=TRUE) # min overlap is 12 as default, but can be adjusted
}

# Most of your reads should successfully merge. If that is not the case upstream parameters may need
# to be revisited: Did you trim away the overlap between your reads?

# Obtain sequence tables
seqtab.perbatch <- list()
for(batch in unique(filt$Batch)){
  seqtab.perbatch[[batch]] <- makeSequenceTable(mergers[[batch]])
  message("Batch ", batch, ": ", dim(seqtab.perbatch[[batch]])[1], " samples and ", dim(seqtab.perbatch[[batch]])[2], " sequences")
}

# Merged all sequence tables into one
seqtab <- seqtab.perbatch[[1]]
for (i in 2:length(seqtab.perbatch)) {
  seqtab <- mergeSequenceTables(seqtab, seqtab.perbatch[[i]])
}

# Remove chimeras from the merged table of all batches
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

dim(seqtab)
dim(seqtab.nochim)

# Percentage of not chimeric reads

sum(seqtab.nochim)/sum(seqtab)

#### End ####

#### Track Pipeline ####
#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))

track <- list()
for(batch in unique(filt$Batch)){
  track[[batch]] <- cbind(sapply(dadaFs[[batch]], getN), sapply(dadaRs[[batch]], getN), sapply(dadaFs.prior[[batch]], getN), sapply(dadaRs.prior[[batch]], getN), sapply(mergers.naive[[batch]], getN), sapply(mergers[[batch]], getN))
}
track <- as.data.frame(do.call(rbind, track))

rownames(out) <- sub("_R1.fastq.gz$", "", rownames(out))
rownames(track) <- sub("_F_filt.fastq$", "", rownames(track))
rownames(seqtab.nochim) <- sub("_F_filt.fastq$", "", rownames(seqtab.nochim))
rownames(seqtab.nochim.naive) <- sub("_F_filt.fastq$", "", rownames(seqtab.nochim.naive))
track <- merge(as.data.frame(out), as.data.frame(track), by = 0, all.x = T) %>% merge(.,rowSums(seqtab.nochim.naive), by.x = "Row.names", by.y = 0) %>% merge(.,rowSums(seqtab.nochim), by.x = "Row.names", by.y = 0)

colnames(track) <- c("sample","input", "filtered", "denoisedF.naive", "denoisedR.naive", "denoisedF.prior", "denoisedR.prior", "merged.naive", "merged.prior", "nonchim.naive", "nonchim.prior")

head(track)

# Save the track

write.csv(track, file = file.path(qc.dir, "Pipeline_Track.csv"))

#### End ####

#### Save RData ####

# Save as .RData to load in following steps or continue later
checkpoint1 <- file.path(R.dir, "1.dada2_track.RData")
if (file.exists(checkpoint1)) {
  load(checkpoint1)
} else {
  save.image(checkpoint1)
}

#### End ####

#### Assign taxonomy for Bacteria ####

# Helper function
assign_taxonomy <- function(seqtab, ref_fasta, ref_fasta_species){
  set.seed(123456789)
  taxa <- assignTaxonomy(seqtab, ref_fasta, tryRC = TRUE) #tryRC = TRUE -> reverse-complement orientation
  taxa_species = addSpecies(taxa, ref_fasta_species, allowMultiple=TRUE)
  return(list(taxa = taxa, taxa_species = taxa_species))
}

#### GTDB
ref_fasta_GTDB = "/mnt/data/databases/dada2/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz"
ref_fasta_species_GTDB = "/mnt/data/databases/dada2/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc53_ssu_r220_species.fa.gz"

taxa_GTDB <- assign_taxonomy(seqtab.nochim, ref_fasta_GTDB, ref_fasta_species_GTDB)
taxa_species_GTDB = taxa_GTDB[["taxa_species"]]

head(unname(taxa_GTDB[["taxa"]])) # Removing sequence rownames for display only
head(unname(taxa_species_GTDB))

#### SILVA

ref_fasta_SILVA = "/mnt/data/databases/dada2/16S_Taxonomy_Databases/SILVA_2024/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
ref_fasta_species_SILVA = "/mnt/data/databases/dada2/16S_Taxonomy_Databases/SILVA_2024/silva_v138.2_assignSpecies.fa.gz"

taxa_SILVA <- assign_taxonomy(seqtab.nochim, ref_fasta_SILVA, ref_fasta_species_SILVA)
taxa_species_SILVA = taxa_SILVA[["taxa_species"]]

head(unname(taxa_SILVA[["taxa"]])) # Removing sequence rownames for display only
head(unname(taxa_species_SILVA))

#### eHOMD

ref_fasta_eHOMD = "/mnt/data/databases/dada2/16S_Taxonomy_Databases/eHOMD/eHOMD_RefSeq_dada2_V15.22.fasta.gz"
ref_fasta_species_eHOMD = "/mnt/data/databases/dada2/16S_Taxonomy_Databases/eHOMD/eHOMD_RefSeq_dada2_assign_species_V15.22.fasta.gz"

taxa_eHOMD<- assign_taxonomy(seqtab.nochim, ref_fasta_eHOMD, ref_fasta_species_eHOMD)
taxa_species_eHOMD = taxa_eHOMD[["taxa_species"]]

head(unname(taxa_eHOMD[["taxa"]])) # Removing sequence rownames for display only
head(unname(taxa_species_eHOMD))

#### End ####

#### Make count table from naive dada ####

# Giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs.n <- colnames(seqtab.nochim.naive)

asv_headers.n <- vector(dim(seqtab.nochim.naive)[2], mode="character")

for (i in 1:dim(seqtab.nochim.naive)[2]) {
  asv_headers.n[i] <- paste(">ASV", i, sep="_")
}

# Making fasta of our final ASV seqs:

fasta_tab.naive <- c(rbind(asv_headers.n, asv_seqs.n))
write(fasta_tab.naive, file.path(res.dir, "fasta_tab_naive.fa"))

# count table:

count_tab.naive <- t(seqtab.nochim.naive)
row.names(count_tab.naive) <- sub(">", "", asv_headers.n)
colnames(count_tab.naive) <- sub("_.*", "", colnames(count_tab.naive)) # Sample names as in metadata (same rule as in "Remove empty samples from list")
write.table(count_tab.naive, file.path(res.dir, "count_tab_naive.tsv"), sep="\t", quote=F, col.names=NA)

#### End ####

#### Make count and taxa tables from dada with priors ####

# Giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochim)

asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Making fasta of our final ASV seqs:

fasta_tab <- c(rbind(asv_headers, asv_seqs))
write(fasta_tab, file.path(res.dir, "fasta_tab.fa"))

# count table:

count_tab <- t(seqtab.nochim)
row.names(count_tab) <- sub(">", "", asv_headers)
colnames(count_tab) <- sub("_.*", "", colnames(count_tab)) # Sample names as in metadata (same rule as in "Remove empty samples from list")
write.table(count_tab, file.path(res.dir, "count_tab.tsv"), sep="\t", quote=F, col.names=NA)

# tax table:
tax_tab_GTDB <- taxa_species_GTDB
row.names(tax_tab_GTDB) <- sub(">", "", asv_headers)
write.xlsx(as.data.frame(tax_tab_GTDB), file.path(tax.dir,"tax_tab_GTDB.xlsx"), rowNames = T)

tax_tab_SILVA <- taxa_species_SILVA
row.names(taxa_species_SILVA) <- sub(">", "", asv_headers)
write.xlsx(as.data.frame(taxa_species_SILVA), file.path(tax.dir,"tax_tab_SILVA.xlsx"), rowNames = T)

tax_tab_eHOMD <- taxa_species_eHOMD
row.names(tax_tab_eHOMD) <- sub(">", "", asv_headers)
write.xlsx(as.data.frame(tax_tab_eHOMD), file.path(tax.dir,"tax_tab_eHOMD.xlsx"), rowNames = T)

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

checkpoint2 <- file.path(R.dir, "2.dada2_taxonomy.RData")
if (file.exists(checkpoint2)) {
  load(checkpoint2)
} else {
  save.image(checkpoint2)
}

#### End ####

#### Blast unassigned species with scripts blast and filter_blast.R ####

# STOP: run blast and filter_blast.R before continuing.
#   blast (run inside res.dir, needs BLASTDB): input fasta_tab.fa; output 2.blast/tax_tab_95_10align_16S
#   filter_blast.R input: 2.blast/tax_tab_95_10align_16S, 2.blast/tax_tab_SILVA.xlsx and fasta_tab.fa;
#   output: 2.blast/tax_tab_SILVA_blast_species.xlsx (read below)

tax_tab <- read.xlsx(file.path(tax.dir,"tax_tab_SILVA_blast_species.xlsx"), rowNames = T)
tax_tab <- tax_tab[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

# Paste Genus and Species name for later visualization
colnames(tax_tab)[colnames(tax_tab) == "Species"] <- "Species.only"
tax_tab$Species <- NA  # Initialize Species column with NA

tax_tab$Species <- ifelse(!is.na(tax_tab$Genus) & !is.na(tax_tab$Species.only),
                          paste(tax_tab$Genus, tax_tab$Species.only),
                          ifelse(!is.na(tax_tab$Genus) & is.na(tax_tab$Species.only),
                                 paste(tax_tab$Genus, "spp."), NA))
tax_tab$Species.only <- NULL

#### End ####

#### Add samples metadata and match to count table####

sample_info <- sample_info[is.na(sample_info$Isolation) | sample_info$Isolation == "Direct_Isolation",]

## Decode sample_info data and convert to factors

# Keep the original numeric codes in new ".code" columns
# (Sex.code also holds the reason why a patient was excluded, see "Remove samples that don´t fit the study conditions")
coded_vars <- c("Sex", "History.of.smoking", "Active.smoker", "Side", "Lobe", "Synchronous.tumor")
sample_info[paste0(coded_vars, ".code")] <- sample_info[coded_vars]

# Labels for each code. The order of the labels is the order of the factor levels.
# Values without a label (e.g. exclusion notes) become NA.
decode <- list(
  Sex                = c("0" = "Female", "1" = "Male"),
  History.of.smoking = c("0" = "None", "1" = "Cigarettes", "2" = "Pipe"),
  Active.smoker      = c("0" = "No", "1" = "Yes"),
  Side               = c("1" = "Left", "0" = "Right"),
  Lobe               = c("1" = "Upper", "2" = "Lower", "3" = "Middle", "4" = "Lingula"),
  Synchronous.tumor  = c("0" = "No", "1" = "Second tumor contralateral", "2" = "Thoracic metastasis of the primary tumor")
)
for (v in names(decode)) {
  sample_info[[v]] <- factor(unname(decode[[v]][as.character(sample_info[[v]])]), levels = unique(decode[[v]]))
}

# History of smoking as yes/no (cigarettes or pipe = yes)
smoking_y_n <- c("0" = "No", "1" = "Yes", "2" = "Yes")
sample_info$History.of.smoking.y.n <- factor(unname(smoking_y_n[as.character(sample_info$History.of.smoking.code)]), levels = c("No", "Yes"))

# Clinical variables already stored as labels: set the order of the levels
#sample_info$Isolation <- factor(sample_info$Isolation, levels = c("Direct_Isolation", "Culture_Enriched"))
sample_info$Lung <- factor(sample_info$Lung, levels = c("Normal", "Diseased"))
sample_info$Histology.NSCLC <- factor(sample_info$Histology.NSCLC, levels = c("Adenocarcinoma", "Squamous cell carcinoma", "Adeno-Squamous cell carcinoma", "Others"))
sample_info$Diagnosis <- factor(sample_info$Diagnosis, levels = c("Benign", "NSCLC", "V.a. NSCLC", "SCLC", "Atypical Carzinoid", "Carcinoid Tumorlet", "Metastasized Tumor"))
sample_info$T <- factor(sample_info$T)
sample_info$N <- factor(sample_info$N)
sample_info$M <- factor(sample_info$M)

## Colors for each variable (saved in the RData files and used by the downstream scripts)

palette2 <- rev(scales::hue_pal()(2))                 # 2 groups: "#00BFC4", "#F8766D" (used by downstream scripts)
palette3 <- c("#00BFC4", "#F8766D", "#00BA38")        # smoking
palette  <- c("#53B400", "#F8766D", "#00C094", "#344cb7", "#00B6EB", "#A58AFF", "#FB61D7") # diagnosis, histology

# Named vector of colors, one per level of a factor in sample_info
level_colors <- function(variable, pal) {
  lv <- levels(sample_info[[variable]])
  setNames(pal[seq_along(lv)], lv)
}

isolation_col <- c(Direct_Isolation = "#1b9e77", Culture_Enriched = "#d95f02")
lung_col      <- c(Normal = "#0072B2", Diseased = "#D55E00")   
diagnosis_col <- level_colors("Diagnosis", palette)
smokerh_col   <- level_colors("History.of.smoking", palette3)
smoker_col    <- level_colors("Active.smoker", palette3)   # also used for History.of.smoking.y.n (same levels No/Yes)
histology_col <- level_colors("Histology.NSCLC", palette)
tumor_col     <- setNames(c("#53B400", "#F8766D", "#00B6EB", "#A58AFF", "#FB61D7"), 0:4) # not used by the downstream scripts

# Examine all tables

head(count_tab)
head(tax_tab)
head(fasta_tab)
head(sample_info)

# Examine consistancy in order between count_tab colnames and coldata rownames 
all(rownames(sample_info) %in% colnames(count_tab))
all(colnames(count_tab) %in% rownames(sample_info))
all(rownames(sample_info) == colnames(count_tab)) 

# Not all samples in sample_info are in the counts table because some samples have no reads after filtering

rownames(sample_info)[!rownames(sample_info) %in% colnames(count_tab)]

noreads <- matrix(0, ncol = length(rownames(sample_info)[!rownames(sample_info) %in% colnames(count_tab)]),
                  nrow = length(rownames(count_tab)),
                  dimnames = list(rownames(count_tab), rownames(sample_info)[!rownames(sample_info) %in% colnames(count_tab)]))
count_tab <- cbind(count_tab, noreads)
noreads.n <- matrix(0, ncol = length(rownames(sample_info)[!rownames(sample_info) %in% colnames(count_tab.naive)]),
                    nrow = length(rownames(count_tab.naive)),
                    dimnames = list(rownames(count_tab.naive), rownames(sample_info)[!rownames(sample_info) %in% colnames(count_tab.naive)]))
count_tab.naive <- cbind(count_tab.naive, noreads.n)

all(rownames(tax_tab) %in% rownames(count_tab))
all(rownames(tax_tab) == rownames(count_tab))

gplots::venn(list(taxonomy=rownames(tax_tab), featuretable=rownames(count_tab)))

#### End ####

#### Phylogenetic tree ####

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")

aligned_asvs <- as(mult, "DNAStringSet")
names(aligned_asvs) <- seqs
writeXStringSet(aligned_asvs, filepath = file.path(res.dir, "aligned_asvs.fasta"), format = "fasta")

# STOP: run the fasttree script before continuing.
#   fasttree (run inside res.dir, needs FastTree): input aligned_asvs.fasta; output tree.nwk
#   FastTree -gtr -gamma -nt aligned_asvs.fasta > tree.nwk
tree <- read.tree(file.path(res.dir, "tree.nwk"))

#### End ####

#### Making our phyloseq object ####
# Examine all tables

head(count_tab)
head(tax_tab)
head(fasta_tab)
head(sample_info)

physeq <- phyloseq(otu_table(count_tab, taxa_are_rows = T),   #taxa_are_rows=F (if your taxa names on the column not the rows)
                   sample_data(sample_info), 
                   tax_table(as.matrix(tax_tab)))

# Adding ASV Fasta sequences and Phylogenetic tree to phyloseq object

dna <- Biostrings::DNAStringSet(asv_seqs)  # Making ASV Fasta sequences 
names(dna) <- sub(">", "", asv_headers)    # ASV name of each sequence (same as rownames(count_tab))

# Rename each tip by matching its sequence
ph_tree = phy_tree(tree)            # Making Phylogenetic tree
ph_tree$tip.label <- names(dna)[match(ph_tree$tip.label, as.character(dna))]
stopifnot(!anyNA(ph_tree$tip.label), setequal(ph_tree$tip.label, taxa_names(physeq)))

# Root the tree at its midpoint. FastTree trees are unrooted, and phyloseq's UniFrac would
# otherwise root them on a randomly chosen ASV, which can differ between analyses.
# The root is kept when ASVs are filtered later (pruned trees stay rooted the same way).
ph_tree <- phangorn::midpoint(ph_tree)
stopifnot(ape::is.rooted(ph_tree))

physeq_dna_tree <- merge_phyloseq(physeq, ph_tree, dna) #Merging  ASV Fasta sequences and Phylogenetic tree to phyloseq object

# Keep the ASV numbers of fasta_tab.fa (ASV_1 -> ASV1); the tree changes the taxa order
taxa_names(physeq_dna_tree) <- sub("^ASV_", "ASV", taxa_names(physeq_dna_tree))
stopifnot(!anyDuplicated(taxa_names(physeq_dna_tree)), all(grepl("^ASV[0-9]+$", taxa_names(physeq_dna_tree))))

physeq
physeq_dna_tree
physeq = physeq_dna_tree

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

checkpoint3 <- file.path(R.dir, "3.physeq.original.RData")
if (file.exists(checkpoint3)) {
  load(checkpoint3)
} else {
  save.image(checkpoint3)
}

#### End ####

## Filter low Abundance Taxa and count table normalization ##

#### Decontamination of ASVs appearing in negative controls ####

# Save initial library size
sample_data(physeq)$LibrarySize <- sample_sums(physeq)

# Select negative controls and delete positive ones

PC <- rownames(sample_info[sample_info$Sample_or_Control == "Positive control sample",])
NC <- rownames(sample_info[sample_info$Sample_or_Control == "Negative control sample",])

ps.decontam <- prune_samples(!sample_names(physeq) %in% PC, physeq)
ps.decontam <- prune_samples(sample_sums(ps.decontam) != 0, ps.decontam)
ps.decontam <- prune_taxa(taxa_sums(ps.decontam) != 0, ps.decontam)

# Decontam per batch, with batches 1-4 pooled into one group
sample_data(ps.decontam)$Batch_pooled <- sample_data(ps.decontam)$Batch
sample_data(ps.decontam)$Batch_pooled[sample_data(ps.decontam)$Batch_pooled %in% c("1-LML_001-006", "2-LML_007-030", "3-LML_031-054", "4-LML_055-077")] <- "1-4 batch"
sample_data(ps.decontam)$is.neg <- sample_data(ps.decontam)$Sample_or_Control == "Negative control sample"

#' The probability threshold below which the null-hypothesis (not a contaminant) should be rejected in favor of
#' the alternate hypothesis (contaminant) was set to 0.5 (default 0.1)
thr = 0.5
contamdf.prev <- isContaminant(ps.decontam, method="prevalence", neg="is.neg", threshold=thr, batch = "Batch_pooled")
table(contamdf.prev$contaminant)

# Supplementary Document 1, Figure 4: distribution of decontam scores across thresholds
hist(contamdf.prev$p, 100, main="P-values for different thresholds in decontam")

svg(file.path(qc.dir, "Decontam_scores_histogram.svg"), width = 10, height = 7)
hist(contamdf.prev$p, 100, main="P-values for different thresholds in decontam")
dev.off()

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps.decontam, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Negative control sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point(size = 2) +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") + 
  theme(text = element_text(size = 30),
        legend.position = "bottom")

# Contaminants0.5 (.tiff/.svg/.xlsx): contaminant ASVs reported in the paper (67 ASVs)
ggsave(file.path(qc.dir, paste0("Contaminants", thr, ".tiff")), width = 10, height = 10, dpi = 300)
ggsave(file.path(qc.dir, paste0("Contaminants", thr, ".svg")), width = 10, height = 10, dpi = 300)

contaminants <- as.data.frame(otu_table(prune_taxa(contamdf.prev$contaminant, ps.decontam)))
contaminants_tax <- merge(tax_tab[rownames(contaminants),], contaminants, by = 0)

write.xlsx(contaminants_tax, file.path(qc.dir, paste0("Contaminants", thr, ".xlsx")))

physeq_decont <- prune_taxa(taxa_sums(physeq) != 0, physeq)
physeq_decont <- prune_taxa(!taxa_names(physeq_decont) %in% rownames(contaminants), physeq_decont)

physeq_decont
physeq
physeq = physeq_decont

# Save library size after decontamination
sample_data(physeq)$LibrarySizeDecontam <- sample_sums(physeq)

#### End ####

#### Make physeq object from naive dada2 to examine library size ####

head(count_tab.naive)
head(fasta_tab.naive)
head(sample_info)

physeq.naive <- phyloseq(otu_table(count_tab.naive, taxa_are_rows = T), sample_data(sample_info))
# physeq.naive = prune_samples(!sample_names(physeq.naive) %in% rownames(samp_ex), physeq.naive)
# physeq.naive <- microViz::ps_reorder(physeq.naive, order(sample_names(physeq.naive)))
# physeq.naive <- prune_taxa(taxa_sums(physeq.naive) != 0, physeq.naive)
sample_data(physeq.naive)$LibrarySize <- sample_sums(physeq.naive)

#### End ####

#### Save sample sums before and after decontamination ####

combined_sums <- data.frame(
  naive = sample_sums(physeq.naive),
  priors = sample_sums(physeq_dna_tree)[match(sample_names(physeq.naive), sample_names(physeq_dna_tree))],
  decont = sample_sums(physeq_decont)[match(sample_names(physeq.naive), sample_names(physeq_decont))]
)

write.xlsx(combined_sums, file.path(qc.dir, "sample_sums.xlsx"), rowNames= T)

sums_to_plot <- combined_sums[
  grepl("LML", rownames(combined_sums)) &
    grepl("A1", rownames(combined_sums)) &
    !grepl("Cul", rownames(combined_sums)), ]

#sums_to_plot <- sums_to_plot[sums_to_plot$decont < 20000, ]

hist(log10(sums_to_plot$decont), 
     main = "Library size after decontamination", 
     xlab = "Log10 Library size", breaks =  200,
     border = "black")

hist(sums_to_plot$decont, 
     main = "Library size after decontamination", 
     xlab = "Library size", breaks =  200,
     border = "black")

plot(density(sums_to_plot$decont), 
     main = "Library size after decontamination", 
     xlab = "Library size", 
     col = "red", 
     lwd = 2)

plot(density(log10(sums_to_plot$decont)), 
     main = "Library size after decontamination", 
     xlab = "Log10 Library size", 
     col = "red", 
     lwd = 2)

#### End ####

#### Remove samples that don´t fit the study conditions ####
samp_ex = sample_data(physeq)[sample_data(physeq)$Exclusion == 1, ]
print(samp_ex)
physeq_study = prune_samples(!sample_names(physeq) %in% rownames(samp_ex), physeq)
physeq_study = prune_taxa(taxa_sums(physeq_study) != 0, physeq_study)

physeq
physeq_study
physeq = physeq_study

# Reorder samples for later paired tests (so Lung A and B have the same position in the vector)
physeq <- microViz::ps_reorder(physeq, order(sample_names(physeq)))

#### End ####

#### Inspect library sizes ####
sample_data(physeq)$LibrarySize <- sample_sums(physeq)

df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() + theme(text = element_text(size = 30))

# Explore library sizes per Diagnosis
breaks <- c(0, 1, 10, 100, 1000, 10000, 50000, Inf)
labels <- c("0", "1-10", "11-100", "101-1000", "1001-10000", "10001-50000", ">50000")
df$LibrarySizeBreaks <- cut(df$LibrarySize, breaks = breaks, labels = labels, include.lowest = TRUE)
df$count = 1
df$LibrarySizeLog10 <- log10(df$LibrarySize)
df <- df[df$Sample_or_Control == "True sample",]

# Box plots
ggplot(data=df, aes(x=Lung, y=LibrarySizeLog10, color=Lung)) +
  geom_boxplot() +
  geom_point() +
  #ggtitle("Library size in directly isolated samples") + 
  facet_wrap(.~Diagnosis, nrow = 2) + theme(text = element_text(size = 30),
                                            axis.ticks.x = element_blank(),
                                            axis.text.x = element_blank())

ggsave(file.path(qc.dir, "LibSizeLog10.svg"), width = 16, height = 10)

#### End ####

#### Calculate and add alpha diversity to phyloseq object for later analysis ####

richness <- estimate_richness(physeq, measures = c("Observed","Chao1", "Shannon","Simpson","InvSimpson","ACE"))
# estimate_richness() changes "-" in sample names to "." (make.names); restore the original names
stopifnot(identical(rownames(richness), make.names(sample_names(physeq))))
rownames(richness) <- sample_names(physeq)
sample_data(physeq) <- merge(sample_data(physeq), richness, by = 0) %>% tibble::column_to_rownames(var = "Row.names")
rm(richness)
sample_data(physeq)$InvSimpson[sample_data(physeq)$InvSimpson == Inf] = NA

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

checkpoint4 <- file.path(R.dir, "4.physeq.decontam.RData")
if (file.exists(checkpoint4)) {
  load(checkpoint4)
} else {
  save.image(checkpoint4)
}
saveRDS(physeq, file.path(R.dir, "physeq_study.RDS"))

#### End ####

#### Rarefaction curves ####

# Supplementary Document 1, Figure 2 (A): rarefaction curves (combined_rarefaction.svg) and
# number of ASVs after each step (asv_number.xlsx): 14,316 naive, 14,954 with priors, 14,887 after decontamination

# Define aesthetics
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")

# Create named list of count tables
count_list <- list(
  ASVs_from_naive_inference = count_tab.naive,
  ASVs_after_inference_with_list_of_priors = count_tab,
  ASVs_after_decontamination = data.frame(otu_table(physeq_decont), check.names = FALSE)  # keep "-" in sample names
)

rare_plots <- list()  # Store ggplot objects
rare_dim <- list()  # Store number of ASVs
                                 
for (name in names(count_list)) {
  count_rar <- count_list[[name]]
  
  # Exclude controls
  count_rar <- count_rar[,colnames(count_rar) %in% rownames(sample_info[sample_info$Sample_or_Control == "True sample",])]
  # Exclude samples with 0 counts
  count_rar <- count_rar[, colSums(count_rar) != 0]
  # Exclude ASVs with 0 counts
  count_rar <- count_rar[rowSums(count_rar) != 0, ]
  
  # Generate rarefaction curves
  result <- rarecurve(t(count_rar), step = 100, label = FALSE)
  
  # Convert to data frame
  rare_df <- map_dfr(seq_along(result), function(i) {
    curve <- result[[i]]
    subsample <- attr(curve, "Subsample")
    
    tibble(
      Sample = paste0("Sample_", i),
      Dataset = name,  # Add dataset identifier
      Individuals = subsample,
      Richness = as.numeric(curve),
      Color = col[(i - 1) %% length(col) + 1],
      LineType = lty[(i - 1) %% length(lty) + 1]
    )
  })
  
  # Create plot
  rare_plot <- ggplot(rare_df, aes(Individuals, Richness)) +
    geom_line(
      aes(group = Sample, color = Color, linetype = LineType),
      linewidth = 0.8,
      alpha = 0.8
    ) +
    scale_color_identity() +
    scale_linetype_identity() +
    labs(
      title = gsub("_", " ", name),
      x = "Sample Size",
      y = "ASVs"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 20),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save individual plot
  ggsave(
    file.path(qc.dir, paste0(name, "_rarefaction.svg")),
    rare_plot,
    width = 8,
    height = 8,
    dpi = 300
  )
  
  rare_plots[[name]] <- rare_plot  # Store for combined plot
  rare_dim[[name]] <- dim(count_rar)  # Store for number of ASVs
}

# Posterior changes for all plots
# for (name in names(count_list)){
#   rare_plots[[name]] <- rare_plots[[name]] +
#     labs(
#       title = NULL,
#     )
# }

# Combine plots 
combined_plot <- plot_grid(plotlist = rare_plots, ncol = 3,
                           labels = "AUTO", label_size = 22, label_fontface = "bold")
combined_plot

ggsave(file.path(qc.dir, "combined_rarefaction.svg"), combined_plot, width = 24, height = 8, dpi = 300)

# Save number of ASVs in the true samples after each step.
asv_number <- do.call(rbind, rare_dim)
colnames(asv_number) <- c("ASVs", "samples")

write.xlsx(asv_number, file.path(qc.dir, "asv_number.xlsx"))

#### End ####

#### Remove NA Phyla####

rank_names(physeq)

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

physeq_oNA <- subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

physeq
physeq_oNA
physeq = physeq_oNA

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

#### End ####

####Define prevalence of each taxa (in how many samples did each taxa appear at least once)####

prev0 = apply(X = otu_table(physeq),
              MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(physeq),
                    tax_table(physeq))

# Save ASV Prevalence and Abundance table before filtering

write.table(prevdf, file.path(res.dir,"asv_prevdf.tsv"), sep="\t", quote=F, col.names=NA)

# Plot Taxa prevalence v. total counts. Each point is a different taxa. 

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(physeq),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#### End ####

#### Remove taxa not seen more than 3 times in at least 1% of the samples #### 
# This protects against an OTU with small mean & trivially large C.V.
# Setting filter parameters :

countperasv = 3
Samplepercentage = 0.01

physeq_filtered = filter_taxa(physeq, function(x) sum(x > countperasv) > (Samplepercentage*length(x)), TRUE)
physeq_filtered
physeq
physeq = physeq_filtered

#### End ####

#### Normalize number of reads in each sample using median sequencing depth ####

total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_mednorm = transform_sample_counts(physeq, standf)

# Transform to relative abundance. Save as new object.
physeq_re = transform_sample_counts(physeq_mednorm, function(x){x / sum(x)})

#### End ####

#### Exploratory plots after filtering and normalization ####
# Check individual phylum Abundance
# Abundance value transformation function

plot_abundance = function(physeq, ylabn = "",
                          Facet = "Class",
                          Color = "Phylum",
                          n = NULL){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = "Lung", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet, nrow = n) + ylab(ylabn) +
    scale_y_log10()
}

# Plot the abundance values before and after transformation

pl_ab_original  = plot_abundance(physeq,"Original Abundances")
pl_ab_original_norm  =plot_abundance(physeq_mednorm,"Normalized to squencing depth Abundances")
pl_ab_original_norm_re  =plot_abundance(physeq_re,"Normalized Relative Abundances")

grid.arrange(pl_ab_original, pl_ab_original_norm, pl_ab_original_norm_re)

#### End ####

#### Save RData ####

# Save as .RData to load in following steps or continue later

checkpoint5 <- file.path(R.dir, "5.phyloseq.filtered.RData")
if (file.exists(checkpoint5)) {
  load(checkpoint5)
} else {
  save.image(checkpoint5)
}
saveRDS(physeq_oNA, file.path(R.dir, "physeq_oNA.RDS"))
saveRDS(physeq, file.path(R.dir, "physeq_abs_filtered.RDS"))
saveRDS(physeq_re, file.path(R.dir, "physeq_re.RDS"))

#### End ####
