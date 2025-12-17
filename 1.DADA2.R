####Introduction####

# This script infers Amplicon Sequence Variants (ASVs) using the DADA2 pipeline.
# Adapter sequences have been pre-trimmed using Cutadapt.
# Taxonomy analysis combines DADA2 and Blast approaches.
# The inferred ASVs are then converted into phyloseq objects.
# The environment is saved at different stages of the workflow in RData files
# for subsequent analyses (e.g. unfiltered phyloseq objects for alpha diversity analysis; filtered phyloseq objects for taxonomy analyses).

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
library(scales)   
library(gplots) 
library(ggstatsplot) 
theme_set(theme_bw())

# Excel
library(openxlsx)

#### End ####

#### Set Directory ####

path <- "Data/trimmed_cutadapt3/" # Change to the directory containing the fastq files
list.files(path)

# Create directories to save results of the analysis 

res.dir <- "Results"
qc.dir <- file.path(res.dir, "1.QC")
tax.dir <- file.path(res.dir, "2.blast")
R.dir <- file.path(res.dir, "RData")

dir.create(qc.dir, recursive = T)
dir.create(tax.dir, recursive = T)
dir.create(R.dir, recursive = T)

#### End ####

#### Quality Control and Filtering####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq

fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE, recursive = T))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE, recursive = T))

# Extract sample names

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles and save as a pdf

pdf(file.path(qc.dir, "quality_profiles_after_cutadapt.pdf"), width = 8, height = 5)
plotQualityProfile(fnFs[c(5,29,103,141,175,226)])
plotQualityProfile(fnRs[c(5,29,103,141,175,226)])
dev.off()

# Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim 

# We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE.
# The maxEE parameter (default maxEE = 2) sets the maximum number of "expected errors" allowed in a read,
# which is a better filter than simply averaging quality scores, was eased to maxEE=c(2,5).
# Watch out with trunclen, reads have to overlap at the end.

# Define your truncation lengths 
truncLenF <- 230  # Forward trunc length
truncLenR <- 230  # Reverse trunc length

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncLenF,truncLenR),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# Check quality after filtering

pdf(file.path(qc.dir, "quality_profiles_after_dada2_filter.pdf"), width = 8, height = 5)
plotQualityProfile(filtFs[c(5,29,103,141,175,226)])
plotQualityProfile(filtRs[c(5,29,103,141,175,226)])
dev.off()

#### End ####

#### Remove empty samples from list ####

# Since we are dealing with low biomass samples, some samples lost all reads after filtering.
filtFs <- sort(list.files(file.path(path, "filtered"), pattern="_F_filt.fastq", full.names = TRUE, recursive = T))
filtRs <- sort(list.files(file.path(path, "filtered"), pattern="_R_filt.fastq", full.names = TRUE, recursive = T))

# Set sample names as in metadata
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) %>% gsub("-", "_", .) 

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

# Remove chimeras
seqtab.nochim.perbatch.naive <- list()
for(batch in unique(filt$Batch)){
  seqtab.nochim.perbatch.naive[[batch]] <- removeBimeraDenovo(seqtab.perbatch.naive[[batch]], method="consensus", multithread=FALSE, verbose=TRUE)
  message("Batch ", batch, ": ", dim(seqtab.nochim.perbatch.naive[[batch]])[1], " samples and ", dim(seqtab.nochim.perbatch.naive[[batch]])[2], " non chimeric sequences")
}

# Merged all sequence tables into one
seqtab.naive <- seqtab.perbatch.naive[[1]]
seqtab.nochim.naive <- seqtab.nochim.perbatch.naive[[1]]
for (i in 2:length(seqtab.nochim.perbatch.naive)) {
  seqtab.naive <- mergeSequenceTables(seqtab.naive, seqtab.perbatch.naive[[i]])
  seqtab.nochim.naive <- mergeSequenceTables(seqtab.nochim.naive, seqtab.nochim.perbatch.naive[[i]])
}

dim(seqtab.naive)
dim(seqtab.nochim.naive)

# Percentage of not chimeric reads

sum(seqtab.nochim.naive)/sum(seqtab.naive)

#### End ####

#### Second ASV inference with list of prior ASVs, merge and final sequence table ####

# Prepare priors for second pass
merged_seqs <- colnames(seqtab.nochim.naive)

# For forward reads: take first truncLenF bases of merged sequence
priorF <- substr(merged_seqs, 1, truncLenF)

# For reverse reads: take last truncLenR bases of merged sequence and reverse-complement
seq_lengths <- nchar(merged_seqs)
priorR <- substr(merged_seqs, seq_lengths - truncLenR + 1, seq_lengths)
priorR_rc <- as.character(reverseComplement(DNAStringSet(priorR)))

# Rerun dada with priors
dadaFs.prior <- list()
dadaRs.prior <- list()
for(batch in unique(filt$Batch)){
  dadaFs.prior[[batch]] <- dada(filt[filt$Batch == batch,"filtFs"], err=errF[[batch]], priors = priorF, multithread=FALSE)
  dadaRs.prior[[batch]] <- dada(filt[filt$Batch == batch,"filtRs"], err=errR[[batch]], priors = priorR_rc, multithread=FALSE)
}

#Merge paired reads
mergers <- list()

for(batch in unique(filt$Batch)){
  mergers[[batch]] <- mergePairs(dadaFs.prior[[batch]], filt[filt$Batch == batch,"filtFs"],
                                 dadaRs.prior[[batch]], filt[filt$Batch == batch,"filtRs"], verbose=TRUE) # min overlap is 12 as default, but can be adjusted
}

# Obtain sequence tables
seqtab.perbatch <- list()
for(batch in unique(filt$Batch)){
  seqtab.perbatch[[batch]] <- makeSequenceTable(mergers[[batch]])
  message("Batch ", batch, ": ", dim(seqtab.perbatch[[batch]])[1], " samples and ", dim(seqtab.perbatch[[batch]])[2], " sequences")
}

# Remove chimeras
seqtab.nochim.perbatch <- list()
for(batch in unique(filt$Batch)){
  seqtab.nochim.perbatch[[batch]] <- removeBimeraDenovo(seqtab.perbatch[[batch]], method="consensus", multithread=FALSE, verbose=TRUE)
  message("Batch ", batch, ": ", dim(seqtab.nochim.perbatch[[batch]])[1], " samples and ", dim(seqtab.nochim.perbatch[[batch]])[2], " non chimeric sequences")
}

# Merged all sequence tables into one
seqtab <- seqtab.perbatch[[1]]
seqtab.nochim <- seqtab.nochim.perbatch[[1]]
for (i in 2:length(seqtab.nochim.perbatch)) {
  seqtab <- mergeSequenceTables(seqtab, seqtab.perbatch[[i]])
  seqtab.nochim <- mergeSequenceTables(seqtab.nochim, seqtab.nochim.perbatch[[i]])
}

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

rownames(out) <- gsub("_.*", "", rownames(out))
rownames(track) <- gsub("_F_filt.fastq", "", rownames(track))
rownames(seqtab.nochim) <- gsub("_F_filt.fastq", "", rownames(seqtab.nochim))
rownames(seqtab.nochim.naive) <- gsub("_F_filt.fastq", "", rownames(seqtab.nochim.naive))
track <- merge(as.data.frame(out), as.data.frame(track), by = 0, all.x = T) %>% merge(.,rowSums(seqtab.nochim.naive), by.x = "Row.names", by.y = 0) %>% merge(.,rowSums(seqtab.nochim), by.x = "Row.names", by.y = 0)

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("sample","input", "filtered", "denoisedF.naive", "denoisedR.naive", "denoisedF.prior", "denoisedR.prior", "merged.naive", "merged.prior", "nonchim.naive", "nonchim.prior")

head(track)

# Save the track

write.csv(track, file = file.path(qc.dir, "Pipeline_Track.csv"))

#### End ####

#### Save RData ####

# Save as .RData to load in following steps or continue later
save.image(file.path(R.dir,"1.dada2_track.RData"))
#load(file.path(R.dir,"1.dada2_track.RData"))

#### End ####

#### Assign taxonomy for Bacteria with ####

# Function
assign_taxonomy <- function(seqtab, ref_fasta, ref_fasta_species){
  set.seed(123456789)
  taxa <- assignTaxonomy(seqtab, ref_fasta, tryRC = TRUE) #tryRC = TRUE -> reverse-complement orientation
  taxa_species = addSpecies(taxa, ref_fasta_species, allowMultiple=TRUE)
  return(list(taxa = taxa, taxa_species = taxa_species))
}

#### GTDB
ref_fasta_GTDB = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz"
ref_fasta_species_GTDB = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc53_ssu_r220_species.fa.gz"

taxa_GTDB <- assign_taxonomy(seqtab.nochim, ref_fasta_GTDB, ref_fasta_species_GTDB)
taxa_species_GTDB = taxa_GTDB[["taxa_species"]]

head(unname(taxa_GTDB[["taxa"]])) # Removing sequence rownames for display only
head(unname(taxa_species_GTDB))

#### SILVA

ref_fasta_SILVA = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/SILVA_2024/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
ref_fasta_species_SILVA = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/SILVA_2024/silva_v138.2_assignSpecies.fa.gz"

taxa_SILVA <- assign_taxonomy(seqtab.nochim, ref_fasta_SILVA, ref_fasta_species_SILVA)
taxa_species_SILVA = taxa_SILVA[["taxa_species"]]

head(unname(taxa_SILVA[["taxa"]])) # Removing sequence rownames for display only
head(unname(taxa_species_SILVA))

#### eHOMD

ref_fasta_eHOMD = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/eHOMD/eHOMD_RefSeq_dada2_V15.22.fasta.gz"
ref_fasta_species_eHOMD = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/eHOMD/eHOMD_RefSeq_dada2_assign_species_V15.22.fasta.gz"

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
colnames(count_tab.naive) <- sub("_F_filt.fastq", "", colnames(count_tab.naive)) %>% gsub("-", "_", .) #%>%
# gsub("2C", "2_Cul", .) %>% gsub("1C", "1_Cul", .)
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
colnames(count_tab) <- sub("_F_filt.fastq", "", colnames(count_tab)) %>% gsub("-", "_", .) #%>%
#gsub("2C", "2_Cul", .) %>% gsub("1C", "1_Cul", .)
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

save.image(file.path(R.dir, "2.dada2_taxonomy.RData"))
#load(file.path(R.dir, "2.dada2_taxonomy.RData"))

#### End ####

#### Blast unassigned species with script 1.1 Blast and 1.2 Filter_blast ####

tax_tab <- read.xlsx(file.path(tax.dir,"tax_tab_SILVA_blast_species.xlsx"), rowNames = T)
tax_tab <- tax_tab[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

# Paste Genus and Species name for later visualization
colnames(tax_tab)[colnames(tax_tab) == "Species"] <- "Species.only"
tax_tab$Species <- NA  # Initialize Species column with NA

tax_tab$Species <- ifelse(!is.na(tax_tab$Genus) & !is.na(tax_tab$Species.only),
                          paste(tax_tab$Genus, tax_tab$Species.only),
                          ifelse(!is.na(tax_tab$Genus) & is.na(tax_tab$Species.only),
                                 paste(tax_tab$Genus, "spp."),
                                 NA))
tax_tab$Species.only <- NULL

#### End ####

#### Add samples metadata and match to count table####

sample_info <- read.xlsx("Data/000-LML_Metadata_clinical_Lab_final_20240114.xlsx", rowNames = T, sheet = "Combined_metadata")
sample_info <- sample_info[is.na(sample_info$Isolation) | sample_info$Isolation == "Direct_Isolation",]

## Decode sample_info data and convert to factors

# Save code in new variables
variables <- c("Sex", "History.of.smoking", "Active.smoker", "Side", "Lobe", "Synchronous.tumor")
variables.code <- paste0(variables, ".code")
for (i in seq_along(variables)) {
  sample_info[[variables.code[i]]] <- sample_info[[variables[i]]]
}

# Convert to factors and decode
#sample_info$Isolation <- factor(sample_info$Isolation, levels = c("Direct_Isolation", "Culture_Enriched"))
sample_info$Lung <- factor(sample_info$Lung, levels = c("Normal", "Diseased"))
sample_info$Histology.NSCLC <- factor(sample_info$Histology.NSCLC, levels = c("Adenocarcinoma", "Squamous cell carcinoma", "Adeno-Squamous cell carcinoma", "Others"))
sample_info$Diagnosis <- factor(sample_info$Diagnosis, levels = c("Benign", "NSCLC", "V.a. NSCLC", "SCLC", "Atypical Carzinoid", "Carcinoid Tumorlet", "Metastasized Tumor"))
sample_info$T <- factor(sample_info$T)
sample_info$N <- factor(sample_info$N)
sample_info$M <- factor(sample_info$M)

sample_info <- sample_info %>%
  mutate(Sex = case_when(
    Sex == 0 ~ "Female",
    Sex == 1 ~ "Male"))
sample_info$Sex <- factor(sample_info$Sex)

sample_info <- sample_info %>%
  mutate(History.of.smoking = case_when(
    History.of.smoking == 0 ~ "None",
    History.of.smoking == 1 ~ "Cigarettes",
    History.of.smoking == 2 ~ "Pipe"))
sample_info$History.of.smoking <- factor(sample_info$History.of.smoking, levels = c("None", "Cigarettes", "Pipe"))

sample_info$History.of.smoking.y.n <- sample_info$History.of.smoking.code
sample_info <- sample_info %>%
  mutate(History.of.smoking.y.n = case_when(
    History.of.smoking.y.n == 0 ~ "No",
    History.of.smoking.y.n == 1 ~ "Yes",
    History.of.smoking.y.n == 2 ~ "Yes"))
sample_info$History.of.smoking.y.n <- factor(sample_info$History.of.smoking.y.n, levels = c("No", "Yes"))

sample_info <- sample_info %>%
  mutate(Active.smoker = case_when(
    Active.smoker == 0 ~ "No",
    Active.smoker == 1 ~ "Yes"))
sample_info$Active.smoker <- factor(sample_info$Active.smoker, levels = c("No", "Yes"))

sample_info <- sample_info %>%
  mutate(Side = case_when(
    Side == 0 ~ "Right",
    Side == 1 ~ "Left"))
sample_info$Side <- factor(sample_info$Side)

sample_info <- sample_info %>%
  mutate(Lobe = case_when(
    Lobe == 1 ~ "Upper",
    Lobe == 2 ~ "Lower",
    Lobe == 3 ~ "Middle",
    Lobe == 4 ~ "Lingula"))
sample_info$Lobe <- factor(sample_info$Lobe, levels = c("Upper", "Lower", "Middle", "Lingula"))

sample_info <- sample_info %>%
  mutate(Synchronous.tumor = case_when(
    Synchronous.tumor == 0 ~ "No",
    Synchronous.tumor == 1 ~ "Second tumor contralateral",
    Synchronous.tumor == 2 ~ "Thoracic metastasis of the primary tumor"))
sample_info$Synchronous.tumor <- factor(sample_info$Synchronous.tumor)

# Select colors for each variable
palette2 <- rev(scales::hue_pal()(2))
palette3 <-  c("#00BFC4", "#F8766D", "#00BA38")
palette <- c("#53B400", "#F8766D", "#00C094", "#344cb7",  "#00B6EB" ,"#A58AFF", "#FB61D7") 

isolation_col <- setNames(palette2[1:length(levels(sample_info$Isolation))], levels(sample_info$Isolation))
lung_col <- setNames(palette2[1:length(levels(sample_info$Lung))], levels(sample_info$Lung))
diagnosis_col <- setNames(palette[1:length(levels(sample_info$Diagnosis))], levels(sample_info$Diagnosis))
smokerh_col <- setNames(palette3[1:length(levels(sample_info$History.of.smoking))], levels(sample_info$History.of.smoking))
smoker_col <- setNames(palette3[1:length(levels(sample_info$Active.smoker))], levels(sample_info$Active.smoker))
histology_col <- setNames(palette[1:length(levels(sample_info$Histology.NSCLC))], levels(sample_info$Histology.NSCLC))
tumor_col <- setNames(c("#53B400", "#F8766D", "#00B6EB" ,"#A58AFF", "#FB61D7"), 0:4)

# Examine all tables

head(count_tab)
head(tax_tab)
head(fasta_tab)
head(sample_info)

# Examine consistancy in order between count_tab colnames and coldata rownames 
# (They have to be in the same order or Deseq2 won't work out)

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

#-----> Fasttree in ubuntu (1.3.Fasttree)
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
names(dna) <- taxa_names(physeq)

ph_tree = phy_tree(tree)            # Making Phylogenetic tree
taxa_names(ph_tree) = taxa_names(dna)

physeq_dna_tree <- merge_phyloseq(physeq, ph_tree, dna) #Merging  ASV Fasta sequences and Phylogenetic tree to phyloseq object

taxa_names(physeq_dna_tree) <- paste0("ASV", seq(ntaxa(physeq_dna_tree)))

physeq
physeq_dna_tree
physeq = physeq_dna_tree

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image(file.path(R.dir,"3.physeq.original.RData"))
#load(file.path(R.dir,"3.physeq.original.RData"))

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

sample_data(ps.decontam)$Batch_pooled <- sample_data(ps.decontam)$Batch
sample_data(ps.decontam)$Batch_pooled[sample_data(ps.decontam)$Batch_pooled %in% c("1-LML_001-006", "2-LML_007-030", "3-LML_031-054", "4-LML_055-077")] <- "1-4 batch"
sample_data(ps.decontam)$is.neg <- sample_data(ps.decontam)$Sample_or_Control == "Negative control sample"

#' The probability threshold below which the null-hypothesis (not a contaminant) should be rejected in favor of
#' the alternate hypothesis (contaminant) was set to 0.4 (default 0.1) 
thr = 0.5
contamdf.prev <- isContaminant(ps.decontam, method="prevalence", neg="is.neg", threshold=thr, batch = "Batch_pooled")
table(contamdf.prev$contaminant)

hist(contamdf.prev$p, 100, main="P-values for different thresholds in decontam")

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
samp_ex = sample_data(physeq)[sample_data(physeq)$Sex.code %in%
                                c("excluded due to history of chemotherapy", "excluded due to history of radiochemotherapy") |
                                sample_data(physeq)$Study_Nr == "LML_112",]
samp_ex
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

ggsave(file.path(qc.dir, "LibSizeLog10_direct.tiff"), width = 16, height = 10, dpi = 300)

# Bar plots library size - Diagnosis relation
ggplot(data=df, aes(x=Lung,  y = count, color=LibrarySizeBreaks, fill = LibrarySizeBreaks)) +
  geom_bar(position="fill", stat="identity") +
  #ggtitle("Library size in directly isolated samples") +
  facet_wrap(.~Diagnosis) + theme(text = element_text(size = 30))
ggsave(file.path(qc.dir, "LibSizeProportion_direct.tiff"), width = 15, height = 12, dpi = 300)
ggplot(data=df, aes(x=Lung,  y = count, color=LibrarySizeBreaks, fill = LibrarySizeBreaks)) +
  geom_bar(stat="identity") +
  #ggtitle("Library size in directly isolated samples") +
  facet_wrap(.~Diagnosis)

ggplot(data=df, aes(x=Lung,  y = count, color=Diagnosis, fill = Diagnosis)) +
  geom_bar(stat="identity") +
  #ggtitle("Library size in directly isolated samples") +
  facet_wrap(.~LibrarySizeBreaks)

ggplot(data=df, aes(x=Lung,  y = count, color=Diagnosis, fill = Diagnosis)) +
  geom_bar(position="fill", stat="identity") +
  #ggtitle("Library size in directly isolated samples") +
  facet_wrap(.~LibrarySizeBreaks)

# Explore correlation of sequencing depth and dna concentration

ggscatterstats(
  data = data.frame(df),
  x = DNA_concentration, 
  y = LibrarySize, 
  xlab = "DNA concentration",
  ylab = "Sequencing depth",
  point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
  # xfill = "#CC79A7", ## fill for marginals on the x-axis
  # yfill = "#009E73", ## fill for marginals on the y-axis
  marginal = F,
  type ="n",
  #title = "Directly isolated samples",
  #results.subtitle = F
) + theme(text = element_text(size = 30),
          plot.title = element_text(size = 28),
          plot.subtitle = element_text(size = 18))
ggsave(file.path(qc.dir, "LibSizeDNA_direct.tiff"), width = 11, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "LibSizeDNA_direct.svg"), width = 11, height = 10, dpi = 300)

#### End ####

#### Calculate and add alpha diversity to phyloseq object for later analysis ####

sample_data(physeq) <- estimate_richness(physeq, measures = c("Observed","Chao1", "Shannon","Simpson","InvSimpson","ACE")) %>%
  merge(sample_data(physeq), ., by = 0) %>% tibble::column_to_rownames(var = "Row.names")
sample_data(physeq)$InvSimpson[sample_data(physeq)$InvSimpson == Inf] = 0

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

#save.image(file.path(R.dir,"4.physeq.decontam.RData"))
load(file.path(R.dir,"4.physeq.decontam.RData"))

#### End ####

#### Rarefaction curves ####

# Define aesthetics
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")

# Create named list of count tables
count_list <- list(
  ASVs_from_naive_inference = count_tab.naive,
  ASVs_after_inference_with_list_of_priors = count_tab,
  ASVs_after_decontamination = data.frame(otu_table(physeq_decont))  # Fixed assignment
)

rare_plots <- list()  # Store ggplot objects

for (name in names(count_list)) {
  count_rar <- count_list[[name]]
  
  # Exclude controls
  count_rar <- count_rar[,colnames(count_rar) %in% rownames(sample_info[sample_info$Sample_or_Control == "True sample",])]
  # Exclude samples with 0 counts
  count_rar <- count_rar[, colSums(count_rar) != 0]
  
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

############ Other options
# Select batch
unique(sample_info$Batch)
batch = "6-LML_101-122"
count_rar <- count_rar[,colnames(count_rar) %in% rownames(sample_info[sample_info$Batch == batch,])]

raremax <- min(rowSums(t(count_rar)))
raremax
rarefy(t(count_rar), raremax)

# Option 1
rarecurve(t(count_rar), step = 100,  col = col, lwd=2, lty = lty, ylab = "ASVs", label = T)
abline(v=(raremax))
rarecurve(t(count_rar), step = 100,  col = col, lwd=2, lty = lty, ylab = "ASVs", label = F)
abline(v=(raremax))

# Option 2
rarecurve(t(count_rar), step = 100, sample = raremax, col = col, lwd=2, lty = lty, ylab = "ASVs", label = T)

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

save.image(file.path(R.dir,"5.phyloseq.filtered.RData"))
#load(file.path(R.dir,"5.phyloseq.filtered.RData"))

#### End ####
