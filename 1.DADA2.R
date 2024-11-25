####Introduction####

#### End ####

#### Package setup ####

Sys.setenv(language = "EN")
library(dada2); packageVersion("dada2")
library(vegan)
library(phyloseq)
library(msa)
library(ape)
library(phangorn)
library(microViz)
library(decontam)
library(ggplot2)
library(gridExtra)
library(openxlsx)
library(dplyr)
theme_set(theme_bw())

#### End ####

#### Set Directory ####

path <- "Data/trimmed_cutadapt3/Samples_NC_PC" # Change to the directory containing the fastq files
list.files(path)

# Create directories to save results of the analysis 

res.dir <- "Results_trimmed_cutadapt3_20241118"
qc.dir <- file.path(res.dir, "1.QC")
R.dir <- file.path(res.dir, "RData")

dir.create(qc.dir, recursive = T)
dir.create(R.dir, recursive = T)

#### End ####

#### Quality Control and Filtering####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE, recursive = T))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE, recursive = T))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles

plotQualityProfile(fnFs[c(4,28,123,219,311,403)])

plotQualityProfile(fnRs[c(4,28,123,219,311,403)])

# Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim 

# We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE
# and  maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read,
# which is a better filter than simply averaging quality scores.
# Watch out with trunclen, reads have to overlap at the end, you have to try out.
# maxEE can be eased maxEE=c(2,5) if too many read are lost because of low quality.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# Check quality after filtering

plotQualityProfile(filtFs[c(4,28,123,219,311,403)])

plotQualityProfile(filtRs[c(4,28,123,219,311,403)])

#### End ####

#### Remove empty samples from list ####
filtFs <- sort(list.files(file.path(path, "filtered"), pattern="_F_filt.fastq", full.names = TRUE, recursive = T))
filtRs <- sort(list.files(file.path(path, "filtered"), pattern="_R_filt.fastq", full.names = TRUE, recursive = T))

# Set sample names as in metadata
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.names <- sub("_F_filt.fastq", "", sample.names) %>% gsub("-", "_", .) %>%
  gsub("2C", "2_Cul", .) %>% gsub("1C", "1_Cul", .)

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#### End ####

#### Learn the Error Rates ####
# In the BigData WF (https://benjjneb.github.io/dada2/bigdata.html),
# it is recommended to learn the error rates and do the sample inference for each run individually.
# Thus we add already the metadata to know the batch of each sample

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

# Merge paired reads and obtain sequence table - just to compare with priors - can be skipped

mergers.naive <- list()

for(batch in unique(filt$Batch)){
mergers.naive[[batch]] <- mergePairs(dadaFs[[batch]], filt[filt$Batch == batch,"filtFs"],
                               dadaRs[[batch]], filt[filt$Batch == batch,"filtRs"], verbose=TRUE) # min overlap is 12 as default, but can be adjusted
}

# Obtain sequence tables
seqtab.naive <- list()
for(batch in unique(filt$Batch)){
  seqtab.naive[[batch]] <- makeSequenceTable(mergers.naive[[batch]])
  print(batch)
  print(dim(seqtab.naive[[batch]]))
}

# Merged all sequence tables into one
seqtab.all.naive <- seqtab.naive[[1]]
for (i in 2:length(seqtab.naive)) {
  seqtab.all.naive <- mergeSequenceTables(seqtab.all.naive, seqtab.naive[[i]])
}

seqtab.nochim.naive <- removeBimeraDenovo(seqtab.all.naive, method="consensus", multithread=FALSE, verbose=TRUE)

# Construct sequence table from forward and reverse for priors
seqtabFs <- list()
seqtabRs <- list()
for(batch in unique(filt$Batch)){
  seqtabFs[[batch]] <- makeSequenceTable(dadaFs[[batch]])
  seqtabRs[[batch]] <- makeSequenceTable(dadaRs[[batch]])
  print(batch)
  print(dim(seqtabFs[[batch]]))
  print(dim(seqtabRs[[batch]]))
}

# Merge sequence tables 
seqtabFs.all <- seqtabFs[[1]]
for (i in 2:length(seqtabFs)) {
  seqtabFs.all <- mergeSequenceTables(seqtabFs.all, seqtabFs[[i]])
}
seqtabRs.all <- seqtabRs[[1]]
for (i in 2:length(seqtabRs)) {
  seqtabRs.all <- mergeSequenceTables(seqtabRs.all, seqtabRs[[i]])
}

#### End ####

#### Second ASV inference with list of prior ASVs, merge and final sequence table ####

# Extract list of prior ASVs
priorF <- colnames(seqtabFs.all)
priorR <- colnames(seqtabRs.all)

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

#Most of your reads should successfully merge. If that is not the case upstream parameters may need
#to be revisited: Did you trim away the overlap between your reads?

seqtab <- list()
for(batch in unique(filt$Batch)){
  seqtab[[batch]] <- makeSequenceTable(mergers[[batch]])
  print(batch)
  print(dim(seqtab[[batch]]))
}

# Merged all sequence tables into one
seqtab.all <- seqtab[[1]]
for (i in 2:length(seqtab)) {
  seqtab.all <- mergeSequenceTables(seqtab.all, seqtab[[i]])
}

# Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=FALSE, verbose=TRUE)

dim(seqtab.nochim)

# Percentage of not chimeric reads

sum(seqtab.nochim)/sum(seqtab.all)

#### End ####

#### Track Pipeline ####
#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))

track <- list()
for(batch in unique(filt$Batch)){
track[[batch]] <- cbind(sapply(dadaFs[[batch]], getN), sapply(dadaRs[[batch]], getN), sapply(mergers.naive[[batch]], getN), sapply(dadaFs.prior[[batch]], getN), sapply(dadaRs.prior[[batch]], getN), sapply(mergers[[batch]], getN))
}
track <- as.data.frame(do.call(rbind, track))

rownames(out) <- gsub("_.*", "", rownames(out))
rownames(track) <- gsub("_F_filt.fastq", "", rownames(track))
rownames(seqtab.nochim) <- gsub("_F_filt.fastq", "", rownames(seqtab.nochim))
rownames(seqtab.nochim.naive) <- gsub("_F_filt.fastq", "", rownames(seqtab.nochim.naive))
track <- merge(as.data.frame(out), as.data.frame(track), by = 0, all.x = T) %>% merge(.,rowSums(seqtab.nochim), by.x = "Row.names", by.y = 0) %>% merge(.,rowSums(seqtab.nochim.naive), by.x = "Row.names", by.y = 0)

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("sample","input", "filtered", "denoisedF.naive", "denoisedR.naive", "merged.naive", "denoisedF.prior", "denoisedR.prior", "merged.prior", "nonchim.prior", "nonchim.naive")


head(track)

# Save the track

write.csv(track, file = file.path(qc.dir, "Pipeline_Track.csv"))

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image(file.path(R.dir,"1.dada2_track.RData"))
#load(file.path(R.dir,"1.dada2_track.RData"))

#### End ####

#### Assign taxonomy for Bacteria with GTDB####

#ref_fasta = "C:/Users/posadas/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc53_ssu_r214_fullTaxo.fa.gz"
#ref_fasta = "C:/Users/posadas/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/GTDB/old/GTDB_bac120_arc53_ssu_r207_fullTaxo.fa.gz"
#ref_fasta = "C:/Users/posadas/Desktop/GitHub/VM/Databases/16S_ribosomal_RNA_blast_dada.format.fasta"
#ref_fasta = "C:/Users/posadas/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/SILVA_2020/silva_nr99_v138_train_set.fa.gz"
ref_fasta = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/SILVA_2024/silva_nr99_v138.2_toSpecies_trainset.fa.gz"

#tryRC = TRUE -> reverse-complement orientation
set.seed(123456789)
taxa <- assignTaxonomy(seqtab.nochim, ref_fasta,tryRC = TRUE )
unname(taxa)

taxa.print  <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#ref_fasta_species = "C:/Users/posadas/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc53_ssu_r214_species.fa.gz"
#ref_fasta_species = "C:/Users/posadas/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/GTDB/old/GTDB_bac120_arc53_ssu_r207_Species.fa.gz"
#ref_fasta_species = "C:/Users/posadas/Desktop/GitHub/VM/Databases/16S_ribosomal_RNA_blast_dada.format.species.fasta"
#ref_fasta_species = "C:/Users/posadas/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/SILVA_2020/silva_species_assignment_v138.fa.gz"
ref_fasta_species = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/SILVA_2024/silva_v138.2_assignSpecies.fa.gz"

taxa_species = addSpecies(taxa,ref_fasta_species, allowMultiple=TRUE)
taxa.print_spp  <- taxa_species  # Removing sequence rownames for display only
rownames(taxa.print_spp) <- NULL
head(taxa.print_spp)

length(taxa.print_spp[,"Species"][!is.na(taxa.print_spp[,"Species"])])
length(taxa.print_spp[,"Genus"][!is.na(taxa.print_spp[,"Genus"])])
length(taxa.print[,"Species"][!is.na(taxa.print[,"Species"])])
#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image(file.path(R.dir, "2.dada2_taxonomy.RData"))
#load("Results/RData/2.dada2_taxonomy.RData")

#### End ####

#### Make Count and Taxa Tables from naive dada ####

# Giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs.n <- colnames(seqtab.nochim.naive)

asv_headers.n <- vector(dim(seqtab.nochim.naive)[2], mode="character")

for (i in 1:dim(seqtab.nochim.naive)[2]) {
  asv_headers.n[i] <- paste(">ASV", i, sep="_")
}

# Making fasta of our final ASV seqs:

fasta_tab.naive <- c(rbind(asv_headers.n, asv_seqs.n))
# write(fasta_tab.naive, file.path(res.dir, "fasta_tab_naive.fa"))

# count table:

count_tab.naive <- t(seqtab.nochim.naive)
row.names(count_tab.naive) <- sub(">", "", asv_headers.n)
colnames(count_tab.naive) <- sub("_F_filt.fastq", "", colnames(count_tab.naive)) %>% gsub("-", "_", .) %>%
  gsub("2C", "2_Cul", .) %>% gsub("1C", "1_Cul", .)
# write.table(count_tab.naive, file.path(res.dir, "count_tab_naive.tsv"), sep="\t", quote=F, col.names=NA)

# tax table:

#tax_tab <- taxa_species
# row.names(tax_tab) <- sub(">", "", asv_headers)
# write.table(tax_tab, "Results/tax_tab_GTDB.tsv", sep="\t", quote=F, col.names=NA)

#### End ####

#### Make Count and Taxa Tables ####

# Giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochim)

asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Making fasta of our final ASV seqs:

fasta_tab <- c(rbind(asv_headers, asv_seqs))
#write(fasta_tab, file.path(res.dir, "fasta_tab.fa"))

# count table:

count_tab <- t(seqtab.nochim)
row.names(count_tab) <- sub(">", "", asv_headers)
colnames(count_tab) <- sub("_F_filt.fastq", "", colnames(count_tab)) %>% gsub("-", "_", .) %>%
  gsub("2C", "2_Cul", .) %>% gsub("1C", "1_Cul", .)
#write.table(count_tab, file.path(res.dir, "count_tab.tsv"), sep="\t", quote=F, col.names=NA)

# tax table:

# tax_tab <- taxa_species
# row.names(tax_tab) <- sub(">", "", asv_headers)
#write.xlsx(as.data.frame(tax_tab), file.path(res.dir,"2.blast/tax_tab_SILVA.xlsx"), rowNames = T)

#### Blast unassigned species with script 1.2 Blast and 1.3 Filter_blast
tax_tab <- read.xlsx(file.path(res.dir,"2.blast/tax_tab_SILVA_blast_species.xlsx"), rowNames = T)

#### End ####

#### Add samples metadata and match to count table####

#sample_info <- read.xlsx("Data/000-LML_Metadata_clinical_Lab_final_20240114.xlsx", rowNames = T, sheet = "Combined_metadata")

## Decode sample_info data and convert to factors

# Save code in new variables
variables <- c("Sex", "History.of.smoking", "Active.smoker", "Side", "Lobe", "Synchronous.tumor")
variables.code <- paste0(variables, ".code")
for (i in seq_along(variables)) {
  sample_info[[variables.code[i]]] <- sample_info[[variables[i]]]
}

# Convert to factors and decode
sample_info$Isolation <- factor(sample_info$Isolation, levels = c("Direct_Isolation", "Culture_Enriched"))
sample_info$Lung <- factor(sample_info$Lung, levels = c("Normal", "Diseased"))
sample_info$Histology.NSCLC <- factor(sample_info$Histology.NSCLC, levels = c("Adenocarcinoma", "Squamous cell carcinoma", "Adeno-Squamous cell carcinoma", "Others"))
sample_info$Diagnosis <- factor(sample_info$diagnosis, levels = c("Benign", "NSCLC", "V.a. NSCLC", "SCLC", "Atypical Carzinoid", "Carcinoid Tumorlet", "Metastasized Tumor"))

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
palette <- c("#53B400", "#F8766D", "#00C094", "#C49A00",  "#00B6EB" ,"#A58AFF", "#FB61D7") 

isolation_col <- setNames(palette2[1:length(levels(sample_info$Isolation))], levels(sample_info$Isolation))
lung_col <- setNames(palette2[1:length(levels(sample_info$Lung))], levels(sample_info$Lung))
diagnosis_col <- setNames(palette[1:length(levels(sample_info$Diagnosis))], levels(sample_info$Diagnosis))
smokerh_col <- setNames(palette3[1:length(levels(sample_info$History.of.smoking))], levels(sample_info$History.of.smoking))
smoker_col <- setNames(palette3[1:length(levels(sample_info$Active.smoker))], levels(sample_info$Active.smoker))
histology_col <- setNames(palette[1:length(levels(sample_info$Histology.NSCLC))], levels(sample_info$Histology.NSCLC))

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

#### Rarefaction curves ####

col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")

# From decontam physeq
count_tab_decontam = data.frame(otu_table(physeq))

# Exclude controls

count_rar <- count_tab_decontam[,!colnames(count_tab_decontam) %in% rownames(sample_info[sample_info$Sample_or_Control %in% c("Positive control sample", "Negative control sample"),])]

count_rar <- count_tab.naive[,!colnames(count_tab.naive) %in% rownames(sample_info[sample_info$Sample_or_Control %in% c("Positive control sample", "Negative control sample"),])]

count_rar <- count_tab[,!colnames(count_tab) %in% rownames(sample_info[sample_info$Sample_or_Control %in% c("Positive control sample", "Negative control sample"),])]

#Exclude samples with 0 counts
#count_rar <- count_rar[, colSums(count_rar) != 0]

# Exclude culture samples

count_rar <- count_rar[,-grep("_Cul", colnames(count_rar))]

count_rar <- count_rar[,-grep("013|035|112", colnames(count_rar))]

# Select batch
unique(sample_info$Batch)
Rarefaction_curves_Batch6_prior_nolab
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

#### Phylogenetic tree ####

 seqs <- getSequences(seqtab.nochim)
 names(seqs) <- seqs # This propagates to the tip labels of the tree
 mult <- msa(seqs, method="ClustalW", type="dna", order="input")
 
 aligned_asvs <- as(mult, "DNAStringSet")
 names(aligned_asvs) <- seqs
 writeXStringSet(aligned_asvs, filepath = file.path(res.dir, "aligned_asvs.fasta"), format = "fasta")

 #-----> Fasttree in ubuntu
 tree <- read.tree(file.path(res.dir, "tree.nwk"))
 
 # The phangorn package is then used to construct a phylogenetic tree. 
 # Here we first construct a neighbor-joining tree, and then fit a GTR+G+I maximum
 # likelihood tree using the neighbor-joining tree as a starting point.
 
# phang.align <- as.phyDat(mult, type="dna", names=getSequence(seqtab.nochim))
# dm <- dist.ml(phang.align)
# treeNJ <- NJ(dm) # Note, tip order != sequence order
# fit = pml(treeNJ, data=phang.align)
# 
# # Negative edges length changed to 0!
# 
# fitGTR <- update(fit, k=4, inv=0.2)
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "stochastic", control = pml.control(trace = 0))
# 
# detach("package:phangorn", unload=TRUE)

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

#### Remove samples that don´t fit the study conditions ####
samp_ex = sample_data(physeq)[sample_data(physeq)$Sex.code %in%
                      c("excluded due to history of chemotherapy", "excluded due to history of radiochemotherapy") |
                        sample_data(physeq)$Study_Nr == "LML_112",]
samp_ex
physeq_study = prune_samples(!sample_names(physeq) %in% rownames(samp_ex), physeq)

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

# Explore library sizes per diagnosis
breaks <- c(0, 1, 10, 100, 1000, 10000, 50000, Inf)
labels <- c("0", "1-10", "11-100", "101-1000", "1001-10000", "10001-50000", ">50000")
df$LibrarySizeBreaks <- cut(df$LibrarySize, breaks = breaks, labels = labels, include.lowest = TRUE)
df$count = 1
df$LibrarySizeLog <- log(df$LibrarySize)
df <- df[df$Sample_or_Control == "True sample",]

direct <- df[df$Isolation == "Direct_Isolation" & df$Study_Nr != "LML_112",]
cult <- df[df$Isolation == "Culture_Enriched" & df$Study_Nr != "LML_112",]

# Box plots
ggplot(data=direct, aes(x=Lung, y=LibrarySizeLog, color=Lung)) +
  geom_boxplot() +
  geom_point() +
  ggtitle("Library size in directly isolated samples") +
  facet_wrap(.~diagnosis) + theme(text = element_text(size = 30))

ggsave(file.path(qc.dir, "LibSizeLog_direct.tiff"), width = 13, height = 10, dpi = 300)

ggplot(data=cult, aes(x=Lung, y=LibrarySizeLog, color=Lung)) +
  geom_boxplot() +
  geom_point() +
  ggtitle("Library size in culture enriched samples") +
  facet_wrap(.~diagnosis) + theme(text = element_text(size = 30))

ggsave(file.path(qc.dir, "LibSizeLog_cult.tiff"), width = 14, height = 10, dpi = 300)

# Bar plots library size - diagnosis relation
ggplot(data=direct, aes(x=Lung,  y = count, color=LibrarySizeBreaks, fill = LibrarySizeBreaks)) +
  geom_bar(position="fill", stat="identity") +
  ggtitle("Library size in directly isolated samples") +
  facet_wrap(.~diagnosis)
ggplot(data=direct, aes(x=Lung,  y = count, color=LibrarySizeBreaks, fill = LibrarySizeBreaks)) +
  geom_bar(stat="identity") +
  ggtitle("Library size in directly isolated samples") +
  facet_wrap(.~diagnosis)

ggplot(data=cult, aes(x=Lung,  y = count, color=LibrarySizeBreaks, fill = LibrarySizeBreaks)) +
  geom_bar(position="fill", stat="identity") +
  ggtitle("Library size in culture enriched samples") +
  facet_wrap(.~diagnosis)
ggplot(data=cult, aes(x=Lung,  y = count, color=LibrarySizeBreaks, fill = LibrarySizeBreaks)) +
  geom_bar(stat="identity") +
  ggtitle("Library size in culture enriched samples") +
  facet_wrap(.~diagnosis)

ggplot(data=direct, aes(x=Lung,  y = count, color=diagnosis, fill = diagnosis)) +
  geom_bar(stat="identity") +
  ggtitle("Library size in directly isolated samples") +
  facet_wrap(.~LibrarySizeBreaks)
ggplot(data=cult, aes(x=Lung,  y = count, color=diagnosis, fill = diagnosis)) +
  geom_bar(stat="identity") +
  ggtitle("Library size in culture enriched samples") +
  facet_wrap(.~LibrarySizeBreaks)

ggplot(data=direct, aes(x=Lung,  y = count, color=diagnosis, fill = diagnosis)) +
  geom_bar(position="fill", stat="identity") +
  ggtitle("Library size in directly isolated samples") +
  facet_wrap(.~LibrarySizeBreaks)
ggplot(data=cult, aes(x=Lung,  y = count, color=diagnosis, fill = diagnosis)) +
  geom_bar(position="fill", stat="identity") +
  ggtitle("Library size in culture enriched samples") +
  facet_wrap(.~LibrarySizeBreaks)

# Explore correlation of sequencing depth and dna concentration
library(ggstatsplot)
ggscatterstats(
  data = data.frame(direct),
  x = DNA_concentration, 
  y = LibrarySize, 
  xlab = "DNA concentration",
  ylab = "Sequencing depth",
  point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
  xfill = "#CC79A7", ## fill for marginals on the x-axis
  yfill = "#009E73", ## fill for marginals on the y-axis
  type ="n",
  title = "Directly isolated samples",
  #results.subtitle = F
) + theme(text = element_text(size = 16),
          plot.title = element_text(size = 18),
          plot.subtitle = element_text(size = 14))
ggscatterstats(
  data = cult,
  x = DNA_concentration, 
  y = LibrarySize, 
  xlab = "DNA concentration",
  ylab = "Sequencing depth",
  point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
  xfill = "#CC79A7", ## fill for marginals on the x-axis
  yfill = "#009E73", ## fill for marginals on the y-axis
  type ="n",
  title = "Cultured enriched samples",
  #results.subtitle = F
) + theme(text = element_text(size = 16),
          plot.title = element_text(size = 18),
          plot.subtitle = element_text(size = 14))

#### End ####

## Phyloseq Object Filtering ## First filtering step for low count samples and NA phyla ##

#### Remove samples with less than 50 total reads ####
sample_sums(physeq) #Nr of reads per Sample

#physeq_above50reads = prune_samples(sample_sums(physeq)>=50, physeq)

#physeq = physeq_above100reads

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

#save.image(file.path(R.dir,"3.physeq.original_2024.10.RData"))
load(file.path(R.dir,"3.physeq.original.RData"))

#### End ####

## Filter low Abundance Taxa and count table normalization ##

#### Remove ASVs appearing in negative controls ####
# Select negative controls and delete positive ones

PC <- rownames(sample_info[sample_info$Sample_or_Control == "Positive control sample",])
NC <- rownames(sample_info[sample_info$Sample_or_Control == "Negative control sample",])

ps.decontam <- prune_samples(!sample_names(physeq) %in% PC, physeq)
ps.decontam <- prune_samples(sample_data(ps.decontam)$LibrarySize != 0, ps.decontam)
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


contaminants <- as.data.frame(otu_table(prune_taxa(contamdf.prev$contaminant, ps.decontam)))
contaminants_tax <- merge(tax_tab[rownames(contaminants),], contaminants, by = 0)

write.xlsx(contaminants_tax, file.path(qc.dir, paste0("Contaminants", thr, ".xlsx")))

physeq_decont <- prune_taxa(taxa_sums(physeq) != 0, physeq)
physeq_decont <- prune_taxa(!taxa_names(physeq_decont) %in% rownames(contaminants), physeq_decont)

physeq_decont
physeq
physeq = physeq_decont

sample_data(physeq)$LibrarySizeDecontam <- sample_sums(physeq)

#### End ####

#### Calculate and add alpha diversity to phyloseq object for later analysis ####

sample_data(physeq) <- estimate_richness(physeq, measures = c("Observed","Chao1", "Shannon","Simpson","InvSimpson","ACE")) %>%
  merge(sample_data(physeq), ., by = 0) %>% tibble::column_to_rownames(var = "Row.names")
sample_data(physeq)$InvSimpson[sample_data(physeq)$InvSimpson == Inf] = 0

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

#save.image(file.path(R.dir,"4.physeq.decontam_2024.10.RData"))
load(file.path(R.dir,"4.physeq.decontam.RData"))

#### End ####

#### Remove NA Phyla####

rank_names(physeq)

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

physeq_oNA <- subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

physeq
physeq_oNA

# > physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 23378 taxa and 496 samples ]
# sample_data() Sample Data:       [ 496 samples by 67 sample variables ]
# tax_table()   Taxonomy Table:    [ 23378 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 23378 tips and 23376 internal nodes ]
# refseq()      DNAStringSet:      [ 23378 reference sequences ]
# > physeq_oNA
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 23131 taxa and 496 samples ]
# sample_data() Sample Data:       [ 496 samples by 67 sample variables ]
# tax_table()   Taxonomy Table:    [ 23131 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 23131 tips and 23129 internal nodes ]
# refseq()      DNAStringSet:      [ 23131 reference sequences ]

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

#save ASV Prevalence and Abundance table before filtering

write.table(prevdf, file.path(res.dir,"asv_prevdf.tsv"), sep="\t", quote=F, col.names=NA)

#Plot Taxa prevalence v. total counts. Each point is a different taxa. 

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(physeq),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#### End ####

#### Remove taxa not seen more than 3 times in at least 5% of the samples #### 
# This protects against an OTU with small mean & trivially large C.V.
# Setting filter parameters :

countperphyla = 3
Samplepercentage = 0.01

physeq_filtered = filter_taxa(physeq, function(x) sum(x > countperphyla) > (Samplepercentage*length(x)), TRUE)
physeq_filtered
physeq

# > physeq_filtered
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2695 taxa and 496 samples ]
# sample_data() Sample Data:       [ 496 samples by 67 sample variables ]
# tax_table()   Taxonomy Table:    [ 2695 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2695 tips and 2694 internal nodes ]
# refseq()      DNAStringSet:      [ 2695 reference sequences ]
# > physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 23131 taxa and 496 samples ]
# sample_data() Sample Data:       [ 496 samples by 67 sample variables ]
# tax_table()   Taxonomy Table:    [ 23131 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 23131 tips and 23129 internal nodes ]
# refseq()      DNAStringSet:      [ 23131 reference sequences ]

physeq = physeq_filtered

#### End ####

#### Normalize number of reads in each sample using median sequencing depth.####

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


#plot the abundance values before and after transformation

pl_ab_original  = plot_abundance(physeq,"Original Abundances")
pl_ab_original_norm  =plot_abundance(physeq_mednorm,"Normalized to squencing depth Abundances")
pl_ab_original_norm_re  =plot_abundance(physeq_re,"Normalized Relative Abundances")

grid.arrange(pl_ab_original, pl_ab_original_norm, pl_ab_original_norm_re)

#### End ####

#### Save RData ####
# Save as .RData to load in following steps or continue later

save.image(file.path(R.dir,"5.phyloseq.filtered.2024.10.RData"))
load(file.path(R.dir,"5.phyloseq.filtered.RData"))

#### End ####
