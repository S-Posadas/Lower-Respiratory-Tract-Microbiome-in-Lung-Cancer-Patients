#### Package setup ####

Sys.setenv(language = "EN")
library(Maaslin2)
library(openxlsx)
library(phyloseq)
library(ggrepel)
theme_set(theme_bw())

#### End ####

#### Load RData with phyloseq object ####

res.dir <- "Results_trimmed_cutadapt3_20241216"
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir,"5.phyloseq.filtered.RData"))

# Create directory for results

maaslin.dir=file.path(res.dir, "Maaslin2")
dir.create(maaslin.dir)

# Change species names for visualization in Maaslin

tax <- as.data.frame(tax_table(physeq_re))   
tax[] <- lapply(tax, gsub, pattern = " ", replacement = "_")
tax$Species[grep("_spp.", tax$Species)] <- NA

tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(physeq_re) <- tax; rm(tax)

#### End ####

#### Agglomerate taxa ####

physeq_ph <- tax_glom(physeq_maas0, "Phylum")
taxa_names(physeq_ph) <- tax_table(physeq_ph)[,"Phylum"]
physeq_g <- tax_glom(physeq_maas0, "Genus")
taxa_names(physeq_g) <- tax_table(physeq_g)[,"Genus"]

#### End ####

#### Select phyloseq object #### 

physeq_maas0 <- physeq_g

# Add LibrarySizeGroup
sample_data(physeq_maas0)$LibrarySizeGroup <- ifelse(sample_data(physeq_maas0)$LibrarySizeDecontam >= 1000, ">=1000", "<1000")

# Remove controls
physeq_maas0 <- subset_samples(physeq_maas0, sample_names(physeq_maas0) %in% c(rownames(sample_info[!sample_info$Sample_or_Control %in% c("Negative control sample", "Positive control sample"),]))) 
physeq_maas0

# Remove samples with no reads
physeq_maas0 <- subset_samples(physeq_maas0, !is.na(colSums(otu_table(physeq_maas0)))) 
physeq_maas0

# Keep only directly isolated samples
physeq_maas0 <- subset_samples(physeq_maas0, sample_data(physeq_maas0)$Isolation == "Direct_Isolation") 
physeq_maas0

# Select main diagnosis
physeq_maas0 = subset_samples(physeq_maas0, sample_data(physeq_maas0)$Diagnosis %in% c("NSCLC", "SCLC", "Benign"))
physeq_maas0

# Select filter for minimum ASV count

n = 0
length(sample_sums(physeq)[sample_sums(physeq) >=n])
less_n <- names(sample_sums(physeq)[sample_sums(physeq) <n])
physeq_maas0 = subset_samples(physeq_maas0, !sample_names(physeq_maas0) %in% less_n) 
x = length(sample_names(physeq_maas0))
x

maas.dir=file.path(maaslin.dir, paste0("more_than_", n, "_counts_", x, "_samples_genus"))
dir.create(maas.dir)

#### Format metadata ####

metadata <- data.frame(sample_data(physeq_maas0))
metadata[] <- lapply(metadata, as.factor)
colnames(metadata)
# vars <- c("DNA_concentration", "Batch", "Lung",                     
#   "Age", "Sex",                    
#   "Active.smoker"  ,  "Pack.years"            ,     "History.of.smoking.y.n" ,
#   "FEV1.l."                , "FEV1."                   ,               
#   "Side"                     , "Lobe"                      ,             
#   "Synchronous.tumor"          ,     "Diagnosis"                   ,           
#   "Histology.NSCLC"                        , "T", "N", "M",
#   "OP"                                     ,  "VATS.RATS",  "Anatomical.resection" , "Postop..pneumonia"                      ,
#   "Microbial.detection"                    ,  "Prolonged.air.leak...7.days."           ,
#   "Bronchial.dehiscence"                   , "Empyema" , "Wound.dehiscence" ,
#   "History.of.smoking.y.n"                 ,
#   "LibrarySizeDecontam"                    ,
#   "Observed"                               ,  "Chao1"  , "ACE"  ,"Shannon"   , "Simpson"   ,  "InvSimpson")
# metadata <- metadata[,vars]

#### Function for volcano plot ####
volcano <- function(maas.dir = maas.dir, name = name){
  maaslin_all <- read.csv(file.path(maas.dir, name, "all_results.tsv"), sep = "\t")
  
  ggplot(maaslin_all, aes(x = coef, y = -log10(qval))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.1), col = "red") +
    annotate("text", x = min_prev, y = 1.5, label = "q-value = 0.1", col = "red", size = 6) +
    geom_text_repel(aes(label = feature)) + theme(axis.text = element_text(size = 18),
                                                  axis.title = element_text(size = 20))
  ggsave(file.path(maas.dir, name, "all_results.svg"), dpi = 300, width = 8, height = 8)
  ggsave(file.path(maas.dir, name, "all_results.tiff"), dpi = 300, width = 8, height = 8)
  
}

#### Cut offs ####
min_abund = 0.001
min_prev = 0.1 

#### Q0: Is there a difference between batches? ####

name = "main_diagnosis_batch"
data <- data.frame(otu_table(physeq_maas0))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "Batch",
                  # reference = "Batch,Normal",
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

#### Q0.1: Is there a difference between small and big library sizes? ####

name = "all_directly_isolated_libsize"
data <- data.frame(otu_table(physeq_maas0))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "LibrarySizeGroup",
                  # reference = "Batch,Normal",
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

name = "all_directly_isolated_libsize_batch"
data <- data.frame(otu_table(physeq_maas0))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "LibrarySizeGroup",
                  random_effects = c("Batch"),
                  reference = "LibrarySizeGroup,<1000",
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

#### Q1: Is there a difference between Diseased vs the parallel Normal Lung? ####

name = "all_directly_isolated_lung"
data <- data.frame(otu_table(physeq_maas0))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
  min_abundance = min_abund,  min_prevalence = min_prev,
  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
  #  max_significance = 0.05,
  fixed_effects = "Lung",
  #  random_effects = c("Batch"),
  reference = "Lung,Normal",
  correction = "BH",
  standardize = FALSE)

name = "all_directly_isolated_lung_batch"
data <- data.frame(otu_table(physeq_maas0))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
  min_abundance = min_abund,  min_prevalence = min_prev,
  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
  #  max_significance = 0.05,
  fixed_effects = "Lung",
  random_effects = c("Batch"),
  reference = "Lung,Normal",
  correction = "BH",
  standardize = FALSE)

#### Q2: Is there a difference between Diseased vs the parallel Normal Lung without samples with synchronous tumor (after decontamination, only directly isolated samples) ####

physeq_maas <- physeq_maas0

# Remove non synchronous tumor
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Synchronous.tumor == "No")

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "all_directly_isolated_non_synchronous_batch"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
  min_abundance = min_abund,  min_prevalence = min_prev,
  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
  #  max_significance = 0.05,
  fixed_effects = "Batch",
  # reference = "Batch,Normal",
  correction = "BH",
  standardize = FALSE)

name = "all_directly_isolated_non_synchronous_lung"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
  min_abundance = min_abund,  min_prevalence = min_prev,
  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
  #  max_significance = 0.05,
  fixed_effects = "Lung",
  #  random_effects = c("Batch"),
  reference = "Lung,Normal",
  correction = "BH",
  standardize = FALSE)

name = "all_directly_isolated_non_synchronous_lung_batch"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
  min_abundance = min_abund,  min_prevalence = min_prev,
  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
  #  max_significance = 0.05,
  fixed_effects = "Lung",
  random_effects = c("Batch"),
  reference = "Lung,Normal",
  correction = "BH",
  standardize = FALSE)


#### Q3: What are the difference between NSCLC vs SCLC vs Bening Tumor? ####

physeq_maas <- physeq_maas0

# Keep only diseased lung
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Lung == "Diseased") 

# Select NSCLC
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Diagnosis %in% c("NSCLC", "SCLC", "Benign"))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "diseased_diagnosis"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
  min_abundance = min_abund,  min_prevalence = min_prev,
  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
  #  max_significance = 0.05,
  fixed_effects = "Diagnosis",
 # reference = "Diagnosis,Benign",
  correction = "BH",
  standardize = FALSE)
volcano(maas.dir, name)

name = "diseased_diagnosis_batch"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
  min_abundance = min_abund,  min_prevalence = min_prev,
  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
  #  max_significance = 0.05,
  fixed_effects = "Diagnosis",
  # reference = "Batch,Normal",
  random_effects = c("Batch"),
  correction = "BH",
  standardize = FALSE)

#### Q4: What are the difference between NSCLC main histologies (Adeno vs Squamous)? ####

physeq_maas <- physeq_maas0

# Keep only diseased lung
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Lung == "Diseased") 

# Select NSCLC
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Diagnosis %in% c("NSCLC"))

# Select histology
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "histology_diagnosis"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "Histology.NSCLC",
                  correction = "BH",
                  standardize = FALSE)

name = "diseased_histology_batch"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "Histology.NSCLC",
                  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)

#### Q5: Does history of smoking have an impact on the lung microbiome? ####
physeq_maas <- physeq_maas0

# Keep only diseased lung
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Lung == "Diseased") 

# Remove NA in history of smoking
physeq_maas = subset_samples(physeq_maas, !is.na(sample_data(physeq_maas)$History.of.smoking.y.n))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "all_directly_isolated_smoking"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "History.of.smoking.y.n",
                  #  random_effects = c("Batch"),
                  plot_heatmap = T,
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

name = "all_directly_isolated_batch_smoking"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "History.of.smoking.y.n",
                  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)

# Select NSCLC
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Diagnosis %in% c("NSCLC"))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "diseased_NSCLC_smoking"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "History.of.smoking.y.n",
                  #  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

name = "diseased_NSCLC__batch_smoking"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "History.of.smoking.y.n",
                  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)

volcano(maas.dir, name)

#### Q6c: Does the M stage have an impact on the lung microbiome? ####
physeq_maas <- physeq_maas0

# Keep only diseased lung
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Lung == "Diseased") 

# Remove NA in M
physeq_maas = subset_samples(physeq_maas, !is.na(sample_data(physeq_maas)$M))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "all_directly_isolated_M"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "M",
                  #  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

name = "all_directly_isolated_batch_M"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "M",
                  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)

# Select NSCLC
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Diagnosis %in% c("NSCLC"))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "diseased_NSCLC_M"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "M",
                  #  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

name = "diseased_NSCLC__batch_M"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "M",
                  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)

volcano(maas.dir, name)

#### Q6b: Does the N stage have an impact on the lung microbiome? ####
physeq_maas <- physeq_maas0

# Keep only diseased lung
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Lung == "Diseased") 

# Remove NA in M
physeq_maas = subset_samples(physeq_maas, !is.na(sample_data(physeq_maas)$N))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "all_directly_isolated_N"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "N",
                  #  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

name = "all_directly_isolated_batch_N"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "N",
                  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)

# Select NSCLC
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Diagnosis %in% c("NSCLC"))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "diseased_NSCLC_N"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "N",
                  #  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

name = "diseased_NSCLC_batch_N"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "N",
                  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)

volcano(maas.dir, name)

#### Q6a: Does the T stage have an impact on the lung microbiome? ####
physeq_maas <- physeq_maas0

# Keep only diseased lung
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Lung == "Diseased") 

# Remove NA in M
physeq_maas = subset_samples(physeq_maas, !is.na(sample_data(physeq_maas)$T))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "all_directly_isolated_T"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "T",
                  #  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

name = "all_directly_isolated_batch_T"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "T",
                  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)

# Select NSCLC
physeq_maas = subset_samples(physeq_maas, sample_data(physeq_maas)$Diagnosis %in% c("NSCLC"))

# Remove not present taxa
physeq_maas <- prune_taxa(taxa_sums(physeq_maas) != 0, physeq_maas)
physeq_maas

# Select metadata
metadata <- data.frame(sample_data(physeq_maas))
metadata[] <- lapply(metadata, as.factor)

name = "diseased_NSCLC_T"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "T",
                  #  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)
volcano(maas.dir, name)

name = "diseased_NSCLC_batch_T"
data <- data.frame(otu_table(physeq_maas))
all(colnames(data) == rownames(metadata))
mas_1 <- Maaslin2(input_data = data,  input_metadata = metadata,  output = file.path(maas.dir, name),
                  min_abundance = min_abund,  min_prevalence = min_prev,
                  normalization = "TSS",  transform = "LOG",  analysis_method = "LM",
                  #  max_significance = 0.05,
                  fixed_effects = "T",
                  random_effects = c("Batch"),
                  correction = "BH",
                  standardize = FALSE)

volcano(maas.dir, name)


############# 


#### Corncob ####

library(corncob)
corn_da <- differentialTest(formula = ~ Stentinsgesamt,
                            phi.formula = ~ 1,
                            formula_null = ~ 1,
                            phi.formula_null = ~ 1,
                            data = physeq_re,
                            test = "Wald", boot = FALSE,
                            fdr_cutoff = 0.05)

fdr_corncob <- corn_da$significant_taxa
dim(data.frame(fdr_corncob))

head(sort(corn_da$p_fdr))  
