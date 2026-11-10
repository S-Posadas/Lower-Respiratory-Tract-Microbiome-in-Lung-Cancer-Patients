# HAllA input: diagnosis and smoking vs ASVs #
# Figure 2E, Supplementary Figure 6.
#
# Samples: directly isolated, diseased lung, NSCLC / SCLC /
# benign, with known smoking status (104 samples).
# X: batch-corrected relative abundances of the ASVs present in
#    more than 2 % of these samples, CLR-transformed (pseudocount: half the
#    smallest non-zero value; microbiome::transform).
# Y: Diagnosis (1 = NSCLC, 2 = SCLC, 3 = Benign) and Smoking (1 = No, 2 = Yes).

#### Package setup ####

suppressPackageStartupMessages({
  library(phyloseq)
  library(microbiome)
})

#### End ####

#### Load batch-corrected phyloseq object ####

res.dir <- "Results_direct"
halla.dir <- file.path(res.dir, "halla", "diagnosis")
dir.create(halla.dir, recursive = T, showWarnings = F)

rdata <- new.env()
load(file.path(res.dir, "RData", "6.phyloseq.corrected.RData"), envir = rdata)
physeq_halla <- rdata$physeq_re_bc
rm(rdata)

prevalence_cutoff <- 0.02  

#### End ####

#### Samples and ASVs ####

sd <- data.frame(sample_data(physeq_halla))
keep <- sd$Sample_or_Control == "True sample" &
  sd$Isolation == "Direct_Isolation" &
  sd$Lung == "Diseased" &
  sd$Diagnosis %in% c("Benign", "NSCLC", "SCLC") &
  !is.na(sd$History.of.smoking.y.n)
physeq_halla <- prune_samples(keep, physeq_halla)
physeq_halla <- prune_taxa(taxa_sums(physeq_halla) > 0, physeq_halla)
cat("HAllA samples:", nsamples(physeq_halla), "\n")

#### End ####

#### X: CLR-transformed relative abundances ####

X <- as(otu_table(physeq_halla), "matrix")
if (!taxa_are_rows(physeq_halla)) X <- t(X)
X <- X[rowMeans(X > 0) > prevalence_cutoff, ]
# ASVs in numeric order (the order of the taxa in the phyloseq object follows
# the tree; HAllA's clustering can depend on the input order)
X <- X[order(as.numeric(sub("^ASV", "", rownames(X)))), ]
X_clr <- as.matrix(microbiome::transform(X, "clr"))
cat("HAllA ASVs:", nrow(X_clr), "\n")

write.table(X_clr, file.path(halla.dir, "matrixX_clr.txt"), sep = "\t", quote = FALSE)

#### End ####

#### Y: diagnosis and smoking (numeric codes) ####

sd <- data.frame(sample_data(physeq_halla))
Y <- rbind(Diagnosis = c(NSCLC = 1, SCLC = 2, Benign = 3)[as.character(sd$Diagnosis)],
           Smoking = c(No = 1, Yes = 2)[as.character(sd$History.of.smoking.y.n)])
colnames(Y) <- rownames(sd)
stopifnot(identical(colnames(Y), colnames(X_clr)), !anyNA(Y))

write.table(Y, file.path(halla.dir, "matrixY.txt"), sep = "\t", quote = FALSE)

#### End ####
