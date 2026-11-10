# PICRUSt2 input #
# ASV sequences (seqs.fna) and raw counts (raw_counts.tsv) of the directly
# isolated true samples with reads left after decontamination, from the
# filtered phyloseq object (before MMUPHin: PICRUSt2 needs raw counts).

#### Package setup ####

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
})

#### End ####

#### Load filtered phyloseq object (1.DADA2_direct.R) ####

res.dir <- "Results_direct"

rdata <- new.env()
load(file.path(res.dir, "RData", "5.phyloseq.filtered.RData"), envir = rdata)
physeq_picrust <- rdata$physeq_filtered
rm(rdata)

picrust.dir <- file.path(res.dir, "picrust2")
dir.create(picrust.dir, recursive = T, showWarnings = F)

#### End ####

#### Samples and ASVs ####

lib_size <- colSums(otu_table(physeq_picrust))
physeq_picrust <- prune_samples(!is.na(lib_size) & lib_size != 0, physeq_picrust)

sd <- data.frame(sample_data(physeq_picrust))
physeq_picrust <- prune_samples(sd$Sample_or_Control == "True sample" & sd$Isolation == "Direct_Isolation" &
                                  sd$LibrarySizeDecontam > 0, physeq_picrust)
physeq_picrust <- prune_taxa(taxa_sums(physeq_picrust) != 0, physeq_picrust)
if (!taxa_are_rows(physeq_picrust)) physeq_picrust <- t(physeq_picrust)

cat("Samples:", nsamples(physeq_picrust), "| ASVs:", ntaxa(physeq_picrust), "\n")

#### End ####

#### Export ####

writeXStringSet(refseq(physeq_picrust), file.path(picrust.dir, "seqs.fna"), format = "fasta")

counts <- as(otu_table(physeq_picrust), "matrix")
write.table(data.frame(ASV = rownames(counts), counts, check.names = FALSE),
            file.path(picrust.dir, "raw_counts.tsv"), sep = "\t", quote = F, row.names = F)

#### End ####
