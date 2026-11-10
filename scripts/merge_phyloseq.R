#### Introduction ####

# Combines the phyloseq objects of the directly isolated samples (1.DADA2_direct.R) and of the
# culture-enriched samples (1.DADA2_culture.R) into one object, for the analyses that compare
# both sample types.
# - ASVs with the same sequence in both objects get the name of the directly isolated object.
#   Culture-only ASVs get new names, continuing the numbering (ASV<max + 1>, ...).
# - For shared ASVs the taxonomy of the directly isolated object is kept (differences are printed).
# - Only true samples are taken from the culture-enriched object: the controls were processed in
#   both runs and are taken from the directly isolated object.
# - The phylogenetic tree of the directly isolated object is discarded
# - A "source" column ("direct" / "cult") is added to the sample data.

#### End ####

#### Package setup ####

Sys.setenv(language = "EN")

library(phyloseq)
library(Biostrings)  # for DNAStringSet manipulation
library(dplyr)       # for data frame operations

#### End ####

#### Set Directory ####

direct.dir <- file.path("Results_direct", "RData")
cult.dir   <- file.path("Results_cult", "RData")
out.dir    <- file.path("Results_merged")
R.dir      <- file.path(out.dir, "RData")
dir.create(R.dir, recursive = T, showWarnings = F)

# Phyloseq objects to merge
objects_to_merge <- list(physeq_re_bc    = c(file = "physeq_re_bc.rds"),
                         physeq_filtered = c(file = "5.phyloseq.filtered.RData", object = "physeq_filtered"))

read_object <- function(dir, src) {
  if (is.na(src["object"])) return(readRDS(file.path(dir, src[["file"]])))
  e <- new.env()
  load(file.path(dir, src[["file"]]), envir = e)
  e[[src[["object"]]]]
}

#### End ####

#### Merge function ####

merge_phyloseq_by_sequence <- function(ps_direct, ps_cult) {

  # 1. Extract components ----------------------------------------------------
  otu_direct <- as(otu_table(ps_direct), "matrix")
  otu_cult   <- as(otu_table(ps_cult), "matrix")
  if (!taxa_are_rows(ps_direct)) otu_direct <- t(otu_direct)
  if (!taxa_are_rows(ps_cult))   otu_cult   <- t(otu_cult)

  tax_direct <- as(tax_table(ps_direct), "matrix")
  tax_cult   <- as(tax_table(ps_cult), "matrix")

  sam_direct <- as(sample_data(ps_direct), "data.frame")
  sam_cult   <- as(sample_data(ps_cult), "data.frame")

  ref_direct <- refseq(ps_direct)
  ref_cult   <- refseq(ps_cult)

  stopifnot(!any(colnames(otu_cult) %in% colnames(otu_direct)))  # Sample names must be unique

  if (!is.null(phy_tree(ps_direct, errorIfNULL = FALSE)) || !is.null(phy_tree(ps_cult, errorIfNULL = FALSE))) {
    message("Phylogenetic trees are discarded")
  }

  # 2. New names for the culture ASVs -----------------------------------------
  # Same sequence as a direct ASV: the direct name; otherwise the next free number
  seqs_direct <- setNames(as.character(ref_direct), names(ref_direct))
  cult_asv_old <- names(ref_cult)
  matches <- match(as.character(ref_cult), seqs_direct)
  matched_idx <- !is.na(matches)

  direct_numbers <- as.integer(sub("^ASV", "", names(seqs_direct)))
  stopifnot(!anyNA(direct_numbers))  # Direct ASV names are ASV1, ASV2, ...

  new_cult_names <- character(length(cult_asv_old))
  new_cult_names[matched_idx] <- names(seqs_direct)[matches[matched_idx]]
  new_cult_names[!matched_idx] <- paste0("ASV", max(direct_numbers) + seq_len(sum(!matched_idx)))

  rename_map <- setNames(new_cult_names, cult_asv_old)
  message(sum(matched_idx), " culture ASVs match a direct ASV; ", sum(!matched_idx), " are new")

  # 3. Check taxonomy for matching ASVs --------------------------------------
  matched_cult_asvs   <- cult_asv_old[matched_idx]
  matched_direct_asvs <- new_cult_names[matched_idx]

  if (length(matched_cult_asvs) > 0) {
    tax_matched_cult   <- tax_cult[matched_cult_asvs, , drop = FALSE]
    tax_matched_direct <- tax_direct[matched_direct_asvs, colnames(tax_cult), drop = FALSE]

    # Rows with a difference (NA in one table and a name in the other counts as equal)
    diff_rows <- !apply(tax_matched_cult == tax_matched_direct, 1, all, na.rm = TRUE)

    if (any(diff_rows)) {
      warning("Taxonomy differs for ", sum(diff_rows), " matching ASVs; the direct taxonomy is kept. See details below.")
      for (i in which(diff_rows)) {
        cat("\nASV:", matched_direct_asvs[i], "\n")
        cat("  Taxonomy in ps_direct:\n")
        print(tax_matched_direct[i, , drop = FALSE])
        cat("  Taxonomy in ps_cult:\n")
        print(tax_matched_cult[i, , drop = FALSE])
      }
    }
  }

  # 4. Harmonize sample data columns -----------------------------------------
  all_cols <- union(colnames(sam_direct), colnames(sam_cult))
  for (col in all_cols) {
    if (!col %in% colnames(sam_direct)) sam_direct[[col]] <- NA
    if (!col %in% colnames(sam_cult))   sam_cult[[col]]   <- NA
  }
  sam_direct <- sam_direct[, all_cols, drop = FALSE]
  sam_cult   <- sam_cult[, all_cols, drop = FALSE]

  # 5. Rename components of ps_cult ------------------------------------------
  rownames(otu_cult) <- rename_map[rownames(otu_cult)]
  rownames(tax_cult) <- rename_map[rownames(tax_cult)]
  names(ref_cult)    <- rename_map[names(ref_cult)]

  # 6. Combine OTU tables (ASVs absent from one sample type get 0) -----------
  all_samples <- c(colnames(otu_direct), colnames(otu_cult))
  all_asvs <- union(rownames(otu_direct), rownames(otu_cult))

  combined_otu <- matrix(0, nrow = length(all_asvs), ncol = length(all_samples),
                         dimnames = list(all_asvs, all_samples))
  combined_otu[rownames(otu_direct), colnames(otu_direct)] <- otu_direct
  combined_otu[rownames(otu_cult), colnames(otu_cult)] <- otu_cult

  # 7. Combine taxonomy tables (direct taxonomy for shared ASVs) -------------
  new_asvs <- setdiff(rownames(tax_cult), rownames(tax_direct))
  tax_combined <- as.matrix(bind_rows(as.data.frame(tax_direct),
                                      as.data.frame(tax_cult)[new_asvs, , drop = FALSE]))
  rownames(tax_combined) <- c(rownames(tax_direct), new_asvs)

  # 8. Combine sample data ----------------------------------------------------
  sam_direct$source <- "direct"
  sam_cult$source   <- "cult"
  sam_combined <- bind_rows(sam_direct, sam_cult)
  rownames(sam_combined) <- c(rownames(sam_direct), rownames(sam_cult))

  # 9. Combine reference sequences --------------------------------------------
  ref_combined <- c(ref_direct, ref_cult[names(ref_cult) %in% new_asvs])

  # 10. Create new phyloseq object --------------------------------------------
  phyloseq(otu_table(combined_otu, taxa_are_rows = TRUE),
           tax_table(tax_combined),
           sample_data(sam_combined),
           refseq(ref_combined))
}

#### End ####

#### Merge and save ####

for (ps in names(objects_to_merge)) {

  ps_direct <- read_object(direct.dir, objects_to_merge[[ps]])
  ps_cult   <- read_object(cult.dir, objects_to_merge[[ps]])

  # Only true samples from the culture-enriched object (the controls are in both objects)
  ps_cult <- subset_samples(ps_cult, Sample_or_Control == "True sample")

  merged_ps <- merge_phyloseq_by_sequence(ps_direct, ps_cult)

  print(ps_direct)
  print(ps_cult)
  print(merged_ps)

  saveRDS(merged_ps, file.path(R.dir, paste0(ps, ".rds")))
}

#### End ####
