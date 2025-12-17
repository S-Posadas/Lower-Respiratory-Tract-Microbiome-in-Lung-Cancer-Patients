#### Package setup ####

Sys.setenv(language = "EN")

library(dplyr)
library(openxlsx)

#### Set Directory ####

res.dir <- "Results"
blast.dir <- file.path(res.dir, "2.blast")

#### End ####

#### Format BLAST table to get a single taxa per ASV ####

# Read BLAST alignment table
blast_data <- read.csv(file.path(blast.dir, "tax_tab_95_10align_16S"), row.names = NULL, sep = "")

# Fix column names 
colnames(blast_data) <- c(colnames(blast_data)[-1], "Species")
colnames(blast_data)[colnames(blast_data) == "stitle"] <- "Genus"

# Calculate total number of ASVs from the last ASV ID
last_asv <- blast_data[nrow(blast_data), "seqid"]
total_asvs <- as.numeric(gsub("ASV_", "", last_asv))

# Ensure numeric columns are properly formatted
blast_data$pident <- as.numeric(blast_data$pident)  # Percentage identity
blast_data$evalue <- as.numeric(blast_data$evalue)  # E-value

dim(blast_data)

# Filter by percentage identity

identity_threshold <- 99
filtered_by_identity <- blast_data[blast_data$pident >= identity_threshold, ]

dim(filtered_by_identity)
unique(filtered_by_identity$Genus)
length(unique(filtered_by_identity$Genus))

# Select best hits per ASV
# For each ASV, keep alignments with:
# 1. Highest percentage identity
# 2. Lowest E-value (for ties)
# 3. Single row if only one genus/species remains

filtered_hits <- list()

for (current_asv in unique(filtered_by_identity$seqid)) {
  # Get all hits for current ASV
  asv_hits <- filtered_by_identity[filtered_by_identity$seqid == current_asv, ]
  
  # Keep hits with maximum identity
  asv_hits <- asv_hits[asv_hits$pident == max(asv_hits$pident), ]
  
  # Among those, keep hits with minimum E-value
  asv_hits <- asv_hits[asv_hits$evalue == min(asv_hits$evalue), ]
  
  # If only one genus and species remains, keep single row
  if (length(unique(asv_hits$Genus)) == 1 & 
      length(unique(asv_hits$Species)) == 1) {
    asv_hits <- asv_hits[1, ]
  }
  
  filtered_hits[[current_asv]] <- asv_hits
}

# Combine filtered hits into dataframe
best_hits <- as.data.frame(do.call(rbind, filtered_hits))

dim(best_hits)
min(best_hits$pident)

# Mark ASVs with multiple different genera as ambiguous
ambiguous_genera <- best_hits %>%
  group_by(seqid) %>%
  filter(n_distinct(Genus) > 1)

length(unique(ambiguous_genera$seqid))

# Set Genus and Species to NA for ambiguous ASVs
best_hits[best_hits$seqid %in% unique(ambiguous_genera$seqid), c("Genus", "Species")] <- NA

# Mark ASVs with multiple species (same genus) as ambiguous at species level
ambiguous_species <- best_hits %>%
  group_by(seqid) %>%
  filter(n_distinct(Species) > 1)

length(unique(ambiguous_species$seqid))

# Set Species to NA for ambiguous ASVs
best_hits[best_hits$seqid %in% unique(ambiguous_species$seqid), "Species"] <- NA

# Keep only one row per ASV
final_hits <- list()
for (current_asv in unique(best_hits$seqid)) {
  asv_data <- best_hits[best_hits$seqid == current_asv, ]
  final_hits[[current_asv]] <- asv_data[1, ]  # Take first row
}

blast_final <- as.data.frame(do.call(rbind, final_hits))

dim(blast_final)
min(as.numeric(blast_final$pident))

# Add unassigned ASVs:

# Extract ASV numbers
blast_final$asv_number <- as.numeric(gsub("ASV_", "", blast_final$seqid))

# Identify missing ASVs
all_asv_numbers <- 1:total_asvs
missing_asv_numbers <- all_asv_numbers[!all_asv_numbers %in% blast_final$asv_number]

if (length(missing_asv_numbers) > 0) {
  # Create placeholder rows for missing ASVs
  missing_asvs <- paste("ASV", missing_asv_numbers, sep = "_")
  placeholder_rows <- cbind(
    seqid = missing_asvs,
    matrix(NA, 
           nrow = length(missing_asvs), 
           ncol = ncol(blast_final) - 1)
  )
  colnames(placeholder_rows) <- colnames(blast_final)
  
  # Add placeholder rows to final table
  blast_final <- rbind(blast_final, placeholder_rows)
}

# Sort by ASV number
blast_final$asv_number <- as.numeric(gsub("ASV_", "", blast_final$seqid))
blast_final <- blast_final[order(blast_final$asv_number), ]
blast_final$asv_number <- NULL  # Remove helper column

#Save filtered BLAST results
write.xlsx(blast_final, file.path(blast.dir,"blast_final_99_16S.xlsx"))

#### End ####

#### Add blast species to unassigned SILVA species ####

# Read taxa table from SILVA
tax_tab_SILVA <- read.xlsx(file.path(blast.dir,"tax_tab_SILVA.xlsx"), rowNames = T)
dim(tax_tab_SILVA)
colnames(tax_tab_SILVA) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species.bayes", "Species.100")

# Read taxa table from blast
blast <- read.xlsx(file.path(blast.dir, "blast_final_99_16S.xlsx"), rowNames = T)
dim(blast)

# Change blast taxonomy to SILVA taxonomy based on known synonyms
correction_list <- list(
  # Format: c("Incorrect Genus (BLAST)", "Species pattern (BLAST)", "Corrected Genus (SILVA)")
  c("Mycolicibacterium", ".*", "Mycobacterium"),
  c("Phocaeicola", "(vulgatus|dorei|massiliensis)", "Bacteroides"),
  c("Hallella", "(multisaccharivorax|mizrahii|colorans)", "Prevotella"),
  c("Ureibacillus", "(massiliensis|chungkukjangi)", "Lysinibacillus"),
  c("[Ruminococcus]", "(faecis|torques)", "Mediterraneibacter"),
  c("Roseateles", "oligotrophus", "Paucibacter"),
  c("Roseateles", "asaccharophilus", "Kinneretia"),
  c("Macrococcoides", ".*", "Macrococcus"),
  c("Metamycoplasma", ".*", "Mycoplasma")
)

# Apply corrections
for (correction in correction_list) {
  old_genus <- correction[1]
  species_pattern <- correction[2]
  new_genus <- correction[3]
  
  # Find rows matching the correction criteria
  match_rows <- !is.na(blast$Genus) & 
    blast$Genus == old_genus &
    grepl(species_pattern, blast$Species)
  
  if (any(match_rows)) {
    blast[match_rows, "Genus"] <- new_genus}
}

# Special case: fix species name spelling
blast[!is.na(blast$Genus) & blast$Genus == "Kinneretia" & blast$Species == "asaccharophilus", "Species"] = "asaccharophila"

# Find ASVs where BLAST can provide species-level resolution
# when SILVA has missing species assignments
samegenus_asvs <- rownames(blast[
  (is.na(tax_tab_SILVA$Species.100) | grepl("/", tax_tab_SILVA$Species.100)) &
    !is.na(blast$Species) &
    !is.na(tax_tab_SILVA$Genus) &
    blast$Genus == tax_tab_SILVA$Genus, ])
length(samegenus_asvs)

# Find ASVs with genus disagreements
conflicting_genera <- rownames(blast[
  !is.na(blast$Genus) &
    !is.na(tax_tab_SILVA$Genus) &
    blast$Genus != tax_tab_SILVA$Genus, ])
length(conflicting_genera)

# Save discrepant genera
comparison_table <- merge(
  tax_tab_SILVA[conflicting_genera, c("Genus", "Species.bayes", "Species.100")],
  blast[conflicting_genera, c("Genus", "Species")],
  by = "row.names")

write.xlsx(comparison_table, file.path(blast.dir, "discrepant_SILVA_blast_genus.xlsx"), rowNames = TRUE)

# Integrate BLAST species into SILVA table
tax_tab_SILVA$Species <- tax_tab_SILVA$Species.100
tax_tab_SILVA[samegenus_asvs, "Species"] <- blast[samegenus_asvs, "Species"]

# Save final taxa table
write.xlsx(tax_tab_SILVA, file.path(blast.dir, "tax_tab_SILVA_blast_species.xlsx"), rowNames = T)

#### End ####
