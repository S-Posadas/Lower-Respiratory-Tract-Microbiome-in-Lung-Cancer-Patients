
# 6.Taxonomy analysis #

#### Package setup ####

library(phyloseq)
library(ggplot2)
library(gtable)
library(ggplotify)
library(openxlsx)
library(dplyr)
library(pgirmess)

theme_set(theme_bw())

#### End ####

#### Load RData with phyloseq object ####

res.dir <- "Results_trimmed_cutadapt3_20241216"
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir, "2.dada2_taxonomy.RData"))

#### End ####

#### Plot theme settings ####
text_size = 30
thm <- theme(text = element_text(size = text_size),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "bottom")

thm.x <- theme(text = element_text(size = text_size),
               axis.text.x = element_text(angle = -30, hjust = 0),
               legend.position = "bottom")

thm.x.2 <- theme(text = element_text(size = text_size),
                 axis.text.x = element_text(angle = -90),
                 legend.position = "bottom")
#### End ####

#### Comparing different databases ####

# Examine all tables

head(count_tab)
head(fasta_tab)
head(sample_info)
head(tax_tab_SILVA)
colnames(tax_tab_SILVA) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species.bayes", "Species")
tax_tab_SILVA <- tax_tab_SILVA[,colnames(tax_tab_SILVA) != "Species.bayes"]
head(tax_tab_GTDB)
colnames(tax_tab_GTDB) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species.bayes", "Species")
tax_tab_GTDB <- tax_tab_GTDB[,colnames(tax_tab_GTDB) != "Species.bayes"]
head(tax_tab_eHOMD)

# Create lists with taxa tables to consider
tax_tables = list(SILVA = tax_tab_SILVA, GTDB = tax_tab_GTDB, eHOMD = tax_tab_eHOMD)

# Create lists to store the results of each database
result_tax_tables = list()
result_physeq_Mock_top = list()
result_physeq_Mock_re_top = list()

for (name in names(tax_tables)) {
  
  # Change Species names for visualization 
  tax_tab = data.frame(tax_tables[[name]])
  colnames(tax_tab)[colnames(tax_tab) == "Species"] = "Species.only"
  
  # Add the Species column based on the conditions
  tax_tab$Species = NA  # Initialize Species column with NA
  
  tax_tab$Species = ifelse(!is.na(tax_tab$Genus) & !is.na(tax_tab$Species.only),
                        paste(tax_tab$Genus, tax_tab$Species.only),
                        ifelse(!is.na(tax_tab$Genus) & is.na(tax_tab$Species.only),
                               paste(tax_tab$Genus, "spp."),
                               NA))
  tax_tab$Species.only = NULL
  tax_tab = as.matrix(tax_tab)
  
  # Create phyloseq object
  physeq <- phyloseq(otu_table(count_tab, taxa_are_rows = T),   #taxa_are_rows=F (if your taxa names on the column not the rows)
                     sample_data(sample_info), 
                     tax_table(as.matrix(tax_tab)))
  
  # Normalize number of reads in each sample using median sequencing depth
  
  total = median(sample_sums(physeq))
  standf = function(x, t=total) round(t * (x / sum(x)))
  physeq_mednorm = transform_sample_counts(physeq, standf)

  # Transform to relative abundance
  physeq_re = transform_sample_counts(physeq_mednorm, function(x){x / sum(x)})

  # Select Mock Communities
  physeq_Mock <- subset_samples(physeq, grepl("Mock", sample_data(physeq)$Study_Nr))
  sample_names(physeq_Mock) <- paste(sample_data(physeq_Mock)$Study_Nr, name, sep = "_")
  sample_data(physeq_Mock)$Study_Nr = paste(name, "database")
  physeq_Mock <- prune_taxa(taxa_sums(physeq_Mock) != 0, physeq_Mock)
  physeq_Mock_re <- subset_samples(physeq_re, grepl("Mock", sample_data(physeq)$Study_Nr))
  sample_names(physeq_Mock_re) <- paste(sample_data(physeq_Mock_re)$Study_Nr, name, sep = "_")
  sample_data(physeq_Mock_re)$Study_Nr = paste(name, "database")
  otu_table(physeq_Mock_re) <- otu_table(physeq_Mock_re) *100 # Transform to percentage
  physeq_Mock_re <- prune_taxa(taxa_sums(physeq_Mock_re) != 0, physeq_Mock_re)
  
  #sample_names(physeq_Mock) <- paste(sample_data(physeq_Mock)$Study_Nr, name, sep = "_")
  taxa_names(physeq_Mock) <- paste(taxa_names(physeq_Mock), name, sep = "_")
  #sample_names(physeq_Mock_re) <- paste(sample_data(physeq_Mock_re)$Study_Nr, name, sep = "_")
  taxa_names(physeq_Mock_re) <- paste(taxa_names(physeq_Mock_re), name, sep = "_")
    
  # Select top n ASVs
  n = 40 # adjust number to wished top ones
  top_Mock <- names(sort(taxa_sums(physeq_Mock), decreasing=TRUE))[1:n]  
  
  physeq_Mock_top <- prune_taxa(top_Mock, physeq_Mock)
  physeq_Mock_re_top <- prune_taxa(top_Mock, physeq_Mock_re)
  
  # Store the result back in the list
  result_tax_tables[[name]] = tax_tab
  result_physeq_Mock_top[[name]] = physeq_Mock_top
  result_physeq_Mock_re_top[[name]] = physeq_Mock_re_top

}

# Access the created phyloseq objects
physeq_Mock_top_SILVA = result_physeq_Mock_top$SILVA
physeq_Mock_re_top_SILVA = result_physeq_Mock_re_top$SILVA
physeq_Mock_top_GTDB = result_physeq_Mock_top$GTDB
physeq_Mock_re_top_GTDB = result_physeq_Mock_re_top$GTDB
physeq_Mock_top_eHOMD = result_physeq_Mock_top$eHOMD
physeq_Mock_re_top_eHOMD = result_physeq_Mock_re_top$eHOMD

# Create phyloseq object with expected taxa
mock_asv <- c(Bacillus_subtilis = 17,
              Enterococcus_faecalis = 10,
              Escherichia_coli = 10,
              Lactobacillus_fermentum = 18,
              Listeria_monocytogenes = 14,
              Salmonella_enterica = 10,
              Pseudomonas_aeruginosa = 4,
              Staphylococcus_aureus = 16) %>% as.data.frame(.) %>% otu_table(., taxa_are_rows = TRUE)

mock_tax <- tax_table(rbind(Bacillus_subtilis = c("","","","","", "Bacillus", "Bacillus subtilis"), 
                            Enterococcus_faecalis = c("","","","","", "Enterococcus", "Enterococcus faecalis"), 
                            Escherichia_coli = c("","","","","", "Escherichia", "Escherichia coli"), 
                            Lactobacillus_fermentum = c("","","","","", "Lactobacillus", "Lactobacillus fermentum"), 
                            Listeria_monocytogenes = c("","","","","", "Listeria", "Listeria monocytogenes"), 
                            Salmonella_enterica = c("","","","","", "Salmonella", "Salmonella enterica"), 
                            Pseudomonas_aeruginosa = c("","","","","", "Pseudomonas", "Pseudomonas aeruginosa"), 
                            Staphylococcus_aureus = c("","","","","", "Staphylococcus", "Staphylococcus aureus")))
colnames(mock_tax) <- c("Kingdom","Phylum","Phylum","Order","Family", "Genus", "Species")

# Create a new sample data frame with sample metadata
new_sample_data <- data.frame(Study_Nr = "Expected")
rownames(new_sample_data) <- "."  # Set row name to match sample ID
new_sample_data <- sample_data(new_sample_data)

# Combine the new OTU table and sample data into a new phyloseq object
physeq_expected_mock <- phyloseq(mock_asv, new_sample_data, mock_tax)

# Combine phyloseq objects from all databases
physeq_combined_abs <- merge_phyloseq(physeq_Mock_top_SILVA, physeq_Mock_top_GTDB, physeq_Mock_top_eHOMD)
sample_data(physeq_combined_abs)$Study_Nr <- factor(sample_data(physeq_combined_abs)$Study_Nr, levels = c("Expected", "SILVA database", "eHOMD database", "GTDB database"))

physeq_combined_re <- merge_phyloseq(physeq_Mock_re_top_SILVA, physeq_Mock_re_top_GTDB, physeq_Mock_re_top_eHOMD, physeq_expected_mock)
sample_data(physeq_combined_re)$Study_Nr <- factor(sample_data(physeq_combined_re)$Study_Nr, levels = c("Expected", "SILVA database", "eHOMD database", "GTDB database"))

# Plot relative Abundance

plot_bar(physeq_combined_re, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 4, title = NULL),
         color = guide_legend(ncol = 4, title = NULL))

ggsave(file.path(qc.dir, "Mock_database_comparison_rel_genus.tiff"), width = 16, height = 8, dpi = 300)
ggsave(file.path(qc.dir, "Mock_database_comparison_rel_genus.svg"), width = 16, height = 8, dpi = 300)

plot_bar(physeq_combined_re, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "Mock_database_comparison_rel_species.tiff"), width = 20, height = 22, dpi = 300)
ggsave(file.path(qc.dir, "Mock_database_comparison_rel_species.svg"), width = 20, height = 22, dpi = 300)

# Plot absolute Abundance
plot_bar(physeq_combined_abs, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm.x.2 +
  guides(fill = guide_legend(ncol = 4, title = NULL),
         color = guide_legend(ncol = 4, title = NULL))

ggsave(file.path(qc.dir, "Mock_database_comparison_abs_genus.tiff"), width = 14, height = 12, dpi = 300)
ggsave(file.path(qc.dir, "Mock_database_comparison_abs_genus.svg"), width = 14, height = 12, dpi = 300)

#### End ####

######################################################################################
################ Everything below here belongs to old script #########################
######################################################################################

#### Select phyloseq object with positive controls ####

#Take physeq study because some contaminants are ASVs from the postive controls
physeq_PC <- physeq

#### End ####

#### Change Species names for visualization ####

tax <- as.data.frame(tax_table(physeq_PC))   
colnames(tax)[colnames(tax) == "Species"] <- "Species.only"

# Add the Species column based on the conditions
tax$Species <- NA  # Initialize Species column with NA

tax$Species <- ifelse(!is.na(tax$Genus) & !is.na(tax$Species.only),
                      paste(tax$Genus, tax$Species.only),
                      ifelse(!is.na(tax$Genus) & is.na(tax$Species.only),
                             paste(tax$Genus, "spp."),
                             NA))
tax$Species.only <- NULL
tax <- as.matrix.data.frame(tax)
tax <- phyloseq::tax_table(tax)
phyloseq::tax_table(physeq_PC) <- tax; rm(tax)

#### End ####

#### Relative abundance ####
# Normalize number of reads in each sample using median sequencing depth

total = median(sample_sums(physeq_PC))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_PC_mednorm = transform_sample_counts(physeq_PC, standf)

# Transform to relative abundance
physeq_PC_re = transform_sample_counts(physeq_PC_mednorm, function(x){x / sum(x)})

#### End ####

#### Plot them settings ####
text_size = 30
thm <- theme(text = element_text(size = text_size),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "bottom")

thm.x <- theme(text = element_text(size = text_size),
               axis.text.x = element_text(angle = -30, hjust = 0),
               legend.position = "bottom")

thm.x2 <- theme(text = element_text(size = text_size),
               #axis.text.x = element_text(hjust = 0),
               legend.position = "bottom")
#### End ####

#### Positive controls plots ####

#### TSBC ####

# Select samples
physeq_TSBC <- subset_samples(physeq_PC, sample_data(physeq_PC)$Study_Nr %in% c("TSBC-3", "TSBC-4"))
physeq_TSBC <- prune_taxa(taxa_sums(physeq_TSBC) != 0, physeq_TSBC)
sample_data(physeq_TSBC)$Study_Nr <- c("TSBC-A", "TSBC-B")
physeq_TSBC_re <- subset_samples(physeq_PC_re, sample_data(physeq_PC_re)$Study_Nr %in% c("TSBC-3", "TSBC-4"))
otu_table(physeq_TSBC_re) <- otu_table(physeq_TSBC_re) *100 # Transform to percentage
physeq_TSBC_re <- prune_taxa(taxa_sums(physeq_TSBC_re) != 0, physeq_TSBC_re)
sample_data(physeq_TSBC_re)$Study_Nr <- c("TSBC-A", "TSBC-B")

# Select top n ASVs
n = 15 # adjust number to wished top ones
top_TSBC <- names(sort(taxa_sums(physeq_TSBC), decreasing=TRUE))[1:n]  

physeq_TSBC_top <- prune_taxa(top_TSBC, physeq_TSBC)
physeq_TSBC_re_top <- prune_taxa(top_TSBC, physeq_TSBC_re)

# Absolute abundance

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + theme(text = element_text(size = text_size),
                                                                                        axis.text.x = element_text(size = 20),
                                                                                        axis.ticks.x = element_blank(),
                                                                                        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC_abs_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC_abs_genus.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC_abs_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC_abs_species.svg"), width = 15, height = 12, dpi = 300)

#Relative Abundance

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))


ggsave(file.path(qc.dir, "TSBC_rel_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC_rel_genus.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "TSBC_rel_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC_rel_species.svg"), width = 15, height = 12, dpi = 300)

#### TSBC-1 ####

# Select samples
physeq_TSBC <- subset_samples(physeq_PC, sample_data(physeq_PC)$Study_Nr %in% c("TSBC-1"))
physeq_TSBC <- prune_taxa(taxa_sums(physeq_TSBC) != 0, physeq_TSBC)
physeq_TSBC_re <- subset_samples(physeq_PC_re, sample_data(physeq_PC_re)$Study_Nr %in% c("TSBC-1"))
otu_table(physeq_TSBC_re) <- otu_table(physeq_TSBC_re) *100 # Transform to percentage
physeq_TSBC_re <- prune_taxa(taxa_sums(physeq_TSBC_re) != 0, physeq_TSBC_re)

# Select top n ASVs
n = 30 # adjust number to wished top ones
top_TSBC <- names(sort(taxa_sums(physeq_TSBC), decreasing=TRUE))[1:n]

physeq_TSBC_top <- prune_taxa(top_TSBC, physeq_TSBC)
physeq_TSBC_re_top <- prune_taxa(top_TSBC, physeq_TSBC_re)

# Absolute abundance

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + theme(text = element_text(size = text_size),
                                                                                        axis.text.x = element_text(size = 20),
                                                                                        axis.ticks.x = element_blank(),
                                                                                        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-1_abs_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-1_abs_genus.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-1_abs_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-1_abs_species.svg"), width = 15, height = 12, dpi = 300)

#Relative Abundance

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-1_rel_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-1_rel_genus.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "TSBC-1_rel_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-1_rel_species.svg"), width = 15, height = 12, dpi = 300)

#### TSBC-2 ####

# Select samples
physeq_TSBC <- subset_samples(physeq_PC, sample_data(physeq_PC)$Study_Nr %in% c("TSBC-2"))
physeq_TSBC <- prune_taxa(taxa_sums(physeq_TSBC) != 0, physeq_TSBC)
physeq_TSBC_re <- subset_samples(physeq_PC_re, sample_data(physeq_PC_re)$Study_Nr %in% c("TSBC-2"))
otu_table(physeq_TSBC_re) <- otu_table(physeq_TSBC_re) *100 # Transform to percentage
physeq_TSBC_re <- prune_taxa(taxa_sums(physeq_TSBC_re) != 0, physeq_TSBC_re)

# Select top n ASVs
n = 30 # adjust number to wished top ones
top_TSBC <- names(sort(taxa_sums(physeq_TSBC), decreasing=TRUE))[1:n] 

physeq_TSBC_top <- prune_taxa(top_TSBC, physeq_TSBC)
physeq_TSBC_re_top <- prune_taxa(top_TSBC, physeq_TSBC_re)

# Absolute abundance

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + theme(text = element_text(size = text_size),
                                                                                        axis.text.x = element_text(size = 20),
                                                                                        axis.ticks.x = element_blank(),
                                                                                        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-2_abs_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-2_abs_genus.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-2_abs_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-2_abs_species.svg"), width = 15, height = 12, dpi = 300)

#Relative Abundance

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-2_rel_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-2_rel_genus.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "TSBC-2_rel_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-2_rel_species.svg"), width = 15, height = 12, dpi = 300)


#### TSBC-3 ####

# Select samples
physeq_TSBC <- subset_samples(physeq_PC, sample_data(physeq_PC)$Study_Nr %in% c("TSBC-3"))
physeq_TSBC <- prune_taxa(taxa_sums(physeq_TSBC) != 0, physeq_TSBC)
sample_data(physeq_TSBC)$Study_Nr <- "TSBC-A"
physeq_TSBC_re <- subset_samples(physeq_PC_re, sample_data(physeq_PC_re)$Study_Nr %in% c("TSBC-3"))
otu_table(physeq_TSBC_re) <- otu_table(physeq_TSBC_re) *100 # Transform to percentage
physeq_TSBC_re <- prune_taxa(taxa_sums(physeq_TSBC_re) != 0, physeq_TSBC_re)
sample_data(physeq_TSBC_re)$Study_Nr <- "TSBC-A"

# Select top n ASVs
n = 5 # adjust number to wished top ones
top_TSBC <- names(sort(taxa_sums(physeq_TSBC), decreasing=TRUE))[1:n]

physeq_TSBC_top <- prune_taxa(top_TSBC, physeq_TSBC)
physeq_TSBC_re_top <- prune_taxa(top_TSBC, physeq_TSBC_re)

# Absolute abundance

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm + theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-3_abs_genus.tiff"), width = 6, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-3_abs_genus.svg"), width = 6, height = 10, dpi = 300)

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm + theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-3_abs_species.tiff"), width = 12, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-3_abs_species.svg"), width = 10, height = 12, dpi = 300)

#Relative Abundance

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm + theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-3_rel_genus.tiff"), width = 6, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-3_rel_genus.svg"), width = 6, height = 10, dpi = 300)

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm + theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "TSBC-3_rel_species.tiff"), width = 12, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-3_rel_species.svg"), width = 10, height = 12, dpi = 300)

#### TSBC-4 ####

# Select samples
physeq_TSBC <- subset_samples(physeq_PC, sample_data(physeq_PC)$Study_Nr %in% c("TSBC-4"))
physeq_TSBC <- prune_taxa(taxa_sums(physeq_TSBC) != 0, physeq_TSBC)
sample_data(physeq_TSBC)$Study_Nr <- "TSBC-B"
physeq_TSBC_re <- subset_samples(physeq_PC_re, sample_data(physeq_PC_re)$Study_Nr %in% c("TSBC-4"))
otu_table(physeq_TSBC_re) <- otu_table(physeq_TSBC_re) *100 # Transform to percentage
physeq_TSBC_re <- prune_taxa(taxa_sums(physeq_TSBC_re) != 0, physeq_TSBC_re)
sample_data(physeq_TSBC_re)$Study_Nr <- "TSBC-B"

# Select top n ASVs
n = 5 # adjust number to wished top ones
top_TSBC <- names(sort(taxa_sums(physeq_TSBC), decreasing=TRUE))[1:n]

physeq_TSBC_top <- prune_taxa(top_TSBC, physeq_TSBC)
physeq_TSBC_re_top <- prune_taxa(top_TSBC, physeq_TSBC_re)

# Absolute abundance

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm + theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-4_abs_genus.tiff"), width = 6, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-4_abs_genus.svg"), width = 6, height = 10, dpi = 300)

plot_bar(physeq_TSBC_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm + theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "TSBC-4_abs_species.tiff"), width = 12, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-4_abs_species.svg"), width = 8, height = 12, dpi = 300)

#Relative Abundance

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm + theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "TSBC-4_rel_genus.tiff"), width = 6, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-4_rel_genus.svg"), width = 6, height = 10, dpi = 300)

plot_bar(physeq_TSBC_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm + theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "TSBC-4_rel_species.tiff"), width = 12, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "TSBC-4_rel_species.svg"), width = 8, height = 12, dpi = 300)

#### PC-Fecal ####

# Select samples
physeq_Fecal <- subset_samples(physeq_PC, sample_data(physeq_PC)$Study_Nr == "PC-Fecal")
physeq_Fecal <- prune_taxa(taxa_sums(physeq_Fecal) != 0, physeq_Fecal)
physeq_Fecal_re <- subset_samples(physeq_PC_re, sample_data(physeq_PC_re)$Study_Nr == "PC-Fecal")
otu_table(physeq_Fecal_re) <- otu_table(physeq_Fecal_re) *100 # Transform to percentage
physeq_Fecal_re <- prune_taxa(taxa_sums(physeq_Fecal_re) != 0, physeq_Fecal_re)

# Select top n ASVs
n = 30 # adjust number to wished top ones
top_Fecal <- names(sort(taxa_sums(physeq_Fecal), decreasing=TRUE))[1:n]

physeq_Fecal_top <- prune_taxa(top_Fecal, physeq_Fecal)
physeq_Fecal_re_top <- prune_taxa(top_Fecal, physeq_Fecal_re)

# Absolute abundance

plot_bar(physeq_Fecal, y =  "Abundance", title = "Abundance based on Genus", fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "Fecal_abs_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Fecal_abs_genus.svg"), width = 8, height = 10, dpi = 300)

#Relative Abundance

plot_bar(physeq_Fecal_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))


ggsave(file.path(qc.dir, "Fecal_rel_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Fecal_rel_genus.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_Fecal_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "Fecal_rel_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "Fecal_rel_species.svg"), width = 15, height = 12, dpi = 300)

#### PC-Gut ####

# Select samples
physeq_Gut <- subset_samples(physeq_PC, sample_data(physeq_PC)$Study_Nr == "PC-Gut")
physeq_Gut <- prune_taxa(taxa_sums(physeq_Gut) != 0, physeq_Gut)
physeq_Gut_re <- subset_samples(physeq_PC_re, sample_data(physeq_PC_re)$Study_Nr == "PC-Gut")
otu_table(physeq_Gut_re) <- otu_table(physeq_Gut_re) *100 # Transform to percentage
physeq_Gut_re <- prune_taxa(taxa_sums(physeq_Gut_re) != 0, physeq_Gut_re)

# Select top n ASVs
n = 30 # adjust number to wished top ones
top_Gut <- names(sort(taxa_sums(physeq_Gut), decreasing=TRUE))[1:n] 

physeq_Gut_top <- prune_taxa(top_Gut, physeq_Gut)
physeq_Gut_re_top <- prune_taxa(top_Gut, physeq_Gut_re)

# Absolute abundance

plot_bar(physeq_Gut_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus")) + facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "Gut_abs_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Gut_abs_genus.svg"), width = 8, height = 10, dpi = 300)

#Relative Abundance

plot_bar(physeq_Gut_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))


ggsave(file.path(qc.dir, "Gut_rel_genus.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Gut_rel_genus.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_Gut_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "Gut_rel_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "Gut_rel_species.svg"), width = 15, height = 12, dpi = 300)

#### Mock ####

# Select samples
physeq_Mock <- subset_samples(physeq_PC, sample_data(physeq_PC)$Study_Nr == "PC-Mock")
physeq_Mock <- prune_taxa(taxa_sums(physeq_Mock) != 0, physeq_Mock)
physeq_Mock_re <- subset_samples(physeq_PC_re, sample_data(physeq_PC_re)$Study_Nr == "PC-Mock")
otu_table(physeq_Mock_re) <- otu_table(physeq_Mock_re) *100 # Transform to percentage
physeq_Mock_re <- prune_taxa(taxa_sums(physeq_Mock_re) != 0, physeq_Mock_re)

# Select top n ASVs
n = 30 # adjust number to wished top ones
top_Mock <- names(sort(taxa_sums(physeq_Mock), decreasing=TRUE))[1:n]  

physeq_Mock_top <- prune_taxa(top_Mock, physeq_Mock)
physeq_Mock_re_top <- prune_taxa(top_Mock, physeq_Mock_re)

# Absolute abundance
plot_bar(physeq_Mock_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "Mock_abs_genus_SILVA.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Mock_abs_genus_SILVA.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_Mock_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

#Relative Abundance

plot_bar(physeq_Mock_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))


ggsave(file.path(qc.dir, "Mock_rel_genus_SILVA.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Mock_rel_genus_SILVA.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_Mock_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "Mock_rel_species_SILVA.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "Mock_rel_species_SILVA.svg"), width = 15, height = 12, dpi = 300)

# Expected Mock
# Specify the abundance for the expected Mock Community
mock_asv <- otu_table(as.data.frame(c(Bacillus_subtilis = 17, 
                Enterococcus_faecalis = 10, 
                Escherichia_coli = 10, 
                Lactobacillus_fermentum = 18, 
                Listeria_monocytogenes = 14, 
                Salmonella_enterica = 10, 
                Pseudomonas_aeruginosa = 4, 
                Staphylococcus_aureus = 16)), taxa_are_rows = TRUE)

mock_asv <- c(Bacillus_subtilis = 17,
              Enterococcus_faecalis = 10,
              Escherichia_coli = 10,
              Lactobacillus_fermentum = 18,
              Listeria_monocytogenes = 14,
              Salmonella_enterica = 10,
              Pseudomonas_aeruginosa = 4,
              Staphylococcus_aureus = 16) %>% as.data.frame(.) %>% otu_table(., taxa_are_rows = TRUE)

mock_tax <- tax_table(rbind(Bacillus_subtilis = c("","","","","", "Bacillus", "Bacillus subtilis"), 
              Enterococcus_faecalis = c("","","","","", "Enterococcus", "Enterococcus faecalis"), 
              Escherichia_coli = c("","","","","", "Escherichia", "Escherichia coli"), 
              Lactobacillus_fermentum = c("","","","","", "Lactobacillus", "Lactobacillus fermentum"), 
              Listeria_monocytogenes = c("","","","","", "Listeria", "Listeria monocytogenes"), 
              Salmonella_enterica = c("","","","","", "Salmonella", "Salmonella enterica"), 
              Pseudomonas_aeruginosa = c("","","","","", "Pseudomonas", "Pseudomonas aeruginosa"), 
              Staphylococcus_aureus = c("","","","","", "Staphylococcus", "Staphylococcus aureus")))
colnames(mock_tax) <- c("Kingdom","Phylum","Phylum","Order","Family", "Genus", "Species")

# Create a new sample data frame with sample metadata
new_sample_data <- data.frame(Study_Nr = "Expected_Mock_Community")
rownames(new_sample_data) <- "."  # Set row name to match sample ID
new_sample_data <- sample_data(new_sample_data)

# Combine the new OTU table and sample data into a new phyloseq object
physeq_expected_mock <- phyloseq(mock_asv, new_sample_data, mock_tax)

physeq_combined <- merge_phyloseq(physeq_Mock_re_top, physeq_expected_mock)

#Relative Abundance

plot_bar(physeq_expected_mock, y =  "Abundance", title = paste("Abundance based\non Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))


ggsave(file.path(qc.dir, "Mock_exp_genus.tiff"), width = 6, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Mock_exp_genus.svg"), width = 6, height = 10, dpi = 300)

plot_bar(physeq_expected_mock, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "Mock_exp_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "Mock_exp_species.svg"), width = 15, height = 12, dpi = 300)

sample_data(physeq_combined)$Study_Nr <- c(rep("Sequenced Mock Community", 7), "Expected Mock Community")

plot_bar(physeq_combined, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "Mock_exp_seq_genus.tiff"), width = 11, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Mock_exp_seq_genus.svg"), width = 11, height = 10, dpi = 300)

plot_bar(physeq_combined, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "Mock_exp_seq_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "Mock_exp_seq_species.svg"), width = 15, height = 12, dpi = 300)
#### End ####

#### Mock -different databases
# Select samples
physeq_Mock <- subset_samples(physeq_PC, sample_data(physeq_PC)$Study_Nr == "PC-Mock")
physeq_Mock <- prune_taxa(taxa_sums(physeq_Mock) != 0, physeq_Mock)
physeq_Mock_re <- subset_samples(physeq_PC_re, sample_data(physeq_PC_re)$Study_Nr == "PC-Mock")
otu_table(physeq_Mock_re) <- otu_table(physeq_Mock_re) *100 # Transform to percentage
physeq_Mock_re <- prune_taxa(taxa_sums(physeq_Mock_re) != 0, physeq_Mock_re)

# Select top n ASVs
n = 30 # adjust number to wished top ones
top_Mock <- names(sort(taxa_sums(physeq_Mock), decreasing=TRUE))[1:n]  

physeq_Mock_top_SILVA <- prune_taxa(top_Mock, physeq_Mock)
physeq_Mock_re_top_SILVA <- prune_taxa(top_Mock, physeq_Mock_re)

# Absolute abundance
plot_bar(physeq_Mock_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "Mock_abs_genus_SILVA.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Mock_abs_genus_SILVA.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_Mock_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

#Relative Abundance

plot_bar(physeq_Mock_re_top_GTDB, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))


ggsave(file.path(qc.dir, "Mock_rel_genus_SILVA.tiff"), width = 8, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Mock_rel_genus_SILVA.svg"), width = 8, height = 10, dpi = 300)

plot_bar(physeq_Mock_re_top, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "Mock_rel_species_SILVA.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "Mock_rel_species_SILVA.svg"), width = 15, height = 12, dpi = 300)

# Expected Mock
# Specify the abundance for the expected Mock Community
mock_asv <- otu_table(as.data.frame(c(Bacillus_subtilis = 17, 
                                      Enterococcus_faecalis = 10, 
                                      Escherichia_coli = 10, 
                                      Lactobacillus_fermentum = 18, 
                                      Listeria_monocytogenes = 14, 
                                      Salmonella_enterica = 10, 
                                      Pseudomonas_aeruginosa = 4, 
                                      Staphylococcus_aureus = 16)), taxa_are_rows = TRUE)

mock_asv <- c(Bacillus_subtilis = 17,
              Enterococcus_faecalis = 10,
              Escherichia_coli = 10,
              Lactobacillus_fermentum = 18,
              Listeria_monocytogenes = 14,
              Salmonella_enterica = 10,
              Pseudomonas_aeruginosa = 4,
              Staphylococcus_aureus = 16) %>% as.data.frame(.) %>% otu_table(., taxa_are_rows = TRUE)

mock_tax <- tax_table(rbind(Bacillus_subtilis = c("","","","","", "Bacillus", "Bacillus subtilis"), 
                            Enterococcus_faecalis = c("","","","","", "Enterococcus", "Enterococcus faecalis"), 
                            Escherichia_coli = c("","","","","", "Escherichia", "Escherichia coli"), 
                            Lactobacillus_fermentum = c("","","","","", "Lactobacillus", "Lactobacillus fermentum"), 
                            Listeria_monocytogenes = c("","","","","", "Listeria", "Listeria monocytogenes"), 
                            Salmonella_enterica = c("","","","","", "Salmonella", "Salmonella enterica"), 
                            Pseudomonas_aeruginosa = c("","","","","", "Pseudomonas", "Pseudomonas aeruginosa"), 
                            Staphylococcus_aureus = c("","","","","", "Staphylococcus", "Staphylococcus aureus")))
colnames(mock_tax) <- c("Kingdom","Phylum","Phylum","Order","Family", "Genus", "Species")

# Create a new sample data frame with sample metadata
new_sample_data <- data.frame(Study_Nr = "Expected_Mock_Community")
rownames(new_sample_data) <- "."  # Set row name to match sample ID
new_sample_data <- sample_data(new_sample_data)

# Combine the new OTU table and sample data into a new phyloseq object
physeq_expected_mock <- phyloseq(mock_asv, new_sample_data, mock_tax)

sample_data(physeq_expected_mock)$Study_Nr <- "Expected"
sample_data(physeq_Mock_re_top_SILVA)$Study_Nr <- "SILVA database"
sample_names(physeq_Mock_re_top_SILVA) <- paste(sample_names(physeq_Mock_re_top_SILVA), "SILVA", sep = "_")
sample_data(physeq_Mock_re_top_GTDB)$Study_Nr <- "GTDB database"
sample_names(physeq_Mock_re_top_GTDB) <- paste(sample_names(physeq_Mock_re_top_GTDB), "GTDB", sep = "_")
sample_data(physeq_Mock_re_top_eHOMD)$Study_Nr <- "eHOMD database"
sample_names(physeq_Mock_re_top_eHOMD) <- paste(sample_names(physeq_Mock_re_top_eHOMD), "eHOMD", sep = "_")

taxa_names(physeq_Mock_re_top_SILVA) <- paste(taxa_names(physeq_Mock_re_top_SILVA), "SILVA", sep = "_")
taxa_names(physeq_Mock_re_top_GTDB) <- paste(taxa_names(physeq_Mock_re_top_GTDB), "GTDB", sep = "_")
taxa_names(physeq_Mock_re_top_eHOMD) <- paste(taxa_names(physeq_Mock_re_top_eHOMD), "eHOMD", sep = "_")

physeq_combined <- merge_phyloseq(physeq_Mock_re_top_SILVA, physeq_Mock_re_top_GTDB, physeq_Mock_re_top_eHOMD, physeq_expected_mock)

#Relative Abundance

plot_bar(physeq_expected_mock, y =  "Abundance", title = paste("Abundance based\non Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))


ggsave(file.path(qc.dir, "Mock_exp_genus.tiff"), width = 6, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Mock_exp_genus.svg"), width = 6, height = 10, dpi = 300)

plot_bar(physeq_expected_mock, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "Mock_exp_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "Mock_exp_species.svg"), width = 15, height = 12, dpi = 300)

sample_data(physeq_combined)$Study_Nr <- c(rep("Sequenced Mock Community", 7), "Expected Mock Community")

plot_bar(physeq_combined, y =  "Abundance", title = paste("Abundance from top", n, "ASVs\nbased on Genus"), fill="Genus")+ facet_wrap(~Study_Nr, scales="free_x", nrow = 1) + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 2, title = NULL),
         color = guide_legend(ncol = 2, title = NULL))

ggsave(file.path(qc.dir, "Mock_exp_seq_genus.tiff"), width = 11, height = 10, dpi = 300)
ggsave(file.path(qc.dir, "Mock_exp_seq_genus.svg"), width = 11, height = 10, dpi = 300)

plot_bar(physeq_combined, y =  "Abundance", title = paste("Abundance from top", n, "ASVs based on Species"), fill="Species")+ facet_wrap(~Study_Nr, scales="free_x") + 
  geom_bar(aes(color=Species, fill=Species), stat = "Identity", position = "stack") + thm +
  guides(fill = guide_legend(ncol = 1, title = NULL),
         color = guide_legend(ncol = 1, title = NULL))

ggsave(file.path(qc.dir, "Mock_exp_seq_species.tiff"), width = 15, height = 15, dpi = 300)
ggsave(file.path(qc.dir, "Mock_exp_seq_species.svg"), width = 15, height = 12, dpi = 300)
#### End ####
