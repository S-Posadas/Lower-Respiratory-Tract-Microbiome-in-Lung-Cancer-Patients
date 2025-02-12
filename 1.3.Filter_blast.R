library(dplyr)
library(openxlsx)

res.dir <- "Results_trimmed_cutadapt3_20241216"

#### Format blast table to get a single taxa per ASV blast 16S ####
# Read alignment table from blast and arrange column names
blastu <- read.csv(file.path(res.dir,"2.blast/tax_tab_95_10align_16S"), row.names=NULL, sep="")
colnames(blastu) <- c(colnames(blastu)[-1],"Species")
colnames(blastu)[colnames(blastu) == "stitle"] <- "Genus"
total_asv <- as.numeric(gsub("ASV_", "", blastu[length(rownames(blastu)),"seqid"]))

# Make sure pident and evalue are numeric
blastu$pident = as.numeric(blastu$pident)
blastu$evalue = as.numeric(blastu$evalue)
head(blastu)
dim(blastu)

# Filter alingments with percentage of identical matches higher than n 
# We already used 95% in blast, but just in case we want to be more strict
# In our case we are not really filtering
blastu = blastu[blastu$pident >= 99,] 
dim(blastu)

# Check how many unique genera were obtained
unique(blastu$Genus)

# Check how many unique genera were obtained
length(unique(blastu$Genus))

# Keep alingments with the highest percentage of identical matches #and lowest e value
# For ASVs with alignment to a single species, keep only one row
filt = list()
for(asv in unique(blastu$seqid)){
  ASV_n = blastu[blastu$seqid == asv,]
  ASV_n = ASV_n[ASV_n$pident == max(ASV_n$pident),] # highest identity
  ASV_n = ASV_n[ASV_n$evalue == min(ASV_n$evalue),] # lowest e value
  if(length(unique(ASV_n$Genus)) == 1 & length(unique(ASV_n$Species)) == 1){
    ASV_n = ASV_n[1,]
  }
  filt[[asv]] = ASV_n
}
filtu <- as.data.frame(do.call(rbind, filt)) 
dim(filtu)
min(filtu$pident)

# Check if Genus is not unique within each ASV
non_unique_genusu <- filtu %>%
  group_by(seqid) %>%
  filter(n_distinct(Genus) > 1)
length(unique(non_unique_genusu$seqid)) 

filtu[filtu$seqid %in% unique(non_unique_genusu$seqid), c("Genus", "Species")] <- NA

# Check if Species is not unique within each ASV
non_unique_speciesu <- filtu %>%
  group_by(seqid) %>%
  filter(n_distinct(Species) > 1)
length(unique(non_unique_speciesu$seqid)) 

filtu[filtu$seqid %in% unique(non_unique_speciesu$seqid), c("Species")] <- NA

# Keep only one row per ASV
length(unique(filtu$seqid)) # Check how many unique ASVs are there

final = list()
for(asv in unique(filtu$seqid)){
  ASV_n = filtu[filtu$seqid == asv,]
  ASV_n = ASV_n[1,]
  final[[asv]] = ASV_n
}
finalu <- as.data.frame(do.call(rbind, final)) 
dim(finalu)

min(as.numeric(finalu$pident))

# Add unassigned ASVs 
finalu$asv <- gsub("ASV_", "", finalu$seqid)
asvNA <- paste("ASV", c(1:total_asv)[!1:total_asv %in% finalu$asv], sep = "_")
asvNA <- cbind(asvNA, matrix(NA, nrow = length(asvNA), ncol = 10))
colnames(asvNA) <- colnames(finalu)
finalu <- as.data.frame(rbind(finalu, asvNA))
finalu$asv <- as.numeric(gsub("ASV_", "", finalu$seqid))
finalu <- finalu[order(finalu$asv),]
finalu$asv <- NULL

#Save data in desired format
write.xlsx(finalu, file.path(res.dir,"2.blast/blast_final_99_16S.xlsx"))

#### End ####

#### Add blast species to unassigned SILVA species ####

# Read taxa table from SILVA
tax_tab_SILVA <- read.xlsx(file.path(res.dir,"2.blast/tax_tab_SILVA.xlsx"), rowNames = T)
dim(tax_tab_SILVA)
colnames(tax_tab_SILVA) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species.bayes", "Species.100")

tax_tab_GTDB <- read.xlsx(file.path(res.dir,"2.blast/tax_tab_GTDB.xlsx"), rowNames = T)
dim(tax_tab_GTDB)
colnames(tax_tab_GTDB) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species.bayes", "Species.100")

tax_tab_eHOMD <- read.xlsx(file.path(res.dir,"2.blast/tax_tab_eHOMD.xlsx"), rowNames = T)
dim(tax_tab_eHOMD)

# Read taxa table from blast
blast <- read.xlsx(file.path(res.dir, "2.blast/blast_final_99_16S.xlsx"), rowNames = T)
dim(blast)

# Change blast taxonomy to SILVA taxonomy
# blast[!is.na(blast$Genus) & blast$Genus == "Schaalia" & blast$Species == "odontolytica", "Species"] = "odontolyticus"
# blast[!is.na(blast$Genus) & blast$Genus == "Schaalia", "Genus"] = "Actinomyces"
# blast[!is.na(blast$Genus) & blast$Genus == "Caldimonas","Genus"] = "Schlegelella"
# blast[!is.na(blast$Genus) & blast$Genus == "Ligilactobacillus","Genus"] = "Lactobacillus"
# blast[!is.na(blast$Genus) & blast$Genus == "Limosilactobacillus","Genus"] = "Lactobacillus"
# blast[!is.na(blast$Genus) & blast$Genus == "Lacticaseibacillus","Genus"] = "Lactobacillus"
# blast[!is.na(blast$Genus) & blast$Genus == "Levilactobacillus","Genus"] = "Lactobacillus"
# tax_tab_SILVA[grep("Clostridium sensu stricto", tax_tab_SILVA$Genus),"Genus"] <- "Clostridium"
# blast[!is.na(blast$Genus) & blast$Genus == "Lancefieldella" & blast$Species %in% c("rimae", "parvula"), "Genus"] = "Atopobium"
# blast[!is.na(blast$Genus) & blast$Genus == "Atopobium" & blast$Species == "parvula", "Species"] = "parvulum"
blast[!is.na(blast$Genus) & blast$Genus == "Mycolicibacterium", "Genus"] = "Mycobacterium"
blast[!is.na(blast$Genus) & blast$Genus == "Phocaeicola" & blast$Species %in% c("vulgatus", "dorei", "massiliensis"), "Genus"] = "Bacteroides"
blast[!is.na(blast$Genus) & blast$Genus == "Hallella" & blast$Species %in% c("multisaccharivorax", "mizrahii", "colorans"), "Genus"] = "Prevotella"
blast[!is.na(blast$Genus) & blast$Genus == "Ureibacillus" & blast$Species %in% c("massiliensis", "chungkukjangi"), "Genus"] = "Lysinibacillus"
blast[!is.na(blast$Genus) & blast$Genus == "[Ruminococcus]" & blast$Species %in% c("faecis", "torques"), "Genus"] = "Mediterraneibacter"
blast[!is.na(blast$Genus) & blast$Genus == "Roseateles" & blast$Species %in% c("oligotrophus"), "Genus"] = "Paucibacter"
blast[!is.na(blast$Genus) & blast$Genus == "Roseateles" & blast$Species %in% c("asaccharophilus"), "Genus"] = "Kinneretia"
blast[!is.na(blast$Genus) & blast$Genus == "Kinneretia" & blast$Species == "asaccharophilus", "Species"] = "asaccharophila"
# blast[!is.na(blast$Genus) & blast$Genus == "Segatella", "Genus"] = "Prevotella"
# blast[!is.na(blast$Genus) & blast$Genus == "Hoylesella", "Genus"] = "Prevotella"
blast[!is.na(blast$Genus) & blast$Genus == "Macrococcoides", "Genus"] = "Macrococcus"
blast[!is.na(blast$Genus) & blast$Genus == "Metamycoplasma", "Genus"] = "Mycoplasma"

# blast["ASV_266",]
#blast[!is.na(tax_tab_SILVA$Genus) & tax_tab_SILVA$Genus == "Phocaeicola",]
# unique(blast[!is.na(blast$Genus) & blast$Genus == "Phocaeicola","Species"])
# tax_tab_SILVA[!is.na(blast$Genus) & blast$Genus == "Caloramator",]
#tax_tab_SILVA["ASV_12072",]

# Explore discrepant genus

all(rownames(tax_tab_SILVA) == rownames(blast))
samegenus <- rownames(blast[(is.na(tax_tab_SILVA$Species.100) | grepl("/", tax_tab_SILVA$Species.100))
                            & !is.na(blast$Species)
                            & !is.na(tax_tab_SILVA$Genus)
                            & blast$Genus == tax_tab_SILVA$Genus,])

#tax_tab_SILVA[samegenus,"Species"] <- blast[samegenus, "Species"]
length(samegenus)

othergenus <- rownames(blast[!is.na(blast$Genus)
                             & !is.na(tax_tab_SILVA$Genus)
                             & blast$Genus != tax_tab_SILVA$Genus,])

length(othergenus)
unique(blast[othergenus,"Genus"])
unique(tax_tab_SILVA[othergenus,"Genus"])

merg <- merge(tax_tab_SILVA[othergenus,c("Genus", "Species.bayes", "Species.100")], blast[othergenus,c("Genus", "Species")], by = 0)

write.xlsx(merg, file.path(res.dir, "2.blast/discrepant_SILVA_blast_genus.xlsx"), rowNames = T)

head(tax_tab_SILVA)
tax_tab_SILVA$Species <- tax_tab_SILVA$Species.100
tax_tab_SILVA[samegenus,"Species"] <- blast[samegenus, "Species"]

# Save final taxa table
write.xlsx(tax_tab_SILVA, file.path(res.dir, "2.blast/tax_tab_SILVA_blast_species.xlsx"), rowNames = T)

#### End ####

#### Explore other databases ####
unique(tax_tab_SILVA$Phylum)
unique(tax_tab_eHOMD$Phylum)

otherphylum <- rownames(tax_tab_eHOMD[!is.na(tax_tab_eHOMD$Phylum)
                                     & !is.na(tax_tab_SILVA$Phylum)
                                     & tax_tab_eHOMD$Phylum != tax_tab_SILVA$Phylum,])
length(otherphylum)
tax_tab_eHOMD[!is.na(tax_tab_eHOMD$Phylum)
              & is.na(tax_tab_SILVA$Phylum),]
# Explore discrepant genus

all(rownames(tax_tab_SILVA) == rownames(tax_tab_eHOMD))
samegenus <- rownames(tax_tab_eHOMD[(is.na(tax_tab_SILVA$Species.100) | grepl("/", tax_tab_SILVA$Species.100))
                            & !is.na(tax_tab_eHOMD$Species)
                            & !is.na(tax_tab_SILVA$Genus)
                            & tax_tab_eHOMD$Genus == tax_tab_SILVA$Genus,])

#tax_tab_SILVA[samegenus,"Species"] <- blast[samegenus, "Species"]
length(samegenus)

othergenus <- rownames(tax_tab_eHOMD[!is.na(tax_tab_eHOMD$Genus)
                             & !is.na(tax_tab_SILVA$Genus)
                             & tax_tab_eHOMD$Genus != tax_tab_SILVA$Genus,])

length(othergenus)
unique(tax_tab_eHOMD[othergenus,"Genus"])
unique(tax_tab_SILVA[othergenus,"Genus"])

merg <- merge(tax_tab_SILVA[othergenus,], tax_tab_eHOMD[othergenus,], by = 0)

write.xlsx(merg, file.path(res.dir, "2.blast/discrepant_SILVA_eHOMD_genus.xlsx"), rowNames = T)

tax_tab_SILVA$Species <- tax_tab_SILVA$Species.100
tax_tab_SILVA[samegenus,"Species"] <- blast[samegenus, "Species"]

write.xlsx(tax_tab_SILVA, file.path(res.dir, "2.blast/tax_tab_SILVA_SILVA_blast_species.xlsx"), rowNames = T)
