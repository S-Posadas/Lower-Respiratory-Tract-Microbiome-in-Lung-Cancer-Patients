# Beta diversity, taxonomy analysis

#### Package Setup ####

Sys.setenv(language = "EN")

# Core utilities
library(plyr)
library(reshape)  # Load before tidyverse (for Rhea)

# Tidyverse (load dplyr AFTER plyr)
library(dplyr)
library(tidyr)    # Better than reshape for most tasks (for Rhea)
library(ggplot2)

# Other packages
library(phyloseq)
library(ggordiplots)
library(microbiomeMarker)
library(Maaslin2)
library(ggrepel)
library(gridExtra)
library(openxlsx)
library(cowplot)
library(SummarizedExperiment)
library(lefser)
source("Results_trimmed_cutadapt3_20241216/2025/Scripts/my_lefser.R")

theme_set(theme_bw())

#### End ####

#### Directory Setup ####

res.dir <- "Results"
load(file.path(res.dir, "RData", "6.phyloseq.corrected.RData"))

# Create main output directories (Set # in front to avoid creating unwanted directories)
b.stats <- file.path(res.dir, "Beta_stats")
b.plots <- file.path(res.dir, "Beta_plots")
bar.dir <- file.path(res.dir, "Most_abundant_plots")
maaslin.dir <- file.path(res.dir, "Maaslin2")
lefser.dir <- file.path(res.dir, "LEfSer")
rhea.dir <- file.path(res.dir, "Rhea")

for(dir in c("b.stats", "b.plots", "maaslin.dir", "bar.dir", "lefse.dir", "rhea.dir")){
  if(exists(dir)){                  # Create directories that were assigned to a variable
    dir <- get(dir) 
    dir.create(dir, recursive = T, showWarnings = FALSE)
  }else{
    message("Directory 'dir' not found.")
  }
}

#### End ####

#### Set theme ####
text_size = 30
thm <- theme(text = element_text(size = text_size),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "bottom")

thm.beta <- theme(text = element_text(size = text_size),
             legend.position = "bottom")

thm.x <- theme(text = element_text(size = text_size),
               axis.text.x = element_text(angle = -30, hjust = 0),
               legend.position = "bottom")

#### End ####

#### Functions - Beta diversity ####

beta_analysis <- function(physeq){
  # Calculate Bray-Curtis, weighted UniFrac and unweighted UniFrac distances
  set.seed(123)
  bray_dist <- phyloseq::distance(physeq, method = "bray")
  set.seed(123)
  wunifrac_dist <- phyloseq::distance(physeq, method = "wunifrac")
  set.seed(123)
  uunifrac_dist <- phyloseq::distance(physeq, method = "uunifrac")
  
  # Perform NMDS ordination
  set.seed(123)
  nmds_bray <- vegan::metaMDS(bray_dist, k = 2, trymax = 100)
  set.seed(123)
  nmds_wunifrac <- vegan::metaMDS(wunifrac_dist, k = 2, trymax = 100)
  set.seed(123)
  nmds_uunifrac <- vegan::metaMDS(uunifrac_dist, k = 2, trymax = 100)
  
  # Perform PCoA ordination
  pcoa_bray <- phyloseq::ordinate(physeq_beta, method = "PCoA", distance = "bray")
  pcoa_wunifrac <- phyloseq::ordinate(physeq_beta, method = "PCoA", distance = "wunifrac")
  pcoa_uunifrac <- phyloseq::ordinate(physeq_beta, method = "PCoA", distance = "uunifrac")
  
  return(list(bray_dist = bray_dist,
              wunifrac_dist = wunifrac_dist,
              uunifrac_dist = uunifrac_dist,
              nmds_bray = nmds_bray,
              nmds_wunifrac = nmds_wunifrac,
              nmds_uunifrac = nmds_uunifrac,
              pcoa_bray = pcoa_bray,
              pcoa_wunifrac = pcoa_wunifrac,
              pcoa_uunifrac = pcoa_uunifrac))
}

beta_plot <- function(physeq, beta, q, variable, col = NULL, prefix = "", legend = NULL){
  for(dist in c("bray", "wUF", "uwUF")){
    NP <- plot(phy = physeq, beta = beta, dist = dist, var = variable, palette = col, legend_title = legend)
    
    #ggsave(plot = NP, file.path(b.plots, paste0(params[[q]]$name, prefix, "_", dist, ".tiff")), width = 15, height = 8, dpi = 300)
    ggsave(plot = NP, file.path(b.plots, paste0(params[[q]]$name, prefix, "_", dist, ".svg")), width = 15, height = 8, dpi = 300)
    
  }
}

plot <- function(phy, beta, dist,var, palette = NULL, legend_title = NULL, nrow = NULL){
  if(dist == "bray"){NM = beta$nmds_bray ; PC = beta$pcoa_bray}
  else if(dist == "wUF"){NM = beta$nmds_wunifrac ; PC = beta$pcoa_wunifrac}
  else if(dist == "uwUF"){NM = beta$nmds_uunifrac ; PC = beta$pcoa_uunifrac}
  
  NMDS = plot_ordination(phy, NM, color=var) + geom_point(size=8, alpha=0.5) +
    stat_ellipse(type = "t") + thm.beta + #scale_color_manual(values = palette) +
    annotate("text", x = max(scores(NM)[,"NMDS1"]), y = max(scores(NM)[,"NMDS2"]),
             label = paste("Stress:", round(NM$stress, 3)), hjust = 1, vjust = 1, size = 7) 
  
  PCOA = plot_ordination(phy, PC, color=var) + geom_point(size=8, alpha=0.5) +
    stat_ellipse(type = "t") + thm.beta #+ scale_color_manual(values = palette)
  
  # Add color palette only if provided
  if (!is.null(palette)) {
    NMDS <- NMDS + scale_color_manual(values = palette)
    PCOA <- PCOA + scale_color_manual(values = palette)
  }
  
  if (!is.null(legend_title)) {
    NMDS <- NMDS + labs(color = legend_title)
    PCOA <- PCOA + labs(color = legend_title)
  } 
  
  if (!is.null(nrow)) {
    NMDS <- NMDS + guides(color = guide_legend(nrow = nrow))
    PCOA <- PCOA + guides(color = guide_legend(nrow = nrow))
  } 
  
  NP <- grid.arrange(NMDS + ggtitle("NMDS") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)),
                     PCOA + ggtitle("PCoA") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)),
                     nrow=1)
  
  NP <- grid.arrange(NP, get_plot_component(PCOA, pattern = "guide-box-bottom"), nrow= 2, heights = c(10,1))
  
}

#### Functions - Maaslin ####

# Clean taxonomy names
clean_tax_names_maaslin <- function(physeq) {
  tax <- as.data.frame(tax_table(physeq))
  tax[] <- lapply(tax, gsub, pattern = " ", replacement = "_")
  tax$Species[grep("_spp.", tax$Species)] <- NA
  tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis","Family"] <- paste(
    tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis","Order"], "Incertae_Sedis", sep = "_")
  tax_table(physeq) <- tax_table(as.matrix(tax))
  return(physeq)
}

# Function to agglomerate taxa
agglomerate_taxa <- function(physeq, level) {
  phy_agg <- tax_glom(physeq, level)
  taxa_names(phy_agg) <- tax_table(phy_agg)[, level]
  return(phy_agg)
}

# Volcano plot
volcano <- function(output_dir) {
  results <- read.csv(file.path(output_dir, "all_results.tsv"), sep = "\t")
  ggplot(results, aes(x = coef, y = -log10(qval))) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), col = "red") +
    geom_vline(xintercept = 0, col = "gray") +
    annotate("text", x = min(results$coef, na.rm = TRUE), y = 1.5, 
             label = "q-value = 0.05", col = "red", size = 6) +
    geom_text_repel(aes(label = feature), size = 6, color = "black") + 
    theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20))
  ggsave(file.path(output_dir, "all_results.svg"), width = 8, height = 8)
  ggsave(file.path(output_dir, "all_results.tiff"), width = 8, height = 8)
}

box_plot_maaslin <- function(output_dir,
                             physeq,
                             condition_var,
                             palette,
                             tax_level,
                             legend_lab = legend_name){
  results <- read.csv(file.path(output_dir, "significant_results.tsv"), sep = "\t")
  nfeatures <- length(unique(results$feature))
  if(nfeatures > 0){
    n = ceiling(nfeatures/5)
    ncol = min(5, nfeatures)
    
    if(is.null(palette)){
      palette <- scales::hue_pal()(nrow(unique(sample_data(physeq)[,condition_var])))
    }
    
    if(tax_level == "ASV"){
      physeq_plot <- prune_taxa(taxa_names(physeq) %in% results$feature, physeq)
      otu_table(physeq_plot) <- otu_table(physeq_plot)*100
      plot_abundance_ASV(physeq_plot, 
                         Color = condition_var, 
                         palette = palette,
                         n = n,
                         ncol = ncol,
                         legend_title = legend_lab) 
    }else{
      physeq_plot <- prune_taxa(tax_table(physeq)[, tax_level] %in% results$feature, physeq)
      otu_table(physeq_plot) <- otu_table(physeq_plot)*100
      plot_abundance(physeq_plot, 
                     Facet = tax_level, 
                     Color = condition_var, 
                     palette = palette,
                     n = n,
                     ncol = ncol,
                     legend_title = legend_lab) 
    }

    if(nfeatures == 1){w = 4}else if(nfeatures < 5){w = 3.2*nfeatures}else{w = 16}
    ggsave(file.path(output_dir, "significant_results.svg"), width = w, height = 6*n, dpi = 300, limitsize = FALSE)
    #ggsave(file.path(output_dir, "significant_results.tiff"), width = w, height = 6*n, dpi = 300)
  }
}

# min_abund = 0.0005
# min_prev = 0.02
# max_sig = 0.05
# tax_level = "Phylum"
# physeq = physeq_sub
# name = p$name
# fixed_effects = p$fixed
# random_effects = p$rand
# reference = p$ref

run_maaslin <- function(physeq, name, fixed_effects, 
                        random_effects = NULL, reference = NULL, 
                        min_abund = 0.0005, min_prev = 0.02,
                        max_sig = 0.05, tax_level = "ASV", legend_name = NULL) {
  
  # Prepare data
  metadata <- data.frame(sample_data(physeq))
  metadata[] <- lapply(metadata, as.factor)
  data <- data.frame(otu_table(physeq))
  stopifnot(all(colnames(data) == rownames(metadata)))
  
  # Create unique output name
  output_name <- paste0(name, "_", tax_level)
  
  # Run MaAsLin
  Maaslin2(
    input_data = data,
    input_metadata = metadata,
    output = file.path(maaslin.dir, output_name),
    min_abundance = min_abund,
    min_prevalence = min_prev,
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM",
    max_significance = max_sig,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    reference = reference,
    correction = "BH",
    standardize = FALSE
  )
  
  write.table(
    metadata, 
    file = file.path(maaslin.dir, output_name, "metadata.tsv"),  # Replace with your desired file path/name
    sep = "\t",             # Use tab delimiter for TSV
    row.names = FALSE,      # Exclude row names (usually unnecessary for metadata)
    col.names = TRUE,       # Include column headers
    quote = FALSE           # Disable quoting of strings
  )
  # Generate volcano plot
  volcano(file.path(maaslin.dir, output_name))
  
  # Generate box plot
  box_plot_maaslin(output_dir = file.path(maaslin.dir, output_name), physeq = physeq, condition_var = p$fixed,
                   palette = p$col, tax_level = tax_level, legend_lab = legend_name)
  
  # Process significant results if needed
  if (!is.null(max_sig)) {
    sig_path <- file.path(maaslin.dir, output_name, "significant_results.tsv")
    if (file.exists(sig_path)) {
      significant_results <- read.delim(sig_path)
      tax_df <- as.data.frame(tax_table(physeq))
      tax_df$feature <- rownames(tax_df)
      significant_results <- merge(significant_results, tax_df, 
                                   by = "feature", all.x = TRUE)
      openxlsx::write.xlsx(significant_results, 
                           file.path(maaslin.dir, output_name, "significant_results.xlsx"))
    }
  }
}

#### Functions - LEfSer ####

# Clean taxonomy names
clean_tax_names_lefser <- function(physeq) {
  tax <- as.data.frame(tax_table(physeq))
  tax[] <- lapply(tax, gsub, pattern = " ", replacement = "_")
  tax[] <- lapply(tax, gsub, pattern = "\\.", replacement = "_") #Very important for lefser
  tax$Species[grep("_spp.", tax$Species)] <- NA
  tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis","Family"] <- paste(
    tax[!is.na(tax$Family) & tax$Family == "Incertae_Sedis","Order"], "Incertae_Sedis", sep = "_")
  #tax[is.na(tax)] <- "unclassified"
  tax$Kingdom <- NULL
  tax_table(physeq) <- tax_table(as.matrix(tax))
  return(physeq)
}

run_lefser_analysis <- function(physeq, group_var, output_dir, 
                               # norm = "CPM", 
                               # taxa_rank = "all",
                               # kw_cutoff = 0.05,
                               # wilcoxon_cutoff = 0.05, 
                               # lda_cutoff = 2,
                               # bootstrap_n = 30,
                               col_palette = NULL) {
  
  # Check that group variable exists
  if (!group_var %in% colnames(sample_data(physeq))) {
    stop("Group variable not found in sample data")
  }
  
  # Check group sizes
  groups <- sample_data(physeq)[[group_var]]
  if (any(table(groups) < 2)) {
    sink(file = file.path(output_dir, "my_analysis_log.txt"))
    cat ("***************************************","\n")
    cat ("All groups must have at least 2 samples","\n")
    cat ("***************************************","\n","\n")
    sink()
    return(NULL)
  }
  
  # Redundant filter
  # physeq <- subset_samples(physeq, !is.na(sample_data(physeq)[[group_var]]))
  # physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  
  if (nsamples(physeq) == 0 || ntaxa(physeq) == 0) {
    message("No samples or taxa left for group: ", group_var)
    return(NULL)
  }
  
  # Extract otu table
  input_lfse <-  unclass(otu_table(physeq))
  
  # Extract coldata
  coldata_lfse <- as(sample_data(physeq), "data.frame")
  
  ## Create a SummarizedExperiment object
  SE_lfse <- SummarizedExperiment(assays = list(counts = input_lfse), colData = coldata_lfse, rowData = tax_table(physeq)[,c("Phylum", "Class", "Order", "Family", "Genus")])
  SE_lfse <- lefser::relativeAb(SE_lfse)
  
  marker <- tryCatch({
    my_lefserClades(SE_lfse, classCol = group_var)
  }, error = function(e) {
    message("LEfSe error: ", e)
    return(NULL)
  })
  
  if (is.null(nrow(marker))) {
    message("No significant markers for group: ", group_var)
    sink(file = file.path(output_dir, "my_analysis_log.txt"))
    cat ("***************************************","\n")
    cat ("No significant markers for group: ", group_var,"\n")
    cat ("***************************************","\n","\n")
    sink()
    return(NULL)
  }
  
  write.xlsx(marker, file.path(output_dir, paste("marker_table", nrow(marker), "m.xlsx", sep = "_")), rowNames = T)
  
  # Add color palette if col_palette is NULL
  if (is.null(col_palette)) {
    group_vals <- sample_data(physeq)[[group_var]]
    if (!is.factor(group_vals)) {
      group_vals <- factor(group_vals)
    }
    n_groups <- nlevels(group_vals)
    col_palette <- scales::hue_pal()(n_groups) # Default ggplot2 palette
  }
  
  create_plots <- function(res_sub, prefix, col_palette) {
    try({
      w = 16
      h = nrow(res_sub)*0.5
      
      p_val <- my_lefserPlot(res_sub, colors = col_palette, label.font.size = 6, other.font.size = 20)
      #ggsave(file.path(output_dir, paste0(prefix, "_lda_bar.tiff")), p_lda, width = w, height = h, dpi = 300)
      ggsave(file.path(output_dir, paste0(prefix, "_dot.svg")), p_val, width = w, height = h, dpi = 300)
      
      p_lda <- lefser::lefserPlot(res_sub, colors = unname(col_palette), label.font.size = 6)
      #ggsave(file.path(output_dir, paste0(prefix, "_lda_bar.tiff")), p_lda, width = w, height = h, dpi = 300)
      ggsave(file.path(output_dir, paste0(prefix, "_lda_bar.svg")), p_lda, width = w, height = h, dpi = 300)
      
     # p_clad <- lefser::lefserPlotClad(res_sub, showTipLabels = T, colors = col_palette)
      #ggsave(file.path(output_dir, paste0(prefix, "_clad_bar.svg")), p_clad, width = 10, height = 10, dpi = 300)
      
    })
  }
  
  # Top by LDA
  for (n in c(20, 50)) {
    marker$absolute <- abs(marker$scores)
    marker_lda <- marker[order(marker$absolute, decreasing = TRUE), ]
    if (nrow(marker_lda) > n) {
      marker_lda <- marker_lda[1:n, ]
      create_plots(marker_lda, paste0("top", n, "lda"), col_palette)
    }else{
      create_plots(marker_lda, paste0("all_markers"), col_palette)}
  }
  
  # Top by adjusted p-value
  for (n in c(20, 50)) {
    marker_padj <- marker[order(marker$kw_p.value), ]
    if (nrow(marker_padj) > n) {
      marker_padj <- marker_padj[1:n, ]
      create_plots(marker_padj, paste0("top", n, "pval"), col_palette)
    }else{create_plots(marker_padj, paste0("all_markers"), col_palette)}
  }
  
  return(marker)
}

#### Functions - Taxonomy plots ####

#Logarithmic scale for plots
marks_no_sci <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = F)

# Calculate taxa sums
taxasum = function(physeq_object, taxa){
  tapply(taxa_sums(physeq_object), tax_table(physeq_object)[, taxa], sum, na.rm=TRUE) %>%
    sort(. , TRUE)
}

# Stacked bar plot
facets_condition <- function(phy, var, top, taxa, nrow = 1){
  #sum = tapply(taxa_sums(phy), tax_table(phy)[, taxa], sum, na.rm=TRUE)
  #top = names(sort(sum, TRUE))[1:top]
  top = prune_taxa((tax_table(phy)[, taxa] %in% top), phy)
  top = tax_glom(top, taxa)
  plot_bar(top, taxa, fill=taxa) + geom_bar(stat = "Identity", position = "stack") +
    facet_wrap(var, nrow = nrow) +
    ylab("Relative abundance %") +
    theme(text = element_text(size=22),
          axis.text = element_text(size = 18),
          #axis.text.x.bottom = element_text(angle = -30),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 22),
          legend.title = element_text(size = 24)) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
     # breaks = c(0, 0.01, 0.1, 1, 10, 100),
      breaks = c(0, 1, 5, 10, 20, 30, 100),
      labels = function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)
    )
}

# Plöt abundance
plot_abundance = function(physeq, ylabn = "Relative abundance (%)",
                          Facet = "Species",
                          Color = "Diagnosis",
                          X_var = Color,
                          palette,
                          n = NULL,
                          ncol = NULL,
                          legend_title = NULL){
  mphyseq = tax_glom(physeq, Facet)
  mphyseq = psmelt(mphyseq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  marks_no_sci <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = F)
  P <- ggplot(data = mphyseq,
              mapping = aes_string(x = X_var, y = "Abundance",
                                   color = Color, fill = Color)) + thm +
    # geom_point(size = 1, alpha = 0.1,
    #            position = position_jitter(width = 0.3)) +
    geom_boxplot(fill = NA) +
    facet_wrap(facets = Facet, nrow = n, ncol = ncol) + ylab(ylabn) +
    # stat_compare_means(method = "wilcox") +
    scale_fill_manual(values = palette, drop = FALSE) + scale_color_manual(values = palette, drop = FALSE) +
    scale_x_discrete(drop = FALSE) + xlab(NULL) +
    #  scale_y_log10()
    scale_y_log10(labels = marks_no_sci)
  
  if (!is.null(legend_title)) {
    P <- P + labs(color = legend_title, fill = legend_title)
  } 
  print(P)
}

plot_abundance_ASV = function(physeq, ylabn = "Relative abundance (%)",
                              Color = "Diagnosis",
                              X_var = Color,
                              palette,
                              n = NULL,
                              ncol = NULL,
                              legend_title = NULL){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  marks_no_sci <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = F)
  P <- ggplot(data = mphyseq,
         mapping = aes_string(x = X_var, y = "Abundance",
                              color = Color, fill = Color)) + thm +
    # geom_point(size = 1, alpha = 0.1,
    #            position = position_jitter(width = 0.3)) +
    geom_boxplot(fill = NA) +
    facet_wrap(facets = "OTU", nrow = n, ncol = ncol) + ylab(ylabn) +
    # stat_compare_means(method = "wilcox") +
    scale_fill_manual(values = palette) + scale_color_manual(values = palette) +
    scale_x_discrete(drop = FALSE) + xlab(NULL) +
    #  scale_y_log10()
    scale_y_log10(labels = marks_no_sci)
  
  if (!is.null(legend_title)) {
    P <- P + labs(color = legend_title, fill = legend_title)
  } 
  print(P)
}

# Make phyloseq according to condition
phy_cond <- function(phy, condition, abundance = "relative"){
  
  # Transform all variables to factors just in case...
  df <- as.data.frame(lapply(sample_data(phy),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
  row.names(df) <- sample_names(phy)
  sample_data(phy) <- sample_data(df)
  
  # Merging samples according to condition
  ps.condition <- merge_samples(phy, condition)  
  sample_data(ps.condition)[,condition] <- rownames(sample_data(ps.condition))
  ps.condition_per <- transform_sample_counts(ps.condition, function(OTU) OTU/sum(OTU)*100) #2nd transformation to make it again in percentage
  # sample_data(ps.condition_per)[,condition] <- rownames(sample_data(ps.condition_per))
  
  if(abundance == "relative"){
    return(ps.condition_per)}else if(abundance == "absolute"){
      return(ps.condition)
    }
  
}

# Plots for questions
plot_most_abundant <- function(q, tax_levels = c("Phylum", "Genus"), n_top = 20, col_palette, legend_title = NULL) {
  # Create output directory for this question
  bar_output_dir <- file.path(bar.dir, params[[q]]$name)
  dir.create(bar_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Subset phyloseq objects for the question
  phy_abs <- subset_physeq(physeq_tax, q)
  phy_rel <- subset_physeq(physeq_tax_re, q)
  
  # Get grouping variable from parameters
  condition_var <- params[[q]]$fixed
  
  # Create merged phyloseq objects (absolute/relative)
  phy_abs_merged <- phy_cond(phy_abs, condition_var, "absolute")
  phy_rel_merged <- phy_cond(phy_rel, condition_var, "relative")
  
  # Generate palette for grouping variable
  group_levels <- levels(sample_data(phy_abs_merged)[[condition_var]])
  palette <- viridis::viridis(length(group_levels))
  names(palette) <- group_levels
  
  for (tax_level in tax_levels) {
    # Get top taxa for this level
    top_taxa <- names(taxasum(phy_rel, tax_level))[1:n_top]
    
    # Plot 1: Merged Absolute Abundance (Stacked Bar)
    p_abs <- plot_bar(prune_taxa(tax_table(phy_abs_merged)[, tax_level] %in% top_taxa, phy_abs_merged), 
                      fill = tax_level) +
      geom_bar(aes_string(color = tax_level, fill = tax_level), stat = "identity") +
      facet_wrap(as.formula(paste("~", condition_var)), scales = "free_x") +
      ylab("Absolute Abundance") + thm.r +
      guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1)) #+
    # scale_fill_manual(values = distinct_palette(n_top))
    
    # Plot 2: Merged Relative Abundance (Stacked Bar)
    p_rel <- plot_bar(prune_taxa(tax_table(phy_rel_merged)[, tax_level] %in% top_taxa, phy_rel_merged), 
                      fill = tax_level) +
      geom_bar(aes_string(color = tax_level, fill = tax_level), stat = "identity") +
      facet_wrap(as.formula(paste("~", condition_var)), scales = "free_x") +
      ylab("Relative Abundance (%)") + thm.r +
      guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1)) #+
    #  scale_fill_manual(values = distinct_palette(n_top))
    
    # Plot 3: Non-Merged Relative Abundance (Barplot)
    p_bar <- facets_condition(phy_rel_merged, condition_var, top_taxa, tax_level)
    
    # Plot 4: Non-Merged Relative Abundance (Boxplot)
    phy_top <- prune_taxa(tax_table(phy_rel)[, tax_level] %in% top_taxa, phy_rel)
    p_box <- plot_abundance(phy_top, 
                            Facet = tax_level, 
                            Color = condition_var, 
                            palette = col_palette,
                            n = 2,
                            legend_title = legend_title)  # Facet rows
    
    # Save plots
    for(format in c("svg")){
      ggsave(file.path(bar_output_dir, paste0(tax_level, "_absolute.", format)), p_abs, width = 12, height = 8, dpi = 300)
      ggsave(file.path(bar_output_dir, paste0(tax_level, "_relative.", format)), p_rel, width = 12, height = 8, dpi = 300)
      ggsave(file.path(bar_output_dir, paste0(tax_level, "_barplot.", format)), p_bar, width = 13, height = 8, dpi = 300)
      ggsave(file.path(bar_output_dir, paste0(tax_level, "_boxplot.", format)), p_box, width = 18, height = 16, dpi = 300)
    #  ggsave(file.path(bar_output_dir, paste0(tax_level, "_boxplot.", format)), p_box, width = 20, height = 9, dpi = 300)
    }
    
  }
}

#### Questions to answer ####

#' Questions
#' Q1: Is there a difference between Diseased vs the parallel Normal Lung?
#' Q2: Is there a difference between Diseased vs the parallel Normal Lung without samples with synchronous tumor?
#' Q3.1: What are the difference between NSCLC vs SCLC vs Bening Tumor?
#' Q3.2: What are the difference between NSCLC vs Bening Tumor in nonsmokers?
#' Q4.1: What are the difference between NSCLC main histologies (Adenocarcinoma vs Squamous)?
#' Q4.2: What are the difference between NSCLC main histologies (Adenocarcinoma vs Squamous) in nonsmokers?
#' Q4:.3 What are the difference between NSCLC main histologies (Adenocarcinoma vs Squamous) in smokers?
#' Q5: Does history of smoking have an impact on the lung microbiome?
#' Q6.1: Does the T stage have an impact on the lung microbiome?
#' Q6.2: Does the N stage have an impact on the lung microbiome?
#' Q6.3: Does the M stage have an impact on the lung microbiome?
#' Q7: Does the microbiome differ between patients who develop post operative pneumonia?

# Maaslin2 analysis parameters
params <- list(
  Q1.1a = list(name = "Q1.1a.all_directly_isolated_lung", fixed = "Lung", rand = "Study_Nr", ref = "Lung,Normal", col = lung_col, paired = T,  pairing_var = "Study_Nr"),
  Q1.2a = list(name = "Q1.2a.main-D_directly_isolated_lung", fixed = "Lung", rand = "Study_Nr", ref = "Lung,Normal", col = lung_col, paired = T,  pairing_var = "Study_Nr"),
  Q1.3 = list(name = "Q1.3a.main-D-individually_directly_isolated_lung"),
  Q2.1a = list(name = "Q2.1a.all_directly_isolated_non_synchronous_lung", fixed = "Lung", rand = "Study_Nr", ref = "Lung,Normal", col = lung_col, paired = T,  pairing_var = "Study_Nr"),
  Q2.2a = list(name = "Q2.2a.main-D_directly_isolated_non_synchronous_lung", fixed = "Lung", rand = "Study_Nr", ref = "Lung,Normal", col = lung_col, paired = T,  pairing_var = "Study_Nr"),
  Q2.3 = list(name = "Q2.3a.main-D-individually_directly_isolated_non_synchronous_lung"),
  Q3.1a = list(name = "Q3.1a.diagnosis", fixed = "Diagnosis", rand = NULL, ref = "Diagnosis,Benign", col = diagnosis_col),
  Q3.2a = list(name = "Q3.2a.diagnosis_nonsmokers", fixed = "Diagnosis", rand = NULL, ref = "Diagnosis,Benign", col = diagnosis_col),
  Q3.3a = list(name = "Q3.3a.diagnosis_smokers", fixed = "Diagnosis", rand = NULL, ref = "Diagnosis,Benign", col = diagnosis_col),
  Q4.1a = list(name = "Q4.1a.histology", fixed = "Histology.NSCLC", rand = NULL, ref = NULL, col = histology_col, lab_plot = "Histology"),
  Q4.2a = list(name = "Q4.2a.histology_nonsmokers", fixed = "Histology.NSCLC", rand = NULL, ref = NULL, col = histology_col, lab_plot = "Histology"),
  Q4.3a = list(name = "Q4.3a.histology_smokers", fixed = "Histology.NSCLC", rand = NULL, ref = NULL, col = histology_col, lab_plot = "Histology"),
  Q5.1a = list(name = "Q5.1a.main-D_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, lab_plot = "History of smoking"),
 Q5.2a = list(name = "Q5.2a.NSCLC_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, lab_plot = "History of smoking"),
 Q5.3a = list(name = "Q5.3a.Benign_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, lab_plot = "History of smoking"),
 Q5.4a = list(name = "Q5.4a.Adenocarcinoma_smoking", fixed = "History.of.smoking.y.n", rand = NULL, ref = "History.of.smoking.y.n,No", col = smoker_col, lab_plot = "History of smoking"),
 Q6.1a = list(name = "Q6.1a.NSCLC_T", fixed = "T", rand = NULL, ref = NULL, col = scales::hue_pal()(5)[c(4,2,5,1)], lab_plot = "T"),
 Q6.2a = list(name = "Q6.2a.NSCLC_N", fixed = "N", rand = NULL, ref = NULL, col = scales::hue_pal()(5)[c(3:4,2,5,1)]),
 Q6.2.1a = list(name = "Q6.2a.NSCLC_N_twolevels", fixed = "Ntwo", rand = NULL, ref = NULL, col = scales::hue_pal()(5)[c(3:4,2,5,1)], lab_plot = "N"),
 Q6.3a = list(name = "Q6.3a.NSCLC_M", fixed = "M", rand = NULL, ref = NULL, col = scales::hue_pal()(5)[3:4]),
 Q7a = list(name = "Q7a.NSCLC_PostopPneumonia", fixed = "Postop..pneumonia", rand = NULL, ref = NULL, col = scales::hue_pal()(2)[2:1], lab_plot = "Pneumonia"),
)

# Function to subset phyloseq object
subset_physeq <- function(physeq, q){
  # Remove samples with no reads
  physeq_0 <- subset_samples(physeq, colSums(otu_table(physeq)) != 0 & !is.na(colSums(otu_table(physeq))))
  
  # Pre-filter 
  physeq_0 <- physeq_0 %>%
    subset_samples(Sample_or_Control == "True sample" &
                     Isolation == "Direct_Isolation" &
                     LibrarySizeDecontam > 0)

  # Custom subsetting for specific questions
  if(grepl("Q1|Q2", q)){
    study_nrs <- sample_data(physeq_0)$Study_Nr
    dup_studies <- study_nrs[duplicated(study_nrs) | duplicated(study_nrs, fromLast = TRUE)]
    # FIX: Use logical vector with prune_samples()
    keep <- sample_data(physeq_0)$Study_Nr %in% dup_studies
    physeq_sub <- prune_samples(keep, physeq_0)  # This replaces subset_samples()
    #physeq_sub <- subset_samples(physeq_0, Study_Nr %in% dup_studies)
   
  
    if (grepl("Q1.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    } else if (grepl("Q2", q)) {
      physeq_sub <- subset_samples(physeq_sub, Synchronous.tumor == "No")
      if (grepl("Q2.2", q)) {
        physeq_sub <- subset_samples(physeq_sub, Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
      }}
  
  } else if (grepl("Q3", q)) {
    physeq_sub <- subset_samples(physeq_0,
                                 Lung == "Diseased" &
                                   Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    if (grepl("Q3.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    }else if (grepl("Q3.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    }
  } else if (grepl("Q4", q)) {
    physeq_sub <- subset_samples(physeq_0,
                                 Lung == "Diseased" & 
                                   Diagnosis == "NSCLC" &
                                   Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma")
    )
    if (grepl("Q4.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    }else if (grepl("Q4.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    }
  } else if (grepl("Q5", q)) {
    physeq_sub <- subset_samples(physeq_0,
                                 Lung == "Diseased" & 
                                   Diagnosis %in% c("Benign", "NSCLC", "SCLC") &
                                   !is.na(History.of.smoking.y.n)
    )
    if (grepl("Q5.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis == "NSCLC")
    }else if (grepl("Q5.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis == "Benign")
    }else if (grepl("Q5.4", q)) {
      physeq_sub <- subset_samples(physeq_sub, Diagnosis == "NSCLC" & Histology.NSCLC == "Adenocarcinoma")
    }
  } else if (grepl("Q6", q)) {
    physeq_sub <- subset_samples(physeq_0,
                                 Lung == "Diseased" & 
                                   Diagnosis == "NSCLC"
    )
    if (grepl("Q6.1", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(T))
    }else if (grepl("Q6.2", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(N))
    }else if (grepl("Q6.3", q)) {
      physeq_sub <- subset_samples(physeq_sub, !is.na(M))
    }
  } else if (grepl("Q7", q)) {
    physeq_sub <- subset_samples(physeq_0,
                                 Lung == "Diseased" & 
                                   Diagnosis == "NSCLC" &
                                   !is.na(Postop..pneumonia)
    )
  } else {
    physeq_sub <- physeq_0
  }
  
  physeq_sub <- prune_taxa(taxa_sums(physeq_sub) != 0, physeq_sub) 
  return(physeq_sub)
}

#### Beta diversity ####

# Select phyloseq object
physeq_beta0 <- physeq_re_bc

# Calculate and plot beta diversity

for(q in names(params)){
  
  print(q)
  
  # Subset samples 
  physeq_beta <- subset_physeq(physeq_beta0, q)
  
  # Calculate distances and ordinate
  beta_list <- beta_analysis(physeq_beta)
  
  # Plot variable
  beta_plot(physeq_beta, beta_list, q = q, variable = params[[q]]$fixed, col = params[[q]]$col, legend = params[[q]]$lab_plot)
  
  # Plot batch
  beta_plot(physeq_beta, beta_list, q = q, variable = "Batch", prefix = "_batch")
  
  # PERMANOVA
  ad <- list()
  for(i in names(beta_list[1:3])){
    set.seed(123)
    test.adonis1 <- adonis2(as.formula(paste("beta_list[[i]] ~", params[[q]]$fixed)), data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
    set.seed(123)
    test.adonis2 <- adonis2(as.formula(paste("beta_list[[i]] ~", params[[q]]$fixed, " + Batch")), data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
    
    ad[[i]] <- rbind(test.adonis1, test.adonis2)
    if(grepl("Q1|Q2", q)){
      set.seed(123)
      test.adonis3 <- adonis2(as.formula(paste("beta_list[[i]] ~", params[[q]]$fixed, " + Study_Nr")), data = data.frame(sample_data(physeq_beta)), permutations = 99999, by = "margin")
      
      ad[[i]] <- rbind(test.adonis1, test.adonis2, test.adonis3)
    }
  }
  ad <- as.data.frame(do.call(rbind, ad))
  ad
  write.xlsx(ad, file.path(b.stats, paste0(q, ".xlsx")), rowNames = T)
  
}

#### Abundant taxonomy plots ####

# Select phyloseq object
physeq_tax <- physeq
physeq_tax_re <- physeq_re_bc
otu_table(physeq_tax_re) <- otu_table(physeq_tax_re)*100
  
# Taxonomic levels to plot
tax_levels <- c("Phylum", "Family", "Genus")

# Plot most abundant taxa
for (q in names(params)) {
  p <- params[[q]]
  
  plot_most_abundant(q, tax_levels = tax_levels, col_palette = p$col, legend_title = p$lab_plot)  
  
}

#### Maaslin2 Analysis ####

# Select phyloseq object
physeq_maas0 <- physeq_re_bc

# Clean taxonomy names
physeq_maas0 <- clean_tax_names_maaslin(physeq_maas0)

# Taxonomic levels to analyze
tax_levels <- c("ASV", "Phylum", "Family", "Genus")

# Run Maaslin2 analyses
for (q in names(params)) {
  p <- params[[q]]
  
  # Custom subsetting for specific questions
  physeq_sub <- subset_physeq(physeq_maas0, q)

  # Run for each taxonomic level
  for (level in tax_levels) {
    if (level == "ASV") {
      # Use original phyloseq object
      run_maaslin(
        physeq = physeq_sub,
        name = p$name,
        fixed_effects = p$fixed,
        random_effects = p$rand,
        reference = p$ref,
        tax_level = level,
        legend_name = p$lab_plot
      )
    } else {
      # Agglomerate to specified taxonomic level
      phy_agg <- agglomerate_taxa(physeq_sub, level)
      run_maaslin(
        physeq = phy_agg,
        name = p$name,
        fixed_effects = p$fixed,
        random_effects = p$rand,
        reference = p$ref,
        tax_level = level,
        legend_name = p$lab_plot
      )
    }
  }
}

#### LEfSeR Analysis ####

# Select phyloseq object
physeq_lefse0 <- physeq_re_bc

# Clean taxonomy names
physeq_lefse0 <- clean_tax_names_lefser(physeq_lefse0)

# Run LEfSe analyses
for (q in names(params)) {
  p <- params[[q]]
  
  # Custom subsetting for specific questions
  physeq_processed <- subset_physeq(physeq_lefse0, q)
  
  # Run LEfSe only for ASV level and non-batch-adjusted questions
  lefse_output_dir <- file.path(lefser.dir, p$name)
  dir.create(lefse_output_dir, showWarnings = FALSE, recursive = TRUE)
  run_lefser_analysis(
    physeq = physeq_processed,
    group_var = p$fixed,
    output_dir = lefse_output_dir,
    col_palette = p$col
  )
}

#### Save matrix for correlations ####

# List of phyloseq objects to use for correlations
  physeq_list <- list(
    #abs = clean_tax_names_maaslin(physeq),
    rel = clean_tax_names_maaslin(physeq_re),
    rel_bc = clean_tax_names_maaslin(physeq_re_bc)
  )


# Loop through each phyloseq object and save matrices
for (name in names(physeq_list)) {
  physeq_current <- physeq_list[[name]]
  physeq_current <- subset_samples(physeq_current, colSums(otu_table(physeq_current)) != 0 & !is.na(colSums(otu_table(physeq_current))))

  physeq_current <- physeq_current %>%
    subset_samples(Sample_or_Control == "True sample" &
                     Isolation == "Direct_Isolation" &
                     LibrarySizeDecontam > 0)

  otu_table(physeq_current) <- otu_table(physeq_current)*100

  matrix <- list()  # Initialize results list
    # Process each taxonomic level
  for (taxon in c("Phylum", "Class", "Order", "Family", "Genus")) {
    physeq_glom <- tax_glom(physeq_current, taxrank = taxon)
    taxa_names(physeq_glom) <- tax_table(physeq_glom)[, taxon]
    matrix[[taxon]] <- merge(sample_data(physeq_glom), t(otu_table(physeq_glom)), by = 0)
  }

  # Add ASV-level data
  matrix[["ASV"]] <- merge(sample_data(physeq_current), t(otu_table(physeq_current)), by = 0)

  # Save to Excel
  write.xlsx(matrix, file = file.path(rhea.dir, paste0("correlation_matrix_", name, "_per.xlsx")))
}

#### End ####

#### Rhea ####
subset_data <- function(data, q){
  
  # Custom subsetting for specific questions
  if (grepl("Q1.2", q)) {
    data_sub <- dplyr::filter(data, Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
  } else if (grepl("Q2", q)) {
    data_sub <- dplyr::filter(data, Synchronous.tumor == "No")
    if (grepl("Q2.2", q)) {
      data_sub <- dplyr::filter(data_sub, Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    }
  } else if (grepl("Q3", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" &
                                Diagnosis %in% c("Benign", "NSCLC", "SCLC"))
    if (grepl("Q3.2", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    }else if (grepl("Q3.3", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    }
  } else if (grepl("Q4", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" & 
                                Diagnosis == "NSCLC" &
                                Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma")
    )
    if (grepl("Q4.2", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "No")
    }else if (grepl("Q4.3", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(History.of.smoking.y.n) & History.of.smoking.y.n == "Yes")
    }
  } else if (grepl("Q5", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" & 
                                Diagnosis %in% c("Benign", "NSCLC", "SCLC") &
                                !is.na(History.of.smoking.y.n)
    )
    if (grepl("Q5.2", q)) {
      data_sub <- dplyr::filter(data_sub, Diagnosis == "NSCLC")
    }else if (grepl("Q5.3", q)) {
      data_sub <- dplyr::filter(data_sub, Diagnosis == "Benign")
    }
  } else if (grepl("Q6", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" & 
                                Diagnosis == "NSCLC"
    )
    if (grepl("Q6.1", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(T))
    }else if (grepl("Q6.2", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(N))
    }else if (grepl("Q6.3", q)) {
      data_sub <- dplyr::filter(data_sub, !is.na(M))
    }
  } else if (grepl("Q7", q)) {
    data_sub <- dplyr::filter(data,
                              Lung == "Diseased" & 
                                Diagnosis == "NSCLC" &
                                !is.na(Postop..pneumonia)
    )
  } else {
    data_sub <- data
  }
  
  return(data_sub)
}

#' Set paramethers for Rhea

#'  Taxonomic levels to analyze
tax_levels <- c("ASV", "Phylum", "Family", "Genus")

#' Please give the input file name 
input_file <-file.path(rhea.dir, "correlation_matrix_rel_bc_per.xlsx")            #<--- CHANGE ACCORDINGLY !!!

#' Please enter the position in the table (column number) where the dependant variable starts (e.g. Richness)
#' Note: the first column containing sample names does not count!
dependant_variables_start <- 68

#' Please enter the order of the group names
#' If no group names are writting groups are ordered automatically
group_order=c("")

#' Please enter the position in the table (column number) where relative abundances of OTUs or taxonomic groups start
#' Note: the first column containing sample names does not count!
taxonomic_variables_start <- 68

#' The cutoff of relative abundance; all values below this cutoff will be zeroed (default cutoff is 0.5 %)
abundance_cutoff <- 0.05 #0.05% = 0.005

#' The prevalence cutoff; at least one group must have a number of samples above the selected treshold
#' for the variable to be tested (default cutoff is 0.3 = 30 % of samples are positive within a given group) 
prevalence_cutoff <- 0.02

#' The minimum median abundance value that must be observed in at least one group before statistical test is performed
max_median_cutoff <- 0

#' Whether the max_median_cutoff should be used
perform_max_median_cutoff <- F

#' Replace 0 Value with NA 
#' YES: Replace zeros with NA (Default)
#' NO: Consider zeros in statistics
ReplaceZero = "NO"

#' Set the graphical output parameter 
#' 1 = without individual values as dots
#' 2 = with individual values as dots
#' 3 = with individual values as dots and with sample names
PlotOption = 1

#' Set the significance cutoff level (default cutoff is 0.05 but it can be set lower)
sig.cutoff <- 0.05

for (q in names(params)) {
  p <- params[[q]]
  
#' The name of the independant variable that the analysis will be performed on
independant_variable_name <- params[[q]]$fixed

### For paired tests only:
if (!is.null(params[[q]]$paired)){
  if(isTRUE(params[[q]]$paired)){
    paired = T
    
    #' The name of the dependent variable that the analysis will be performed on (different time points)
    dependant_variable_name <- independant_variable_name
    
    #' The name of the ID variable that the analysis will be performed on (the id of the subject that was followed overtime)
    id_name = params[[q]]$pairing_var
    
    #' Missing values are removed by default ("NO").
    #' If missing values should be replaced by Skillings-Mack method please change this parameter to "YES"
    ReplaceMissingValues ="NO"
  }else(paired = F)
}else(paired = F)

### End of settings for paired tests

for (taxon in tax_levels) {
  
#' Taxon for analysis. Needs to be set if for loop is turned off 
taxon <- taxon

# Load the correlation matrix
original_table <- read.xlsx(input_file, sheet = taxon)
original_table <- subset_data(original_table, q)

# Run Rhea
source(file.path(res.dir, "Scripts/Rhea.R"))
  }
}
