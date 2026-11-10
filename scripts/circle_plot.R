# Circular plot of the NSCLC cohort #
# Supplementary Figure 7: patient characteristics (history of smoking,
# histology, T, N and M stage) and the relative abundance of the ten most
# abundant phyla per sample (tumour-bearing lung, directly isolated), with the
# samples clustered by Bray-Curtis distance at phylum level (hierarchical
# clustering, complete linkage).

#### Package setup ####

suppressPackageStartupMessages({
  library(phyloseq)
  library(circlize)
})

#### End ####

#### Load batch-corrected phyloseq object (2.MMUPHin.R) ####

res.dir <- "Results_direct"

rdata <- new.env()
load(file.path(res.dir, "RData", "6.phyloseq.corrected.RData"), envir = rdata)
physeq_circ <- rdata$physeq_re_bc
histology_col <- rdata$histology_col
smoker_col <- rdata$smoker_col
rm(rdata)

circ.dir <- file.path(res.dir, "circle_plot")
dir.create(circ.dir, recursive = T, showWarnings = F)

#### End ####

#### Samples, phyla and clustering ####

sd <- data.frame(sample_data(physeq_circ))
physeq_circ <- prune_samples(sd$Sample_or_Control == "True sample" & sd$Isolation == "Direct_Isolation" &
                               sd$Lung == "Diseased" & sd$Diagnosis == "NSCLC", physeq_circ)
cat("Samples:", nsamples(physeq_circ), "\n")

# Phylum-level relative abundances
p.phylum <- tax_glom(physeq_circ, taxrank = "Phylum")

# Ten most abundant phyla (summed over all samples)
phylum_sums <- tapply(taxa_sums(physeq_circ), tax_table(physeq_circ)[, "Phylum"], sum, na.rm = TRUE)
top <- names(sort(phylum_sums, TRUE))[1:10]
top.ps <- prune_taxa(tax_table(p.phylum)[, "Phylum"] %in% top, p.phylum)

data.p <- as(otu_table(top.ps), "matrix")
if (!taxa_are_rows(top.ps)) data.p <- t(data.p)
rownames(data.p) <- tax_table(top.ps)[rownames(data.p), "Phylum"]
data.p <- as.data.frame(t(data.p))  # samples x phyla

# Order of the phyla (stacking order and colours) as in the published figure;
# by abundance if the ten phyla were different
phylum_order <- c("Pseudomonadota", "Patescibacteria", "Campylobacterota", "Actinomycetota", "Bdellovibrionota",
                  "Bacillota", "Chloroflexota", "Cyanobacteriota", "Bacteroidota", "Fusobacteriota")
if (!setequal(phylum_order, colnames(data.p))) phylum_order <- top
data.p <- data.p[, phylum_order]

# Cluster the samples by Bray-Curtis distance (all phyla)
hc <- hclust(phyloseq::distance(p.phylum, method = "bray"))
labels <- hc$labels[hc$order]
dend <- as.dendrogram(hc)  # drawn in black by circos.dendrogram()

# Samples in the order of the clustering
data.p <- data.p[labels, ]
sd <- data.frame(sample_data(p.phylum))[labels, ]

#### End ####

#### Colours ####

barcolor <- scales::hue_pal()(ncol(data.p))  # phyla

smoker_col["Yes"] <- "#990000"
histology_col["Others"] <- "#FFD700"
sd$Histology.NSCLC[sd$Histology.NSCLC %in% "Adeno-Squamous cell carcinoma"] <- "Others"

# Colour scale from light to dark over the stages
stage_colors <- function(stages, low_color, high_color) {
  stages <- sort(as.numeric(stages))
  normalized <- (stages - min(stages)) / (max(stages) - min(stages))
  rgb_colors <- colorRamp(c(low_color, high_color))(normalized)
  setNames(apply(rgb_colors, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255)), stages)
}
stages_present <- function(x) as.character(sort(unique(as.numeric(as.character(x[!is.na(x)])))))
colors_T <- stage_colors(stages_present(sd$T), "#E6F0FF", "#0066CC")  # blue
colors_N <- stage_colors(stages_present(sd$N), "#FFEFDB", "#FF7700")  # orange
colors_M <- stage_colors(stages_present(sd$M), "#E5FFE5", "#008000")  # green

#### End ####

#### Circle plot ####

# One ring with one coloured cell per sample (white if the value is missing)
category_track <- function(values, colors) {
  cols <- colors[as.character(values)]
  cols[is.na(cols)] <- "white"
  circos.track(ylim = c(0, 1), track.margin = c(0, 0), track.height = .05, bg.lty = 0,
               panel.fun = function(x, y) {
                 pos <- seq_along(values) - .5
                 circos.rect(pos - .5, 0, pos + .5, 1, border = NA, col = cols)
               })
}

draw_circle_plot <- function(file, sample_labels = FALSE) {
  svg(file, width = 12, height = 10)
  circos.clear()
  circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 2, points.overflow.warning = F)
  circos.initialize("a", xlim = c(0, nrow(data.p)))

  # Sample names (outermost ring, optional)
  if (sample_labels) {
    circos.track(ylim = c(0, 1), bg.border = NA, track.height = .01, track.margin = c(.01, 0),
                 panel.fun = function(x, y) {
                   for (i in seq_along(labels)) {
                     circos.text(i - .5, 0, labels[i], adj = c(0, .5), facing = "clockwise", niceFacing = T, cex = .8)
                   }
                 })
  }

  # Patient characteristics
  category_track(sd$History.of.smoking.y.n, smoker_col)
  category_track(sd$Histology.NSCLC, histology_col)
  category_track(sd$T, colors_T)
  category_track(sd$N, colors_N)
  category_track(sd$M, colors_M)

  # Stacked relative abundance of the ten phyla
  circos.track(ylim = c(0, 1), track.margin = c(0, 0), track.height = .35, bg.lty = 0,
               panel.fun = function(x, y) {
                 mat <- as.matrix(data.p)
                 pos <- seq_len(nrow(mat)) - .5
                 for (i in seq_len(ncol(mat))) {
                   seq1 <- rowSums(mat[, seq_len(i - 1), drop = F])
                   seq2 <- rowSums(mat[, seq_len(i), drop = F])
                   circos.rect(pos - .5, seq1, pos + .5, seq2, border = NA, lwd = .1, col = barcolor[i])
                 }
               })

  # Dendrogram of the Bray-Curtis clustering
  circos.track(ylim = c(0, attr(dend, "height")), bg.border = NA, track.margin = c(0, .0015), track.height = .35,
               panel.fun = function(x, y) circos.dendrogram(dend))

  # Legends
  legend("topleft", legend = c("Yes", "No"), fill = smoker_col[c("Yes", "No")],
         title = "History of smoking", cex = 1.3, title.cex = 1.4, text.font = 2)
  hist_levels <- intersect(names(histology_col), unique(sd$Histology.NSCLC))
  legend("bottomleft", legend = hist_levels, fill = histology_col[hist_levels],
         title = "Histology", cex = 1.3, title.cex = 1.4, text.font = 2)
  legend("bottomright", legend = colnames(data.p), fill = barcolor,
         title = "Phylum", cex = 1.3, title.cex = 1.4, text.font = 2)
  legend("topright", legend = c(names(colors_T), names(colors_N), names(colors_M)),
         fill = c(colors_T, colors_N, colors_M), title = "T/N/M", ncol = 3,
         cex = 1.3, title.cex = 1.4, text.font = 2)

  circos.clear()
  dev.off()
}

draw_circle_plot(file.path(circ.dir, "circle_plot_bray.svg"))
draw_circle_plot(file.path(circ.dir, "circle_plot_bray_labels.svg"), sample_labels = TRUE)

# Source data: phylum relative abundances and patient characteristics, in plot order
write.csv(cbind(sd[, c("Study_Nr", "History.of.smoking.y.n", "Histology.NSCLC", "T", "N", "M")], data.p),
          file.path(circ.dir, "circle_plot_data.csv"))

#### End ####
