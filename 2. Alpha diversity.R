
# Alpha diversity #

#### Package setup ####

Sys.setenv(language = "EN")
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(ggpattern)
library(grid)
library(gridExtra)
library(openxlsx)
library(ggstatsplot)
library(tidyr)
library(cowplot)
theme_set(theme_bw())

#### End ####

#### Load output from DADA2 pipeline with original phyloseq object ####

res.dir <- "Results_trimmed_cutadapt3_20241216"
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir,"4.physeq.decontam.RData"))

#### End ####

#### Create directories for results ####

a.stats <- file.path(res.dir, "2.Alpha_stats")
a.plots <- file.path(res.dir, "3.Alpha_plots")

dir.create(a.stats, recursive = T)
dir.create(a.plots, recursive = T)

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
#### End ####

#### Select samples ####

physeq_alpha = physeq

physeq_alpha = subset_samples(physeq_alpha, sample_names(physeq_alpha) %in% c(rownames(sample_info[!sample_info$Sample_or_Control %in% c("Negative control sample", "Positive control sample"),]))) 
physeq_alpha = subset_samples(physeq_alpha, sample_data(physeq_alpha)$Isolation == "Direct_Isolation") 
physeq_alpha = subset_samples(physeq_alpha, sample_data(physeq_alpha)$Diagnosis %in% c("Benign", "SCLC", "NSCLC")) 

sample_data(physeq_alpha)$LibrarySizeGroup <- ifelse(sample_data(physeq_alpha)$LibrarySizeDecontam >= 1000, ">=1000", "<1000")

#### End ####

#### Convert sample data to data frame ####

rich <- data.frame(sample_data(physeq_alpha))
dim(rich)

# Convert to long data frame for the facets
rich_long <- rich %>%
  pivot_longer(cols = c(Shannon, InvSimpson),
               names_to = "alpha_measure",
               values_to = "alpha_value")

rich_long$alpha_measure <- factor(rich_long$alpha_measure, levels = c("Shannon", "InvSimpson"))

#### End ####

#### Explore correlation of sequencing depth and alpha diversity ####

for(index in c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")){

ggscatterstats(
  data = rich,
  x = !!sym(index), 
  y = LibrarySizeDecontam, 
  xlab = index,
  ylab = "Library Size",
  point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
  marginal = F,
  type ="n",
  title = "Directly isolated samples",
) + theme(text = element_text(size = 32),
          plot.title = element_text(size = 28),
          plot.subtitle = element_text(size = 18))
  
ggsave(file.path(qc.dir, paste0("LibSize_", index, "_direct.tiff")), width = 11, height = 10, dpi = 300)
ggsave(file.path(qc.dir, paste0("LibSize_", index, "_direct.svg")), width = 11, height = 10, dpi = 300)

}

#### End ####

#### Q0: Is there a difference between Batches ####
kw <- list()

for (n in c("Shannon", "InvSimpson")) {
  kruskal <- kruskal.test(rich[,n] ~ rich[,"Batch"])
  eta_squared <- (kruskal$statistic - 3 + 1)/(nrow(rich) - 3)
  f <- sqrt(eta_squared/(1-eta_squared))
  
  tab <- c(kruskal$method, n, "Batch", kruskal$statistic, kruskal$parameter, kruskal$p.value, f)
  kw[[n]] <- tab
  
}

kw <- as.data.frame(do.call(rbind, kw))
colnames(kw) <- c("Test", "Variable1", "Variable2", "Chi2", "df", "p-value", "Effect size")

kw$p.adj.BH <- p.adjust(kw$`p-value`, method = "BH")
kw$p.adj.bonferroni <- p.adjust(kw$`p-value`, method = "bonferroni")
kw

write.xlsx(list(kw,rich), file.path(a.stats, "Q0.xlsx"), rowNames = T)

ph <- list()
for (n in c("Shannon", "InvSimpson")) {
  posthocw <- pairwise.wilcox.test(rich[,n],rich[,"Batch"], paired = F, p.adjust= "BH")
  df <- rbind(posthocw$p.value, c(posthocw$method, posthocw$p.adjust), c(n,"Batch"))
  ph[[n]] <- df
}
ph

write.xlsx(ph, file = file.path(a.stats, "Q0_posthoc.MW.xlsx"), rowNames = T)

P.Q0 <- ggplot(rich_long, aes(x = Batch, y = alpha_value, color = Batch, fill = Batch)) +
  geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "Batch", fill = "Batch") + thm

P.Q0

ggsave(plot = P.Q0, file.path(a.plots, "P.Q0.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q0, file.path(a.plots, "P.Q0.svg"), width = 10, height = 8, dpi = 300)

#### Q0.1: Is there a difference between <1000 and >= 1000 ASV counts? ####

wt <- list()

for (n in c("Shannon", "InvSimpson")) {
  W_test <- wilcox.test(rich[, n] ~ rich[, "LibrarySizeGroup"], data = rich, exact = T)
  z <- abs(qnorm(W_test$p.value/2))
  r <- z/sqrt(nrow(rich))
  
  tab <- c(W_test$method, n, "Lung", W_test$statistic, W_test$p.value, r)
  wt[[n]] <- tab
}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

wt$p.adj.BH <- p.adjust(wt$`p-value`, method = "BH")
wt$p.adj.bonferroni <- p.adjust(wt$`p-value`, method = "bonferroni")
wt

write.xlsx(list(wt,rich), file.path(a.stats, "Q0.1.xlsx"), rowNames = T)

P.Q0.1 <- ggplot(rich_long, aes(x = LibrarySizeGroup, y = alpha_value, color = LibrarySizeGroup, fill = LibrarySizeGroup)) +
  geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL) + thm

P.Q0.1

ggsave(plot = P.Q0.1, file.path(a.plots, "P.Q0.1.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q0.1, file.path(a.plots, "P.Q0.1.svg"), width = 10, height = 8, dpi = 300)

# In diseased lung
wt <- list()

for (n in c("Shannon", "InvSimpson")) {
  for (i in c("Diseased", "Normal")) {  
    W_test <- wilcox.test(rich[rich$Lung == i, n] ~ rich[rich$Lung == i, "LibrarySizeGroup"], data = rich, exact = T)
    z <- abs(qnorm(W_test$p.value/2))
    r <- z/sqrt(nrow(rich[rich$LibrarySizeGroup == i,]))
    
    tab <- c(W_test$method, n, "Lung", i, W_test$statistic, W_test$p.value, r)
    wt[[paste0(n, "_", i)]] <- tab
  }}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variable1", "Variable2", "Lung", "Statistic", "p-value", "Effect size")

wt$p.adj.BH.same.lung <- NA
wt$p.adj.bonferroni.same.lung <- NA

for (i in c("Diseased", "Normal")){
  wt[wt$Lung == i,]$p.adj.BH.same.lung <- p.adjust(wt[wt$Lung == i,]$`p-value`, method = "BH")
  wt[wt$Lung == i,]$p.adj.bonferroni.same.lung <- p.adjust(wt[wt$Lung == i,]$`p-value`, method = "bonferroni")
}
wt
write.xlsx(list(wt,rich), file.path(a.stats, "Q0.1b.xlsx"), rowNames = T)

Q0.1b <- ggplot(rich_long[rich_long$Lung == "Diseased",], aes(x = LibrarySizeGroup, y = alpha_value, color = LibrarySizeGroup, fill = LibrarySizeGroup)) +
  geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "Lung", color = "Diagnosis", fill = "Diagnosis") + thm.x

Q0.1b

ggsave(plot = Q0.1b, file.path(a.plots, "P.Q0.1b.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = Q0.1b, file.path(a.plots, "P.Q0.1b.svg"), width = 10, height = 8, dpi = 300)

#### Q1: Is there a difference between Diseased vs the parallel Normal Lung? ####

# In all diagnosis
wt <- list()

for (n in c("Shannon", "InvSimpson")) {
    W_test <- wilcox.test(rich[rich$Lung == "Diseased", n], rich[rich$Lung == "Normal", n], data = rich, paired = T, exact = T)
    z <- abs(qnorm(W_test$p.value/2))
    r <- z/sqrt(nrow(rich))
    
    tab <- c(W_test$method, n, "Lung", W_test$statistic, W_test$p.value, r)
    wt[[n]] <- tab
}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

wt$p.adj.BH <- p.adjust(wt$`p-value`, method = "BH")
wt$p.adj.bonferroni <- p.adjust(wt$`p-value`, method = "bonferroni")
wt

write.xlsx(list(wt,rich), file.path(a.stats, "Q1a.xlsx"), rowNames = T)

P.Q1a <- ggplot(rich_long, aes(x = Lung, y = alpha_value, color = Lung, fill = Lung)) +
  geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = lung_col) + scale_color_manual(values = lung_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL) + thm

P.Q1a

ggsave(plot = P.Q1a, file.path(a.plots, "P.Q1a.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q1a, file.path(a.plots, "P.Q1a.svg"), width = 10, height = 8, dpi = 300)

# In each of the three main diagnoses
wt <- list()

for (n in c("Shannon", "InvSimpson")) {
  for (i in c("NSCLC", "SCLC", "Benign")) {  
  W_test <- wilcox.test(rich[rich$Diagnosis == i & rich$Lung == "Diseased", n],
                          rich[rich$Diagnosis == i & rich$Lung == "Normal", n], data = rich, paired = T, exact = T)
  z <- abs(qnorm(W_test$p.value/2))
  r <- z/sqrt(nrow(rich[rich$Diagnosis == i,]))
  
  tab <- c(W_test$method, n, "Lung", i, W_test$statistic, W_test$p.value, r)
  wt[[paste0(n, "_", i)]] <- tab
}}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variable1", "Variable2", "Diagnosis", "Statistic", "p-value", "Effect size")

wt$p.adj.BH.same.diag <- NA
wt$p.adj.bonferroni.same.diag <- NA

for (i in c("NSCLC", "SCLC", "Benign")){
  wt[wt$Diagnosis == i,]$p.adj.BH.same.diag <- p.adjust(wt[wt$Diagnosis == i,]$`p-value`, method = "BH")
  wt[wt$Diagnosis == i,]$p.adj.bonferroni.same.diag <- p.adjust(wt[wt$Diagnosis == i,]$`p-value`, method = "bonferroni")
}
wt
write.xlsx(list(wt,rich), file.path(a.stats, "Q1b.xlsx"), rowNames = T)

P.Q1b <- ggplot(rich_long[rich_long$Diagnosis %in% c("NSCLC", "SCLC", "Benign"),], aes(x = Lung, y = alpha_value, color = Diagnosis, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = diagnosis_col) + scale_color_manual(values = diagnosis_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "Lung", color = "Diagnosis", fill = "Diagnosis") + thm.x

P.Q1b

ggsave(plot = P.Q1b, file.path(a.plots, "P.Q1b.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q1b, file.path(a.plots, "P.Q1b.svg"), width = 10, height = 8, dpi = 300)

#### Q2: Is there a difference between Diseased vs the parallel Normal Lung without samples with synchronous tumor (after decontamination, only directly isolated samples) ####

# In all diagnosis
wt <- list()

for (n in c("Shannon", "InvSimpson")) {
  W_test <- wilcox.test(rich[!is.na(rich$Synchronous.tumor) & rich$Synchronous.tumor == "No" & rich$Lung == "Diseased", n],
                        rich[!is.na(rich$Synchronous.tumor) & rich$Synchronous.tumor == "No" & rich$Lung == "Normal", n],
                        data = rich, paired = T, exact = T)
  z <- abs(qnorm(W_test$p.value/2))
  r <- z/sqrt(nrow(rich[rich$Synchronous.tumor == "No", ]))
  
  tab <- c(W_test$method, n, "Lung", W_test$statistic, W_test$p.value, r)
  wt[[n]] <- tab
}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

wt$p.adj.BH <- p.adjust(wt$`p-value`, method = "BH")
wt$p.adj.bonferroni <- p.adjust(wt$`p-value`, method = "bonferroni")
wt

write.xlsx(list(wt,rich), file.path(a.stats, "Q2a.xlsx"), rowNames = T)

P.Q2a <- ggplot(rich_long[!is.na(rich_long$Synchronous.tumor) & rich_long$Synchronous.tumor == "No",], aes(x = Lung, y = alpha_value, color = Lung, fill = Lung)) +
  geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = lung_col) + scale_color_manual(values = lung_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL) + thm

P.Q2a

ggsave(plot = P.Q2a, file.path(a.plots, "P.Q2a.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q2a, file.path(a.plots, "P.Q2a.svg"), width = 10, height = 8, dpi = 300)

# In each of the three main diagnoses
wt <- list()

for (n in c("Shannon", "InvSimpson")) {
  for (i in c("NSCLC", "SCLC")) {  
    W_test <- wilcox.test(rich[rich$Diagnosis == i & !is.na(rich$Synchronous.tumor) & 
                                 rich$Synchronous.tumor == "No" & rich$Lung == "Diseased", n],
                            rich[rich$Diagnosis == i & !is.na(rich$Synchronous.tumor) &
                                   rich$Synchronous.tumor == "No" & rich$Lung == "Normal", n], data = rich, paired = T, exact = T)
    z <- abs(qnorm(W_test$p.value/2))
    r <- z/sqrt(nrow(rich[rich$Diagnosis == i & rich$Synchronous.tumor == "No",]))
    
    tab <- c(W_test$method, n, "Lung", i, W_test$statistic, W_test$p.value, r)
    wt[[paste0(n, "_", i)]] <- tab
  }}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variable1", "Variable2", "Diagnosis", "Statistic", "p-value", "Effect size")

wt$p.adj.BH.same.diag <- NA
wt$p.adj.bonferroni.same.diag <- NA

for (i in c("NSCLC", "SCLC", "Benign")){
  wt[wt$Diagnosis == i,]$p.adj.BH.same.diag <- p.adjust(wt[wt$Diagnosis == i,]$`p-value`, method = "BH")
  wt[wt$Diagnosis == i,]$p.adj.bonferroni.same.diag <- p.adjust(wt[wt$Diagnosis == i,]$`p-value`, method = "bonferroni")
}

write.xlsx(list(wt,rich), file.path(a.stats, "Q2b.xlsx"), rowNames = T)

P.Q2b <- ggplot(rich_long[rich_long$Diagnosis %in% c("NSCLC", "SCLC") & rich_long$Synchronous.tumor == "No",], aes(x = Lung, y = alpha_value, color = Diagnosis, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = diagnosis_col) + scale_color_manual(values = diagnosis_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "Lung", color = "Diagnosis", fill = "Diagnosis") + thm.x

P.Q2b

ggsave(plot = P.Q2b, file.path(a.plots, "P.Q2b.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q2b, file.path(a.plots, "P.Q2b.svg"), width = 10, height = 8, dpi = 300)

#### Conclusion ####
#' There is no significant difference between diseased and normal lung, thus from now on
#' we will only consider the diseased lung

#### Q3: What are the difference between NSCLC vs SCLC vs Bening Tumor? ####
# Create data frames with main diagnoses and diseased lung
rich.diag <- rich[rich$Diagnosis %in% c("NSCLC", "SCLC", "Benign") & rich$Lung == "Diseased",]
rich_long.diag <- rich_long[rich_long$Diagnosis %in% c("NSCLC", "SCLC", "Benign") & rich_long$Lung == "Diseased",]

kw <- list()

for (n in c("Shannon", "InvSimpson")) {
    kruskal <- kruskal.test(rich.diag[,n] ~ rich.diag[,"Diagnosis"])
    eta_squared <- (kruskal$statistic - 3 + 1)/(nrow(rich.diag) - 3)
    f <- sqrt(eta_squared/(1-eta_squared))
    
    tab <- c(kruskal$method, n, "Diagnosis", kruskal$statistic, kruskal$parameter, kruskal$p.value, f)
    kw[[n]] <- tab
    
  }

kw <- as.data.frame(do.call(rbind, kw))
colnames(kw) <- c("Test", "Variable1", "Variable2", "Chi2", "df", "p-value", "Effect size")

kw$p.adj.BH <- p.adjust(kw$`p-value`, method = "BH")
kw$p.adj.bonferroni <- p.adjust(kw$`p-value`, method = "bonferroni")

write.xlsx(list(kw,rich.diag), file.path(a.stats, "Q3.xlsx"), rowNames = T)

ph <- list()
for (n in c("Shannon", "InvSimpson")) {
    posthocw <- pairwise.wilcox.test(rich.diag[,n],rich.diag[,"Diagnosis"], paired = F, p.adjust= "BH")
    df <- rbind(posthocw$p.value, c(posthocw$method, posthocw$p.adjust), c(n,"Diagnosis"))
    ph[[n]] <- df
}
ph

write.xlsx(ph, file = file.path(a.stats, "Q3_posthoc.MW.xlsx"), rowNames = T)

P.Q3 <- ggplot(rich_long.diag, aes(x = Diagnosis, y = alpha_value, color = Diagnosis, fill = Diagnosis)) +
  geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = diagnosis_col) + scale_color_manual(values = diagnosis_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "Diagnosis", fill = "Diagnosis") + thm

P.Q3

ggsave(plot = P.Q3, file.path(a.plots, "P.Q3.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q3, file.path(a.plots, "P.Q3.svg"), width = 10, height = 8, dpi = 300)

#### Q4: What are the difference between NSCLC main histologies (Adeno vs Squamous)? ####
# Create data frames with NSCLC, main histologies and diseased lung
rich.NSCLC <- rich[rich$Diagnosis == "NSCLC"  & rich$Lung == "Diseased" & 
                     rich$Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"),]
rich_long.NSCLC <- rich_long[rich_long$Diagnosis == "NSCLC" & rich_long$Lung == "Diseased" & 
                               rich_long$Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"),]

# Previous smokers and non smokers
mwut <- list()

for (n in c("Shannon", "InvSimpson")){
    if(nrow(rich.NSCLC) <51){exact = T}else{exact = F}
    U_test <- wilcox.test(rich.NSCLC[,n] ~ rich.NSCLC[,"Histology.NSCLC"], data = rich.NSCLC, exact = exact)
    z <- abs(qnorm(U_test$p.value/2))
    r <- z/sqrt(nrow(rich.NSCLC))           #Effect size calculation
    tab <- c(U_test$method, n, "Histology.NSCLC", U_test$statistic, U_test$p.value, r)
    mwut[[n]] <- tab
  }

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

mwut$p.adj.BH <- p.adjust(mwut$`p-value`, method = "BH")
mwut$p.adj.bonferroni <- p.adjust(mwut$`p-value`, method = "bonferroni")

mwut

write.xlsx(list(mwut,rich), file.path(a.stats, "Q4a.xlsx"), rowNames = T)

P.Q4a <- ggplot(rich_long.NSCLC, aes(x = Histology.NSCLC, y = alpha_value, color = Histology.NSCLC, fill = Histology.NSCLC)) +
  geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = histology_col) + scale_color_manual(values = histology_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "Histology", fill = "Histology") + thm +
  guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))

P.Q4a

ggsave(plot = P.Q4a, file.path(a.plots, "P.Q4a.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q4a, file.path(a.plots, "P.Q4a.svg"), width = 10, height = 8, dpi = 300)

# Only previous smokers
mwut <- list()

for (n in c("Shannon", "InvSimpson")){
  if(nrow(rich.NSCLC[rich.NSCLC$History.of.smoking.y.n == "Yes",]) <51){exact = T}else{exact = F}
  U_test <- wilcox.test(rich.NSCLC[rich.NSCLC$History.of.smoking.y.n == "Yes",n] ~
                          rich.NSCLC[rich.NSCLC$History.of.smoking.y.n == "Yes","Histology.NSCLC"], data = rich.NSCLC, exact = exact)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(rich.NSCLC[rich.NSCLC$History.of.smoking.y.n == "Yes",]))           #Effect size calculation
  tab <- c(U_test$method, n, "Histology.NSCLC", U_test$statistic, U_test$p.value, r)
  mwut[[n]] <- tab
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

mwut$p.adj.BH <- p.adjust(mwut$`p-value`, method = "BH")
mwut$p.adj.bonferroni <- p.adjust(mwut$`p-value`, method = "bonferroni")

mwut

write.xlsx(list(mwut,rich), file.path(a.stats, "Q4b.xlsx"), rowNames = T)

P.Q4b <- ggplot(rich_long.NSCLC[!is.na(rich_long.NSCLC$History.of.smoking.y.n) & rich_long.NSCLC$History.of.smoking.y.n == "Yes",], aes(x = Histology.NSCLC, y = alpha_value, color = Histology.NSCLC, fill = Histology.NSCLC)) +
  geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = histology_col) + scale_color_manual(values = histology_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "Histology", fill = "Histology") + thm +
  guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))

P.Q4b

ggsave(plot = P.Q4b, file.path(a.plots, "P.Q4b.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q4b, file.path(a.plots, "P.Q4b.svg"), width = 10, height = 8, dpi = 300)

#### Q5: Does history of smoking have an impact on the lung microbiome? ####

# NSCLC
mwut <- list()

for (n in c("Shannon", "InvSimpson")){
  if(nrow(rich.NSCLC) <51){exact = T}else{exact = F}
  U_test <- wilcox.test(rich.NSCLC[,n] ~ rich.NSCLC[,"History.of.smoking.y.n"], data = rich.NSCLC, exact = exact)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(rich.NSCLC))           #Effect size calculation
  tab <- c(U_test$method, n, "History.of.smoking.y.n", U_test$statistic, U_test$p.value, r)
  mwut[[n]] <- tab
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

mwut$p.adj.BH <- p.adjust(mwut$`p-value`, method = "BH")
mwut$p.adj.bonferroni <- p.adjust(mwut$`p-value`, method = "bonferroni")

mwut

write.xlsx(list(mwut,rich), file.path(a.stats, "Q5a.xlsx"), rowNames = T)

P.Q5a <- ggplot(subset(rich_long.NSCLC, !is.na(History.of.smoking.y.n)), aes(x = History.of.smoking.y.n, y = alpha_value, color = History.of.smoking.y.n, fill = History.of.smoking.y.n)) +
  geom_point() + geom_boxplot(alpha = 0.5, na.rm = T) +
  scale_fill_manual(values = smoker_col) + scale_color_manual(values = smoker_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "History of smoking", fill = "History of smoking") + thm 

P.Q5a

ggsave(plot = P.Q5a, file.path(a.plots, "P.Q5a.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q5a, file.path(a.plots, "P.Q5a.svg"), width = 10, height = 8, dpi = 300)

# NSCLC Adenocarcinoma
# Create data frames with Adenocarcinoma and diseased lung
rich.adeno <- rich.NSCLC[rich.NSCLC$Histology.NSCLC == "Adenocarcinoma",]
rich_long.adeno <- rich_long.NSCLC[rich_long.NSCLC$Histology.NSCLC == "Adenocarcinoma",]

mwut <- list()

for (n in c("Shannon", "InvSimpson")){
  if(nrow(rich.NSCLC) <51){exact = T}else{exact = F}
  U_test <- wilcox.test(rich.adeno[,n] ~ rich.adeno[,"History.of.smoking.y.n"], data = rich.adeno, exact = exact)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(rich.adeno))           #Effect size calculation
  tab <- c(U_test$method, n, "History.of.smoking.y.n", U_test$statistic, U_test$p.value, r)
  mwut[[n]] <- tab
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

mwut$p.adj.BH <- p.adjust(mwut$`p-value`, method = "BH")
mwut$p.adj.bonferroni <- p.adjust(mwut$`p-value`, method = "bonferroni")

mwut

write.xlsx(list(mwut,rich), file.path(a.stats, "Q5aa.xlsx"), rowNames = T)

P.Q5aa <- ggplot(subset(rich_long.adeno, !is.na(History.of.smoking.y.n)), aes(x = History.of.smoking.y.n, y = alpha_value, color = History.of.smoking.y.n, fill = History.of.smoking.y.n)) +
  geom_point() + geom_boxplot(alpha = 0.5, na.rm = T) +
  scale_fill_manual(values = smoker_col) + scale_color_manual(values = smoker_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "History of smoking", fill = "History of smoking") + thm 

P.Q5aa

ggsave(plot = P.Q5aa, file.path(a.plots, "P.Q5aa.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q5aa, file.path(a.plots, "P.Q5aa.svg"), width = 10, height = 8, dpi = 300)

# Benign
# Create data frame with benign and diseased lung
rich.ben <- rich[rich$Diagnosis == "Benign"  & rich$Lung == "Diseased",]
rich_long.ben <- rich_long[rich_long$Diagnosis == "Benign" & rich_long$Lung == "Diseased",]

mwut <- list()

for (n in c("Shannon", "InvSimpson")){
  if(nrow(rich.ben) <51){exact = T}else{exact = F}
  U_test <- wilcox.test(rich.ben[,n] ~ rich.ben[,"History.of.smoking.y.n"], data = rich.ben, exact = exact)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(rich.ben))           #Effect size calculation
  tab <- c(U_test$method, n, "History.of.smoking.y.n", U_test$statistic, U_test$p.value, r)
  mwut[[n]] <- tab
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

mwut$p.adj.BH <- p.adjust(mwut$`p-value`, method = "BH")
mwut$p.adj.bonferroni <- p.adjust(mwut$`p-value`, method = "bonferroni")

mwut

write.xlsx(list(mwut,rich), file.path(a.stats, "Q5b.xlsx"), rowNames = T)

P.Q5b <- ggplot(subset(rich_long.ben, !is.na(History.of.smoking.y.n)), aes(x = History.of.smoking.y.n, y = alpha_value, color = History.of.smoking.y.n, fill = History.of.smoking.y.n)) +
  geom_point() + geom_boxplot(alpha = 0.5, na.rm = T) +
  scale_fill_manual(values = smoker_col) + scale_color_manual(values = smoker_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "History of smoking", fill = "History of smoking") + thm 

P.Q5b

ggsave(plot = P.Q5b, file.path(a.plots, "P.Q5b.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q5b, file.path(a.plots, "P.Q5b.svg"), width = 10, height = 8, dpi = 300)

# Correlation between packyears and alpha diversity
pack.years.NSCLC <- list()
for(index in c("Shannon", "InvSimpson")){
  
  pack.years.NSCLC[[index]] <- ggscatterstats(
    data = rich.NSCLC,
    y = !!sym(index), 
    x = Pack.years, 
    ylab = index,
    xlab = "Pack years",
    point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
    marginal = F,
    type ="n",
    title = NULL
  ) + theme(text = element_text(size = 32),
            plot.subtitle = element_text(size = 16))
}

plot_grid(top = textGrob("Relationship between alpha diversity and pack years in NSCLC",gp=gpar(fontsize=30, fontface="bold")),
          nrow= 2, rel_heights = c(1,8),
          plot_grid(pack.years.NSCLC[[1]], pack.years.NSCLC[[2]]))

ggsave(file.path(a.plots, "Q5c_NSCLC.tiff"), width = 18, height = 10, dpi = 300)
ggsave(file.path(a.plots, "Q5c_NSCLC.svg"), width = 18, height = 10, dpi = 300)

pack.years.ben <- list()
for(index in c("Shannon", "InvSimpson")){
  
  pack.years.ben[[index]] <- ggscatterstats(
    data = rich.ben,
    y = !!sym(index), 
    x = Pack.years, 
    ylab = index,
    xlab = "Pack years",
    point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
    marginal = F,
    type ="n",
    title = NULL
  ) + theme(text = element_text(size = 32),
            plot.subtitle = element_text(size = 16))
}

plot_grid(top = textGrob("Relationship between alpha diversity and pack years in benign tumor",gp=gpar(fontsize=30, fontface="bold")),
          nrow= 2, rel_heights = c(1,8),
          plot_grid(pack.years.ben[[1]], pack.years.ben[[2]]))

ggsave(file.path(a.plots, "Q5d_benign.tiff"), width = 18, height = 10, dpi = 300)
ggsave(file.path(a.plots, "Q5d_benign.svg"), width = 18, height = 10, dpi = 300)


#### Q6: Does the N and M stage of NSCLC have an impact on the lung microbiome ####

# T
kw <- list()

for (n in c("Shannon", "InvSimpson")) {
  kruskal <- kruskal.test(rich.NSCLC[,n] ~ rich.NSCLC[,"T"])
  eta_squared <- (kruskal$statistic - 3 + 1)/(nrow(rich.NSCLC) - 3)
  f <- sqrt(eta_squared/(1-eta_squared))
  
  tab <- c(kruskal$method, n, "T", kruskal$statistic, kruskal$parameter, kruskal$p.value, f)
  kw[[n]] <- tab
  
}

kw <- as.data.frame(do.call(rbind, kw))
colnames(kw) <- c("Test", "Variable1", "Variable2", "Chi2", "df", "p-value", "Effect size")

kw$p.adj.BH <- p.adjust(kw$`p-value`, method = "BH")
kw$p.adj.bonferroni <- p.adjust(kw$`p-value`, method = "bonferroni")
kw

write.xlsx(list(kw,rich), file.path(a.stats, "Q6a.xlsx"), rowNames = T)

ph <- list()
for (n in c("Shannon", "InvSimpson")) {
  posthocw <- pairwise.wilcox.test(rich.NSCLC[,n],rich.NSCLC[,"T"], paired = F, p.adjust= "BH")
  df <- rbind(posthocw$p.value, c(posthocw$method, posthocw$p.adjust), c(n,"T"))
  ph[[n]] <- df
}
ph

write.xlsx(ph, file = file.path(a.stats, "Q6a_posthoc.MW.xlsx"), rowNames = T)

P.Q6a <- ggplot(subset(rich_long.NSCLC, !is.na(T)), aes(x = T, y = alpha_value)) + 
  geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "T") + theme(text = element_text(size = text_size))

P.Q6a

ggsave(plot = P.Q6a, file.path(a.plots, "P.Q6a.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q6a, file.path(a.plots, "P.Q6a.svg"), width = 10, height = 8, dpi = 300)

# N 
kw <- list()

for (n in c("Shannon", "InvSimpson")) {
  kruskal <- kruskal.test(rich.NSCLC[,n] ~ rich.NSCLC[,"N"])
  eta_squared <- (kruskal$statistic - 3 + 1)/(nrow(rich.NSCLC) - 3)
  f <- sqrt(eta_squared/(1-eta_squared))
  
  tab <- c(kruskal$method, n, "N", kruskal$statistic, kruskal$parameter, kruskal$p.value, f)
  kw[[n]] <- tab
  
}

kw <- as.data.frame(do.call(rbind, kw))
colnames(kw) <- c("Test", "Variable1", "Variable2", "Chi2", "df", "p-value", "Effect size")

kw$p.adj.BH <- p.adjust(kw$`p-value`, method = "BH")
kw$p.adj.bonferroni <- p.adjust(kw$`p-value`, method = "bonferroni")
kw

write.xlsx(list(kw,rich), file.path(a.stats, "Q6b.xlsx"), rowNames = T)

ph <- list()
for (n in c("Shannon", "InvSimpson")) {
  posthocw <- pairwise.wilcox.test(rich.NSCLC[,n],rich.NSCLC[,"N"], paired = F, p.adjust= "BH")
  df <- rbind(posthocw$p.value, c(posthocw$method, posthocw$p.adjust), c(n,"N"))
  ph[[n]] <- df
}
ph

write.xlsx(ph, file = file.path(a.stats, "Q6b_posthoc.MW.xlsx"), rowNames = T)

P.Q6b <- ggplot(subset(rich_long.NSCLC, !is.na(N)), aes(x = N, y = alpha_value)) + 
  geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "N") + theme(text = element_text(size = text_size))

P.Q6b

ggsave(plot = P.Q6b, file.path(a.plots, "P.Q6b.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q6b, file.path(a.plots, "P.Q6b.svg"), width = 10, height = 8, dpi = 300)

# M 

mwut <- list()

for (n in c("Shannon", "InvSimpson")){
  if(nrow(rich.NSCLC) <51){exact = T}else{exact = F}
  U_test <- wilcox.test(rich.NSCLC[,n] ~ rich.NSCLC[,"M"], data = rich.NSCLC, exact = exact)
  z <- abs(qnorm(U_test$p.value/2))
  r <- z/sqrt(nrow(rich.NSCLC))           #Effect size calculation
  tab <- c(U_test$method, n, "M", U_test$statistic, U_test$p.value, r)
  mwut[[n]] <- tab
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

mwut$p.adj.BH <- p.adjust(mwut$`p-value`, method = "BH")
mwut$p.adj.bonferroni <- p.adjust(mwut$`p-value`, method = "bonferroni")

mwut

write.xlsx(list(mwut,rich), file.path(a.stats, "Q6c.xlsx"), rowNames = T)

P.Q6c <- ggplot(subset(rich_long.NSCLC, !is.na(M)), aes(x = M, y = alpha_value)) + 
  geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "M") + theme(text = element_text(size = text_size))

P.Q6c

ggsave(plot = P.Q6c, file.path(a.plots, "P.Q6c.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P.Q6c, file.path(a.plots, "P.Q6c.svg"), width = 10, height = 8, dpi = 300)

#### End ####
