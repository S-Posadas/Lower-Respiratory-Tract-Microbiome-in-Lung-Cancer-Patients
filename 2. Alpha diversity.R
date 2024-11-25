
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
theme_set(theme_bw())

#### End ####

#### Load output from DADA2 pipeline with original phyloseq object ####

res.dir <- "Results_trimmed_cutadapt3_20241118"
R.dir <- file.path(res.dir, "RData")
load(file.path(R.dir,"4.physeq.decontam.RData"))
physeq_o <- physeq
colnames(physeq_o)
#### End ####



#### Automatic statistics and plots for exploration of alpha diversity ####
#### Define the variables and their possible values for a loop ####
only_tumor <- c(TRUE, FALSE)
only_non_synchronous <- c(TRUE, FALSE)
only_direct <- c(TRUE, FALSE)

# Generate all combinations
combinations <- expand.grid(
  #only_direct = only_direct,
  only_tumor = only_tumor,
  only_non_synchronous = only_non_synchronous
  
)

for (i in 1:nrow(combinations)) {
  only_tumor <- combinations$only_tumor[i]
  only_non_synchronous <- combinations$only_non_synchronous[i]
 # only_direct <- combinations$only_direct[i]
  
#### Select samples of interest ####

only_tumor = F
only_non_synchronous = F
all_Diagnosis = T
selected_Diagnosis = c("NSCLC", "SCLC", "Benign")
only_main_histology = F
only_smoker = F
only_non_smoker = F
only_direct = T
only_culture = F

prefix = "2024.10_"

physeq = physeq_o
physeq = subset_samples(physeq, sample_names(physeq) %in% c(rownames(sample_info[!sample_info$Sample_or_Control %in% c("Negative control sample", "Positive control sample"),]))) 


if (only_direct == T){
  physeq = subset_samples(physeq, sample_data(physeq)$Isolation == "Direct_Isolation") 
  prefix = paste0(prefix, "direct_isolation_")
}else if (only_culture == T){
  physeq = subset_samples(physeq, sample_data(physeq)$Isolation == "Culture_Enriched")
  prefix = paste0(prefix, "culture_enriched_")
}

if (only_tumor == T){
  physeq = subset_samples(physeq, sample_data(physeq)$Lung == "Diseased")
  prefix = paste0(prefix, "diseased_lung_")
}

if (only_non_synchronous == T & (all_Diagnosis == T | selected_Diagnosis[1] != "Benign")){
  physeq = subset_samples(physeq, sample_data(physeq)$Synchronous.tumor == "No")
  prefix = paste0(prefix, "no_synchronous_")
}

if (only_main_histology == T & all_Diagnosis == F & length(selected_Diagnosis) == 1){if (selected_Diagnosis == "NSCLC"){
  physeq = subset_samples(physeq, sample_data(physeq)$Histology.NSCLC %in% c("Adenocarcinoma", "Squamous cell carcinoma"))
  prefix = paste0(prefix, "main_histology_")
}}

if (only_smoker == T){
  physeq = subset_samples(physeq, sample_data(physeq)$History.of.smoking %in% c("Cigarettes", "Pipe"))
  prefix = paste0(prefix, "smokers_")
}else if (only_non_smoker == T & all_Diagnosis == F & length(selected_Diagnosis) == 1){if (selected_Diagnosis == "NSCLC"){
  physeq = subset_samples(physeq, sample_data(physeq)$History.of.smoking %in% c("None"))
  prefix = paste0(prefix, "nonsmokers_") }}

if (all_Diagnosis == T){
  physeq = physeq
  prefix = paste0(prefix, "all_diag_")
} else {
  physeq = subset_samples(physeq, sample_data(physeq)$Diagnosis %in% selected_Diagnosis)
  prefix = paste0(prefix, paste(selected_Diagnosis, collapse = "."))
  
}

#### End ####

#### Create directories for results ####

a.stats <- file.path(res.dir, "2.Alpha_stats_2024.10", prefix)
a.plots <- file.path(res.dir, "3.Alpha_plots_2024.10", prefix)

dir.create(a.stats, recursive = T)
dir.create(a.plots, recursive = T)

#### End ####

#### Convert sample data to data frame ####

rich <- data.frame(sample_data(physeq))
dim(rich)

# Convert to long data frame for the facets
rich_long <- rich %>%
  pivot_longer(cols = c(Observed, Shannon, InvSimpson),
               names_to = "alpha_measure",
               values_to = "alpha_value")

rich_long$alpha_measure <- factor(rich_long$alpha_measure, levels = c("Observed", "Shannon", "InvSimpson"))

#### End ####

#### Statistics ####
# Lung and Isolation are paired -> Wilcoxon
# Rest is unpaired -> Mann-Whitney or Kruskal-Wallis

# Mann-Whitney-U-Test 
if(length(selected_Diagnosis) == 1 & all_Diagnosis == F){if (selected_Diagnosis == "SCLC"){
  test.mw = c("Sex", "Active.smoker", "Side")}
  else if(selected_Diagnosis == "Benign"){ test.mw = c("Sex", "Side", "History.of.smoking.y.n")}
  else if(selected_Diagnosis == "NSCLC"){
    if(only_main_histology == T & only_non_smoker == F & only_smoker == F){
      test.mw = c("Sex", "Active.smoker", "Side", "Histology.NSCLC", "History.of.smoking.y.n")
    }else if(only_main_histology == T & (only_non_smoker == T | only_smoker == T)){
      test.mw = c("Sex", "Active.smoker", "Side", "Histology.NSCLC")
    }else {test.mw = c("Sex", "Active.smoker", "Side", "History.of.smoking.y.n")}}
}else{ test.mw = c("Sex", "Active.smoker", "Side", "History.of.smoking.y.n")}

mwut <- list()

for (n in c("Shannon", "InvSimpson")){
  for (i in test.mw){
    if(nrow(rich) <51){exact = T}else{exact = F}
    U_test <- wilcox.test(rich[,n] ~ rich[,i], data = rich, paired = F, exact = exact)
    z <- abs(qnorm(U_test$p.value/2))
    r <- z/sqrt(nrow(rich))           #Effect size calculation
    tab <- c(U_test$method, n, i, U_test$statistic, U_test$p.value, r)
    mwut[[paste0(i, "_", n)]] <- tab
  }
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

mwut$p.adj.BH.same.var <- NA
mwut$p.adj.bonferroni.same.var <- NA
  
for (i in test.mw){
mwut[mwut$Variable2 == i,]$p.adj.BH.same.var <- p.adjust(mwut[mwut$Variable2 == i,]$`p-value`, method = "BH")
mwut[mwut$Variable2 == i,]$p.adj.bonferroni.same.var <- p.adjust(mwut[mwut$Variable2 == i,]$`p-value`, method = "bonferroni")
}

mwut$p.adj.BH <- p.adjust(mwut$`p-value`, method = "BH")
mwut$p.adj.bonferroni <- p.adjust(mwut$`p-value`, method = "bonferroni")

write.xlsx(list(mwut,rich), file.path(a.stats, "Mann_Whitney_u_test.xlsx"), rowNames = T)

# Wilcoxon-Test 
# Wilcoxon needs paired data
# Inspect the Study_Nr that do not have one of each

table(rich$Lung, rich$Study_Nr) #"LML_017","LML_086"

if (only_tumor == F | (only_direct == F & only_culture == F)){ if (only_tumor == F & only_direct == F & only_culture == F){ test.w = c("Lung", "Isolation")}
  else if (only_tumor == F & (only_direct == T | only_culture == T)){ test.w = "Lung" }else if(only_tumor == T & only_direct == F & only_culture == F){ test.w = "Isolation"}

wt <- list()

for (n in c("Shannon", "InvSimpson")) {
  for (i in test.w) {
    #direct.wt <- direct
   # direct.wt <- direct[rich$Synchronous.tumor == 0,]
    W_test <- wilcox.test(rich[, n] ~ rich[, i], data = rich, paired = T, exact = T)
    z <- abs(qnorm(W_test$p.value/2))
    r <- z/sqrt(nrow(rich))
    
    tab <- c(W_test$method, n, i, W_test$statistic, W_test$p.value, r)
    wt[[paste0(n, "_", i)]] <- tab
  }
}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

wt$p.adj.BH.same.var <- NA
wt$p.adj.bonferroni.same.var <- NA

for (i in test.w){
  wt[wt$Variable2 == i,]$p.adj.BH.same.var <- p.adjust(wt[wt$Variable2 == i,]$`p-value`, method = "BH")
  wt[wt$Variable2 == i,]$p.adj.bonferroni.same.var <- p.adjust(wt[wt$Variable2 == i,]$`p-value`, method = "bonferroni")
}

wt$p.adj.BH <- p.adjust(wt$`p-value`, method = "BH")
wt$p.adj.bonferroni <- p.adjust(wt$`p-value`, method = "bonferroni")

write.xlsx(list(wt,rich), file.path(a.stats, "Wilcoxon_test.xlsx"), rowNames = T)
}

# Kruskall-Wallis
#https://bjoernwalther.com/kruskal-wallis-test-in-r-rechnen/
if (length(selected_Diagnosis) > 1 & "NSCLC" %in% selected_Diagnosis | all_Diagnosis == T) {test.kw = c("Diagnosis", "History.of.smoking", "Histology.NSCLC", "Lobe")
}else if (selected_Diagnosis == "NSCLC" & all_Diagnosis == F) {test.kw = c("History.of.smoking", "Histology.NSCLC", "Lobe", "T", "N", "M")
}else if (selected_Diagnosis == "SCLC" & all_Diagnosis == F) {test.kw = c("Lobe")
}else if (selected_Diagnosis == "Benign" & all_Diagnosis == F)  {test.kw = F}

if (test.kw[1] != F){
kw <- list()

for (n in c("Shannon", "InvSimpson")) {
  for(i in test.kw){
    if(i == "Diagnosis"){rich.kw <- rich[rich$Diagnosis != "V.a. NSCLC",]}else{
      rich.kw <- rich}
    kruskal <- kruskal.test(rich.kw[,n] ~ rich.kw[,i])
    eta_squared <- (kruskal$statistic - 3 + 1)/(nrow(rich.kw) - 3)
    f <- sqrt(eta_squared/(1-eta_squared))
    
    tab <- c(kruskal$method, n, i, kruskal$statistic, kruskal$parameter, kruskal$p.value, f)
    kw[[paste0(n, "_", i)]] <- tab
    
  }}

kw <- as.data.frame(do.call(rbind, kw))
colnames(kw) <- c("Test", "Variable1", "Variable2", "Chi2", "df", "p-value", "Effect size")

kw$p.adj.BH.same.var <- NA
kw$p.adj.bonferroni.same.var <- NA

for (i in test.kw){
  kw[kw$Variable2 == i,]$p.adj.BH.same.var <- p.adjust(kw[kw$Variable2 == i,]$`p-value`, method = "BH")
  kw[kw$Variable2 == i,]$p.adj.bonferroni.same.var <- p.adjust(kw[kw$Variable2 == i,]$`p-value`, method = "bonferroni")
}

kw$p.adj.BH <- p.adjust(kw$`p-value`, method = "BH")
kw$p.adj.bonferroni <- p.adjust(kw$`p-value`, method = "bonferroni")

write.xlsx(list(kw, rich.kw), file = file.path(a.stats, "Kruskal.wallis.xlsx"))

# Post hoc analysis

ph <- list()
for (n in c("Shannon", "InvSimpson")) {
  for(i in test.kw){
    if(i == "Diagnosis"){rich.kw <- rich[rich$Diagnosis != "V.a. NSCLC",]}else{
      rich.kw <- rich}
    p.adj = "BH"
    posthocw <- pairwise.wilcox.test(rich.kw[,n],rich.kw[,i], paired = F, p.adjust=p.adj)
    df <- rbind(posthocw$p.value, c(posthocw$method, posthocw$p.adjust), c(n,i))
    ph[[paste0(n, "_", i)]] <- df
  }}
write.xlsx(ph, file = file.path(a.stats, paste0("Kruskal.wallis_posthoc.MW.",p.adj,".xlsx")), rowNames = T)

# Post hoc analysis
ph <- list()
for (n in c("Shannon", "InvSimpson")) {
  for(i in test.kw){
    if(i == "Diagnosis"){rich.kw <- rich[rich$Diagnosis != "V.a. NSCLC",]}else{
      rich.kw <- rich}
    p.adj = "BH"
    rich.kw <- rich.kw[,c(n,i)]
    colnames(rich.kw) <- c("Variable", "Group")
    posthocd <- rstatix::dunn_test(rich.kw, Variable ~ Group, p.adjust.method = p.adj)
    posthocd$.y. <- paste(n, i, p.adj)
    ph[[paste0(n, "_", i)]] <- posthocd
  }}

write.xlsx(ph, file = file.path(a.stats, paste0("Kruskal.wallis_posthoc.dunn.",p.adj,".xlsx")), rowNames = T)
}
#### End ####

#### Convert sample data to data frame ####

# Convert to long data frame for the facets
rich_long <- rich %>%
  pivot_longer(cols = c(Shannon, InvSimpson),
               names_to = "alpha_measure",
               values_to = "alpha_value")

rich_long$alpha_measure <- factor(rich_long$alpha_measure, levels = c("Shannon", "InvSimpson"))

#### End ####

#### Visualization ####
text_size = 30
thm <- theme(text = element_text(size = text_size),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             legend.position = "bottom")

thm.x <- theme(text = element_text(size = text_size),
               axis.text.x = element_text(angle = -30, hjust = 0),
               legend.position = "bottom")
#### Plots ####

# Lung
if (only_tumor == F){
P1 <- ggplot(rich_long, aes(x = Lung, y = alpha_value, color = Lung, fill = Lung)) 

P1.1 <- P1 + geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = lung_col) + scale_color_manual(values = lung_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL) + thm

P1.1

ggsave(plot = P1.1, file.path(a.plots, "01-Lung.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P1.1, file.path(a.plots, "01-Lung.svg"), width = 10, height = 8, dpi = 300)
}

# Diagnosis
if (all_Diagnosis == T | length(selected_Diagnosis) > 1){
P2 <- ggplot(rich_long, aes(x = Diagnosis, y = alpha_value, color = Diagnosis, fill = Diagnosis)) 

P2.1 <- P2 + geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = diagnosis_col) + scale_color_manual(values = diagnosis_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "Diagnosis", fill = "Diagnosis") + thm

P2.1

ggsave(plot = P2.1, file.path(a.plots, "02-Diagnosis.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P2.1, file.path(a.plots, "02-Diagnosis.svg"), width = 10, height = 8, dpi = 300)
}

# Smoking
P4 <- ggplot(subset(rich_long, !is.na(History.of.smoking)), aes(x = History.of.smoking, y = alpha_value, color = History.of.smoking, fill = History.of.smoking)) 

P4.1 <- P4 + geom_point() + geom_boxplot(alpha = 0.5, na.rm = T) +
  scale_fill_manual(values = smokerh_col) + scale_color_manual(values = smokerh_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "History of smoking", fill = "History of smoking") + thm 

P4.1

ggsave(plot = P4.1, file.path(a.plots, "04-History.of.smoking.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P4.1, file.path(a.plots, "04-History.of.smoking.svg"), width = 10, height = 8, dpi = 300)

P41 <- ggplot(subset(rich_long, !is.na(History.of.smoking.y.n)), aes(x = History.of.smoking.y.n, y = alpha_value, color = History.of.smoking.y.n, fill = History.of.smoking.y.n)) 

P41.1 <- P41 + geom_point() + geom_boxplot(alpha = 0.5, na.rm = T) +
  scale_fill_manual(values = smoker_col) + scale_color_manual(values = smoker_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "History of smoking", fill = "History of smoking") + thm 

P41.1

ggsave(plot = P41.1, file.path(a.plots, "04-History.of.smoking.y.nes.no.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P41.1, file.path(a.plots, "04-History.of.smoking.y.nes.no.svg"), width = 10, height = 8, dpi = 300)

P5 <- ggplot(subset(rich_long, !is.na(Active.smoker)), aes(x = Active.smoker, y = alpha_value, color = Active.smoker, fill = Active.smoker)) 

P5.1 <- P5 + geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = smoker_col) + scale_color_manual(values = smoker_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "Active smoker", fill = "Active smoker") + thm

P5.1

ggsave(plot = P5.1, file.path(a.plots, "05-Active.smoker.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P5.1, file.path(a.plots, "05-Active.smoker.svg"), width = 10, height = 8, dpi = 300)

# Histology
if (all_Diagnosis == F & length(selected_Diagnosis) == 1){if (selected_Diagnosis == "NSCLC") {
P3 <- ggplot(subset(rich_long, !is.na(Histology.NSCLC)), aes(x = Histology.NSCLC, y = alpha_value, color = Histology.NSCLC, fill = Histology.NSCLC)) 

P3.1 <- P3 + geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = histology_col) + scale_color_manual(values = histology_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL, color = "Histology", fill = "Histology") + thm +
  guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))

P3.1

ggsave(plot = P3.1, file.path(a.plots, "03-Histology.NSCLC.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P3.1, file.path(a.plots, "03-Histology.NSCLC.svg"), width = 10, height = 8, dpi = 300)
}}

# Lung + Diagnosis
if (only_tumor == F){
P12 <- ggplot(rich_long, aes(x = Lung, y = alpha_value, color = Diagnosis, fill = Diagnosis)) 

P12.1 <- P12 + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = diagnosis_col) + scale_color_manual(values = diagnosis_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "Lung", color = "Diagnosis", fill = "Diagnosis") + thm.x

P12.1

ggsave(plot = P12.1, file.path(a.plots, "012-Diagnosis_lung.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P12.1, file.path(a.plots, "012-Diagnosis_lung.svg"), width = 10, height = 8, dpi = 300)

P13 <- ggplot(rich_long, aes(x = Diagnosis, y = alpha_value, color = Lung, fill = Lung)) 

P13.1 <- P13 + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = lung_col) + scale_color_manual(values = lung_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "Diagnosis", color = "Lung", fill = "Lung") + thm.x

P13.1

ggsave(plot = P13.1, file.path(a.plots, "013-Lung_Diagnosis.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P13.1, file.path(a.plots, "013-Lung_Diagnosis.svg"), width = 10, height = 8, dpi = 300)
}

# Side
P14 <- ggplot(rich_long, aes(x = Side, y = alpha_value, color = Side, fill = Side)) 

P14.1 <- P14 + geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL) + thm

P14.1

ggsave(plot = P14.1, file.path(a.plots, "014-Side.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P14.1, file.path(a.plots, "014-Side.svg"), width = 10, height = 8, dpi = 300)

# Lobe
P15 <- ggplot(subset(rich_long, !is.na(Lobe)), aes(x = Lobe, y = alpha_value, color = Lobe, fill = Lobe)) 

P15.1 <- P15 + geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL) + thm

P15.1

ggsave(plot = P15.1, file.path(a.plots, "014-Lobe.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P15.1, file.path(a.plots, "014-Lobe.svg"), width = 10, height = 8, dpi = 300)

# Isolation
if (only_direct == F & only_culture == F){
P6 <- ggplot(rich_long, aes(x = Isolation, y = alpha_value, color = Isolation, fill = Isolation)) 

P6.1 <- P6 + geom_point() + geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = isolation_col) + scale_color_manual(values = isolation_col) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL) + thm +
  guides(color = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))

P6.1

ggsave(plot = P6.1, file.path(a.plots, "06-Isolation.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P6.1, file.path(a.plots, "06-Isolation.svg"), width = 10, height = 8, dpi = 300)

P61 <- ggplot(rich_long, aes(x = Lung, y = alpha_value, shape = Isolation, color = Diagnosis)) 

P61.1 <- P61 + geom_boxplot_pattern(aes(pattern = Isolation, fill = Diagnosis), pattern_fill = "white", alpha = 0.3, pattern_spacing = 0.02, outlier.shape=NA) +
  scale_fill_manual(values = diagnosis_col) + scale_color_manual(values = diagnosis_col) +
  scale_pattern_manual(values = c("stripe","none"), name= "Isolation") +
  labs(y = NULL, x = NULL, color = "Diagnosis", fill = "Diagnosis", shape = "Isolation") + thm.x +
  geom_point(position=position_jitterdodge(jitter.width = 0),alpha=0.6, size = 3) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  guides(pattern = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))

P61.1

ggsave(plot = P61.1, file.path(a.plots, "061-Lung_Diagnosis_isolation.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P61.1, file.path(a.plots, "061-Lung_Diagnosis_isolation.svg"), width = 10, height = 8, dpi = 300)
}

# Sex
P7 <- ggplot(rich_long, aes(x = Sex, y = alpha_value, color = Sex, fill = Sex)) 

P7.1 <- P7 + geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL) + thm

P7.1

ggsave(plot = P7.1, file.path(a.plots, "07-Sex.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P7.1, file.path(a.plots, "07-Sex.svg"), width = 10, height = 8, dpi = 300)

# Synchronous tumor
if (only_non_synchronous == F){
P16 <- ggplot(rich_long, aes(x = Synchronous.tumor, y = alpha_value, color = Synchronous.tumor, fill = Synchronous.tumor)) 

P16.1 <- P16 + geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = NULL) + thm +
  guides(color = guide_legend(nrow = 3), fill = guide_legend(nrow = 3))

P16.1

ggsave(plot = P16.1, file.path(a.plots, "016-Synchronous.tumor.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P16.1, file.path(a.plots, "016-Synchronous.tumor.svg"), width = 10, height = 8, dpi = 300)
}

# T
P8 <- ggplot(subset(rich_long, !is.na(T)), aes(x = T, y = alpha_value)) 
  
P8.1 <- P8 + geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "T") + theme(text = element_text(size = text_size))

P8.1
  
ggsave(plot = P8.1, file.path(a.plots, "08.1-T.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P8.1, file.path(a.plots, "08.1-T.svg"), width = 10, height = 8, dpi = 300)

# N
P8.2 <- ggplot(subset(rich_long, !is.na(N)), aes(x = N, y = alpha_value)) 

P8.2 <- P8.2 + geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "N") + theme(text = element_text(size = text_size))

P8.2

ggsave(plot = P8.2, file.path(a.plots, "08.2-N.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P8.2, file.path(a.plots, "08.2-N.svg"), width = 10, height = 8, dpi = 300)

# M
P8.3 <- ggplot(subset(rich_long, !is.na(M)), aes(x = M, y = alpha_value)) 

P8.3 <- P8.3 + geom_point() + geom_boxplot(alpha = 0.5) +
  facet_wrap(. ~ alpha_measure, scales = "free_y") +
  labs(y = NULL, x = "M") + theme(text = element_text(size = text_size))

P8.3

ggsave(plot = P8.3, file.path(a.plots, "08.1-M.tiff"), width = 10, height = 8, dpi = 300)
ggsave(plot = P8.3, file.path(a.plots, "08.1-M.svg"), width = 10, height = 8, dpi = 300)

#### End ####
}


#### Alpha diversity following structure of the questions ####


#### Directed questions ####

physeq = physeq_o
physeq = subset_samples(physeq, sample_names(physeq) %in% c(rownames(sample_info[!sample_info$Sample_or_Control %in% c("Negative control sample", "Positive control sample"),]))) 
physeq = subset_samples(physeq, sample_data(physeq)$Isolation == "Direct_Isolation") 

#### Convert sample data to data frame ####

rich <- data.frame(sample_data(physeq))
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

#### Create directories for results ####

a.stats <- file.path(res.dir, "2.Alpha_stats", "Directed_questions")
a.plots <- file.path(res.dir, "3.Alpha_plots", "Directed_questions")

dir.create(a.stats, recursive = T)
dir.create(a.plots, recursive = T)

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
#### End ####

#### Q1: Is there a difference between Diseased vs the parallel Normal Lung? ####

# In all diagnosis
wt <- list()

for (n in c("Shannon", "InvSimpson")) {
    W_test <- wilcox.test(rich[, n] ~ rich[, "Lung"], data = rich, paired = T, exact = T)
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
  W_test <- wilcox.test(rich[rich$Diagnosis == i, n] ~ rich[rich$Diagnosis == i, "Lung"], data = rich, paired = T, exact = T)
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
  W_test <- wilcox.test(rich[rich$Synchronous.tumor == "No", n] ~ rich[rich$Synchronous.tumor == "No", "Lung"], data = rich, paired = T, exact = T)
  z <- abs(qnorm(W_test$p.value/2))
  r <- z/sqrt(nrow(rich[rich$Synchronous.tumor == "No", ]))
  
  tab <- c(W_test$method, n, "Lung", W_test$statistic, W_test$p.value, r)
  wt[[n]] <- tab
}

wt <- as.data.frame(do.call(rbind, wt))
colnames(wt) <- c("Test", "Variable1", "Variable2", "Statistic", "p-value", "Effect size")

wt$p.adj.BH <- p.adjust(wt$`p-value`, method = "BH")
wt$p.adj.bonferroni <- p.adjust(wt$`p-value`, method = "bonferroni")

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
    W_test <- wilcox.test(rich[rich$Diagnosis == i & rich$Synchronous.tumor == "No", n] ~
                            rich[rich$Diagnosis == i & rich$Synchronous.tumor == "No", "Lung"], data = rich, paired = T, exact = T)
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

write.xlsx(list(kw,rich), file.path(a.stats, "Q3.xlsx"), rowNames = T)

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
    U_test <- wilcox.test(rich.NSCLC[,n] ~ rich.NSCLC[,"Histology.NSCLC"], data = rich.NSCLC, paired = F, exact = exact)
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
                          rich.NSCLC[rich.NSCLC$History.of.smoking.y.n == "Yes","Histology.NSCLC"], data = rich.NSCLC, paired = F, exact = exact)
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
  U_test <- wilcox.test(rich.NSCLC[,n] ~ rich.NSCLC[,"History.of.smoking.y.n"], data = rich.NSCLC, paired = F, exact = exact)
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

# Benign
# Create data frame with benign and diseased lung
rich.ben <- rich[rich$Diagnosis == "Benign"  & rich$Lung == "Diseased",]
rich_long.ben <- rich_long[rich_long$Diagnosis == "Benign" & rich_long$Lung == "Diseased",]

mwut <- list()

for (n in c("Shannon", "InvSimpson")){
  if(nrow(rich.ben) <51){exact = T}else{exact = F}
  U_test <- wilcox.test(rich.ben[,n] ~ rich.ben[,"History.of.smoking.y.n"], data = rich.ben, paired = F, exact = exact)
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
  U_test <- wilcox.test(rich.NSCLC[,n] ~ rich.NSCLC[,"M"], data = rich.NSCLC, paired = F, exact = exact)
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

