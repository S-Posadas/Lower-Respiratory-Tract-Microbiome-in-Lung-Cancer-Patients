# Lower Respiratory Tract Microbiome in Lung Cancer Patients

Analysis code for the 16S rRNA gene sequencing (V1–V2) of paired bronchial
lavage samples from patients with suspected lung cancer. The samples were
collected from the tumour-bearing and the contralateral lung during diagnostic
bronchoscopy, in a single-centre prospective study. Each sample was sequenced
directly and after culture enrichment.
>
> **Raw reads:** European Nucleotide Archive (ENA), project
> [PRJEB94648](https://www.ebi.ac.uk/ena/browser/view/PRJEB94648).

## Workflow

| Step | Script | What it does |
|---|---|---|
| Quality control | `QC` | FastQC and MultiQC on the raw and trimmed reads |
| Trimming | `cutadapt` | Removes the 27F/338R primers and Nextera adapters |
| ASV inference | `1.DADA2_direct.R`, `1.DADA2_culture.R` | DADA2 per sequencing batch in two passes (naive pass, then all-batch priors), SILVA taxonomy assignment, contaminants removal with decontam, phyloseq objects as output, ASV alignment for phylogenetic tree (ClustalW). |
| Phylogenetic tree | `fasttree` | Approximate maximum-likelihood tree with FastTree |
| Species refinement | `blast`, `filter_blast.R` | BLASTn against NCBI 16S rRNA, hits with ≥ 99 % identity and ≥ 90 % query coverage are kept, species are added to the SILVA taxonomy |
| Batch correction | `2.MMUPHin.R` | Sequencing-batch effects correction with MMUPHin |
| Combined object | `merge_phyloseq.R` | Combination of the direct and culture-enriched phyloseq objects |
| Downstream analyses | see below | Diversity, composition, differential abundance, functional prediction |

### Downstream analyses

| Script | Analysis |
|---|---|
| `paired_lungs.R` | Tumour-bearing vs contralateral lung: within- vs between-patient distances, percentile ranks and composition |
| `alpha.R` | Alpha diversity (Shannon, inverse Simpson). Paired Wilcoxon tests for lung side, and permutation ANCOVA adjusted for sex, age and smoking for the other questions |
| `beta.R` | Beta diversity (Bray–Curtis, weighted and unweighted UniFrac): NMDS/PCoA, PERMANOVA (adonis2) and betadisper |
| `taxa_plots.R` | Box, bar, ternary and mirrored composition plots of the most abundant taxa |
| `maaslin.R` | MaAsLin2 differential abundance, leave-one-out sensitivity analyses, forest and Venn plots |
| `lefse.R`, `lefser_functions.R` | LEfSe (lefser, modified to return the Kruskal–Wallis p-values) |
| `halla`, `halla_input.R`, `halla_plots.R` | HAllA association between ASVs and diagnosis or smoking |
| `picrust2`, `picrust2_input.R`, `picrust2.R` | PICRUSt2 MetaCyc pathways, compared with LinDA (ggpicrust2) |
| `core_and_culture.R` | Core genera, and direct vs culture-enriched samples |
| `delta.R` | Shift in abundance after culture enrichment vs N stage, M stage and smoking |
| `circle_plot.R` | Circular overview of the NSCLC cohort |
| `mock_database_comparison.R` | Mock community: SILVA vs eHOMD vs GTDB |

