# Lower-Respiratory-Tract-Microbiome-in-Lung-Cancer-Patients

This repository contains code and documentation for the analysis of 16S rRNA gene sequencing data generated from paired bronchial lavage samples collected from patients with suspected lung cancer. The data originate from a single-centre, prospective study including 119 patients, where paired samples were obtained from both lungs during diagnostic bronchoscopy.

## Analysis Workflow

The analysis pipeline consists of the following steps:

### 1. Quality Control
- FastQC and MultiQC
- Adapter trimming with cutadapt

### 2. Sequence Processing
- Quality filtering and trimming, and denoising of reads using the **DADA2** pipeine
- Amplicon sequence variant (ASV) inference using a predefined list of expected sequences to improve sensitivity

### 3. Taxonomic Assignment
- Database selection based on optimal reconstruction of a sequenced mock community
- Taxonomic classification using the **SILVA** reference database
- Species-level refinement via **BLASTn** against the **NCBI 16S rRNA** database

### 4. Contamination and Batch Correction
- Identification and removal of contaminant features using **Decontam** (prevalence-based method)
- Correction for sequencing run–associated batch effects using **MMUPHin**

### 5. Downstream Analysis
- Alpha and beta diversity analysis using **phyloseq**
- Differential abundance testing with **MaAsLin2** and a modified **Rhea** pipeline
- Identification of discriminant taxa using **lefser**

### 6. Visualization
- Data visualization in **R** using **phyloseq**, **ggplot2**, **ggordiplots**
