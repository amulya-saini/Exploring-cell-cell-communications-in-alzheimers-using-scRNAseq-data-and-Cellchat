# Exploring cell-cell communications in Alzheimer's using scRNA seq data and Cellchat

## Project Overview
Alzheimer's Disease (AD) stands as a formidable challenge in neuroscience due to its intricate pathophysiology and devastating impact on cognitive function. To decode the complexities of this neurodegenerative disorder, this project harnesses the power of single-cell RNA sequencing (scRNA-seq) technology coupled with computational tools for intercellular communication analysis.

## Background
AD is marked by progressive cognitive decline and is notorious for its heterogeneous cellular landscape within the brain. Traditional bulk RNA sequencing struggles to capture the nuanced interactions between various cell types, necessitating the adoption of single-cell technologies. Leveraging recent advancements in scRNA-seq and computational methodologies like CellChat, I delve into the intricate web of cell-to-cell interactions within the AD brain.

## Objectives
The primary goal is to decipher the transcriptomic alterations and intercellular communication patterns among neuronal and non-neuronal cells, including microglia and astrocytes, in the context of AD pathology. By mining publicly available scRNA-seq data (GEO Accession Number: GSE227157), I aim to uncover novel insights into the mechanisms driving disease onset and progression.

## Data Description
1. Organism: Mus musculus
2. Experiment Type: Expression profiling by high throughput sequencing
3. Summary: Daily intraperitoneal (i.p.) injections of exogenous APC were administered to wild-type (WT) and 5xFAD mice for 5 months, with daily saline injections given as controls to WT and 5xFAD mice.
4. Web Link: [DOI: 10.18632/aging.205624](https://doi.org/10.18632/aging.205624)
5. Citation: Fatmi MK, Wang H, Slotabec L, Wen C et al. Single-Cell RNA-seq reveals transcriptomic modulation of Alzheimer's disease by activated protein C. Aging (Albany NY) 2024 Feb 21;16(4):3137-3159. PMID: 38385967
6. Submission Date: Mar 12, 2023
7. Last Update Date: Feb 25, 2024

# Installation Instructions

## R Packages

# Install the following R packages:

```
install.packages(c("devtools", "patchwork", "Matrix", "Seurat", "dplyr", "tidyr", "circlize", "digest", "NMF", "presto", "BiocManager"))
```

# Install additional packages from GitHub:

```
devtools::install_github("jinworks/CellChat")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jokergoo/circlize")
devtools::install_github("immunogenomics/presto")
```

# Analysis Pipeline

## Data Preprocessing:
- Preprocess scRNA-seq data, including quality control, normalization, and dimensionality reduction.
- Assess quality control metrics, including the number of expressed features, total counts, and mitochondrial gene content.

## Dimensionality Reduction:
- Use Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP) for dimensionality reduction and visualization.

## Clustering:
- Cluster cells based on transcriptomic profiles to identify distinct cell populations.

## Differential Expression Analysis:
- Perform differential expression analysis to identify marker genes for each cell cluster.

## Cell Type Annotation:
- Identify cell types using clustering and marker gene expression analysis.

## Cell-Cell Communication Analysis:
- Utilize CellChat to investigate cell-cell communication networks and signaling pathways affected by APC treatment.

## Visualization:
- Visualize transcriptomic changes, cell clusters, and communication networks using dimensionality reduction techniques, heatmaps, and network plots.

## Interpretation:
- Interpret the results to gain insights into the molecular mechanisms underlying APC-mediated modulation of AD.


# Usage
To replicate the analysis or explore the data further, follow the steps outlined in the provided R scripts. Ensure that the required R packages are installed and the data files are correctly sourced.

# References
- Fatmi MK, Wang H, Slotabec L, Wen C et al. Single-Cell RNA-seq reveals transcriptomic modulation of Alzheimer's disease by activated protein C. Aging (Albany NY) 2024 Feb 21;16(4):3137-3159. PMID: 38385967
- Jin, S., Guerrero‐Juarez, C. F., Zhang, L., Chang, I., Ramos, R., Kuan, C., Myung, P., & Plikus, M. V. (2021). Inference and analysis of cell-cell communication using CellChat. Nature Communications, 12(1). https://doi.org/10.1038/s41467-021-21246-9
- [CellChat V2](https://github.com/jinworks/CellChat) 
