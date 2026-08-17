# Multiomics Reveals Compartmentalized Immune Responses and Tissue-Vascular Signatures in Lyme Disease

Preprint: [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.05.27.656061v1)  
Latest data release: [Zenodo](https://doi.org/10.5281/zenodo.15395686) (concept DOI)

## Summary

Lyme disease (LD) is growing in incidence, with nearly 500,000 cases diagnosed annually in the United States. Despite treatment, some patients experience persistent symptoms. The immune mechanisms underlying LD remain poorly understood. We conducted a multiomic longitudinal analysis of 49 LD patients and matched controls, integrating plasma proteomics, metabolomics, PBMC immunophenotyping, and a meta-analysis of skin lesions, with key longitudinal analyses performed in an antibiotic-naïve subset. We identified compartmentalized immune responses in acute LD, with coordinated alterations in circulating plasma proteins and metabolites linked to endothelial barrier stability, metabolic reprogramming, and symptom severity, predominantly reflecting tissue and vascular immune processes at the site of infection. In contrast, PBMC responses were limited, revealing a disconnect between localized tissue responses and systemic immunity. These findings provide novel insights into LD pathophysiology and highlight the potential for diagnostics leveraging tissue and vascular immune markers detectable in blood. They also provide a resource for biomarker discovery and predictive modeling to improve LD management. 

Keywords: Lyme disease, *Borrelia burgdorferi*, PBMC, skin, endothelium, multiomics, immune response.

## About

This repository contains the scripts for generating the manuscript figures and the preprocessing workflows. The companion Zenodo release supplies the objects used to generate figures, selected intermediate caches, and the metadata required by those scripts, prepared in accordance with the HIPAA Safe Harbor standard.

Together, these resources are intended to reproduce the analyses and figures reported in this study and are not intended as a general purpose analysis package for arbitrary user-supplied datasets.

## Overall categories

All the scripts and data are organized into four repeating categories:

1. Proteomics & metabolomics
2. Internal PBMC scRNA-seq
3. Public skin RNA
4. Integrative modeling

## Asset allocation and directory hierarchy

The assets are allocated as follows:

- GitHub hosts the scripts for analysis and  rendering the figures.
  - Figure rendering scripts are stored as Quarto markdown files (`.qmd`) in the repository root.
  - Supporting analysis pipelines are organized in the `scripts/` subdirectory.
- Zenodo hosts three downloadable bundles.
  - `Proteomics_Metabolomics_scRNAseq_Raw.zip` contains one top-level `raw/` directory.
  - `Proteomics_Metabolomics_scRNAseq_Processed.zip` contains the top-level `processed/`, `intermediate/`, and `metadata/` directories.
  - `Flow_Data_and_Preprocessing.zip` (89 GB) contains the raw flow cytometry data and its preprocessing workflow. It is self-contained and is not extracted into the hierarchy below.

The two archives other than the flow cytometry archive are complementary. Extract both directly into the repository's `data/` directory and not into a subdirectory or into an automatically named directory created by the default extraction behavior of your operating system. The directory tree shown below shows the paths and folder names compatible with the latest release of the data and scripts. Older paths and folder names from previous releases are not supported.

```
Lyme/
├── README.md                             ----On GitHub----
├── fig_1.qmd ... fig_7.qmd               # Main figure scripts
├── fig_s*.qmd                            # Supplementary figure scripts
│
├── scripts/                              ----On GitHub----
│   ├── 02_internal_pbmc_scrna/           # PBMC scRNA-seq of our cohort
│   │   └── 01_pbmc_preprocessing.Rmd     # PBMC preprocessing and downstream analysis
│   └── 03_public_skin_rna/               # Public skin RNA preprocessing
│       ├── 01a_prepare_public_expression_data.R
│       ├── 01b_integrate_and_annotate_GSE169440.R
│       └── 02_test_GSE169440_subpopulations.R
│
├── data/                               
│   ├── raw/                              ----On Zenodo----
│   │   ├── 01_proteomics_metabolomics/   # OlinkPreprocessed.RData (cohort multiomics input)
│   │   ├── 02_internal_pbmc_scrna/       # 10x/<sample>/ triplets
│   │   ├── 03_public_skin_rna/           # GSE169440/<sample>_GEX_HHT_cellranger/filtered_feature_bc_matrix
│   │   └── 04_integrative_modeling/      # Published tables used for comparison panels
│   │
│   ├── processed/                        ----On Zenodo----
│   │   ├── 01_proteomics_metabolomics/   # Data.RData and matrices used to generate figures
│   │   ├── 02_internal_pbmc_scrna/       # Integrated and annotated Seurat objects; cache of differential expression results over time
│   │   ├── 03_public_skin_rna/           # Preprocessed public skin/PBMC transcriptomic data and abundance tables used to generate figures
│   │   └── 04_integrative_modeling/      # Cached model evaluation tables
│   │
│   ├── metadata/                         ----On Zenodo----
│   │   ├── sampleData.csv                # Primary cohort clinical metadata
│   │   └── sc_sampleData.csv             # scRNA-seq map between libraries and patients
│   │
│   └── intermediate/                     ----On Zenodo----
│       ├── 01_proteomics_metabolomics/   # FELLA graph and PageRank caches
│       ├── 02_internal_pbmc_scrna/       # SCE, Monaco prediction, and pseudobulk caches
│       ├── 03_public_skin_rna/           # ncells.RData, SubCluster_Labels.RData, Public Data Limma.RData
│       ├── 04_integrative_modeling/      # T1_T3_Multiomics_Community_Results.xlsx and cached model/figure inputs
│       └── flow_gating/                  # bcell, tcell, monocyte, dcnk.RData
```

### Preprocessing scripts and their scope

The files in the `scripts/` directory document the preprocessing workflows and analytical choices used to create intermediate and processed objects. They are not intended to regenerate every deposited cache from raw data.

The script
`scripts/02_internal_pbmc_scrna/01_pbmc_preprocessing.Rmd`
documents the import, integration, annotation, and downstream analyses of PBMC sequencing data from our cohort.

The three scripts inside `scripts/03_public_skin_rna/` prepare public transcriptomic data and integrate GSE169440 to annotate cell subpopulations in unaffected skin and erythema migrans lesions.

## System requirements

The following environment was tested:

- Windows 11 x64 (build 22631).
- R 4.3.
- Quarto 1.9.37.

The R packages used across the figure and preprocessing scripts are listed below. Individual scripts use subsets of this list.

- Core and reporting: `pacman`, `here`, `fs`, `knitr`, `ragg`, `systemfonts`, `rsvg`, `rasterpdf`, `DiagrammeR`, `DiagrammeRsvg`, `datapasta`, `openxlsx`, `officer`, `flextable`, `DT`.
- Data handling: `tidyverse`, `dplyr`, `tidyr`, `tibble`, `readr`, `readxl`, `purrr`, `forcats`, `stringr`, `rlang`, `magrittr`, `janitor`, `plyr`, `broom`, `matrixStats`, `conflicted`.
- Plotting and networks: `ggplot2`, `ggrepel`, `ggstance`, `gghighlight`, `ggnewscale`, `geomtextpath`, `ggbump`, `ggforce`, `ggtext`, `ggthemes`, `ggprism`, `ggpubr`, `ggplotify`, `ggsignif`, `patchwork`, `cowplot`, `camcorder`, `colorspace`, `RColorBrewer`, `circlize`, `ComplexHeatmap`, `pheatmap`, `corrplot`, `factoextra`, `ggraph`, `scales`, `igraph`, `UpSetR`, `VennDiagram`.
- Statistics and modeling: `Hmisc`, `psych`, `emmeans`, `lmtest`, `lme4`, `lmerTest`, `MASS`, `mice`, `impute`, `caret`, `e1071`, `glmnet`, `glmmLasso`, `pROC`, `stringdist`, `Matrix.utils`, `WGCNA`, `progress`, `future`, `furrr`, `progressr`.
- Bioinformatics: `DESeq2`, `limma`, `edgeR`, `sva`, `BiocParallel`, `AnnotationDbi`, `org.Hs.eg.db`, `clusterProfiler`, `enrichplot`, `enrichR`, `DOSE`, `msigdbr`, `FELLA`, `KEGGREST`, `GEOquery`, `DropletUtils`, `multtest`, `metap`, `scCATCH`, `SingleR`, `celldex`, `ROntoTools`, `lemur`, `Seurat`, `SingleCellExperiment`, `SummarizedExperiment`, `flowCore`, `flowStats`, `xCell`, `venn`.
- External lookups: `rentrez`, `webchem`.

Exact package versions are recorded by `sessionInfo()` at the end of each figure script and preprocessing workflow.

### Hardware and resource requirements

No specialized hardware is required for routine figure rendering or for the representative Figure 2 demo.
The main exceptions are the full flow cytometry analysis, whose archive is approximately 89 GB before extraction, and the optional full refit in `fig_s5.qmd`, which can use up to 24 CPU cores. Routine rendering uses the supplied cached results unless the parameter that enables a full refit is selected.

## Quick start

No standalone software package was developed for this study. To reproduce the analysis, combine the code from GitHub with the data from Zenodo:

1. Install R 4.3, Quarto 1.9.37, and the required R packages listed above. Use CRAN or Bioconductor where available and the package's official source repository otherwise. RStudio is optional. Typical installation time is approximately 30 minutes on a normal desktop computer, excluding data download time.

2. Clone the repository:

   ```text
   git clone https://github.com/subramanian-group/Lyme.git
   cd Lyme
   ```

3. Download `Proteomics_Metabolomics_scRNAseq_Raw.zip` and `Proteomics_Metabolomics_scRNAseq_Processed.zip` from the [latest Zenodo data release](https://doi.org/10.5281/zenodo.15395686). Download `Flow_Data_and_Preprocessing.zip` only when rendering Figure 6d, e and Supplementary Figure 7, which require the raw flow cytometry data and its standalone preprocessing workflow.

4. Extract `Proteomics_Metabolomics_scRNAseq_Raw.zip` and `Proteomics_Metabolomics_scRNAseq_Processed.zip` directly into `data/`. After extraction, the repository must contain `data/raw/`, `data/processed/`, `data/intermediate/`, and `data/metadata/`; an extra wrapper directory indicates an incorrect extraction. Keep the flow cytometry archive separate. The four small derived gating objects used by the figure scripts, except `fig_6.qmd` and `fig_s7.qmd`, are already supplied under `data/intermediate/flow_gating/`.

   `fig_6.qmd` and `fig_s7.qmd` also use assets from the standalone flow cytometry archive. When rendering Figure 6d, e and Supplementary Figure 7, place the files extracted from `Flow_Data_and_Preprocessing.zip` into a `Flow/` directory under `data/`.

   Download and extraction of the Zenodo archives takes approximately 30 minutes on a 1-Gbps connection.

5. From a system terminal opened in the repository root, run:

   ```text
   quarto render fig_2.qmd
   ```

   Use PowerShell or Command Prompt on Windows, or Bash or Zsh on macOS/Linux. The Terminal tab in RStudio may also be used. Alternatively, open `fig_2.qmd` in RStudio and select **Render**.

### Representative demo

Rendering Figure 2 from `fig_2.qmd` uses only one small deposited input, `data/raw/01_proteomics_metabolomics/OlinkPreprocessed.RData` (approximately 1.3 MB), and does not require the large archives containing single-cell or flow cytometry data.

Expected Figure 2 panel outputs include:

- `fig_2_files/figure_02a_pairwise_protein_differential_expression_heatmap.png`
- `fig_2_files/figure_02b_protein_gsea_pathway_enrichment_heatmap.png`
- `fig_2_files/figure_02c_differential_protein_zscore_heatmap.png`
- `fig_2_files/figure_02d_angiopoietin_longitudinal_boxplots.png`
- `fig_2_files/figure_02e_protein_clinical_correlation_edge_bundle.png`

Typical demo runtime is under 20 minutes.

## License and versioning

The code in this repository is licensed under the [MIT License](LICENSE). The accompanying Zenodo dataset is licensed separately under CC BY 4.0.

## Citation

Rostomily C et al. Multiomics reveals compartmentalized immune responses and tissue-vascular signatures in Lyme disease. bioRxiv. 2025. https://doi.org/10.1101/2025.05.27.656061
