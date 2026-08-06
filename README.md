# Multiomics Reveals Compartmentalized Immune Responses and Tissue-Vascular Signatures in Lyme Disease

Preprint: [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.05.27.656061v1)  
Latest data release: [Zenodo](https://doi.org/10.5281/zenodo.15395686) (concept DOI)

## Summary

Lyme disease (LD) is growing in incidence, with nearly 500,000 cases diagnosed annually in the United States. Despite treatment, some patients experience persistent symptoms. The immune mechanisms underlying LD remain poorly understood. We conducted a multiomic longitudinal analysis of 49 LD patients and matched controls, integrating plasma proteomics, metabolomics, PBMC immunophenotyping, and a meta-analysis of skin lesions, with key longitudinal analyses performed in an antibiotic-naïve subset. We identified compartmentalized immune responses in acute LD, with coordinated alterations in circulating plasma proteins and metabolites linked to endothelial barrier stability, metabolic reprogramming, and symptom severity, predominantly reflecting tissue and vascular immune processes at the site of infection. In contrast, PBMC responses were limited, revealing a disconnect between localized tissue responses and systemic immunity. These findings provide novel insights into LD pathophysiology and highlight the potential for diagnostics leveraging tissue and vascular immune markers detectable in blood. They also provide a resource for biomarker discovery and predictive modeling to improve LD management. 

Keywords: Lyme disease, *Borrelia burgdorferi*, PBMC, skin, endothelium, multiomics, immune response.

## About

This repository contains the scripts for generating the manuscript figures and the preprocessing workflows. The companion Zenodo release supplies the figure-ready objects, selected intermediate caches, and the HIPAA safe-harbor compliant metadata required by those scripts.

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

The two non-flow cytometry archives are complementary. Extract both directly into the repository's `data/` directory and not into a subdirectory or into an automatically-named directory created by the default extraction behavior of your operating system. The directory tree shown below shows the paths and folder names compatible with the latest release of the data and scripts. Older paths and folder names from previous releases are not supported.

```
Lyme/
├── README.md                             ----On GitHub----
├── fig_1.qmd ... fig_7.qmd               # Main figure scripts
├── fig_s*.qmd                            # Supplementary figure scripts
│
├── scripts/                              ----On GitHub----
│   ├── 02_internal_pbmc_scrna/           # PBMC scRNA-seq of our cohort
│   │   └── lyme.Rmd                      # PBMC preprocessing and downstream-analysis source
│   └── 03_public_skin_rna/               # Public skin RNA preprocessing
│       └── scRNAseq EM.R                 # GSE169440 integration and broad annotation
│
├── data/                               
│   ├── raw/                              ----On Zenodo----
│   │   ├── 01_proteomics_metabolomics/   # OlinkPreprocessed.RData (cohort multiomics input)
│   │   ├── 02_internal_pbmc_scrna/       # 10x/<sample>/ triplets
│   │   ├── 03_public_skin_rna/           # GSE169440/<sample>_GEX_HHT_cellranger/filtered_feature_bc_matrix
│   │   └── 04_integrative_modeling/      # Published comparison-panel tables
│   │
│   ├── processed/                        ----On Zenodo----
│   │   ├── 01_proteomics_metabolomics/   # Data.RData and figure-ready matrices
│   │   ├── 02_internal_pbmc_scrna/       # Integrated and annotated Seurat objects; time-DE cache
│   │   ├── 03_public_skin_rna/           # Preprocessed public skin/PBMC transcriptomic data and figure-ready abundance tables
│   │   └── 04_integrative_modeling/      # Cached model evaluation tables
│   │
│   ├── metadata/                         ----On Zenodo----
│   │   ├── sampleData.csv                # Primary cohort clinical metadata
│   │   └── sc_sampleData.csv             # scRNA-seq library-to-patient map
│   │
│   └── intermediate/                     ----On Zenodo----
│       ├── 01_proteomics_metabolomics/   # FELLA graph and PageRank caches
│       ├── 02_internal_pbmc_scrna/       # SCE, Monaco-prediction, and pseudobulk caches
│       ├── 03_public_skin_rna/           # ncells.RData, SubCluster_Labels.RData, Public Data Limma.RData
│       ├── 04_integrative_modeling/      # T1_T3_Multiomics_Community_Results.xlsx and cached model/figure inputs
│       └── flow_gating/                  # bcell, tcell, monocyte, dcnk.RData
```

### Preprocessing scope

Figure S12 uses `publicDataClean.RData`, `xCell_Labels.csv`, and `GSE169440_subpopulation_abundance.csv` under `data/processed/03_public_skin_rna/`. The compact abundance table contains the figure-ready public-skin cell counts without the much larger single-cell object.

The PBMC preprocessing script, `Lyme.Rmd`, documents the 10x-to-Seurat workflow and downstream analyses. Zenodo hosts intermediate and processed objects generated by this preprocessing script and consumed by the downstream figure scripts. It is not a one-command rebuild of every cache: in particular, the  `Lyme.Rmd` does not currently serialize the final annotated Seurat object, so that object is supplied directly in the Zenodo data release. 

The skin preprocessing script generates `combined3.RData` and `ncells.RData`. `SubCluster_Labels.RData` and `Public Data Limma.RData` are supplied downstream objects from the original analyses.

## Quick start

To reproduce the analysis, the code from GitHub must be combined with the data from Zenodo.

1. Clone the repository
```
git clone https://github.com/subramanian-group/Lyme.git
cd Lyme
```
2. Download `Proteomics_Metabolomics_scRNAseq_Raw.zip` and `Proteomics_Metabolomics_scRNAseq_Processed.zip` from the [latest Zenodo data release](https://doi.org/10.5281/zenodo.15395686). Download `Flow_Data_and_Preprocessing.zip` only when rendering Figure 6d, e and Supplementary Figure 7, which require the raw flow cytometry data and its standalone preprocessing workflow. 

3. Extract `Proteomics_Metabolomics_scRNAseq_Raw.zip` and `Proteomics_Metabolomics_scRNAseq_Processed.zip` directly into `data/`. After extraction, the repository must contain `data/raw/`, `data/processed/`, `data/intermediate/`, and `data/metadata/`; an extra wrapper directory indicates an incorrect extraction. Keep the flow cytometry archive separate. The four small derived gating objects used by the figure scripts, except  `fig_6.qmd` and `fig_s7.qmd`,   are already supplied under `data/intermediate/flow_gating/`.

   `fig_6.qmd` and `fig_s7.qmd` also use assets from the standalone flow cytometry archive. When rendering Figure 6d, e and Supplementary Figure 7, place the files extracted from `Flow_Data_and_Preprocessing.zip` into a `Flow/` directory under `data/`.

4. To render a figure, open the corresponding `.qmd` file (for example, `fig_3.qmd`) in RStudio and click `Render`. Paths used by the figure scripts in the latest release are relative to the repository root.

## Software requirements

- R (>= 4.1.0).
- Key R packages: `Seurat`, `DESeq2`, `limma`, `glmmLasso`, `igraph`, `tidyverse`.
- Detailed session information is included at the end of each script.

## Citation

Rostomily C et al. Multiomics reveals compartmentalized immune responses and tissue-vascular signatures in Lyme disease. bioRxiv. 2025. https://doi.org/10.1101/2025.05.27.656061
