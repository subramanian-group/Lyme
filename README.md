# Lyme Multiomics Study

This repository contains all the programing scripts used to generate the figures for our Lyme-disease multi-omics study.  

**Multiomics Reveals Compartmentalized Immune Responses and Tissue-Vascular Signatures in Lyme Disease**  
Clifford Rostomily, Ameek Bhalla, Henry Hampton, Lance Pflieger, Mary Brunkow, Denise Cooper, Ajay Akhade, Brett Smith, Kelsey Scherler, Pamela Troisch, Kai Wang, Chris Lausted, Gary Wormser, Leroy Hood, Noa Rappaport, and Naeha Subramanian.  

Correspondence to Naeha Subramanian (nsubrama@systemsbiology.org).  

---

## Code overview  
Each Quarto file corresponds to a single figure in the manuscript and is designed to be self-contained: running the file will generate all panels for that figure, as well as intermediate processed data outputs needed for figure construction. Figures follow the manuscript numbering.

| File | Description | Outputs |
|------|-------------|---------|
| **`fig_1.qmd`** | Cohort overview, symptom burden, and clinical correlations. | *Figure 1a–e* |
| **`fig_2.qmd`** | Plasma proteomics: fast- and slow‑resolving proteins, pathway enrichment, and protein–symptom correlations. | *Figure 2a–e* |
| **`fig_3.qmd`** | Integrated protein–metabolite community analysis and LASSO diagnostic models. | *Figure 3a–d* |
| **`fig_4.qmd`** | Metabolite levels in key pathways, and pathway enrichment. | *Figure 4a–d* |
| **`fig_5.qmd`** | PBMC flow‑cytometry and scRNA‑seq analyses showing minimal peripheral activation. | *Figure 5a–e* |
| **`fig_6.qmd`** | Outlier case with severe disseminated disease: altered immune, metabolic, and cell profiles. | *Figure 6a–e* |
| **`fig_7.qmd`** | Skin‑lesion scRNA‑seq linking tissue signals to circulating markers. | *Figure 7a–e* |
| **`fig_s1.qmd`** | Circulating immune‑protein trajectories in acute disease. | *Figure S1a–f* |
| **`fig_s2.qmd`** | Fever and EM patterns overlaid on PCA map of Olink proteins. | *Figure S2a–c* |
| **`fig_s3.qmd`** | Community sizes and protein–metabolite–symptom correlation matrix. | *Figure S3a–b* |
| **`fig_s4.qmd`** | Flow‑cytometry analysis pipeline. | *Figure S4* |
| **`fig_s5.qmd`** | Representative flow density plots of major PBMC populations. | *Figure S5* |
| **`fig_s6.qmd`** | PBMC composition changes and correlations with symptoms and circulating proteins. | *Figure S6a–c* |
| **`fig_s7.qmd`** | UMAPs of PBMC scRNA‑seq samples showing cell‑type clusters and sample variation. | *Figure S7a–b* |
| **`fig_s8.qmd`** | Differentially expressed genes in classical monocytes overlaid on PBMC UMAP. | *Figure S8* |
| **`fig_s9.qmd`** | Pronounced PBMC changes in the patient with severe disseminated disease (flow‑cytometry plots). | *Figure S9* |
| **`fig_s10.qmd`** | Immune and stromal cell abundance analyses of public skin and PBMC transcriptomes. | *Figure S10* |
| **`fig_s11.qmd`** | Skin scRNA‑seq:gross-clustering of cell populations and marker‑gene bubble plot. | *Figure S11a–b* |
| **`fig_s12.qmd`** | Skin scRNA‑seq:sub‑clustering of cell populations and marker-gene bubbble plots (myeloid, myocyte, and T‑cells) | *Figure S12a–c* |
| **`fig_s13.qmd`** | Skin scRNA‑seq:sub‑clustering of cell populations and marker-gene bubbble plots  (endothelial, fibroblast, and keratinocytes) | *Figure S13a–c* |
