if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required to resolve paths from the repository root.")
}

skin_raw_dir <- here::here("data", "raw", "03_public_skin_rna", "GSE169440")
skin_processed_dir <- here::here("data", "processed", "03_public_skin_rna")
skin_intermediate_dir <- here::here("data", "intermediate", "03_public_skin_rna")
skin_subset_dir <- file.path(skin_intermediate_dir, "unlabeled_subsets")
combined_output <- file.path(skin_processed_dir, "combined3.RData")
cell_count_output <- file.path(skin_intermediate_dir, "ncells.RData")
force_rebuild <- tolower(Sys.getenv("LYME_REBUILD_PREPROCESSING", "false")) %in%
  c("1", "true", "yes")

cached_outputs <- c(combined_output, cell_count_output)

if (all(file.exists(cached_outputs)) && !force_rebuild) {
  message(
    "Using deposited caches; raw GSE169440 integration was not rerun:\n- ",
    paste(cached_outputs, collapse = "\n- ")
  )
} else {
  if (!dir.exists(skin_raw_dir)) {
    stop(
      "Cannot rebuild the GSE169440 integration because the raw 10x directory is missing: ",
      skin_raw_dir,
      "\nUse the deposited combined3.RData and ncells.RData caches or restore the raw inputs."
    )
  }

  required_packages <- c(
    "Seurat", "GEOquery", "ggplot2", "cowplot", "patchwork", "multtest",
    "metap", "dplyr", "scCATCH", "SingleR", "celldex"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
  }

  invisible(lapply(required_packages, library, character.only = TRUE))
  dir.create(skin_processed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(skin_intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(skin_subset_dir, recursive = TRUE, showWarnings = FALSE)

geo = getGEO("GSE169440", destdir = skin_raw_dir)

dirs = list.dirs(skin_raw_dir, full.names = T, recursive = T)
dirs = grep("HHT.*filter",dirs,value=T)
dirs = grep("SKL|SKN",dirs,value=T)
names = sub("_GEX_HHT_cellranger$", "", basename(dirname(dirs)))
data = mapply(function(x,name){
  d = Read10X(
    data.dir=x,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
  )
  d = CreateSeuratObject(counts = d, project = name)
  d[["percent.mt"]] = PercentageFeatureSet(d, pattern = "^MT-")
  d = subset(d, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 15)
  d = NormalizeData(d,normalization.method = "LogNormalize",scale.factor = 10000)
  d = FindVariableFeatures(d, selection.method = "vst", nfeatures = 2000)
},dirs,names)
names(data) = names
features = SelectIntegrationFeatures(object.list = data)
anchors = FindIntegrationAnchors(object.list = data, anchor.features = features)
combined = IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"
combined = ScaleData(combined, verbose = FALSE)
combined = RunPCA(combined, npcs = 30, verbose = FALSE)

ElbowPlot(combined,ndims = 30)

combined = RunUMAP(combined, reduction = "pca", dims = 1:30)
combined = FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined = FindClusters(combined, resolution = .6)

DimPlot(combined, reduction = "umap",label = T)

md = geo$`GSE169440-GPL24676_series_matrix.txt.gz`@phenoData@data
md = md[md$title%in%combined@meta.data$orig.ident,]

metaNew = md[match(combined@meta.data$orig.ident,md$title),]
rownames(metaNew) = rownames(combined@meta.data)
metaNew$condition = gsub(".*[1-9] ","",metaNew$source_name_ch1)
metaNew$condition[metaNew$condition=="Unaffected Lesion"] = "Unaffected Skin"
metaNew$id = gsub("Subject (.*[1-9]) .*","\\1",metaNew$source_name_ch1)

combined = AddMetaData(object = combined,metadata = metaNew)

DimPlot(combined, reduction = "umap", group.by = "orig.ident",label = T)
DimPlot(combined, reduction = "umap",label = T,)

DefaultAssay(combined) <- "RNA"
# markers = lapply(as.numeric(levels(combined$seurat_clusters)),function(x)FindConservedMarkers(combined, ident.1 = x, grouping.var = "condition", verbose = FALSE))
# names(markers) = levels(combined$seurat_clusters)
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_small = markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

### scCATCH
scCATCH_markers = findmarkergenes(combined,
                                  species = 'Human',
                                  match_CellMatch = T,
                                  tissue = c('Skin','Dermis','Epithelium',
                                             'Blood','Peripheral blood',
                                             'Blood vessel','Artery','Adventitia','Antecubital vein',
                                             'Hair follicle','Muscle','Skeletal muscle',
                                             'Abdominal adipose tissue','Adipose tissue','Brown adipose tissue',
                                             'Fat pad','Subcutaneous adipose tissue','Visceral adipose tissue',
                                             'White adipose tissue'))
scCATCH_ann <- scCATCH(object = scCATCH_markers$clu_markers,
                       species = 'Human',
                       tissue = c('Skin','Dermis','Epithelium',
                                  'Blood','Peripheral blood',
                                  'Blood vessel','Artery','Adventitia','Antecubital vein',
                                  'Hair follicle','Muscle','Skeletal muscle',
                                  'Abdominal adipose tissue','Adipose tissue','Brown adipose tissue',
                                  'Fat pad','Subcutaneous adipose tissue','Visceral adipose tissue',
                                  'White adipose tissue'))

### SingleR
hpca.se <- celldex::HumanPrimaryCellAtlasData()
SingleR_ann = SingleR(test = GetAssayData(combined[['RNA']],slot='counts'), 
                      ref = hpca.se,
                      labels = hpca.se$label.main,
                      clusters = combined$seurat_clusters)

labels = data.frame(cluster = 0:31,
                    scCATCH = scCATCH_ann$cell_type[match(0:31,scCATCH_ann$cluster)],
                    SingleR = SingleR_ann$labels)

# Optional intermediate checkpoint:
# save(
#   list = c('combined','dirs','names','markers','labels','scCATCH_ann','SingleR_ann'),
#   file = file.path(skin_intermediate_dir, "combined2.RData")
# )

FeaturePlot(object = combined, features = "MS4A1",label = T)

labels_broad = c(
  `0` = "T Cell",
  `1` = "T Cell",
  `2` = "Monocyte",
  `3` = "Keratinocyte",
  `4` = "Fibroblast",
  `5` = "Fibroblast",
  `6` = "Endothelial",
  `7` = "T Cell",
  `8` = "Monocyte",
  `9` = "NK Cells",
  `10` = "Fibroblast",
  `11` = "Keratinocyte",
  `12` = "T Cell",
  `13` = "Keratinocyte",
  `14` = "Myocyte",
  `15` = "DC",
  `16` = "Melanocyte",
  `17` = "Fibroblast",
  `18` = "Myocyte",
  `19` = "B Cell",
  `20` = "DC",
  `21` = "T Cell",
  `22` = "Keratinocyte",
  `23` = "B Cell",
  `24` = "Endothelial",
  `25` = "Keratinocyte",
  `26` = "Unknown",
  `27` = "Endothelial",
  `28` = "Fibroblast",
  `29` = "Neuronal",
  `30` = "Keratinocyte",
  `31` = "Endothelial")

markers_broad = list(
  `0` = c("CD3E"),
  `1` = c("CD3E"),
  `2` = c("RNASE1","CD14","CD68"),
  `3` = c("KRT5","LY6D"),
  `4` = c("COL1A1","DCN"),
  `5` = c("COL1A1","DCN"),
  `6` = c("PECAM1","VWF"),
  `7` = c("CD3E"),
  `8` = c("RNASE1","CD14","CD68"),
  `9` = c("GNLY"),
  `10` = c("COL1A1","DCN"),
  `11` = c("KRT5","LY6D"),
  `12` = c("CD3E"),
  `13` = c("KRT5","LY6D"),
  `14` = c("TAGLN"),
  `15` = c("FCER1A","CD1A"),
  `16` = c("MLANA"),
  `17` = c("COL1A1","DCN"),
  `18` = c("TAGLN"),
  `19` = c("MS4A1"),
  `20` = c("CLEC9A","CD1A"),
  `21` = c("CD3E"),
  `22` = c("KRT5","LY6D"),
  `23` = c("IGKC"),
  `24` = c("PECAM1","VWF"),
  `25` = c("KRT5","LY6D"),
  `26` = c(""),
  `27` = c("PECAM1","VWF"),
  `28` = c("COL1A1","DCN"),
  `29` = c("NRXN1"),
  `30` = c("KRT5","LY6D"),
  `31` = c("PECAM1","VWF"))

cell_class = c(
  `0` = "T Cell",
  `1` = "T Cell",
  `2` = "Myeloid",
  `3` = "Keratinocyte",
  `4` = "Fibroblast",
  `5` = "Fibroblast",
  `6` = "Endothelial",
  `7` = "T Cell",
  `8` = "Myeloid",
  `9` = "NK Cells",
  `10` = "Fibroblast",
  `11` = "Keratinocyte",
  `12` = "T Cell",
  `13` = "Keratinocyte",
  `14` = "Myocyte",
  `15` = "Myeloid",
  `16` = "Melanocyte",
  `17` = "Fibroblast",
  `18` = "Myocyte",
  `19` = "B Cell",
  `20` = "Myeloid",
  `21` = "T Cell",
  `22` = "Keratinocyte",
  `23` = "B Cell",
  `24` = "Endothelial",
  `25` = "Keratinocyte",
  `26` = "Unknown",
  `27` = "Endothelial",
  `28` = "Fibroblast",
  `29` = "Neuronal",
  `30` = "Keratinocyte",
  `31` = "Endothelial")

labels = cbind(labels,`Manual` = labels_broad,`subcluster_class` = cell_class)

combined[["broad.labels"]] = labels$Manual[match(combined$seurat_clusters, labels$cluster)]
combined[["subcluster_class"]] = labels$subcluster_class[match(combined$seurat_clusters, labels$cluster)]


save(list=c('combined','dirs','names','markers','labels','scCATCH_ann',
            'SingleR_ann','markers_broad'),
     file = combined_output)


DimPlot(combined,group.by = 'broad.labels',label = T,split.by = 'condition') + NoLegend()
DimPlot(combined,group.by = 'broad.labels',label = T) + NoLegend()

FeaturePlot(combined,'CXCL9',split.by = 'condition')
FeaturePlot(combined,'S100A12',split.by = 'condition')

n_cells = combined@meta.data%>%group_by(orig.ident,broad.labels)%>%summarise(ncells = n())
n_cells = n_cells%>%group_by(orig.ident)%>%mutate(total_cells = sum(ncells))
n_cells = as.data.frame(n_cells)

save(n_cells, file = cell_count_output)

### Save Subsets
classes = unique(combined$subcluster_class)
for(i in 1:length(classes)){
  sub = subset(combined, subset = subcluster_class == classes[i])
  subset_name <- gsub("[^A-Za-z0-9]+", "_", classes[i])
  save(
    sub,
    file = file.path(skin_subset_dir, paste0("Unlabeled_", subset_name, ".RData"))
  )
}
rm(sub)





}

print(sessionInfo())
