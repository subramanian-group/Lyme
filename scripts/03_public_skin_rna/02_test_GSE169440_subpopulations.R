if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required to resolve paths from the repository root.")
}

skin_processed_dir <- here::here("data", "processed", "03_public_skin_rna")
skin_intermediate_dir <- here::here("data", "intermediate", "03_public_skin_rna")
skin_plot_dir <- file.path(skin_intermediate_dir, "differential_testing_plots")

labeled_files <- file.path(
  skin_intermediate_dir,
  c(
    "Labeled_Endothelial.RData", "Labeled_Fibroblast.RData",
    "Labeled_Keratinocyte.RData", "Labeled_Myeloid.RData",
    "Labeled_Myocytes.RData", "Labeled_TCell.RData"
  )
)
required_inputs <- c(
  file.path(skin_intermediate_dir, "ncells.RData"),
  file.path(skin_processed_dir, "combined3.RData"),
  labeled_files
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop(
    "Missing required deposited input(s):\n- ",
    paste(missing_inputs, collapse = "\n- ")
  )
}

subpopulation_output <- file.path(skin_intermediate_dir, "SubPop_DA_and_DEGs.RData")
top_population_output <- file.path(skin_intermediate_dir, "TopPop_DA_and_DEGs.RData")
analysis_outputs <- c(subpopulation_output, top_population_output)
force_rebuild <- tolower(Sys.getenv("LYME_REBUILD_PREPROCESSING", "false")) %in%
  c("1", "true", "yes")
smoke_test <- tolower(Sys.getenv("LYME_SMOKE_TEST", "false")) %in%
  c("1", "true", "yes")

if (smoke_test) {
  message(
    "Smoke test passed: all deposited inputs for skin differential testing are present. ",
    "The DESeq2 models were not refit."
  )
} else if (all(file.exists(analysis_outputs)) && !force_rebuild) {
  load(subpopulation_output)
  load(top_population_output)
  message(
    "Using locally generated differential-testing caches:\n- ",
    paste(analysis_outputs, collapse = "\n- ")
  )
} else {
  required_packages <- c(
    "DESeq2", "Seurat", "GEOquery", "ggplot2", "cowplot", "patchwork",
    "multtest", "metap", "dplyr", "scCATCH", "Matrix.utils", "magrittr",
    "purrr", "stringr", "tibble", "ggsignif"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
  }

  invisible(lapply(required_packages, library, character.only = TRUE))
  dir.create(skin_intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(skin_plot_dir, recursive = TRUE, showWarnings = FALSE)

# Functions ---------------------------------------------------------------

# DESeq Wrapper

DESeq_wrapper = function(sub){
  # Subset metadata to only include the cluster and sample IDs to aggregate across
  groups <- sub@meta.data[, c("Identity", "orig.ident")]
  # Aggregate across cluster-sample groups
  pb <- aggregate.Matrix(t(sub@assays$RNA@counts),
                         groupings = groups, fun = "sum")
  splitf <- sapply(stringr::str_split(rownames(pb),
                                      pattern = "_",
                                      n = 2),`[`, 1)
  pb <- split.data.frame(pb, factor(splitf)) %>% lapply(function(u)
    set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+")))
  get_sample_ids <- function(x){
    pb[[x]] %>%
      colnames()
  }
  de_samples <- map(1:length(unique(sub$Identity)), get_sample_ids) %>%
    unlist()
  samples_list <- map(1:length(unique(sub$Identity)), get_sample_ids)
  
  get_cluster_ids <- function(x){
    rep(names(pb)[x],
        each = length(samples_list[[x]]))
  }
  de_cluster_ids <- map(1:length(unique(sub$Identity)), get_cluster_ids) %>%
    unlist()
  gg_df <- data.frame(cluster_id = de_cluster_ids,
                      sample_id = de_samples)
  ## Turn named vector into a numeric vector of number of cells per sample
  n_cells <- table(sub$orig.ident)
  sids = names(n_cells)
  n_cells = as.numeric(n_cells)
  ## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
  m <- match(sids, sub$orig.ident)
  ## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
  ei <- data.frame(sub@meta.data[m, ],
                   n_cells, row.names = NULL)
  gg_df <- cbind(gg_df, ei[match(gg_df$sample_id,ei$orig.ident), ])
  metadata <- gg_df %>% dplyr::select(cluster_id, sample_id, condition, id)
  metadata$cluster_id = as.factor(metadata$cluster_id)
  clusters <- levels(metadata$cluster_id)
  de_samples <- map(1:length(unique(sub$Identity)), get_sample_ids) %>%
    unlist()
  
  resList = list()
  ddsList = list()
  for(i in 1:length(clusters)){
    message(paste0(i,': DE Testing on population: ',clusters[i]))
    if(min(table(metadata$condition[metadata$cluster_id==clusters[i]]))<2){
      message(i,': DEGs not computed for ',clusters[i],'. Cluster does not have enough cells.')
      ddsList[[i]] = NULL
      resList[[i]] = NULL
    }else if(length(table(metadata$condition[metadata$cluster_id==clusters[i]]))!=2){
      message(i,': DEGs not computed for ',clusters[i],'. Cluster only exists for one group.')
      ddsList[[i]] = NULL
      resList[[i]] = NULL
    }else{
      # Subset the metadata
      cluster_metadata <- metadata[which(metadata$cluster_id == clusters[i]), ]
      cluster_metadata$condition = factor(cluster_metadata$condition,levels=c('Unaffected Skin','EM Lesion'))
      # Assign the rownames of the metadata to be the sample IDs
      rownames(cluster_metadata) <- cluster_metadata$sample_id
      # Subset the counts
      counts <- pb[[clusters[i]]]
      cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
      dds <- DESeqDataSetFromMatrix(cluster_counts,
                                    colData = cluster_metadata,
                                    design = ~ id + condition) ##NOTE the id + condition this is the big change
      message('conducting degs for cluster ',clusters[i])
      dds <- DESeq(dds,quiet = T)
      ddsList[[i]] = dds
      # resultsNames(dds) # lists the coefficients
      resList[[i]] <- results(dds, name="condition_EM.Lesion_vs_Unaffected.Skin")
    }
  }
  names(resList) = clusters
  names(ddsList) = clusters
  return(list(res = resList,dds=ddsList))
}


# SubPopulations ------------------------------------------------------------------

load(file.path(skin_intermediate_dir, "ncells.RData"))
files = labeled_files
pops = gsub(".*led_(.*)\\.RD.*","\\1",files)

resList = list()
daList = list()
for(i in 1:length(files)){
  print(paste0('Working on SubPopulation: ',pops[i]))
  
  # Differential Expression
  print('Conducting Differential Expression')
  load(files[i])
  DefaultAssay(sub) = 'RNA'
  resList[[i]] = DESeq_wrapper(sub)
  
  #Sub-population Summary Plots
  DefaultAssay(sub) = 'integrated'
  print('Plotting Markers')
  
  png(file.path(skin_plot_dir, paste0("SubPops_", pops[i], ".png")), width = 1000, height = 800)
  dim_plot = DimPlot(sub,group.by = 'Identity',label = T,pt.size = 2,label.size = 7) + NoLegend() + ggtitle(pops[i])
  print(dim_plot)
  dev.off()
  
  png(file.path(skin_plot_dir, paste0("DotPlot_", pops[i], ".png")), width = 1000*2, height = 1200*2)
  markers_small = markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)
  feat = Map(c,markers_used,split(markers_small$gene,markers_small$cluster))
  feat = lapply(feat,unique)
  feat = unique(unlist(feat))
  dp = DotPlot(sub,features = feat,group.by = 'seurat_clusters',cluster.idents = T,dot.scale = 30,)
  dp$data$cluster = labels[as.numeric(as.character(dp$data$id))+1]
  dp = dp + 
    facet_grid(rows = vars(cluster), 
               scales = "free_y", space = "free_y", switch = "y") +
    theme_minimal(base_size = 35) + 
    theme(axis.text.x = element_text(angle=45))
  print(dp)
  dev.off()
  
  # Subpopulation differential abundance
  print('Calculating Differential Abundance')
  meta = sub@meta.data %>% group_by(orig.ident,Identity) %>% mutate(count = n())
  meta = meta %>% group_by(orig.ident) %>% mutate(total_cells_sub = n())
  meta_sub = meta[!duplicated(paste0(meta$orig.ident,meta$Identity)),c('orig.ident','id','Identity','condition','count','total_cells_sub')]
  meta_sub$total_cells_all = n_cells$total_cells[match(meta_sub$orig.ident,n_cells$orig.ident)]
  meta_sub$pop = meta_sub$count/meta_sub$total_cells_sub
  meta_sub$pot = meta_sub$count/meta_sub$total_cells_all
  meta_split = split(meta_sub,meta_sub$Identity)
  wilcox_parent = list()
  wilcox_total = list()
  for(j in 1:length(meta_split)){
    ids = names(table(meta_split[[j]]$id)==2)[table(meta_split[[j]]$id)==2]
    meta_split[[j]] = meta_split[[j]][meta_split[[j]]$id%in%ids,]
    meta_split[[j]] = meta_split[[j]][order(meta_split[[j]]$id),]
    if(nrow(meta_split[[j]])>1){
      meta_split_split = split(meta_split[[j]],meta_split[[j]]$condition)
      x = meta_split_split$`EM Lesion`$pop
      y = meta_split_split$`Unaffected Skin`$pop
      wilcox_parent[[j]] = data.frame(p = wilcox.test(x,y,
                                                      alternative = 'two.sided',
                                                      paired = T)$p.value,
                                      `Fold Change` = mean(x/y))
      x = meta_split_split$`EM Lesion`$pot
      y = meta_split_split$`Unaffected Skin`$pot
      wilcox_total[[j]] = data.frame(p = wilcox.test(x,y,
                                                     alternative = 'two.sided',
                                                     paired = T)$p.value,
                                     `Fold Change` = mean(x/y))
    }else{
      wilcox_parent[[j]] = NA
      wilcox_total[[j]] = NA
    }
  }
  names(wilcox_parent) = names(meta_split)
  names(wilcox_total) = names(meta_split)
  meta_sub = do.call(rbind,meta_split)
  daList[[i]] = list(parent = wilcox_parent, total = wilcox_total, plotting_data = meta_sub)
  
  png(file.path(skin_plot_dir, paste0("POP_Boxplot_", pops[i], ".png")), width = 1100*2, height = 800*2)
  p = ggplot(meta_sub,aes(x = condition,y = pop, fill = condition))+
    geom_boxplot(outlier.shape = NULL,size=2)+
    geom_point(size=5)+
    geom_line(aes(group=id))+
    facet_wrap(vars(Identity),scales = 'free')+
    ylab('Percent of Parent')+
    xlab('')+
    scale_y_continuous(expand = expansion(mult = c(.1, .1))) +
    ggsignif::geom_signif(test = 'wilcox.test',
                          comparisons = list(c('Unaffected Skin','EM Lesion')),
                          test.args = list(paired=T),textsize = 10)+
    theme_minimal(base_size = 40)
  print(p)
  dev.off()
  
  png(file.path(skin_plot_dir, paste0("POT_Boxplot_", pops[i], ".png")), width = 1000*2, height = 800*2)
  p = ggplot(meta_sub,aes(x = condition,y = pot,fill = condition))+
    geom_boxplot(outlier.shape = NULL,size=2)+
    geom_point(size=5)+
    geom_line(aes(group=id))+
    facet_wrap(vars(Identity),scales = 'free')+
    ylab('Percent of Total')+
    xlab('')+
    scale_y_continuous(expand = expansion(mult = c(.1, .1))) +
    ggsignif::geom_signif(test = 'wilcox.test',
                          comparisons = list(c('Unaffected Skin','EM Lesion')),
                          test.args = list(paired=T),textsize = 10)+
    theme_minimal(base_size = 40)
  print(p)
  dev.off()
}
names(resList) = pops
names(daList) = pops

save(list = c('resList','daList'),
     file = subpopulation_output)

# SubPop Dotplot

# load("Data/GEO/GSE169440/SubPop_DA_and_DEGs.RData")



# Top Clusters ------------------------------------------------------------

load(file.path(skin_intermediate_dir, 'ncells.RData'))

# Differential Expression
load(file.path(skin_processed_dir, 'combined3.RData'))
labs = labels$Manual
names(labs) = labels$cluster
combined <- RenameIdents(combined, labs)
combined$Identity = Idents(combined)
DefaultAssay(combined) = 'RNA'
res = DESeq_wrapper(combined)

#Sub-population Summary Plots
DefaultAssay(combined) = 'integrated'

png(file.path(skin_plot_dir, "AllPops.png"), width = 1000, height = 800)
dim_plot = DimPlot(combined,group.by = 'Identity',label = T,pt.size = 2,label.size = 7) + NoLegend() + ggtitle('All Clusters')
print(dim_plot)
dev.off()

png(file.path(skin_plot_dir, "DotPlot_AllPops.png"), width = 1000*3, height = 1200*3)
markers_small = markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
feat = Map(c,markers_broad,split(markers_small$gene,markers_small$cluster))
feat = lapply(feat,unique)
feat = unique(unlist(feat))
dp = DotPlot(combined,features = feat,group.by = 'seurat_clusters',cluster.idents = T,dot.scale = 20,)
dp$data$cluster = labels$Manual[as.numeric(as.character(dp$data$id))+1]
dp = dp + 
  facet_grid(rows = vars(cluster), 
             scales = "free_y", space = "free_y", switch = "y") +
  theme_minimal(base_size = 35) + 
  theme(axis.text.x = element_text(angle=45))
print(dp)
dev.off()

# Subpopulation differential abundance
meta = combined@meta.data %>% group_by(orig.ident,Identity) %>% mutate(count = n())
meta = meta %>% group_by(orig.ident) %>% mutate(total_cells_sub = n())
meta_sub = meta[!duplicated(paste0(meta$orig.ident,meta$Identity)),c('orig.ident','id','Identity','condition','count','total_cells_sub')]
meta_sub$total_cells_all = n_cells$total_cells[match(meta_sub$orig.ident,n_cells$orig.ident)]
meta_sub$pop = meta_sub$count/meta_sub$total_cells_sub
meta_sub$pot = meta_sub$count/meta_sub$total_cells_all
meta_split = split(meta_sub,meta_sub$Identity)
wilcox_parent = list()
wilcox_total = list()
for(j in 1:length(meta_split)){
  ids = names(table(meta_split[[j]]$id)==2)[table(meta_split[[j]]$id)==2]
  meta_split[[j]] = meta_split[[j]][meta_split[[j]]$id%in%ids,]
  meta_split[[j]] = meta_split[[j]][order(meta_split[[j]]$id),]
  if(nrow(meta_split[[j]])>1){
    meta_split_split = split(meta_split[[j]],meta_split[[j]]$condition)
    x = meta_split_split$`EM Lesion`$pop
    y = meta_split_split$`Unaffected Skin`$pop
    wilcox_parent[[j]] = data.frame(p = wilcox.test(x,y,
                                                    alternative = 'two.sided',
                                                    paired = T)$p.value,
                                    `Fold Change` = mean(x/y))
    x = meta_split_split$`EM Lesion`$pot
    y = meta_split_split$`Unaffected Skin`$pot
    wilcox_total[[j]] = data.frame(p = wilcox.test(x,y,
                                                   alternative = 'two.sided',
                                                   paired = T)$p.value,
                                   `Fold Change` = mean(x/y))
  }else{
    wilcox_parent[[j]] = NA
    wilcox_total[[j]] = NA
  }
}
names(wilcox_parent) = names(meta_split)
names(wilcox_total) = names(meta_split)
meta_sub = do.call(rbind,meta_split)
da = list(parent = wilcox_parent, total = wilcox_total, plotting_data = meta_sub)

png(file.path(skin_plot_dir, "POP_Boxplot_AllPops.png"), width = 1100*3, height = 800*3)
p = ggplot(meta_sub,aes(x = condition,y = pop, fill = condition))+
  geom_boxplot(outlier.shape = NULL,size=2)+
  geom_point(size=5)+
  geom_line(aes(group=id))+
  facet_wrap(vars(Identity),scales = 'free')+
  ylab('Percent of Parent')+
  xlab('')+
  scale_y_continuous(expand = expansion(mult = c(.1, .1))) +
  ggsignif::geom_signif(test = 'wilcox.test',
                        comparisons = list(c('Unaffected Skin','EM Lesion')),
                        test.args = list(paired=T),textsize = 10)+
  theme_minimal(base_size = 40)
print(p)
dev.off()

png(file.path(skin_plot_dir, "POT_Boxplot_AllPops.png"), width = 1000*3, height = 800*3)
p = ggplot(meta_sub,aes(x = condition,y = pot,fill = condition))+
  geom_boxplot(outlier.shape = NULL,size=2)+
  geom_point(size=5)+
  geom_line(aes(group=id))+
  facet_wrap(vars(Identity),scales = 'free')+
  ylab('Percent of Total')+
  xlab('')+
  scale_y_continuous(expand = expansion(mult = c(.1, .1))) +
  ggsignif::geom_signif(test = 'wilcox.test',
                        comparisons = list(c('Unaffected Skin','EM Lesion')),
                        test.args = list(paired=T),textsize = 10)+
  theme_minimal(base_size = 40)
print(p)
dev.off()

save(list = c('res','da'),
     file = top_population_output)
}

print(sessionInfo())
