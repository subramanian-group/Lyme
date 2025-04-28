pacman::p_load(DESeq2,Seurat,GEOquery,ggplot2,cowplot,patchwork,
               multtest,metap,dplyr,scCATCH,Matrix.utils,magrittr,purrr)

# Functions ---------------------------------------------------------------

# DESeq Wrapper
remotes::install_github("cvarrichio/Matrix.utils") 

DESeq_wrapper = function(sub){
  # Subset metadata to only include the cluster and sample IDs to aggregate across
  groups <- sub@meta.data[, c("Identity", "orig.ident")]
  # Aggregate across cluster-sample groups
  pb <- Matrix.utils::aggregate.Matrix(t(sub@assays$RNA@counts),
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

load('data/GEO/GSE169440/ncells.RData')
files = list.files('Data/GEO/GSE169440',full.names = T,pattern = "Label")
pops = gsub(".*led_(.*)\\.RD.*","\\1",files)

files <- files |> stringr::str_subset("SubCluster_Labels", negate = TRUE)

resList = list()
daList = list()
for(i in 1:length(files)){
  print(paste0('Working on SubPopulation: ',pops[i]))
  
  # Differential Expression
  print('Conducting Differential Expression')
  load(files[i])
  DefaultAssay(sub) = 'RNA'
  # resList[[i]] = DESeq_wrapper(sub)
  
  #Sub-population Summary Plots
  # DefaultAssay(combined) = 'integrated'
  DefaultAssay(sub) = 'integrated'
  print('Plotting Markers')
  
  # png(paste0('figures_final/scRNAseq EM/AllPops.png'),width = 1000,height = 800)
  dim_plot = DimPlot(sub,group.by = 'Identity',label = T,pt.size = 2,label.size = 7) + NoLegend() + ggtitle(pops[i])
  print(dim_plot)
  # dev.off()
  
  # # png(paste0('AllPops_Split.png'),width = 1200,height = 800)
  # dim_plot = DimPlot(sub,group.by = 'Identity',label = T,pt.size = 2,label.size = 7, split.by = 'condition') + NoLegend() + ggtitle('All Clusters')
  # print(dim_plot)
  # # dev.off()
  
    # png(paste0("figures_final/scRNAseq EM/DotPlot_",pops[i],'.png'),width = 1000*2,height = 1200*2)
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
    # dev.off()
  