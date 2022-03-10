# make a slum executable script using functions from package scUtils above that for loops over all seurats:
# (load parameters ?)
# load data
# run basic qc (cell removal)
# Run basic pre-processing on all samples - merged
# Run batch-detection
# For each batch:
#   re-run pre-processing
#   run standard clustering
#   Detect Doublets:
#   select cluster with high doublet scores and remove
# Export final seurat objects


## example with Chen et al

##########
### load data
##########

# must be loaded from params:
raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/ChenDropseq/"
seurat_raw_name = "ChenDropseq_seurat_raw.rds"
n_cores = 30
id_column = "Cell_ID"
dataset_name = "ChenDropseq"
# other params
global_seed = 123456
minUMI = 1000
minFeatures = 500
maxUMI = Inf
maxUMI_dynamic = 20
max_mt = 10
genes_to_exclude = c()
exclude_author = FALSE
min_cells_sample = 100
sample_column = "SRA_ID"
nfeatures_vst_prelim = 1000
npcs_PCA = 60
target_cluster_number = 60

##
features_exclude_list=jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/features_exclude_list.json")
features_exclude_list=sapply(features_exclude_list,unlist)
gene_to_remove = features_exclude_list$hvgs_exclude_long_extra

seurat_raw = readRDS(paste0(raw_data_path,seurat_raw_name))
seurat_raw_meta = seurat_raw@meta.data

##########
### basic qc
##########

maxUMI_use = min(median(seurat_raw@meta.data$nCount_RNA)*maxUMI_dynamic,maxUMI)
seurat_raw = subset(seurat_raw,subset = nCount_RNA < maxUMI_use & nCount_RNA > minUMI & nFeature_RNA > minFeatures & percent_mt < max_mt)

# if specified: exclude based on author annotation
if(exclude_author & "Author_Exclude" %in% colnames(seurat_raw@meta.data)){
  if(all(!is.na(seurat_raw@meta.data$Author_Exclude))){
    seurat_raw = subset(seurat_raw,subset = Author_Exclude == "no")
  }
}

# if there are samples with very few cells: remove
keep_samples=names(table(seurat_raw@meta.data[,sample_column]))[table(seurat_raw@meta.data[,sample_column])>min_cells_sample]
seurat_raw@meta.data$tmp_id = seurat_raw@meta.data[,sample_column]
seurat_raw = subset(seurat_raw, subset = tmp_id %in% keep_samples)


##########
### basic pre-processing
##########

seurat_processed = scUtils::seurat_recipe(seurat_raw,nfeatures_vst = nfeatures_vst_prelim,sample_column = "SRA_ID",
                           normalize_data = TRUE,remove_hvgs = TRUE,genes_to_remove = gene_to_remove,
                           calcUMAP = TRUE,
                           findClusters = TRUE,
                           npcs_PCA = npcs_PCA,
                           clusterRes = 1,
                           seed = global_seed)

# cluster resolution:
seurat_processed = scUtils::determine_cluster_resolution(seurat_object = seurat_processed,
                                      target_cluster_number = target_cluster_number,
                                      resolutions = c(0.5, 0.75, 1, 1.5, 2:10),
                                      min_cells = 5,
                                      graph_name = "RNA_snn",
                                      cluster_col_name = "preliminary_clusters",
                                      return_seurat = TRUE,
                                      seed = 1234)

DimPlot(seurat_processed,group.by = "preliminary_clusters")

##########
### Detect batches
##########

# todo: move params to param list
# TODO: make even more stringent
# TODO: move function to scUtils or to hypoMap_datasets ?
# TODO: rework function ?
# TODO: don't plot in actual script

source("R/functions.R")
source("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/hypoMap_datasets/R/functions.R")

seurat_newbatches= identifyBatches(seurat_object = seurat_processed, sample_variable = sample_column, cell_id = id_column,
                                         max_entropy = 0.9 ,trees= 2000,n_cores = n_cores,n_dim = npcs_PCA,embedding_name = "pca",plotClustering =TRUE,
                                         plotEntropy=TRUE,seed=1237)
batch_ids = paste0(dataset_name,seurat_newbatches$new_batch)
# or different cutoff:
batch_ids = cut_batches(all_entropies_unsorted = seurat_newbatches$all_entropies_unsorted,
            max_entropy=0.8,
            hc=seurat_newbatches$clustering,
            seurat_object=seurat_processed,
            sample_variable=sample_column)
batch_ids =  paste0(dataset_name,batch_ids)

# add new batches
seurat_processed$Batch_ID = batch_ids
table(seurat_processed@meta.data[,sample_column],seurat_processed$Batch_ID)
DimPlot(seurat_processed,group.by = "Batch_ID")


##########
### For each detected batch: run processing & doublet detection
##########

# TODO: add wrapper fordoublet detection to scUtils

# ...

seurat_object_batch_list <- Seurat::SplitObject(seurat_processed, split.by = "Batch_ID")

for(i in 1:length(seurat_object_batch_list)){
  current_seurat = seurat_object_batch_list[[i]]
  keep_genes = rownames(current_seurat)[!rownames(current_seurat) %in% gene_to_remove]
  current_seurat_new = CreateSeuratObject(current_seurat@assays$RNA@counts[keep_genes,],meta.data = current_seurat@meta.data)
  current_seurat_new = scUtils::seurat_recipe(current_seurat_new,nfeatures_vst = nfeatures_vst_prelim,
                                            normalize_data = TRUE,remove_hvgs = TRUE,genes_to_remove = gene_to_remove,
                                            calcUMAP = TRUE,
                                            findClusters = FALSE,
                                            npcs_PCA = npcs_PCA,
                                            clusterRes = 1,
                                            seed = global_seed)
  current_seurat_new = scUtils::apply_DoubletFinder( # scUtils::
    seurat_object = current_seurat_new,
    npcs_PCA = npcs_PCA,
    pN_fixed = 0.25,
    pK_max = 0.1,
    doublet_formation_rate = 0.05,
    adjust_nExp = FALSE,
    doublet_cluster_tresh = 0.75,
    cluster_column = "preliminary_clusters",
    return_seurat = TRUE
  )
  #current_seurat_new$Doublet = doublet_anno
}

DimPlot(current_seurat_new,group.by = "Doublet")

##########
### Merge batches again
##########

##########
### TODO overview
##########

## load table with dataset specifications

## load json with all other attributes

## for each dataset start job with pr-processing script

## also start jobs for Doublet detection script with dependencies on pre-processing script

## save seurat objects as pre-pocessed



## next script: load all and merge to one object

## remove Doublets etc.

## run signature enrichment ?

##
