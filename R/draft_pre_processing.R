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

seurat_raw = readRDS(paste0(ChenDropseq_raw_data_path,seurat_raw_name))
seurat_raw_meta = ChenDropseq_seurat_raw@meta.data

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

# TODO: need to chaneg function so that it can take a list of genes for exclusion in hvg search
# TODO: HVG detection should be sample-aware ---> use identify_variable_features and add optionally to seurat recipe
# TODO: use function from scUTIls
# TODO: add a function that finds cluster using a certain threshold/res

source("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scHarmonize/processing/processing functions.R")

seurat_processed = scUtils::seurat_recipe(seurat_raw,nfeatures_vst = nfeatures_vst_prelim,
                           clean_hvg = TRUE,
                           normalize_data = TRUE,
                           calcUMAP = TRUE,
                           findClusters = TRUE,
                           npcs_PCA = npcs_PCA,
                           clusterRes = 2.5,
                           seed = global_seed)

DimPlot(seurat_processed,group.by = "Sample_ID")

##########
### Detect batches
##########

# todo: move params to param list
# TODO: make even more stringent
# TODO: move function to scUtils or to hypoMap_datasets ?
# TODO: rework function ?
# TODO: don't plot in actual script

seurat_newbatches= identifyBatches(seurat_object = seurat_processed, sample_variable = sample_column, cell_id = id_column,
                                         max_entropy = 0.9 ,trees= 5000,n_cores = n_cores,n_dim = npcs_PCA,embedding_name = "pca",plotClustering =TRUE,
                                         plotEntropy=TRUE,seed=1237)
batch_ids = paste0(dataset_name,seurat_newbatches$new_batch)
# or different cutoff:
batch_ids = cut_batches(all_entropies_unsorted = seurat_newbatches$all_entropies_unsorted,
            max_entropy=0.85,
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

# TODO: add wrapper forr doublet detection to scUtils

# ...

seurat_object_batch_list <- Seurat::SplitObject(seurat_processed, split.by = "Batch_ID")


##########
### Merge batches again
##########


