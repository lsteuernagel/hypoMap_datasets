
# This script runs the basic processing of raw seurat object including batch detection.

##########
### Load parameters and packages
##########
message(" Load parameters and packages ")

library(magrittr)
library(scUtils)
library(Seurat)

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# read features to exlude
features_exclude_list= unlist(jsonlite::read_json(parameter_list$genes_to_exclude_file))

# load seurat
seurat_raw = readRDS(paste0(parameter_list$data_path,parameter_list$raw_file))
seurat_raw_meta = seurat_raw@meta.data

##########
### basic qc
##########
message(" Basic QC ")

# add exclude features
#seurat_raw[["percent_exclude_features"]] <- Seurat::PercentageFeatureSet(seurat_raw,features=features_exclude_list[features_exclude_list %in% rownames(seurat_raw)])
seurat_raw[["percent_exclude_features"]] <- Matrix::colSums(seurat_raw@assays$RNA@counts[features_exclude_list[features_exclude_list %in% rownames(seurat_raw)],]) / seurat_raw@meta.data$nCount_RNA


# max umi
maxUMI_use = as.numeric(min(median(seurat_raw@meta.data$nCount_RNA)*parameter_list$maxUMI_dynamic,parameter_list$maxUMI))
# subset seurat:
seurat_raw = subset(seurat_raw,subset = nCount_RNA <= maxUMI_use & nCount_RNA >= parameter_list$minUMI & nFeature_RNA >= parameter_list$minFeatures & percent_mt <= parameter_list$max_mt)

# if specified: exclude based on author annotation
if(parameter_list$exclude_author & "Author_Exclude" %in% colnames(seurat_raw@meta.data)){
  if(all(!is.na(seurat_raw@meta.data$Author_Exclude))){
    seurat_raw = subset(seurat_raw,subset = Author_Exclude == "no")
  }
}

# if there are samples with very few cells: remove
keep_samples=names(table(seurat_raw@meta.data[,parameter_list$sample_column]))[table(seurat_raw@meta.data[,parameter_list$sample_column]) > parameter_list$min_cells_sample]
seurat_raw@meta.data$tmp_id = seurat_raw@meta.data[,parameter_list$sample_column]
seurat_raw = subset(seurat_raw, subset = tmp_id %in% keep_samples)
seurat_raw@meta.data = seurat_raw@meta.data[,!grepl("tmp_id",colnames(seurat_raw@meta.data))]

##########
### basic pre-processing
##########
message(" Basic pre-processing ")

seurat_processed = scUtils::seurat_recipe(seurat_raw,
                                          nfeatures_vst = parameter_list$nfeatures_vst_prelim,
                                          sample_column = parameter_list$sample_column,
                                          normalize_data = TRUE,
                                          remove_hvgs = TRUE,
                                          genes_to_remove = features_exclude_list,
                                          calcUMAP = FALSE,
                                          findClusters = TRUE,
                                          npcs_PCA = parameter_list$npcs_PCA,
                                          clusterRes = 1,
                                          k.param = parameter_list$k_param,
                                          seed = parameter_list$global_seed)

# find best cluster resolution:
seurat_processed = scUtils::determine_cluster_resolution(seurat_object = seurat_processed,
                                                         target_cluster_number = parameter_list$target_cluster_number,
                                                         resolutions =  parameter_list$resolutions_to_try,
                                                         min_cells = 5,
                                                         graph_name = "RNA_snn",
                                                         cluster_col_name = "preliminary_clusters",
                                                         return_seurat = TRUE,
                                                         seed = parameter_list$global_seed)
# remove cluster resolutions and other unwanted columns:
seurat_processed@meta.data = seurat_processed@meta.data[,!grepl("RNA_snn",colnames(seurat_processed@meta.data))]
seurat_processed@meta.data = seurat_processed@meta.data[,!grepl("seurat_clusters",colnames(seurat_processed@meta.data))]

## Mark clusters with very high features to exclude
percent_exclude_features_per_cluster = data.frame(preliminary_clusters = names(tapply(seurat_processed@meta.data[,"percent_exclude_features"],INDEX = seurat_processed@meta.data[,"preliminary_clusters"],FUN = median)),pct_median = tapply(seurat_processed@meta.data[,"percent_exclude_features"],INDEX = seurat_processed@meta.data[,"preliminary_clusters"],FUN = median))
clusters_with_high_pct = percent_exclude_features_per_cluster$preliminary_clusters[percent_exclude_features_per_cluster$pct_median >= parameter_list$max_pctExclude]
seurat_processed@meta.data$Process_Exclude = "no"
seurat_processed@meta.data$Process_Exclude[seurat_processed@meta.data$preliminary_clusters %in% clusters_with_high_pct] = "yes"

## check cell cycle scoring:
seurat_processed = Seurat::CellCycleScoring(seurat_processed,s.features = stringr::str_to_title(cc.genes$s.genes), g2m.features = stringr::str_to_title(cc.genes$g2m.genes))

##########
### Detect batches
##########
message(" Detect batches within Dataset ")

# source from functions:
source("R/functions.R")

# TODO: add downsampling to not run into memory problems:
if(ncol(seurat_processed) > parameter_list$downsampling_max){
  downsampled_meta = scUtils::downsample_balanced_iterative(metadata = seurat_processed@meta.data,
                                                            target_sample_size = parameter_list$downsampling_max,
                                                            predictor_var = "preliminary_clusters",
                                                            stepsize = 100,
                                                            global_seed = parameter_list$global_seed)
  seurat_processed_downsampled = subset(seurat_processed,cells = downsampled_meta[,1])
}else{
  seurat_processed_downsampled = seurat_processed
}
# run RF on downsampled data
seurat_newbatches= identifyBatches(seurat_object = seurat_processed_downsampled,
                                   sample_variable = parameter_list$sample_column,
                                   cell_id = parameter_list$id_column,
                                   max_entropy = parameter_list$max_entropy_batch_detection,
                                   trees= min(ncol(seurat_processed_downsampled),parameter_list$trees_rf),
                                   n_cores = parameter_list$n_cores,
                                   n_dim = parameter_list$npcs_PCA,
                                   embedding_name = "pca",
                                   plotClustering =FALSE,
                                   plotEntropy=FALSE,
                                   seed=parameter_list$global_seed)

# get per sample_id table to map back onto full processed object:
batch_id_to_sample_id = as.data.frame(table(seurat_processed_downsampled@meta.data[,parameter_list$sample_column],seurat_newbatches$new_batch))
colnames(batch_id_to_sample_id)[1:2] = c(parameter_list$sample_column,"Batch_ID")[1:2]
batch_id_to_sample_id$Batch_ID = paste0(parameter_list$dataset_name,"_",batch_id_to_sample_id$Batch_ID)
batch_id_to_sample_id = batch_id_to_sample_id[batch_id_to_sample_id$Freq > 0,]
batch_id_to_sample_id = batch_id_to_sample_id %>% dplyr::distinct(!!rlang::sym(parameter_list$sample_column),Batch_ID)
#batch_ids = paste0(parameter_list$dataset_name,seurat_newbatches$new_batch)

# add new batches
temp_meta = dplyr::left_join(seurat_processed@meta.data,batch_id_to_sample_id)
rownames(temp_meta) = temp_meta[,parameter_list$id_column]
seurat_processed@meta.data = temp_meta

##########
### save as processed file
##########
message(" Save processed file ")

# get folder from raw file
split_by_slash = strsplit(parameter_list$raw_file,split = "/")[[1]]
target_dir = paste0(split_by_slash[1:min(1,(length(split_by_slash)-1))],collapse = "/")

# save
output_file = paste0(parameter_list$data_path,target_dir,"/",parameter_list$dataset_name,parameter_list$processed_suffix,".rds")
saveRDS(seurat_processed,file = output_file)

message(" Complete ")

