
# This script runs the doublet detection per batch after the basic preprocessing

##########
### Load parameters and packages
##########

message(" Load parameters and packages ")

library(magrittr)
library(scUtils)
library(scUtils)

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
#parameter_list = readRDS(file=paramfile_name)
parameter_list = jsonlite::read_json(param_file)

# read features to exlude
features_exclude_list= unlist(jsonlite::read_json(parameter_list$genes_to_exclude_file))

# get filepath to processed_file
split_by_slash = strsplit(parameter_list$raw_file,split = "/")[[1]]
target_dir = paste0(split_by_slash[1:min(1,(length(split_by_slash)-1))],collapse = "/")
processed_file = paste0(parameter_list$data_path,target_dir,"/",parameter_list$dataset_name,parameter_list$processed_suffix,".rds")

# load seurat
seurat_processed = readRDS(processed_file)
seurat_processed_meta = seurat_processed@meta.data


##########
### Doublet detection
##########
message(" Doublet detection ")

seurat_object_batch_list <- Seurat::SplitObject(seurat_processed, split.by = "Batch_ID")

for(i in 1:length(seurat_object_batch_list)){

  current_seurat = seurat_object_batch_list[[i]]
  # remove features that should not be HVGs
  keep_genes = rownames(current_seurat)[!rownames(current_seurat) %in% features_exclude_list]
  current_seurat_new = CreateSeuratObject(current_seurat@assays$RNA@counts[keep_genes,],meta.data = current_seurat@meta.data)

  # run preprocessing for the current batch
  seurat_processed = scUtils::seurat_recipe(current_seurat_new,
                                            nfeatures_vst = parameter_list$nfeatures_vst_prelim,
                                            sample_column = parameter_list$sample_column,
                                            normalize_data = TRUE,
                                            remove_hvgs = TRUE,
                                            genes_to_remove = features_exclude_list,
                                            calcUMAP = FALSE,
                                            findClusters = TRUE,
                                            npcs_PCA = parameter_list$npcs_PCA,
                                            clusterRes = 1,
                                            k.param = parameter_list$k_param
                                            seed = parameter_list$global_seed)

  # run doublet detection
  current_seurat_new = scUtils::apply_DoubletFinder( # scUtils::
    seurat_object = current_seurat_new,
    npcs_PCA = parameter_list$npcs_PCA,
    pN_fixed = parameter_list$pN_fixed,
    pK_max = parameter_list$pK_max,
    doublet_formation_rate = parameter_list$doublet_formation_rate,
    adjust_nExp = FALSE,
    doublet_cluster_tresh = parameter_list$doublet_cluster_tresh,
    cluster_column = "preliminary_clusters",
    return_seurat = TRUE
  )
  #current_seurat_new$Doublet = doublet_anno
  seurat_object_batch_list[[i]] = current_seurat_new
}

##########
### add doublet information back into processed seurat
##########

message(" Add Doublet annotation ")


##########
### save as processed file
##########

message(" Save processed file ")

saveRDS(seurat_processed,file = processed_file)

message(" Complete ")


