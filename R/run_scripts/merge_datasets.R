
# This script merges all individual datasets to on object that can be used for downstream analysis and integration

##########
### Load parameters and packages
##########

message(" Load parameters and packages ")

library(magrittr)
library(scUtils)

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# read features to exlude
features_exclude_list= unlist(jsonlite::read_json(parameter_list$genes_to_exclude_file))

##########
### Load all datasets into list
##########

message(" Load datasets ")

# get all files
all_processed_files = list.files(parameter_list$data_path,pattern = parameter_list$processed_suffix,full.names = TRUE,recursive = TRUE)

# load all files
all_processed_seurats = lapply(all_processed_files, readRDS)

message(" Merge datasets ")

# merge all files:
merged_seurat = merge(x = all_processed_seurats[[1]],y = all_processed_seurats[2:length(all_processed_seurats)])

##########
### Normalized and run feature detection
##########

message(" Add variable features ")

# normalize data
seurat_object <- Seurat::NormalizeData(object = seurat_object,  verbose = F, assay = "RNA")

# find HVGs
seurat_object = scUtils::identify_variable_features(seurat_object,
                                          n_hvgs_sizes = parameter_list$feature_set_sizes,
                                          batch_var = parameter_list$sample_column,
                                          assay_name = "RNA",
                                          method = "vst",
                                          ignore_genes_vector = features_exclude_list,
                                          returnSeurat = TRUE,
                                          seed = parameter_list$global_seed)

##########
### save merged object
##########

message(" Save merged file ")

# make file name
merged_file_name = paste0(parameter_list$data_path,parameter_list$merged_name)

# save to rds
saveRDS(seurat_object,file = paste0(merged_file_name,".rds"))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = seurat_object,filename = paste0(merged_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(merged_file_name,".h5seurat"), dest =  paste0(merged_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)

message(" Complete ")





