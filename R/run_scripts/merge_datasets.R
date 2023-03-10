
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
merged_seurat <- Seurat::NormalizeData(object = merged_seurat,  verbose = F, assay = "RNA")

# find HVGs
merged_seurat = scUtils::identify_variable_features(merged_seurat,
                                                    n_hvgs_sizes = parameter_list$feature_set_sizes,
                                                    batch_var = parameter_list$sample_column,
                                                    assay_name = "RNA",
                                                    method = "vst",
                                                    ignore_genes_vector = features_exclude_list,
                                                    returnSeurat = TRUE,
                                                    seed = parameter_list$global_seed)
##########
### Run basic processing (without integration)
##########

message(" Add basic processing ")

merged_seurat = scUtils::seurat_recipe(merged_seurat,
                                       nfeatures_vst = parameter_list$nfeatures_vst_prelim,
                                       sample_column = "Batch_ID",
                                       normalize_data = TRUE,
                                       remove_hvgs = TRUE,
                                       genes_to_remove = features_exclude_list,
                                       calcUMAP = TRUE,
                                       findClusters = TRUE,
                                       npcs_PCA = parameter_list$npcs_PCA,
                                       clusterRes = 2,
                                       k.param = parameter_list$k_param,
                                       seed = parameter_list$global_seed)




##########
### save merged object
##########

message(" Save merged file ")

# make file name
merged_file_name = paste0(parameter_list$data_path,parameter_list$merged_name)

# save to rds
saveRDS(merged_seurat,file = paste0(merged_file_name,".rds"))

# save h5seurat
SeuratDisk::SaveH5Seurat(object = merged_seurat,filename = paste0(merged_file_name,".h5seurat"), overwrite = TRUE, verbose = TRUE)

# save to anndata
SeuratDisk::Convert( paste0(merged_file_name,".h5seurat"), dest =  paste0(merged_file_name,".h5ad"),assay="RNA",verbose=TRUE,overwrite=TRUE)

message(" Complete ")


##########
### save QC plots
##########

library(ggplot2)

qc_dir = paste0(parameter_list$qc_path,"/")
system(paste0("mkdir -p ",paste0(qc_dir)))
columns_to_plot = c("Doublet","Author_Exclude","Batch_ID","seurat_clusters","Dataset")

for(column_to_plot in columns_to_plot){
  if(column_to_plot %in% colnames(merged_seurat@meta.data)){
    if(length(unique(merged_seurat@meta.data[,column_to_plot]))>10){
      p1 = Seurat::DimPlot(object = merged_seurat,group.by = column_to_plot,label=TRUE,raster = FALSE)+NoLegend()
    }else{
      p1 = Seurat::DimPlot(object = merged_seurat,group.by = column_to_plot,raster = FALSE)
    }
    p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 2048)
    ggplot2::ggsave(filename = paste0(qc_dir,"merged","_",column_to_plot,".png"),
                    plot = p1r, "png",dpi=600,width=300,height = 300,units="mm")
    ggplot2::ggsave(filename = paste0(qc_dir,"merged","_",column_to_plot,".pdf"),
                    plot = p1r, "pdf",dpi=600,width=300,height = 300,units="mm")
  }
}

message(" Complete with plots")
