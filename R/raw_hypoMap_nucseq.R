
##########
### Load count data into seurat objects and save as raw
##########

library(Seurat)
library(dplyr)
library(scUtils)

## this script creates the raw seurat object for our ow in-hoise nucseq data (no sra tables available..)

###
nucseq_data_dir = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Dowsett10xnuc/"
feature_column_idx = 2
matrix_type_regex = "filtered_feature_bc"

# get relevant files
dataset_grep = "Dowsett10xnuc"
dataset_dir = list.dirs(hypomap_data_dir,recursive = FALSE)[grepl(dataset_grep,list.dirs(hypomap_data_dir,recursive = FALSE))]
dataset_grep = strsplit(dataset_dir,"/")[[1]][length(strsplit(dataset_dir,"/")[[1]])]
sample_files = list.files(path = dataset_dir,recursive = TRUE,full.names = TRUE)
message("Processing dataset ",dataset_grep)

# get samples
sample_overview = list.files(dataset_dir,recursive = FALSE,full.names = FALSE)

# init list with temp seurats
all_run_seurats =list()
# for all runs:
for(sample_run in sample_overview){
  skip_sample =FALSE
  message("  Reading ",sample_run)
  # find all files for current run
  sample_run_files = sample_files[grepl(sample_run,sample_files)]
  # load matrix
  if(any(grepl(".mtx",sample_run_files))){
    sample_run_files = sample_run_files[grepl(matrix_type_regex,sample_run_files)]
    sra_run_counts <- scUtils::Read10xFormat(mtx = sample_run_files[grepl("matrix.mtx",sample_run_files)], cells = sample_run_files[grepl("barcodes",sample_run_files)], features = sample_run_files[grepl("features",sample_run_files)], feature.column = feature_column_idx)
  }else if(any(grepl(".dge.txt",sample_run_files))){
    sra_run_counts <- scUtils::ReadDGEFormat(dge =  sample_run_files[grepl(".dge.txt",sample_run_files)], feature.column = 1)
  }else{
    skip_sample =TRUE
    message("Warning: cannot find file to load count matrix for run ",sample_run)
  }
  if(!skip_sample){
    # add run name to column names
    name_to_use = tolower(gsub(" ","",sample_run))
    colnames(sra_run_counts) = paste0(colnames(sra_run_counts),"_",name_to_use,"_","Dowsett10xnuc")

    # make seurat object
    sra_run_seurat = SeuratObject::CreateSeuratObject(sra_run_counts,project = name_to_use,min.cells = 0,  min.features = 0 )
    sra_run_seurat@meta.data$Cell_ID = colnames(sra_run_seurat)
    sra_run_seurat@meta.data$Run_ID = name_to_use
    all_run_seurats[[name_to_use]] = sra_run_seurat
  }
}
# merge seurat objects:
dataset_seurat <- merge(all_run_seurats[[1]], y = all_run_seurats[2:length(all_run_seurats)], project = dataset_grep)
# add mt
dataset_seurat[["percent.mt"]] <- PercentageFeatureSet(dataset_seurat, pattern = "^mt-")

# save seurat object
message("Saving dataset to",paste0(dataset_dir,"/",dataset_grep,"_seurat_raw.rds"))
saveRDS(dataset_seurat,file = paste0(dataset_dir,"/",dataset_grep,"_seurat_raw.rds"))





