
##########
### Load count data into seurat objects and save as raw
##########

library(Seurat)
library(dplyr)
library(scUtils)

###
hypomap_data_dir = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/"
feature_column_idx = 2
matrix_type_regex = "filtered_feature_bc"
sra_table_dir = paste0(hypomap_data_dir,"SRAtables/")
ignore_cols = c("Assay Type", "AvgSpotLen","Bases","BioProject","BioSample", "Bytes","Center Name","Consent","DATASTORE filetype","DATASTORE provider","DATASTORE region", "Experiment","GEO_Accession (exp)","LibraryLayout","LibrarySelection" ,"LibrarySource","Organism","Platform")

# get all tables
all_sra_tables = list.files(sra_table_dir)

# init list with dataset seurats
# for each dataset:
for(dataset_table in all_sra_tables){

  # read sra table of current dataset
  sra_table = data.table::fread(paste0(sra_table_dir,dataset_table),data.table = F)
  # get relevant files
  dataset_grep = gsub("SraRunTable_","",gsub(".txt","",dataset_table))
  dataset_dir = list.dirs(hypomap_data_dir,recursive = FALSE)[grepl(dataset_grep,list.dirs(hypomap_data_dir,recursive = FALSE))]
  dataset_grep = strsplit(dataset_dir,"/")[[1]][length(strsplit(dataset_dir,"/")[[1]])]
  sra_files = list.files(path = dataset_dir,recursive = TRUE,full.names = TRUE)
  message("Processing dataset ",dataset_grep)
  # init list with temp seurats
  all_run_seurats =list()
  # for all runs:
  for(sra_run in sra_table$Run){
    skip_sample =FALSE
    message("  Reading ",sra_run)
    # find all files for current run
    sra_run_files = sra_files[grepl(sra_run,sra_files)]
    # load matrix
    if(any(grepl(".mtx",sra_run_files))){
      sra_run_files = sra_run_files[grepl(matrix_type_regex,sra_run_files)]
      sra_run_counts <- scUtils::Read10xFormat(mtx = sra_run_files[grepl("matrix.mtx",sra_run_files)], cells = sra_run_files[grepl("barcodes",sra_run_files)], features = sra_run_files[grepl("features",sra_run_files)], feature.column = feature_column_idx)
    }else if(any(grepl(".dge.txt",sra_run_files))){
      sra_run_counts <- scUtils::ReadDGEFormat(dge =  sra_run_files[grepl(".dge.txt",sra_run_files)], feature.column = 1)
    }else{
      skip_sample =TRUE
      message("Warning: cannot find file to load count matrix for run ",sra_run)
    }
    if(!skip_sample){
      # add run name to column names
      colnames(sra_run_counts) = paste0(colnames(sra_run_counts),"_",sra_run)

      # make seurat object
      sra_run_seurat = SeuratObject::CreateSeuratObject(sra_run_counts,project = sra_run,min.cells = 0,  min.features = 0 )
      sra_run_seurat@meta.data$Cell_ID = colnames(sra_run_seurat)
      sra_run_seurat@meta.data$Run_ID = sra_run
      all_run_seurats[[sra_run]] = sra_run_seurat
    }
  }
  # merge seurat objects:
  dataset_seurat <- merge(all_run_seurats[[1]], y = all_run_seurats[2:length(all_run_seurats)], project = dataset_grep)
  # add mt
  dataset_seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(dataset_seurat, pattern = "^mt-")
  ## add metadata from sra run table
  sra_table_toadd = sra_table[,!colnames(sra_table) %in% ignore_cols]
  meta_temp = dplyr::left_join(dataset_seurat@meta.data,sra_table_toadd,by=c("Run_ID"="Run")) %>% dplyr::select(-orig.ident) %>% as.data.frame()
  rownames(meta_temp) <- meta_temp$Cell_ID
  dataset_seurat@meta.data = meta_temp

  # save seurat object
  message("Saving dataset to",paste0(dataset_dir,"/",dataset_grep,"_seurat_raw.rds"))
  saveRDS(dataset_seurat,file = paste0(dataset_dir,"/",dataset_grep,"_seurat_raw.rds"))
}





