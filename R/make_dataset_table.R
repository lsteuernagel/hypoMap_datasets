## This script creates the table that will be used to pre-process datasets
# requires some manualy decisions which are added here.

##load
library(magrittr)
#hypothalamus_neurons_map = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/hypothalamus_neurons_map.rds")
# count clusters
per_dataset_clusters = hypothalamus_neurons_map@meta.data %>% dplyr::group_by(Dataset,K98_pruned) %>% dplyr::count() %>% dplyr::filter(n > 20) %>% dplyr::group_by(Dataset) %>% dplyr::count(name = "clusters_per_dataset")
per_dataset_clusters$estimated_total_clusters = per_dataset_clusters$clusters_per_dataset + 15
per_dataset_clusters$Dataset = paste0(per_dataset_clusters$Dataset,"10x")
per_dataset_clusters$Dataset[per_dataset_clusters$Dataset=="Campbell10x"] = "CampbellDropseq"
per_dataset_clusters$Dataset[per_dataset_clusters$Dataset=="Chen10x"] = "ChenDropseq"
per_dataset_clusters$Dataset[per_dataset_clusters$Dataset=="Rossi10x"] = "RossiDropseq"
per_dataset_clusters$Dataset[per_dataset_clusters$Dataset=="Lee_Idol10x"] = "LeeDropseq"
per_dataset_clusters$Dataset[per_dataset_clusters$Dataset=="wenDropSeq10x"] = "wenDropseq"
per_dataset_clusters$Dataset[per_dataset_clusters$Dataset=="wen10x10x"] = "wen10x"
#per_dataset_clusters$Dataset[per_dataset_clusters$Dataset=="kimDev10x"] = "KimDev10x"
# add new ones
per_dataset_clusters_new = data.frame(Dataset=c("Dowsett10xnuc","Affinati10x","Anderson10x","Rupp10x","Morris10x"),estimated_total_clusters=c(90,90,70,90,50))
per_dataset_clusters=dplyr::bind_rows(per_dataset_clusters,per_dataset_clusters_new)

## load sra tables
hypomap_data_dir = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/"
all_seurat_files = list.files(hypomap_data_dir,"seurat_raw.rds",recursive = TRUE)

data_set_table = data.frame(seurat_file = all_seurat_files)
data_set_table$Dataset = stringr::str_extract(data_set_table$seurat_file,pattern = "/[a-zA-Z\\_0-9]+_raw.rds")
data_set_table$Dataset = data_set_table$Dataset %>% stringr::str_remove(pattern = "/")  %>% stringr::str_remove(pattern = "_seurat_raw.rds")
data_set_table = dplyr::left_join(data_set_table,per_dataset_clusters %>% dplyr::select(Dataset,estimated_total_clusters),by=c("Dataset"="Dataset"))

# manually set author exclude
data_set_table$exclude_author = FALSE
data_set_table$exclude_author[data_set_table$Dataset == "Dowsett10xnuc"] = TRUE
data_set_table$exclude_author[data_set_table$Dataset == "Kim10x"] = TRUE
data_set_table$exclude_author[data_set_table$Dataset == "Moffit10x"] = TRUE

# manually set target doublet rate:
data_set_table$doublet_formation_rate = 0.05
data_set_table$doublet_formation_rate[data_set_table$Dataset == "Dowsett10xnuc"] = 0.02
data_set_table$doublet_formation_rate[data_set_table$Dataset == "Kim10x"] = 0.02
#also some of the maller ones:
data_set_table$doublet_formation_rate[data_set_table$Dataset == "kimDev10x"] = 0.02
data_set_table$doublet_formation_rate[data_set_table$Dataset == "LeeDropseq"] = 0.02
#data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/")

# save
data.table::fwrite(data_set_table,"data/dataset_overview.tsv",sep="\t")


