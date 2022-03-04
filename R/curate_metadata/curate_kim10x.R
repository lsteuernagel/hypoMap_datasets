##########
### Info
##########

# data: 10.17632/ypx3sw2f7c.3
# https://www.ncbi.nlm.nih.gov/pubmed/31626771

# Special case: we cannot find the SRA but it is a dataset with good quality that covers many cell types and quite similar to other datasets (in terms of processing) so we will also include

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
Kim10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Kim10x/"
Kim10x_seurat_raw = readRDS(paste0(Kim10x_raw_data_path,"Kim10x_seurat_raw_original.rds"))
Kim10x_seurat_raw_meta = Kim10x_seurat_raw@meta.data

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)
table(hypoMap_neurons_v1_curated_metadata$Strain[hypoMap_neurons_v1_curated_metadata$Dataset=="Kim"])

##########
### Curate metadata
##########

# set relevant columns
Kim10x_seurat_raw_meta$Dataset = "Kim10x"
Kim10x_seurat_raw_meta$Sample_ID = paste0(Kim10x_seurat_raw_meta$Dataset,"_",Kim10x_seurat_raw_meta$orig.ident.x)
Kim10x_seurat_raw_meta$Author_Region = "Ventromedial hypothalamus"
Kim10x_seurat_raw_meta$Technology = "10xv2"
Kim10x_seurat_raw_meta$Sex = Kim10x_seurat_raw_meta$sex_label
Kim10x_seurat_raw_meta$Age = "6+ weeks"
Kim10x_seurat_raw_meta$Diet = "Normal chow"
Kim10x_seurat_raw_meta$Strain = "C57Bl/6N"

# author cell type
Kim10x_seurat_raw_meta$Author_CellType = Kim10x_seurat_raw_meta$cell_cluster_label
# author class
Kim10x_seurat_raw_meta$Author_Class=NA
Kim10x_seurat_raw_meta$Author_Class[!is.na(Kim10x_seurat_raw_meta$Author_CellType)]="Neurons"
#Kim10x_seurat_raw_meta$Author_Class[grepl("Ependy",Kim10x_seurat_raw_meta$Author_CellType)]="Ependymal"
Kim10x_seurat_raw_meta$Author_Class[grepl("Macro|Micro",Kim10x_seurat_raw_meta$Author_CellType)]="Immune"
#Kim10x_seurat_raw_meta$Author_Class[grepl("Tany",Kim10x_seurat_raw_meta$Author_CellType)]="Tanycyte"
Kim10x_seurat_raw_meta$Author_Class[grepl("Oligo",Kim10x_seurat_raw_meta$Author_CellType)]="Oligodendrocytes"
Kim10x_seurat_raw_meta$Author_Class[grepl("OPC",Kim10x_seurat_raw_meta$Author_CellType)]="NG/OPC"
Kim10x_seurat_raw_meta$Author_Class[grepl("Astro",Kim10x_seurat_raw_meta$Author_CellType)]="Astrocytes"
Kim10x_seurat_raw_meta$Author_Class[grepl("Mural",Kim10x_seurat_raw_meta$Author_CellType)]="Mural"
Kim10x_seurat_raw_meta$Author_Class[grepl("Endo",Kim10x_seurat_raw_meta$Author_CellType)]="Vascular"
Kim10x_seurat_raw_meta$Author_Class[grepl("Others",Kim10x_seurat_raw_meta$Author_CellType)]="NA"
table(Kim10x_seurat_raw_meta$Author_Class)

# add additional neuro labels
Kim10x_seurat_raw_meta$Author_CellType[Kim10x_seurat_raw_meta$Author_Class == "Neurons" & !is.na(Kim10x_seurat_raw_meta$neuron_cluster_label)]=
  paste0(Kim10x_seurat_raw_meta$neuron_cluster_label[Kim10x_seurat_raw_meta$Author_Class == "Neurons" & !is.na(Kim10x_seurat_raw_meta$neuron_cluster_label)])

# set "NA" to actual NA
Kim10x_seurat_raw_meta$Author_Class[Kim10x_seurat_raw_meta$Author_Class=="NA"]=NA

# author condition
Kim10x_seurat_raw_meta$Author_Condition = paste0(Kim10x_seurat_raw_meta$behavior_label,"_",Kim10x_seurat_raw_meta$sex_exp_label,"_",Kim10x_seurat_raw_meta$housing_label)
table(Kim10x_seurat_raw_meta$Author_Condition)

# author exclusion
Kim10x_seurat_raw_meta$Author_exclude = "no"
Kim10x_seurat_raw_meta$Author_exclude[grepl("Doublet",Kim10x_seurat_raw_meta$Author_CellType) | is.na(Kim10x_seurat_raw_meta$Author_CellType)] = "yes"







# use only the ones not drom reisdent intruder assay :
# Control_Exp_Single, Control_F_Naive_Group, Plain_F_Naive_Group , Plain_Naive_Group


