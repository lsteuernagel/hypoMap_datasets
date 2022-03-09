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
Kim10x_seurat_raw_meta$Pooled = "yes"

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
Kim10x_seurat_raw_meta$Author_Exclude = "no"
Kim10x_seurat_raw_meta$Author_Exclude[grepl("Doublet",Kim10x_seurat_raw_meta$Author_CellType) | is.na(Kim10x_seurat_raw_meta$Author_CellType)] = "yes"

Kim10x_seurat_meta = Kim10x_seurat_raw_meta

# infer sex
Kim10x_seurat_meta$inferred_sex = scUtils::infer_sex(Kim10x_seurat_raw,sample_column="orig.ident.x",min_xist_female = 0.7,max_xist_male = 0.1)

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Kim10x_seurat_meta %>% dplyr::distinct(Sample_ID,.keep_all = TRUE) %>% dplyr::select(-umis_label,-genes_label,-Author_exclude,-barcode,-sample_name,-Cell_ID,-orig.ident.x,-orig.ident.y,-cell_cluster_id,-cell_cluster_label,-cell_cluster_color,-neuron_cluster_id,
                                                                                                          -neuron_cluster_label,-neuron_cluster_color,-VMH_neuron_cluster_id,-VMH_neuron_cluster_label,-VMH_neuron_cluster_color,-VMH_neuron_cca_cluster_id,-VMH_neuron_cca_cluster_label,-VMH_neuron_cca_cluster_color,-Author_CellType,-Cell_ID,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(Kim10x_raw_data_path,"Kim10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Kim10x_seurat_meta_final = Kim10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset, Sample_ID, Technology,Strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex,Author_Condition, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Class,Author_CellType,Author_Exclude) %>%
  as.data.frame()
rownames(Kim10x_seurat_meta_final) = Kim10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(Kim10x_seurat_meta_final) == ncol(Kim10x_seurat_raw)){
  Kim10x_seurat_raw@meta.data = Kim10x_seurat_meta_final
  saveRDS(Kim10x_seurat_raw,paste0(Kim10x_raw_data_path,"Kim10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}

##########
### Save subset version
##########

# use only the ones not drom reisdent intruder assay :
# Control_Exp_Single, Control_F_Naive_Group, Plain_F_Naive_Group , Plain_Naive_Group#

# also remove low quality and doublets here to reduce size

# I manually copied the non-subsetted objects to a backup : Kim10x_seurat_raw_original.rds

##### Subset !

Kim10x_seurat_subset = subset(Kim10x_seurat_raw, subset = Author_Exclude == "no" & Author_Condition %in% c("Control_Exp_Single","Control_F_Naive_Group","Plain_F_Naive_Group","Plain_Naive_Group"))
Kim10x_seurat_subset

saveRDS(Kim10x_seurat_subset,paste0(Kim10x_raw_data_path,"Kim10x_seurat_raw.rds"))




