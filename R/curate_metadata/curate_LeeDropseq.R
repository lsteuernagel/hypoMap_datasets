##########
### Info
##########

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119960

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
LeeDropseq_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/LeeDropseq/"
LeeDropseq_seurat_raw = readRDS(paste0(LeeDropseq_raw_data_path,"LeeDropseq_seurat_raw.rds"))
LeeDropseq_seurat_raw_meta = LeeDropseq_seurat_raw@meta.data

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)
table(hypoMap_neurons_v1_curated_metadata$Strain[hypoMap_neurons_v1_curated_metadata$Dataset=="Lee_Idol"])

##########
### Curate metadata
##########

LeeDropseq_seurat_meta = LeeDropseq_seurat_raw_meta

LeeDropseq_seurat_meta$Sample_ID = paste0("LeeDropseq_",LeeDropseq_seurat_raw_meta$Run_ID)
LeeDropseq_seurat_meta$Dataset = "LeeDropseq"
LeeDropseq_seurat_meta$Technology ="Dropseq"
LeeDropseq_seurat_meta$Author_Region = "Hypothalamus"
LeeDropseq_seurat_meta$Diet ="HFHC_2weeks"
LeeDropseq_seurat_meta$Sex = "M"
LeeDropseq_seurat_meta$Age = "6+ weeks"
LeeDropseq_seurat_meta$Pooled = "yes"
LeeDropseq_seurat_meta$Strain =LeeDropseq_seurat_meta$strain
LeeDropseq_seurat_meta$Author_Class = NA
LeeDropseq_seurat_meta$Author_CellType = NA

# infer sex
LeeDropseq_seurat_meta$inferred_sex = scUtils::infer_sex(LeeDropseq_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

rownames(LeeDropseq_seurat_meta) = LeeDropseq_seurat_meta$Cell_ID

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = LeeDropseq_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Cell_ID,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(LeeDropseq_raw_data_path,"LeeDropseq_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
LeeDropseq_seurat_meta_final = LeeDropseq_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(LeeDropseq_seurat_meta_final) = LeeDropseq_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(LeeDropseq_seurat_meta_final) == ncol(LeeDropseq_seurat_raw)){
  LeeDropseq_seurat_raw@meta.data = LeeDropseq_seurat_meta_final
  saveRDS(LeeDropseq_seurat_raw,paste0(LeeDropseq_raw_data_path,"LeeDropseq_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}

