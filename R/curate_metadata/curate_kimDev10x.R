##########
### Info
##########

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132355
# https://www.nature.com/articles/s41467-020-18231-z

# no author clusters or other metadata found

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
KimDev10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/kimDev10x/"
KimDev10x_seurat_raw = readRDS(paste0(KimDev10x_raw_data_path,"kimDev10x_seurat_raw.rds"))
KimDev10x_seurat_raw_meta = KimDev10x_seurat_raw@meta.data

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)

##########
### Curate metadata
##########

KimDev10x_seurat_meta = KimDev10x_seurat_raw_meta

KimDev10x_seurat_meta$Dataset = "KimDev10x"
KimDev10x_seurat_meta$Tissue = "Hypothalamus"
KimDev10x_seurat_meta$Technology ="10x v2"
KimDev10x_seurat_meta$Author_Region = "Hypothalamus"
KimDev10x_seurat_meta$Diet ="Normal chow"
KimDev10x_seurat_meta$Sex = "M"
KimDev10x_seurat_meta$Age = "6+ weeks"
KimDev10x_seurat_meta$Pooled = "yes"
KimDev10x_seurat_meta$Author_Class = NA
KimDev10x_seurat_meta$Author_CellType = NA

# infer sex
KimDev10x_seurat_meta$inferred_sex = scUtils::infer_sex(KimDev10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

# curate some sample ids
KimDev10x_seurat_meta$Sample_ID = paste0("KimDev10x_",KimDev10x_seurat_meta$Run_ID)

rownames(KimDev10x_seurat_meta) = KimDev10x_seurat_meta$Cell_ID


##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = KimDev10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Cell_ID,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(KimDev10x_raw_data_path,"KimDev10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
KimDev10x_seurat_meta_final = KimDev10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(KimDev10x_seurat_meta_final) = KimDev10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(KimDev10x_seurat_meta_final) == ncol(KimDev10x_seurat_raw)){
  KimDev10x_seurat_raw@meta.data = KimDev10x_seurat_meta_final
  saveRDS(KimDev10x_seurat_raw,paste0(KimDev10x_raw_data_path,"kimDev10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}



