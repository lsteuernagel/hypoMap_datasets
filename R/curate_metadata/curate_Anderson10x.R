
##########
### Info
##########

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193921

# does not look like there is any available metdata

# 10x Chromium 3
#  	strain: C57BL/6N
# ventrolateral subdivision of the ventromedial hypothalamus (VMHvl)

##########
### Load data
##########

library(dplyr)
library(Seurat)

## read raw file and extract metadata
Anderson10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Anderson10x/"
Anderson10x_seurat_raw = readRDS(paste0(Anderson10x_raw_data_path,"Anderson10x_seurat_raw.rds"))
Anderson10x_seurat_raw_meta = Anderson10x_seurat_raw@meta.data

##########
### Curate metadata
##########

Anderson10x_seurat_meta = Anderson10x_seurat_raw_meta
Anderson10x_seurat_meta$Cell_ID

### set other columns
Anderson10x_seurat_meta$Dataset = "Anderson10x"
Anderson10x_seurat_meta$Technology ="10xv3"
Anderson10x_seurat_meta$Diet ="Normal chow"
Anderson10x_seurat_meta$Pooled = NA
Anderson10x_seurat_meta$Age = NA
Anderson10x_seurat_meta$Author_CellType = NA
Anderson10x_seurat_meta$Author_Class = NA
# ventrolateral subdivision of the ventromedial hypothalamus (VMHvl)
Anderson10x_seurat_meta$Author_Region = "Ventromedial hypothalamus (ventrolateral)"

# infer sex
Anderson10x_seurat_meta$inferred_sex = scUtils::infer_sex(Anderson10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)
Anderson10x_seurat_meta$Sex = "F"

# curate some sample ids
Anderson10x_seurat_meta$Sample_ID = paste0("Anderson10x_",Anderson10x_seurat_meta$Run_ID)

# mark cells that authors excluded (or did not annotate)
Anderson10x_seurat_meta$Author_Exclude = "no"
##
Anderson10x_seurat_meta$Author_Condition = paste0(Anderson10x_seurat_meta$sex,"_",Anderson10x_seurat_meta$Developmental_stage,"_",Anderson10x_seurat_meta$reproductive_state,"_",Anderson10x_seurat_meta$social_behavior)

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Anderson10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(Anderson10x_raw_data_path,"Anderson10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Anderson10x_seurat_meta_final = Anderson10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`, Technology,Author_Region,Strain = strain,Diet,Pooled,Age,Author_Condition,inferred_sex,Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Exclude,Author_Class,Author_CellType  ) %>%
  as.data.frame()
rownames(Anderson10x_seurat_meta_final) = Anderson10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(Anderson10x_seurat_meta_final) == ncol(Anderson10x_seurat_raw)){
  Anderson10x_seurat_raw@meta.data = Anderson10x_seurat_meta_final
  saveRDS(Anderson10x_seurat_raw,paste0(Anderson10x_raw_data_path,"Anderson10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}






