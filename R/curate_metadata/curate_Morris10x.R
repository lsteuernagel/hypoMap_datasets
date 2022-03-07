##########
### Info
##########

# GSE167927

# no metadata found
# strain: C57BL/6J
# age: Postnatal day 10-12

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
Morris10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Morris10x/"
Morris10x_seurat_raw = readRDS(paste0(Morris10x_raw_data_path,"Morris10x_seurat_raw.rds"))
Morris10x_seurat_raw_meta = Morris10x_seurat_raw@meta.data

##########
### Curate metadata
##########

Morris10x_seurat_meta = Morris10x_seurat_raw_meta

Morris10x_seurat_meta$Dataset = "Morris10x"
Morris10x_seurat_meta$Sample_ID = paste0(Morris10x_seurat_meta$Dataset,"_",Morris10x_seurat_meta$Run_ID)
Morris10x_seurat_meta$Technology ="10xv2"
Morris10x_seurat_meta$Author_Region = "Suprachiasmatic nucleus"
Morris10x_seurat_meta$Diet ="Normal chow"
Morris10x_seurat_meta$Age = "0-3 weeks"
Morris10x_seurat_meta$Pooled = "yes"
Morris10x_seurat_meta$Author_Class = NA
Morris10x_seurat_meta$Author_CellType = NA

# infer sex
Morris10x_seurat_meta$inferred_sex = scUtils::infer_sex(Morris10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.8,max_xist_male = 0.1)

# curate sex
Morris10x_seurat_meta$Sex = "U"


##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Morris10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Cell_ID,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(Morris10x_raw_data_path,"Morris10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Morris10x_seurat_meta_final = Morris10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(Morris10x_seurat_meta_final) = Morris10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(Morris10x_seurat_meta_final) == ncol(Morris10x_seurat_raw)){
  Morris10x_seurat_raw@meta.data = Morris10x_seurat_meta_final
  saveRDS(Morris10x_seurat_raw,paste0(Morris10x_raw_data_path,"Morris10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}


