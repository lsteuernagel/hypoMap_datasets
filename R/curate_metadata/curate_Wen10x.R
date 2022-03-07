##########
### Info
##########

# GSE132608

# https://www.ncbi.nlm.nih.gov/pubmed/32066983

# no metadata found
#  	strain: C57BL/6, age: 8 weeks, sampling enviroment: constant dark, tissue: suprachiasmatic nucleus

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
wen10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/wen10x/"
wen10x_seurat_raw = readRDS(paste0(wen10x_raw_data_path,"wen10x_seurat_raw.rds"))
wen10x_seurat_raw_meta = wen10x_seurat_raw@meta.data

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)
table(hypoMap_neurons_v1_curated_metadata$Diet[hypoMap_neurons_v1_curated_metadata$Dataset=="Wen_10X"])

##########
### Curate metadata
##########

wen10x_seurat_meta = wen10x_seurat_raw_meta

wen10x_seurat_meta$Dataset = "Wen10x"
wen10x_seurat_meta$Sample_ID = paste0(wen10x_seurat_meta$Dataset,"_",wen10x_seurat_meta$Run_ID)
wen10x_seurat_meta$Technology ="10xv2"
wen10x_seurat_meta$Author_Region = "Suprachiasmatic nucleus"
wen10x_seurat_meta$Diet =  "Normal chow"
#wen10x_seurat_meta$Author_Condition = wen10x_seurat_meta$sampling_enviroment
wen10x_seurat_meta$Age = "6+ weeks"
wen10x_seurat_meta$Pooled = "yes"
wen10x_seurat_meta$Author_Class = NA
wen10x_seurat_meta$Author_CellType = NA
wen10x_seurat_meta$Sex = "M"

# infer sex
wen10x_seurat_meta$inferred_sex = scUtils::infer_sex(wen10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

# author ex set to no
wen10x_seurat_meta$Author_Exclude = NA


##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = wen10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Cell_ID,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(wen10x_raw_data_path,"wen10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
wen10x_seurat_meta_final = wen10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(wen10x_seurat_meta_final) = wen10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(wen10x_seurat_meta_final) == ncol(wen10x_seurat_raw)){
  wen10x_seurat_raw@meta.data = wen10x_seurat_meta_final
  saveRDS(wen10x_seurat_raw,paste0(wen10x_raw_data_path,"wen10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}


