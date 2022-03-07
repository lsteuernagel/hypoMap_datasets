##########
### Info
##########

# GSE130597

# https://www.ncbi.nlm.nih.gov/pubmed/31249056

# no metadata found
# tissue: lateral hypothalamus, strain: C57BL/6J

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
RossiDropseq_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/RossiDropseq/"
RossiDropseq_seurat_raw = readRDS(paste0(RossiDropseq_raw_data_path,"RossiDropseq_seurat_raw.rds"))
RossiDropseq_seurat_raw_meta = RossiDropseq_seurat_raw@meta.data

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)
table(hypoMap_neurons_v1_curated_metadata$Diet[hypoMap_neurons_v1_curated_metadata$Dataset=="Rossi"])

##########
### Curate metadata
##########

RossiDropseq_seurat_meta = RossiDropseq_seurat_raw_meta

RossiDropseq_seurat_meta$Dataset = "RossiDropseq"
RossiDropseq_seurat_meta$Sample_ID = paste0(RossiDropseq_seurat_meta$Dataset,"_",RossiDropseq_seurat_meta$Run_ID)
RossiDropseq_seurat_meta$Technology ="Dropseq"
RossiDropseq_seurat_meta$Author_Region = "Lateral hypothalamus"
RossiDropseq_seurat_meta$Diet = RossiDropseq_seurat_meta$Group
RossiDropseq_seurat_meta$Diet[RossiDropseq_seurat_meta$Diet=="high fat diet"] = "HFD_9-16weeks"
RossiDropseq_seurat_meta$Diet[RossiDropseq_seurat_meta$Diet=="control"] = "Normal chow"
RossiDropseq_seurat_meta$Age = "6+ weeks"
RossiDropseq_seurat_meta$Pooled = "no"
RossiDropseq_seurat_meta$Author_Class = NA
RossiDropseq_seurat_meta$Author_CellType = NA
RossiDropseq_seurat_meta$Sex = "M"

# infer sex
RossiDropseq_seurat_meta$inferred_sex = scUtils::infer_sex(RossiDropseq_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

# author ex set to no
RossiDropseq_seurat_meta$Author_Exclude = "no"


##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = RossiDropseq_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Cell_ID,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(RossiDropseq_raw_data_path,"RossiDropseq_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
RossiDropseq_seurat_meta_final = RossiDropseq_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(RossiDropseq_seurat_meta_final) = RossiDropseq_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(RossiDropseq_seurat_meta_final) == ncol(RossiDropseq_seurat_raw)){
  RossiDropseq_seurat_raw@meta.data = RossiDropseq_seurat_meta_final
  saveRDS(RossiDropseq_seurat_raw,paste0(RossiDropseq_raw_data_path,"RossiDropseq_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}


