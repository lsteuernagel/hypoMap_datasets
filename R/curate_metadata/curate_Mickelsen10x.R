##########
### Info
##########

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125065

# no author clusters or other metadata found

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
Mickelsen10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Mickelsen10x/"
Mickelsen10x_seurat_raw = readRDS(paste0(Mickelsen10x_raw_data_path,"Mickelsen10x_seurat_raw.rds"))
Mickelsen10x_seurat_raw_meta = Mickelsen10x_seurat_raw@meta.data

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)
table(hypoMap_neurons_v1_curated_metadata$Diet[hypoMap_neurons_v1_curated_metadata$Dataset=="Mickelsen"])


##########
### Curate metadata
##########

Mickelsen10x_seurat_meta = Mickelsen10x_seurat_raw_meta

Mickelsen10x_seurat_meta$Dataset = "Mickelsen10x"
Mickelsen10x_seurat_meta$Sample_ID = paste0(Mickelsen10x_seurat_meta$Dataset,"_",Mickelsen10x_seurat_meta$Run_ID)
Mickelsen10x_seurat_meta$Tissue = "Hypothalamus"
Mickelsen10x_seurat_meta$Technology ="10xv2"
Mickelsen10x_seurat_meta$Author_Region = "Lateral hypothalamus"
Mickelsen10x_seurat_meta$Diet ="Normal chow"
Mickelsen10x_seurat_meta$Age = "3-6 weeks"
Mickelsen10x_seurat_meta$Pooled = "yes"
Mickelsen10x_seurat_meta$Author_Class = NA
Mickelsen10x_seurat_meta$Author_CellType = NA

# infer sex
Mickelsen10x_seurat_meta$inferred_sex = scUtils::infer_sex(Mickelsen10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.5,max_xist_male = 0.1)

# curate sex (different from author SRA anno!)
Mickelsen10x_seurat_meta$Sex = "M"
Mickelsen10x_seurat_meta$Sex[Mickelsen10x_seurat_meta$Sample_ID=="Mickelsen10x_SRR8441476"] = "F"

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Mickelsen10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Cell_ID,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(Mickelsen10x_raw_data_path,"Mickelsen10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Mickelsen10x_seurat_meta_final = Mickelsen10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(Mickelsen10x_seurat_meta_final) = Mickelsen10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(Mickelsen10x_seurat_meta_final) == ncol(Mickelsen10x_seurat_raw)){
  Mickelsen10x_seurat_raw@meta.data = Mickelsen10x_seurat_meta_final
  saveRDS(Mickelsen10x_seurat_raw,paste0(Mickelsen10x_raw_data_path,"Mickelsen10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}





#Mickelsen10x
