
##########
### Info
##########

# GSE113576

# get metdata
# Moffit10x_seurat_author_metadata_tmp = readxl::read_xlsx("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypothalamus/moffit/NIHMS1024025-supplement-Table_S1.xlsx",skip = 1)
# data.table::fwrite(Moffit10x_seurat_author_metadata_tmp,"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/moffit10x_metadata.tsv",sep = "\t")
# # Moffit10x

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
Moffit10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Moffit10x/"
Moffit10x_seurat_raw = readRDS(paste0(Moffit10x_raw_data_path,"Moffit10x_seurat_raw.rds"))
Moffit10x_seurat_raw_meta = Moffit10x_seurat_raw@meta.data

# load author metadata
Moffit10x_seurat_author_meta= data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/moffit10x_metadata.tsv",data.table = F,header = TRUE)

# add barcode as exra column
Moffit10x_seurat_raw_meta$Barcode = stringr::str_extract(Moffit10x_seurat_raw_meta$Cell_ID,pattern = "[ACGT]{5,}") # leave out -1
Moffit10x_seurat_author_meta$Barcode = stringr::str_extract(Moffit10x_seurat_author_meta$`Cell name`,pattern = "[ACGT]{5,}")

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)
table(hypoMap_neurons_v1_curated_metadata$Age[hypoMap_neurons_v1_curated_metadata$Dataset=="Moffit"])

##########
### Curate metadata
##########

### set other columns
Moffit10x_seurat_raw_meta$Dataset = "Moffit10x"
Moffit10x_seurat_raw_meta$Sample_ID = paste0(Moffit10x_seurat_raw_meta$Dataset,"_",Moffit10x_seurat_raw_meta$Run_ID)
Moffit10x_seurat_raw_meta$Diet ="Normal chow"
Moffit10x_seurat_raw_meta$Author_Region = "Preoptic region"
Moffit10x_seurat_raw_meta$Diet ="Normal chow"
Moffit10x_seurat_raw_meta$Pooled ="Yes"
Moffit10x_seurat_raw_meta$Sex[Moffit10x_seurat_raw_meta$sex=="male"] = "M"
Moffit10x_seurat_raw_meta$Sex[Moffit10x_seurat_raw_meta$sex=="female"] = "F"
Moffit10x_seurat_raw_meta$Age = "6+ weeks"
Moffit10x_seurat_raw_meta$Strain = "C57Bl6/J"

##########
### Curate metadata cell level
##########

Moffit10x_seurat_author_meta$Author_Sample = paste0(Moffit10x_seurat_author_meta$`Replicate number`,"_",Moffit10x_seurat_author_meta$Sex)

# match samples and join
matched_Sample_names = scUtils::match_sample_names(table_A = Moffit10x_seurat_author_meta,table_B = Moffit10x_seurat_raw_meta,sample_col_A = "Author_Sample",sample_col_B = "Run_ID",barcode_col = "Barcode")
# based on this barcodes match 1:1 !!!

# add run id
Moffit10x_seurat_author_meta = dplyr::left_join(Moffit10x_seurat_author_meta,matched_Sample_names,by=c("Author_Sample"="Samples_A")) %>% dplyr::rename(Run_ID = Samples_B)

# curate celltypes
Moffit10x_seurat_author_meta$Author_CellType = Moffit10x_seurat_author_meta$`Neuronal cluster (determined from clustering of inhibitory or excitatory neurons)`
Moffit10x_seurat_author_meta$Author_CellType[Moffit10x_seurat_author_meta$Author_CellType==""] = Moffit10x_seurat_author_meta$`Non-neuronal cluster (determined from clustering of all cells)`[Moffit10x_seurat_author_meta$Author_CellType==""]
Moffit10x_seurat_author_meta$Author_CellType = gsub(" ","_",Moffit10x_seurat_author_meta$Author_CellType)

# rework cell classes to one schema
Moffit10x_seurat_author_meta$Author_Class = Moffit10x_seurat_author_meta$`Cell class (determined from clustering of all cells)`
Moffit10x_seurat_author_meta$Author_Class[grepl("Excitatory|Inhibitory",Moffit10x_seurat_author_meta$Author_Class)]="Neurons"
Moffit10x_seurat_author_meta$Author_Class[Moffit10x_seurat_author_meta$Author_Class=="Astrocyte"]="Astrocytes"
Moffit10x_seurat_author_meta$Author_Class[grepl("Endothelial",Moffit10x_seurat_author_meta$Author_Class)]="Endothelial"
Moffit10x_seurat_author_meta$Author_Class[grepl("Mural",Moffit10x_seurat_author_meta$Author_Class)]="Mural"
Moffit10x_seurat_author_meta$Author_Class[grepl("Oligo|oligo",Moffit10x_seurat_author_meta$Author_Class)]="Oligodendrocytes"
Moffit10x_seurat_author_meta$Author_Class[grepl("Ependym",Moffit10x_seurat_author_meta$Author_Class)]="Ependymal"
Moffit10x_seurat_author_meta$Author_Class[Moffit10x_seurat_author_meta$Author_Class=="Macrophage"]="Immune"
Moffit10x_seurat_author_meta$Author_Class[Moffit10x_seurat_author_meta$Author_Class=="Microglia"]="Immune"
Moffit10x_seurat_author_meta$Author_Class[Moffit10x_seurat_author_meta$Author_Class=="Tanycyte"]="Tanycytes"
table(Moffit10x_seurat_author_meta$Author_Class)

# join
Moffit10x_seurat_meta = dplyr::left_join(Moffit10x_seurat_raw_meta,Moffit10x_seurat_author_meta %>% dplyr::select(Run_ID,Barcode,Author_CellType,Author_Class) ,by=c("Run_ID"="Run_ID","Barcode"="Barcode"))

# exlude cell_ids
Moffit10x_seurat_meta$Author_Exclude = "no"
Moffit10x_seurat_meta$Author_Exclude[is.na(Moffit10x_seurat_meta$Author_Class) | Moffit10x_seurat_meta$Author_Class %in% c("Ambiguous","Unstable") ] = "yes"
# An entry of ‘ambiguous’ indicates that the cell was marked as a putative doublet and was not further
#analyzed. An entry of ‘unstable’ indicates that the cell was not part of a stable cluster and was not further analyzed

##########
### other metadata
##########


Moffit10x_seurat_meta$Technology = paste0("10xv2")

# infer sex
Moffit10x_seurat_meta$inferred_sex = scUtils::infer_sex(Moffit10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)


##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Moffit10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Author_Exclude,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(Moffit10x_raw_data_path,"Moffit10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Moffit10x_seurat_meta_final = Moffit10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Exclude,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(Moffit10x_seurat_meta_final) = Moffit10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(Moffit10x_seurat_meta_final) == ncol(Moffit10x_seurat_raw)){
  Moffit10x_seurat_raw@meta.data = Moffit10x_seurat_meta_final
  saveRDS(Moffit10x_seurat_raw,paste0(Moffit10x_raw_data_path,"Moffit10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}
