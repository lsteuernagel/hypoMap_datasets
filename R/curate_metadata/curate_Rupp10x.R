
##########
### Info
##########

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172461
#https://www.biorxiv.org/content/10.1101/2021.12.10.472115v1

# cannot find age and strain info

# see also Affinati dataset for another one from Rupp Lab

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
Rupp10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Rupp10x/"
Rupp10x_seurat_raw = readRDS(paste0(Rupp10x_raw_data_path,"Rupp10x_seurat_raw.rds"))
Rupp10x_seurat_raw_meta = Rupp10x_seurat_raw@meta.data

# load author metadata
Rupp10x_seurat_author_meta = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE172461_metadata.csv",data.table = FALSE)

# add barcode as exra column
Rupp10x_seurat_raw_meta$Barcode = stringr::str_extract(Rupp10x_seurat_raw_meta$Cell_ID,pattern = "[ACGT]{5,}")
Rupp10x_seurat_author_meta$Barcode = stringr::str_extract(Rupp10x_seurat_author_meta$ID,pattern = "[ACGT]{5,}")

##########
### Curate metadata sample level
##########

# match samples and join
matched_Sample_names = scUtils::match_sample_names(table_A = Rupp10x_seurat_author_meta,table_B = Rupp10x_seurat_raw_meta,sample_col_A = "Sample",sample_col_B = "Run_ID",barcode_col = "Barcode")

Rupp10x_seurat_author_meta = dplyr::left_join(Rupp10x_seurat_author_meta,matched_Sample_names,by=c("Sample"="Samples_A")) %>% dplyr::rename(Run_ID = Samples_B)
# reduce to per sample
Rupp10x_seurat_author_meta_perSample = Rupp10x_seurat_author_meta %>% dplyr::select(Sample,Run_ID,Run,Treatment) %>%
  dplyr::distinct(.keep_all = TRUE)
# add to metadata
Rupp10x_seurat_meta = dplyr::left_join(Rupp10x_seurat_raw_meta,Rupp10x_seurat_author_meta_perSample,by="Run_ID")

### set other columns
Rupp10x_seurat_meta$Dataset = "Rupp10x"
Rupp10x_seurat_meta$Technology ="10xv3"
Rupp10x_seurat_meta$Strain = NA#"C57Bl6/J"
Rupp10x_seurat_meta$Diet ="Normal chow"
Rupp10x_seurat_meta$Pooled ="Yes"
Rupp10x_seurat_meta$Age = "6+ weeks"
Rupp10x_seurat_meta$Author_Region = "Mediobasal hypothalamus"

##########
### Curate metadata cell level
##########

# Curate Author Cell types
Rupp10x_seurat_author_meta$Author_Class = gsub("_NA","",paste0(Rupp10x_seurat_author_meta$Cell_type))
Rupp10x_seurat_author_meta$Author_CellType = Rupp10x_seurat_author_meta$Author_Class
Rupp10x_seurat_author_meta$Author_CellType[Rupp10x_seurat_author_meta$Author_CellType=="Neuron"] = Rupp10x_seurat_author_meta$Neuron_cluster[Rupp10x_seurat_author_meta$Author_CellType=="Neuron"]
# rework cell classes to one schema
table(Rupp10x_seurat_author_meta$Author_Class)
Rupp10x_seurat_author_meta$Author_Class[Rupp10x_seurat_author_meta$Author_Class=="Neuron"]="Neurons"
Rupp10x_seurat_author_meta$Author_Class[Rupp10x_seurat_author_meta$Author_Class=="Astro"]="Astrocytes"
Rupp10x_seurat_author_meta$Author_Class[grepl("Endothelial",Rupp10x_seurat_author_meta$Author_Class)]="Endothelial"
Rupp10x_seurat_author_meta$Author_Class[grepl("Mural",Rupp10x_seurat_author_meta$Author_Class)]="Mural"
Rupp10x_seurat_author_meta$Author_Class[grepl("Oligo",Rupp10x_seurat_author_meta$Author_Class)]="Oligodendrocytes"
Rupp10x_seurat_author_meta$Author_Class[grepl("Ependym",Rupp10x_seurat_author_meta$Author_Class)]="Ependymal"
Rupp10x_seurat_author_meta$Author_Class[Rupp10x_seurat_author_meta$Author_Class=="Micro"]="Immune"
Rupp10x_seurat_author_meta$Author_Class[Rupp10x_seurat_author_meta$Author_Class=="OPC"]="NG/OPC"

# this dataset is clean so duplicated barcodes within samples are not a problem and we can join directly
Rupp10x_seurat_meta = dplyr::left_join(Rupp10x_seurat_meta,Rupp10x_seurat_author_meta %>% dplyr::select(Run_ID,Barcode,Author_CellType,Author_Class) ,by=c("Run_ID"="Run_ID","Barcode"="Barcode"))

##########
### other metadata
##########

# infer sex
Rupp10x_seurat_meta$inferred_sex = scUtils::infer_sex(Rupp10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

# curate some sample ids
Rupp10x_seurat_meta$Sample_ID = gsub("Sample_","",gsub("-","_",paste0("Rupp10x_",Rupp10x_seurat_meta$Sample)))

# add 10x run
Rupp10x_seurat_meta$Run10x = Rupp10x_seurat_meta$Run

# add treament
Rupp10x_seurat_meta$Author_Condition = Rupp10x_seurat_meta$treatment

# mark cells that authors excluded (or did not annotate)
Rupp10x_seurat_meta$Author_Exclude = "no"
Rupp10x_seurat_meta$Author_Exclude[Rupp10x_seurat_meta$Author_Class=="NA" | is.na(Rupp10x_seurat_meta$Author_Class)] = "yes"
Rupp10x_seurat_meta$Author_Exclude[Rupp10x_seurat_meta$Author_Class=="Doublets"] = "yes"

# I don't know the gender --> probaly later use inferred
Rupp10x_seurat_meta$Sex = NA

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Rupp10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA,-Sex,-Author_Exclude)
data.table::fwrite(per_sample_summary,file = paste0(Rupp10x_raw_data_path,"Rupp10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Rupp10x_seurat_meta_final = Rupp10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`, Run10x,Technology,Strain,Diet,Pooled,Age,Author_Condition,Author_Region,inferred_sex,Sex,nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Exclude,Author_Class,Author_CellType  ) %>%
  as.data.frame()
rownames(Rupp10x_seurat_meta_final) = Rupp10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(Rupp10x_seurat_meta_final) == ncol(Rupp10x_seurat_raw)){
  Rupp10x_seurat_raw@meta.data = Rupp10x_seurat_meta_final
  saveRDS(Rupp10x_seurat_raw,paste0(Rupp10x_raw_data_path,"Rupp10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}




