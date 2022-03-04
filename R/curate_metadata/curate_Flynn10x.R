##########
### Info
##########

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146692
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7595735/

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
Flynn10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Flynn10x/"
Flynn10x_seurat_raw = readRDS(paste0(Flynn10x_raw_data_path,"Flynn10x_seurat_raw.rds"))
Flynn10x_seurat_raw_meta = Flynn10x_seurat_raw@meta.data

# load author metadata
Flynn10x_seurat_author_meta= data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/vph-2020_global_metadata.txt",data.table = F,header = TRUE)

# add barcode as exra column
Flynn10x_seurat_raw_meta$Barcode = stringr::str_extract(Flynn10x_seurat_raw_meta$Cell_ID,pattern = "[ACGT]{5,}\\-?[0-9]?")
Flynn10x_seurat_author_meta$Barcode = stringr::str_extract(Flynn10x_seurat_author_meta$index,pattern = "[ACGT]{5,}\\-?[0-9]?")

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)

##########
### Curate metadata sample level
##########

# match samples and join
matched_Sample_names = scUtils::match_sample_names(table_A = Flynn10x_seurat_author_meta,table_B = Flynn10x_seurat_raw_meta,sample_col_A = "sampleid",sample_col_B = "Run_ID",barcode_col = "Barcode")

Flynn10x_seurat_author_meta = dplyr::left_join(Flynn10x_seurat_author_meta,matched_Sample_names,by=c("sampleid"="Samples_A")) %>% dplyr::rename(Run_ID = Samples_B)
# reduce to per sample
Flynn10x_seurat_author_meta_perSample = Flynn10x_seurat_author_meta %>% dplyr::select(sampleid,Run_ID,batch,sample_name,chemistry) %>% dplyr::distinct(.keep_all = TRUE)
# add to metadata
Flynn10x_seurat_meta = dplyr::left_join(Flynn10x_seurat_raw_meta,Flynn10x_seurat_author_meta_perSample,by="Run_ID")

### set other columns
Flynn10x_seurat_meta$Dataset = "Flynn10x"
Flynn10x_seurat_meta$Diet ="Normal chow"
Flynn10x_seurat_meta$Pooled ="Yes"
Flynn10x_seurat_meta$Age = "3-6 weeks"

##########
### Curate metadata cell level
##########

Flynn10x_seurat_author_meta$Author_Class = Flynn10x_seurat_author_meta$celltype_broad
# rework cell classes to one schema
Flynn10x_seurat_author_meta$Author_Class[Flynn10x_seurat_author_meta$Author_Class=="neuron"]="Neurons"
Flynn10x_seurat_author_meta$Author_Class[Flynn10x_seurat_author_meta$Author_Class=="astrocyte"]="Astrocytes"
Flynn10x_seurat_author_meta$Author_Class[Flynn10x_seurat_author_meta$Author_Class=="endothelial"]="Endothelial"
Flynn10x_seurat_author_meta$Author_Class[Flynn10x_seurat_author_meta$Author_Class=="fibroblast-like"]="Fibroblasts"
Flynn10x_seurat_author_meta$Author_Class[Flynn10x_seurat_author_meta$Author_Class=="oligodendrocyte"]="Oligodendrocytes"
Flynn10x_seurat_author_meta$Author_Class[Flynn10x_seurat_author_meta$Author_Class=="pericyte"]="Pericytes"
Flynn10x_seurat_author_meta$Author_Class[Flynn10x_seurat_author_meta$Author_Class=="macrophage"]="Immune"
Flynn10x_seurat_author_meta$Author_Class[Flynn10x_seurat_author_meta$Author_Class=="microglia"]="Immune"
Flynn10x_seurat_author_meta$Author_Class[Flynn10x_seurat_author_meta$Author_Class=="tancycte"]="Tanycytes"
table(Flynn10x_seurat_author_meta$Author_Class)
# add cellclusters
Flynn10x_seurat_author_meta$Author_CellType = paste0(Flynn10x_seurat_author_meta$Author_Class,"_",Flynn10x_seurat_author_meta$cluster_revised_R1)
Flynn10x_seurat_author_meta$Author_CellType[Flynn10x_seurat_author_meta$Author_CellType=="NA_NA"]=NA

# join
Flynn10x_seurat_meta = dplyr::left_join(Flynn10x_seurat_meta,Flynn10x_seurat_author_meta %>% dplyr::select(Run_ID,Barcode,Author_CellType,Author_Class) ,by=c("Run_ID"="Run_ID","Barcode"="Barcode"))

##########
### other metadata
##########

Flynn10x_seurat_meta$Technology = paste0("10x",Flynn10x_seurat_meta$`10x_chemistry`)

# infer sex
Flynn10x_seurat_meta$inferred_sex = scUtils::infer_sex(Flynn10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

# curate some sample ids
Flynn10x_seurat_meta$Sample_ID = paste0("Flynn10x_",Flynn10x_seurat_meta$sampleid)

# mark cells that authors excluded (or did not annotate)
Flynn10x_seurat_meta$Author_Exclude = "no"
Flynn10x_seurat_meta$Author_Exclude[Flynn10x_seurat_meta$Author_Class=="NA" | is.na(Flynn10x_seurat_meta$Author_Class)] = "yes"

# sex
Flynn10x_seurat_meta$Sex =Flynn10x_seurat_meta$sex
Flynn10x_seurat_meta$Sex[Flynn10x_seurat_meta$Sex =="male"]="M"
Flynn10x_seurat_meta$Sex[Flynn10x_seurat_meta$Sex =="female"]="F"

#otehr
Flynn10x_seurat_meta$Author_Region = Flynn10x_seurat_meta$source_name
Flynn10x_seurat_meta$Age

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Flynn10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Author_Exclude,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(Flynn10x_raw_data_path,"Flynn10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Flynn10x_seurat_meta_final = Flynn10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Exclude,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(Flynn10x_seurat_meta_final) = Flynn10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(Flynn10x_seurat_meta_final) == ncol(Flynn10x_seurat_raw)){
  Flynn10x_seurat_raw@meta.data = Flynn10x_seurat_meta_final
  saveRDS(Flynn10x_seurat_raw,paste0(Flynn10x_raw_data_path,"Flynn10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}






