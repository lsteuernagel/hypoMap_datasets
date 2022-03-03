
##########
### Info
##########

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172204
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8184210/

# cannot find age and strain info

##########
### Load data
##########

library(dplyr)
library(Seurat)

## read raw file and extract metadata
Affinati10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Affinati10x/"
Affinati10x_seurat_raw = readRDS(paste0(Affinati10x_raw_data_path,"Affinati10x_seurat_raw.rds"))
Affinati10x_seurat_raw_meta = Affinati10x_seurat_raw@meta.data
# load author metadata
Affinati10x_seurat_author_meta = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE172204_metadata.csv",data.table = FALSE)

# add barcode as exra column
Affinati10x_seurat_raw_meta$Barcode = stringr::str_extract(Affinati10x_seurat_raw_meta$Cell_ID,pattern = "[ACGT]+\\-[0-9]")
Affinati10x_seurat_author_meta$Barcode = stringr::str_extract(Affinati10x_seurat_author_meta$ID,pattern = "[ACGT]+\\-[0-9]")
# match samples and join
matched_Sample_names = match_sample_names(table_A = Affinati10x_seurat_author_meta,table_B = Affinati10x_seurat_raw_meta,sample_col_A = "Sample",sample_col_B = "Run_ID",barcode_col = "Barcode")

##########
### Curate metadata sample level
##########

Affinati10x_seurat_author_meta = dplyr::left_join(Affinati10x_seurat_author_meta,matched_Sample_names,by=c("Sample"="Samples_A")) %>% dplyr::rename(Run_ID = Samples_B)
# reduce to per sample
Affinati10x_seurat_author_meta_perSample = Affinati10x_seurat_author_meta %>% dplyr::select(Sample,Run_ID,Run) %>%
  dplyr::distinct(.keep_all = TRUE)
# add to metadata
Affinati10x_seurat_meta = dplyr::left_join(Affinati10x_seurat_raw_meta,Affinati10x_seurat_author_meta_perSample,by="Run_ID")

### set other columns
Affinati10x_seurat_meta$Dataset = "Affinati10x"
Affinati10x_seurat_meta$Technology ="10xv3"
Affinati10x_seurat_meta$Strain = NA#"C57Bl6/J"
Affinati10x_seurat_meta$Diet ="Normal chow"
Affinati10x_seurat_meta$Pooled ="Yes"
Affinati10x_seurat_meta$Age = NA

##########
### Curate metadata cell level
##########

# Curate Author Cell types
Affinati10x_seurat_author_meta$Author_Class = gsub("_NA","",paste0(Affinati10x_seurat_author_meta$Cell_type))
Affinati10x_seurat_author_meta$Author_CellType = gsub("_NA","",paste0(Affinati10x_seurat_author_meta$Cell_type,"_",Affinati10x_seurat_author_meta$Cell_cluster,"_",Affinati10x_seurat_author_meta$Neuron_cluster))
# rework cell classes to one schema
table(Affinati10x_seurat_author_meta$Author_Class)
Affinati10x_seurat_author_meta$Author_Class[Affinati10x_seurat_author_meta$Author_Class=="Neuron"]="Neurons"
Affinati10x_seurat_author_meta$Author_Class[Affinati10x_seurat_author_meta$Author_Class=="Astro"]="Astrocytes"
Affinati10x_seurat_author_meta$Author_Class[grepl("Endo",Affinati10x_seurat_author_meta$Author_Class)]="Endothelial"
Affinati10x_seurat_author_meta$Author_Class[grepl("Mural",Affinati10x_seurat_author_meta$Author_Class)]="Mural"
Affinati10x_seurat_author_meta$Author_Class[grepl("Oligo",Affinati10x_seurat_author_meta$Author_Class)]="Oligodendrocytes"
Affinati10x_seurat_author_meta$Author_Class[grepl("Ependym",Affinati10x_seurat_author_meta$Author_Class)]="Ependymal"
Affinati10x_seurat_author_meta$Author_Class[Affinati10x_seurat_author_meta$Author_Class=="Micro"]="Immune"
Affinati10x_seurat_author_meta$Author_Class[Affinati10x_seurat_author_meta$Author_Class=="Tany"]="Tanycytes"
Affinati10x_seurat_author_meta$Author_Class[Affinati10x_seurat_author_meta$Author_Class=="OPC"]="NG/OPC"

Affinati10x_seurat_meta = dplyr::left_join(Affinati10x_seurat_meta,Affinati10x_seurat_author_meta %>% dplyr::select(Run_ID,Barcode,Author_CellType,Author_Class) ,by=c("Run_ID"="Run_ID","Barcode"="Barcode"))

##########
### other metadata
##########

# infer sex
Affinati10x_seurat_meta$inferred_sex = infer_sex(Affinati10x_seurat_raw,sample_column="Run_ID",id_column="Cell_ID",min_xist_female = 0.7,max_xist_male = 0.1)

# curate some sample ids
Affinati10x_seurat_meta$Sample_ID = paste0("Affinati10x_",Affinati10x_seurat_meta$Sample)

# mark cells that authors excluded (or did not annotate)
Affinati10x_seurat_meta$Author_Exclude = "no"
Affinati10x_seurat_meta$Author_Exclude[Affinati10x_seurat_meta$Author_Class=="NA" | is.na(Affinati10x_seurat_meta$Author_Class)] = "yes"

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Affinati10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(Affinati10x_raw_data_path,"Affinati10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Affinati10x_seurat_meta_final = Affinati10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`, Run10x = Run,Technology,Strain,Diet,Pooled,Age,Author_Region=tissue,inferred_sex, nCount_RNA, nFeature_RNA,Author_Exclude,Author_Class,Author_CellType  ) %>%
  as.data.frame()
rownames(Affinati10x_seurat_meta_final) = Affinati10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(Affinati10x_seurat_meta_final) == ncol(Affinati10x_seurat_raw)){
  Affinati10x_seurat_raw@meta.data = Affinati10x_seurat_meta_final
  saveRDS(Affinati10x_seurat_raw,paste0(Affinati10x_raw_data_path,"Affinati10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}




