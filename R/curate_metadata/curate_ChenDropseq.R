##########
### Info
##########

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5782816/
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87544
# four normal fed mice and three food-deprived (24 hr) mice were used
# All single-cell RNA-seq experiments were performed in five batches within 3 months
# Samples do not add up with metadata, I assume that a mouse was not split across batches --> then there are 5 normal mice and 3 FD, some with multiple samples per batch

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
ChenDropseq_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/ChenDropseq/"
ChenDropseq_seurat_raw = readRDS(paste0(ChenDropseq_raw_data_path,"ChenDropseq_seurat_raw.rds"))
ChenDropseq_seurat_raw_meta = ChenDropseq_seurat_raw@meta.data

# load author metadata
ChenDropseq_seurat_author_meta= data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv",data.table = F,header = TRUE)

# add barcode as exra column
ChenDropseq_seurat_raw_meta$Barcode = stringr::str_extract(ChenDropseq_seurat_raw_meta$Cell_ID,pattern = "[ACGT]{5,}\\-?[0-9]?")
ChenDropseq_seurat_author_meta$Barcode = stringr::str_extract(ChenDropseq_seurat_author_meta$V1,pattern = "[ACGT]{5,}\\-?[0-9]?")

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)

##########
### set general metadata
##########

### set other columns
ChenDropseq_seurat_raw_meta$Dataset = "ChenDropseq"
ChenDropseq_seurat_raw_meta$Technology ="Dropseq"
ChenDropseq_seurat_raw_meta$Strain = "C57BL6 X DBA2 "
ChenDropseq_seurat_raw_meta$Pooled = NA#"Yes"
ChenDropseq_seurat_raw_meta$Author_Region = "Hypothalamus"

## check v1 metadata for diet and age:
#matched_Sample_names_v1 = scUtils::match_sample_names(table_A = ChenDropseq_seurat_raw_meta,table_B = hypoMap_neurons_v1_curated_metadata,sample_col_A = "Run_ID",sample_col_B = "Age",barcode_col = "Barcode",min_pct = 0.1)
ChenDropseq_seurat_raw_meta$Age = "6+ weeks"
ChenDropseq_seurat_raw_meta$Sex = "F"

##########
### Clean up cell annotations
##########

# set Author_CellType
ChenDropseq_seurat_author_meta$Author_CellType = ChenDropseq_seurat_author_meta$SVM_clusterID

## build author class
ChenDropseq_seurat_author_meta$Author_Class=NA
ChenDropseq_seurat_author_meta$Author_Class[!is.na(ChenDropseq_seurat_author_meta$Author_CellType)]="Neurons"
ChenDropseq_seurat_author_meta$Author_Class[grepl("Ependy",ChenDropseq_seurat_author_meta$Author_CellType)]="Ependymal"
ChenDropseq_seurat_author_meta$Author_Class[grepl("Macro|Micro",ChenDropseq_seurat_author_meta$Author_CellType)]="Immune"
ChenDropseq_seurat_author_meta$Author_Class[grepl("Tany",ChenDropseq_seurat_author_meta$Author_CellType)]="Tanycytes"
ChenDropseq_seurat_author_meta$Author_Class[grepl("IMO|MO",ChenDropseq_seurat_author_meta$Author_CellType)]="Oligodendrocytes"
ChenDropseq_seurat_author_meta$Author_Class[grepl("OPC|POPC",ChenDropseq_seurat_author_meta$Author_CellType)]="NG/OPC"
ChenDropseq_seurat_author_meta$Author_Class[grepl("Astro",ChenDropseq_seurat_author_meta$Author_CellType)]="Astrocytes"
ChenDropseq_seurat_author_meta$Author_Class[grepl("Epith",ChenDropseq_seurat_author_meta$Author_CellType)]="Endothelial"
ChenDropseq_seurat_author_meta$Author_Class[grepl("SCO|zothers",ChenDropseq_seurat_author_meta$Author_CellType)]=NA # set them to NA, need to propagate good labels
table(ChenDropseq_seurat_author_meta$Author_Class)

##########
### Match data on cell level
##########

# "[ACGT]{5,}\\-?[0-9]?"
ChenDropseq_seurat_author_meta$Author_Sample = stringr::str_extract(ChenDropseq_seurat_author_meta$V1,pattern = "B[0-9]{1,}")
ChenDropseq_seurat_author_meta$Author_Sample = paste0(ChenDropseq_seurat_author_meta$Author_Sample,stringr::str_extract(ChenDropseq_seurat_author_meta$V1,pattern = "_[a-zA-Z]+"))

# these samples do not match properly with the SRR ids:
# table(ChenDropseq_seurat_raw_meta$Run_ID)
# matched_Sample_names = scUtils::match_sample_names(table_A = ChenDropseq_seurat_author_meta,table_B = ChenDropseq_seurat_raw_meta,sample_col_A = "Author_Sample",sample_col_B = "Run_ID",barcode_col = "Barcode",min_pct = 0.1)

# I solve this by joining just on barcode --> but then we need to take care of duplicated barcodes
duplicated_barcodes = ChenDropseq_seurat_author_meta$Barcode[duplicated(ChenDropseq_seurat_author_meta$Barcode)]
message("length duplicated barcodes: ",length(duplicated_barcodes))

# join non-duplicated barcodes
ChenDropseq_seurat_meta_all = dplyr::left_join(ChenDropseq_seurat_raw_meta %>% dplyr::filter(! Barcode %in% duplicated_barcodes) ,ChenDropseq_seurat_author_meta %>% dplyr::select(Author_Class,Author_CellType,Barcode) %>% dplyr::filter(! Barcode %in% duplicated_barcodes),by=c("Barcode"="Barcode"))
# cannot join duplicated barcodes without knowing the ids -- just use srr metadata
ChenDropseq_seurat_meta_dup = ChenDropseq_seurat_raw_meta %>% dplyr::filter( Barcode %in% duplicated_barcodes)
# rbind:
ChenDropseq_seurat_meta = dplyr::bind_rows(ChenDropseq_seurat_meta_all,ChenDropseq_seurat_meta_dup)
# reorder
rownames(ChenDropseq_seurat_meta) = ChenDropseq_seurat_meta$Cell_ID
ChenDropseq_seurat_meta = ChenDropseq_seurat_meta[match(rownames(ChenDropseq_seurat_raw_meta),rownames(ChenDropseq_seurat_meta)),]
message("nrow ChenDropseq_seurat_meta: ",nrow(ChenDropseq_seurat_meta)," , should be: ",nrow(ChenDropseq_seurat_raw_meta))

##########
### Curate remainder of metadata
##########

# infer sex
ChenDropseq_seurat_meta$inferred_sex = scUtils::infer_sex(ChenDropseq_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

# curate some sample ids
ChenDropseq_seurat_meta$Sample_ID = paste0("ChenDropseq_",ChenDropseq_seurat_meta$Run_ID)

# mark cells that authors excluded (or did not annotate)
ChenDropseq_seurat_meta$Author_Exclude = "no"
ChenDropseq_seurat_meta$Author_Exclude[is.na(ChenDropseq_seurat_meta$Author_CellType)] = "yes"

## add conditions
ChenDropseq_seurat_meta$Diet = ChenDropseq_seurat_meta$mouse_status
ChenDropseq_seurat_meta$Diet[ChenDropseq_seurat_meta$Diet=="normal"] = "Normal chow"
ChenDropseq_seurat_meta$Diet[ChenDropseq_seurat_meta$Diet=="food deprived"] = "Fasted"
ChenDropseq_seurat_meta$Sex = "F"

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = ChenDropseq_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(ChenDropseq_raw_data_path,"ChenDropseq_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
ChenDropseq_seurat_meta_final = ChenDropseq_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=Strain,Diet,Pooled,Age,Author_Region,inferred_sex,Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Exclude,Author_Class,Author_CellType  ) %>%
  as.data.frame()
rownames(ChenDropseq_seurat_meta_final) = ChenDropseq_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(ChenDropseq_seurat_meta_final) == ncol(ChenDropseq_seurat_raw)){
  ChenDropseq_seurat_raw@meta.data = ChenDropseq_seurat_meta_final
  saveRDS(ChenDropseq_seurat_raw,paste0(ChenDropseq_raw_data_path,"ChenDropseq_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}



