
##########
### Info
##########

# GSE132730


# metadata ix extracted from this rds: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3890941
# and then saved as a table
# romanovDev_folder_static = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypothalamus/romanovDev/"
# GSM3890941_Hypothalamus_traject_P23_b2 = readRDS(paste0(romanovDev_folder_static,"geo_seurat/","GSM3890941_Hypothalamus_traject_P23_srt_annotated_wo_blood.rds"))
# hypothalamus_romanovDev_seurat <- GSM3890941_Hypothalamus_traject_P23_b2
# hypothalamus_romanovDev_meta = hypothalamus_romanovDev_seurat@meta.data
#data.table::fwrite(hypothalamus_romanovDev_meta,"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE132730_GSM3890941_romanovDev_metadata.tsv",sep="\t")

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
RomanovDev10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/RomanovDev10x/"
RomanovDev10x_seurat_raw = readRDS(paste0(RomanovDev10x_raw_data_path,"RomanovDev10x_seurat_raw.rds"))
RomanovDev10x_seurat_raw_meta = RomanovDev10x_seurat_raw@meta.data

# load author metadata
RomanovDev10x_seurat_author_meta= data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE132730_GSM3890941_romanovDev_metadata.tsv",data.table = F,header = TRUE)

# add barcode as exra column
RomanovDev10x_seurat_raw_meta$Barcode = stringr::str_extract(RomanovDev10x_seurat_raw_meta$Cell_ID,pattern = "[ACGT]{5,}") # leave out -1
RomanovDev10x_seurat_author_meta$Barcode = stringr::str_extract(RomanovDev10x_seurat_author_meta$CellID,pattern = "[ACGT]{5,}")

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)
table(hypoMap_neurons_v1_curated_metadata$Technology[hypoMap_neurons_v1_curated_metadata$Dataset=="Romanov_Dev"])

##########
### Curate metadata
##########

### set other columns
RomanovDev10x_seurat_raw_meta$Dataset = "RomanovDev10x"
RomanovDev10x_seurat_raw_meta$Diet ="Normal chow"
RomanovDev10x_seurat_raw_meta$Author_Region = "Hypothalamus"
RomanovDev10x_seurat_raw_meta$Pooled ="Yes"
RomanovDev10x_seurat_raw_meta$Sex = "U"
RomanovDev10x_seurat_raw_meta$Age = "3-6 weeks"
RomanovDev10x_seurat_raw_meta$Strain =RomanovDev10x_seurat_raw_meta$strain

# sample id based on SRR because it does not fit with the author sampleid
RomanovDev10x_seurat_raw_meta$Sample_ID = paste0(RomanovDev10x_seurat_raw_meta$Dataset,"_",RomanovDev10x_seurat_raw_meta$Run_ID)

##########
### Curate metadata cell level
##########

# match samples and join
# matched_Sample_names = scUtils::match_sample_names(table_A = RomanovDev10x_seurat_author_meta,table_B = RomanovDev10x_seurat_raw_meta,sample_col_A = "SampleID",sample_col_B = "Run_ID",barcode_col = "Barcode")
# this did not work well: one sample but to SRR numbers :/

# curate celltypes
RomanovDev10x_seurat_author_meta$Author_CellType = NA
RomanovDev10x_seurat_author_meta$Author_Class = RomanovDev10x_seurat_author_meta$Class

# rework cell classes to one schema --> pretty good already!
table(RomanovDev10x_seurat_author_meta$Author_Class )
RomanovDev10x_seurat_author_meta$Author_Class[grepl("Vascular",RomanovDev10x_seurat_author_meta$Author_Class)]="Endothelial"
RomanovDev10x_seurat_author_meta$Author_Class[grepl("Oligo",RomanovDev10x_seurat_author_meta$Author_Class)]="Oligodendrocytes"

# join (there seem to be no duplicated barcodes !)
RomanovDev10x_seurat_meta = dplyr::left_join(RomanovDev10x_seurat_raw_meta,RomanovDev10x_seurat_author_meta %>% dplyr::select(Barcode,Author_CellType,Author_Class) ,by=c("Barcode"="Barcode"))

# exlude cell_ids
RomanovDev10x_seurat_meta$Author_Exclude = "no"
RomanovDev10x_seurat_meta$Author_Exclude[is.na(RomanovDev10x_seurat_meta$Author_Class) | RomanovDev10x_seurat_meta$Author_Class %in% c("Excluded") ] = "yes"

##########
### other metadata
##########

RomanovDev10x_seurat_meta$Technology = paste0("10xv2")

# infer sex
RomanovDev10x_seurat_meta$inferred_sex = scUtils::infer_sex(RomanovDev10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = RomanovDev10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Cell_ID,-Barcode,-Author_CellType,-Author_Exclude,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(RomanovDev10x_raw_data_path,"RomanovDev10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
RomanovDev10x_seurat_meta_final = RomanovDev10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Exclude,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(RomanovDev10x_seurat_meta_final) = RomanovDev10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(RomanovDev10x_seurat_meta_final) == ncol(RomanovDev10x_seurat_raw)){
  RomanovDev10x_seurat_raw@meta.data = RomanovDev10x_seurat_meta_final
  saveRDS(RomanovDev10x_seurat_raw,paste0(RomanovDev10x_raw_data_path,"RomanovDev10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}
