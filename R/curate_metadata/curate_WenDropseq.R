##########
### Info
##########

# GSE117295

# https://www.ncbi.nlm.nih.gov/pubmed/32066983

# no metadata found

# see also Wen10x!

## sample metadata via geo
# gds <- GEOquery::getGEO("GSE117295")
# wenDropseq_phenodata=gds$GSE117295_series_matrix.txt.gz@phenoData@data
# data.table::fwrite(wenDropseq_phenodata %>% dplyr::select(zeitgeber = title, GEO_ID = geo_accession),"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE117295_wenDropseq_sample_metadata.tsv",sep="\t")

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
wenDropseq_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/wenDropseq/"
wenDropseq_seurat_raw = readRDS(paste0(wenDropseq_raw_data_path,"wenDropseq_seurat_raw.rds"))
wenDropseq_seurat_raw_meta = wenDropseq_seurat_raw@meta.data

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)
table(hypoMap_neurons_v1_curated_metadata$Diet[hypoMap_neurons_v1_curated_metadata$Dataset=="Wen_10X"])

# sample metadata from geo
wenDropseq_author_sample_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE117295_wenDropseq_sample_metadata.tsv",data.table = FALSE)


##########
### Curate metadata
##########

# add zeitgeber
wenDropseq_seurat_meta = dplyr::left_join(wenDropseq_seurat_raw_meta,wenDropseq_author_sample_metadata,by=c("Sample Name"="GEO_ID"))
wenDropseq_seurat_meta$Author_Condition = wenDropseq_seurat_meta$zeitgeber

# set others
wenDropseq_seurat_meta$Dataset = "wenDropseq"
wenDropseq_seurat_meta$Sample_ID = paste0(wenDropseq_seurat_meta$Dataset,"_",wenDropseq_seurat_meta$Run_ID)
wenDropseq_seurat_meta$Technology ="Dropseq"
wenDropseq_seurat_meta$Author_Region = "Suprachiasmatic nucleus"
wenDropseq_seurat_meta$Diet =  "Normal chow"
#wenDropseq_seurat_meta$Author_Condition = wenDropseq_seurat_meta$sampling_enviroment
wenDropseq_seurat_meta$Age = "6+ weeks"
wenDropseq_seurat_meta$Pooled = "yes"
wenDropseq_seurat_meta$Author_Class = NA
wenDropseq_seurat_meta$Author_CellType = NA
wenDropseq_seurat_meta$Sex = "M"

# infer sex
wenDropseq_seurat_meta$inferred_sex = scUtils::infer_sex(wenDropseq_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

# author ex set to no
wenDropseq_seurat_meta$Author_Exclude = NA


##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = wenDropseq_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Cell_ID,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(wenDropseq_raw_data_path,"wenDropseq_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
wenDropseq_seurat_meta_final = wenDropseq_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=strain,Diet,Pooled,Age,Author_Condition,Author_Region, inferred_sex, Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(wenDropseq_seurat_meta_final) = wenDropseq_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(wenDropseq_seurat_meta_final) == ncol(wenDropseq_seurat_raw)){
  wenDropseq_seurat_raw@meta.data = wenDropseq_seurat_meta_final
  saveRDS(wenDropseq_seurat_raw,paste0(wenDropseq_raw_data_path,"wenDropseq_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}


