##########
### Info
##########

# in house dataset


##########
### Load preprocessed object
##########

# Because this datset was already processed with the same pipeline as the others, but does have contain SRA numbers (at the time of processing)
# we just load a preprocseed version and use this as the 'raw' object. this also ensure backwards compatibility because we use the same cells as in earlier HypoMap paper versions

library(dplyr)
library(Seurat)
library(scUtils)

## temp: load original object and save metadata into tables folder:
# snuc_hypo_master = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/yeo_data/hypothalamus_nucSeq/snuc_hypo_master_211101.RDS")
# Dowsett10xnuc_seurat_author_meta = snuc_hypo_master@meta.data
# Dowsett10xnuc_seurat_author_meta$Author_Sample = stringr::str_extract(Dowsett10xnuc_seurat_author_meta$Sample,"[a-zA-Z]+[0-9]+")
# Dowsett10xnuc_seurat_author_meta$Author_Sample = gsub("fast","fasted",Dowsett10xnuc_seurat_author_meta$Author_Sample)
# Dowsett10xnuc_seurat_author_meta$Cell_ID = paste0(gsub("_[0-9]+","",rownames(Dowsett10xnuc_seurat_author_meta)),"_",Dowsett10xnuc_seurat_author_meta$Author_Sample,"_","Dowsett10xnuc")
# data.table::fwrite(Dowsett10xnuc_seurat_author_meta,"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/Dowsett10xnuc_metadata.tsv",sep="\t")

## read raw file and extract metadata
Dowsett10xnuc_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Dowsett10xnuc/"
Dowsett10xnuc_seurat_raw = readRDS(paste0(Dowsett10xnuc_raw_data_path,"Dowsett10xnuc_seurat_raw.rds"))
Dowsett10xnuc_seurat_raw_meta = Dowsett10xnuc_seurat_raw@meta.data

## load original metadata
Dowsett10xnuc_seurat_author_meta = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/Dowsett10xnuc_metadata.tsv",data.table = F)

##########
### Curate author metadata
##########

### set other columns
Dowsett10xnuc_seurat_author_meta$Dataset = "Dowsett10xnuc"
Dowsett10xnuc_seurat_author_meta$Technology ="10xv3"
Dowsett10xnuc_seurat_author_meta$Strain ="C57BL/6J"
Dowsett10xnuc_seurat_author_meta$Pooled = "Yes"
Dowsett10xnuc_seurat_author_meta$Age = "6+ weeks"
Dowsett10xnuc_seurat_author_meta$Author_CellType = Dowsett10xnuc_seurat_author_meta$Cluster_IDs
Dowsett10xnuc_seurat_author_meta$Author_Class = NA
# ventrolateral subdivision of the ventromedial hypothalamus (VMHvl)
Dowsett10xnuc_seurat_author_meta$Author_Region = "Hypothalamus"

# curate some sample ids
Dowsett10xnuc_seurat_author_meta$Sample_ID = paste0("Dowsett10xnuc_",Dowsett10xnuc_seurat_author_meta$Author_Sample)

# mark cells that authors excluded (or did not annotate)
Dowsett10xnuc_seurat_author_meta$Author_Exclude = "no"

## diet:
Dowsett10xnuc_seurat_author_meta$Diet = Dowsett10xnuc_seurat_author_meta$nutr.cond
Dowsett10xnuc_seurat_author_meta$Diet[Dowsett10xnuc_seurat_author_meta$Diet=="adlib"]="Normal chow"
Dowsett10xnuc_seurat_author_meta$Diet[Dowsett10xnuc_seurat_author_meta$Diet=="fast"]="Fasted"


##########
### subset raw to author
##########

Dowsett10xnuc_seurat_meta = dplyr::left_join(Dowsett10xnuc_seurat_raw_meta,Dowsett10xnuc_seurat_author_meta %>% dplyr::select(-orig.ident,-nCount_RNA,-nFeature_RNA),by=c("Cell_ID"))

##########
### Curate other metadata
##########

# infer sex
Dowsett10xnuc_seurat_meta$inferred_sex = scUtils::infer_sex(Dowsett10xnuc_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)
Dowsett10xnuc_seurat_meta$Sex = "M"

Dowsett10xnuc_seurat_meta$Author_Exclude[is.na(Dowsett10xnuc_seurat_meta$Author_Exclude)] = "yes"

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Dowsett10xnuc_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Author_Class,-nCount_RNA,-percent.mt,-Author_Exclude,
                                                                                                              -nFeature_RNA,-Barcode,-sum,-detected,-percent_top_50,-percent_top_100,-percent_top_200,-percent_top_500,-subsets_Mito_detected,-subsets_Mito_percent,-total,-scDblFinder.weighted,-scDblFinder.ratio,-scDblFinder.score,
                                                                                                              -scDblFinder.class,-cxds_score,-cxds_call,-bcds_score,-bcds_call,-hybrid_score,-hybrid_call,-nCount_SCT,-nFeature_SCT,-integrated_snn_res.0.8,-seurat_clusters,-integrated_snn_res.1.5,-Cluster_IDs)
data.table::fwrite(per_sample_summary,file = paste0(Dowsett10xnuc_raw_data_path,"Dowsett10xnuc_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Dowsett10xnuc_seurat_meta_final = Dowsett10xnuc_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,Sample_ID, Technology,Author_Region,Strain ,Diet,Pooled,Age,inferred_sex,Sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Exclude,Author_Class,Author_CellType  ) %>%
  as.data.frame()
rownames(Dowsett10xnuc_seurat_meta_final) = Dowsett10xnuc_seurat_meta_final$Cell_ID

# overwrite metadata
if(nrow(Dowsett10xnuc_seurat_meta_final) == ncol(Dowsett10xnuc_seurat_raw)){
  Dowsett10xnuc_seurat_raw@meta.data = Dowsett10xnuc_seurat_meta_final
  #saveRDS(Dowsett10xnuc_seurat_raw,paste0(Dowsett10xnuc_raw_data_path,"Dowsett10xnuc_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}

######## SUBSET ! ########
Dowsett10xnuc_seurat_raw = subset(Dowsett10xnuc_seurat_raw,subset = Author_Exclude == "no")

# save
saveRDS(Dowsett10xnuc_seurat_raw,paste0(Dowsett10xnuc_raw_data_path,"Dowsett10xnuc_seurat_raw.rds"))



