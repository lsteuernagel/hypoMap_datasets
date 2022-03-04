##########
### Info
##########

# GSE93374
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5323293/

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
CampbellDropseq_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/CampbellDropseq/"
CampbellDropseq_seurat_raw = readRDS(paste0(CampbellDropseq_raw_data_path,"CampbellDropseq_seurat_raw.rds"))
CampbellDropseq_seurat_raw_meta = CampbellDropseq_seurat_raw@meta.data

# load author metadata
campbell_metadata= data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE93374_cell_metadata.txt",data.table = F,header = TRUE,select=1:11)
campbell_mapping_all_a = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE93374_cell_metadata.txt",data.table = F,header = TRUE,select=13)
campbell_mapping_all_a = campbell_mapping_all_a%>% filter(`All Cell Clusters`!="") %>% tidyr::separate(`All Cell Clusters`, into=c("all_cluster_id","all_cluster"),sep="\\.")
campbell_mapping_sub_s = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE93374_cell_metadata.txt",data.table = F,header = TRUE,select=14)
campbell_mapping_sub_s = campbell_mapping_sub_s%>% filter(`All Cell Subclusters`!="") %>% tidyr::separate(`All Cell Subclusters`, into=c("sub_cluster_id","sub_cluster"),sep="\\.")
campbell_mapping_neurons_n = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/GSE93374_cell_metadata.txt",data.table = F,header = TRUE,select=15)
campbell_mapping_neurons_n = campbell_mapping_neurons_n%>% filter(`Neuron Subclusters`!="") %>% tidyr::separate(`Neuron Subclusters`, into=c("neuron_cluster_id","neuron_cluster"),sep="\\.")

campbell_metadata = left_join(campbell_metadata,campbell_mapping_all_a,by=c("7.clust_all"="all_cluster_id"))
campbell_metadata = left_join(campbell_metadata,campbell_mapping_sub_s,by=c("9.clust_all_micro"="sub_cluster_id"))
CampbellDropseq_seurat_author_meta = left_join(campbell_metadata,campbell_mapping_neurons_n,by=c("10.clust_neurons"="neuron_cluster_id"))
# add a sample ID
CampbellDropseq_seurat_author_meta = CampbellDropseq_seurat_author_meta %>% dplyr::mutate(Author_Sample = paste0(`3.batches`,"_",`5.Diet`,"_",`4.sex`))
# add barcode as exra column
CampbellDropseq_seurat_raw_meta$Barcode = stringr::str_extract(CampbellDropseq_seurat_raw_meta$Cell_ID,pattern = "[ACGT]{5,}\\-?[0-9]?")
CampbellDropseq_seurat_author_meta$Barcode = stringr::str_extract(CampbellDropseq_seurat_author_meta$`1.ID`,pattern = "[ACGT]{5,}\\-?[0-9]?")

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)

##########
### set general metadata
##########

### set other columns
CampbellDropseq_seurat_author_meta$Dataset = "CampbellDropseq"
CampbellDropseq_seurat_author_meta$Technology ="Dropseq"
CampbellDropseq_seurat_author_meta$Strain = "C57BL6/J"
CampbellDropseq_seurat_author_meta$Pooled ="Yes"
CampbellDropseq_seurat_author_meta$Author_Region =  "Arcuate Nucleus/Median Eminence"

## check v1 metadata for diet and age:
matched_Sample_names_v1 = scUtils::match_sample_names(table_A = CampbellDropseq_seurat_author_meta,table_B = hypoMap_neurons_v1_curated_metadata,sample_col_A = "Author_Sample",sample_col_B = "Diet",barcode_col = "Barcode",min_pct = 0.1)
CampbellDropseq_seurat_author_meta$Diet = CampbellDropseq_seurat_author_meta$`5.Diet`
CampbellDropseq_seurat_author_meta$Diet[CampbellDropseq_seurat_author_meta$Diet=="Chow"] = "Normal chow"
CampbellDropseq_seurat_author_meta$Diet[CampbellDropseq_seurat_author_meta$Diet=="Ch10"] = "Normal chow"
CampbellDropseq_seurat_author_meta$Diet[CampbellDropseq_seurat_author_meta$Diet=="Fast"] = "Fasted"
CampbellDropseq_seurat_author_meta$Diet[CampbellDropseq_seurat_author_meta$Diet=="Refed"] = "Refed"
CampbellDropseq_seurat_author_meta$Diet[CampbellDropseq_seurat_author_meta$Diet=="HFD"] = "HFD_1week"

## check v1 metadata for age:
matched_Sample_names_v1 = scUtils::match_sample_names(table_A = CampbellDropseq_seurat_author_meta,table_B = hypoMap_neurons_v1_curated_metadata,sample_col_A = "Author_Sample",sample_col_B = "Age",barcode_col = "Barcode",min_pct = 0.1)
CampbellDropseq_seurat_author_meta$Age= "6+ weeks"
CampbellDropseq_seurat_author_meta$Age[CampbellDropseq_seurat_author_meta$Author_Sample %in% c("b1_Chow_M","b2_Chow_M","b3_Chow_M")] = "3-6 weeks"

##########
### Clean up cell annotations
##########

# Curate Author Cell types
CampbellDropseq_seurat_author_meta$Author_Class = gsub("[0-9]","",CampbellDropseq_seurat_author_meta$all_cluster)
CampbellDropseq_seurat_author_meta$Author_Class[is.na(CampbellDropseq_seurat_author_meta$Author_Class)] = "NA"
CampbellDropseq_seurat_author_meta$neuron_cluster[is.na(CampbellDropseq_seurat_author_meta$neuron_cluster)] = "NA"
CampbellDropseq_seurat_author_meta$Author_CellType = CampbellDropseq_seurat_author_meta$sub_cluster
CampbellDropseq_seurat_author_meta$Author_CellType[CampbellDropseq_seurat_author_meta$Author_Class =="Neurons"] = CampbellDropseq_seurat_author_meta$neuron_cluster[CampbellDropseq_seurat_author_meta$Author_Class =="Neurons"]
# rework cell classes to one schema
table(CampbellDropseq_seurat_author_meta$Author_Class)
CampbellDropseq_seurat_author_meta$Author_Class[CampbellDropseq_seurat_author_meta$Author_Class=="Neuron"]="Neurons"
CampbellDropseq_seurat_author_meta$Author_Class[CampbellDropseq_seurat_author_meta$Author_Class=="Astrocyte"]="Astrocytes"
CampbellDropseq_seurat_author_meta$Author_Class[grepl("Endothelial",CampbellDropseq_seurat_author_meta$Author_Class)]="Endothelial"
CampbellDropseq_seurat_author_meta$Author_Class[grepl("Mural",CampbellDropseq_seurat_author_meta$Author_Class)]="Mural"
CampbellDropseq_seurat_author_meta$Author_Class[grepl("Oligo",CampbellDropseq_seurat_author_meta$Author_Class)]="Oligodendrocytes"
CampbellDropseq_seurat_author_meta$Author_Class[grepl("Ependym",CampbellDropseq_seurat_author_meta$Author_Class)]="Ependymal"
CampbellDropseq_seurat_author_meta$Author_Class[CampbellDropseq_seurat_author_meta$Author_Class=="PVMMicro"]="Immune"
CampbellDropseq_seurat_author_meta$Author_Class[CampbellDropseq_seurat_author_meta$Author_Class=="Tanycyte"]="Tanycytes"
table(CampbellDropseq_seurat_author_meta$Author_Class)

##########
### Match data on sample level
##########

# match samples and join
matched_Sample_names = scUtils::match_sample_names(table_A = CampbellDropseq_seurat_author_meta,table_B = CampbellDropseq_seurat_raw_meta,sample_col_A = "Author_Sample",sample_col_B = "Run_ID",barcode_col = "Barcode",min_pct = 0.1)
# add RUN_IDs via matched samples
CampbellDropseq_seurat_author_meta = dplyr::left_join(CampbellDropseq_seurat_author_meta,matched_Sample_names,by=c("Author_Sample"="Samples_A")) %>% dplyr::rename(Run_ID = Samples_B)
# reduce to per sample
CampbellDropseq_seurat_author_meta_perSample = CampbellDropseq_seurat_author_meta %>% dplyr::select(Dataset,Technology,Strain,Pooled,Author_Region,Diet,Age,Run_ID,Author_Sample,Author_Batch = `3.batches`) %>%
  dplyr::distinct(.keep_all = TRUE)
# add to metadata
CampbellDropseq_seurat_meta = dplyr::left_join(CampbellDropseq_seurat_raw_meta,CampbellDropseq_seurat_author_meta_perSample,by="Run_ID")


##########
### Match data on cell level
##########

# matching does not work very well --> probably due to excluded cells
# I solve this by joining just on barcode --> but then we need to take care of duplicated barcodes
duplicated_barcodes = CampbellDropseq_seurat_author_meta$Barcode[duplicated(CampbellDropseq_seurat_author_meta$Barcode)]
message("length duplicated barcodes: ",length(duplicated_barcodes))

# join non-duplicated barcodes
CampbellDropseq_seurat_meta_all = dplyr::left_join(CampbellDropseq_seurat_meta %>% dplyr::filter(! Barcode %in% duplicated_barcodes) ,CampbellDropseq_seurat_author_meta %>% dplyr::select(Author_Class,Author_CellType,Run_ID,Barcode) %>% dplyr::filter(! Barcode %in% duplicated_barcodes) %>% dplyr::select(-Run_ID),by=c("Barcode"="Barcode"))
# join duplicated barcodes
CampbellDropseq_seurat_meta_dup = dplyr::left_join(CampbellDropseq_seurat_meta %>% dplyr::filter( Barcode %in% duplicated_barcodes) ,CampbellDropseq_seurat_author_meta %>% dplyr::select(Author_Class,Author_CellType,Run_ID,Barcode) %>% dplyr::filter(Barcode %in% duplicated_barcodes),by=c("Barcode"="Barcode","Run_ID"="Run_ID"))
# rbind:
CampbellDropseq_seurat_meta = dplyr::bind_rows(CampbellDropseq_seurat_meta_all,CampbellDropseq_seurat_meta_dup)
# reorder
rownames(CampbellDropseq_seurat_meta) = CampbellDropseq_seurat_meta$Cell_ID
CampbellDropseq_seurat_meta = CampbellDropseq_seurat_meta[match(rownames(CampbellDropseq_seurat_raw_meta),rownames(CampbellDropseq_seurat_meta)),]
message("nrow CampbellDropseq_seurat_meta: ",nrow(CampbellDropseq_seurat_meta)," , should be: ",nrow(CampbellDropseq_seurat_raw_meta))

##########
### Curate remainder of metadata
##########

# infer sex
CampbellDropseq_seurat_meta$inferred_sex = scUtils::infer_sex(CampbellDropseq_seurat_raw,sample_column="Run_ID",min_xist_female = 0.7,max_xist_male = 0.1)

# curate some sample ids
CampbellDropseq_seurat_meta$Sample_ID = paste0("CampbellDropseq_",CampbellDropseq_seurat_meta$Author_Sample)

# mark cells that authors excluded (or did not annotate)
CampbellDropseq_seurat_meta$Author_Exclude = "no"
CampbellDropseq_seurat_meta$Author_Exclude[is.na(CampbellDropseq_seurat_meta$Author_Class)] = "yes"

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = CampbellDropseq_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA)
data.table::fwrite(per_sample_summary,file = paste0(CampbellDropseq_raw_data_path,"CampbellDropseq_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
CampbellDropseq_seurat_meta_final = CampbellDropseq_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,Author_Batch,SRA_ID = Run_ID, Sample_ID, GEO_ID = `Sample Name`,Technology,Strain=Strain.x,Diet,Pooled,Age,Author_Region,inferred_sex, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Exclude,Author_Class,Author_CellType  ) %>%
  as.data.frame()
rownames(CampbellDropseq_seurat_meta_final) = CampbellDropseq_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(CampbellDropseq_seurat_meta_final) == ncol(CampbellDropseq_seurat_raw)){
  CampbellDropseq_seurat_raw@meta.data = CampbellDropseq_seurat_meta_final
  saveRDS(CampbellDropseq_seurat_raw,paste0(CampbellDropseq_raw_data_path,"CampbellDropseq_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}

# CampbellDropseq_seurat_PROCESS = scUtils::seurat_recipe(CampbellDropseq_seurat_raw,calcUMAP = TRUE,nfeatures_vst=1000)
# DimPlot(CampbellDropseq_seurat_PROCESS,group.by = "Author_CellType")
