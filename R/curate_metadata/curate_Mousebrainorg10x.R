

##########
### Info
##########

# SRP135960
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6086934/
# http://mousebrain.org/adolescent/downloads.html

# with this curation we retain much more cells than in the original loom file from mousebrain:
# mapping onto the hypomap v1 shows that this makes sense and man of these cells indeed are true hypothalamic cell types!

##########
### Load data
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
Mousebrainorg10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Mousebrainorg10x/"
Mousebrainorg10x_seurat_raw = readRDS(paste0(Mousebrainorg10x_raw_data_path,"Mousebrainorg10x_seurat_raw.rds"))
Mousebrainorg10x_seurat_raw_meta = Mousebrainorg10x_seurat_raw@meta.data

## read sample metdata ---> nearly all this information is already part of the SRP metadata so I won't join it
Mousebrainorg10x_seurat_cell_meta= data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/mousebrain_org_cellmetdata.tsv",data.table = F,header = TRUE)
# read cell metdata --> see below for how to extract this from the loom files at: http://mousebrain.org/adolescent/downloads.html
Mousebrainorg10x_seurat_sample_meta= data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/mousebrain_org_sample_metadata.txt",data.table = F,header = TRUE)

# add barcode as extra column
Mousebrainorg10x_seurat_raw_meta$Barcode = stringr::str_extract(Mousebrainorg10x_seurat_raw_meta$Cell_ID,pattern = "[ACGT]{5,}") # leave out -1
Mousebrainorg10x_seurat_cell_meta$Barcode = stringr::str_extract(Mousebrainorg10x_seurat_cell_meta$CellID,pattern = "[ACGT]{5,}")
# for some reason the mousebrain cell ids have the first two bases from the barcode missing --> cut these from the srp data as well to allow matching:
Mousebrainorg10x_seurat_raw_meta$Barcode = substring(Mousebrainorg10x_seurat_raw_meta$Barcode,first = 3)

# load hypoMpa v1 annotations
hypoMap_neurons_v1_curated_metadata = data.table::fread("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/hypoMap_neurons_v1_curated_metadata.tsv",data.table = FALSE)
table(hypoMap_neurons_v1_curated_metadata$Strain[hypoMap_neurons_v1_curated_metadata$Dataset=="Mousebrainorg"])

##########
### Add metadata sample level
##########

# not necessary ..
# Mousebrainorg10x_seurat_meta = dplyr::left_join(Mousebrainorg10x_seurat_raw_meta,Mousebrainorg10x_seurat_sample_meta,by=c("Sample Name"="SampleID"))

##########
### Add metadata cell level
##########

# match samples and join ---> not necessary: SampleID is available in both
# matched_Sample_names = scUtils::match_sample_names(table_A = Mousebrainorg10x_seurat_cell_meta,table_B = Mousebrainorg10x_seurat_raw_meta,sample_col_A = "SampleID",sample_col_B = "Run_ID",barcode_col = "Barcode")

# there will be duplicates:
Mousebrainorg10x_seurat_raw_meta$pasted_barcode = paste0(Mousebrainorg10x_seurat_raw_meta$Barcode,"_",Mousebrainorg10x_seurat_raw_meta$`Sample Name`)
Mousebrainorg10x_seurat_cell_meta$pasted_barcode = paste0(Mousebrainorg10x_seurat_cell_meta$Barcode,"_",Mousebrainorg10x_seurat_cell_meta$SampleID)
duplicated_barcodes_samples = unique(c(Mousebrainorg10x_seurat_cell_meta_all$pasted_barcode[duplicated(Mousebrainorg10x_seurat_cell_meta_all$pasted_barcode)],
                                     Mousebrainorg10x_seurat_raw_meta$pasted_barcode[duplicated(Mousebrainorg10x_seurat_raw_meta$pasted_barcode)]))
message("length duplicated barcodes: ",length(duplicated_barcodes_samples))

# join non-duplicated barcodes
Mousebrainorg10x_seurat_cell_meta_all = dplyr::left_join(Mousebrainorg10x_seurat_raw_meta %>% dplyr::filter(! pasted_barcode %in% duplicated_barcodes_samples) ,Mousebrainorg10x_seurat_cell_meta %>% dplyr::filter(! pasted_barcode %in% duplicated_barcodes_samples),by=c("pasted_barcode"="pasted_barcode"))
# cannot join duplicated barcodes without knowing the ids -- just use srr metadata
Mousebrainorg10x_seurat_cell_meta_dup = Mousebrainorg10x_seurat_raw_meta %>% dplyr::filter( pasted_barcode %in% duplicated_barcodes_samples)
# rbind:
Mousebrainorg10x_seurat_meta = dplyr::bind_rows(Mousebrainorg10x_seurat_cell_meta_all,Mousebrainorg10x_seurat_cell_meta_dup)
# reorder
rownames(Mousebrainorg10x_seurat_meta) = Mousebrainorg10x_seurat_meta$Cell_ID
Mousebrainorg10x_seurat_meta = Mousebrainorg10x_seurat_meta[match(rownames(Mousebrainorg10x_seurat_raw_meta),rownames(Mousebrainorg10x_seurat_meta)),]
message("nrow Mousebrainorg10x_seurat_meta: ",nrow(Mousebrainorg10x_seurat_meta)," , should be: ",nrow(Mousebrainorg10x_seurat_raw_meta))

# set author exclude:
Mousebrainorg10x_seurat_meta$Author_Exclude = "no"
Mousebrainorg10x_seurat_meta$Author_Exclude[is.na(Mousebrainorg10x_seurat_meta$Class)] = "yes"

##########
### Curate metadata
##########
# cell class
Mousebrainorg10x_seurat_meta$Author_Class = Mousebrainorg10x_seurat_meta$Class
Mousebrainorg10x_seurat_meta$Author_Class[Mousebrainorg10x_seurat_meta$Author_Class=="Neuron"]="Neurons"
Mousebrainorg10x_seurat_meta$Author_Class[Mousebrainorg10x_seurat_meta$Author_Class=="Astrocyte"]="Astrocytes"
Mousebrainorg10x_seurat_meta$Author_Class[grepl("Vascular",Mousebrainorg10x_seurat_meta$Author_Class)]="Endothelial"
Mousebrainorg10x_seurat_meta$Author_Class[grepl("Oligo",Mousebrainorg10x_seurat_meta$Author_Class)]="Oligodendrocytes"
Mousebrainorg10x_seurat_meta$Author_Class[grepl("Ependym",Mousebrainorg10x_seurat_meta$Author_Class)]="Ependymal"
# cell type
Mousebrainorg10x_seurat_meta$Author_CellType = paste0(Mousebrainorg10x_seurat_meta$ClusterName,"_",Mousebrainorg10x_seurat_meta$Taxonomy_group)
less_than_5 = names(table(Mousebrainorg10x_seurat_meta$Author_CellType))[which(table(Mousebrainorg10x_seurat_meta$Author_CellType)<5)]
Mousebrainorg10x_seurat_meta$Author_CellType[Mousebrainorg10x_seurat_meta$Author_CellType  %in% less_than_5] = NA
Mousebrainorg10x_seurat_meta$Author_CellType[Mousebrainorg10x_seurat_meta$Author_CellType == "NA_NA"] = NA
# age
Mousebrainorg10x_seurat_meta$Age[Mousebrainorg10x_seurat_meta$Age != "p19"] = "3-6 weeks"
Mousebrainorg10x_seurat_meta$Age[Mousebrainorg10x_seurat_meta$Age == "p19"] = "0-3 weeks"
# sex
table(Mousebrainorg10x_seurat_meta$sex)
Mousebrainorg10x_seurat_meta$Sex = NA
Mousebrainorg10x_seurat_meta$Sex[Mousebrainorg10x_seurat_meta$sex == "female"] = "F"
Mousebrainorg10x_seurat_meta$Sex[Mousebrainorg10x_seurat_meta$sex == "male"] = "M"

### set other columns
Mousebrainorg10x_seurat_meta$Dataset = "Mousebrainorg10x"
Mousebrainorg10x_seurat_meta$Sample_ID = Mousebrainorg10x_seurat_meta$`Sample Name` # use their original sample IDs
Mousebrainorg10x_seurat_meta$Diet ="Normal chow"
Mousebrainorg10x_seurat_meta$Author_Region = "Hypothalamus"
Mousebrainorg10x_seurat_meta$Diet ="Normal chow"
Mousebrainorg10x_seurat_meta$Pooled ="no"
Mousebrainorg10x_seurat_meta$Strain = Mousebrainorg10x_seurat_meta$strain
Mousebrainorg10x_seurat_meta$Author_Condition = Mousebrainorg10x_seurat_meta$Isolate # use isolate (mouse donor as 'condition')
Mousebrainorg10x_seurat_meta$GEO_ID = NA
Mousebrainorg10x_seurat_meta$Run10x = Mousebrainorg10x_seurat_meta$cell_index # add cell index here
Mousebrainorg10x_seurat_meta$Technology = "10xv2"

# inferred sex
Mousebrainorg10x_seurat_meta$inferred_sex = scUtils::infer_sex(Mousebrainorg10x_seurat_raw,sample_column="Run_ID",min_xist_female = 0.5,max_xist_male = 0.1)

##########
### Save updated version
##########

## save per sample summary of metadata
per_sample_summary = Mousebrainorg10x_seurat_meta %>% dplyr::distinct(Run_ID,.keep_all = TRUE) %>% dplyr::select(-Author_CellType,-Author_Exclude,-Author_Class,-nCount_RNA,-percent.mt,-nFeature_RNA,-Barcode.x,-pasted_barcode,-CellID,-Class,-ClusterName,-Clusters,-Description,-Developmental_compartment,-LeafOrder,-MitoRiboRatio,-Neurotransmitter,-OriginalClusters,-Probable_location,-Region,-Subclass, -TaxonomyRank1, -TaxonomyRank2,-TaxonomyRank3,-TaxonomyRank4,-TaxonomySymbol,-Taxonomy_group,-Barcode.y,-Barcode)#,-Clusters,-Description,-Developmental_compartment,-LeafOrder,-MitoRiboRatio,-Neurotransmitter,-OriginalClusters,-Probable_location,-Region,-SampleID,-Subclass, -TaxonomyRank1, -TaxonomyRank2,-TaxonomyRank3,-TaxonomyRank4,-TaxonomySymbol,-Taxonomy_group,-_NGenes,-_tSNE1,-_tSNE2,-Barcode.y,-Barcode)
data.table::fwrite(per_sample_summary,file = paste0(Mousebrainorg10x_raw_data_path,"Mousebrainorg10x_per_sample_summary.tsv"),sep = "\t")

## make final sorting and selection of columns
Mousebrainorg10x_seurat_meta_final = Mousebrainorg10x_seurat_meta %>%
  dplyr::select(Cell_ID,Dataset,SRA_ID = Run_ID, Sample_ID, GEO_ID,Technology,Run10x,Strain=strain,Diet,Pooled,Age,Author_Region, inferred_sex, Sex,Author_Condition, nCount_RNA, nFeature_RNA,percent_mt = percent.mt,Author_Exclude,Author_Class,Author_CellType) %>%
  as.data.frame()
rownames(Mousebrainorg10x_seurat_meta_final) = Mousebrainorg10x_seurat_meta_final$Cell_ID

## overwrite raw rds object
if(nrow(Mousebrainorg10x_seurat_meta_final) == ncol(Mousebrainorg10x_seurat_raw)){
  Mousebrainorg10x_seurat_raw@meta.data = Mousebrainorg10x_seurat_meta_final
  saveRDS(Mousebrainorg10x_seurat_raw,paste0(Mousebrainorg10x_raw_data_path,"Mousebrainorg10x_seurat_raw.rds"))
}else{
  message("metdata does not contain the same number of cells as raw seurat !")
}

######## Code for metadata extracation from older loom file:
#
# # https://satijalab.org/loomR/loomR_tutorial.html
# # Connect to the loom file in read/write mode
# library(loomR)
#
# loom_level5_mousebrain <- connect(filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/public_scseq/mousebrain_org/all_level5/","l5_all.loom"), mode = "r+")
#
# # Viewing a dataset in the 'row_attrs' group with S3 $ chaining
# all_col_attributes = names(loom_level5_mousebrain$col.attrs)
# # since get.attribute.df is broken (https://github.com/mojaveazure/loomR/issues/34) I iterate over all columns and manually build the metadata
# columndata_list = list()
# for(i in 1:length(all_col_attributes)){
#   columndata_list[[i]]=loom_level5_mousebrain[[paste0("col_attrs/",all_col_attributes[i])]][]
# }
# metadata_mousebrain_level5 = as.data.frame(do.call(cbind,columndata_list),stringsAsFactors=F)
# colnames(metadata_mousebrain_level5) = all_col_attributes
# # manually select interesting columns
# metadata_relev_mousebrain_level5 = metadata_mousebrain_level5[,c(6,9,41,42,47,48,54,59,61,65,72,82,83,98,101:106,115,122,123)]
# metadata_relev_mousebrain_level5$`_tSNE1`=as.numeric(metadata_relev_mousebrain_level5$`_tSNE1`)
# metadata_relev_mousebrain_level5$`_tSNE2`=as.numeric(metadata_relev_mousebrain_level5$`_tSNE2`)
# # save cell metadata
# data.table::fwrite(metadata_relev_mousebrain_level5,"/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/metadata_tables/mousebrain_org_cellmetdata.tsv",sep="\t")
#
# # read and add sample metadata (optionally)
# # mousebrain_org_sample_metadata = data.table::fread(paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/public_scseq/mousebrain_org/all_level5/","mousebrain_org_sample_metadata.txt"),data.table = F)
# # metadata_relev_mousebrain_level5 = left_join(metadata_relev_mousebrain_level5,mousebrain_org_sample_metadata,by=c("SampleID"="SampleID"))
# # metadata_relev_mousebrain_level5$row_index = 1:nrow(metadata_relev_mousebrain_level5)
# # colnames(metadata_relev_mousebrain_level5)


