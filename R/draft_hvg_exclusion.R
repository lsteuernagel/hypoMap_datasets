
# Kim et al . https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7534821/#S15
# - removing sex-specific genes (Ddx3y, Eif2s3y, Uty, Kdm5d, Xist, Tsix),
# -immediate early genes (e.g., Fos, Fosl2, Junb, Egr1, Arc, Homer1; 139 genes in total from (Wu et al., 2017)),
# - 30 retro-virus-induced genes (e.g., B2m, Bst2, Oasl2, Ifit1)
# - and 1,000 noise-sensitive genes (high abundance genes sensitive to technical noise; see also Table S2).

# - additionally I remove some more rpl and cox genes
# - additionally I compile a list of genes that are expressed only in 4 or less of the original datasets in hypoMap v1 and exclude those as well

library(magrittr)
library(Seurat)

# example seurat:
raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/ChenDropseq/"
seurat_raw_name = "ChenDropseq_seurat_raw.rds"
seurat_raw = readRDS(paste0(raw_data_path,seurat_raw_name))

# load genes from kim et al: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7534821/#S15
hvg_exclude_kim_et_al = readxl::read_excel("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/NIHMS1629366-supplement-Table_S2.xlsx",sheet = 1)

# make a grep for various problematic genes
ignore_genes_regex="Xist|Uty|Utx|Fos|Fosb|Jun|Junb|Jund|Malat1|Heph|Rpl|Rps|Cox|rpl|rps|mt-" # which genes to skip in feature selection
regex_based_genes = rownames(seurat_raw)[grepl(ignore_genes_regex,rownames(seurat_raw))]

# make a shortlist from kim et al and above genes
hvgs_exclude_small = c(hvg_exclude_kim_et_al$`Sex-specific`,hvg_exclude_kim_et_al$IEGs,hvg_exclude_kim_et_al$`retro-virus-induced`) %>% na.omit() %>% as.character()
hvgs_exclude_small = hvgs_exclude_small[hvgs_exclude_small %in% rownames(seurat_raw)]
# make longer list with ethri 1000 genes
hvgs_exclude_long= unique(c(hvgs_exclude_small,hvg_exclude_kim_et_al$`noise-sensitive`),regex_based_genes)


## load hypothalamus and check genes that are very specific for one dataset
hypothalamus_neurons_map = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/hypothalamus_neurons_map.rds")
gene_expr_dataset = hypothalamus_neurons_map@assays$RNA@counts
gene_sums = Matrix::rowSums(gene_expr_dataset)
gene_expr_dataset = gene_expr_dataset[gene_sums > 50 ,]
gene_expr_dataset[gene_expr_dataset!=0] = 1
gene_expr_dataset_dense = as.matrix(gene_expr_dataset)

dataset_factor = hypothalamus_neurons_map@meta.data$Dataset
per_dataset_occ = data.frame(t(rowsum(t(gene_expr_dataset_dense),dataset_factor)))
per_dataset_occ = data.frame(t(apply(per_dataset_occ,1,function(x){return(round((x/sum(x)),5))})))

thresh = 0.0001
per_dataset_occ_binary = as.matrix(per_dataset_occ)
per_dataset_occ_binary[per_dataset_occ_binary>thresh] = 1
per_dataset_occ_binary[per_dataset_occ_binary<1] = 0
per_dataset_occ_min = data.frame(occ=rowSums(per_dataset_occ_binary),gene = rownames(per_dataset_occ_binary))

# get a list
exclude_genes_occurence = per_dataset_occ_min$gene[per_dataset_occ_min$occ<=5]

#make an extra long list
hvgs_exclude_long_extra = unique(c(hvgs_exclude_long,exclude_genes_occurence))

# make lists
features_exclude_list = list(
  hvgs_exclude_small = hvgs_exclude_small,
  hvgs_exclude_long = hvgs_exclude_long,
  hvgs_exclude_long_extra = hvgs_exclude_long_extra
)
scUtils::writeList_to_JSON(features_exclude_list,filename = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/features_exclude_list.json")

##
features_exclude_list=jsonlite::read_json("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/features_exclude_list.json")
features_exclude_list=sapply(features_exclude_list,unlist)

