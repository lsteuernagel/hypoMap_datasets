## This script creates the json with general parameters --> make otehr jsons or edit this script if other params should be used
# requires some manually decisions which are added here.

param_list = list()

# must be loaded from params:
param_list$data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/"
param_list$processed_suffix = "_seurat_processed"
#param_list$processed_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_processed/"

# for final merged object:
param_list$merged_name = "hypoMap_merged_raw"
param_list$feature_set_sizes = c(750,1000,1250,1500,2000,2500,3000,4000) # TODO

# processing
param_list$n_cores = 50
param_list$id_column = "Cell_ID"
param_list$global_seed = 123456
param_list$minUMI = 1000
param_list$minFeatures = 500
param_list$maxUMI = 1000000
param_list$maxUMI_dynamic = 20
param_list$max_mt = 10
param_list$genes_to_exclude_file = "data/features_exclude_list.json"
param_list$max_pctExclude = 40 # preliminary cluster with more than this percentage of genes_to_exclude (median) will be flagged
param_list$min_cells_sample = 100
param_list$sample_column = "Sample_ID"
param_list$nfeatures_vst_prelim = 1000
param_list$npcs_PCA = 70
param_list$k_param = 30
param_list$resolutions_to_try = c(0.5, 0.75, 1, 1.5, 2:10)
param_list$downsampling_max = 20000 # downsample to maximum of x cells to avoid mem overflow or very long run times during rf training
param_list$max_entropy_batch_detection = 0.8 ## important !
param_list$trees_rf = 20000

# doublet
param_list$pN_fixed = 0.25 # this can be the same globally
param_list$pK_max = 0.1 #  this can be the same globally
param_list$doublet_cluster_tresh = 0.7 # this can be the same globally
param_list$min_avg_cells_cluster = 200 # to avoid accientally overclustering small datasets all cluster should have on average this number of cells, so if target_cluster_number = 50 --> 10000 cells required, else target_cluster_number will be scaled down.

# save
scUtils::writeList_to_JSON(param_list,filename = "data/parameters_pre_processing_v2_1.json")

