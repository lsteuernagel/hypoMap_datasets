##########
### Overview
##########

## we want to do some final clean up
# - remove Doublets
# - remove other Process_Exclude
# - remove Author_Exclude for selected datasets
# - Find out which are the problemetic neuron clusters in hypoMap v1 and check whetehr to remove them
# - Find out what the Tany/astrocyte author class clusters that dock onto the neuron clusters are and remove

##########
### Load
##########
library(Seurat)
## load pre-processed dataset
global_parameters = jsonlite::read_json("data/parameters_pre_processing_v2_1.json")
hypoMap_merged_raw = readRDS(paste0(global_parameters$data_path,"hypoMap_merged_raw.rds"))

## ad column for final removal
hypoMap_merged_raw@meta.data$Final_Exclude = hypoMap_merged_raw@meta.data$Process_Exclude # Exclude all problemetics marked during pre-processing
hypoMap_merged_raw@meta.data$Final_Exclude[hypoMap_merged_raw@meta.data$Doublet == "Doublet"] = "yes"  # Exclude all Doublets marked during pre-processing

# something is wrong with exclude features
features_exclude_list= unlist(jsonlite::read_json("data/features_exclude_list.json"))
hypoMap_merged_raw@meta.data$percent_exclude_features2 <- Matrix::colSums(hypoMap_merged_raw@assays$RNA@counts[features_exclude_list[features_exclude_list %in% rownames(hypoMap_merged_raw)],]) / Matrix::colSums(hypoMap_merged_raw@assays$RNA@counts)

# repeat clustering with high res
hypoMap_merged_raw = FindClusters(hypoMap_merged_raw,resolution = 5)
hypoMap_merged_raw@meta.data$Processing_clusters = hypoMap_merged_raw@meta.data$seurat_clusters

# neuron hypoMap:
hypothalamus_neurons_map = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap/hypoMap_largeFiles/hypothalamus_neurons_map.rds")
DimPlot(hypothalamus_neurons_map,group.by = "K98_pruned",label = TRUE,label.size = 2)+NoLegend()

##########
### visualize
##########

# - remove Doublets
p1= DimPlot(hypoMap_merged_raw,group.by = "seurat_clusters",raster = F,shuffle = TRUE,label=TRUE,repel = TRUE)+NoLegend()
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r

# - remove Doublets
p1= DimPlot(hypoMap_merged_raw,group.by = "Doublet",raster = F,shuffle = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r
# - remove other Process_Exclude
p1= DimPlot(hypoMap_merged_raw,group.by = "Process_Exclude",raster = F,shuffle = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r
# - remove Author_Exclude for selected datasets
p1= DimPlot(hypoMap_merged_raw,group.by = "Author_Exclude",raster = F,shuffle = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r

p1= DimPlot(hypoMap_merged_raw,group.by = "Phase",raster = F,shuffle = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r
# general QC:
p1= FeaturePlot(hypoMap_merged_raw,features =  "S.Score",raster = F,order = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r
p1= FeaturePlot(hypoMap_merged_raw,features =  "G2M.Score",raster = F,order = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r
p1= FeaturePlot(hypoMap_merged_raw,features =  "nCount_RNA",raster = F,order = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r
p1= FeaturePlot(hypoMap_merged_raw,features =  "percent_exclude_features2",raster = F,order = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r
p1= FeaturePlot(hypoMap_merged_raw,features =  "nFeature_RNA",raster = F,order = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r
p1= FeaturePlot(hypoMap_merged_raw,features =  "percent_mt",raster = F,order = TRUE)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r

##########
### Other Doublets
##########

# - mark some clusters that are containing Doublets for full removal
treshold_auc = 0.03
doublet_summary = hypoMap_merged_raw@meta.data %>% dplyr::group_by(Processing_clusters) %>% dplyr::add_count(name="n_cells") %>%
  dplyr::group_by(Doublet,Processing_clusters,n_cells) %>% dplyr::count(name="n_doublets") %>% dplyr::ungroup() %>%
  dplyr::mutate(pct_doublets = n_doublets / n_cells) %>% dplyr::filter(Doublet =="Doublet") %>% dplyr::select(-Doublet)
# which clusters to remove:
suggested_Doublet_clusters = as.character(doublet_summary$Processing_clusters[doublet_summary$pct_doublets > 0.3])
### manual curation:
## don't remove fully: 47, 122, 83, 109,43,130,10, 85
suggested_Doublet_clusters = suggested_Doublet_clusters[!suggested_Doublet_clusters %in% c("47","122","83","109","43","130","10","85")]
# maybe remove fully?: 166
## remove fully: 99, 49, 155, 168, 165, 150, 154, 114, 152, 163,82,
suggested_Doublet_clusters = unique(c(suggested_Doublet_clusters,
                                      c("99","49","155","168","165","150","154","114","152","163","82","113","161","158","161","135","90")))
## make final doublet anno:
hypoMap_merged_raw@meta.data$Final_Doublet = hypoMap_merged_raw@meta.data$Doublet
hypoMap_merged_raw@meta.data$Final_Doublet[hypoMap_merged_raw@meta.data$Processing_clusters %in% suggested_Doublet_clusters] = "Doublet"
## protect: 163!, 118?
hypoMap_merged_raw@meta.data$Final_Doublet[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("163","118")] = "Singlet"

p1= DimPlot(hypoMap_merged_raw,group.by = "Final_Doublet",raster = F,shuffle = TRUE,label=F)
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r


##########
### Problematic neuron cluster from v1
##########

# - Find out which are the problemetic neuron clusters in hypoMap v1 and check whetehr to remove them
# K169-51 and K169-86
# add barcode as exra column
hypoMap_merged_raw@meta.data$Barcode = stringr::str_extract(hypoMap_merged_raw@meta.data$Cell_ID,pattern = "[ACGT]{5,}\\-?[0-9]?")
hypothalamus_neurons_map@meta.data$Barcode = stringr::str_extract(hypothalamus_neurons_map@meta.data$Cell_ID,pattern = "[ACGT]{5,}\\-?[0-9]?")
barcode_list_169_51 = hypothalamus_neurons_map@meta.data$Barcode[hypothalamus_neurons_map@meta.data$K169_pruned == "K169-51"]
barcode_list_169_86 = hypothalamus_neurons_map@meta.data$Barcode[hypothalamus_neurons_map@meta.data$K169_pruned == "K169-86"]
# get all pontetial cells --> some fals positives, but should not be too bad
cells_list_169_51 = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Barcode %in% barcode_list_169_51 &
                                                           !hypoMap_merged_raw@meta.data$Dataset %in% c('Affinati10x' , 'Anderson10x' , 'Dowsett10xnuc' , 'Morris10x', 'RossiDropseq' , 'Rupp10x' )]
cells_list_169_86 = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Barcode %in% barcode_list_169_86 &
                                                           !hypoMap_merged_raw@meta.data$Dataset %in% c('Affinati10x' , 'Anderson10x' , 'Dowsett10xnuc' , 'Morris10x', 'RossiDropseq' , 'Rupp10x' )]

# make final exclusion
hypoMap_merged_raw@meta.data$Final_Exclude = hypoMap_merged_raw@meta.data$Process_Exclude
#hypoMap_merged_raw@meta.data$Final_Exclude[hypoMap_merged_raw@meta.data$Process_Exclude == "yes"] = "yes"
hypoMap_merged_raw@meta.data$Final_Exclude[hypoMap_merged_raw@meta.data$Final_Doublet == "Doublet"] = "yes"
hypoMap_merged_raw@meta.data$Final_Exclude[hypoMap_merged_raw@meta.data$Cell_ID %in% cells_list_169_51] = "yes"
hypoMap_merged_raw@meta.data$Final_Exclude[hypoMap_merged_raw@meta.data$Cell_ID %in% cells_list_169_86] = "yes"

##########
### Problematic Tany/astrocyte
##########

# - Find out what the Tany/astrocyte author class clusters that dock onto the neuron clusters are and remove
# --> Affinati, probably Doublets
# 67 # part of 96 --> umap1 > 0
hypoMap_merged_raw@meta.data$Final_Exclude[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("67","45","96","69") & hypoMap_merged_raw@reductions$umap@cell.embeddings[,1] > 0] = "yes"

# also problematic: 35 and 31 (something strange in the affinati data), 160 (hbbs)
hypoMap_merged_raw@meta.data$Final_Exclude[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("35","31","160") ] = "yes"

# cluster 27:
hypoMap_merged_raw@meta.data$Final_Exclude[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("27") ] = "yes"


##########
### Add author exclude from selected
##########

# already subset and QC to author exclude: Kim10x, Dowsett10x,Moffit10x -
#
# also exclude: Flynn10x
hypoMap_merged_raw@meta.data$Final_Exclude[hypoMap_merged_raw@meta.data$Dataset == "Flynn10x" & hypoMap_merged_raw@meta.data$Author_Exclude == "yes"] = "yes"

##########
### Subset
##########

hypoMap_merged_raw@meta.data$Final_Exclude = factor(hypoMap_merged_raw@meta.data$Final_Exclude,levels = c("yes","no"))
p1= DimPlot(hypoMap_merged_raw,group.by = "Final_Exclude",raster = F,shuffle = TRUE,label=F)#+NoLegend()
p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
p1r

hypoMap_merged_raw@meta.data$Final_Exclude = as.character(hypoMap_merged_raw@meta.data$Final_Exclude)

hypoMap_merged_raw_subset = subset(hypoMap_merged_raw,subset = Final_Exclude == "no")

##########
## save final object:
##########

# also save original with new columns:
saveRDS(hypoMap_merged_raw,file = paste0(global_parameters$data_path,"hypoMap_merged_raw.rds"))

# save filtered
saveRDS(hypoMap_merged_raw_subset,file = paste0(global_parameters$data_path,"hypoMap_merged_filtered.rds"))


##########
## old:
##########
#
# # 161 yes, 135yes, 111 no, 95 no, 98 no
# # 90 , 64 no
# clusters=c("49")
# for(c in clusters){
#   print(c)
#   temp_markers = FindMarkers(hypoMap_merged_raw,ident.1 = c,max.cells.per.ident = 3000,min.diff.pct = 0.1,logfc.threshold = 0.3,only.pos = TRUE)
#   temp_markers$gene = rownames(temp_markers)
#   assign(x = paste0("markers_",c) ,value = temp_markers)
# }
#
# Idents(hypoMap_merged_raw) = "Processing_clusters"
# temp_markers = FindMarkers(hypoMap_merged_raw,ident.1 = "27",max.cells.per.ident = 2000,min.diff.pct = 0.1,logfc.threshold = 0.3)
# temp_markers$gene = rownames(temp_markers)
#
# p1= FeaturePlot(hypoMap_merged_raw,features =  "Sst",raster = F,order = TRUE)
# p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
# p1r
#
# cellsh = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Processing_clusters=="90"]
# #cellsh = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Dataset=="Anderson10x"  & hypoMap_merged_raw@meta.data$Doublet =="Doublet"]
# p1= DimPlot(hypoMap_merged_raw,group.by = "Phase",raster = F,shuffle = TRUE,cells.highlight = cellsh,sizes.highlight = 0.1)
# p1r = scUtils::rasterize_ggplot(p1,pixel_raster = 1596)
# p1r
#
# table(hypoMap_merged_raw@meta.data$Dataset[hypoMap_merged_raw@meta.data$Processing_clusters=="49"])
