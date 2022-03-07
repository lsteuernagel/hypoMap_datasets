##########
### start
##########

library(dplyr)
library(Seurat)
library(scUtils)

## read raw file and extract metadata
Mousebrainorg10x_raw_data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Mousebrainorg10x/"
Mousebrainorg10x_seurat = readRDS(paste0(Mousebrainorg10x_raw_data_path,"Mousebrainorg10x_seurat_raw.rds"))

Mousebrainorg10x_seurat_PROCESSED = subset(Mousebrainorg10x_seurat,subset = nCount_RNA > 1000 & percent_mt < 10)
Mousebrainorg10x_seurat_PROCESSED = seurat_recipe(Mousebrainorg10x_seurat_PROCESSED,nfeatures_vst = 1200,calcUMAP = TRUE,findClusters = TRUE,clusterRes = 3)
DimPlot(Mousebrainorg10x_seurat_PROCESSED,group.by = "seurat_clusters")
DimPlot(Mousebrainorg10x_seurat_PROCESSED,group.by = "Author_CellType",label = TRUE,label.size = 2.4,order = TRUE)+NoLegend()
FeaturePlot(Mousebrainorg10x_seurat_PROCESSED,"Meis2")

##########
### Map onto full hypoMap
##########

### map onto full map
library(mapscvi)
Mousebrainorg10x_seurat_PROCESSED_mapped_full = mapscvi::map_new_seurat_hypoMap(Mousebrainorg10x_seurat_PROCESSED,reference_mode = "hypoMap_full", suffix="query_morris_mmSCN11_full",max_epochs=20)
Mousebrainorg10x_seurat_PROCESSED_mapped_full@meta.data$predicted_Curated_Class = Mousebrainorg10x_seurat_PROCESSED_mapped_full@meta.data$predicted
Mousebrainorg10x_seurat_PROCESSED_mapped_full@meta.data$prediction_probability_Curated_Class = Mousebrainorg10x_seurat_PROCESSED_mapped_full@meta.data$prediction_probability
#plot
plot_query_labels(query_seura_object=Mousebrainorg10x_seurat_PROCESSED_mapped_full,reference_seurat=mapscvi::reference_hypoMap_full,label_col="Curated_Class",
                  label_col_query = "predicted_Curated_Class",overlay = TRUE,query_pt_size = 0.4,labelonplot = TRUE,label.size=2)
FeaturePlot(Mousebrainorg10x_seurat_PROCESSED_mapped_full,features = "prediction_probability",reduction = "umap_scvi")

##########
### Map onto neuron hypoMap
##########

neurons_with_high_prob = Mousebrainorg10x_seurat_PROCESSED_mapped_full@meta.data$Cell_ID[Mousebrainorg10x_seurat_PROCESSED_mapped_full$predicted_Curated_Class == "Neurons" & Mousebrainorg10x_seurat_PROCESSED_mapped_full$prediction_probability_Curated_Class > 0.75]

library(mapscvi)
Mousebrainorg10x_seurat_PROCESSED_mapped_neurons = mapscvi::map_new_seurat_hypoMap(Mousebrainorg10x_seurat_PROCESSED_mapped_full,reference_mode = "hypoMap_neurons", suffix="query_morris_mmSCN11_neurons",subset_col = "Cell_ID",
                                                                       subset_values = neurons_with_high_prob,label_col = "K98_pruned",max_epochs=20)
Mousebrainorg10x_seurat_PROCESSED_mapped_neurons@meta.data$predicted_K98_pruned = Mousebrainorg10x_seurat_PROCESSED_mapped_neurons@meta.data$predicted
Mousebrainorg10x_seurat_PROCESSED_mapped_neurons@meta.data$predicted_K98_named = mapscvi::add_paired_annotation(input_annotation = Mousebrainorg10x_seurat_PROCESSED_mapped_neurons@meta.data$predicted_K98_pruned,
                                                                                                    reference_annotations = mapscvi::reference_hypoMap_neurons@meta.data[,c("K98_pruned","K98_named")])

#plot
plot_query_labels(query_seura_object=Mousebrainorg10x_seurat_PROCESSED_mapped_neurons,reference_seurat=mapscvi::reference_hypoMap_neurons,label_col="K98_named",
                  label_col_query = "predicted_K98_named",overlay = TRUE,query_pt_size = 0.4,labelonplot = TRUE,label.size=2)
FeaturePlot(Mousebrainorg10x_seurat_PROCESSED_mapped_neurons,features = "prediction_probability",reduction = "umap_scvi")

FeaturePlot(Mousebrainorg10x_seurat_PROCESSED_mapped_neurons,features = "Ghrh",reduction = "umap_scvi")

DimPlot(Mousebrainorg10x_seurat_PROCESSED_mapped_neurons,group.by = "Author_CellType",label = TRUE,label.size = 2.4,order = TRUE,reduction = "umap_scvi")+NoLegend()

