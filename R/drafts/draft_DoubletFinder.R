#https://github.com/chris-mcginnis-ucsf/DoubletFinder
seurat_object_doublet = seurat_object_batch_list$ChenDropseqbatch_1
pN_fixed = 0.25 # this can be the same globally
pK_max = 0.1 #  this can be the same globally
doublet_cluster_tresh = 0.5 # this can be the same globally
doublet_formation_rate = 0.05 # or 0.075  --->for nExp_poi, this one might differ between datasets
target_cluster_number = 50 # need to make a good estimate
min_avg_cells_cluster = 200

seurat_object_doublet = scUtils::seurat_recipe(seurat_object_doublet,nfeatures_vst = nfeatures_vst_prelim,
                                          clean_hvg = TRUE,
                                          normalize_data = TRUE,
                                          calcUMAP = TRUE,
                                          findClusters = TRUE,
                                          npcs_PCA = npcs_PCA,
                                          clusterRes = 0.5,
                                          seed = global_seed)
# determine cluster number:
# expected target_cluster_number should eitehr be constant or adjusted per dataset (manually based on known heterogeneity)
# additionally to avid too many small clusters after louvain clustering, the min_avg_cells_cluster sets a maximum based on the number of cells in the batch (e.g. to avoid separating a dataset with 3000 cells into 100 different clusters)
target_cluster_number_adjusted = floor(min(target_cluster_number, ncol(seurat_object_doublet) / min_avg_cells_cluster))

# execute: determine_cluster_resolution()

######## determine_cluster_resolution function for clusters:
target_cluster_number = 50 # need to make a good estimate
graph_name = "RNA_snn"
cluster_col_name = "seurat_clusters"
resolutions = c(0.5,0.75,1,1.5,2:10)
seed = 1234
return_seurat =TRUE)
# other args to Seurat::FindClusters
# run clustering (louvain)
seurat_object_doublet = Seurat::FindClusters(seurat_object_doublet,resolution = resolutions,graph.name = graph_name,random.seed =seed)
# identify res clostest to target
n_clusters_per_res = apply(seurat_object_doublet@meta.data[,paste0(graph_name,"_res.",resolutions)],2,function(x,min_cells){length(table(x)[table(x)>min_cells])},min_cells=5)
res_with_target_n = which(abs(n_clusters_per_res-target_cluster_number)== min(abs(n_clusters_per_res-target_cluster_number)))[1]
# return
if(return_seurat){return(seurat_object_doublet@meta.data[,cluster_col_name] = seurat_object_doublet@meta.data[,res_with_target_n])}else{return(res_with_target_n)}

#### doublet finder
# execute: apply_DoubletFinder()

library(DoubletFinder)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.seurat_object_doublet <- paramSweep_v3(seurat_object_doublet, PCs = 1:npcs_PCA, sct = FALSE,num.cores=n_cores)
sweep.stats_seurat_object_doublet <- summarizeSweep(sweep.res.seurat_object_doublet, GT = FALSE)
bcmvn_seurat_object_doublet <- find.pK(sweep.stats_seurat_object_doublet)
pK_at_BCmetrci_max = as.numeric(as.character(bcmvn_seurat_object_doublet$pK[bcmvn_seurat_object_doublet$BCmetric == max(bcmvn_seurat_object_doublet$BCmetric)]))
if(pK_at_BCmetrci_max > pK_max){pK_at_BCmetrci_max = pK_max}

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
# gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
# sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
# bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seurat_object_doublet@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(seurat_object_doublet@meta.data$Author_Class)
nExp_poi <- round(doublet_formation_rate*nrow(seurat_object_doublet@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
# nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seurat_object_doublet <- DoubletFinder::doubletFinder_v3(seurat_object_doublet, PCs = 1:npcs_PCA, pN = pN_fixed, pK = pK_at_BCmetrci_max, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
reuse.pANN = paste0("pANN_",pN_fixed,"_",pK_at_BCmetrci_max,"_",nExp_poi)
#seurat_object_doublet <- DoubletFinder::doubletFinder_v3(seurat_object_doublet, PCs = 1:npcs_PCA, pN = pN_fixed, pK = pK_at_BCmetrci_max, nExp = nExp_poi.adj, reuse.pANN = reuse.pANN, sct = FALSE)

DimPlot(seurat_object_doublet,group.by = paste0("DF.classifications_",gsub("pANN_","",reuse.pANN)))#"DF.classifications_0.25_0.11_736")#paste0("DF.classifications_",reuse.pANN))
FeaturePlot(seurat_object_doublet,features = reuse.pANN)
DimPlot(seurat_object_doublet,group.by = "seurat_clusters",label = TRUE,label.size = 2.5)+NoLegend()
DimPlot(seurat_object_doublet,group.by = "Author_Class",label = TRUE,label.size = 2.5)+NoLegend()
sort(tapply(seurat_object_doublet@meta.data$pANN_0.25_0.01_736,seurat_object_doublet@meta.data$seurat_clusters,FUN = "mean"),decreasing = TRUE)

# calculate pct of doublets per cluster
doublet_stats_per_cluster = seurat_object_doublet@meta.data %>% dplyr::select(cluster = !!rlang::sym("seurat_clusters"),doublet_finder_classification = !!rlang::sym(paste0("DF.classifications_",gsub("pANN_","",reuse.pANN)))) %>%
  dplyr::group_by(cluster) %>% dplyr::add_count(name="cells_per_cluster") %>% dplyr::filter(doublet_finder_classification=="Doublet") %>% dplyr::add_count(name="doublet_per_cluster") %>%
  dplyr::distinct(cluster,cells_per_cluster,doublet_per_cluster) %>% dplyr::mutate(doublet_pct = doublet_per_cluster / cells_per_cluster)
# also filter out full cluster with certain pct!
cells_in_doublet_clusters = seurat_object_doublet@meta.data$Cell_ID[ seurat_object_doublet@meta.data[,"seurat_clusters"] %in% doublet_stats_per_cluster$cluster[doublet_stats_per_cluster$doublet_pct >= doublet_cluster_tresh]]
DimPlot(seurat_object_doublet,group.by = paste0("DF.classifications_",gsub("pANN_","",reuse.pANN)),cells.highlight = cells_in_doublet_clusters,sizes.highlight = 0.2)#"DF.classifications_0.25_0.11_736")#paste0("DF.classifications_",reuse.pANN))

# final Doublet column:
seurat_object_doublet@meta.data$Doublet = seurat_object_doublet@meta.data[,paste0("DF.classifications_",gsub("pANN_","",reuse.pANN))]
seurat_object_doublet@meta.data$Doublet[seurat_object_doublet@meta.data$Cell_ID %in% cells_in_doublet_clusters] = "Doublet"



#########testing
DimPlot(seurat_object_doublet,group.by ="Doublet")

FeaturePlot(seurat_object_doublet,features = "Tgfb3",order = TRUE)

tmp_markers = FindMarkers(seurat_object_doublet,ident.1 = "35",min.diff.pct = 0.1,max.cells.per.ident = 4000)
tmp_markers$gene = rownames(tmp_markers)


tmp_markers_btw = FindMarkers(seurat_object_doublet,ident.1 = "35",ident.2 = "19",min.diff.pct = 0.1,max.cells.per.ident = 4000)
tmp_markers_btw$gene = rownames(tmp_markers_btw)
