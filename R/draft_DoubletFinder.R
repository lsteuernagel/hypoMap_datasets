#https://github.com/chris-mcginnis-ucsf/DoubletFinder
seurat_object_doublet = seurat_object_batch_list$ChenDropseqbatch_1
pN_fixed = 0.25

seurat_object_doublet = scUtils::seurat_recipe(seurat_object_doublet,nfeatures_vst = nfeatures_vst_prelim,
                                          clean_hvg = TRUE,
                                          normalize_data = TRUE,
                                          calcUMAP = TRUE,
                                          findClusters = TRUE,
                                          npcs_PCA = npcs_PCA,
                                          clusterRes = 2.5,
                                          seed = global_seed)

library(DoubletFinder)

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.seurat_object_doublet <- paramSweep_v3(seurat_object_doublet, PCs = 1:npcs_PCA, sct = FALSE,num.cores=n_cores)
sweep.stats_seurat_object_doublet <- summarizeSweep(sweep.res.seurat_object_doublet, GT = FALSE)
bcmvn_seurat_object_doublet <- find.pK(sweep.stats_seurat_object_doublet)
pK_at_BCmetrci_max = as.numeric(as.character(bcmvn_seurat_object_doublet$pK[bcmvn_seurat_object_doublet$BCmetric == max(bcmvn_seurat_object_doublet$BCmetric)]))

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
# sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
# gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results
# sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
# bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(seurat_object_doublet@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(seurat_object_doublet@meta.data$Author_Class)
nExp_poi <- round(0.075*nrow(seurat_object_doublet@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seurat_object_doublet <- DoubletFinder::doubletFinder_v3(seurat_object_doublet, PCs = 1:npcs_PCA, pN = pN_fixed, pK = pK_at_BCmetrci_max, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
reuse.pANN = paste0("pANN_",pN_fixed,"_",pK_at_BCmetrci_max,"_",nExp_poi)
seurat_object_doublet <- DoubletFinder::doubletFinder_v3(seurat_object_doublet, PCs = 1:npcs_PCA, pN = pN_fixed, pK = pK_at_BCmetrci_max, nExp = nExp_poi.adj, reuse.pANN = reuse.pANN, sct = FALSE)

DimPlot(seurat_object_doublet,group.by = paste0("DF.classifications_",gsub("pANN_","",reuse.pANN)))#"DF.classifications_0.25_0.11_736")#paste0("DF.classifications_",reuse.pANN))
FeaturePlot(seurat_object_doublet,features = reuse.pANN)
DimPlot(seurat_object_doublet,group.by = "seurat_clusters",label = TRUE,label.size = 2.5)+NoLegend()
DimPlot(seurat_object_doublet,group.by = "Author_Class",label = TRUE,label.size = 2.5)+NoLegend()
sort(tapply(seurat_object_doublet@meta.data$pANN_0.25_0.01_736,seurat_object_doublet@meta.data$seurat_clusters,FUN = "mean"),decreasing = TRUE)

# todo calculate pct of doublets per cluster
# also filter out full cluster with certain pct!
