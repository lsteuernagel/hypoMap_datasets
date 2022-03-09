seurat_object.list <- Seurat::SplitObject(seurat_raw, split.by = sample_column)
nfeatures_comparison = 3000
verbose =FALSE
assay="RNA"
sort_var_col = "vst.variance.standardized"
if(is.null(sort_var_col)){sort_var_col=4}
message("Finding shared variable features: ",nfeatures_comparison)
hvg_list = list()
for(i in 1:length(seurat_object.list)){
  message(names(seurat_object.list)[i])
  seurat_object.list[names(seurat_object.list)[i]] <- FindVariableFeatures(object = seurat_object.list[[i]],assay = assay,selection.method = "vst", nfeatures = nfeatures_comparison, verbose = verbose)
  #hvg_list[[names(seurat_object.list)[i]]] = seurat_object.list[[i]]@assays[[assay]]@var.features
  hvg_list[[names(seurat_object.list)[i]]]= as.data.frame(t(seurat_object.list[[i]]@assays[[assay]]@meta.features[,sort_var_col,drop=FALSE]))
}
# get matrix
hvg_matrix = as.data.frame(t(data.table::rbindlist(hvg_list,use.names = TRUE,fill=TRUE)))
# make ranks to only workj with the top x1000
hvg_matrix_rank = apply(hvg_matrix,MARGIN = 2,FUN = rank)
hvg_matrix_rank_median = apply(hvg_matrix_rank,1,median)
top_n_features_per_rank = names(sort(hvg_matrix_rank_median,decreasing = TRUE)[1:nfeatures_comparison])
hvg_matrix_comparison = hvg_matrix[top_n_features_per_rank,]

# Use correlation distances
distmat = as.dist(correlation_based_distance(hvg_matrix_comparison,method = "pearson"))
hc <- hclust(distmat, method = "ward.D")
plotClustering =TRUE
if(plotClustering){
  print("Plotting clustering as basis for agglomerative merging.")
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  heatmap.2( as.matrix(distmat), Rowv=as.dendrogram(hc),
             symm=TRUE, trace="none", col=colors,
             margins=c(2,10), labCol=FALSE )
}
hc$height

# hvg_list_ranked = list()
# for(i in 1:length(hvg_list)){
#   temp = data.frame(rank = 1:length(hvg_list[[names(hvg_list)[i]]]))
#   rownames(temp) = hvg_list[[names(hvg_list)[i]]]
#   temp = as.data.frame(t(temp))
#   hvg_list_ranked[[names(hvg_list)[i]]] = temp
# }





hvg_ranks = do.call(dplyr::bind_cols,hvg_list_ranked)
hvg_ranks = as.data.frame(t(data.table::rbindlist(hvg_list_ranked,use.names = TRUE,fill=TRUE)))
colnames(hvg_ranks) = names(hvg_list_ranked)

index=10
rownames(hvg_list_ranked$SRR4340023)[index]
rownames(hvg_list_ranked$SRR4340041)[index]

a1=(simplify2array(hvg_list_ranked))
dim(a1)
