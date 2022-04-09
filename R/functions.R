
writeList_to_JSON = function (list_with_rows, filename)
{
  jsonfile = jsonlite::toJSON(list_with_rows, pretty = TRUE,
                              auto_unbox = TRUE, digits = NA)
  writeLines(jsonfile, con = paste0(filename))
}


##########
### function: identifyBatches
##########

#' Use random forests to identify batches within datasets by merging similar samples!
#' @param seurat_object
#' @param sample_variable
#' @param cell_id
#' @param max_entropy
#' @param trees
#' @param sampsize_pct
#' @param n_cores
#' @param n_dim
#' @param embedding_name
#' @param plotClustering
#' @param plotEntropy
#' @param seed
#' @return list with new_batch_labels, clustering and entrop results

## be sure that Cell_ID and assay rownames are identical !!!!
identifyBatches = function(seurat_object ,  sample_variable = "Sample_ID", cell_id = "Cell_ID", max_entropy = 0.9 ,trees= 1000,sampsize_pct=0.632, n_cores = 10,n_dim = 50,embedding_name = "pca",plotClustering =FALSE,plotEntropy=TRUE,seed=1234){

  require(tidyverse)
  require(stringr)
  require(data.table)
  require(ggplot2)
  require(gplots)
  require(RColorBrewer)
  require(Seurat)

  set.seed(seed)
  ### Use a random forest to find similar batch labels
  print("Running initial random forest to use class probabilities as input to correlation-based clustering")
  # based on uncorrected pca
  if(ncol(seurat_object@reductions[[embedding_name]]@cell.embeddings)>n_dim){n_dim = ncol(seurat_object@reductions[[embedding_name]]@cell.embeddings)}
  train_predictors_all = seurat_object@reductions[[embedding_name]]@cell.embeddings[,1:n_dim]
  train_response_all = as.factor(seurat_object@meta.data[,sample_variable])
  # length of current_labels
  n_batches = length(unique(train_response_all))
  # run random forest
  votes = rf_classProbabilities(train_predictors=train_predictors_all,train_response=train_response_all,trees=trees,n_cores=n_cores,seed1=seed)
  # Use correlation distances
  distmat = as.dist(correlation_based_distance(votes))
  # Cluster and merge samples (iteratively)
  hc <- hclust(distmat, method = "ward.D")
  if(plotClustering){
    print("Plotting clustering as basis for agglomerative merging.")
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    heatmap.2( as.matrix(distmat), Rowv=as.dendrogram(hc),
               symm=TRUE, trace="none", col=colors,
               margins=c(2,10), labCol=FALSE )
  }

  ### Define new_batch based on RF

  # - Start with all samples as their own cluster

  # get embedding for rf
  all_features = seurat_object@reductions[[embedding_name]]@cell.embeddings[,1:n_dim]
  # result list
  all_entropies = list()
  # artificially inserted column for no merging
  all_entropies[[1]] = rep(1,nrow(seurat_object@meta.data))
  # temorary rf results
  stored_cluster_votes =list()
  # - go through cutoffs to merge
  print("Iterate over potential clusters: ")
  for(j in 1:length(hc$height)){
    print(paste0(">>>>>j: ",j))
    cluster_labels = cutree(hc,h = round(hc$height[j],5)+0.0001)
    # go through all clusters
    unique_cluster_labels = unique(cluster_labels)
    all_labels_entropy = list()
    for(k in (unique_cluster_labels)){
      # - subset to cells of the given cluster
      labels = names(cluster_labels)[cluster_labels==k]
      all_current_cells = as.character(seurat_object@meta.data[,cell_id] [seurat_object@meta.data[,sample_variable] %in% labels])
      all_current_labels = seurat_object@meta.data[,sample_variable][seurat_object@meta.data[,sample_variable] %in% labels]
      name = paste0(labels,collapse = "__")
      if(length(labels)==1){
        #(if only 1 sample per cluster: set entropy to 1
        cluster_entropy = rep(1,length(all_current_cells))
        names(cluster_entropy)=all_current_cells
      }else{
        # if this combination was already run in rf used saved result:
        if(name %in% names(stored_cluster_votes)){
          print(paste0("Using previously calculated probabilities for: ",name))
          cluster_votes = stored_cluster_votes[[name]]
          cluster_entropy = apply(cluster_votes,1,entropy,logfun="log2") / log2(length(labels))
        }else{
          print(paste0("Running random forest for: ",paste0(unique(all_current_labels),collapse = "__")))
          # - train random forest only on these cells and only with cluster label
          cluster_votes = rf_classProbabilities(train_predictors=all_features[all_current_cells,],train_response=as.factor(all_current_labels),trees=trees,n_cores=n_cores)
          # save
          stored_cluster_votes[[name]] = cluster_votes
          # - calculate entropy on class probabilities and normalize by total samples in cluster
          cluster_entropy = apply(cluster_votes,1,entropy,logfun="log2") / log2(length(labels))

        }
      }
      # save result
      all_labels_entropy[[k]] = cluster_entropy
    }
    # - fuse all cells per cutoff level
    all_labels_entropy_fused = unlist(all_labels_entropy)
    all_entropies[[j+1]] = all_labels_entropy_fused
  }
  # - compare between cutoffs
  all_entropies_unsorted = as.data.frame(do.call(cbind,all_entropies))
  median_entropy = apply(all_entropies_unsorted,2,median)

  if(plotEntropy){
    print("Plotting median entropy of different merges.")
    print(median_entropy)
    plot(length(median_entropy):1,median_entropy)
    abline(h = max_entropy)
    # plot(11:1,c(1,apply(all_entropies_unsorted,2,median)))
  }

  # >>> cut
  # get mean and find level where mean is > threshold
  best_merge_level = max(as.numeric(gsub("V","",names(median_entropy)[median_entropy >= max_entropy])))-1 # -1 because of the artificially inserted column for no merging
  print(paste0("Merging into ",length(all_entropies)+1-best_merge_level," batches"))
  # make vector with new labels
  new_label_vec=cutree(hc,h = hc$height[best_merge_level]+0.00001)

  ### assign new labels
  new_batch = rep(NA,nrow(seurat_object@meta.data))
  for(i in 1:length(unique(new_label_vec))){
    names_i = names(new_label_vec)[new_label_vec==i]
    new_batch[seurat_object@meta.data[,sample_variable] %in% names_i] = paste0("batch_",i)
  }
  print("Finalized batch merging. New labels:")
  print(table(new_batch))

  # return both the new batches and also all entropies in acse the threshold has to be changed.
  results = list(new_batch = new_batch,all_entropies_unsorted=all_entropies_unsorted, clustering = hc)
  return(results)

}

### Function to plot and cut entropy

cut_batches = function(all_entropies_unsorted,max_entropy=0.9,hc,seurat_object,sample_variable){

  print("Plotting median entropy of different merges.")
  median_entropy = apply(all_entropies_unsorted,2,median)
  print(median_entropy)
  plot(length(median_entropy):1,median_entropy)
  abline(h = max_entropy)

  # >>> cut
  # get mean and find level where mean is > threshold
  median_entropy = apply(all_entropies_unsorted,2,median)
  best_merge_level = max(as.numeric(gsub("V","",names(median_entropy)[median_entropy >= max_entropy])))-1 # -1 because of the artificially inserted column for no merging
  # make vector with new labels
  new_label_vec=cutree(hc,h = hc$height[best_merge_level]+0.00001)

  ### assign new labels
  new_batch = rep(NA,nrow(seurat_object@meta.data))
  for(i in 1:length(unique(new_label_vec))){
    names_i = names(new_label_vec)[new_label_vec==i]
    new_batch[seurat_object@meta.data[,sample_variable] %in% names_i] = paste0("batch_",i)
  }
  return(new_batch)

}


### one RF round function
# input:
# PCA emebedding
# current labels
# trees
# ncores
rf_classProbabilities = function(train_predictors,train_response,trees=500,sampsize_pct=0.632,n_cores=4,seed1=123){

  require(tidyverse)
  require(randomForest)
  require(foreach)
  require(doParallel)
  # length of current_labels
  n_batches = length(unique(train_response))
  # run random forest
  registerDoParallel(cores=n_cores)
  rf_res <- foreach(ntree=rep(trees/n_cores, n_cores), .combine=randomForest::combine,.multicombine=TRUE, .packages='randomForest') %dopar% {
    set.seed(seed1)
    randomForest(x=train_predictors, y=as.factor(train_response), ntree=ntree,sampsize=ceiling(sampsize_pct*nrow(train_predictors)))
  }
  # obtain class probabilities
  votes = rf_res$votes / n_cores
  return(votes)
}
# get correlation based distances
correlation_based_distance = function(x,method="spearman"){
  result_mat = matrix(nrow = ncol(x),ncol = ncol(x))
  for(i in 1:(ncol(x)-1)){
    for(j in (i):ncol(x)){
      # see http://www.econ.upf.edu/~michael/stanford/maeb6.pdf ,
      # also: https://stats.stackexchange.com/questions/232758/distance-and-correlations  & https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering
      # this is the same:  mycor = function(x1,x2){sqrt(2*(1-cor(x1,x2,method = "spearman")))}
      mycor = function(x1,x2){sqrt(2 - (2*cor(x1,x2,method = "spearman")))}
      dists = mycor(x[,i],x[,j])
      result_mat[i,j] = dists
      result_mat[j,i] = dists
    }
  }
  result_mat[ncol(x),ncol(x)]=0
  colnames(result_mat) = colnames(x)
  rownames(result_mat) = colnames(x)
  return(result_mat)
}

### helper functions
entropy= function(x,logfun ="log2"){
  log_vec = do.call(logfun,list(x))
  log_vec[is.infinite(log_vec)] = 0
  return(-sum(x * log_vec))
}
