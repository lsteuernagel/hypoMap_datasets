###

# make a slum executable script using functions from package scUtils above that for loops over all seurats:
# Run basic pre-processing on all samples - merged
# Run batch-detection
# For each batch:
# re-run pre-processing
# run standard clustering
# Detect Doublets:
# (https://github.com/chris-mcginnis-ucsf/DoubletFinder)
# select cluster with high doublet scores and remove
# Export final seurat objects

# function to match samples:

match_sample_names = function(table_A,table_B,sample_col_A,sample_col_B,barcode_col = "Barcode",min_pct=0.8){

  A_barcode_list = split(table_A[,barcode_col],table_A[,sample_col_A])
  B_barcode_list = split(table_B[,barcode_col],table_B[,sample_col_B])

  resvec=vector()
  for(i in 1:length(A_barcode_list)){
    current_bcs = A_barcode_list[[i]]
    overlap_bcs = sapply(B_barcode_list,function(x,el){length(intersect(x,el))},el=current_bcs)
    likely_sample=names(which(overlap_bcs==max(overlap_bcs) & max(overlap_bcs) > min_pct*length(current_bcs)))
    if(length(likely_sample)>0){
      message("Associating ",names(A_barcode_list)[i]," with ",likely_sample, " based on an overlap of ",max(overlap_bcs)," out of ",length(current_bcs)," barcodes.")
      resvec[names(A_barcode_list)[i]] =  likely_sample
    }else{
      resvec[names(A_barcode_list)[i]] =  NA
    }
  }
  resdf = data.frame(Samples_A = names(resvec),Samples_B = resvec)
  return(resdf)
}


## function to infer sex per sample:
infer_sex = function(seurat_object,sample_column,id_column,min_xist_female = 0.7,max_xist_male = 0.1){
  xist_expr = Seurat::FetchData(seurat_object,"Xist")[,1]
  xist_expr[xist_expr>0]=1
  tmp_df = data.frame(Cell_ID = seurat_object@meta.data[,id_column], Sample_ID =  seurat_object@meta.data[,sample_column],xist_expr=xist_expr) %>% group_by(Sample_ID) %>%
    dplyr::add_count(name="n_cells")  %>% dplyr::mutate(n_xist = sum(xist_expr)) %>% dplyr::mutate(xist_pct = n_xist / n_cells) %>%
    dplyr::mutate(inferred_sex = case_when(xist_pct > min_xist_female ~ "F", xist_pct < max_xist_male ~ "M", TRUE ~ "U"))
  return(tmp_df$inferred_sex)
}





