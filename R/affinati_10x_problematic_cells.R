
##########
### Check on Affinati10x_seurat_processed
##########

Affinati10x_seurat_processed = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Affinati10x/Affinati10x_seurat_processed.rds")

###
## there is this large set of neurons in the middle that are weird and also there are many problematic cells marked as mixed
DimPlot(Affinati10x_seurat_processed,group.by = "Author_CellType",label=TRUE,label.size = 2.5)+NoLegend()
# the processing identified some clusters as doublets, most make sense
DimPlot(Affinati10x_seurat_processed,group.by = "Doublet",label=F,label.size = 2)
# the processing identified some clusters as doublets
DimPlot(Affinati10x_seurat_processed,group.by = "preliminary_clusters",label=TRUE,label.size = 2)+NoLegend()

# the cluster 35
cellsh = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("35") & hypoMap_merged_raw@meta.data$Dataset=="Affinati10x"]
DimPlot(Affinati10x_seurat_processed,group.by = "Author_Class",label = TRUE,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()
specifc_markers35 = FoldChange_pct(object = Affinati10x_seurat_processed@assays$RNA,cells.1 = cellsh,cells.2 = Affinati10x_seurat_processed@meta.data$Cell_ID[!Affinati10x_seurat_processed@meta.data$Cell_ID %in% cellsh])
specifc_markers35$gene = rownames(specifc_markers35)
FeaturePlot(Affinati10x_seurat_processed,"Avp")

# also: 67
cellsh = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("67") & hypoMap_merged_raw@meta.data$Dataset=="Affinati10x"]
DimPlot(Affinati10x_seurat_processed,group.by = "Author_Class",label = TRUE,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()
specifc_markers67 = FoldChange_pct(object = Affinati10x_seurat_processed@assays$RNA,cells.1 = cellsh,cells.2 = Affinati10x_seurat_processed@meta.data$Cell_ID[!Affinati10x_seurat_processed@meta.data$Cell_ID %in% cellsh])
specifc_markers67$gene = rownames(specifc_markers67)
FeaturePlot(Affinati10x_seurat_processed,"Avp")
#,"45","96" are tany and ependymo, 69 astrcytes

# the cluster 31
cellsh = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("31")]
DimPlot(Affinati10x_seurat_processed,group.by = "Author_Class",label = TRUE,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()
# check agrp
FeaturePlot(Affinati10x_seurat_processed,"Agrp")
# get agrp clusters to compare: 68 vs 44
Idents(Affinati10x_seurat_processed) = "preliminary_clusters"
within_agrp_markers = FindMarkers(Affinati10x_seurat_processed,ident.1 = "44",ident.2 = c("68"),max.cells.per.ident = 2000,min.diff.pct = 0.1,logfc.threshold = 0.3)
within_agrp_markers$gene = rownames(within_agrp_markers)

# different genes:
FeaturePlot(Affinati10x_seurat_processed,"Xist")

cellsi = Affinati10x_seurat_processed@meta.data$Cell_ID[Affinati10x_seurat_processed@meta.data$preliminary_clusters=="10"]
a1=hypoMap_merged_raw@meta.data[hypoMap_merged_raw@meta.data$Cell_ID %in% cellsi, ]
# the cluster 27
cellsh = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("27")]
DimPlot(Affinati10x_seurat_processed,group.by = "Author_Class",label = TRUE,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()

# the cluster 69
cellsh = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("69")]
DimPlot(Affinati10x_seurat_processed,group.by = "Author_Class",label = TRUE,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()

# the cluster 69
cellsh = hypoMap_merged_raw@meta.data$Cell_ID[hypoMap_merged_raw@meta.data$Processing_clusters %in% c("96")]
DimPlot(Affinati10x_seurat_processed,group.by = "Author_Class",label = TRUE,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()


s##########
### Check on Anderson10x_seurat_processed
##########

Anderson10x_seurat_processed = readRDS("/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/Anderson10x/Anderson10x_seurat_processed.rds")
DimPlot(Anderson10x_seurat_processed,group.by = "Doublet")
DimPlot(Anderson10x_seurat_processed,group.by = "preliminary_clusters",label = TRUE)+NoLegend()
DimPlot(Anderson10x_seurat_processed,group.by = "Author_Class",label = TRUE)+NoLegend()


Idents(Anderson10x_seurat_processed) = "preliminary_clusters"
temp_markers = FindMarkers(Anderson10x_seurat_processed,ident.1 = "44",max.cells.per.ident = 2000,min.diff.pct = 0.1,logfc.threshold = 0.3)
temp_markers$gene = rownames(temp_markers)

FeaturePlot(Anderson10x_seurat_processed,"Plp1")
# find cluster 60 cells

temp_markers2 = FindMarkers(Anderson10x_seurat_processed,ident.1 = "44",ident.2 = c("34","6","25","24","29"),max.cells.per.ident = 2000,min.diff.pct = 0.1,logfc.threshold = 0.3)
temp_markers2$gene = rownames(temp_markers2)

cellsh = Anderson10x_seurat_processed@meta.data$Cell_ID[Anderson10x_seurat_processed@meta.data$preliminary_clusters=="44"]

DimPlot(Anderson10x_seurat_processed,group.by = "Author_Class",label = TRUE,cells.highlight = cellsh,sizes.highlight = 0.1)+NoLegend()


target_cluster = "35"
cells1=Anderson10x_seurat_processed@meta.data$Cell_ID[Anderson10x_seurat_processed@meta.data$preliminary_clusters==target_cluster & Anderson10x_seurat_processed$Doublet=="Doublet"]
cells2=Anderson10x_seurat_processed@meta.data$Cell_ID[Anderson10x_seurat_processed@meta.data$preliminary_clusters==target_cluster & Anderson10x_seurat_processed$Doublet=="Singlet"]

Idents(Anderson10x_seurat_processed) = "Cell_ID"
specifc_markers = FoldChange_pct(object = Anderson10x_seurat_processed@assays$RNA,cells.1 = cells1,cells.2 = cells2)
specifc_markers$gene = rownames(specifc_markers)



