# hypoMap_datasets
Processing and curation of all datasets that are part of HypoMap v2

# Overview 

This repository contains all code to create the basic merged hypoMap object prior ton integration and harmonization/annotation.
Key steps:

- Loading and conversion of matrix files and SRA metdata to Seurat objects
- curation of each dataset with Author metadata and information from publications where possible
- basic QC and preprocessing
- additional batch detection within each dataset to group highly similar samples and avoid fragmentation into too many samples.
- Doublet detection
- Merging to final Seurat object and export to h5 files

It is not a fully automatized R package due to the manual work involved in curating the individual datasets. It relies on the following R packages:
[scUtils](https://github.sf.mpg.de/lsteuernagel/scUtils) and [Seurat](https://satijalab.org/seurat/) as well as some other minor dependencies that should be loaded via the other two.

Most code will be executed via slurm jobs using a singularity image with all required dependencies that could fro example be pulled from a suitable [Docker image](https://hub.docker.com/r/lsteuernagel/r_scvi/tags).

The execute and queue management scripts are writte in R as well and require some packages like digest (for hashing), magrittr and dplyr.

# Step-by step 

This section intends to explain the different scripts involved and should allow to reproduce the above steps.

## Loading and conversion of matrix files and SRA metdata to Seurat objects

The input are the raw cellranger and dropseq tools files as well as tables from the SRA that come out of our alignment pipeline and should be stored someweher on Disk. In our case:

```
/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/
```

This directory contains subdirectories for each dataset. Based on the table stored in

```
/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/
```

which looks like this:

| seurat_file                                      | Dataset          | estimated_total_clusters | exclude_author | doublet_formation_rate |
|--------------------------------------------------|------------------|--------------------------|----------------|------------------------|
| Affinati10x/Affinati10x_seurat_raw.rds           | Affinati10x      | 90                       | FALSE          | 0.05                   |
| Anderson10x/Anderson10x_seurat_raw.rds           | Anderson10x      | 70                       | FALSE          | 0.05                   |
| CampbellDropseq/CampbellDropseq_seurat_raw.rds   | CampbellDropseq  | 67                       | FALSE          | 0.05                   |
| ChenDropseq/ChenDropseq_seurat_raw.rds           | ChenDropseq      | 76                       | FALSE          | 0.05                   |
| Dowsett10xnuc/Dowsett10xnuc_seurat_raw.rds       | Dowsett10xnuc    | 90                       | TRUE           | 0.02                   |
| Flynn10x/Flynn10x_seurat_raw.rds                 | Flynn10x         | 41                       | FALSE          | 0.05                   |
| Kim10x/Kim10x_seurat_raw.rds                     | Kim10x           | 72                       | TRUE           | 0.02                   |
| kimDev10x/kimDev10x_seurat_raw.rds               | kimDev10x        | 19                       | FALSE          | 0.02                   |
| LeeDropseq/LeeDropseq_seurat_raw.rds             | LeeDropseq       | 48                       | FALSE          | 0.02                   |
| Mickelsen10x/Mickelsen10x_seurat_raw.rds         | Mickelsen10x     | 43                       | FALSE          | 0.05                   |
| Moffit10x/Moffit10x_seurat_raw.rds               | Moffit10x        | 64                       | TRUE           | 0.05                   |
| Morris10x/Morris10x_seurat_raw.rds               | Morris10x        | 50                       | FALSE          | 0.05                   |
| Mousebrainorg10x/Mousebrainorg10x_seurat_raw.rds | Mousebrainorg10x | 37                       | FALSE          | 0.05                   |
| RomanovDev10x/RomanovDev10x_seurat_raw.rds       | RomanovDev10x    | 18                       | FALSE          | 0.05                   |
| RossiDropseq/RossiDropseq_seurat_raw.rds         | RossiDropseq     | 51                       | FALSE          | 0.05                   |
| Rupp10x/Rupp10x_seurat_raw.rds                   | Rupp10x          | 90                       | FALSE          | 0.05                   |
| wen10x/wen10x_seurat_raw.rds                     | wen10x           | 21                       | FALSE          | 0.05                   |
| wenDropseq/wenDropseq_seurat_raw.rds             | wenDropseq       | 61                       | FALSE          | 0.05                   |
