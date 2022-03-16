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

This directory contains subdirectories for each dataset. 

The script "R/raw_hypoMap_datasets.R" loads all data and creates a "_seurat_raw.rds" file in each subdirectory containing all cells and metadata from SRA.
Our own nucseq data is processed in an additional script in a similar way.

## Curation of each dataset with Author metadata and information from publications where possible

After creation of the seurat objects further metadata curation is required. For each dataset we add (if available) author annotations regarding both the samples as well as cells (cluster, passed-qc etc.). We curate on high lvele annotation: "Author_Class" using the same names across all datasets.
After curation the "_seurat_raw.rds" file is overwritten with the updated data.

See "R/curate_metadata/" for a script per dataset (requires adjustement of hardcoded filepaths!)

TODO: Add more details on added columns like AUthor_Class

## Execute processing pipeline

The core script to execute the slurm jobs is "R/execute_hypoMap_datasets.R" which consists of 4 steps:

- Load input data overview and parameters
- Preprocessing
- Doublet-detection
- Merging

Step 2-4 are called as slurm jobs that are dependent on the previous step. For each step the main script exports required parameters as a json file that will then be loaded by an R script which is called in an bash script ("R/run_scripts/run_Rscript_slurm.sh") that is executed by sbatch.

### Input datasets

Based on the table stored in [data/dataset_overview.tsv] datasets are loaded and processed. .

It looks like this:

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

The (not-reproducible) script: [R/make_dataset_table.R] was used to build this table. We manually added estimations of expected clusters for the processing based on prior work.

### Parameters

The parameters are loaded from a json files that can be adjusted (or copied and edited). [R/make_parameter_json.R] can be used to create the json.

For the hypoMap pipeline we used this parameter file: [data/parameters_pre_processing_v2_1.json] 

It looks like this:

```yaml
{
  "data_path": "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_rawdata/",
  "processed_suffix": "_seurat_processed",
  "merged_name": "hypoMap_merged_raw",
  "feature_set_sizes": [750, 1000, 1250, 1500, 2000, 2500, 3000, 4000],
  "n_cores": 50,
  "id_column": "Cell_ID",
  "global_seed": 123456,
  "minUMI": 1000,
  "minFeatures": 500,
  "maxUMI": 1000000,
  "maxUMI_dynamic": 20,
  "max_mt": 10,
  "genes_to_exclude_file": "data/features_exclude_list.json",
  "max_pctExclude": 40,
  "min_cells_sample": 100,
  "sample_column": "Sample_ID",
  "nfeatures_vst_prelim": 1000,
  "npcs_PCA": 70,
  "k_param": 30,
  "resolutions_to_try": [0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10],
  "downsampling_max": 20000,
  "max_entropy_batch_detection": 0.8,
  "trees_rf": 20000,
  "pN_fixed": 0.25,
  "pK_max": 0.1,
  "doublet_cluster_tresh": 0.7,
  "min_avg_cells_cluster": 200
}
```

TODO: Explain parameters.

### Preprocessing

lorem ipsum

