##########
### [0] Load
##########

# direct output and logs to some files on the local filesystem:
# where to store temporary json's with params for jobs:
param_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_params/"
# where to save log files --> use this path in the slurm.sh files!
log_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/hypoMap_v2_slurmlogs/"

# load overview table with per dataset information
dataset_table = data.table::fread("data/dataset_overview.tsv",data.table = FALSE)

# load json file with all other information
params_pre_processing = jsonlite::read_json("data/parameters_pre_processing_v2_1.json")

### try to creat dir if necessary:
system(paste0("mkdir -p ",paste0(param_path)))
system(paste0("mkdir -p ",paste0(log_path)))

##########
### [1] Run pre-processing and batch id job for each dataset
##########

slurm_id_1_per_dataset = vector()

for(i in 1:nrow(dataset_table)){
  # set parameters for current dataset
  param_set = params_pre_processing
  param_set$dataset_name = dataset_table$Dataset[i]
  param_set$raw_file = dataset_table$seurat_file[i]
  param_set$target_cluster_number = dataset_table$estimated_total_clusters[i] # need to make a good estimate
  param_set$exclude_author = dataset_table$exclude_author[i]
  # make unique id:
  job_id=digest::digest(param_set)
  # write to JSON as transfer file
  param_file = paste0(param_path,"preprocessing_params_",job_id,".json")
  scUtils::writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "run_scripts/preprocessing.R"
  # set sbatch params:
  jobname = paste0("preprocessing_",dataset_name,"_",job_id)
  outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
  errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
  #run job
  slurm_id_1_per_dataset[dataset_name] = system(paste0("sbatch -j ",jobid," -o ",outputfile," -e ",errorfile," bash/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_path),intern = TRUE)
}


##########
### [2] Run doublet detection job for each dataset ---> depends on [1] for each dataset
##########

slurm_id_2_per_dataset = vector()

for(i in 1:nrow(dataset_table)){
  # set parameters for current dataset
  param_set = params_pre_processing
  param_set$dataset_name = dataset_table$Dataset[i]
  param_set$raw_file = dataset_table$seurat_file[i]
  param_set$doublet_formation_rate = dataset_table$doublet_formation_rate[i] # or 0.075  --->for nExp_poi, this one might differ between datasets
  # make unique id:
  job_id=digest(param_set)
  # write to JSON as transfer file
  param_file = paste0(param_path,"doublet_detection_params_",job_id,".json")
  scUtils::writeList_to_JSON(list_with_rows = param_set,filename = param_file)
  # execute job
  script_path = "run_scripts/doublet_detection.R"
  # set sbatch params:
  jobname = paste0("doublet_detection_",dataset_name,"_",job_id)
  outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
  errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
  dependency_ids = slurm_id_1_per_dataset[dataset_name]
  slurm_id_2_per_dataset[dataset_name] = system(paste0("sbatch -j ",jobid," -o ",outputfile," -e ",errorfile," --dependency=afterok:",dependency_ids," bash/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_path),intern = TRUE)
}


##########
### [3] Merge final result and run pre-processing+HVG detection, save as merged dataset ---> depends on all [2]
##########

# TODO
# .... param_set$
# need: which feature sizes to add
#   exclude_author = dataset_table$Dataset[i]

# write to JSON as transfer file
param_file = paste0(param_path,"merge_datasets_params_",job_id,".json")
scUtils::writeList_to_JSON(list_with_rows = param_set,filename = param_file)
# not a loop
script_path = "run_scripts/merge_datasets.R"
# set sbatch params:
jobname = paste0("merge_datasets_",dataset_name,"_",job_id)
outputfile = paste0(log_path,jobname,"_","slurm-%j.out")
errorfile = paste0(log_path,jobname,"_","slurm-%j.err")
dependency_ids = paste0(slurm_id_2_per_dataset,collapse = ":")
slurm_id_3 = system(paste0("sbatch -j ",jobid," -o ",outputfile," -e ",errorfile," --dependency=afterok:",dependency_ids," bash/run_Rscript_slurm.sh ",singularity_path," ",script_path," ",param_path),intern = TRUE)


