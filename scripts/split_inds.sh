#!/bin/bash
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=split
#SBATCH --time=24:00:00
#SBATCH --partition=med
#SBATCH --array=1-15

file=$(ls /home/mlacava/scripts/split_tasks | sed -n ${SLURM_ARRAY_TASK_ID}p)
echo array task $SLURM_ARRAY_TASK_ID
echo $file
bash /home/mlacava/scripts/split_tasks/$file
