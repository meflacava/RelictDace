#!/bin/bash
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=admix
#SBATCH --time=7-00:00:00
#SBATCH --partition=high
#SBATCH --output=admix-%A-%a.out
#SBATCH --error=admix-%A-%a.err
#SBATCH --array=1-20

#Submit job with: sbatch NGSadmix.sh [path/to/beagle] [path/to/desired/output]

#Modify number of array jobs to test desired range of k values (e.g., --array=1-10 tests k=1-10)

#Example:
#sbatch /home/mlacava/filtering/rwpopgen/rwpopgen_PCA.beagle.gz /home/mlacava/filtering/rwpopgen/admix/ngsadmix_rwpopgen

#$1 #your path to beagle file
#$2 #your path to desired output file name (will add k tested to end)

module load angsd
module load ngsTools

NGSadmix -likes $1 -K ${SLURM_ARRAY_TASK_ID} -minMaf 0.05 -o $2_k${SLURM_ARRAY_TASK_ID} -P 10
