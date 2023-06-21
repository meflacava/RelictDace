#!/bin/bash
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=admix
#SBATCH --time=7-00:00:00
#SBATCH --partition=med
#SBATCH --output=/home/mlacava/slurm_logs/admix-%A-%a.out
#SBATCH --error=/home/mlacava/slurm_logs/admix-%A-%a.err
#SBATCH --array=1-20

#Submit job with: sbatch NGSadmix.sh [path/to/input/beagle] [path/to/desired/output]
#EDIT: number of K values you want to test above in --array=1-n (e.g., 1-10)

#Example:
#sbatch NGSadmix.sh /home/mlacava/filtering/rwpopgen/rwpopgen_PCA.beagle.gz /home/mlacava/filtering/rwpopgen/admix/ngsadmix_rwpopgen

#$1 #your path to input beagle file path/name (produced by PCAngsd)
#$2 #your desired output file path/name (will add k tested to end)

module load angsd
module load ngsTools

for i in {1..10};
    do NGSadmix -likes $1 -K ${SLURM_ARRAY_TASK_ID} -minMaf 0.05 -o $2_k${SLURM_ARRAY_TASK_ID}_run$i -P 10;
done
