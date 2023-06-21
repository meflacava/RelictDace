#!/bin/bash
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=index
#SBATCH --time=01:00:00
#SBATCH --partition=med
cd /home/mlacava/align/
module load bio3
bwa index GCA_022829085.1_ASM2282908v1_genomic.fna
