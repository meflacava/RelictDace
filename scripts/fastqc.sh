#!/bin/bash
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=fastqc
#SBATCH --time=12:00:00
#SBATCH --partition=med
#SBATCH --output=fastqc.out

# Run this script with sbatch fastqc.sh

#load module
module load fastqc/0.11.9

cd /home/mlacava/raw_data/RelictDace_BMAG070

fastqc /home/mlacava/raw_data/RelictDace_BMAG070/*.fastq
