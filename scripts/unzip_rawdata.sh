#!/bin/bash
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=unzip
#SBATCH --time=12:00:00
#SBATCH --partition=med

# Run this script with sbatch unzip_rawdata.sh

cd /home/mlacava/raw_data/RelictDace_BMAG070

gzip -drk /home/mlacava/raw_data/RelictDace_BMAG070
