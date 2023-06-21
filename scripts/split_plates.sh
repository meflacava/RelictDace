#!/bin/bash
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=splitp8
#SBATCH --time=1-00:00:00
#SBATCH --partition=med
#SBATCH --output=splitp8.out

# Run this script with sbatch split_plates.sh

mkdir /scratch/$SLURM_JOBID
cd /scratch/$SLURM_JOBID

/home/maccamp/miniconda2/bin/fastq-multx -m 0 -B /home/mlacava/raw_data/RelictDace_BMAG070/barcodes_RD.tsv \
/home/mlacava/raw_data/RelictDace_BMAG070/Undetermined_S0_L008_I1_001.fastq.gz \
/home/mlacava/raw_data/RelictDace_BMAG070/Undetermined_S0_L008_R1_001.fastq.gz \
/home/mlacava/raw_data/RelictDace_BMAG070/Undetermined_S0_L008_R2_001.fastq.gz \
-o n/a \
-o ./%_R1.fastq \
-o ./%_R2.fastq

mv ./*.fastq /home/mlacava/raw_data/RelictDace_BMAG070/demulti8/

cd ..

rm -r $SLURM_JOBID
