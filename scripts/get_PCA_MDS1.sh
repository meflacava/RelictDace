#!/bin/bash
#SBATCH --time=2-00
#SBATCH --mem=50G# Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=high # Partition to submit to
#SBATCH --output=pca-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=pca-%A-%a.err # File to which STDERR will be written
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mlacava@ucdavis.edu

#Submit job with sbatch get_PCA_MDS1.sh [path/to/bamlist] [path/to/chromosome/list.txt] [path/to/desired/output]

bamlist=$1 #provide name of bamlist with full path
ref=$2 #provide name of reference genome file with full path
out=$3
nInd=$(wc $bamlist | awk '{print $1}')
#minInd=$[$nInd*9/10]
minInd=1313

#mfreq=$4
module load angsd

angsd -bam ${bamlist} -out ${out} -doIBS 1 -rf ${ref} -doGLF 2 -doCounts 1 -doMajorMinor 1 -minMapQ 20 -minQ 20 -SNP_pval 1e-3 -minMaf 0.05 -makeMatrix 1 -doCov 1 -GL 1 -doMaf 1 -doPost 2 -minInd $minInd
