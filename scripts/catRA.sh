#!/bin/bash
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=catA
#SBATCH --time=24:00:00
#SBATCH --partition=med

#concatenate reads across sequencing lanes 
#this script is for combining forward reads across lanes 1, 2, and 8 
for i in `cat ~/raw_data/RelictDace_BMAG070/names`; do cat ~/raw_data/RelictDace_BMAG070/demulti1/inds/$i\_RA* ~/raw_data/RelictDace_BMAG070/demulti2/inds/$i\_RA* ~/raw_data/RelictDace_BMAG070/demulti8/inds/$i\_RA* > ~/raw_data/RelictDace_BMAG070/compiled/$i\_RA.fastq; done
