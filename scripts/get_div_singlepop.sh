#!/bin/bash

#SBATCH -t 0-02:00:00
#SBATCH --job-name=div
#SBATCH -p high
#SBATCH --cpus-per-task=4
#SBATCH --mem=60GB
#SBATCH --output=div-%A.out # File to which STDOUT will be written
#SBATCH --error=div-%A.err # File to which STDERR will be written
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mlacava@ucdavis.edu

#Submit this job with: sbatch get_div_singlepop.sh [pop bamlist] [pop code]

#MODIFY your reference genome path and name below, otherwise run as is

#SUBMIT this job with: sbatch get_div_singlepop.sh [pop bamlist] [pop code]
#Example: sbatch get_div_singlepop.sh BTOC.bamlist BTOC

module load angsd

bamlist=$1
out=$2
ref='/home/mlacava/align/GCA_022829085.1_ASM2282908v1_genomic.fna'


## Step 1: find a global estimate of the SFS

#Estimate the site allele frequency likelihood
# P=number of threads
# Use whatever quality scores you want (probably the same you used for PCAngsd)
angsd -bam ${bamlist} -doSaf 1 -anc ${ref} -GL 1 -P 4 -minMapQ 20 -minQ 20 -out ${out}
#output: .arg, .saf.gz, .saf.idx, .saf.pos.gz

#Obtain ML estimate of the SFS
# we are calculating a folded sfs because we do not know the ancestral state, and instead we are using the reference genome (if you know the ancestral state, provide it and exclude -fold 1 here and in next step)
realSFS ${out}.saf.idx -maxIter 100 -P 4 -fold 1 > ${out}.sfs


##Step 2: calculate thetas for each site

realSFS saf2theta ${out}.saf.idx -sfs ${out}.sfs -outname ${out} -fold 1
#output: out.thetas.gz and out.thetas.idx

#extract logscale persite thetas
thetaStat print ${out}.thetas.idx > ${out}.per-site.thetas.idx.txt


##Step 3: estimate Tajima's D and other stats

thetaStat do_stat ${out}.thetas.idx
#output: .thetas.idx.pestPG
