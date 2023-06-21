#!/bin/bash

#SBATCH -t 0-02:00:00
#SBATCH --job-name=div
#SBATCH -p high
#SBATCH --cpus-per-task=4
#SBATCH --mem=60GB
#SBATCH --output=div-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=div-%A-%a.err # File to which STDERR will be written
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --array=1-30

#Set SBATCH array range above based on number of pops in all.pop.list [e.g., 1-30 for 30 pops to test]
# wc -l all.pop.list #check how many sites are in your pop list

#Files to have in folder you run this job:
# - bamlist for each pop you want div stats for (named [pop code].bamlist, e.g., BTOC.bamlist)
# - all.pop.list = list of pop names, one per line (e.g., BTOC) that correspond to the bamlists in the same folder

#MODIFY your reference genome path and name below, otherwise run as is

#SUBMIT this job with: SBATCH get_div_all.sh

module load angsd

pop=$(sed -n ${SLURM_ARRAY_TASK_ID}p all.pop.list)
ref='/home/mlacava/align/GCA_022829085.1_ASM2282908v1_genomic.fna'


## Step 1: find a global estimate of the SFS

#Estimate the site allele frequency likelihood
# P=number of threads
# Use whatever quality scores you want (probably the same you used for PCAngsd)
angsd -bam ${pop}.bamlist -doSaf 1 -anc ${ref} -GL 1 -P 4 -minMapQ 20 -minQ 20 -out ${pop}
#output: .arg, .saf.gz, .saf.idx, .saf.pos.gz

#Obtain ML estimate of the SFS
# we are calculating a folded sfs because we do not know the ancestral state, and instead we are using the reference genome (if you know the ancestral state, provide it and exclude -fold 1 here and in next step)
realSFS ${pop}.saf.idx -maxIter 100 -P 4 -fold 1 > ${pop}.sfs


##Step 2: calculate thetas for each site

realSFS saf2theta ${pop}.saf.idx -sfs ${pop}.sfs -outname ${pop} -fold 1
#output: out.thetas.gz and out.thetas.idx

#extract logscale persite thetas
thetaStat print ${pop}.thetas.idx > ${pop}.per-site.thetas.idx.txt


##Step 3: estimate Tajima's D and other stats

thetaStat do_stat ${pop}.thetas.idx
#output: .thetas.idx.pestPG
