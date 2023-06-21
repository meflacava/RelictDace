#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH --job-name=fst
#SBATCH -p high
#SBATCH --output=fst-%A-%a.out # File to which STDOUT will be written
#SBATCH --error=fst-%A-%a.err # File to which STDERR will be written
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mlacava@ucdavis.edu

#RUN get_div_all.sh or get_div_singlepop.sh before running this script

#MODIFY code below in 2 ways:
# - provide path to your existing .saf.idx files from diversity stats
# - insert your username to submit each job under your account

#SUBMIT with: sbatch get_fst.sh [pop list]
# - pop list is the list of site codes you want to calc Fst among (this script will loop throgh all possible pairs of these sites and create and submit an SBATCH job to run Fst for them) (e.g., all.pop.list)

module load angsd

infile=$1 #your pop list (e.g., all.pop.list)
n=$(wc -l $infile | awk '{print $1}')

x=1 
while [ $x -le $n ]
do
    y=$(( $x + 1 ))
    while [ $y -le $n ]
    do
    
    pop1=$( (sed -n ${x}p $infile) )
    pop2=$( (sed -n ${y}p $infile) )

        echo "#!/bin/bash" > ${pop1}.${pop2}.sh
        echo "" >> ${pop1}.${pop2}.sh
        #change relative or absolute path to your existing .saf.idx files below
        echo "realSFS ../diversity/${pop1}.saf.idx ../diversity/${pop2}.saf.idx > ${pop1}.${pop2}.2dsfs" >> ${pop1}.${pop2}.sh
        echo "" >> ${pop1}.${pop2}.sh
        echo "realSFS fst index ../diversity/${pop1}.saf.idx ../diversity/${pop2}.saf.idx -sfs ${pop1}.${pop2}.2dsfs -fstout ${pop1}.${pop2}" >> ${pop1}.${pop2}.sh
        echo "" >> ${pop1}.${pop2}.sh
        echo "realSFS fst stats ${pop1}.${pop2}.fst.idx 2> ${pop1}.${pop2}_global.fst" >> ${pop1}.${pop2}.sh

        #insert your FARM username below to submit sbatch job, increase time requested if your jobs don't finish
        sbatch -J mlacava -p med --mem=32G -t 06:00:00 -c 1 ${pop1}.${pop2}.sh

    y=$(( $y + 1 ))
    
    done

x=$(( $x + 1 ))

done
