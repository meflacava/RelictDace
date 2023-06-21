#!/bin/bash -l
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=align
#SBATCH --error=/home/mlacava/slurm_logs/align.%j.err
#SBATCH --output=/home/mlacava/slurm_logs/align.%j.out
#SBATCH --cpus-per-task=4
#SBATCH --partition=med
#SBATCH --time=1:00:00

##This script was written by Shannon Joslin, modified by Shannon Kieran and Melanie LaCava##

##On command line, specify 1.Source fastqs directory (full path), 2.Alignment directory (full path of where reference genome is and where you want alignments to go), 3.Reference filename (no path)

#BELOW: can edit amount of time to request for individual SBATCH scripts (set at 24 hours)

set -e # exits upon failing command
set -v # verbose -- all lines
#set -x # trace of all commands after expansion before execution


# set up directories ##
data_dir=$1 ##where you keep your demultiplexed, renamed reads##
align_dir=$2 ##where you want to keep your aligned reads directory##
ref=$3 ##name of reference genome (that has already been BWA indexed)

## NOTE: Need list of all sequence file names in a .list file.
#     Need ID'd loci

####################
###  makeRADseqID .list  ###
##this makes a list of all your RAs and RBs and puts it in your alignment directory##
####################
cd ${data_dir}

ls *_RA.fastq > id1
paste id1 id1 > id2
sed -i 's/_RA.fastq//g' id2
mv id2  ${align_dir}/RADseqID.list ##rename this file, and be sure to use sed or vim's %s/oldname/newname/g to change it in all the spots it pops up later##
rm id1


#######################
###  align RAD seq  ###
#######################

cd $align_dir
mkdir -p RAD_alignments
##this parses the list of individuals and creates an alignment script for each one, and starts it running on medium.
wc=$(wc -l RADseqID.list | awk '{print $1}')
x=1
while [[ $x -le $wc ]]
do
    string="sed -n ${x}p RADseqID.list"
    str=$($string)

        var=$(echo $str | awk -F"\t" '{print $1}')
        set -- $var
        c1=$1
    c2=$2
    echo "name is $1 path is $2"
        echo "#!/bin/bash -l" >> aln_${c1}.sh #make a script to align each sample
        echo "#SBATCH -o RAD_alignments/${c1}-%j.out" >> aln_${c1}.sh
    echo "#SBATCH -p med" >> aln_${c1}.sh
    echo "#SBATCH --time=24:00:00" >> aln_${c1}.sh
    echo "making aln_${c1}.sh"
        echo "" >> aln_${c1}.sh
    echo "cd ${data_dir}" >> aln_${c1}.sh
    echo "bwa mem ${align_dir}/${ref} ${data_dir}/${c2}_RA.fastq ${data_dir}/${c2}_RB.fastq > ${align_dir}/RAD_alignments/${c1}.sam"  >> aln_${c1}.sh #aligns the reads and outputs a sam file##
    echo "samtools view -bS ${align_dir}/RAD_alignments/${c1}.sam > ${align_dir}/RAD_alignments/${c1}.bam" >> aln_${c1}.sh #outputs a compressed (bam)file
        echo "samtools sort ${align_dir}/RAD_alignments/${c1}.bam -o ${align_dir}/RAD_alignments/${c1}_sorted.bam" >> aln_${c1}.sh #sorts the bamfile
    echo "samtools view -b -f 0x2 ${align_dir}/RAD_alignments/${c1}_sorted.bam > ${align_dir}/RAD_alignments/${c1}_sorted_proper.bam" >> aln_${c1}.sh #filters reads that aren't properly paired (F/R)
        echo "samtools rmdup ${align_dir}/RAD_alignments/${c1}_sorted_proper.bam ${align_dir}/RAD_alignments/${c1}_sorted_proper_rmdup.bam" >> aln_${c1}.sh #removes PCR duplicates#
        echo "sleep 2m" >> aln_${c1}.sh #waits for this to finish
        echo "samtools index ${align_dir}/RAD_alignments/${c1}_sorted_proper_rmdup.bam ${align_dir}/RAD_alignments/${c1}_sorted_proper_rmdup.bam.bai" >> aln_${c1}.sh #indexes bamfil
        echo "reads=\$(samtools view -c ${align_dir}/RAD_alignments/${c1}_sorted.bam)" >> aln_${c1}.sh #calculates reads for the aligned file
        echo "ppalign=\$(samtools view -c ${align_dir}/RAD_alignments/${c1}_sorted_proper.bam)" >> aln_${c1}.sh #calculates reads for the paired-filtered file
        echo "rmdup=\$(samtools view -c ${align_dir}/RAD_alignments/${c1}_sorted_proper_rmdup.bam)" >> aln_${c1}.sh #calculates reads for the PCR-dup-filtered file.
        echo "echo \"\${reads},\${ppalign},\${rmdup}\" > ${align_dir}/RAD_alignments/${c1}.stats" >> aln_${c1}.sh #outputs a stats file
        sbatch -J ${c1}_radaln aln_${c1}.sh #starts the job running
        rm aln_${c1}.sh #removes the alignment script, because they are clutter. Keep them if you think you might want them for troubleshooting, but it's almost always easier to make the scripts again then to change them all one by one.

        x=$(( $x + 1 ))

done
