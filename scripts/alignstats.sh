#!/bin/bash -l
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=stats
#SBATCH --error=/home/mlacava/slurm_logs/align.%j.err
#SBATCH --output=/home/mlacava/slurm_logs/align.%j.out
#SBATCH --cpus-per-task=4
#SBATCH --partition=high
#SBATCH --time=01:00:00

###This script was written by Shannon Joslin and modified by Shannon Kieran and Melanie LaCava##
##command line arguments: 1. align directory name, 2. output stats file name
set -e # exits upon failing command
set -v # verbose -- all lines
#set -x # trace of all commands after expansion before execution


# set up directories
RAD_dir=$1

#############################
###  create master stats  ###
#############################

##This just concatonates your stats and calculates some percentages. Use it to pick individuals with more than 50K reads to use in downstream analysis##
cd ${RAD_dir}
mkdir -p stats #makes a subdirectory called "stats"
cp *.stats stats/.
cd ${RAD_dir}/stats
touch STATS.all
touch STATNAMES.all
for i in *.stats
do
echo "$(head -n 1 $i)" >> STATS.all
echo "$i" >> STATNAMES.all
done

pr -mts, STATNAMES.all STATS.all > statsummary.all

sed -i 's/\.stats//g' statsummary.all
file="${RAD_dir}/stats/statsummary.all"
touch $2
n=$(wc -l ${file} | awk '{print $1}')
x=1

while [ $x -le $n ]
do
    string="sed -n ${x}p $file"
        str=$($string)
    nam=$(echo $str | awk -F"," '{print $1}')
    aln=$(echo $str | awk -F"," '{print $2}')
    ppd=$(echo $str | awk -F"," '{print $3}')
    rmd=$(echo $str | awk -F"," '{print $4}')
    p2=$(echo "scale=3 ; $ppd / $aln" | bc)
    p3=$(echo "scale=3 ; $rmd / $ppd" | bc)
    p4=$(echo "scale=3 ; $rmd / $aln" | bc)
    echo "${nam},${aln},${ppd},${rmd},${p2},${p3},${p4}" >> $2
    x=$(( $x + 1 ))
done
