#!/bin/bash
#SBATCH --mail-user=mlacava@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=download
#SBATCH --time=3-00:00:00
#SBATCH --partition=med
#SBATCH --output=download_data.out
#SBATCH --error=download_data.err

# Run with sbatch download_data.sh

cd /home/mlacava/raw_data/RelictDace_BMAG070
echo "starting at " date

#download
wget -r -nH -nc -np -R index.html*    http://slimsdata.genomecenter.ucdavis.edu/Data/qic3hmnk9t/Unaligned/Project_ASMB_BMAG070_Run1/
date

wget -r -nH -nc -np -R index.html*    http://slimsdata.genomecenter.ucdavis.edu/Data/zn6isc6xo3/Unaligned/Project_ASMB_BMAG070_Run2/
date

wget -r -nH -nc -np -R index.html*    http://slimsdata.genomecenter.ucdavis.edu/Data/spfehf5o6h/Unaligned/Project_ASMB_BMAG070_Run3/
date
