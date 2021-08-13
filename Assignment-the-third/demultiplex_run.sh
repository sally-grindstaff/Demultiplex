#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --time=1-0:00:00
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mail-user='sgrindst@uoregon.edu'
#SBATCH --mail-type=BEGIN,END,FAIL

conda activate bgmp_py39

/usr/bin/time -v ./demultiplex.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -b /projects/bgmp/shared/2017_sequencing/indexes.txt -o ../output/

cd ../output

files=$(ls -1 | grep -v 'stats.txt' | grep -v 'stats.tsv' | grep -v '.gz')
for file in $files;
do
gzip -f $file;
done

exit