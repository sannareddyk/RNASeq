#!/bin/bash

#SBATCH -J fastqc
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-24:00
#SBATCH -p xxx
#SBATCH --mail-user=


FASTQS_PATH="/mnt/home1/bioinfo/ks765/xxx”
OUT_PATH="/mnt/home1/bioinfo/ks765/xxx/fastqc"

for file in $FASTQS_PATH/*_R1_001.fastq.gz
do

fastqc -o $OUT_PATH --noextract -f fastq -q $file

done
