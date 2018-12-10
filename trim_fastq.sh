#!/bin/bash

#SBATCH -J trim_fastq
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 0-24:00
#SBATCH -p IACT
#SBATCH --mail-user=xxx

#FASTQ_MERGED_PATH='/mnt/home1/bioinfo/ks765/xxx/xxx’
FASTQ_MERGED_PATH='/mnt/home1/bioinfo/ks765/xxx/fastq'
TRIMMED_DATA='/mnt/home1/bioinfo/ks765/xxx/fastq_trimmed'

#for file in $FASTQ_MERGED_PATH/*.sequence.txt.gz
for file in $FASTQ_MERGED_PATH/*R1_001.fastq.gz
do
	outfile=$(basename "$file")
	outfile=${outfile%_R1_001.fastq.gz}

	echo "Running sickle on merged fastq $file" >${outfile}.log
	sickle se -f $file -t sanger -n -g -o $TRIMMED_DATA/${outfile}_sickletrimmed.fq.gz &>>${outfile}.log

	echo "Now running trim_galore on $TRIMMED_DATA/${outfile}_sickletrimmed.fq.gz" >>${outfile}.log
	trim_galore -q 20 --length 20 --gzip --trim-n --stringency 3 --fastqc $TRIMMED_DATA/${outfile}_sickletrimmed.fq.gz -o $TRIMMED_DATA &>>${outfile}.log

done
