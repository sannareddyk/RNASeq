#!/bin/bash

REFERENCE_PATH="/mnt/home1/bioinfo/ks765/mm10_reference"
TRIMMED_DATA="/mnt/home1/bioinfo/ks765/xxx/fastq_trimmed"
ALIGNED_DATA="/mnt/home1/bioinfo/ks765/xxx/aligned_data"

for file in $TRIMMED_DATA/*_sickletrimmed_trimmed.fq.gz 
do
	outfile=$(basename "$file")
        outfile=${outfile%_sickletrimmed_trimmed.fq.gz}

	echo export SBATCH_CMD=\"STAR --genomeDir $REFERENCE_PATH --runThreadN 8 --readFilesIn $file --readFilesCommand zcat --quantMode GeneCounts --sjdbGTFfile $REFERENCE_PATH/Sabine_mm10_genes.gtf --outFileNamePrefix $ALIGNED_DATA/${outfile}. --outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outSAMunmapped Within\" >>~/cmd/SUBMIT_4

	echo sbatch ~/cmd/sbatch.sh >>~/cmd/SUBMIT_4

done
