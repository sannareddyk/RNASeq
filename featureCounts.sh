#!/bin/bash

#SBATCH -J feature_counts
#SBATCH --mail-type=END
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 0-24:00
#SBATCH -p xxx
#SBATCH --mail-user=xxx


BAMS_PATH='/mnt/home1/bioinfo/ks765/xxx/aligned_data'
OUT_PATH='/mnt/home1/bioinfo/ks765/xxx/featureCounts'
GTF='/mnt/b2/home1/bioinfo/ks765/hg19_BAM/Homo_sapiens.GRCh37.87.gtf'

for file in $BAMS_PATH/*.Aligned.sortedByCoord.out.bam
do
outfile=$(basename "$file")
outfile=${outfile%.Aligned.sortedByCoord.out.bam}
#echo $outfile

	featureCounts -T 2 -t exon -g gene_id -M -a $GTF -o $OUT_PATH/${outfile}.featureCounts $file

	tail -n +3 $OUT_PATH/${outfile}.featureCounts | cut -f1,6,7 | sort>$OUT_PATH/${outfile}.fcounts.txt

done
