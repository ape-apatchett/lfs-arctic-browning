#!/bin/bash -l

#SBATCH -A xxx
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -t 02:00:00


# Load required modules
module load bowtie2/2.5.4 samtools/1.20

fastqdir=/xxx/trimmomatic/trimmed_fastq
HG38_INDEX="/sw/data/bowtie_data/prebuilt/GRCh38_noalt_as/GRCh38_noalt_as"
OUTDIR=/xxx/trimmomatic/human_filtered_browning
mkdir -p $OUTDIR

#Loop over all R1 fastq files
for R1 in $fastqdir/*_combined_R1.fastq; do

SAMPLE=$(basename $R1 | sed 's/_combined_R1.*//')
R2="${fastqdir}/${SAMPLE}_combined_R2.fastq"
BAM="$OUTDIR/${SAMPLE}_hg38.bam"
UNMAPPED_R1="$OUTDIR/${SAMPLE}_R1_nohuman.fastq"
UNMAPPED_R2="$OUTDIR/${SAMPLE}_R2_nohuman.fastq"

echo "Processing sample: $SAMPLE"

bowtie2 -x $HG38_INDEX \
 -1 $R1 -2 $R2 \
 --threads $SLURM_CPUS_PER_TASK \
 | samtools sort -@ $SLURM_CPUS_PER_TASK -o $BAM

samtools index $BAM

MAPPED=$(samtools view -c -F 4 $BAM)
echo "$SAMPLE mapped reads: $MAPPED"

samtools view -b -f 12 $BAM \
 | samtools sort -n -@ $SLURM_CPUS_PER_TASK -o $OUTDIR/${SAMPLE}_unmapped.bam

samtools fastq -1 $UNMAPPED_R1 -2 $UNMAPPED_R2 -0 /dev/null -s /dev/null \
 -n $OUTDIR/${SAMPLE}_unmapped.bam

done


