#!/bin/bash -l

#SBATCH -A xxx
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G
#SBATCH -t 16:00:00

module load bioinfo-tools eggNOG-mapper/2.1.9 eggNOG_data/5.0.0 diamond/2.1.6 seqtk/1.4

# Define paths

fastqdir=/xxx/trimmomatic/human_filtered_browning
outdir=/xxx/trimmomatic/eggnog_out_nohuman_browning
logdir=/xxx/trimmomatic/eggnog_logs_browning
tmpdir=$SNIC_TMP


mkdir -p $outdir
mkdir -p $logdir
mkdir -p $tmpdir
 

# Loop over all R1 files
for R1 in $(ls $fastqdir/*_R1_nohuman.fastq | sort); do
s=$(basename "$R1" _R1_nohuman.fastq)
R2="${R1/_R1_/_R2_}"
echo "[$(date)] Running EggNOG-mapper for $s"

# Convert fastq to fasta
F1="$tmpdir/${s}_R1.fasta"
F2="$tmpdir/${s}_R2.fasta"
seqtk seq -a "$R1" > "$F1"
seqtk seq -a "$R2" > "$F2"

# Merge paired-end fasta temporarily
merged_fasta="$tmpdir/${s}_merged.fasta"
cat "$F1" "$F2" > "$merged_fasta"

# Run EggNOG
emapper.py \
 -i "$merged_fasta" \
 --itype metagenome \
 --translate \
 -m diamond \
 --cpu $SLURM_CPUS_PER_TASK \
 --output $outdir/${s}_eggnog \
 --data_dir $EGGNOG_DATA_ROOT

# Clean up merged fastq
rm -f "$F1" "$F2" "$merged_fasta"

echo "[$(date)] Done $s"
done
