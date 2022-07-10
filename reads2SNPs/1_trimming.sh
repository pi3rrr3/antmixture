# 1_trimming.sh

#!/bin/bash -l
#SBATCH -J trim
#SBATCH -o logs/trimming/trim_%j.out
#SBATCH -e logs/trimming/trim_%j.err
#SBATCH --account=project_2001099
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-59
#SBATCH --ntasks 8
#SBATCH --mem=4G

cd /scratch/project_2001099/nouhaudp/reseq/fastq

module load trimmomatic/0.39

# Get file name
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p name.list)

# Get sample ID, for adapter removal
sample=${file%_BDS*}

# Trim reads
trimmomatic PE -threads 8 -phred33 raw/$file"_1.fq.gz" raw/$file"_2.fq.gz" \
-trimlog trim/logs/$file".trim.log" \
trim/$file"_1.trim.pair.fq.gz" trim/unpair/$file"_1.trim.unpair.fq.gz" trim/$file"_2.trim.pair.fq.gz" trim/unpair/$file"_2.trim.unpair.fq.gz" \
ILLUMINACLIP:illumina_raw_data/adapter/$sample/${sample}_all_adapters.fa:2:30:10:2:keepBothReads LEADING:20 TRAILING:20 MINLEN:50


