# 2_mapping.sh

#!/bin/bash -l
#SBATCH -J map
#SBATCH -o /scratch/project_2001099/nouhaudp/reseq/bam/logs/map_%j.out
#SBATCH -e /scratch/project_2001099/nouhaudp/reseq/bam/logs/map_%j.err
#SBATCH --account=project_2001099
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-60
#SBATCH --ntasks 8
#SBATCH --mem=12G


# load modules and set wd

module load biokit
module load picard/2.21.4
cd /scratch/project_2001099/nouhaudp/reseq

# Define paths

REFPATH=/scratch/project_2001099/nouhaudp/reseq/ref
FASTQPATH=/scratch/project_2001099/nouhaudp/reseq/fastq/trim
TMPBAMPATH=/scratch/project_2002161/nouhaudp/tmp
MYBAMPATH=/scratch/project_2001099/nouhaudp/reseq/bam


# Get file name and sample ID

file=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/project_2001099/nouhaudp/reseq/fastq/name.list)
shortfile=${file%_BD*}

echo "###### STARTING file: $file - sample ID: $shortfile"


### 1. Map, sort & index

bwa mem -t8 $REFPATH/Formica_hybrid_v1_wFhyb_Sapis.fa $FASTQPATH/$file"_1.trim.pair.fq.gz" $FASTQPATH/$file"_2.trim.pair.fq.gz" | samtools sort -@8 -m512M -o $TMPBAMPATH/${shortfile}".bam" -

samtools index -@4 $TMPBAMPATH/${shortfile}".bam"


### 2. Compute insert size distribution

java -Xmx4G -jar /appl/soft/bio/picard/picard-tools-2.21.4/picard.jar CollectInsertSizeMetrics \
I=$TMPBAMPATH/${shortfile}".bam" \
O=$MYBAMPATH/stats/${shortfile}"_insert_size_metrics.txt" \
H=$MYBAMPATH/stats/${shortfile}"_insert_size_hist.pdf"


### 3. Duplicate filtering

java -Xmx4G -jar /appl/soft/bio/picard/picard-tools-2.21.4/picard.jar MarkDuplicates \
I=$TMPBAMPATH/${shortfile}".bam" \
O=$MYBAMPATH/nodupl/${shortfile}"_nodupl.bam" \
M=$MYBAMPATH/nodupl/stats/${shortfile}"_dupl_metrics.txt" \
REMOVE_DUPLICATES=T TMP_DIR=$TMPBAMPATH


### 4. Index and compute stats

samtools index -@4 $MYBAMPATH/nodupl/${shortfile}"_nodupl.bam"
samtools flagstat -@4 $MYBAMPATH/nodupl/${shortfile}"_nodupl.bam" > $MYBAMPATH/nodupl/stats/${shortfile}"_nodupl.flagstat"

echo "###### DONE! file: $file - sample ID: $shortfile"

# samtools idxstats -@4 $BAMPATH/nodupl/${shortfile}"_nodupl.bam" > $BAMPATH/nodupl/stats/${shortfile}"_nodupl.idxstat"

# End of slurm script



# Collect insert size metrics
cd /scratch/project_2001099/nouhaudp/reseq/bam/stats
for i in *txt ; do echo $i ;  head -n8 $i | tail -n1 ; done


# Collect duplicate metrics
cd /scratch/project_2001099/nouhaudp/reseq/bam/stats
for i in *dupl_metrics.txt ; do echo $i ;  head -n8 $i | tail -n1 ; done
