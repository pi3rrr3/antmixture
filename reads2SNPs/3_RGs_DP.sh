# 3_RGs_DP.sh

###
### Add read groups -------------------------------------------------------------------------------------------------------
###


#!/bin/bash -l
#SBATCH -J add_RGs
#SBATCH -o /scratch/project_2001099/nouhaudp/reseq/bam/logs/add_RGs_%j.out
#SBATCH -e /scratch/project_2001099/nouhaudp/reseq/bam/logs/add_RGs_%j.err
#SBATCH --account=project_2001099
#SBATCH -t 2:00:00
#SBATCH -p small
#SBATCH --array=1-60
#SBATCH --ntasks 1
#SBATCH --mem-per-cpu=4G

module load biokit

FINALDIR=/scratch/project_2001099/nouhaudp/reseq/bam/nodupl_RG

cd /scratch/project_2001099/nouhaudp/reseq/bam/nodupl

# Get file
file=$(sed -n "$SLURM_ARRAY_TASK_ID"p /scratch/project_2001099/nouhaudp/reseq/fastq/name.list)

# Get sample ID
sample=${file%_BDSW*}

# Get library info (eg BDSW190006966-1a)
plcmd='s/.+(BDS.+\da).+/$1/'
lib=$(echo "$file" | perl -pe "$plcmd") ; echo $lib

# Get flowcell/lane info (eg HJNHWDSXX_L4)
plcmd='s/.+a_(.+)/$1/'
fcl=$(echo "$file" | perl -pe "$plcmd") ; echo $fcl

# Add RGs
java -Xmx4G -jar /appl/soft/bio/picard/picard-tools-2.21.4/picard.jar AddOrReplaceReadGroups \
    I=$sample"_nodupl.bam" \
    O=$FINALDIR/$sample"_nodupl_wRG.bam" \
    RGID=$lib"_"$fcl \
    RGPL=illumina \
    RGLB=$lib \
    RGSM=$sample \
    RGPU=$fcl \
    TMP_DIR=/scratch/project_2001099/nouhaudp/tmpdir



### End




###
### Compute sequencing depth per sample --------------------------------------------------------------------------------------------
###

cd /scratch/project_2001099/nouhaudp/reseq/bam/nodupl_RG
ls *bam > ../bam.list ; cd ..

#!/bin/bash -l
#SBATCH -J blckstr
#SBATCH -o /scratch/project_2001099/nouhaudp/reseq/bam/logs/mosdepth_%j.out
#SBATCH -e /scratch/project_2001099/nouhaudp/reseq/bam/logs/mosdepth_%j.err
#SBATCH --account=project_2001099
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-60
#SBATCH --ntasks 1
#SBATCH --mem=1G

cd /scratch/project_2001099/nouhaudp/reseq/bam/nodupl_RG

module load samtools
module load bioconda/3
source activate analyse

file=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../bam.list)
myprefix=${file%.bam}
samtools index $file
mosdepth -t 1 -n -x ../coverage/$myprefix $file



# Get results

cd coverage
rm -rf *global*

for i in *nodupl*.mosdepth.summary.txt
 do echo $i ; grep 'Scaffold' $i | awk '{sum+=$4;} END{print sum/NR;}'
done > all.mosdepth.summary.tmp

grep -v 'txt' all.mosdepth.summary.tmp > tmp1
grep 'txt' all.mosdepth.summary.tmp | perl -npe 's/_nodupl_wRG.mosdepth.summary.txt//' > tmp2
paste tmp2 tmp1 > all.mosdepth.summary.txt ; rm *tmp*

scp puhti:/scratch/project_2001099/nouhaudp/reseq/bam/coverage/all.mosdepth.summary.txt .

# in R
# options(stringsAsFactors=F)
# tt = read.table("all_mosdepth_summary_Fhyb_v1.txt", h = F)
# summary(tt$V2)
# overall mean: 14.2 (9.79 - 18.88)
