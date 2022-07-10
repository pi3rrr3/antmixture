# 6_SNP_calling.sh


cd /scratch/project_2001099/nouhaudp/reseq/bam/nodupl_RG


# Split reference genome in 50 kb regions for parallelization purposes
source activate freebayes
fasta_generate_regions.py /scratch/project_2001099/nouhaudp/reseq/ref/Formica_hybrid_v1_wFhyb_Sapis.fa.fai 50000 > /scratch/project_2001099/nouhaudp/reseq/ref/Formica_hybrid_v1_wFhyb_Sapis_50kb_regions.txt



###
### Full SNP calling
###

module load bioconda
source activate freebayes # version:  v1.3.1

cd /scratch/project_2001099/nouhaudp/reseq/snp
RSQ=/scratch/project_2001099/nouhaudp/reseq

ls $RSQ/bam/nodupl_RG/*bam > bam.list

/appl/soft/bio/bioconda/miniconda3/envs/freebayes/bin/freebayes-puhti \
  -time 72 \
  -regions $RSQ/ref/Formica_hybrid_v1_wFhyb_Sapis_50kb_regions.txt \
  -f $RSQ/ref/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -L bam.list \
  -k --genotype-qualities \
  -out $RSQ/snp/vcf/all_samples_raw.vcf



module load biokit

# Compress & sort

#!/bin/bash -l
#SBATCH -J cmprsrt
#SBATCH -o /scratch/project_2001099/nouhaudp/reseq/snp/vcf/compress_sort.out
#SBATCH -e /scratch/project_2001099/nouhaudp/reseq/snp/vcf/compress_sort.err
#SBATCH --account=project_2001099
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 1
#SBATCH --mem=8G

module load biokit
cd /scratch/project_2001099/nouhaudp/reseq/snp/vcf

# compress
/appl/soft/bio/samtools/htslib-1.9/bgzip all_samples_raw.vcf

# sort
bcftools sort -m 1G -O z -o all_samples.vcf.gz -T /scratch/project_2001099/nouhaudp/reseq/snp/vcf/tmp_sort all_samples_raw.vcf.gz

# index
/appl/soft/bio/samtools/htslib-1.9/tabix -p vcf all_samples.vcf.gz
bcftools index -n all_samples.vcf.gz

### END
