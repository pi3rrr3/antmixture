# 6_phasing.sh


cd /scratch/project_2001099/nouhaudp/reseq/snp/phasing
BAMDIR=/scratch/project_2001099/nouhaudp/reseq/bam/data/
VCF=/scratch/project_2001099/nouhaudp/reseq/snp/vcf/SNP_A.FEMALES.minDP8.AN15percMiss.vcf.gz



###
#### Whatshap -------------------------------------------------------------------------------------
###

{
#!/bin/bash -l
#SBATCH -J whap
#SBATCH -o /scratch/project_2001099/nouhaudp/reseq/snp/phasing/logs/whatshap_%a.out
#SBATCH -e /scratch/project_2001099/nouhaudp/reseq/snp/phasing/logs/whatshap_%a.err
#SBATCH --account=project_2001099
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-60
#SBATCH --ntasks 1
#SBATCH --mem=2G

# conda mambo
module load biokit
module load bioconda/3
source activate /scratch/project_2001099/nouhaudp/conda/gnmcs

cd /scratch/project_2001099/nouhaudp/reseq/snp/phasing

# Define paths
BAMDIR=/scratch/project_2001099/nouhaudp/reseq/bam/data
masterVCF=/scratch/project_2001099/nouhaudp/reseq/snp/vcf/SNP_A.FEMALES.minDP8.AN15percMiss.vcf.gz
REF=/scratch/project_2001099/nouhaudp/reseq/ref/Formica_hybrid_v1_wFhyb_Sapis.fa

# Define focal ind
ind=$(sed -n "$SLURM_ARRAY_TASK_ID"p female.list)

# 1. Extract single individual from VCF & fix floats in GQ field (or issue later when parsing with whatshap)
bcftools view -s $ind $masterVCF -Ov | perl -npe 's/(.\/.)\:.+\:(.+\:.+\:.+\:.+\:.+\:.+\:.+)/$1\:99\:$2/' | grep -v '^Scaffold00' | bgzip > whatshap/${ind}.GQfixed.vcf.gz &&\
bcftools index -t whatshap/${ind}.GQfixed.vcf.gz

# 2. Phase
whatshap phase -o whatshap/${ind}.phased.vcf.gz --reference $REF \
               whatshap/${ind}.GQfixed.vcf.gz $BAMDIR/${ind}_nodupl_wRG.bam &&\
whatshap stats whatshap/${ind}.phased.vcf.gz --tsv=whatshap/${ind}.phased.tsv

### end

### Combine individual VCFs (with a fix for a GL issue)
#for i in *.phased.vcf.gz ; do echo $i ; perl -pi -e 's/(\#\#FORMAT=<ID=GL,Number=)G/$1\./' $i ; bcftools index -t $i ; done
for i in *.phased.vcf.gz ; do
  echo $i
  mv $i ${i}.tmp
  bcftools view ${i}.tmp -Ov | perl -pi -e 's/(\#\#FORMAT=<ID=GL,Number=)G/$1\./' | bgzip > $i
  bcftools index -t $i
done

bcftools merge -l female_phased_vcf.list -Oz > SNP_A.FEMALES.minDP8.AN15percMiss.whap.vcf.gz &&\
bcftools index -t SNP_A.FEMALES.minDP8.AN15percMiss.whap.vcf.gz
}





###
#### ShapeIt --------------------------------------------------------------------------------------
###

{
cd /scratch/project_2001099/nouhaudp/reseq/snp/phasing

#!/bin/bash -l
#SBATCH -J shpt4
#SBATCH -o /scratch/project_2001099/nouhaudp/reseq/snp/phasing/logs/shpt4_Scaffold%a.out
#SBATCH -e /scratch/project_2001099/nouhaudp/reseq/snp/phasing/logs/shpt4_Scaffold%a.err
#SBATCH --account=project_2001099
#SBATCH -t 12:00:00
#SBATCH -p small
#SBATCH --array=1-27
#SBATCH --ntasks 4
#SBATCH --mem=2G

cd /scratch/project_2001099/nouhaudp/reseq/snp/phasing/shapeit

mychr=$(sed -n "$SLURM_ARRAY_TASK_ID"p scaffold.list)

module load bioconda/3
source activate /scratch/project_2001099/nouhaudp/conda/gnmcs

shapeit4 --input ../whatshap/SNP_A.FEMALES.minDP8.AN15percMiss.whap.vcf.gz \
  --output "SNP_A.FEMALES.minDP8.AN15percMiss.whap.shapeit."$mychr".vcf.gz" \
  --region $mychr --sequencing --thread 4 \
  --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,1b,1p,1b,1p,10m --pbwt-depth 8 --use-PS 0.0001

### End of bash script

}


## Combine all chromosomes

ls SNP_A.*.shapeit*vcf.gz > file.list
bcftools concat -Oz -f file.list -o SNP_A.FEMALES.minDP8.AN15percMiss.whap.shapeit.allScafs.vcf.gz

# Remove temp files
rm -rf *it.Scaffold*vcf*
rm -rf ../whatshap/*vcf*


#### This is the end.
