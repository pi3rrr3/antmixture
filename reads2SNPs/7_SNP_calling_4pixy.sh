# 7_SNP_calling_4pixy.sh


###
### SNP CALLING & FILTERING -----------------------------------------------------------------------
###

{
#!/bin/bash -l
#SBATCH -J allSites2
#SBATCH -o /scratch/project_2001099/nouhaudp/reseq/snp/allSitesVCF/raw/logs/allSites-2_%a.out
#SBATCH -e /scratch/project_2001099/nouhaudp/reseq/snp/allSitesVCF/raw/logs/allSites-2_%a.err
#SBATCH --account=project_2001099
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-27
#SBATCH --ntasks 1
#SBATCH --mem=8G

# Load module and set wd
module load biokit
cd /scratch/project_2001099/nouhaudp/reseq/bam/data

# Get scaffold ID
REF=/scratch/project_2001099/nouhaudp/reseq/ref
scaffold=$(sed -n "$SLURM_ARRAY_TASK_ID"p $REF/scaffold.list)

# Write gvcf
bcftools mpileup -f $REF/Formica_hybrid_v1_wFhyb_Sapis.fa \
  -b all_fems.list \
 -r ${scaffold} | bcftools call -m -Oz -f GQ -o /scratch/project_2001099/nouhaudp/reseq/snp/allSitesVCF/raw/${scaffold}_allFems.vcf.gz

cd /scratch/project_2001099/nouhaudp/reseq/snp/allSitesVCF/raw

# Set DP thresholds & filter on missing data
mindp=354 # 6*59
maxdp=1770 # 30*59
bcftools view -Ou ${scaffold}_allFems.vcf.gz | bcftools filter -e "DP < ${mindp} | DP > ${maxdp}" --set-GTs . -Ou | bcftools filter -e "F_MISSING>0.5" -Oz > ${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz

# Create a filtered VCF containing only invariant sites
vcftools --gzvcf ${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz \
--max-maf 0 --recode --recode-INFO-all --stdout | bgzip -c > ${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz

# Create a filtered VCF containing only variant sites, keep SNPqual >= 30 and filter out HW excess
vcftools --gzvcf ${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc.vcf.gz --mac 1 --hwe 0.001 --minQ 30 --recode --recode-INFO-all --stdout | bgzip -c > ${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz

# Index both vcfs using tabix
tabix ${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz
tabix ${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz

# Combine the two VCFs using bcftools concat
bcftools concat \
--allow-overlaps \
${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc_invar.vcf.gz \
${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc_var.vcf.gz \
-Oz > ${scaffold}_allFems_filtered.vcf.gz

# Index for concat afterwards
tabix ${scaffold}_allFems_filtered.vcf.gz

# Remove temp files
rm ${scaffold}_allFems_meanDP${mindp}-${maxdp}_maxNA50perc*

# End bash script
}


# Combine VCFs
for file in *allFems_filtered.vcf.gz ; do tabix $file ; done

ls *allFems_filtered.vcf.gz > all_filtered_vcfs.list

bcftools concat \
--allow-overlaps \
-f all_filtered_vcfs.list \
-Oz > allFems_filtered.vcf.gz

tabix allFems_filtered.vcf.gz





###
### PIXY COMPUTATION ------------------------------------------------------------------------------
###


# create population list
cd /scratch/project_2001099/nouhaudp/reseq/snp/analyses/pixy
awk '$6 == "f" {print $0}' ../../../sample_table_R_ordered_LanFix.tab | cut -f1,15 | grep -v 'FA07_5q' > pop.list


#!/bin/bash -l
#SBATCH -J pixy
#SBATCH -o /scratch/project_2001099/nouhaudp/reseq/snp/analyses/pixy/pixy2.out
#SBATCH -e /scratch/project_2001099/nouhaudp/reseq/snp/analyses/pixy/pixy2.err
#SBATCH --account=project_2001099
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=4G

module load biokit
module load bioconda/3
source activate /scratch/project_2001099/nouhaudp/conda/gnmcs

cd /scratch/project_2001099/nouhaudp/reseq/snp/analyses/pixy

pixy --stats pi fst dxy \
--vcf /scratch/project_2001099/nouhaudp/reseq/snp/allSitesVCF/raw/allFems_filtered.vcf.gz \
--populations pop.list \
--bed_file windows_20kb.bed \
--output_prefix allFems \
--n_cores 8


# sort output - FST
head -n1 allFems_fst.txt > tmp
grep -v 'pop' allFems_fst.txt | sort -k3,3 -k4,4n > tmp2
cat tmp tmp2 > allFems_fst_sort.txt

# sort output - Dxy
head -n1 allFems_dxy.txt > tmp
grep -v 'pop' allFems_dxy.txt | sort -k3,3 -k4,4n > tmp2
cat tmp tmp2 > allFems_dxy_sort.txt

# sort output - Pi
head -n1 allFems_pi.txt > tmp
grep -v 'pop' allFems_pi.txt | sort -k2,2 -k3,3n > tmp2
cat tmp tmp2 > allFems_pi_sort.txt

rm tmp tmp2
