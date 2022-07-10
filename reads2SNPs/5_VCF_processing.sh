# 5_VCF_processing.sh

#!/bin/bash -l
#SBATCH -J proc
#SBATCH --account=project_2001099
#SBATCH -t 48:00:00
#SBATCH -p small
#SBATCH --ntasks 8
#SBATCH --mem=12G


module load biokit
module load bioconda/3
source activate /scratch/project_2001099/nouhaudp/conda/gnmcs

cd /scratch/project_2001099/nouhaudp/reseq/snp/vcf/filtering

### VCF normalization
vt normalize -n -r ../../../ref/Formica_hybrid_v1_wFhyb_Sapis.fa ../all_samples.vcf.gz | bgzip -c > all_samples.normalized.vcf.gz


##
## 0. Pre-filtering: indels, non-SNPs, read imbalance, decomposition ---------------
##

bcftools filter --threads 8 -Oz -s+ --SnpGap 2 all_samples.normalized.vcf.gz > all_samples.normalized.SnpGap_2.vcf.gz && \

bcftools filter --threads 8 -Oz -e 'TYPE!="snp"' -s NonSnp -m+ all_samples.normalized.SnpGap_2.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.vcf.gz && \

bcftools filter --threads 8 -Oz -s Balance -m+ -i 'RPL>=1 && RPR>=1 & SAF>=1 && SAR>=1' all_samples.normalized.SnpGap_2.NonSNP.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.vcf.gz && \

bcftools view --threads 8 -O z -f PASS all_samples.normalized.SnpGap_2.NonSNP.Balance.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.vcf.gz && \

bcftools view --threads 8 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.vcf.gz | vcfallelicprimitives --keep-info --keep-geno -t decomposed | sed '/^##/! s/|/\//g' | sed 's/\.:\.:\.:\.:\.:\.:\.:\./\.\/\.:\.:\.:\.:\.:\.:\.:\./g' | bcftools sort --temp-dir /scratch/project_2001099 --max-mem 4G -O z > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz && \

bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz
echo "bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz"
bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz



##
## 1. SNP QUAL >= 30, biallelic -------------------------------------------------------------------
##

bcftools filter --threads 8 --include 'QUAL >= 30 && TYPE="snp"' -Oz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz
echo "gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'"
gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz | grep -vc '#'

bcftools view --threads 18 --min-alleles 2 --max-alleles 2 all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.vcf.gz -Oz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz
echo "gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'"
gunzip -c all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz | grep -vc '#'
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz



##
## 2. Nuclear genome only -------------------------------------------------------------------------
##

bcftools view -R scaffold.tab all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.vcf.gz -Oz > all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.vcf.gz && \
bcftools index -t all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.vcf.gz

bcftools index -n all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.vcf.gz



##
## 3. DP filtering --------------------------------------------------------------------------------
##

### First, correct field IDs (MNV coding issue not replaced in header)
# Extract header from VCF
module load biokit
bcftools view -h all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.vcf.gz > header.vcf

# Fix fields
perl -npe 's/<ID=AO,Number=A/<ID=AO,Number=\./' header.vcf | perl -npe 's/<ID=AD,Number=R/<ID=AD,Number=\./' | perl -npe 's/<ID=QA,Number=A/<ID=QA,Number=\./' | perl -npe 's/<ID=GL,Number=G/<ID=GL,Number=\./' > header_AO_AD_QA_GL.vcf

# Replace corrected header
bcftools reheader -h header_AO_AD_QA_GL.vcf -o all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.vcf.gz all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.vcf.gz


### Compute individual DP thresholds

vcftools --depth --gzvcf all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.vcf.gz --out all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader

# in R
{

  dat = read.table("all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.idepth", h = T)
  dat$MIN_DEPTH = round((dat$MEAN_DEPTH-.5)/2) # min depth is 50% of mean depth
  dat$MAX_DEPTH = round((dat$MEAN_DEPTH)*2) # max depth is twice mean depth
  write.table(dat, "all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.idepth.MinMax", col.names = T, row.names = F, quote = F, sep = "\t")

}



#!/bin/bash -l
#SBATCH -J splitDP
#SBATCH --account=project_2001099
#SBATCH -o logs/splitDP_%a.out
#SBATCH -e logs/splitDP_%a.err
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-60
#SBATCH --ntasks 1
#SBATCH --mem=2G

module load biokit
cd /scratch/project_2001099/nouhaudp/reseq/snp/vcf/filtering

# Get file name and sample ID
ind=$(sed -n "$SLURM_ARRAY_TASK_ID"p ind.list)
mkdir individualDP/${ind}
cd individualDP/${ind}

# Extract single individual from VCF and apply DP filters
dpmin=2 # not relevant yet, will be increased afterwards
dpmax=$(grep $ind ../../all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.idepth.MinMax | cut -f5)
echo "Processing individual $ind with $dpmin < DP <= $dpmax"
bcftools view -s ${ind} -Ou ../../all_samples.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.vcf.gz | bcftools filter -e "FORMAT/DP < ${dpmin} | FORMAT/DP >= ${dpmax}" --set-GTs . -Oz > ${ind}.vcf.gz &&
bcftools index -t ${ind}.vcf.gz





##
## 4. HWE testing ---------------------------------------------------------------------------------
##


# Merge all individual VCFs
ls */*gz | grep -v 'm.' > female.list
bcftools merge -l female.list -Oz --threads 1 > ../females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.indDP.vcf.gz


### Test HWE in females
vcftools --gzvcf females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.indDP.vcf.gz --hardy --out FEMALES_hwe

### Extract bad sites, in R
Rscript hwe_filtering.R
{
  # # hwe_filtering.R
  # # get the filename to read the hwe test
  # filename <- "FEMALES_hwe" # file name of reference table
  # # print the command line options
  # print(paste("file name with hwe output:", filename))
  #
  # # Read the HW results
  # hwvalues <- read.table(paste(filename,".hwe",sep=""), header=T, stringsAsFactors = F)
  # # Get distribution of the p-values
  # print(paste("he excess", sum(hwvalues$P_HET_EXCESS<0.01)))
  # print(paste("he deficit", sum(hwvalues$P_HET_DEFICIT<0.01)))
  # print(paste("he overall", sum(hwvalues$P_HWE<0.01)))
  #
  # # remove only the sites with excess of heterozygotes
  # indextoremove <- which(hwvalues$P_HET_EXCESS<0.01)
  # # Create BED file with the sites that fail the HW test for excess of heterozygotes
  # # BED file has three entries:
  # # chromosome, position-1, position (it is position-1 because it is assumed that the first base is 0)
  # position <- hwvalues[indextoremove,2]
  # bedmatrix <- matrix(c(hwvalues[indextoremove,1],format(position-1, scientific = F),format(position, scientific = F)), ncol=3, byrow=F)
  # write("#Sites with heterozygosity excess", file=paste("nohwe_excess_",filename,".bed",sep=""))
  # write(t(bedmatrix), file=paste("nohwe_excess_",filename,".bed",sep=""), ncolumns = 3, append=T)
}

### Remove bad sites in both sex-specific VCFs

vcftools --gzvcf females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.indDP.vcf.gz --exclude-bed nohwe_excess_FEMALES_hwe.bed --recode --recode-INFO-all --stdout | bgzip >  females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.indDP.hwe.vcf.gz





##
## 5. Min DP filtering  ---------------------------------------------------------------------------
##

#!/bin/bash -l
#SBATCH -J DP_miss
#SBATCH --account=project_2001099
#SBATCH -o logs/DP_missingData.out
#SBATCH -e logs/DP_missingData.err
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --ntasks 2
#SBATCH --mem=4G

module load biokit

cd /scratch/project_2001099/nouhaudp/reseq/snp/vcf/filtering/

### Females

# DP 5
dpmin=5
bcftools filter --threads 1 -i "FORMAT/DP>=${dpmin}" --set-GTs . -Oz females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.indDP.hwe.vcf.gz > females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.minDP${dpmin}.hwe.vcf.gz
bcftools index -t females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.minDP${dpmin}.hwe.vcf.gz

# DP 8
dpmin=8
bcftools filter --threads 1 -i "FORMAT/DP>=${dpmin}" --set-GTs . -Oz females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.indDP.hwe.vcf.gz > females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.minDP${dpmin}.hwe.vcf.gz
bcftools index -t females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.minDP${dpmin}.hwe.vcf.gz





##
## 6. Remove sites with missing data --------------------------------------------------------------
##

# min DP 8 for females, 15% missing data (60 individuals, 120 chromosomes)
inputfile=females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.minDP8.hwe.vcf.gz
bcftools view --threads 1 -Oz -i 'AN >= 102' ${inputfile} > ${inputfile%.vcf.gz}.AN15percMiss.vcf.gz
bcftools index -t ${inputfile%.vcf.gz}.AN15percMiss.vcf.gz
bcftools index -n ${inputfile%.vcf.gz}.AN15percMiss.vcf.gz
vcftools --gzvcf ${inputfile%.vcf.gz}.AN15percMiss.vcf.gz --missing-indv --out ${inputfile%.vcf.gz}.AN15percMiss.MISSING_DATA

 cp temp_files/females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.minDP8.hwe.AN15percMiss.vcf.gz ../SNP_A.FEMALES.minDP8.AN15percMiss.vcf.gz


 # min DP 8 for females, 25% missing data (60 individuals, 120 chromosomes)
 inputfile=females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.minDP8.hwe.vcf.gz
 bcftools view --threads 1 -Oz -i 'AN >= 90' ${inputfile} > ${inputfile%.vcf.gz}.AN25percMiss.vcf.gz
 bcftools index -t ${inputfile%.vcf.gz}.AN25percMiss.vcf.gz
 bcftools index -n ${inputfile%.vcf.gz}.AN25percMiss.vcf.gz
 vcftools --gzvcf ${inputfile%.vcf.gz}.AN25percMiss.vcf.gz --missing-indv --out ${inputfile%.vcf.gz}.AN25percMiss.MISSING_DATA

  cp females.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.SNPQ30.biall.scaffolds.fixedHeader.minDP8.hwe.AN25percMiss.vcf.gz ../SNP_A.FEMALES.minDP8.AN25percMiss.vcf.gz
