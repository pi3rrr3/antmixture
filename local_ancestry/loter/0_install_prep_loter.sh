# local_ancestry_Loter.sh


###
#### Install Loter --------------------------------------------------------------------------------
###

cd /projappl/project_2001099
git clone https://github.com/bcm-uga/Loter.git
module load bioconda/3 ; source activate /scratch/project_2001099/nouhaudp/conda/gnmcs
cd Loter ; python setup.py install



###
#### Create one VCF per parental pool or hybrid population ----------------------------------------
###

cd /scratch/project_2001099/nouhaudp/reseq/snp/phasing/shapeit

VCF=SNP_A.FEMALES.minDP8.AN15percMiss.whap.shapeit.allScafs.vcf.gz

bcftools view -S aq.list $VCF -Oz > SNP_A.aquilonia.females.phased.allScafs.vcf.gz
bcftools view -S pol.list $VCF -Oz > SNP_A.polyctena.females.phased.allScafs.vcf.gz
bcftools view -S aq_nonFin.list $VCF -Oz > SNP_A.aquilonia_nonFin.females.phased.allScafs.vcf.gz
bcftools view -S pol_nonFin.list $VCF -Oz > SNP_A.polyctena_nonFin.females.phased.allScafs.vcf.gz
bcftools view -S LanW_female.list $VCF -Oz > SNP_A.LanW.females.phased.allScafs.vcf.gz
bcftools view -S LanR_female.list $VCF -Oz > SNP_A.LanR.females.phased.allScafs.vcf.gz
bcftools view -S Bun_female.list $VCF -Oz > SNP_A.Bun.females.phased.allScafs.vcf.gz
bcftools view -S Pik_female.list $VCF -Oz > SNP_A.Pik.females.phased.allScafs.vcf.gz

bcftools view -H /scratch/project_2001099/nouhaudp/reseq/snp/phasing/shapeit/SNP_A.aquilonia.females.phased.allScafs.vcf.gz



###
#### Run ------------------------------------------------------------------------------------------
###

cd /scratch/project_2001099/nouhaudp/reseq/snp/analyses/local_ancestry/loter

# Extract positions
gunzip -c ../../../phasing/shapeit/SNP_A.FEMALES.minDP8.AN15percMiss.whap.shapeit.allScafs.vcf.gz | grep -v '^#' | cut -f1,2 > positions.tab

# Following https://github.com/bcm-uga/Loter/blob/eb8970d5b6abddf9e6f4665212769f36e62d1d73/python-package/Local_Ancestry_Example.ipynb
module load bioconda/3 ; source activate /scratch/project_2001099/nouhaudp/conda/gnmcs
python run_loter.py


#
##
### This is the end.
