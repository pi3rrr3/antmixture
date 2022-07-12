# twisst_run.sh
 

###
### Write geno files ------------------------------------------------------------------------------
###


# Create ploidy table
cut -f1,6 sample_table_R_ordered_LanFix.tab | grep -v 'id' | awk '$2 == "f" {print $0}' | sed 's/     f/      2/' > ploidy.tab

# sample list
gunzip -c SNP_A.FEMALES.minDP8.AN15percMiss.whap.shapeit.allScafs.vcf.gz | grep 'CHROM' | sed 's/        /,/g' | perl -npe 's/\#CHROM.+FORMAT,//'

# Phased dataset
{
python $PROJAPPL/genomics_general/VCF_processing/parseVCF.py \
  -i SNP_A.FEMALES.minDP8.AN15percMiss.whap.shapeit.allScafs.vcf.gz \
  --skipIndels --excludeDuplicates \
  --ploidyFile ploidy.tab \
  -s CBAQ1_1w,CBAQ2_2w,CBAQ3_1w,Lai_1w,Lai_2w,Loa_1w,CF14a_1w,CF4b_1w,CF8b_1w,Pus2_1w,Bun24_10w,Bun24_11w,Bun24_12w,Bun24_8w,Bun24_9w,Bun26_4w,Bun26_5w,Bun26_6w,Bun26_7w,Bun26_8w,FA07_4q,FA15_4q,FA15_5q,FAu14_2q,FAu14_5q,FA12_10q,FA12_11q,FA12_12q,FA12_8q,FA12_9q,FA07_5q,FA07_6q,FA16_4q,FA16_5q,FA16_6q,FA04_4q,FA04_5q,FA04_6q,FA04_7q,FA04_8q,Mar1_15q,Mar1_16q,Mar1_6q,Mar1_7q,Mar1_8q,Mar2_10q,Mar2_6q,Mar2_7q,Mar2_8q,Mar2_9q,Att1_1w,Fis2_1w,Jar6_1w,Lok3_1w,CBCH1_1w,CBCH2_2w,CBCH3_1w,CAGa_1w,NAZa_1w,VDa_1w \
  | bgzip > phased.geno.gz &
}





###
### Get trees for variable window sizes -----------------------------------------------------------
###

{
#!/bin/bash -l
#SBATCH -J twisst
#SBATCH --account=project_2001099
#SBATCH -t 1-00:00:00
#SBATCH -p small
#SBATCH --ntasks=4
#SBATCH --mem=8G

cd /scratch/project_2001099/nouhaudp/reseq/snp/analyses/twisst

DATADIR=/scratch/project_2001099/nouhaudp/reseq/snp/geno/
RESDIR=/scratch/project_2001099/nouhaudp/reseq/snp/analyses/twisst/trees

module load bioconda/3
source activate /scratch/project_2001099/nouhaudp/conda/gnmcs
export PYTHONPATH="$PROJAPPL/genomics_general:$PYTHONPATH"

# Fraction of missing data allowed
fractMiss=10

for x in 200 100 50 20 ; do

# compute max missing data per individual allowed (set to 15%)
(( Mi = x - x*fractMiss/100 ))

echo "Inferring trees with window size $x SNPs, at least $Mi sites per individual"

python3 $PROJAPPL/genomics_general/phylo/phyml_sliding_windows.py -T 4 \
  -g $DATADIR/phased_out.geno.gz \
  --prefix $RESDIR/phased_out_max${fractMiss}percMi.phyml_bionj.w$x \
  -w $x -Mi $Mi --windType sites --model GTR --optimise n

done

}





###
### Run TWISST ------------------------------------------------------------------------------------
###

### Prepare group (hierarchy) file
awk '$6 == "f"' sample_table_R_ordered_LanFix.tab | cut -f1,15 > tmp
perl -npe 's/      /_A     /' tmp > tmpA
perl -npe 's/      /_B     /' tmp > tmpB
cat tmpA tmpB | sort > group.tab
rm tmp tmpA tmpB


### Run Twisst

{

cd /scratch/project_2001099/nouhaudp/reseq/snp/analyses/twisst

for x in 20 50 100 200 ; do

echo "Running Twisst for window size $x"

  for pop in Bun LanR LanW Pik ; do

    echo "Analysing $pop"

    python3 $PROJAPPL/twisst/twisst.py \
    -t "trees/phased_out_max10percMi.phyml_bionj.w"$x".trees.gz" \
    -w "results/"$pop".w"$x".weights.csv.gz" \
    -g Aq \
    -g Pol \
    -g $pop \
    -g Out \
    --outgroup Out \
    --groupsFile group.tab

  done

done

}



#
##
### This is the end.
