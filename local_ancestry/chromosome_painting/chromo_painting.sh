#!/bin/bash -l

# this is a wrapper to ID ancestry-informative markers and genotype them in hybrids (painting).
# Each job will be run on one scaffold independently.

#SBATCH -J chromopaint
#SBATCH -o logs/chromopaint_80perc_%a.out
#SBATCH -e logs/chromopaint_80perc_%a.err
#SBATCH --account=project_2001099
#SBATCH -t 72:00:00
#SBATCH -p small
#SBATCH --array=1-27
#SBATCH --ntasks 1
#SBATCH --mem=8G

cd /scratch/project_2001099/nouhaudp/reseq/snp/analyses/local_ancestry/chromopaint

# Load r-env-singularity
module load r-env-singularity

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2001099" >> ~/.Renviron

# Get scaffold name
scaffold=$(sed -n "$SLURM_ARRAY_TASK_ID"p scaffold.list)
mkdir $scaffold
cp chromo_painting.R $scaffold
cd $scaffold

# Run the R script
srun singularity_wrapper exec Rscript --no-save chromo_painting.R $scaffold
