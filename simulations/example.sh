#!/bin/bash -l
#SBATCH -J bunpik
#SBATCH -o /scratch/project_2001099/nouhaudp/reseq/snp/analyses/twisst/sims/logs/LanWLanR_SO/%a.out
#SBATCH --account=project_2001099
#SBATCH -t 24:00:00
#SBATCH -p small
#SBATCH --array=1-100
#SBATCH --ntasks 1
#SBATCH --mem=2G

module load bioconda/3
source activate $PROJCONDA/gnmcs

mkdir /scratch/project_2001099/nouhaudp/reseq/snp/analyses/twisst/sims/LanWLanR_SO/sim_$SLURM_ARRAY_TASK_ID
cd  /scratch/project_2001099/nouhaudp/reseq/snp/analyses/twisst/sims/LanWLanR_SO/sim_$SLURM_ARRAY_TASK_ID

cp -r /projappl/project_2001099/twisst .
cp -r /projappl/project_2001099/genomics_general .
cp -r /scratch/project_2001099/nouhaudp/reseq/snp/analyses/twisst/sims/LanWLanR_SO/LanWLanR_SO.py twisst/.

cd twisst
export PYTHONPATH="$PROJAPPL/twisst:$PYTHONPATH"
export PYTHONPATH="$PROJAPPL/genomics_general:$PYTHONPATH"

/scratch/project_2001099/nouhaudp/conda/gnmcs/bin/python LanWLanR_SO.py
