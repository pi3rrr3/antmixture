# gnmcs environment
conda create -p /scratch/project_2001099/nouhaudp/conda/gnmcs
source activate /scratch/project_2001099/nouhaudp/conda/gnmcs
conda install -c bioconda eigensoft tabix shapeit4 mummer4 minimap2 phyml vt bamutil -p /scratch/project_2001099/nouhaudp/conda/gnmcs
conda install mosdepth vcflib -p /scratch/project_2001099/nouhaudp/conda/gnmcs
conda install -c etetoolkit ete3 -p /scratch/project_2001099/nouhaudp/conda/gnmcs
conda install whatshap #nomkl
conda install -c numpy pandas scikit-learn scikit-allel scipy -p /scratch/project_2001099/nouhaudp/conda/gnmcs
conda install sniffles -p /scratch/project_2001099/nouhaudp/conda/gnmcs
conda install ngmlr -p /scratch/project_2001099/nouhaudp/conda/gnmcs
conda install samtools=1.9 -p /scratch/project_2001099/nouhaudp/conda/gnmcs
conda install -c conda-forge msprime=1.0.2 -p /scratch/project_2001099/nouhaudp/conda/gnmcs
conda install -c conda-forge pixy -p /scratch/project_2001099/nouhaudp/conda/gnmcs
