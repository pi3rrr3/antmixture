import allel
import numpy as np

def vcf2npy(vcfpath):
  callset = allel.read_vcf(vcfpath)
  haplotypes_1 = callset['calldata/GT'][:,:,0]
  haplotypes_2 = callset['calldata/GT'][:,:,1]
  m, n = haplotypes_1.shape
  mat_haplo = np.empty((2*n, m))
  mat_haplo[::2] = haplotypes_1.T
  mat_haplo[1::2] = haplotypes_2.T
  return mat_haplo.astype(np.uint8)


import os

aq_nonFin = vcf2npy(os.path.join(os.pardir, '/scratch/project_2001099/nouhaudp/reseq/snp/phasing/shapeit', 'SNP_A.aquilonia_nonFin.females.phased.allScafs.vcf.gz'))
pol_nonFin = vcf2npy(os.path.join(os.pardir, '/scratch/project_2001099/nouhaudp/reseq/snp/phasing/shapeit', 'SNP_A.polyctena_nonFin.females.phased.allScafs.vcf.gz'))
pik = vcf2npy(os.path.join(os.pardir, '/scratch/project_2001099/nouhaudp/reseq/snp/phasing/shapeit', 'SNP_A.Pik.females.phased.allScafs.vcf.gz'))
bun = vcf2npy(os.path.join(os.pardir, '/scratch/project_2001099/nouhaudp/reseq/snp/phasing/shapeit', 'SNP_A.Bun.females.phased.allScafs.vcf.gz'))
lanw = vcf2npy(os.path.join(os.pardir, '/scratch/project_2001099/nouhaudp/reseq/snp/phasing/shapeit', 'SNP_A.LanW.females.phased.allScafs.vcf.gz'))
lanr = vcf2npy(os.path.join(os.pardir, '/scratch/project_2001099/nouhaudp/reseq/snp/phasing/shapeit', 'SNP_A.LanR.females.phased.allScafs.vcf.gz'))


import loter.locanc.local_ancestry as lc

# Non-Finnish reference panel

res_loter = lc.loter_smooth(l_H=[aq_nonFin, pol_nonFin], h_adm=pik, num_threads=1)
np.savetxt("loter_nonFin.pik.females.txt", res_loter, fmt="%i")

res_loter = lc.loter_smooth(l_H=[aq_nonFin, pol_nonFin], h_adm=bun, num_threads=1)
np.savetxt("loter_nonFin.bun.females.txt", res_loter, fmt="%i")

res_loter = lc.loter_smooth(l_H=[aq_nonFin, pol_nonFin], h_adm=lanw, num_threads=1)
np.savetxt("loter_nonFin.lanw.females.txt", res_loter, fmt="%i")

res_loter = lc.loter_smooth(l_H=[aq_nonFin, pol_nonFin], h_adm=lanr, num_threads=1)
np.savetxt("loter_nonFin.lanr.females.txt", res_loter, fmt="%i")
