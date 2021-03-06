# antmixture


This repo contains scripts related to the preprint entitled "Rapid and repeatable genome evolution across three hybrid ant populations" available on [bioRxiv](https://doi.org/10.1101/2022.01.16.476493).

## Bioinformatic pipeline 
The folder `reads2SNPs` contains scripts for trimming, mapping and deduplicating reads, and then calling, filtering and phasing variant sites (check the preprint for software versions). Note the script 7_SNP_calling_4pixy.sh performs variant and invariant site calling to estimate diversity & divergence metrics in 20kb windows with [pixy](https://doi.org/10.1111/1755-0998.13326).

## Data analysis
The folder `local_ancestry` contains scripts to perform local ancestry inference using [Loter](https://doi.org/10.1093/molbev/msy126), topology weighting with [Twisst](https://doi.org/10.1534/genetics.116.194720) and good old chromosome painting.
