# Selection in Latin Americans
This repository contains the scripts used to detect and classify signals of selection in Latin Americans published in Mendoza-Revilla Javier et al. 2021

## Preparing files

For this tutorial we are going to use genomic data from Peruvians (`PEL`) from the 1000 Genomes Project (1KGP) as our target admixed population, and `CHB`, `IBS`, and `YRI` as our reference populations. Note that we are using `CHB` as a proxy for the Native American reference population.

Process the VCF for chr22 only:

```
bcftools view --force-samples ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
-S inds/AMR_YRI_tokeep.txt -Ou | bcftools view -m2 -M2 -v snps -Ou | \
bcftools norm -d snps -Ov > PEL_REFs_clean.chr22.vcf

```

Compute match rate introgression:

```
Rscript compute_mr_sprime.R test.sprime.score Denisova_chr22.gtformat Vindija_chr22.gtformat > test.sprime.matchrates.txt
```
## Run AdaptMix

```
Rscript compute_mr_sprime.R test.sprime.score Denisova_chr22.gtformat Vindija_chr22.gtformat > test.sprime.matchrates.txt
```

## Citation
Mendoza-Revilla Javier et al. "Disentangling signatures of selection before and after European colonization in Latin Americans." 
