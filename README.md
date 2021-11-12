# Selection in Latin Americans
This repository contains the scripts used to detect and classify signals of selection in Latin Americans published in Mendoza-Revilla Javier et al. 2021

## Modules needed
* module load R/4.0.2
* module load plink/1.90b6.16
* module load tabix/0.2.6
* module load samtools/1.10
* module load vcftools/0.1.16

## Preparing files

For this tutorial we are going to use genomic data from Peruvians (`PEL`) from the 1000 Genomes Project (1KGP) as our target admixed population, and `CHB`, `IBS`, and `YRI` as our reference populations. Note that we are using `CHB` as a proxy for the Native American reference population.

Process the 1KP VCF for chr22 only:

```
bcftools view ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
-S inds/PEL_REFs_tokeep.txt -Ou | bcftools view -m2 -M2 -v snps -Ou | \
bcftools norm -d snps -Ov > PEL_REFs.chr22.vcf
```

Convert the VCF file to Chromopainter (CP) format:

```
VCF2CP.pl PEL_REFs_clean.chr22.vcf
```

```
Rscript compute_mr_sprime.R test.sprime.score Denisova_chr22.gtformat Vindija_chr22.gtformat > test.sprime.matchrates.txt
```
## Run AdaptMix

```
Rscript compute_mr_sprime.R test.sprime.score Denisova_chr22.gtformat Vindija_chr22.gtformat > test.sprime.matchrates.txt
```

## Citation
Mendoza-Revilla Javier et al. "Disentangling signatures of selection before and after European colonization in Latin Americans." 
