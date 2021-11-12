# Selection in Latin Americans
This repository contains the scripts used to detect and classify signals of selection in Latin Americans published in Mendoza-Revilla Javier et al. 2021

## Modules needed
* R/4.0.2
* plink/1.90b6.16
* tabix/0.2.6
* samtools/1.10
* vcftools/0.1.16
* plink2/2.00a2

## Preparing files

For this tutorial we are going to use genomic data from Peruvians (`PEL`) from the 1000 Genomes Project (1KGP) as our target admixed population, and `CHB`, `IBS`, and `YRI` as our reference populations. Note that we are using `CHB` as a proxy for the Native American reference population.

Process the 1KGP VCF for chr22 only:

```
bcftools view ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz \
-S PEL_REFs_tokeep.txt -Ou | bcftools view -m2 -M2 -v snps -Ou | \
bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' -Ou | \
bcftools norm -d snps -Ov > PEL_REFs.chr22.vcf
```

Grab only 100K SNPs to make this tutorial run quicker:

```
bcftools query -f '%CHROM\t%POS\n' PEL_REFs.chr22.vcf | shuf -n 100000 | sort -n > PEL_REFs_100KSNPs.chr22.txt
bcftools view PEL_REFs.chr22.vcf -T PEL_REFs_100KSNPs.chr22.txt -Ov > PEL_REFs_small.chr22.vcf
```

Convert the VCF file to Chromopainter (CP) format to run AdaptMix:

```
plink2 --vcf PEL_REFs_small.chr22.vcf --export hapslegend --out PEL_REFs_small.chr22
```

Convert the VCF file to PLINK format to run ADMIXTURE:

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
