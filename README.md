# Selection in Latin Americans
This repository contains the scripts used to detect and classify signals of selection in Latin Americans published in Mendoza-Revilla Javier et al. 2021

## Estimate ...

Get genotype information from archaic genomes i.e. Vindija Neanderthal and Denisovan Altai. Using `bcftools` do:

```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' Denisova_chr22.vcf.gz > Denisova_chr22.gtformat
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' Vindija_chr22.vcf.gz > Vindija_chr22.gtformat
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
