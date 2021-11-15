# Selection in Latin Americans
This repository contains the scripts used to detect and classify signals of selection in Latin Americans published in Mendoza-Revilla Javier et al. 2021

## Modules needed
* R/4.0.2
* plink/1.90b6.16
* tabix/0.2.6
* samtools/1.10
* vcftools/0.1.16
* plink2/2.00a2
* perl/5.30.1

## Preparing files
Our new software, `AdaptMix`, uses genotype data in ChromoPainter (CP) format. For a description of the file format see: https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromopainter_info.html.

For this tutorial we are going to use genomic data from Peruvians (`PEL`) from the 1000 Genomes Project (1KGP) as our target admixed population, and `CHB`, `IBS`, and `YRI` as our reference populations. Note that we are using `CHB` as a proxy for the Native American reference population.

Process the 1KGP VCF (merged file with ALL chromsomes) and extract 20K random SNPs:

```
bcftools query -f '%CHROM\t%POS\n' ALL.chrALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes2.vcf.gz \
| shuf -n 20000 | sort -n > 20KSNPs_toextract.txt

bcftools view ALL.chrALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes2.vcf.gz \
-S PEL_REFs_tokeep.txt -Ou | bcftools view -m2 -M2 -v snps -Ou | \
bcftools view -T 20KSNPs_toextract.txt -Ou | \
bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' -Ou | \
bcftools norm -d snps -Ov > PEL_REFs_ALLCHR_20K.vcf
```

Convert the VCF file to CP format to run AdaptMix:

```
plink2 --vcf PEL_REFs_small.chr22.vcf --export haps --out PEL_REFs_small.chr22
perl impute2chromopainter2.pl PEL_REFs_small.chr22.haps genetic_map_chr22_combined_b37.20140701.txt PEL_REFs_small.chr22.chromopainter
gzip PEL_REFs_small.chr22.chromopainter.haps
```

Convert the VCF file to PLINK format to run ADMIXTURE:

```
plink --vcf PEL_REFs_small.chr22.vcf --make-bed --out PEL_REFs_small.chr22
```

## Run ADMIXTURE 
Note that we are using ADMIXTURE to estimate ancestry proportions in PEL, but you can estimate this using other approaches e.g. SOURCEFIND (https://github.com/sahwa/sourcefindV2). 

```
plink --bfile PEL_REFs_small.chr22 --indep-pairwise 50 5 0.2 --out PEL_REFs_small.chr22
plink --bfile PEL_REFs_small.chr22 --extract PEL_REFs_small.chr22.prune.in --make-bed --out PEL_REFs_small_pruned.chr22
./admixture PEL_REFs_small_pruned.chr22.bed 3
```

## Run AdaptMix

```
Rscript compute_mr_sprime.R test.sprime.score Denisova_chr22.gtformat Vindija_chr22.gtformat > test.sprime.matchrates.txt
```

## Citation
Mendoza-Revilla Javier et al. "Disentangling signatures of selection before and after European colonization in Latin Americans." 
