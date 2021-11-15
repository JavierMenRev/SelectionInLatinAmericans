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
`AdaptMix` uses the ChromoPainter (CP) file format as input. We provide a script to prepare your data starting from a VCF file. 

The example VCF file contains data from Peruvians (`PEL`) from the 1000 Genomes Project as the target admixed population, and `CHB`, `IBS`, and `YRI` as the reference populations. Note that `CHB` is used as a proxy for the Native American reference population.

We fist split the VCF files by chromosomes, convert the VCF to haps/sample format, and then to CP format:

```
for chr in {1..22}
do
  ./plink2 --vcf PEL_REFs_ALLCHR_20K.vcf --chr ${chr} --export haps --out PEL_REFs_ALLCHR_20K_chr${chr}

  perl impute2chromopainter2.pl PEL_REFs_ALLCHR_20K_chr${chr}.haps genetic_map_chr${chr}_combined_b37.20140701.txt PEL_REFs_ALLCHR_20K_chr${chr}.chromopainter
  gzip PEL_REFs_ALLCHR_20K_chr${chr}.chromopainter.haps

done
```

We will also need an id file to run `AdaptMix`. This file has three columns:
* 1 - Population ID
* 2 - Individual ID
* 3 - Flag 1 or 0 whether to include the individual in the analysis or not

```
Rscript make_id_file.R PEL_REFs_ALLCHR_20K_chr1.sample PEL_REFs_ids.txt > PEL_REFs_ALLCHR_20K.ids.txt
```

We will lastly convert the VCF file to PLINK format to run ADMIXTURE (output will be needed to run `AdaptMix`):

```
plink --vcf PEL_REFs_ALLCHR_20K.vcf --make-bed --out PEL_REFs_ALLCHR_20K
```

## Run ADMIXTURE 
Note that we are using ADMIXTURE to estimate ancestry proportions in PEL, but you can estimate these proportions using other approaches e.g. SOURCEFIND (https://github.com/sahwa/sourcefindV2).

```
./admixture PEL_REFs_ALLCHR_20K.bed 3
```

## Run AdaptMix

```
Rscript run_AdaptMix.R PEL PEL_REFs_ALLCHR_20K_chr .chromopainter.haps.gz PEL_REFs_ALLCHR_20K.ids.txt PEL_REFs_ALLCHR_20K.3.Q PEL_REFs_ALLCHR_20K.fam CHB,IBS,YRI PEL_REFs_ALLCHR_20K_adaptmix.txt 
```

Input arguments:
* 1 - Population ID of the targed admixed population
* 2 - Prefix of CP file (before "chr" string in file)
* 3 - Postfix of CP file (after "chr" string in file)
* 4 - ID file
* 5 - Q file from ADMIXTURE (or any other software) 
* 6 - Fam file from PLINK used to run ADMIXTURE (ids should be the same order as the Q file)
* 7 - Reference population IDs separated by a comma (the order of the reference population should match the ancestries i.e. columns in the ADMIXTURE file)
* 8 - Output file

Output:
The first two lines indicate the population ID of the target admixed population and the drift estimate.
Line 3 contains the header of the output file:

* 1 - chrom (chromosome ID)
* 2 - pos (position in base pair)
* 3 - log10.pval.target.1 (-log10 P-value)
* 4 - obs.freq.target.1 (observed allele frequency in the targed admixed population)
* 5 - exp.freq.target.1 (expected allele frequency in the targed admixed population)
* 6 - AIC.neutral.target.1 (AIC for SNP being neutral)
* 7 - AIC.postadmix.target.1 (AIC for SNP being selected post-admixture)
* 8 - AIC.insurr.source1target.1 (AIC for SNP being selected in reference population 1)
* 9 - AIC.insurr.source2target.1 (AIC for SNP being selected in reference population 2)
* 10 - AIC.insurr.source3target.1 (AIC for SNP being selected in reference population 3)
* 11 - sel.postadmix.target.1 (selection coefficient in selection strenght times number of generation units assuming SNP is selected post-admixture)
* 12 - sel.insurr.source1target.1 (selection coefficient in selection strenght times number of generation units assuming SNP is selected in reference population 1)
* 13 - sel.insurr.source2target.1 (selection coefficient in selection strenght times number of generation units assuming SNP is selected in reference population 2)
* 14 - sel.insurr.source3target.1 (selection coefficient in selection strenght times number of generation units assuming SNP is selected in reference population 3)
* 15 - LRT.pval.postadmix.target.1 (selection coefficient in selection strenght times number of generation units assuming SNP is selected post-admixture)
* 16 - LRT.pval.insurr.source1.target.1 
* 17 - LRT.pval.insurr.source2.target.1 
* 18 LRT.pval.insurr.source3.target.1










## Citation
Mendoza-Revilla Javier et al. "Disentangling signatures of selection before and after European colonization in Latin Americans." 
