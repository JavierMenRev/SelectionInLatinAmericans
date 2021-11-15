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

The example VCF file contains data from Peruvians (`PEL`) from the 1000 Genomes Project (1KGP) as our target admixed population, and `CHB`, `IBS`, and `YRI` as our reference populations. Note that we are using `CHB` as a proxy for the Native American reference population.

We fist split the VCF files by chromosomes, convert the VCF to haps/sample format, and lastly to CP format:

```
for chr in {1..22}
do
  ./plink2 --vcf PEL_REFs_ALLCHR_20K.vcf --chr ${chr} --export haps --out PEL_REFs_ALLCHR_20K_chr${chr}

  perl impute2chromopainter2.pl \
  PEL_REFs_ALLCHR_20K_chr${chr}.haps genetic_map_chr${chr}_combined_b37.20140701.txt PEL_REFs_ALLCHR_20K_chr${chr}.chromopainter
  gzip PEL_REFs_ALLCHR_20K_chr${chr}.chromopainter.haps

done
```

We will also need an id file to run `AdaptMix`. This file simply contains the individual id (column 1), population id (column 2), and a 1 or 0 flag that indicates whether to include the individual in the analysis or not (column 3). We can process the sample file which has the individuals in the same order as our CP files with a small `R` script:

```
Rscript make_id_file.R PEL_REFs_ALLCHR_20K_chr1.sample PEL_REFs_ids.txt > PEL_REFs_ALLCHR_20K.ids.txt
```

The PEL_REFs_ids.txt contains the individual id in column 1 and the population id in column 2. The script also assumes all individuals will be used in the analysis. We then convert the VCF file to PLINK format to run ADMIXTURE (output will be needed to run AdaptMix later):

```
plink --vcf PEL_REFs_ALLCHR_20K.vcf --make-bed --out PEL_REFs_ALLCHR_20K
```

## Run ADMIXTURE 
Note that we are using ADMIXTURE to estimate ancestry proportions in PEL, but you can estimate this using other approaches e.g. SOURCEFIND (https://github.com/sahwa/sourcefindV2).

```
./admixture PEL_REFs_ALLCHR_20K.bed 3
```

## Run AdaptMix

```
Rscript run_AdaptMix.R PEL PEL_REFs_ALLCHR_20K_chr .chromopainter.haps.gz PEL_REFs_ALLCHR_20K.ids.txt PEL_REFs_ALLCHR_20K.3.Q PEL_REFs_ALLCHR_20K.fam CHB,IBS,YRI PEL_REFs_ALLCHR_20K_adaptmix.txt 
```

Input arguments:

1 - Population ID of the targed admixed population
2 - Prefix of CP file
3 - Postfix of CP file
4 - ID file
5 - Q file from ADMIXTURE (or any other software) 
6 - Fam file from PLINK used to run ADMIXTURE (ids should be the same order as the Q file)
7 - Reference population IDs separated by a comma (the order of the reference population should match the ancestries i.e. columns in the ADMIXTURE file)
8 - Output file


## Citation
Mendoza-Revilla Javier et al. "Disentangling signatures of selection before and after European colonization in Latin Americans." 
