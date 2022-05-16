# Selection in Latin Americans
This repository contains the scripts used to detect and classify signals of selection in Latin Americans published in Mendoza-Revilla Javier et al. 2022 

## Modules used
* R/4.0.2
* plink/1.90b6.16
* vcftools/0.1.16
* plink2/2.00a2
* perl/5.30.1

## Running AdaptMix

`AdaptMix` takes four files as input, and is run from a command line in the following way, providing four file names: 

R < run_AdaptMix.R parameter.input.file genotypes.input.filenames id.file output.file --no-save > screenoutput.out

The first three file names are input files, described below. The last is the output file name. "screenoutput.out" saves any technical output from AdaptMix (with no results) and can be discarded.

###################################
## Input File1:  parameter.input.file

The first input file contains the following 4 rows:
pop.vec: [pop1 pop2 ... pop_K]
surrogate.vec: [pop_1 pop_2 ... pop_S]
drift.maf.bins: [0.0 ... 0.5]
min.allele.freq.shift: [0.0,...,1.0]

The first row (pop.vec) lists the target populations you wish to test for admixture in. 

The second row (surrogate.vec) lists the surrogate populations that are to be tested for selection against the corresponding ancestral source populations that the target groups are assumed to descend from. These should be listed here in order of the ancestry proportion columns of "id.file" (see below) that provide the ancestral source populations' contributions to each target individual. For each pop.vec population, AdaptMix will test for (presumed pre-admixture) selection between each source, i.e. as defined by the ancestry proportion columns, and its corresponding surrogate.vec whose individuals are defined by the second column of the "id.file" (see below), in addition to testing for selection post-admixture in each pop.vec.

The third row specifies the bins of expected minor allele frequency (MAF) to use when binning SNPs to calculate drift under a neutral model. I.e. for each target population, drift will be calculated separately for SNPs within each bin. The first value must be 0, the last value 0.5, and values must be sequential inbetween these. The program will exit with an error message if any target population has no SNPs within one of the bins. For the paper "Disentangling Signatures of Selection Before and After European Colonization in Latin Americans", we set this as "0 0.5", which infers only a single drift value for all SNPs, regardless of their expected MAF. However, we have since updated the program, and we recommend an additional AdaptMix run that makes this as fine as possible (e.g. bins of size .01), so long as you have enough data, to eliminate (potentially false) small p-values.

The fourth row specifies the minimum drift value, translating to the expected (minimum) difference between observed and expected allele frequencies under neutrality. For example, a value of 0.01 indicates that it is normal, assuming neutrality, to have a frequency difference of 0.01 between the observed and expected frequencies. This also helps eliminate false positives for SNPs with low minor allele frequency. Such low MAF SNPs may have low inferred drift values, so that even very small differences between observed and expected frequencies are interpreted as evidence of selection unless a mininum value is set. (Use a very small non-zero value, e.g. 0.00000001, to essentially ignore this.)

#######################################
## Input File2: genotypes.input.filenames

This file contains X rows listing the X files (and their directory locations) that contain the genotype data for each target and surrogate individual. For example, this might contain 22 rows corresponding to the 22 chromosomes. 

Each of the files pointed to by "genotypes.input.filenames" should be in ChromoPainter (CP) file format as input, where individuals' haploid genomes are in rows and the columns are SNPs. The first row gives the number of haplotypes in the file, the second row gives the number of SNPs, and the third row a "P" in column 1 with remaining columns the basepair positions of each SNP. The remaining rows are individuals' haploid genomes (with no spaces), with each two consecutive rows an individual. The allowed value for each person's haploid at each SNP is {0,1,?}, where "?" denotes missing data. Though each individual is represented by two rows (i.e. two haplotypes) in CP format, haplotype information is ignored. I.e. you can randomly assign the two alleles for heterozygotes to either row of an individual.

We provide a script to prepare your data starting from a VCF. The example VCF file contains data from Peruvians (`PEL`) from the 1000 Genomes Project as the target admixed population, and `CHB`, `IBS`, and `YRI` as the reference populations. Note that `CHB` is used as a proxy for the Native American reference population.

We fist split the VCF files by chromosomes and convert to CP format using vcf_to_chromopainter_main.R (from https://github.com/sahwa/vcf_to_chromopainter):

```
for chr in {1..22}
do
  vcftools --vcf PEL_REFs_ALLCHR_20K.vcf --chr ${chr} --recode --stdout | gzip -c > PEL_REFs_chr${chr}_20K.vcf.gz 
  vcf_to_chromopainter_main.R -g PEL_REFs_chr${chr}_20K.vcf.gz -o PEL_REFs_chr${chr}_20K.chromopainter.haps -u F
done
```

Alternatively, if the data is in haps/sample (SHAPEIT) format, we can convert to CP format using `impute2chromopainter2.pl`. We show this by first converting the VCF file to haps/sample file using `plink2`.

```
for chr in {1..22}
do
  ./plink2 --vcf PEL_REFs_ALLCHR_20K.vcf --chr ${chr} --export haps --out PEL_REFs_ALLCHR_20K_chr${chr}

  perl impute2chromopainter2.pl PEL_REFs_ALLCHR_20K_chr${chr}.haps genetic_map_chr${chr}_combined_b37.20140701.txt PEL_REFs_ALLCHR_20K_chr${chr}.chromopainter
  gzip PEL_REFs_ALLCHR_20K_chr${chr}.chromopainter.haps

done
```

################################
## Input File 3: id.file

Finally, AdaptMix also requires an id file. This file contains 2+S columns, where S is the number of ancestral populations specified in "surrogate.vec" of the "parameter.input.file" input file above. (E.g. S may be the number of clusters K when running ADMIXTURE.) 

Each row of "id.file" is a person, with column 1 containing their individual ID and column 2 containing their population ID label. These ID labels must match those provided in "pop.vec" and "surrogate.vec" of "parameter.input.file" (though not all labels in column 2 need be among "pop.vec" or "surrogate.vec"). The remaining S columns give the inferred ancestry proportions for each individual (e.g. as inferred by ADMIXTURE), with each column corresponding to the populations -- as ordered -- in surrogate.vec.

#################################
###########################
## RUNNING ADMIXTURE TO GET ANCESTRY PROPORTIONS 

We can convert the VCF file to PLINK format to run ADMIXTURE (output will be needed to run `AdaptMix`):

```
plink --vcf PEL_REFs_ALLCHR_20K.vcf --make-bed --out PEL_REFs_ALLCHR_20K
```

## Run ADMIXTURE 
Note that we are using ADMIXTURE to estimate ancestry proportions in PEL, but you can estimate these proportions using other approaches e.g. SOURCEFIND (https://github.com/sahwa/sourcefindV2).

```
./admixture PEL_REFs_ALLCHR_20K.bed 3
```

#########################################
########################################3
## RUN AdaptMix EXAMPLE:

```
R < run_AdaptMix.R example/PEL_analysis_paramfile.txt example/PEL_REFs_ALLCHR_chr.txt example/PEL_REFs.ids.txt example/PEL_REFs_ALLCHR_adaptmix.txt --no-save > screenoutput.out
```

Input arguments:
* 1 - Parameter file specifying the target and ancestral populations
* 2 - File listing each of the CP genotype files to be used for analysis
* 3 - ID file
* 4 - Output file name

########################################
#######################################
## AdaptMix OUTPUT FILE:

The first lines give the drift estimates for each target admixed population (row) for SNPs within each bin (column) as specified by "drift.maf.bin".

Following this, the next line contains the header of the output file for the remaining lines, of which there is one line per SNP:

* 1 - chrom (chromosome ID)
* 2 - pos (position in base pair)
* 3 - log10.pval.target.1 (-log10 P-value, testing for any evidence of selection)
* 4 - obs.freq.target.1 (observed allele frequency in the target admixed population)
* 5 - exp.freq.target.1 (expected allele frequency in the target admixed population)
* 6 - AIC.neutral.target.1 (AIC for SNP being neutral, = -2*log-likelihood(neutral))
* 7 - AIC.postadmix.target.1 (AIC for SNP being selected post-admixture, = 2-2*log-likelihood(post-admixture))
* 8 - AIC.insurr.source1target.1 (AIC for SNP being selected in surrogate population 1, = 2-2*log-likelihood(sel.in.pop1))
* 9 - AIC.insurr.source2target.1 (AIC for SNP being selected in surrogate population 2, = 2-2*log-likelihood(sel.in.pop2))
* 10 - AIC.insurr.source3target.1 (AIC for SNP being selected in surrogate population 3, = 2-2*log-likelihood(sel.in.pop3))
* 11 - sel.postadmix.target.1 (selection coefficient for SNP being selected post-admixture)
* 12 - sel.insurr.source1target.1 (selection coefficient for SNP being selected in surrogate population 1)
* 13 - sel.insurr.source2target.1 (selection coefficient for SNP being selected in surrogate population 2)
* 14 - sel.insurr.source3target.1 (selection coefficient for SNP being selected in surrogate population 3)

Note that depending on the number of target/surrogate population used the number of columns will be different from that shown here. The above assumes 3 surrogate ("surrogate.vec") populations and one target ("pop.vec") population. 

## Citation
Mendoza-Revilla, Javier, et al. "Disentangling signatures of selection before and after European colonization in Latin Americans." Molecular biology and evolution 39.4 (2022): msac076.
