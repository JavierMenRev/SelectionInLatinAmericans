# Selection in Latin Americans
This repository contains the scripts used to detect and classify signals of selection in Latin Americans published in Mendoza-Revilla Javier et al. 2022 

This project was directed by Garrett Hellenthal (University College London) and Andres Ruiz-Linares (Université d'Aix-Marseille). Code developed by Garrett Hellenthal and Javier Mendoza Revilla.

## Modules used
* R/4.0.2
* tabix/0.2.6
* vcftools/0.1.16
* plink/1.90b6.16
* plink2/2.00a2
* perl/5.30.1

## Running AdaptMix

`AdaptMix` takes four files as input, and is run from a command line in the following way, providing four file names: 

```
R < run_AdaptMix.R parameter.input.file genotypes.input.filenames id.file output.file --no-save > screenoutput.out
```

The first three file names are input files, described below. The last is the output file name. "screenoutput.out" saves any technical output from AdaptMix (with no results) and can be discarded.

### Input File 1:  parameter.input.file

The first input file contains the following 4 rows:
* pop.vec: [pop1 pop2 ... pop_K]
* surrogate.vec: [pop_1 pop_2 ... pop_S]
* drift.maf.bins: [0.0 ... 0.5]
* min.allele.freq.shift: [0.0,...,1.0]

The first row (pop.vec) lists the target populations you wish to test for admixture in. 

The second row (surrogate.vec) lists the surrogate populations that are to be tested for selection against the corresponding ancestral source populations that the target groups are assumed to descend from. These should be listed here in order of the ancestry proportion columns of "id.file" (see below) that provide the ancestral source populations' contributions to each target individual. For each pop.vec population, AdaptMix will test for (presumed pre-admixture) selection between each source, i.e. as defined by the ancestry proportion columns, and its corresponding surrogate.vec whose individuals are defined by the second column of the "id.file" (see below), in addition to testing for selection post-admixture in each pop.vec.

The third row specifies the bins of expected minor allele frequency (MAF) to use when binning SNPs to calculate drift under a neutral model. I.e. for each target population, drift will be calculated separately for SNPs within each bin. The first value must be 0, the last value 0.5, and values must be sequential in between these. The program will exit with an error message if any target population has no SNPs within one of the bins. For the paper "Disentangling Signatures of Selection Before and After European Colonization in Latin Americans", we set this as "0 0.5", which infers only a single drift value for all SNPs, regardless of their expected MAF. However, we have since updated the program, and we recommend an additional AdaptMix run that makes this as fine as possible (e.g. bins of size 0.01), so long as you have enough data, to eliminate (potentially false) small P-values.

The fourth row specifies the minimum drift value, translating to the expected (minimum) difference between observed and expected allele frequencies under neutrality. For example, a value of 0.01 indicates that it is normal, assuming neutrality, to have a frequency difference of 0.01 between the observed and expected frequencies. This also helps eliminate false positives for SNPs with low minor allele frequency. Such low MAF SNPs may have low inferred drift values, so that even very small differences between observed and expected frequencies are interpreted as evidence of selection unless a mininum value is set. (Use a very small non-zero value, e.g. 0.00000001, to essentially ignore this.)

### Input File 2: genotypes.input.filenames

This file contains n rows listing the n files (and their directory locations) that contain the genotype data for each target and surrogate individual. For example, this might contain 22 rows corresponding to the 22 chromosomes. 

Each of the files pointed to by "genotypes.input.filenames" should be in ChromoPainter (CP) file format as input, where individuals' haploid genomes are in rows and the columns are SNPs. The first row gives the number of haplotypes in the file, the second row gives the number of SNPs, and the third row a "P" in column 1 with remaining columns the basepair positions of each SNP. The remaining rows are individuals' haploid genomes (with no spaces), with each two consecutive rows an individual. The allowed value for each person's haploid at each SNP is {0,1,?}, where "?" denotes missing data. Though each individual is represented by two rows (i.e. two haplotypes) in CP format, haplotype information is ignored. I.e. you can randomly assign the two alleles for heterozygotes to either row of an individual.

### Input File 3: id.file

This file contains 2+S columns, where S is the number of ancestral populations specified in "surrogate.vec" of the "parameter.input.file" input file above. (E.g. S may be the number of clusters K when running ADMIXTURE.) 

Each row of "id.file" is a person, with column 1 containing their individual ID and column 2 containing their population ID label. These ID labels must match those provided in "pop.vec" and "surrogate.vec" of "parameter.input.file" (though not all labels in column 2 need be among "pop.vec" or "surrogate.vec"). The remaining S columns give the inferred ancestry proportions for each individual (e.g. as inferred by ADMIXTURE), with each column corresponding to the populations -- as ordered -- in surrogate.vec.

## Workflow

We provide a script to prepare your data starting from a VCF. The example VCF file contains data from Peruvians (PEL) from the 1000 Genomes Project as the target admixed population, and Han Chinese (CHB), Spanish (IBS), and Yoruba (YRI) as the reference populations. Note that CHB is used as a proxy for the Native American reference population.

We fist split the VCF files by chromosomes and convert to CP format using `vcf_to_chrompainter_AdaptMix.R`:

```
for chr in {1..22}
do
  vcftools --vcf ./example/PEL_REFs_ALLCHR_8KSNPs.vcf --chr ${chr} --recode --stdout > ./example/PEL_REFs_chr${chr}.vcf 
  Rscript vcf_to_chrompainter_AdaptMix.R ./example/PEL_REFs_chr${chr} ./example/PEL_REFs_chr${chr}
  gzip ./example/PEL_REFs_chr${chr}.chromopainter.haps
done
```

Alternatively, if the data is in haps/sample (SHAPEIT) format, we can convert to CP format using `impute2chromopainter2.pl`. We show this by first converting the VCF file to haps/sample file using `plink2`.

```
for chr in {1..22}
do
  plink2 --vcf ./example/PEL_REFs_ALLCHR_8KSNPs.vcf --chr ${chr} --export haps --out ./example/PEL_REFs_chr${chr}
  
  perl impute2chromopainter2.pl ./example/PEL_REFs_chr${chr}.haps ./example/genetic_map/genetic_map_chr${chr}_combined_b37.txt ./example/PEL_REFs_chr${chr}.chromopainter
  gzip ./example/PEL_REFs_chr${chr}.chromopainter.haps
done
```

The next step is to obtain ancestry proportions for the individuals in the target population. We will illustrate this using ADMIXTURE, but you can estimate these proportions using other approaches e.g. SOURCEFIND (See: https://github.com/sahwa/sourcefindV2).

We first convert the VCF files to PLINK format and then run ADMIXTURE using K=3:

```
plink --vcf ./example/PEL_REFs_ALLCHR_8KSNPs.vcf --make-bed --out ./example/PEL_REFs_ALLCHR_8KSNPs
admixture ./example/PEL_REFs_ALLCHR_8KSNPs.bed 3
```

We can then use ADMIXTURE's Q output file (ancestry proportions) to make the "id.file" (see above). We note again that the ancestry proportions for each individual should correspond to the reference populations -- as ordered -- in surrogate.vec.

Ordering the ancestry proportions as listed in surrogate.vec (see below) we get the following estimates for the first 5 individuals:

```
HG01565 PEL 0.699085 0.300905 0.00001
HG01566 PEL 0.398061 0.601929 0.00001
HG01571 PEL 0.713799 0.286191 0.00001
HG01572 PEL 0.944167 0.055823 0.00001
HG01577 PEL 0.301734 0.650348 0.047919
```

Run AdaptMix:

```
R < run_AdaptMix.R ./example/PEL_REFs_paramfile.txt ./example/PEL_REFs_genotypes_files.txt ./example/PEL_REFs_ids.txt ./example/PEL_REFs_ALLCHR_adaptmix.txt --no-save > screenoutput.out
```

Input arguments:
* 1 - Parameter file specifying the target and reference populations
* 2 - File listing each of the CP genotype files to be used for analysis
* 3 - ID file
* 4 - Output file name

AdaptMix output file:

The first lines give the drift estimates for each target admixed population (row) for SNPs within each bin (column) as specified by "drift.maf.bin".

Following this, the next line contains the header of the output file for the remaining lines, of which there is one line per SNP:

```
1 - chrom (chromosome ID)
2 - pos (position in base pair)
3 - log10.pval.target.1 (-log10 P-value, testing for any evidence of selection)
4 - obs.freq.target.1 (observed allele frequency in the target admixed population)
5 - exp.freq.target.1 (expected allele frequency in the target admixed population)
6 - AIC.neutral.target.1 (AIC for SNP being neutral, = -2*log-likelihood(neutral))
7 - AIC.postadmix.target.1 (AIC for SNP being selected post-admixture, = 2-2*log-likelihood(post-admixture))
8 - AIC.insurr.source1target.1 (AIC for SNP being selected in surrogate population 1, = 2-2*log-likelihood(sel.in.pop1))
9 - AIC.insurr.source2target.1 (AIC for SNP being selected in surrogate population 2, = 2-2*log-likelihood(sel.in.pop2))
10 - AIC.insurr.source3target.1 (AIC for SNP being selected in surrogate population 3, = 2-2*log-likelihood(sel.in.pop3))
11 - sel.postadmix.target.1 (selection coefficient for SNP being selected post-admixture)
12 - sel.insurr.source1target.1 (selection coefficient for SNP being selected in surrogate population 1)
13 - sel.insurr.source2target.1 (selection coefficient for SNP being selected in surrogate population 2)
14 - sel.insurr.source3target.1 (selection coefficient for SNP being selected in surrogate population 3)
```

Note that depending on the number of target/surrogate population used the number of columns will be different from that shown here. The above assumes 3 surrogate ("surrogate.vec") populations and one target ("pop.vec") population. 

## AdaptMixSimulator

`AdaptMixSimulator` simulates allele frequencies of a single admixed population using a forward simulator, incorporating selection post-admixture or pre-admixture in one or more of the source populations at a single SNP (randomly chosen from the supplied data -- see below). As with `AdaptMix`, you feed the program real data from a target population (i.e. the population assumed to have experienced selection) and surrogates to each of the sources of that target population. You also specify how well the target population's real sources match each surrogate (i.e. "drift").

The simulator prints out:
      (a) Genotpye data for each surrogate population (i.e. the real data).
      (b) Simulated genotype data for the target population.

The program allows for missingness in surrogates and targets, but at least some target individuals and individuals from each surrogate population should be non-missing at each SNP.

### Simulation protocol: 

(A) selection pre-admixture:
* (A1) sample simulated source population's allele frequencies via beta according to surrogate allele frequencies plus user-defined drift
* (A2) for simulated source population(s) `k` under selection, for selected SNP forward simulate randomly-mating population of size `N_k` diploids for `G_k` generations with selection (incorporating dominance, etc)
* (A3) sample admixed individual's data, at selected SNPs + X (user-specified) neutral SNPs, via binomial model according to target individuals' admixture proportions * simulated source population's allele frequencies

(B) selection post-admixture:
(NOTE: Currently this simulates an equal number of generations `G` of selection in each source pop, and then combines the final frequencies to make the admixed target population. This is likely not ideal, but is what was used in the paper. The challenge is adding selection to an admixed target pop while keeping the admixture proportions of each target individual the same. Another possibility is to simulate one admixed population per target individual, matching that individual's admixture proportions, adding selection to that population and sampling one individual from it, but this is computationally demanding.)
* (B1) sample simulated source populations' allele frequencies via beta according to surrogate allele frequencies plus drift
* (B2) at X neutral SNPs, sample admixed individual's data via binomial according to target individuals' admixture proportions * simulated source population's allele frequencies
* (B3) at selected SNP:
  * (B3i) forward simulate `K` randomly-mating (source) pops of size `N_k` diploids for `G_k` generations with selection (incorporating dominance, etc) 
  * (B3ii) sample admixed inds' data via binomial model according to target individuals' admixture proportions * simulated source population's allele frequencies (with no additional added drift)


### Running AdaptMixSimulator

`AdaptMixSimulator.R` takes three files as input, and is run from a command line in the following way, providing four file names: 

```
R < AdaptMixSimulator.R parameter.input.file genotypes.input.filenames id.file output.filePREFIX --no-save > screenoutput.out
```

E.g. using the provided example files: 

```
R < AdaptMixSimulator.R example/PEL_simulation_paramfile.txt example/PEL_REFs_ALLCHR_chr.txt example/PEL_REFs.ids.txt example/PEL_SIM_ALLCHR --no-save > screenoutput.out
```

### Input files

The first three file names are input files, described below. The last is the prefix for the output file name. "screenoutput.out" saves any technical output from AdaptMixSimulator (with no results), which may contain some helpful information (see "Strategies" below).

"genotypes.input.filenames" and "id.file" are in the same format as described in "AdaptMix", so we do not describe them here. 

The new file, "parameter.input.file" has the following format:

```
selection.post-admixture?: [0,1]
sel.coeff: [0.0,....]
sel.type: [additive,dominant,multiplicative,recessive]
target.pop: [name]
surrogate.pops: [name1,name2,...]
sources.with.selection.preadmixture: [0,1]
generations.selection.each.source: [1,...]
pop.size.sources: [1,...]
drift.btwn.surrogates.and.sources: [0.0,...,1.0]
num.neutral.snps: [1,...]
max.startfrequency.selected.snp: [0.0,...,1.0]
infer.source.freq.using.target.data: [0,1]
divide.into.runs.ofXX.inds(to.reduce.RAM): [1,...]
```

ALL FIELDS must be entered. 

The first line asks whether post-admixture selection should be simulated (1=YES, 0=simulate pre-admixture selection).

The second line specifies the selection coefficient (s) for the single SNP under selection. Under the various models, this is the increased probability of having offspring per each genotype class (count of selected alleles carried):

``` 
           	0	1	  2
additive   	0	s	  2s
dominant	0	s	  s
recessive	0	0	  s
multiplicative	0	s	  s^2+2s
```

**If you want to simulate all SNPs to have NO selection, choose 0 for `sel.coeff`.**

The third line specifies the model of selection (above) -- choose one. 

The fourth line specifies the target population, which will be matched for individual-specific admixture proportions.  

The fifth line specifies the surrogate populations, whose data will be used to simulate the sources of the target population. 

The sixth line allows you to specify which sources are undergoing pre-admixture selection (`1`=YES, `0`=NO; specify for each surrogate population). NOTE that all sources will be "1" if `selection.post-admixture?: 1`.

The seventh line specifies the number of generations of selection for each source. If selection is post-admixture, the mean of these values is used for each source.

The eighth line specifies the effective population size of each source over time (i.e. `N_k`).

The ninth line specifies the "drift" between each real surrogate population and its corresponding source, i.e. how well the source data likely reflects the surrogate data. (See Strategies below for advice on how to specify this.)

The tenth line specifies the number of neutral SNPs to simulate -- these are based on randomly choosing SNPs across the supplied data. Specify `"all"` to use all the data; in this case all but one SNP (randomly chosen to undergo selection) will be simulated as neutral.

The eleventh line specifies that the SNP chosen to be selected will have starting frequency <= this value for the selected allele. The starting frequency (i.e. before selection begins) is determined by the simulated source allele frequencies.

The twelfth line specifies whether you want to infer source allele frequencies using the target genotype data and admixture proportions. This may take longer to run, but can be useful for determining whether the drift values you specified, i.e. measuring how well each surrogate matches its corresponding source, are appropriate. (See Strategies below.)

The thirteenth (last) line specifies whether to only use a certain number of individuals at a time when simulating. This will slow the program down, but may be necessary if the memory used is high. NOTE: if specifiying `"infer.source.freq.using.target.data: 1"`, memory used may still be high after lowering this number. If that is the case, you can go into the R code and try changing `"snps.perrun"` to a smaller value.

### Output files

Three output files will be produced:

1. <output.filePREFIX>.haps -- new genotype input files for running in AdaptMix

2. <output.filePREFIX>.id -- a new ID input file for running in AdaptMix

3. <output.filePREFIX>.truth -- notes parameters used to simulate

### Strategies

One issue is trying to simulate sources that appropriately match their corresponding surrogates, as this can notably affect the power and interpretation of `AdaptMix`. The parameter to toggle this correlation is `"drift.btwn.surrogates.and.sources"`. 

One way to determine the best values to use is to specify `"infer.source.freq.using.target.data: 1"`, which will use only target data (genotypes and admixture proportions) to infer allele frequencies of each source, using a simple binomial model. You can then find (A) the correlation between each source's inferred allele frequency and its corresponding surrogate. You can compare this correlation to (B) that between each simulated source's allele frequencies and its corresponding surrogate.

After `AdaptMixSimulator` is finished running with `"infer.source.freq.using.target.data: 1"`, the bottom of <screenoutput.out> (e.g. `"pic screenoutput.out"`) will provide these correlations, with `"cor.sims"` corresponding to (B) and `"cor.truth"` corresponding to (A). One strategy then, is to toggle `"drift.btwn.surrogates.and.sources"` until these values align. (Though note that for sources where the admixture fractions overall are small in the admixed target individuals, e.g. "YRI" in the provided example, these values may never align well, as there is too little data to reliably calculate `"cor.truth"`. In such cases, the value of drift is likely less important, as there is little influence -- i.e. admixture -- from this source anyway.)

However, an issue with this is that the simulator adds additional drift onto the targets by simulating their genotypes from a binomial that uses the simulated source frequencies. Therefore, you may want to also ensure that the inferred drift from `AdapMix` is similar between your real data and your simulated data. See below.

## Worflow

Below, we outline the workflow for generating a set of neutral SNPs using the forward simulation setting implemented in `AdaptMixSimulator`. In Mendoza-Revilla et al., we used this strategy to estimate a false positive rate (FPR) cutoff, which helped determine whether a given SNP was selected (either pre- or post-admixture). In the `"parameter.input.file"`, we set the `num.neutral.snps parameter` to `"all"` and the `drift.btwn.surrogates.and.sources` parameter to `0.26`, `0.16`, and `0.3`. The example parameter file is shown below.

```
selection.post-admixture?: 1
sel.coeff: 0.5
sel.type: additive
target.pop: PEL
surrogate.pops: CHB IBS YRI
sources.with.selection.preadmixture: 1 0 0
generations.selection.each.source: 10 0 0
pop.size.sources: 10000 10000 10000
drift.btwn.surrogates.and.sources: 0.26 0.16 0.3
num.neutral.snps: all
max.startfrequency.selected.snp: 0.5
infer.source.freq.using.target.data: 1
divide.into.runs.ofXX.inds(to.reduce.RAM): 700
```

Then, running `AdaptMixSimulator` with these commands

```
R < AdaptMixSimulator.R example/PEL_simulation_paramfile.txt example/PEL_REFs_genotypes_files.txt example/PEL_REFs_ids.txt example/PEL_SIM_ALLCHR --no-save > screenoutput.out
```
will give you at the bottom of `screenoutput.out` estimates similar to these:
 
 ```
  source drift.val cor.sims cor.truth
 CHB    0.26      0.9116   0.9114   
 IBS    0.16      0.9289   0.9271   
 YRI    0.3       0.864    0.6994   
```

In this case, the simulated and observed correlations given the drift values are quite similar, except for YRI. As noted above, YRI has contributed little to the admixture process, so modeling the drift may be less important. Note that the specific `drift.btwn.surrogates.and.sources` values may differ greatly from those observed here with your specific dataset, and you may need to adjust these values to find the appropriate settings.

Next, as noted above, we also check that the inferred drift from `AdapMix` is similar between the real data and the simulated data. We first prepare the data and parameters files to run `AdaptMix`:

```
gzip example/PEL_SIM_ALLCHR.haps
echo example/PEL_SIM_ALLCHR.haps.gz > example/PEL_SIM_ALLCHR_genotypes_files.txt
```

Then, we run `AdaptMix`:

```
R < run_AdaptMix.R example/PEL_REFs_paramfile.txt example/PEL_SIM_ALLCHR_genotypes_files.txt example/PEL_SIM_ALLCHR.idfile.txt example/PEL_SIM_ALLCHR_adaptmix.txt --no-save > screenoutput.out
```

The header of example/PEL_SIM_ALLCHR_adaptmix.txt shows the following drift estimates:

```
drift.est [0,0.05) [0.05,0.1) [0.1,0.15) [0.15,0.2) [0.2,0.25) [0.25,0.3) [0.3,0.35) [0.35,0.4) [0.4,0.45) [0.45,0.51)
PEL 0.000497 0.00875 0.015658 0.019456 0.025934 0.02948 0.03389 0.038473 0.037685 0.041418
```

which were close to those observed when running `AdaptMix` on the real data:

```
drift.est [0,0.05) [0.05,0.1) [0.1,0.15) [0.15,0.2) [0.2,0.25) [0.25,0.3) [0.3,0.35) [0.35,0.4) [0.4,0.45) [0.45,0.51)
PEL 0.0001 0.007029 0.012254 0.016864 0.020187 0.024984 0.024494 0.02731 0.033527 0.035545
```

We note that the simulated drifts are larger, which will likely make the cutoff more conservative.

Finally, we can estimate a cutoff for the `AdaptMix` scores by reading the output and calculating a quantile.

```
df <- read.table("example/PEL_SIM_ALLCHR_adaptmix.txt", skip=2, header=TRUE)
> quantile(df$log10.pval.target.1,probs=c(0.99,0.995,0.999,0.9999))
     99%    99.5%    99.9%   99.99% 
1.824360 2.078520 4.877608 9.477474 

```

In Mendoza-Revilla et al., we used a conservative FPR cutoff of `5e-5` (i.e., probs=0.99995), which we deemed appropriate for our dataset, considering the surrogate populations and number of SNPs tested.

As in our study, based on the estimated drifts and a given FPR for your specific dataset, it is possible to estimate the power to detect selection given a specified selection coefficient, as well as to correctly assign a selection scenario (i.e., selection pre- or post-admixture). We note that in our study, only SNPs with an estimated FPR cutoff of `5e-5` were considered selected. We classified these SNPs as either pre- or post-admixture only if the I-score (see above) was lower than `0.001`.

## Citation
Mendoza-Revilla, Javier, et al. "Disentangling signatures of selection before and after European colonization in Latin Americans." Molecular biology and evolution 39.4 (2022): msac076.
