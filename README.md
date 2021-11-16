# Selection in Latin Americans
This repository contains the scripts used to detect and classify signals of selection in Latin Americans published in Mendoza-Revilla Javier et al. 2021

## Modules needed
* R/4.0.2

## Preparing files
`AdaptMix` uses several files as input. 

ChromoPainter's PHASE file:
The first line of the file contains the number of reference (donors) haplotypes. The second line of the file contains the total number of donor and recipient
individuals. The third line contains the letter "P" followed by a vector of the basepair positions of each SNP. The remaining lines of the file contain the genetic variation information of each donor and recipient haplotype, with each row a new haplotype and each column the allelic type at each biallelic SNP, in the same order as the "positions" line. NO misssing data is allowed.

ChromoPainter's ID file:
Each row is ordered to match the rows of the PHASE input file. There are three columns per row, with the first row giving the individual identifier, the second column giving the individual's population label and the third column an indicator for whether the individual should not be included in the analysis (use "0" to specify NOT to include the given individual).

ADMIXTURE Q file:
Each row is an individual and columns represent the proportions for each ancestry.

PLINK FAM file:
Sample information file (accompanying a .bed and .bim file) used to run ADMIXTURE.

## Run AdaptMix

The `AdaptMix` command line is as follows:

```
Rscript run_AdaptMix.R [target_admixed_pop_id] [prefix_HAPS_file] [postfix_HAPS_file] [ID_file] [ADMIXTURE_Q_file] [PLINK_FAM_FILE] [reference_pop_id_1,...,reference_pop_id_n] [output_file] 
```

Input arguments:
* 1: target_admixed_pop_id - Population ID of the target admixed population
* 2: prefix_HAPS_file Prefix of CP file (before "chr" string in file)
* 3: postfix_HAPS_file - Postfix of CP file (after "chr" string in file)
* 4: ID_file - ID file
* 5: ADMIXTURE_Q_file - Q file from ADMIXTURE (or any other software) 
* 6: PLINK_FAM_FILE - Fam file from PLINK used to run ADMIXTURE (ids should be the same order as the Q file)
* 7: reference_pop_id_1,...,reference_pop_id_n - Reference population IDs separated by a comma (the order of the reference population should MATCH the ancestries i.e. columns in the ADMIXTURE Q file)
* 8: output_file - Output file

Output:

The first two lines indicate the ID of the target admixed population and the drift estimate.
Line 3 contains the header of the output file:

* 1 - chrom (chromosome ID)
* 2 - pos (position in base pair)
* 3 - log10.pval.target.1 (-log10 P-value)
* 4 - obs.freq.target.1 (observed allele frequency in the target admixed population)
* 5 - exp.freq.target.1 (expected allele frequency in the target admixed population)
* 6 - AIC.neutral.target.1 (AIC for SNP being neutral)
* 7 - AIC.postadmix.target.1 (AIC for SNP being selected post-admixture)
* 8 - AIC.insurr.source1target.1 (AIC for SNP being selected in reference population 1)
* 9 - AIC.insurr.source2target.1 (AIC for SNP being selected in reference population 2)
* 10 - AIC.insurr.source3target.1 (AIC for SNP being selected in reference population 3)
* 11 - sel.postadmix.target.1 (selection coefficient for SNP being selected post-admixture)
* 12 - sel.insurr.source1target.1 (selection coefficient for SNP being selected in reference population 1)
* 13 - sel.insurr.source2target.1 (selection coefficient for SNP being selected in reference population 2)
* 14 - sel.insurr.source3target.1 (selection coefficient for SNP being selected in reference population 3)
* 15 - LRT.pval.postadmix.target.1 (LRT for SNP being selected post-admixture)
* 16 - LRT.pval.insurr.source1.target.1 (LRT for SNP being selected in reference population 1)
* 17 - LRT.pval.insurr.source2.target.1 (LRT for SNP being selected in reference population 2)
* 18 - LRT.pval.insurr.source3.target.1 (LRT for SNP being selected in reference population 3)

Note that depending on the number of target/reference population used the number of columns will be different from that shown here.

## Citation
Mendoza-Revilla Javier et al. "Disentangling signatures of selection before and after European colonization in Latin Americans." 
