# The Effect of Ancestry in Forensic Record Matching

This repository contains all necessary files for running the analysis presented in the paper *The Effect of Ancestry in Forensic Record Matching*. The included record matching scripts compute a match score between a single SNP profile and a single STR profile, using allele frequency, genotype probability, and true genotype files as inputs.

## Overview

### Data
This analysis utilizes data from 2504 individuals typed at 18 CODIS markers, along with STRs and SNPs present within a 1 megabase window centered at each CODIS locus. Phased and unphased VCFs are available within the `data/` directory.

- SNPs/Imputed STRs: [Gymrek Lab - SNPSTR Imputation](https://gymreklab.com/2018/03/05/snpstr_imputation.html)

### Pipeline Steps

1. **Run STRUCTURE**: Estimate ancestry and allele frequencies.
2. **Run CLUMPP**: Combine replicate STRUCTURE runs.
3. **Generate Reference Panels & Test Sets**: Use `genotype_imputation/imputation_functions.R`.
4. **Phase & Impute Genotypes**: Use BEAGLE and `imputation_functions.R`.
5. **Extract Allele Frequencies**: Use `record_matching/af_extraction.R`.
6. **Extract Genotype Probabilities**: Use `record_matching/gp_extraction.R`.
7. **Compute Record Matching Scores**: Use `record_matching/rm_functions.R`.
8. **Calculate Match Accuracies**: Use `record_matching/rm_functions.R`.

## Running STRUCTURE

STRUCTURE clusters individuals based on inferred ancestry. This is performed using 3245 STR loci, as well as the 18 CODIS loci from 1000 Genomes. Each run is replicated 100 times for consistency.

**Download STRUCTURE:**  
[STRUCTURE v2.3.4](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html)

**Run STRUCTURE using provided parameter files:**
```bash
for i in {1..100}; do
    ./structure -m structure/codis/mainparams.boot -e structure/codis/extraparams30${i}.boot -K 3 -i structure/codis/codis.str -o structure/codis/output/codis_{i}.str
    ./structure -m structure/str/mainparams.boot -e structure/str/extraparams30${i}.boot -K 3 -i structure/str/all_str.str -o structure/str/output/str_{i}.str
done
```

## Running CLUMPP

CLUMPP resolves label switching issues across multiple STRUCTURE replicates.

**Download CLUMPP:**  
[CLUMPP](https://rosenberglab.stanford.edu/clumpp.html)

**Generate CLUMPP input files:**
```bash
bash clumpp/get_codis_ind_file.sh
bash clumpp/get_codis_pop_file.sh
bash clumpp/get_str_ind_file.sh
bash clumpp/get_str_pop_file.sh
```

**Run CLUMPP using provided parameter files:**
```bash
./CLUMPP clumpp/codis/codis_ind_params.txt
./CLUMPP clumpp/codis/codis_pop_params.txt
./CLUMPP clumpp/str/str_ind_params.txt
./CLUMPP clumpp/str/str_pop_params.txt
```

Inspect the CLUMPP output files manually to determine which cluster labels correspond to each superpopulation. Then, create permutation files by copying the permutation table from each `*miscfile.txt` file and adding a header with the correct superpopulation order.

## Generating Reference Panels

Locate the `ref_gen` function in `genotype_imputation/imputation_functions.R`. This function generates a list of individuals for the reference panel in one iteration of genotype imputation and allele frequency calculation. The remaining 1000 Genomes profiles become the test set.

Use the following functions to create reference/test set VCFs:
```r
# Generate VCF files for reference panel
generate_ref_vcfs()

# Generate VCF files for test set
generate_testset_vcfs()
```

## Running Genotype Imputation

Run the imputation using:
```r
impute_codis_gts()
```
This function uses BEAGLE to phase SNP profiles and impute masked CODIS alleles.

Evaluate imputation accuracy using:
```r
imputation_accuracy()
```
This function compares imputed CODIS alleles with true genotypes and returns an accuracy score.

## Extracting Allele Frequencies

For the ancestry estimation scheme, generate raw allele frequency files:
```r
get_raw_af_ancestry_codis()
get_raw_af_ancestry_str()
```
To extract allele frequencies for a specific matching scheme, use:
```r
get_af_*scheme*()
```
This writes allele frequency files to the designated output directory.

## Extracting Genotype Probabilities

Use the `get_gp()` function in `record_matching/gp_extraction.R` to extract genotype probability files from BEAGLE output:
```r
get_gp()
```
This function writes genotype probability files to the specified output directory.

## Running Record Matching

To compute match scores between SNP and STR profiles, use:
```r
rm()
```
This function calculates match scores using allele frequency, genotype probability, and true CODIS genotype files.

After computing match scores for all profile pairs, construct a match score matrix and calculate match accuracy:
```r
match_rm()
```
This function applies four matching algorithms and reports match accuracy.

## Contact
For any questions, contact the authors.

