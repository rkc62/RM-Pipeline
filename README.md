# The Effect of Ancestry in Forensic Record Matching

This repository contains all necessary files for running the analysis presented in the paper *The Effect of Ancestry in Forensic Record Matching*. The included record matching scripts compute a match score between a single SNP profile and a single STR profile, using allele frequency, genotype probability, and true genotype files as inputs.

### Data
This analysis utilizes data from 2504 individuals typed at 18 CODIS markers, along with STRs and SNPs present within a 1 megabase window centered at each CODIS locus. Phased and unphased VCFs are available within the `data/` directory.

- SNPs/Imputed STRs: [Gymrek Lab - SNPSTR Imputation](https://gymreklab.com/2018/03/05/snpstr_imputation.html)

### Pipeline Steps

1. **Run STRUCTURE**: Use STRUCTURE to get ancestry and allele frequency estimates for the ancestry estimation record match scheme.
2. **Run CLUMPP**: Combine replicate STRUCTURE runs. Run distruct to visualize output using the instructions below.
3. **Generate Reference Panels & Test Sets**: Use `genotype_imputation/imputation_functions.R` to generate reference panels and test sets for the null, known superpopulation, empirical, ancestry estimation (CODIS), and ancestry estimation (all STR) schemes.
4. **Phase & Impute Genotypes**: Use BEAGLE and `genotype_imputation/imputation_functions.R` to phase test set SNP profiles and impute test set masked STRs for each record matching scheme.
5. **Extract Allele Frequencies**: Use `record_matching/af_extraction.R` to generate allele frequency files for each individual in the test set.
6. **Extract Genotype Probabilities**: Use `record_matching/gp_extraction.R` to generate genotype probability files for each individual in the test set.
7. **Compute Record Matching Scores**: Use `record_matching/rm_functions.R` to compute the record matching scores for each SNP profile/STR profile pairing in the test set.
8. **Downstream Analyses**: Use `downstream_analysis/calc_match_accuracy.R` to compute match accuracies across four matching methods. Use `downstream_analyis/calc_AUCs` to evaluate record match performance under each setting. Additional scripts for plotting figures and are described below.

## Running STRUCTURE
STRUCTURE clusters individuals based on inferred ancestry. This is performed using 3245 STR loci, as well as 18 CODIS loci from 1000 Genomes. Each run is replicated 100 times for consistency.

**Download STRUCTURE:**  
[STRUCTURE v2.3.4](https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html)

**Run STRUCTURE using provided parameter files:**
```bash
for i in {1..10}; do
    ./structure -m structure/input/mainparams_codisK3.boot -e structure/input/extraparams_codis/extraparams${i}.boot -K 3 -i structure/input/codis.str -o structure/output/k3/codis_{i}.str_f
    ./structure -m structure/input/mainparams_allstrK3.boot -e structure/input/extraparams_allstr/extraparams${i}.boot -K 3 -i structure/input/all_str.str -o structure/output/k3/allstr_{i}.str_f
    ./structure -m structure/input/mainparams_allsnpK3.boot -e structure/input/extraparams_allsnp/extraparams${i}.boot -K 3 -i structure/input/all_snp.str -o structure/output/k3/allstr_{i}.str_f
done
```

## Running CLUMPP
CLUMPP resolves label switching issues across multiple STRUCTURE replicates.

**Download CLUMPP:**  
[CLUMPP](https://rosenberglab.stanford.edu/clumpp.html)

**Generate CLUMPP input files:**
```bash
bash structure_output_manip.sh
```

**Run CLUMPP using provided parameter files:**
```bash
./CLUMPP clumpp/param_files/codis_ind_params.txt
./CLUMPP clumpp/param_files/codis_pop_params.txt
./CLUMPP clumpp/param_files/allstr_ind_params.txt
./CLUMPP clumpp/param_files/allstr_pop_params.txt
./CLUMPP clumpp/param_files/allsnp_ind_params.txt
./CLUMPP clumpp/param_files/allsnp_pop_params.txt
```

Inspect the CLUMPP output files manually to determine which cluster labels correspond to each superpopulation. Then, create permutation files by copying the permutation table from each `*miscfile.txt` file and adding a header with the correct superpopulation order.

**Extract and format individual allele frequencies from CLUMPP average as input for ancestry-aware record matching**
```bash
Rscript structure/extract_structure_afs.R
Rscript structure/format_rm_afs.R
```

## Running distruct
distruct is used to produce visualizations of STRUCTURE outputs.

Download distruct:
https://rosenberglab.stanford.edu/distruct.html

Run distruct:
./distructLinux1.1 -d distruct/drawparams_codis.txt
./distructLinux1.1 -d distruct/drawparams_allstr.txt
./distructLinux1.1 -d distruct/drawparams_allsnp.txt

## Generating Reference Panels
Locate the `ref_gen()` function in `genotype_imputation/imputation_functions.R`. This function generates a list of individuals for the reference panel in one iteration of genotype imputation and allele frequency calculation. The remaining 1000 Genomes profiles become the test set.

Then, use the following functions to create reference/test set VCFs:
```r
# Generate VCF files for reference panel
generate_ref_vcfs()

# Generate VCF files for test set
generate_testset_vcfs()
```
These functions will take in the id lists and a data file path to generate VCF files subset to include only individuals in the reference panel or only individuals in the test set. We omit the total collection of VCF files to reduce overall repo size.

## Running Genotype Imputation
Run the imputation by using the `impute_codis_gts()` function in `genotype_imputation/imputation_functions.R`. This function will use BEAGLE to phase the test set SNP profiles and impute the masked CODIS alleles on each haplotype.

Run the imputation using:
```r
impute_codis_gts()
```
The output appears in data/imputed_genotype_probabilities.

Analyze the accuracy of the imputation by using the `imputation_accuracy()` function in `genotype_imputation/imputation_functions.R`. This function will compare the imptued CODIS alleles to one true CODIS genotype and return a score that indicates imputation accuracy.

Evaluate imputation accuracy using:
```r
imputation_accuracy()
```
The output appears in data/imputation_accuracy.

Plot the results of imputation accuracy using:
```bash
Rscript downstream_analysis/accuracy_plotter.R
Rscript downstream_analysis/locus_accuracy_plotter.R
python downstream_analysis/mixture_panel_accuracy_plotter.py
```

## Extracting Allele Frequencies
If the record matching scheme is ancestry estimation, first generate the raw allele frequency files from the CLUMPP output using either of the following functions based on the desired STRUCTURE input data (18 CODIS loci or all 3245 STRs):
```r
get_raw_af_ancestry_codis()
get_raw_af_ancestry_str()
```
To extract allele frequencies for a specific matching scheme, use the following command to writes allele frequency files to the designated output directory:
```r
get_af_*scheme*()
```
The output appears in data/allele_frequencies.

## Extracting Genotype Probabilities
Use the `get_gp()` function in `record_matching/gp_extraction.R` to extract genotype probability files from BEAGLE output:
```r
get_gp()
```
This function writes genotype probability files to the specified output directory.

## Running Record Matching
To compute match scores between SNP and STR profiles, locate and run the `rm()` function in `record_matching/rm_functions.R` to calcuate the record match score using one allele frequency file, one genotype probability file, and the true CODIS genotype file.
```r
rm()
```
The results for each record matching scheme appear in rm_scores/<scheme>.

## Calculate record match accuracies
After computing match scores for all profile pairs, construct a match score matrix and calculate match accuracy:
```bash
Rscript downstream_analysis/calc_match_acc.R
```
The results appear in data/rm_scores/match_accuracy

## Generate ROC curves and calculate AUC values
Use the total match score matrices to generate ROC curves and evaluate performance with AUC. Replicates with are selected for ROC curve and KDE visualization. AUC is used as a tie breaker when 
```bash
Rscript calc_AUCs.R
Rscript plot_scheme_ROCs.R
Rscript plot_scheme_KDEs.R
Rscript plot_admpop_ROCs.R
Rscript plot_admpop_KDEs.R
```

## Contact
For any questions, contact the authors.

