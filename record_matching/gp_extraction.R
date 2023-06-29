get_gp = function(vcf.f, marker, out.dir){
    # This function extracts the genotype probabilities from a BEAGLE-imputed
    # VCF recoded with vcftools.
    # Input:
    #   vcf.f: a VCF file in .recode.vcf format with BEAGLE-imputed CODIS 
    #       genotypes. See genotype_imputation/imputation_functions.R
    #   marker: a string with the name of the CODIS marker represented in vcf.f
    #   out.dir: a directory for the genotype probability files
    # Output:
    #   This function directly writes genotype probability files to the out
    #   directory

    out.pre = file.path(out.dir, paste('imp_str_', marker, sep=''))

    vcftools.str = paste("vcftools --vcf ",
                            vcf.f,
                            " --extract-FORMAT-info GP",
                            " --out ", out.pre, sep = "")
    system(vcftools.str)
}