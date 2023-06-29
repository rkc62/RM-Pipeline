nested_sampling_ref_gen = function(ref.pool.f){
    # This function generates nested reference panels from size 25 to 375 with 
    # step size 25 using an input file that lists the names of all individuals 
    # to be sampled from. Samples are chosed uniformly at random from the input 
    # pool and all larger sample sizes contain all of the individuals at each 
    # smaller sample size, making this a nested sampling method. 100 iterations 
    # of nested sampling are done at once.
    # Input:
    #   ref.pool.f: file with ids of pool of individuals to be sampled from with
    #       one name on each line. ref.pool.f should comprise all individuals
    #       belonging each of AFR, EAS, and EUR superpopulations
    # Output:
    #   ref.ppl: a list of lists of lists. The innermost list contains 100 
    #       vectors with the ids of the people in the reference panel for a 
    #       certain sample size. Each sample size is a named list.

    ppl = read.table(ref.pool.f, header = FALSE)[,1]
    n.iter = 100
    n.max = 375
    n.list = c(seq(375, 25, by=-25))

    ref.ppl = list()

    ## Iterate through the list of sample sizes and execute nested sampling
    for (n in n.list){
        
        ## When n is n.max, do the initial sample
        if (n == n.max) {
            min.ind = 1
            max.ind = length(ppl)
            samples = replicate(n.iter, sample(min.ind:max.ind, size = n, 
                                replace = FALSE))
        } else {
            tmp.mat = matrix(NA, nrow = n, ncol = n.iter)
            for (j in 1:n.iter){
                tmp.mat[, j] = sample(samples[, j], size = n, replace = FALSE)
            }
            if(NA %in% tmp.mat) stop("something went wrong")
            samples = tmp.mat
            rm(tmp.mat)
        }
        # Add current n list to the output
        sample.size.list = list()
        for (j in 1:n.iter) {
            sampled.ppl = ppl[samples[,as.numeric(j)]]
            if (length(unique(sampled.ppl))<length(sampled.ppl)){
                stop("duplicates detected")
                }
            sample.size.list = append(sample.size.list, list(sampled.ppl))
        }

        ref.ppl = append(ref.ppl, list(sample.size.list))
    }

    names(ref.ppl) = as.character(n.list)
    return(ref.ppl)
}

null_ref_gen = function(n = 375, ref.pool.f){
    # This function generates a list of n reference panel individuals sampled 
    # from those present in the file ref.pool.f. The method for sampling is 
    # 'null' i.e. sampling uniformly at random 
    # Input:
    #   n: number of individuals to be sampled
    #   ref.pool.f: file with ids of pool of individuals to be sampled from with
    #       one name on each line. ref.pool.f should comprise all individuals
    #       belonging each of AFR, EAS, and EUR superpopulations
    # Output:
    #   ref.ppl: a vector with the ids of sampled indivdiduals constituting the 
    #       reference panel

    ref.pool = read.table(ref.pool.f, header = FALSE)[,1]
    ref.ppl = sample(ref.pool, size = n, replace = FALSE)
    return(ref.ppl)
}

known_superpop_ref_gen = function(n = 375, ref.pool.f){
    # This function generates a list of n reference panel individuals sampled 
    # from those present in the file ref.pool.f. The method for sampling is 
    # 'known_superpopulaion' i.e. sampling uniformly at random from a single 
    # representative superpopulation
    # Input:
    #   n: number of individuals to be sampled
    #   ref.pool.f: file with ids of pool of individuals to be sampled from with
    #       one name on each line. ref.pool.f should comprise all individuals
    #       belonging to a single superpopulation
    # Output:
    #   ref.ppl: a vector with the ids of sampled indivdiduals constituting the 
    #       reference panel

    ref.pool = read.table(ref.pool.f, header = FALSE)[,1]
    ref.ppl = sample(ref.pool, size = n, replace = FALSE)
    return(ref.ppl)
}

empirical_ref_gen = function(ref.pool.f, empirical.accuracy.vector, pop.info.f){
    # This function generates a list of n reference panel individuals sampled 
    # from those present in the file ref.pool.f. The method for sampling is 
    # 'empirical' i.e. sampling such that the reference panel comprises a  
    # specific contribution ratio of the representative superpopulations 
    # Input:
    #   ref.pool.f: file with ids of pool of individuals to be sampled from with
    #       one name on each line. ref.pool.f should comprise all individuals
    #       belonging to a single superpopulation
    #   empirical.accuracy.vector: a numerical vector representing how many 
    #       individuals to sample from each superpopulation in the order
    #       c(AFR, EAS, EUR)
    #   pop.info.f: file with 1kgp sample info. See data/popinfo_1kgp.txt

    # Output:
    #   ref.ppl: a vector with the ids of sampled indivdiduals constituting the 
    #       reference panel

    ref.pool = read.table(ref.pool.f, header = FALSE)[,1]
    pop.info = read.table(pop.info.f, header = TRUE)

    ## Define subpopulations in each representative superpopulation
    AFR.subpops = c('ESN', 'GWD', 'LWK', 'MSL', 'YRI')
    EAS.subpops = c('CDX', 'CHB', 'CHS', 'JPT', 'KHV')
    EUR.subpops = c('CEU', 'FIN', 'GBR', 'IBS', 'TSI')

    ## Extract number of individuals to sample from each superpopulation
    n.AFR.ppl = empirical.accuracy.vector[1]
    n.EAS.ppl = empirical.accuracy.vector[2]
    n.EUR.ppl = empirical.accuracy.vector[3]

    ## Generate the superpopulation-specific reference pools from the 1kgp data
    AFR.pool = pop.info$samples[pop.info$pop %in% AFR.subpops]
    EAS.pool = pop.info$samples[pop.info$pop %in% EAS.subpops]
    EUR.pool = pop.info$samples[pop.info$pop %in% EUR.subpops]

    ## Sample individuals according to the empirical.accuracy.vector
    AFR.ref.ppl = sample(AFR.pool, size = n.AFR.ppl, replace = FALSE)
    EAS.ref.ppl = sample(EAS.pool, size = n.EAS.ppl, replace = FALSE)
    EUR.ref.ppl = sample(EUR.pool, size = n.EUR.ppl, replace = FALSE)

    ref.ppl = c(AFR.ref.ppl, EAS.ref.ppl, EUR.ref.ppl)
    return(ref.ppl)
}

ancestry_ref_gen = function(n = 375, ref.pool.f, q.hat, pop.info.f){
    # This function generates a list of n reference panel individuals sampled 
    # from those present in the file ref.pool.f. The method for sampling is 
    # 'ancestry' i.e. sampling such that the proportion of  
    # representative superpopulations in the reference panel match the ancestry
    # vector q.hat extracted from STRUCTURE. 
    # Input:
    #   n: number of individuals to be sampled
    #   ref.pool.f: file with ids of pool of individuals to be sampled from with
    #       one name on each line. ref.pool.f should comprise all individuals
    #       belonging to a single superpopulation
    #   q.hat: a numerical vector summing to one representing the proportions 
    #       each superpopulation should have in the reference panel. 
    #       Proportions should be in the order c(AFR, EAS, EUR)
    #   pop.info.f: file with 1kgp sample info. See data/popinfo_1kgp.txt

    # Output:
    #   ref.ppl: a vector with the ids of sampled indivdiduals constituting the 
    #       reference panel


    ref.pool = read.table(ref.pool.f, header = FALSE)[,1]
    pop.info = read.table(pop.info.f, header = TRUE)

    ## Define subpopulations in each representative superpopulation
    AFR.subpops = c('ESN', 'GWD', 'LWK', 'MSL', 'YRI')
    EAS.subpops = c('CDX', 'CHB', 'CHS', 'JPT', 'KHV')
    EUR.subpops = c('CEU', 'FIN', 'GBR', 'IBS', 'TSI')

    ## Extract the proportions of each superpopulation to be included in the
    ## reference panel
    AFR.prop = q.hat[1]
    EAS.prop = q.hat[2]
    EUR.prop = q.hat[3]

    ## Turn proportions into integers and ensure that the proportions are as 
    ## close as possible to the desired proportions
    total.to.sample = n * c(AFR.prop, EAS.prop, EUR.prop)
    floor.to.sample = floor(total.to.sample)
    while (sum(floor.to.sample) < n){
        most.loss.index = which((total.to.sample - floor.to.sample) == max(total.to.sample - floor.to.sample))
        floor.to.sample[most.loss.index] = floor.to.sample[most.loss.index] + 1
    }

    ## Generate the superpopulation-specific reference pools from the 1kgp data
    AFR.pool = pop.info$samples[pop.info$pop %in% AFR.subpops]
    EAS.pool = pop.info$samples[pop.info$pop %in% EAS.subpops]
    EUR.pool = pop.info$samples[pop.info$pop %in% EUR.subpops]

    ## Sample individuals according to q.hat
    AFR.ref.ppl = sample(AFR.pool, size = floor.to.sample[1], replace = FALSE)
    EAS.ref.ppl = sample(EAS.pool, size = floor.to.sample[2], replace = FALSE)
    EUR.ref.ppl = sample(EUR.pool, size = floor.to.sample[3], replace = FALSE)

    ref.ppl = c(AFR.ref.ppl, EAS.ref.ppl, EUR.ref.ppl)
    return(ref.ppl)
}

generate_ref_vcfs = function(ref.ppl.f, codis.vcf.dir, out.dir){
    # This function takes in a list of ids and generates compressed VCFs 
    # (.vcf.gz format) for the 18 CODIS loci and surrounding SNPs with samples
    # corresponding to the input id list using bcftools. Effectively subsets the
    # master VCFs to include only the reference panel individuals.
    # Input:
    #   ref.ppl.f: file with the ids of reference panel individuals, one name 
    #       per line
    #   codis.vcf.dir: file path to directory with *phased* CODIS vcf files.  
    #       This directory must have .vcf.gz files, tabix-indexed files .tbi, 
    #       and a marker information file (see data/marker_info_1kgp.txt).
    #   out.dir: file path to the output directory
    # Output:
    #   This function writes ouput vcfs directly to the output directory

    marker.info.f = file.path(codis.vcf.dir, 'marker_info_1kgp.txt')
    marker.info = read.table(marker.info.f, header = TRUE)

    codis.stump = '_halfwindow500000WithCODIS.vcf.gz'

    for (m.ind in 1:length(marker.info$str)){
        m = marker.info$str[m.ind]
        vcf.f = file.path(codis.vcf.dir, paste0(m, codis.stump))
        out.f = file.path(out.dir, paste0(m, codis.stump))

        bcftools.str = paste('bcftools view --samples-file', ref.ppl.f,
                                '--output', out.f, vcf.f,
                                sep = ' ')
        system(bcftools.str)
    }
}

generate_testset_vcfs = function(test.ppl.f, codis.vcf.dir, out.dir){
    # This function takes in a list of ids and generates compressed VCFs 
    # (.vcf.gz format) for the surrounding SNPs without CODIS genotypes with 
    # samples corresponding to the input id list using bcftools. Effectively 
    # subsets the master VCFs to include only the reference panel individuals.
    # Input:
    #   test.ppl.f: file with the ids of test set individuals, one name 
    #       per line
    #   codis.vcf.dir: file path to directory with *unphased* CODIS vcf files. 
    #       This directory must have .vcf.gz files, tabix-indexed files .tbi,
    #       and a marker information file (see data/marker_info_1kgp.txt).
    #   out.dir: file path to the output directory
    # Output:
    #   This function writes ouput vcfs directly to the output directory

    marker.info.f = file.path(codis.vcf.dir, 'marker_info_1kgp.txt')
    marker.info = read.table(marker.info.f, header = TRUE)

    snp.stump = '_halfwindow500000SNPsOnly.vcf.gz'

    for (m.ind in 1:length(marker.info$str)){
        m = marker.info$str[m.ind]
        vcf.f = file.path(codis.vcf.dir, paste0(m, snp.stump))
        out.f = file.path(out.dir, paste0(m, snp.stump))

        bcftools.str = paste('bcftools view --samples-file', ref.ppl.f,
                                '--output', out.f, vcf.f,
                                sep = ' ')
        system(bcftools.str)
    }
}

impute_codis_gts = function(ref.vcf.dir, test.vcf.dir, bgl.dir, out.dir){
    # This function takes in a reference panel VCF and a test set VCF and uses
    # BEAGLE to impute the missing CODIS genotypes for all samples in the test 
    # set. It also extracts the .GT format files for imputed CODIS genotypes.
    # Input:
    #   ref.vcf.dir: file path to directory with reference panel vcfs
    #   test.vcf.dir: file path to directory with test set vcfs (no CODIS)
    #   bgl.dir: file path to directory containing BEAGLE software, 
    #       plink-formatted chromosome map, and a marker information file (see
    #       data/marker_info_1kgp.txt)
    #   out.dir: file path to the output directory
    # Output:
    #   This function writes ouput files directly to the output directory.
    #   Output files are raw imputed vcfs (marker_imptued.vcf.gz), imputed CODIS
    #   only vcfs (STR_marker.recode.vcf), and imputed CODIS genotype files 
    #   (STR_marker.GT.FORMAT)

    map.dir = file.path(bgl.dir, 'plink.GRCh37.map')
    beagle = file.path(bgl.dir, 'beagle.22Jul22.46e.jar')
    marker.info.f = file.path(bgl.dir, 'marker_info_1kgp.txt')
    marker.info = read.table(marker.info.f, header = TRUE)

    snp.stump = '_halfwindow500000SNPsOnly.vcf.gz'
    codis.stump = '_halfwindow500000WithCODIS.vcf.gz'

    log.dir = file.path(out.dir, 'log')
        if (!dir.exists(log.dir)) {
            dir.create(log.dir)
        }

    for (m.ind in 1:length(marker.info$str)){
        m = marker.info$str[m.ind]
        out.vcf = file.path(out.dir, paste0(m, '_imputed'))
        ref.vcf = file.path(ref.vcf.dir, paste0(m, codis.stump))
        to.imp.vcf = file.path(test.vcf.dir, paste0(m, snp.stump))

        beag.str = paste0("java -Xmx4g -jar ", beagle, 
                            " gt=", to.imp.vcf,
                            " out=", out.vcf,
                            " ref=", ref.vcf, 
                            " gp=true impute=true iterations=14", 
                            " map=", map.dir, '/', "plink.chr",
                            as.character(marker.info$chr[m.ind]), 
                            ".GRCh37.map", " nthreads=",
                            as.character(bgl.threads))
        system(beag.str)

        vcftools.str.1 <- paste("vcftools --gzvcf ", out.dir, '/',
                                m,  '_imputed.vcf.gz',
                                " --snp ", m, " --recode --recode-INFO-all", 
                                " --out ", out.dir, '/',
                                "STR_", m, sep = "")
        system(vcftools.str.1)
                    
        vcftools.str.2 <- paste("vcftools --gzvcf ", out.dir, '/',
                                "STR_", m,  
                                ".recode.vcf --extract-FORMAT-info GT",
                                " --out ", out.dir, '/',
                                "STR_", m, sep = "")
        system(vcftools.str.2)
    }
    system(paste('mv', file.path(out.dir, '*.log'), log.dir))
}

imputation_accuracy = function(imp.gt.f, true.gt.f){
    # This function calculates the imputation accuracy for a given CODIS locus
    # by comparing the imputed gt to the true gt for all individuals in the
    # test set.
    # Input:
    #   imp.gt.f: file (.GT.FORMAT) with imputed CODIS gts for all individuals
    #       in the test set. Contains only one CODIS locus
    #   true.gt.f: file (.GT.FORMAT) with imputed CODIS gts for all individuals
    #       in the 1kgp data set. Contains only one CODIS locus
    # Output:
    #   score: a number in (0,1) representing the imputation accuracy at the
    #       input locus averaged across all individuals in the test set

    imp.gt.table = read.table(imp.gt.f, header = TRUE)
    true.gt.table = read.table(true.gt.f, header = TRUE)

    ## Extract test ids from the imputed gt file by taking the column names and
    ## ignoring the CHROM and POS columns
    test.ids = names(imp.gt.table)[!(names(imp.gt.table) == 'CHROM' | names(imp.gt.table) == 'POS')]
    
    score = 0
    for (test.id in test.ids){
        test.col = which(names(imp.gt.table) == test.id)
        true.col = which(names(true.gt.table) == test.id)
        
        test.gt = imp.gt.table[1, test.col]
        true.gt = true.gt.table[1, true.col]

        GT.test = strsplit(as.character(test.gt), split = "\\|")[[1]]
        GT.true = strsplit(as.character(true.gt), split = "\\|")[[1]]

        ## Score the accuracy by calulation proportion of overlapping alleles
        score = score + (length(intersect(GT.test, GT.true)) / 2)
    }

    score = score / length(test.ids)
    return(score)
}