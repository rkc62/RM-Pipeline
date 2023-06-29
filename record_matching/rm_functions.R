library(plyr)

true.greater.than.false = function(mat){
    trues = diag(mat)
    falses = c(as.numeric(mat[lower.tri(mat, diag = FALSE)]), as.numeric(mat[upper.tri(mat, diag = FALSE)]))
    maxfalse = max(falses)
    return(trues > maxfalse)
}


ind.to.gt = function(index) {
    # Genotype index to genotype. The output is sorted by a1 first and then a2.
    # For each genotype, a1 >= a2. 
    # It's in the order or (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), (3,0), ...
    a1 = floor((sqrt(8 * index - 7) -1) / 2)
    a2 = index - a1 * (a1 + 1) / 2 - 1
    return(data.frame('a1' = a1, 'a2' = a2))
}


gt.to.ind = function(gt){
    # Genotype to genotpe index. 
    # Input is a vector (or data frame) of length two (or dim 1x2).
    # Output is an index (scalar) that correspond to the input genotype.
    gt = sort(gt, decreasing=T, na.last=T) 
    ind = gt[1] * (gt[1] + 1) / 2 + 1 + gt[2]
    return(ind)
}

read.anc.al.freq = function(filename, allele.count, correct = T){
    # This function reads files containing ancestry allele frequencies and 
    # returns allele frequencies in training set allele frequency files.
    # Input:
    #   filename: filename containing allele frequencies from training set.
    #             the extension must be 'frq'.
    #   correct: T/F, default = T
    #            routine from Doc. correcting for zero allele frequencies
    # Output:
    #   result: list with the following items:
    #       $al.freq: allele frequency vector
    #       $min.af: minimum nonzero allele frequency value BEFORE correction
        
    al.freq.tab = read.table(file=filename, header=TRUE,)
    filesplit = strsplit(filename, split = '/')
    marker = strsplit(filesplit[[1]][length(filesplit[[1]])], split = '_')[[1]][1]
    max.alleles = allele.count[allele.count$marker == marker,]$alleles
    alleles = al.freq.tab$allele
    freq.char = c()
    for (i in 1:max.alleles){
        if (!(i-1) %in% alleles){
            freq.char[i] = 0
        } else {
            freq.char[i] = al.freq.tab[al.freq.tab$allele == (i-1),]$frequency
        }
    }
    al.freq = freq.char
    

    smallest.af = min(al.freq[al.freq > 0])
    n.zeros = sum(al.freq == 0)
    zero.indices = which(al.freq == 0)
    al.freq[al.freq == 0] = (smallest.af ^ 2) / 1000 
    al.freq = al.freq / sum(al.freq)
    
    result = list('al.freq'=al.freq, 'min.af'=smallest.af)
    return(result)
}

read.al.freq = function(filename, correct=T) {
    # This function reads files containing allele frequencies and return
    # allele frequencies in training set allele frequency files.
    # Input:
    #   filename: filename containing allele frequencies from training set.
    #             the extension must be 'frq'.
    #   correct: T/F, default = T
    #            routine from Doc. correcting for zero allele frequencies
    # Output:
    #   result: list with the following items:
    #       $al.freq: allele frequency vector
    #       $min.af: minimum nonzero allele frequency value BEFORE correction
    
    al.freq.tab = read.table(file=filename, head=FALSE, as.is=TRUE, skip=1)
    freq.char = as.character(al.freq.tab[-(1:4)])
    al.freq = as.numeric(matrix(unlist(strsplit(freq.char, ":")), nrow = 2)[2,])
    n.al = al.freq.tab[1, 3]
    stopifnot(n.al == length(al.freq)) # check to see if we extracted everything
    

    smallest.af = min(al.freq[al.freq > 0])
    n.zeros = sum(al.freq == 0)
    zero.indices = which(al.freq == 0)
    al.freq[al.freq == 0] = (smallest.af ^ 2) / 1000 
    al.freq = al.freq / sum(al.freq)
    
    result = list('al.freq'=al.freq, 'min.af'=smallest.af)
    return(result)
}


read.gt = function(filename, phased=T) {
    # This function reads files containing genotypes from BEAGLE.
    # Note that the allele indexing starts from "0" from BEAGLE.
    # Input: 
    #   filename: name of the file containing genotype data
    #   phased: Whether genotype in the filename has been phased or not.
    #           Default is "T". If phased, separator is "|", and if not,
    #           the separator is "/".
    # Output:
    #   Genotypes of individuals in the file. 
    #   row: individuals, 
    #   columns: allele 1 and allele 2, NOT sorted
    
    # stopifnot(strsplit(filename, "\\.")[[1]][2] == 'GT') # check input format
    
    if (phased) {
        split.char = '|'
    } else {
        split.char = '/'
    }
    
    gt.str = read.table(filename, head=T, as.is=T, na.strings=c("NA","-9"))
    gt.data = matrix(as.numeric(unlist(sapply(gt.str[-c(1,2,3)], strsplit, 
                                               split=split.char, fixed=T))),
                      ncol=2, byrow=T)
    return(gt.data)
}

read.gp = function(filename) {
    # This function reads files containing genotype probabilities from BEAGLE.
    # Input: 
    #   filename: name of the file containing genotype probability data
    # Output:
    #   Genotypes probabilities of individuals in the file. 
    #   row: individuals, 
    #   columns: genotype probabilities of each genotype
    #            the ordering of genotype probabilities follows the ones in
    #            the genotype file.
    # Note: The first and the second element of gp.str are chromosom number 
    #       and position, respectively.
    
    # stopifnot(strsplit(filename, "\\.")[[1]][2] == 'GP') # check input format
    gp.str = read.table(filename, head=T, as.is=T, na.strings = c("NA","-9"))
    gp.data = matrix(as.numeric(unlist(sapply(gp.str[-c(1,2)], strsplit, 
                                               split=",", fixed=T))),
                      nrow=dim(gp.str)[2] - 2, byrow=T)
    rownames(gp.data) = names(gp.str)[-c(1,2)]
    return(gp.data)
}

rm = function(gp.f, af.f, gt.f, ancestry = FALSE){
    # This function computes the record matching score between one SNP profile
    # and one STR profile for one CODIS locus. The test hypothesis is that the
    # profiles are drawn from the same individual. The null hypothesis is that 
    # the profiles are drawn from unrelated individuals. Formally, this function
    # calculates log(P(R_A | S_B)) - log(P(R_A)), where R_A is the CODIS STR 
    # profile of individual A and S_B is the flanking SNP profile of individual 
    # B. The first log is the conditional term that can be extracted from the 
    # BEAGLE genotype probability and the second log is the unconditional term
    # that can be extracted from the reference panel or STRUCTURE depending on 
    # the record matching scheme.
    # Input:
    #   gp.f: file with the genotype probabilities from BEAGLE for the CODIS 
    #       locus
    #   af.f: file with the allele frequencies for the CODIS locus based on 
    #       the reference panel or STRUCTURE output
    #   gt.f: file with the true genotype for the STR profile
    #   ancestry: boolean that represents whether the record match scheme is 
    #       ancestry estimation, which has differently formatted allele
    #       frequency files
    # Output:
    #   lambda: the record match score for the input profiles



    #===================
    ## TRUE genotype of A, i.e. known STR genotype, NOT from imputation.
    ## The missing STR genotype is represented with NA

    gen.str.true = read.gt(gt.f, phased = TRUE) # make sure the GT is from phased; if not, change to FALSE
    str.true.ind = apply(gen.str.true, 1, gt.to.ind)

    #===================
    ## Load allele frequencies from reference panel for A
    if (ancestry){
        al.freq = read.anc.al.freq(af.f)
    } else {
        al.freq = read.al.freq(af.f)
    }
    
    al.freq.vec = al.freq$al.freq
    smallest.af = al.freq$min.af
    n.al = length(al.freq.vec)
    n.gen = n.al * (n.al + 1) / 2

    ## Construct HW genotype frequency matrix of A based on 
    ## the allele frequencies for A
    ## This is for the 2nd term in the match score computation. P(R_A)
    outer.mat = outer(al.freq.vec, al.freq.vec)
    hw.gf.mat = 2 * outer.mat - diag(diag(outer.mat))
    hw.gf.vec = hw.gf.mat[!lower.tri(hw.gf.mat)] 
    hw.gf.min = apply(cbind(hw.gf.vec, rep(0.0005, length(hw.gf.vec))), 
                    1, min) #minimum values to be taken for imputation probabilities

    unobserved.mat = hw.gf.mat
    unobserved.mat[hw.gf.mat >= (smallest.af^2)/10] = 0
    unobserved.mat[hw.gf.mat < (smallest.af^2)/10] = 2       
    unobserved.vec = unobserved.mat[!lower.tri(unobserved.mat)] 
    #vector with "2" everywhere with a genotype involving 
    #an allele unobserved in the training data.

    #===================
    ## Step 1: Load BEAGLE imputed genotype probabilities of B, 
    #  P(R_B | S_B) and set imputation probabilities for anything 
    #  unobserved in training set to be 0, which will make the algorithm 
    #  reset them to the small value in gf.min

    val.gp = read.gp(gp.f)
    val.gp[ , unobserved.vec == 2] = 0

    ## Step 2: Set imputation probabilities of genotypes unobserved in 
    #  the training set to be min(training set GF, 0.005) and normalize
    #  imputation probability vector to sum to 1. 
    hw.gf.min.mat = matrix(rep(hw.gf.min, 1), nrow=1, byrow=TRUE)
    scale.mat = diag((1 - rowSums(hw.gf.min.mat * (val.gp == 0))) 
                / rowSums(val.gp))
    val.gp.scale = scale.mat %*% val.gp
    val.gp.scale[val.gp == 0] = hw.gf.min.mat[val.gp == 0]

    match.prob.mat = val.gp.scale
    if (any(is.na(match.prob.mat))) {
        print(paste("marker:", m.name))
        stop("something went wrong in processing const.match.prob.mat function.")
    }

    ## Compute log-likelihood ratio
    conditional = log(apply(match.prob.mat, 1, "[", i=str.true.ind))
        
    train.gf = matrix(rep(hw.gf.vec, 1), nrow=1, byrow=T) 
    unconditional = log(apply(train.gf, 1, "[", i=str.true.ind))

    lambda = conditional - unconditional
    return(lambda)
}

match_rm = function(match.matrix.f){
    # This function performs four match methods on a square matrix of record
    # matching outputs. The methods are one-to-one matching, SNP query matching,
    # STR query matching, and one pair matching. The values returned are 
    # accuracies that represent the success rate fo each matching approach.
    # Input:
    #   match.matrix: a square matrix where rows represent the SNP profiles of
    #       the test set and columns represent the STR profiles of the test set
    # Output:
    #   accuracies: a vector with entries corresponding to one-to-one matching
    #       (Hungarian method), SNP query, STR query, and one pair matching,
    #       in that order

    library(clue)

    mat = read.table(match.matrix.f, header = TRUE, row.names = 1)
    if (!nrow(mat) == ncol(mat)){
        stop('RM score matrix is not square')
    }

    #Remove anyone with missing data at all loci.
    mat = as.matrix(mat[(rowSums(mat^2) != 0),(rowSums(mat^2) != 0)])
    matched.LSAP = solve_LSAP(mat - min(mat) + 1, maximum = TRUE)

    hungarian.acc = mean(matched.LSAP == 1:length(matched.LSAP))

    #  One-to-many accuracy computation
    pickSTR.acc = mean(apply(mat, 2, which.max) == 1:length(mat[,1]))
    pickSNP.acc = mean(apply(mat, 1, which.max) == 1:length(mat[,1]))

    #  Needle-in-haystack (one match) accuracy computation
    onematch.acc = mean(true.greater.than.false(mat))

    # Write accuracy vector
    accuracies = c(hungarian.acc, pickSTR.acc, pickSNP.acc, onematch.acc)

    return(accuracies)
}