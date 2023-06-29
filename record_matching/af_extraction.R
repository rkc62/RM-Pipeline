get_af = function(ref.f, marker, out.dir){
    # This function extracts the allele frequency file (.frq) from a 
    # compressed VCF (.vcf.gz) for a single CODIS marker.
    # Input:
    #   ref.f: compressed vcf (.vcf.gz) file containing the CODIS marker for all
    #       samples in the reference panel for a given test individual
    #   marker: a string representing the CODIS marker of interest
    #   out.dir: the desired directory for the output .frq file
    # Output:
    #   This function writes the .frq file directly to the output directory

    out.pre <- file.path(out.dir, paste('ref_', marker, sep=''))

    vcftools.str <- paste("vcftools --gzvcf ", ref.f,
                          " --snp ", marker,
                          " --freq --out ", out.pre, sep = "")
    system(vcftools.str)
}

get_raw_af_ancestry_codis = function(codis.struc.out.dir, codis.permute.f, n.iter = 100, marker.info.f, allele.data.f, out.dir){
    # This function generates the raw allele frequency files based on n.iter 
    # STRUCTURE runs with CODIS-based input and subsequent CLUMPP aligning.
    # Output files are used to generate the weighted p-hat vectors for each 
    # individual in the test set.
    # Input:
    #   codis.struc.out.dir: a file path to the directory containing all the 
    #       codis STRUCTURE run outs. CLUMPP does not average allele 
    #       frequencies, so this function will do so manually. Out files should
    #       be named codis_$iteration$.str_f, where $iteration$ ranges from 1 
    #       to n.iter
    #   codis.permute.f: a file with cluster permutations for each STRUCTURE run
    #       with columns labeled by superpopulation (AFR, EAS, or EUR) and row
    #       i corresponding to the ith STRUCTURE run
    #   n.iter: the number of STRUCTURE runs to average over (default 100)
    #   marker.info.f: a file with CODIS marker information (see 
    #   data/marker_info_1kgp.txt)
    #   allele.data.f: a file with a row per CODIS marker where the first column
    #       is the name of the marker and the subsequent columns are all the 
    #       alleles for that marker present in 1kgp (see data/allele_info.txt) 
    #   out.dir: a directory for the raw allele frequency files
    # Output:
    #   This function writes allele frequency files directly to the desired out
    #   directory

    permute.table = read.table(codis.permute.f, header = TRUE)
    marker.info = read.table(marker.info.f, header = TRUE)
    allele.info = read.table(allele.info.f, header = FALSE)

    # Define the lines from the STRUCTURE output corresponding to each CODIS 
    # marker
    marker.lines = list('CSF1PO' = c(2594, 2604), 'D10S1248' = c(2609, 2618), 
        'D12S391' = c(2623, 2690), 'D13S317' = c(2695, 2709), 'D18S51' = c(2714, 2738),
        'D19S433' = c(2743, 2765), 'D1S1656' = c(2770, 2801), 'D22S1045' = c(2806, 2816),
        'D2S1338' = c(2821, 2864), 'D2S441' = c(2869, 2885), 'D3S1358' = c(2890, 2912),
        'D5S818' = c(2917, 2933), 'D7S820' = c(2938, 2946), 'D8S1179' = c(2951, 2976),
        'FGA' = c(2981, 3002), 'TH01' = c(3007, 3013), 'TPOX' = c(3018, 3024),
        'vWA' = c(3029, 3051))


    # Extract allele frequencies from each STRUCTURE iteration separately
    for (iter in 1:n.iter){
        print(paste0('Extracting allele frequencies from STRUCTURE run ', iter))
        iter.dir = file.path(out.dir, iter)
        if (!dir.exists(iter.dir)){
            dir.create(iter.dir)
        }

        struc.f = file.path(codis.struc.out.dir, paste0('codis_', iter, '.str_f'))

        # Use the permutation table to set the correct cluster index labels for
        # each iteration
        AFR.cluster.index = as.numeric(permute.table$AFR[iter])
        EAS.cluster.index = as.numeric(permute.table$EAS[iter])
        EUR.cluster.index = as.numeric(permute.table$EUR[iter])

        for (m.ind in 1:length(marker.info$str)){
            m = marker.info$str[m.ind]
            allele.info.row = which(allele.info[,1] == m)

            # Extract the allele frequency information from the STRUCTURE run 
            # for the current locus. This will only have data for alleles 
            # present in 1kgp, so some alleles may be missing
            pull.af.cmd.str = paste0('awk \'NR==', marker.lines[[m]][1], ',NR==', marker.lines[[m]][2], '\' ', struc.file)
            marker.frequencies = system(pull.af.cmd.str, intern = TRUE)
            n.alleles = length(marker.frequencies)
            marker.af.df = data.frame('allele' = rep(NA, n.alleles), 
                                      'sequence' = rep(NA, n.alleles), 
                                      'AFR' = rep(NA, n.alleles), 
                                      'EAS' = rep(NA, n.alleles), 
                                      'EUR' = rep(NA, n.alleles)
                                      )

            # For each represented allele, extract the cluster frequencies from
            # the STRUCTURE output
            for (allele.ind in 1:n.alleles){
                allele.data = strsplit(marker.frequencies[allele.ind], split = " ")
                allele.data = allele.data[[1]][!(allele.data[[1]] == "")]
                allele = allele.data[1]
                allele.info.col = as.numeric(allele) + 2
                allele.seq = allele.info[allele.info.row, allele.info.col]
                AFR.af = allele.data[AFR.cluster.index + 1]
                EAS.af = allele.data[EAS.cluster.index + 1]
                EUR.af = allele.data[EUR.cluster.index + 1]
                marker.af.df[allele.ind, ] = c(as.numeric(allele), allele.seq, AFR.af, EAS.af, EUR.af)
            }
            # Sort the allele data frame to be in numeric order by allele
            sorted.marker.af.df = marker.af.df[order(as.numeric(marker.af.df$allele)),]
            write.table(sorted.marker.af.df, file = file.path(out.dir, paste0(m, '_', iter, '.txt')), col.names = TRUE, row.names = FALSE, quote = FALSE)
        }
    }

    # Average the allele frequencies over all STRUCTURE runs
    for (m.ind in 1:length(marker.info$str)){
        m = marker.info$str[m.ind]

        for (iter in 1:n.iter){
            raw.af.f = file.path(out.dir, paste0(m, '_', iter, '.txt'))
            raw.af = read.table(raw.af.f, header = TRUE)
            if (iter == 1){
                AFR.afs = raw.af$AFR
                EAS.afs = raw.af$EAS
                EUR.afs = raw.af$EUR
            } else {
                AFR.afs = AFR.afs + raw.af$AFR 
                EAS.afs = EAS.afs + raw.af$EAS
                EUR.afs = EUR.afs + raw.af$EUR
            }
        }

        AFR.afs = AFR.afs / n.iter
        EAS.afs = EAS.afs / n.iter
        EUR.afs = EUR.afs / n.iter

        n.alleles = nrow(raw.af)
        marker.af.df = data.frame('allele' = rep(NA, n.alleles), 
                                  'sequence' = rep(NA, n.alleles), 
                                  'AFR' = rep(NA, n.alleles), 
                                  'EAS' = rep(NA, n.alleles), 
                                  'EUR' = rep(NA, n.alleles)
                                 )

        for (i in 1:n.alleles){
            marker.af.df[i,] = c(raw.af$allele[i], raw.af$sequence[i], AFR.afs[i], EAS.afs[i], EUR.afs[i])
        }
        write.table(marker.af.df, file = file.path(out.dir, paste0(m, '_af.txt')), col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
}

get_raw_af_ancestry_str = function(str.struc.out.dir, str.permute.file){
    # This function generates the raw allele frequency files based on n.iter 
    # STRUCTURE runs with STR-based input and subsequent CLUMPP aligning.
    # Output files are used to generate the weighted p-hat vectors for each 
    # individual in the test set.
    # Input:
    #   str.struc.out.dir: a file path to the directory containing all the 
    #       STR STRUCTURE run outs. CLUMPP does not average allele 
    #       frequencies, so this function will do so manually. Out files should
    #       be named str_$iteration$.str_f, where $iteration$ ranges from 1 
    #       to n.iter
    #   str.permute.f: a file with cluster permutations for each STRUCTURE run
    #       with columns labeled by superpopulation (AFR, EAS, or EUR) and row
    #       i corresponding to the ith STRUCTURE run
    #   n.iter: the number of STRUCTURE runs to average over (default 100)
    #   marker.info.f: a file with CODIS marker information (see 
    #   data/marker_info_1kgp.txt)
    #   allele.data.f: a file with a row per CODIS marker where the first column
    #       is the name of the marker and the subsequent columns are all the 
    #       alleles for that marker present in 1kgp (see data/allele_info.txt) 
    #   out.dir: a directory for the raw allele frequency files
    # Output:
    #   This function writes allele frequency files directly to the desired out
    #   directory

    permute.table = read.table(str.permute.f, header = TRUE)
    marker.info = read.table(marker.info.f, header = TRUE)
    allele.info = read.table(allele.info.f, header = FALSE)

    # Define the lines from the STRUCTURE output corresponding to each CODIS 
    # marker
    marker.lines = list('CSF1PO' = c(40124, 40134), 'D10S1248' = c(40139, 40148), 
        'D12S391' = c(40153, 40220), 'D13S317' = c(40225, 40239), 'D18S51' = c(40244, 40268),
        'D19S433' = c(40273, 40295), 'D1S1656' = c(40300, 40331), 'D22S1045' = c(40336, 40346),
        'D2S1338' = c(40351, 40394), 'D2S441' = c(40399, 40415), 'D3S1358' = c(40420, 40442),
        'D5S818' = c(40447, 40463), 'D7S820' = c(40468, 40476), 'D8S1179' = c(40481, 40506),
        'FGA' = c(40511, 40532), 'TH01' = c(40537, 40543), 'TPOX' = c(40548, 40554),
        'vWA' = c(40559, 40581))



    # Extract allele frequencies from each STRUCTURE iteration separately
    for (iter in 1:n.iter){
        print(paste0('Extracting allele frequencies from STRUCTURE run ', iter))
        iter.dir = file.path(out.dir, iter)
        if (!dir.exists(iter.dir)){
            dir.create(iter.dir)
        }

        struc.f = file.path(str.struc.out.dir, paste0('str_', iter, '.str_f'))

        # Use the permutation table to set the correct cluster index labels for
        # each iteration
        AFR.cluster.index = as.numeric(permute.table$AFR[iter])
        EAS.cluster.index = as.numeric(permute.table$EAS[iter])
        EUR.cluster.index = as.numeric(permute.table$EUR[iter])

        for (m.ind in 1:length(marker.info$str)){
            m = marker.info$str[m.ind]
            allele.info.row = which(allele.info[,1] == m)

            # Extract the allele frequency information from the STRUCTURE run 
            # for the current locus. This will only have data for alleles 
            # present in 1kgp, so some alleles may be missing
            pull.af.cmd.str = paste0('awk \'NR==', marker.lines[[m]][1], ',NR==', marker.lines[[m]][2], '\' ', struc.file)
            marker.frequencies = system(pull.af.cmd.str, intern = TRUE)
            n.alleles = length(marker.frequencies)
            marker.af.df = data.frame('allele' = rep(NA, n.alleles), 
                                      'sequence' = rep(NA, n.alleles), 
                                      'AFR' = rep(NA, n.alleles), 
                                      'EAS' = rep(NA, n.alleles), 
                                      'EUR' = rep(NA, n.alleles)
                                      )

            # For each represented allele, extract the cluster frequencies from
            # the STRUCTURE output
            for (allele.ind in 1:n.alleles){
                allele.data = strsplit(marker.frequencies[allele.ind], split = " ")
                allele.data = allele.data[[1]][!(allele.data[[1]] == "")]
                allele = allele.data[1]
                allele.info.col = as.numeric(allele) + 2
                allele.seq = allele.info[allele.info.row, allele.info.col]
                AFR.af = allele.data[AFR.cluster.index + 1]
                EAS.af = allele.data[EAS.cluster.index + 1]
                EUR.af = allele.data[EUR.cluster.index + 1]
                marker.af.df[allele.ind, ] = c(as.numeric(allele), allele.seq, AFR.af, EAS.af, EUR.af)
            }
            # Sort the allele data frame to be in numeric order by allele
            sorted.marker.af.df = marker.af.df[order(as.numeric(marker.af.df$allele)),]
            write.table(sorted.marker.af.df, file = file.path(out.dir, paste0(m, '_', iter, '.txt')), col.names = TRUE, row.names = FALSE, quote = FALSE)
        }
    }

    # Average the allele frequencies over all STRUCTURE runs
    for (m.ind in 1:length(marker.info$str)){
        m = marker.info$str[m.ind]

        for (iter in 1:n.iter){
            raw.af.f = file.path(out.dir, paste0(m, '_', iter, '.txt'))
            raw.af = read.table(raw.af.f, header = TRUE)
            if (iter == 1){
                AFR.afs = raw.af$AFR
                EAS.afs = raw.af$EAS
                EUR.afs = raw.af$EUR
            } else {
                AFR.afs = AFR.afs + raw.af$AFR 
                EAS.afs = EAS.afs + raw.af$EAS
                EUR.afs = EUR.afs + raw.af$EUR
            }
        }

        AFR.afs = AFR.afs / n.iter
        EAS.afs = EAS.afs / n.iter
        EUR.afs = EUR.afs / n.iter

        n.alleles = nrow(raw.af)
        marker.af.df = data.frame('allele' = rep(NA, n.alleles), 
                                  'sequence' = rep(NA, n.alleles), 
                                  'AFR' = rep(NA, n.alleles), 
                                  'EAS' = rep(NA, n.alleles), 
                                  'EUR' = rep(NA, n.alleles)
                                 )

        for (i in 1:n.alleles){
            marker.af.df[i,] = c(raw.af$allele[i], raw.af$sequence[i], AFR.afs[i], EAS.afs[i], EUR.afs[i])
        }
        write.table(marker.af.df, file = file.path(out.dir, paste0(m, '_af.txt')), col.names = TRUE, row.names = FALSE, quote = FALSE)
    }
}

get_af_ancestry = function(q.hat, raw.marker.dir, marker, out.dir){
    # This function generates an allele frequency file for one individual under 
    # the ancestry estimation record match scheme i.e. using a linear 
    # combination of each allele frequency cluster vector weighted by the
    # ancestry clustering vector q-hat to calculate the allele frequency.
    # Input:
    #   q.hat: a numerical vector with entries corresponding to ancestry 
    #       proportion distributed across the three representative 
    #       superpopulations in order AFR, EAS, EUR. Output from CLUMPP
    #   raw.marker.dir: a directory with the raw af files for each CODIS locus
    #       as generated in get_raw_af_ancestry_codis or get_raw_af_ancestry_str
    #   marker: a string with the name of the CODIS marker of interest
    #   out.dir: a directory for the weighted allele frequencies of the 
    #       individual described by q-hat
    # Output:
    #   This function directly writes allele frequency files to the out
    #   directory


    raw.af.table = read.table(file.path(raw.marker.dir, paste0(marker, '_af.txt')), header = TRUE)
    n.alleles = length(raw.af.table$allele)
    ind.af.df = data.frame('allele' = rep(NA, n.alleles), 'sequence' = rep(NA, n.alleles), 'frequency' = rep(NA, n.alleles))
    
    # Elements of q.hat must be in the order AFR, EAS, EUR
    names(q.hat) = c('AFR', 'EAS', 'EUR')
    
    for (allele.ind in 1:n.alleles){
        allele.data = raw.af.table[allele.ind,]
        AFR.af = allele.data$AFR
        EAS.af = allele.data$EAS
        EUR.af = allele.data$EUR
        allele = allele.data$allele
        allele.seq = allele.data$sequence[allele.ind]

        #Calculate the linear combination
        linear.combo = AFR.af * q.hat['AFR'] + EAS.af * q.hat['EAS'] + EUR.af * q.hat['EUR']
        ind.af.df[allele.ind, ] = c(allele, allele.seq, linear.combo)
    }

    write.table(ind.af.df, file = file.path(out.dir, paste0(marker, '_af.txt')), col.names = TRUE, row.names = FALSE, quote = FALSE)
}