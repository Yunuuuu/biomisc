# multiplication works

    Code
      suppressWarnings(run_motif_fisher(maf[, c("Tumor_Sample_Barcode", "Chromosome",
        "Start_position", "Reference_Allele", "Tumor_Seq_Allele2")], ref_genome = "hg38"))
    Message
      
      Attaching package: 'BiocGenerics'
      
      The following objects are masked from 'package:stats':
      
          IQR, mad, sd, var, xtabs
      
      The following objects are masked from 'package:base':
      
          Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
          as.data.frame, basename, cbind, colnames, dirname, do.call,
          duplicated, eval, evalq, get, grep, grepl, intersect, is.unsorted,
          lapply, mapply, match, mget, order, paste, pmax, pmax.int, pmin,
          pmin.int, rank, rbind, rownames, sapply, setdiff, sort, table,
          tapply, union, unique, unsplit, which.max, which.min
      
      
      Attaching package: 'S4Vectors'
      
      The following object is masked from 'package:utils':
      
          findMatches
      
      The following objects are masked from 'package:base':
      
          I, expand.grid, unname
      
      
      Attaching package: 'Biostrings'
      
      The following object is masked from 'package:base':
      
          strsplit
      
      Mapping seqnames of `mut_data` to `ref_genome` ("UCSC")
      Defining signature motif frequency and background context
      i Using "C>T" and "C>G" to define background_mut_freq
      i Using "TCA", "TCT", "TGA", and "AGA" to define signature_mut_freq
      i Using "C" and "G" to define background_context
      i Using "TCA", "TCT", "TGA", and "AGA" to define signature_context
      Performing one-way Fisher's test
    Output
                  sample signature_mut_freq background_mut_freq signature_context
                  <char>              <int>               <int>             <num>
      1: TCGA-DK-A1A6-01                 33                 312               988
      2: TCGA-DK-A1A6-06                 28                 266               845
         background_context non_signature_mut_freq n_mutations enrichment
                      <num>                  <int>       <int>      <num>
      1:               5520                    279         312  0.5909374
      2:               4770                    238         266  0.5942074
         fisher.pvalue        or    ci.low ci.up
                 <num>     <num>     <num> <num>
      1:     0.9998297 0.5426123 0.3881218   Inf
      2:     0.9994671 0.5465323 0.3789297   Inf

