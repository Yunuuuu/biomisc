# multiplication works

    Code
      run_motif_fisher(maf[, c("Tumor_Sample_Barcode", "Chromosome", "Start_position",
        "Reference_Allele", "Tumor_Seq_Allele1")], ref_genome = "hg38")
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
      
      i Mapping seqnames of `mut_data` to `ref_genome` (UCSC)
    Condition
      Warning in `valid.GenomicRanges.seqinfo()`:
      GRanges object contains 8 out-of-bound ranges located on sequences
        chr9, chr21, chr22, chr6, chr19, and chr15. Note that ranges located on
        a sequence whose length is unknown (NA) or on a circular sequence are
        not considered out-of-bound (use seqlengths() and isCircular() to get
        the lengths and circularity flags of the underlying sequences). You can
        use trim() to trim these ranges. See ?`trim,GenomicRanges-method` for
        more information.
    Message
      i Defining signature motif frequency and background context
      i Performing one-way Fisher's test
    Output
                  sample     A     C     G     T   AAA   AAC   AAG   AAT   ACA   ACC
                  <char> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int>
      1: TCGA-DK-A1A6-01  4663  3665  3650  4627   559   215   313   336   306   201
      2: TCGA-DK-A1A6-06  3827  3022  3045  3841   482   172   273   267   228   167
           ACG   ACT   AGA   AGC   AGG   AGT   ATA   ATC   ATG   ATT   CAA   CAC
         <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int>
      1:    58   253   347   243   296   222   244   228   272   348   283   274
      2:    54   220   265   200   258   187   182   182   218   288   241   224
           CAG   CAT   CCA   CCC   CCG   CCT   CGA   CGC   CGG   CGT   CTA   CTC
         <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int>
      1:   338   317   315   274    66   295    45    67    62    44   206   281
      2:   280   238   264   222    57   237    39    61    53    39   163   232
           CTG   CTT   GAA   GAC   GAG   GAT   GCA   GCC   GCG   GCT   GGA   GGC
         <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int>
      1:   349   292   288   165   289   186   241   219    52   259   262   220
      2:   282   255   230   147   232   160   209   179    48   203   226   191
           GGG   GGT   GTA   GTC   GTG   GTT   TAA   TAC   TAG   TAT   TCA   TCC
         <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int> <int>
      1:   273   208   159   168   261   209   295   166   178   255   337   255
      2:   233   168   129   133   216   176   240   131   140   206   267   212
           TCG   TCT   TGA   TGC   TGG   TGT   TTA   TTC   TTG   TTT
         <int> <int> <int> <int> <int> <int> <int> <int> <int> <int>
      1:    45   317   269   237   333   323   279   276   280   542
      2:    36   274   235   186   274   256   235   237   242   484
         background_mut_freq signature_mut_freq background_conetxt signature_context
                       <int>              <int>              <num>             <num>
      1:                 347                 37               7315              1270
      2:                 297                 31               6067              1041
         n_mutations non_signature_mut_freq enrichment fisher_pvalue        OR
               <int>                  <int>      <num>         <num>     <num>
      1:         406                    369  0.6141619     0.9997628 0.5681535
      2:         336                    305  0.6083150     0.9994654 0.5627179
             ci.up ci.low
             <num>  <num>
      1: 0.4149178    Inf
      2: 0.3982515    Inf

