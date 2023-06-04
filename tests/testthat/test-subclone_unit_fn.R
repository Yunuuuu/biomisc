test_that("subclone unit function works well", {
    # read test data ---------------------------------------
    mut.table <- data.table::fread("../testdata/clonal_dissection/SNV/LMS025.Exome.SNV.xls")
    segRData <- load("../testdata/clonal_dissection/ASCAT/D_LMS025_T1.ascat.seg.RData")
    seg.mat.copy <- get(segRData)
    # we read all mutatio data and cnv data of every region for a patient


    mut.table[, SampleID := vapply(
        strsplit(V1, ":", fixed = TRUE), "[[", character(1L), 1L
    )]
    mut.table[, V1 := NULL]
    mut.table[, mutation_id := NULL]
    data.table::setnames(
        mut.table, c(
            "D_LMS025_T1.ref_count",
            "D_LMS025_T1.var_count", "D_LMS025_T1.VAF"
        ),
        c("ref_counts", "alt_counts", "vaf")
    )
    data.table::setcolorder(
        mut.table,
        c("SampleID", "chr", "start", "stop", "ref", "var")
    )
    data.table::setDT(seg.mat.copy)
    data.table::setnames(seg.mat.copy, "sample", "SampleID")
    seg.mat.copy[, SampleID := vapply(
        strsplit(SampleID, "_", fixed = TRUE), "[[", character(1L), 2L
    )]
    data.table::setcolorder(
        seg.mat.copy, c("SampleID", "chr", "startpos", "endpos", "n.het", "cnTotal", "nMajor", "nMinor", "nAraw", "nBraw", "Ploidy", "ACF")
    )

    # create mutation ID ----------------------------------------
    mut.table[, mutation_id := paste(SampleID, chr, start, ref, sep = ":")]

    # test define_subclone_cn works -----------------------------
    # define subclone copy number
    seg.mat.phylo <- define_subclone_cn(
        seg.mat.copy, min_subclonal = 0.01
    )[order(SampleID, chr, startpos)]

    create.subclonal.copy.number <- function(seg.mat.copy,
                                             min.subclonal = 0.1) {
        seg.out <- seg.mat.copy
        seg.out$nMaj1 <- floor(as.numeric(seg.out$nAraw))
        seg.out$nMaj2 <- ceiling(as.numeric(seg.out$nAraw))
        seg.out$nMin1 <- floor(as.numeric(seg.out$nBraw))
        seg.out$nMin2 <- ceiling(as.numeric(seg.out$nBraw))


        seg.out$fracMaj1 <- as.numeric(seg.out$nMaj2) - as.numeric(seg.out$nAraw)
        seg.out$fracMaj2 <- 1 - as.numeric(seg.out$fracMaj1)
        seg.out$fracMin1 <- as.numeric(seg.out$nMin2) - as.numeric(seg.out$nBraw)
        seg.out$fracMin2 <- 1 - as.numeric(seg.out$fracMin1)


        # next, let's deal with the minimum subclonal
        seg.out$fracMaj1 <- ifelse(seg.out$fracMaj1 < min.subclonal, 0,
            ifelse(seg.out$fracMaj1 > (1 - min.subclonal), 1, seg.out$fracMaj1)
        )

        seg.out$fracMaj2 <- ifelse(seg.out$fracMaj2 < min.subclonal, 0,
            ifelse(seg.out$fracMaj2 > (1 - min.subclonal), 1, seg.out$fracMaj2)
        )

        seg.out$fracMin1 <- ifelse(seg.out$fracMin1 < min.subclonal, 0,
            ifelse(seg.out$fracMin1 > (1 - min.subclonal), 1, seg.out$fracMin1)
        )

        seg.out$fracMin2 <- ifelse(seg.out$fracMin2 < min.subclonal, 0,
            ifelse(seg.out$fracMin2 > (1 - min.subclonal), 1, seg.out$fracMin2)
        )

        # how much of the genome for each tumour region is subject to subclonal copy number
        for (region in unique(seg.mat.copy$SampleID))
        {
            seg.region <- seg.out[seg.out$SampleID == region, , drop = FALSE]
            prop.aber <- 1 - length(which((seg.region$fracMaj1 == 1 | seg.region$fracMaj1 == 0) & (seg.region$fracMin1 == 1 | seg.region$fracMin1 == 0))) / nrow(seg.region)
            # print(prop.aber)
        }

        # which ones work?

        seg.sorted <- seg.out[(seg.out$fracMaj1 == seg.out$fracMin1 | seg.out$fracMaj1 == seg.out$fracMin2) | (seg.out$fracMaj1 == 1) | (seg.out$fracMaj1 == 0) | (seg.out$fracMin1 == 1) | (seg.out$fracMin1 == 0), , drop = FALSE]
        seg.problem <- seg.out[which(!rownames(seg.out) %in% rownames(seg.sorted)), ]



        # let's divide the problem further
        seg.out.final <- c()

        for (a in 1:nrow(seg.problem))
        {
            seg.info <- seg.problem[a, , drop = FALSE]
            seg.info$fracA <- seg.info$fracMaj1 * seg.info$fracMin1
            seg.info$fracB <- seg.info$fracMaj1 - seg.info$fracA
            seg.info$fracC <- seg.info$fracMaj2 * seg.info$fracMin2
            seg.info$fracD <- seg.info$fracMaj2 - seg.info$fracC

            # let's see how this works.
            seg.info$nMaj_A <- seg.info$nMaj1
            # seg.info$nMaj_A <- seg.info$nMaj1
            # seg.info$nMin_A <- seg.info$nMin1
            seg.info$nMin_A <- seg.info$nMin1

            seg.info$nMaj_B <- seg.info$nMaj1
            # seg.info$nMaj2_B <- seg.info$nMaj1
            # seg.info$nMin1_B <- seg.info$nMin2
            seg.info$nMin_B <- seg.info$nMin2

            seg.info$nMaj_C <- seg.info$nMaj2
            # seg.info$nMaj2_C <- seg.info$nMaj2
            seg.info$nMin_C <- seg.info$nMin2
            # seg.info$nMin2_C <- seg.info$nMin2

            seg.info$nMaj_D <- seg.info$nMaj2
            # seg.info$nMaj2_D <- seg.info$nMaj2
            seg.info$nMin_D <- seg.info$nMin1
            # seg.info$nMin2_D <- seg.info$nMin1



            seg.out.final <- rbind(seg.out.final, seg.info)
        }

        # let's make the sorted easier to read
        seg.sorted$fracA <- NA
        seg.sorted$fracB <- NA
        seg.sorted$fracC <- NA
        seg.sorted$fracD <- NA
        seg.sorted$nMaj_A <- NA
        seg.sorted$nMaj_B <- NA
        seg.sorted$nMaj_C <- NA
        seg.sorted$nMaj_D <- NA
        seg.sorted$nMin_A <- NA
        seg.sorted$nMin_B <- NA
        seg.sorted$nMin_C <- NA
        seg.sorted$nMin_D <- NA

        # put the columns we do need in the right order
        cols.maj <- apply(cbind(seg.sorted$fracMaj1, seg.sorted$fracMaj2), 1, which.max)
        cols.min <- apply(cbind(seg.sorted$fracMin1, seg.sorted$fracMin2), 1, which.max)
        for (i in 1:nrow(seg.sorted))
        {
            poss.frac.major <- c(seg.sorted$fracMaj1[i], seg.sorted$fracMaj2[i])
            poss.frac.minor <- c(seg.sorted$fracMin1[i], seg.sorted$fracMin2[i])
            poss.cpn.major <- c(seg.sorted$nMaj1[i], seg.sorted$nMaj2[i])
            poss.cpn.minor <- c(seg.sorted$nMin1[i], seg.sorted$nMin2[i])
            cols.major <- cols.maj[i]
            cols.minor <- cols.min[i]

            fracA.major <- poss.frac.major[cols.major]
            fracA.minor <- poss.frac.minor[cols.minor]
            fracB.major <- poss.frac.major[which(!1:2 %in% cols.major)]
            fracB.minor <- poss.frac.minor[which(!1:2 %in% cols.minor)]

            seg.sorted$nMaj_A[i] <- poss.cpn.major[cols.major]
            seg.sorted$nMin_A[i] <- poss.cpn.minor[cols.minor]
            seg.sorted$nMaj_B[i] <- poss.cpn.major[which(!1:2 %in% cols.major)]
            seg.sorted$nMin_B[i] <- poss.cpn.minor[which(!1:2 %in% cols.minor)]

            seg.sorted$fracA[i] <- fracA.major
            seg.sorted$fracB[i] <- fracB.major
            if (fracA.major == 1) {
                seg.sorted$fracA[i] <- fracA.minor
                seg.sorted$fracB[i] <- fracB.minor
                seg.sorted$nMaj_B[i] <- poss.cpn.major[cols.major]
            }

            if (fracA.minor == 1) {
                seg.sorted$fracA[i] <- fracA.major
                seg.sorted$fracB[i] <- fracB.major
                seg.sorted$nMin_B[i] <- poss.cpn.minor[cols.minor]
            }
        }
        seg.final <- rbind(seg.out.final, seg.sorted)

        # let's order this correctly
        seg.final <- seg.final[order(seg.final$SampleID, seg.final$chr, seg.final$startpos), , drop = FALSE]
        # finally, let's choose the columns we want and the order we want
        colnames(seg.final)
        col.names <- c(
            "SampleID", "chr", "startpos", "endpos", "n.het", "cnTotal", "nMajor", "nMinor", "Ploidy",
            "ACF", "nAraw", "nBraw", "fracA", "nMaj_A", "nMin_A", "fracB", "nMaj_B", "nMin_B",
            "fracC", "nMaj_C", "nMin_C", "fracD", "nMaj_D", "nMin_D"
        )
        seg.final <- seg.final[, col.names]
        return(seg.final)
    }
    seg.mat.phylo2 <- create.subclonal.copy.number(
        as.data.frame(seg.mat.copy),
        min.subclonal = 0.01
    )
    data.table::setDT(seg.mat.phylo2)
    seg.mat.phylo2 <- seg.mat.phylo2[order(SampleID, chr, startpos)]
    testthat::expect_equal(seg.mat.phylo[
        , .SD, .SDcols = names(seg.mat.phylo2)
    ], seg.mat.phylo2)

    # test mut_match_cn works well ------------------------------
    # let's run PyClone, but correct for copy number before running
    # prepare pyclone table
    sample <- region <- seg.mat.phylo$SampleID[1]
    region.mut.table <- mut.table
    region.seg.copy <- seg.mat.copy[SampleID %in% region]
    region.seg.phylo <- seg.mat.phylo[SampleID %in% region]
    # just extract the segmented CNV for this sample
    # region.seg.copy <- seg.mat.phylo[SampleID %in% region]
    pyclone.table <- mut_match_cn(region.mut.table, seg.mat.phylo,
        on_sample = "SampleID", on_chr = "chr",
        mut_pos = "start", start_field = "startpos",
        end_field = "endpos"
    )
    data.table::setDT(pyclone.table)
    pyclone.table <- pyclone.table[, list(
        SampleID,
        mutation_id,
        ref_counts = ref_counts,
        alt_counts = alt_counts,
        normal_cn = 2L,
        major_cn = pmax(nMinor, nMajor),
        minor_cn = pmin(nMinor, nMajor),
        major_raw = nAraw, minor_raw = nBraw,
        region = region,
        Reference_Base = ref,
        Alternate_Base = var,
        fracA, nMaj_A, nMin_A,
        fracB, nMaj_B, nMin_B,
        fracC, nMaj_C, nMin_C,
        fracD, nMaj_D, nMin_D
    )]
    identify.subclonal.mut.copy.number.ascat <- function(x, sub.mat.mut, sub.mat.copy, region, sample, sex = "male") {
        mut <- sub.mat.mut[x, , drop = FALSE]
        ww <- which(as.numeric(sub.mat.copy$chr) == as.numeric(mut$chr) &
            as.numeric(sub.mat.copy$startpos) <= as.numeric(mut$start) &
            as.numeric(sub.mat.copy$endpos) >= as.numeric(mut$stop))
        copy <- sub.mat.copy[ww, , drop = FALSE]


        mutation_id <- paste(sample, mut$chr, mut$start, mut$ref, sep = ":")
        ref_counts <- mut[, "ref_counts"]
        alt_counts <- mut[, "alt_counts"]


        normal_cn <- 2
        region <- region
        Reference_Base <- mut$ref
        Alternate_Base <- mut$var

        if (nrow(copy) != 1) {
            minor_cn <- NA
            major_cn <- NA
            major_raw <- NA
            minor_raw <- NA
            fracA <- NA
            nMaj_A <- NA
            nMin_A <- NA
            fracB <- NA
            nMaj_B <- NA
            nMin_B <- NA
            fracC <- NA
            nMaj_C <- NA
            nMin_C <- NA
            fracD <- NA
            nMaj_D <- NA
            nMin_D <- NA

            output <- data.frame(mutation_id,
                ref_counts,
                alt_counts,
                normal_cn,
                minor_cn,
                major_cn,
                major_raw,
                minor_raw,
                region,
                Reference_Base,
                Alternate_Base,
                fracA,
                nMaj_A,
                nMin_A,
                fracB,
                nMaj_B,
                nMin_B,
                fracC,
                nMaj_C,
                nMin_C,
                fracD,
                nMaj_D,
                nMin_D,
                stringsAsFactors = FALSE
            )
            return(output)
        }

        minor_cn <- min(c(copy$nMajor, copy$nMinor))
        major_cn <- max(c(copy$nMajor, copy$nMinor))
        major_raw <- c(copy$nAraw)
        minor_raw <- c(copy$nBraw)
        fracA <- as.numeric(copy$fracA)
        nMaj_A <- as.numeric(copy$nMaj_A)
        nMin_A <- as.numeric(copy$nMin_A)
        fracB <- as.numeric(copy$fracB)
        nMaj_B <- as.numeric(copy$nMaj_B)
        nMin_B <- as.numeric(copy$nMin_B)
        fracC <- as.numeric(copy$fracC)
        nMaj_C <- as.numeric(copy$nMaj_C)
        nMin_C <- as.numeric(copy$nMin_C)
        fracD <- as.numeric(copy$fracD)
        nMaj_D <- as.numeric(copy$nMaj_D)
        nMin_D <- as.numeric(copy$nMin_D)

        output <- data.frame(mutation_id,
            ref_counts,
            alt_counts,
            normal_cn,
            minor_cn,
            major_cn,
            major_raw,
            minor_raw,
            region,
            Reference_Base,
            Alternate_Base,
            fracA,
            nMaj_A,
            nMin_A,
            fracB,
            nMaj_B,
            nMin_B,
            fracC,
            nMaj_C,
            nMin_C,
            fracD,
            nMaj_D,
            nMin_D,
            stringsAsFactors = FALSE
        )

        return(output)
    }
    pyclone.table2 <- data.frame(
        t(sapply(
            1:nrow(region.mut.table),
            identify.subclonal.mut.copy.number.ascat,
            as.data.frame(region.mut.table),
            as.data.frame(region.seg.phylo),
            region, sample
        )),
        stringsAsFactors = FALSE
    )
    data.table::setDT(pyclone.table2)
    pyclone.table2[, SampleID := region]
    data.table::setcolorder(pyclone.table2, names(pyclone.table))
    testthat::expect_true(all(vapply(
        pyclone.table2, function(x) all(lengths(x) == 1L), logical(1L)
    )))
    pyclone.table2[, names(pyclone.table2) := lapply(.SD, unlist),
        .SDcols = names(pyclone.table2)
    ]
    testthat::expect_equal(pyclone.table, pyclone.table2)

    # filter pyclone data -------------------------------------------
    pyclone.table <- pyclone.table[!is.na(pyclone.table$minor_cn)]
    pyclone.table <- pyclone.table[!is.na(pyclone.table$ref_counts)]
    pyclone.table <- pyclone.table[!duplicated(pyclone.table$mutation_id)]
    pyclone.table <- pyclone.table[ref_counts + alt_counts >= 1L]

    # let's load the purity estimate from VAF purity
    sample.purity <- region.seg.copy$ACF[1]

    # test run_subclonal_dissection works fine -----------------------------
    f.function <- function(c, purity, local.copy.number) {
        return(pmin(
            c((purity * c) / (2 * (1 - purity) + purity * local.copy.number)), 1
        ))
    }
    get.mut.mult <- function(CNt, Vaf, cellularity, CNn) {
        return((Vaf * 1 / cellularity) * ((cellularity * CNt) + CNn * (1 - cellularity)))
    }
    get.mutCopyNum <- function(i) {
        mutCopyNumber <- get.mut.mult(
            CNt = (as.numeric(pyClone.tsv$major_raw[i]) + as.numeric(pyClone.tsv$minor_raw[i])),
            Vaf = pyClone.tsv$obsVAF[i],
            cellularity = cellularity,
            CNn = unlist(pyClone.tsv$normal_cn[i])
        )
        return(mutCopyNumber)
    }
    get.absolute.ccf <- function(i) {
        # print(i)
        absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number) {
            f.function <- function(c, purity, local.copy.number) {
                return(min(c((purity * c) / (2 * (1 - purity) + purity * local.copy.number), 1)))
            }
            x <- dbinom(n.alt, depth, prob = sapply(seq(0.01, 1, length.out = 100), f.function, purity, local.copy.number))
            if (min(x) == 0) {
                x[length(x)] <- 1
            }

            names(x) <- seq(0.01, 1, length.out = 100)
            sub.cint <- function(x, prob = 0.95, n.alt, depth) {
                xnorm <- x / sum(x)
                xsort <- sort(xnorm, decreasing = TRUE)
                xcumLik <- cumsum(xsort)
                n <- sum(xcumLik < prob) + 1
                LikThresh <- xsort[n]
                cint <- x[xnorm >= LikThresh]
                all <- as.numeric(names(x))
                cellu <- as.numeric(names(cint))
                l.t <- cellu[1]
                r.t <- cellu[length(cellu)]
                m <- cellu[which.max(cint)]

                prob.subclonal <- sum(xnorm[1:90]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
                prob.clonal <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val

                data.frame(left = l.t, est = m, right = r.t, prob.subclonal = prob.subclonal, prob.clonal = prob.clonal)
            }
            return(sub.cint(x, n.alt = n.alt, depth = depth))
        }
        absolute.calc <- absolute.cancer.cell.fraction(n.alt = unlist(pyClone.tsv$alt_counts)[i], depth = (as.numeric(pyClone.tsv$ref_counts[i]) + as.numeric(pyClone.tsv$alt_counts[i])), purity = cellularity, local.copy.number = (as.numeric(pyClone.tsv$major_raw[i]) + as.numeric(pyClone.tsv$minor_raw[i])))
        absolute.ccf.0.05 <- absolute.calc[1]
        absolute.ccf.0.95 <- absolute.calc[3]
        absolute.ccf <- absolute.calc[2]
        prob.subclonal <- absolute.calc[4]
        prob.clonal <- absolute.calc[5]
        return(cbind(absolute.ccf, absolute.ccf.0.05, absolute.ccf.0.95, prob.subclonal, prob.clonal))
    }

    complete.mutation.table <- pyclone.table
    cellularity <- purity <- sample.purity
    order.by.pos <- TRUE
    pyClone.tsv <- complete.mutation.table[region == region]
    pyClone.tsv$expVAF <- calculate_vaf(
        1L, purity, pyClone.tsv$major_raw + pyClone.tsv$minor_raw
    )
    ## test calculate_vaf --------------------------------
    testthat::expect_equal(
        pyClone.tsv$expVAF,
        f.function(
            rep(1, nrow(pyClone.tsv)), purity,
            unlist(pyClone.tsv$major_raw) + unlist(pyClone.tsv$minor_raw)
        )
    )

    pyClone.tsv$obsVAF <- as.integer(pyClone.tsv$alt_counts) /
        (as.integer(pyClone.tsv$alt_counts) +
            as.integer(pyClone.tsv$ref_counts))

    ## test calculate_abs_ccf --------------------------------
    absolute_ccfs <- calculate_abs_ccf(
        pyClone.tsv$alt_counts, pyClone.tsv$ref_counts,
        CNts = pyClone.tsv$major_raw + pyClone.tsv$minor_raw,
        purity = purity
    )
    absolute_ccfs2 <- t(sapply(1:nrow(pyClone.tsv), get.absolute.ccf)) # nolint
    absolute_ccfs2 <- data.table::as.data.table(absolute_ccfs2)
    data.table::setnames(
        absolute_ccfs2,
        c("left", "right"), c("lower", "higher")
    )
    data.table::setcolorder(absolute_ccfs2, names(absolute_ccfs))
    testthat::expect_true(all(vapply(
        absolute_ccfs2, function(x) all(lengths(x) == 1L), logical(1L)
    )))
    absolute_ccfs2[, names(absolute_ccfs2) := lapply(.SD, unlist),
        .SDcols = names(absolute_ccfs2)
    ]
    testthat::expect_equal(absolute_ccfs, absolute_ccfs2)
    pyClone.tsv$absCCF <- absolute_ccfs$est
    pyClone.tsv$absCCF_lower <- absolute_ccfs$lower
    pyClone.tsv$absCCF_higher <- absolute_ccfs$higher

    ## test calculate_abs_ccf --------------------------------
    phylo.ccfs <- suppressWarnings(calculate_phylo_ccf(
        pyClone.tsv$alt_counts, pyClone.tsv$ref_counts,
        CNts = pyClone.tsv$major_raw + pyClone.tsv$minor_raw,
        purity, observed_vafs = pyClone.tsv$obsVAF,
        expected_vafs = pyClone.tsv$expVAF
    ))
    data.table::setDT(phylo.ccfs)
    pyClone.tsv$mutCopyNum <- phylo.ccfs$phyloCCF

    testthat::expect_equal(
        pyClone.tsv$mutCopyNum,
        sapply(seq_len(nrow(pyClone.tsv)), get.mutCopyNum)
    )
    pyClone.tsv <- cbind(pyClone.tsv, phylo.ccfs)
    testthat::expect_equal(
        pyClone.tsv$phyloCCF_lower,
        sapply(seq_len(nrow(pyClone.tsv)), function(var.count, depth, e, CNtumor, cellularity, CNn, i) {
            VAF <- suppressWarnings(stats::prop.test(
                var.count[i], depth[i], e[i]
            )$conf.int[1])
            return(get.mut.mult(
                CNt = CNtumor[i],
                Vaf = VAF, cellularity = cellularity, CNn = CNn
            ))
        },
        var.count = as.numeric(pyClone.tsv$alt_counts),
        depth = as.numeric(pyClone.tsv$alt_counts) + as.numeric(pyClone.tsv$ref_counts),
        e = as.numeric(pyClone.tsv$expVAF),
        CNtumor = as.numeric(pyClone.tsv$major_raw) + as.numeric(pyClone.tsv$minor_raw),
        cellularity = cellularity,
        CNn = 2
        )
    )
    testthat::expect_equal(
        pyClone.tsv$phyloCCF_higher,
        sapply(seq_len(nrow(pyClone.tsv)),
            function(var.count, depth, e, CNtumor, cellularity, CNn, i) {
                VAF <- suppressWarnings(prop.test(
                    var.count[i],
                    depth[i],
                    e[i]
                )$conf.int[2])
                return(get.mut.mult(CNt = CNtumor[i], Vaf = VAF, cellularity = cellularity, CNn = CNn))
            },
            var.count = as.numeric(pyClone.tsv$alt_counts),
            depth = as.numeric(pyClone.tsv$alt_counts) + as.numeric(pyClone.tsv$ref_counts),
            e = as.numeric(pyClone.tsv$expVAF),
            CNtumor = as.numeric(pyClone.tsv$major_raw) + as.numeric(pyClone.tsv$minor_raw),
            cellularity = cellularity,
            CNn = 2
        )
    )

    # test
    # pre-allocate column for downstream analysis
    pyClone.tsv$no.chrs.bearing.mut <- 1
    pyClone.tsv$whichFrac <- NA_character_
    # pyClone.tsv$CPNChange <- 0L
    pyClone.tsv2 <- as.data.frame(pyClone.tsv)

    # Identification of subclonal mutations  -----------------------
    pyClone.tsv[, is_subclone := mutCopyNum > 0.01 & # nolint
        # pyClone.tsv$absCCF_higher < 1L
        absolute_ccfs$prob.subclonal > 0.5]
    pyClone.tsv[is_subclone & fracA == 1L, whichFrac := "A,B"] # nolint

    # define best CN -----------------------------------------------
    pyClone.tsv[
        is_subclone & fracA != 1L & is.na(fracC), # nolint
        best_cn := calculate_best_cn_for_loss_mut( # nolint
            nMaj_A, nMaj_B, nMin_A, nMin_B, # nolint
            fracA, fracB, fracA, fracB, mutCopyNum # nolint
        )
    ]
    pyClone.tsv[
        is_subclone & fracA != 1L & !is.na(fracC), # nolint
        best_cn := calculate_best_cn_for_loss_mut( # nolint
            nMaj_A + nMaj_B, nMaj_C + nMaj_D, # nolint
            nMin_A + nMin_D, nMin_C + nMin_B, # nolint
            fracA + fracB, fracC + fracD, # nolint
            fracA + fracD, fracC + fracB, mutCopyNum # nolint
        )
    ]

    pyClone.tsv[
        best_cn == 1L, # nolint
        whichFrac := data.table::fcase( # nolint
            is.na(fracC), "A,B", !is.na(fracC), "A,B,C,D" # nolint
        )
    ] # nolint

    # check whether subclonal CN results in clonal mutation
    # otherwise subclonal CN doesn't explain subclonal MCN
    pyClone.tsv[
        best_cn != 1L, # nolint
        explained_by_cn_pvalue := prop_test_pvalues( # nolint
            alt_counts, (alt_counts + ref_counts) * purity, # nolint
            prop = expVAF * best_cn / purity, # nolint
            alternative = "less",
            correction = NULL
        )
    ]

    pyClone.tsv[
        explained_by_cn_pvalue > 0.01, # nolint
        c("phyloCCF", "phyloCCF_lower", "phyloCCF_higher", "no.chrs.bearing.mut", "expVAF") := {
            tmp_expected_prop <- expVAF * best_cn # nolint
            tmp_vafs <- prop_test_ci(
                alt_counts, alt_counts + ref_counts, # nolint
                tmp_expected_prop
            ) # nolint
            tmp_CNts <- major_raw + minor_raw # nolint
            list(
                mutCopyNum / best_cn, # nolint
                calculate_obs_mut(
                    CNts = tmp_CNts,
                    vafs = tmp_vafs[[1L]],
                    purity = purity
                ) / best_cn,
                calculate_obs_mut(
                    CNts = tmp_CNts,
                    vafs = tmp_vafs[[2]],
                    purity = purity,
                ) / best_cn,
                best_cn, tmp_expected_prop
            )
        }
    ]
    pyClone.tsv[, explained_by_cn_pvalue := NULL] # nolint
    pyClone.tsv[, is_subclone := NULL] # nolint
    pyClone.tsv[, best_cn := NULL] # nolint

    # official function
    subclonal.mutations <- which(absolute_ccfs$prob.subclonal > 0.5 &
        pyClone.tsv2$mutCopyNum > 0.01)

    if (length(subclonal.mutations) > 0) {
        # assume subclonal muts are on one chromosome copy
        # therefore mutation copy number must be subclonal fraction
        # of the higher CN subclone (i.e. lost in lower CN subclone)
        # or 1 (i.e. present in both subclones)
        for (a in subclonal.mutations)
        {
            # print(a)
            mut.info <- pyClone.tsv2[a, , drop = FALSE]
            # check whether the fracA is hundred%
            if (mut.info$fracA == 1) {
                pyClone.tsv2$whichFrac[a] <- c("A,B")
                next
            }

            # let's see what happens when we only have A and B
            if (is.na(mut.info$fracC)) {
                # determine which fraction has a copy number loss
                possible.subclonal.fractions <- c(1)
                if (unlist(mut.info$nMaj_A) > unlist(mut.info$nMaj_B)) {
                    possible.subclonal.fractions <- c(possible.subclonal.fractions, unlist(mut.info$fracA))
                }
                if (unlist(mut.info$nMaj_B) > unlist(mut.info$nMaj_A)) {
                    possible.subclonal.fractions <- c(possible.subclonal.fractions, unlist(mut.info$fracB))
                }
                if (unlist(mut.info$nMin_B) > unlist(mut.info$nMin_A)) {
                    possible.subclonal.fractions <- c(possible.subclonal.fractions, unlist(mut.info$fracB))
                }
                if (unlist(mut.info$nMin_A) > unlist(mut.info$nMin_B)) {
                    possible.subclonal.fractions <- c(possible.subclonal.fractions, unlist(mut.info$fracA))
                }


                best.CN <- possible.subclonal.fractions[which.min(abs(mut.info$mutCopyNum / possible.subclonal.fractions - 1))]
                if (best.CN == 1) {
                    pyClone.tsv2$whichFrac[a] <- c("A,B")
                    next
                }
                var.count <- as.numeric(mut.info$alt_counts)
                depth.count <- as.numeric(mut.info$alt_counts) + as.numeric(mut.info$ref_counts)
                expected.prop <- pyClone.tsv2$expVAF[a] * best.CN

                # check whether subclonal CN results in clonal mutation
                # otherwise subclonal CN doesn't explain subclonal MCN
                if (best.CN != 1 & prop.test(var.count, depth.count * purity, expected.prop / purity, alternative = "less")$p.value > 0.01) {
                    pyClone.tsv2$phyloCCF[a] <- mut.info$mutCopyNum / best.CN
                    pyClone.tsv2$phyloCCF_lower[a] <- get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[1], cellularity = cellularity, CNn = 2) / best.CN
                    pyClone.tsv2$phyloCCF_higher[a] <- get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[2], cellularity = cellularity, CNn = 2) / best.CN
                    pyClone.tsv2$no.chrs.bearing.mut[a] <- best.CN
                    # pyClone.tsv2$whichFrac[a] = c('A')
                    pyClone.tsv2$expVAF[a] <- expected.prop
                }
            }
            if (!is.na(mut.info$fracC)) {
                possible.subclonal.fractions <- c(1)
                # there are going to be quite a few options here, let's list them
                if ((unlist(mut.info$nMaj_A) + unlist(mut.info$nMaj_B)) > (unlist(mut.info$nMaj_C) + unlist(mut.info$nMaj_D))) {
                    fracAB <- unlist(mut.info$fracA) + unlist(mut.info$fracB)
                    possible.subclonal.fractions <- c(possible.subclonal.fractions, fracAB)
                }
                if ((unlist(mut.info$nMaj_A) + unlist(mut.info$nMaj_B)) < (unlist(mut.info$nMaj_C) + unlist(mut.info$nMaj_D))) {
                    fracCD <- unlist(mut.info$fracC) + unlist(mut.info$fracD)
                    possible.subclonal.fractions <- c(possible.subclonal.fractions, fracCD)
                }
                if ((unlist(mut.info$nMin_A) + unlist(mut.info$nMin_D)) > (unlist(mut.info$nMin_C) + unlist(mut.info$nMin_B))) {
                    fracAD <- unlist(mut.info$fracA) + unlist(mut.info$fracD)
                    possible.subclonal.fractions <- c(possible.subclonal.fractions, fracAD)
                }
                if ((unlist(mut.info$nMin_A) + unlist(mut.info$nMin_D)) < (unlist(mut.info$nMin_C) + unlist(mut.info$nMin_B))) {
                    fracBC <- unlist(mut.info$fracB) + unlist(mut.info$fracC)
                    possible.subclonal.fractions <- c(possible.subclonal.fractions, fracBC)
                }

                best.CN <- possible.subclonal.fractions[
                    which.min(abs(mut.info$mutCopyNum / possible.subclonal.fractions - 1))
                ]
                if (best.CN == 1) {
                    pyClone.tsv2$whichFrac[a] <- c("A,B,C,D")
                    next
                }
                var.count <- as.numeric(mut.info$alt_counts)
                depth.count <- as.numeric(mut.info$alt_counts) + as.numeric(mut.info$ref_counts)
                expected.prop <- pyClone.tsv2$expVAF[a] * best.CN

                # check whether subclonal CN results in clonal mutation
                # otherwise subclonal CN doesn't explain subclonal MCN
                if (best.CN != 1 & prop.test(var.count, depth.count * purity, expected.prop / purity, alternative = "less")$p.value > 0.01) {
                    pyClone.tsv2$phyloCCF[a] <- mut.info$mutCopyNum / best.CN
                    pyClone.tsv2$phyloCCF_lower[a] <- get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[1], cellularity = cellularity, CNn = 2) / best.CN
                    pyClone.tsv2$phyloCCF_higher[a] <- get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[2], cellularity = cellularity, CNn = 2) / best.CN
                    pyClone.tsv2$no.chrs.bearing.mut[a] <- best.CN
                    # pyClone.tsv2$min.or.maj[a] = min.or.maj
                    pyClone.tsv2$expVAF[a] <- expected.prop
                }
            }
        }
    }
    data.table::setDT(pyClone.tsv2)
    testthat::expect_equal(pyClone.tsv, pyClone.tsv2)

    # Next, let's deal with potentially amplified mutations
    # convert MCN to subclonal fraction - tricky for amplified mutations test
    # for mutations in more than 1 copy
    pyClone.tsv[
        ,
        amp_mut_pvalue := suppressWarnings(prop_test_pvalues( # nolint
            alt_counts, alt_counts + ref_counts, expVAF, # nolint
            alternative = "greater", # nolint
            correction = NULL
        ))
    ]
    # official function
    p.vals <- sapply(seq_len(nrow(pyClone.tsv2)),
        function(var.count, depth, e, i) {
            suppressWarnings(prop.test(var.count[i],
                depth[i],
                e[i],
                alternative = "greater"
            )$p.value)
        },
        var.count = as.numeric(pyClone.tsv2$alt_counts),
        depth = pyClone.tsv2$alt_counts + pyClone.tsv2$ref_counts,
        e = as.numeric(pyClone.tsv2$expVAF)
    )
    data.table::setDF(pyClone.tsv2)
    testthat::expect_equal(pyClone.tsv$amp_mut_pvalue, p.vals)

    # copy numbers of subclones can only differ by 1 or 0 (as assumed when
    # calling subclones)
    # nolint start
    pyClone.tsv[amp_mut_pvalue < 0.05 & mutCopyNum > 1L & is.na(fracC), c("no.chrs.bearing.mut", "best_cn") := {
        suppressWarnings(calculate_best_cn_for_amp_mut(
            nMaj_A, nMaj_B,
            fracA = fracA, fracB = fracB,
            alt_counts = alt_counts, mutCopyNum = mutCopyNum,
            subclone_correction = FALSE
        ))
    }]
    pyClone.tsv[amp_mut_pvalue < 0.05 & mutCopyNum > 1L & !is.na(fracC), c("no.chrs.bearing.mut", "best_cn") := {
        suppressWarnings(calculate_best_cn_for_amp_mut(
            nMaj_A, nMaj_C, nMin_A, nMin_B,
            fracA + fracB, fracC + fracD,
            fracA + fracD, fracC + fracB,
            alt_counts = alt_counts, mutCopyNum = mutCopyNum,
            subclone_correction = FALSE
        ))
    }]
    pyClone.tsv[amp_mut_pvalue < 0.05 & mutCopyNum > 1L, c("phyloCCF", "phyloCCF_lower", "phyloCCF_higher", "expVAF") := {
        list(
            mutCopyNum / best_cn,
            phyloCCF_lower / best_cn,
            phyloCCF_higher / best_cn,
            expVAF * best_cn
        )
    }]
    # nolint end

    # official function
    amplified.muts <- which(p.vals <= 0.05 & pyClone.tsv2$mutCopyNum > 1)
    # copy numbers of subclones can only differ by 1 or 0 (as assumed when calling subclones)
    if (length(amplified.muts) > 0) {
        for (a in amplified.muts)
        {
            mut.info <- pyClone.tsv2[a, , drop = FALSE]
            if (is.na(mut.info$fracC)) {
                max.CN2 <- 0
                max.CN1 <- unlist(mut.info$nMaj_A)
                frac1.mut <- unlist(mut.info$fracA)
                frac2.mut <- 0


                # swap subclones, so that the one with the higher CN is first
                if (unlist(mut.info$nMaj_B) > max.CN1) {
                    max.CN2 <- max.CN1
                    max.CN1 <- unlist(mut.info$nMaj_B)
                    frac2.mut <- frac1.mut
                    frac1.mut <- unlist(mut.info$fracB)
                } else {
                    max.CN2 <- unlist(mut.info$nMaj_B)
                    frac2.mut <- unlist(mut.info$fracB)
                }


                best.err <- mut.info$mutCopyNum - 1
                best.CN <- 1
                for (j in 1:max.CN1) {
                    for (k in (j - 1):min(j, max.CN2)) {
                        potential.CN <- j * frac1.mut + k * frac2.mut
                        err <- abs(mut.info$mutCopyNum / potential.CN - 1)
                        if (err < best.err) {
                            pyClone.tsv2$no.chrs.bearing.mut[a] <- potential.CN
                            best.err <- err
                            best.CN <- potential.CN
                        }
                    }
                }
                pyClone.tsv2$phyloCCF[a] <- pyClone.tsv2$mutCopyNum[a] / best.CN
                pyClone.tsv2$expVAF[a] <- pyClone.tsv2$expVAF[a] * best.CN
                pyClone.tsv2$phyloCCF_lower[a] <- pyClone.tsv2$phyloCCF_lower[a] / best.CN
                pyClone.tsv2$phyloCCF_higher[a] <- pyClone.tsv2$phyloCCF_higher[a] / best.CN
            }
            if (!is.na(mut.info$fracC)) {
                max.CN2 <- 0
                max.CN1 <- unlist(mut.info$nMaj_A)
                frac1.mut <- unlist(mut.info$fracA) + unlist(mut.info$fracB)
                frac2.mut <- 0


                # swap subclones, so that the one with the higher CN is first
                if (mut.info$nMaj_C > max.CN1) {
                    max.CN2 <- max.CN1
                    max.CN1 <- unlist(mut.info$nMaj_C)
                    frac2.mut <- frac1.mut
                    frac1.mut <- unlist(mut.info$fracC) + unlist(mut.info$fracD)
                } else {
                    max.CN2 <- unlist(mut.info$nMaj_C)
                    frac2.mut <- unlist(mut.info$fracC) + unlist(mut.info$fracD)
                }


                best.err <- mut.info$mutCopyNum - 1
                best.CN <- 1
                for (j in 1:max.CN1) {
                    for (k in (j - 1):min(j, max.CN2)) {
                        potential.CN <- j * frac1.mut + k * frac2.mut
                        err <- abs(mut.info$mutCopyNum / potential.CN - 1)


                        if (err < best.err) {
                            pyClone.tsv2$no.chrs.bearing.mut[a] <- potential.CN
                            best.err <- err
                            best.CN <- potential.CN
                        }
                    }
                }

                # next, let's see whether we should also look at the minor allele
                max.CN2 <- 0
                max.CN1 <- unlist(mut.info$nMin_A)
                frac1.mut <- unlist(mut.info$fracA) + unlist(mut.info$fracD)
                frac2.mut <- 0


                # swap subclones, so that the one with the higher CN is first
                if (unlist(mut.info$nMin_B) > max.CN1) {
                    max.CN2 <- max.CN1
                    max.CN1 <- unlist(mut.info$nMin_B)
                    frac2.mut <- frac1.mut
                    frac1.mut <- unlist(mut.info$fracB) + unlist(mut.info$fracC)
                } else {
                    max.CN2 <- unlist(mut.info$nMin_B)
                    frac2.mut <- unlist(mut.info$fracB) + unlist(mut.info$fracC)
                }


                for (j in 1:max.CN1) {
                    for (k in (j - 1):min(j, max.CN2)) {
                        potential.CN <- j * frac1.mut + k * frac2.mut
                        err <- abs(mut.info$mutCopyNum / potential.CN - 1)


                        if (err < best.err) {
                            pyClone.tsv2$no.chrs.bearing.mut[a] <- potential.CN
                            best.err <- err
                            best.CN <- potential.CN
                        }
                    }
                }
                pyClone.tsv2$phyloCCF[a] <- pyClone.tsv2$mutCopyNum[a] / best.CN
                pyClone.tsv2$expVAF[a] <- pyClone.tsv2$expVAF[a] * best.CN
                pyClone.tsv2$phyloCCF_lower[a] <- pyClone.tsv2$phyloCCF_lower[a] / best.CN
                pyClone.tsv2$phyloCCF_higher[a] <- pyClone.tsv2$phyloCCF_higher[a] / best.CN
            }
        }
    }
    data.table::setDT(pyClone.tsv2)
    testthat::expect_equal(
        pyClone.tsv[, !c("best_cn", "amp_mut_pvalue")],
        pyClone.tsv2
    )
    # testthat::expect_equal(
    #     pyClone.tsv[
    #         amp_mut_pvalue <= 0.05 & mutCopyNum > 1L & !is.na(fracC),
    #         !c("best_cn", "amp_mut_pvalue")
    #     ],
    #     pyClone.tsv2[
    #         p.vals <= 0.05   & mutCopyNum > 1L & !is.na(fracC)
    #     ]
    # )
})
