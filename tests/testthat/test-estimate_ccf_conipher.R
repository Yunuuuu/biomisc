test_that("estimate_ccf function works well", {
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

    # define subclone copy number
    seg.mat.phylo <- define_subclone_cn(
        seg.mat.copy,
        min_subclonal = 0.01
    )

    # let's run PyClone, but correct for copy number before running
    # prepare pyclone table --------------------------
    sample <- region <- seg.mat.phylo$SampleID[1]
    region.mut.table <- mut.table
    region.seg.copy <- seg.mat.copy[SampleID %in% region]
    region.seg.phylo <- seg.mat.phylo[SampleID %in% region]
    # just extract the segmented CNV for this sample
    # region.seg.copy <- seg.mat.phylo[SampleID %in% region]
    pyclone.table <- suppressWarnings(mut_match_cn(region.mut.table, seg.mat.phylo,
        on_sample = "SampleID", on_chr = "chr",
        mut_pos = "start", start_field = "startpos",
        end_field = "endpos"
    ))
    pyclone.table <- pyclone.table[, list(
        SampleID,
        gender = "male", mutation_id, chr,
        pos = startpos,
        ref_counts = ref_counts,
        alt_counts = alt_counts,
        normal_cn = 2L,
        major_cn = pmax(nMinor, nMajor),
        minor_cn = pmin(nMinor, nMajor),
        major_raw = nAraw, minor_raw = nBraw,
        purity = ACF,
        ref = ref, alt = var, fracA, nMaj_A, nMin_A,
        fracB, nMaj_B, nMin_B,
        fracC, nMaj_C, nMin_C,
        fracD, nMaj_D, nMin_D
    )]

    # filter pyclone data -------------------------------------------
    pyclone.table <- pyclone.table[!is.na(pyclone.table$minor_cn)]
    pyclone.table <- pyclone.table[!is.na(pyclone.table$ref_counts)]
    pyclone.table <- pyclone.table[!duplicated(pyclone.table$mutation_id)]
    pyclone.table <- pyclone.table[ref_counts + alt_counts >= 1L]

    # let's load the purity estimate from VAF purity
    sample.purity <- region.seg.copy$ACF[1]

    pyclone.table[, normal_cn := define_normal_cn(gender, chr)]
    out <- suppressWarnings(estimate_ccf(pyclone.table,
        conipher = TRUE, min_subclonal = 0.05,
        min_vaf_to_explain = 0.05
    ))
    calculate_phylo_ccf_withBH <- function(region,
                                           complete.mutation.table,
                                           purity,
                                           min.subclonal = 0.05,
                                           minVAFtoEXPLAIN = 0.05) {
        # complete.mutation.table <- pyclone.table
        # min.subclonal <- 0.05
        # minVAFtoEXPLAIN <- 0.05
        # purity <- sample.purity
        pyClone.tsv <- data.frame(complete.mutation.table, stringsAsFactors = FALSE)
        row.names <- pyClone.tsv$mutation_id
        cellularity <- as.numeric(purity)
        f.function <- function(c, purity, local.copy.number, normal.copy.number) {
            return(pmin(1, c((purity * c) / (normal.copy.number * (1 - purity) + purity * local.copy.number))))
        }

        get.mut.mult <- function(CNt, Vaf, cellularity, CNn) {
            return((Vaf * 1 / cellularity) * ((cellularity * CNt) + CNn * (1 - cellularity)))
        }

        pyClone.tsv$expected.VAF <- f.function(rep(1, nrow(pyClone.tsv)), cellularity, unlist(pyClone.tsv$major_raw) + unlist(pyClone.tsv$minor_raw), normal.copy.number = unlist(pyClone.tsv$normal_cn))

        # let's put it back together
        rownames(pyClone.tsv) <- pyClone.tsv$mutation_id
        pyClone.tsv$obs.VAF <- as.numeric(pyClone.tsv$alt_counts) / (as.numeric(pyClone.tsv$alt_counts) + as.numeric(pyClone.tsv$ref_counts))

        get.mutCopyNum <- function(i) {
            # let's simplify things, and assume there are max two subclonal populations
            mutCopyNumber <- get.mut.mult(CNt = (as.numeric(pyClone.tsv$major_raw[i]) + as.numeric(pyClone.tsv$minor_raw[i])), Vaf = pyClone.tsv$obs.VAF[i], cellularity = cellularity, CNn = unlist(pyClone.tsv$normal_cn[i]))
            return(mutCopyNumber)
        }
        get.absolute.ccf <- function(i) {
            # print(i)
            absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number, normal.copy.number) {
                f.function <- function(c, purity, local.copy.number, normal.copy.number) {
                    return(min(1, c((purity * c) / (normal.copy.number * (1 - purity) + purity * local.copy.number))))
                }
                x <- dbinom(n.alt, depth, prob = sapply(seq(0.01, 1, length.out = 100), f.function, purity, local.copy.number, normal.copy.number))
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

            absolute.calc <- absolute.cancer.cell.fraction(
                n.alt = unlist(pyClone.tsv$alt_counts)[i],
                depth = (as.numeric(pyClone.tsv$ref_counts[i]) + as.numeric(pyClone.tsv$alt_counts[i])),
                purity = cellularity,
                local.copy.number = (as.numeric(pyClone.tsv$major_raw[i]) + as.numeric(pyClone.tsv$minor_raw[i])),
                normal.copy.number = as.numeric(unlist(pyClone.tsv$normal_cn))[i]
            )
            absolute.ccf.0.05 <- absolute.calc[1]
            absolute.ccf.0.95 <- absolute.calc[3]
            absolute.ccf <- absolute.calc[2]
            prob.subclonal <- absolute.calc[4]
            prob.clonal <- absolute.calc[5]
            return(cbind(absolute.ccf, absolute.ccf.0.05, absolute.ccf.0.95, absolute.ccf, prob.subclonal, prob.clonal))
        }

        pyClone.tsv$mutCopyNum <- sapply(1:nrow(pyClone.tsv), get.mutCopyNum)
        absolute.ccfs <- t(sapply(1:nrow(pyClone.tsv), get.absolute.ccf))
        absolute.ccfs <- data.frame(absolute.ccfs, stringsAsFactors = FALSE)

        pyClone.tsv$absolute.ccf <- as.numeric(absolute.ccfs$est)
        pyClone.tsv$absolute.ccf.0.05 <- as.numeric(absolute.ccfs$left)
        pyClone.tsv$absolute.ccf.0.95 <- as.numeric(absolute.ccfs$right)
        pyClone.tsv$phyloCCF <- as.numeric(pyClone.tsv$mutCopyNum)
        pyClone.tsv$phyloCCF.0.05 <- sapply(1:nrow(pyClone.tsv), function(var.count, depth, e, CNtumor, cellularity, CNn, i) {
            VAF <- prop.test(
                var.count[i],
                depth[i],
                e[i]
            )$conf.int[1]
            return(get.mut.mult(CNt = CNtumor[i], Vaf = VAF, cellularity = cellularity, CNn = CNn[i]))
        },
        var.count = as.numeric(pyClone.tsv$alt_counts),
        depth = as.numeric(pyClone.tsv$alt_counts) + as.numeric(pyClone.tsv$ref_counts),
        e = as.numeric(pyClone.tsv$expected.VAF),
        CNtumor = as.numeric(pyClone.tsv$major_raw) + as.numeric(pyClone.tsv$minor_raw),
        cellularity = cellularity,
        CNn = as.numeric(unlist(pyClone.tsv$normal_cn))
        )
        pyClone.tsv$phyloCCF.0.95 <- sapply(1:nrow(pyClone.tsv), function(var.count, depth, e, CNtumor, cellularity, CNn, i) {
            VAF <- prop.test(
                var.count[i],
                depth[i],
                e[i]
            )$conf.int[2]
            return(get.mut.mult(CNt = CNtumor[i], Vaf = VAF, cellularity = cellularity, CNn = CNn[i]))
        },
        var.count = as.numeric(pyClone.tsv$alt_counts),
        depth = as.numeric(pyClone.tsv$alt_counts) + as.numeric(pyClone.tsv$ref_counts),
        e = as.numeric(pyClone.tsv$expected.VAF),
        CNtumor = as.numeric(pyClone.tsv$major_raw) + as.numeric(pyClone.tsv$minor_raw),
        cellularity = cellularity,
        CNn = as.numeric(unlist(pyClone.tsv$normal_cn))
        )


        pyClone.tsv$no.chrs.bearing.mut <- 1

        # let's check which subclonal mutations can be explained by copy number
        # pyClone.tsv.subclonal         <- pyClone.tsv[which(pyClone.tsv$absolute.ccf.0.95<1&pyClone.tsv$mutCopyNum>0.01),]
        subclonal.mutations <- which(pyClone.tsv$absolute.ccf.0.95 < 1 & pyClone.tsv$mutCopyNum > 0.01)
        subclonal.mutations <- which(absolute.ccfs$prob.subclonal > 0.5 & pyClone.tsv$mutCopyNum > 0.01 & pyClone.tsv$obs.VAF >= minVAFtoEXPLAIN)
        pyClone.tsv$whichFrac <- NA
        pyClone.tsv$CPNChange <- 0


        p.vals <- sapply(1:nrow(pyClone.tsv), function(var.count, depth, e, i) {
            prop.test(var.count[i],
                depth[i],
                e[i],
                alternative = "less"
            )$p.value
        },
        var.count = as.numeric(pyClone.tsv$alt_counts),
        depth = as.numeric(pyClone.tsv$alt_counts) + as.numeric(pyClone.tsv$ref_counts),
        e = as.numeric(pyClone.tsv$expected.VAF)
        )

        p.vals.adj <- p.adjust(p.vals, method = "BH")
        subclonal.mutations <- which(pyClone.tsv$absolute.ccf.0.95 < 1 & pyClone.tsv$mutCopyNum > 0.01 & p.vals.adj < 0.01)

        if (length(subclonal.mutations) > 0) {
            # assume subclonal muts are on one chromosome copy
            # therefore mutation copy number must be subclonal fraction
            # of the higher CN subclone (i.e. lost in lower CN subclone)
            # or 1 (i.e. present in both subclones)
            for (a in subclonal.mutations)
            {
                # print(a)
                mut.info <- pyClone.tsv[a, , drop = FALSE]
                # check whether the fracA is hundred%
                if (mut.info$fracA == 1) {
                    pyClone.tsv$whichFrac[a] <- c("A,B")
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
                    if (best.CN == 1 | best.CN >= (1 - min.subclonal)) {
                        pyClone.tsv$whichFrac[a] <- c("A,B")
                        next
                    }
                    var.count <- as.numeric(mut.info$alt_counts)
                    depth.count <- as.numeric(mut.info$alt_counts) + as.numeric(mut.info$ref_counts)
                    expected.prop <- pyClone.tsv$expected.VAF[a] * best.CN

                    # check whether subclonal CN results in clonal mutation
                    # otherwise subclonal CN doesn't explain subclonal MCN
                    if (best.CN != 1 & prop.test(var.count, depth.count * purity, expected.prop / purity, alternative = "less")$p.value > 0.01) {
                        if ((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[1], cellularity = cellularity, CNn = mut.info$normal_cn) / best.CN) > 1) {
                            next
                        }

                        if ((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[1], cellularity = cellularity, CNn = mut.info$normal_cn) / best.CN) <= 1) {
                            pyClone.tsv$phyloCCF[a] <- mut.info$mutCopyNum / best.CN
                            pyClone.tsv$phyloCCF.0.05[a] <- get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[1], cellularity = cellularity, CNn = mut.info$normal_cn) / best.CN
                            pyClone.tsv$phyloCCF.0.95[a] <- get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[2], cellularity = cellularity, CNn = mut.info$normal_cn) / best.CN
                            pyClone.tsv$no.chrs.bearing.mut[a] <- best.CN
                            # pyClone.tsv$whichFrac[a] = c('A')
                            pyClone.tsv$expected.VAF[a] <- expected.prop
                            pyClone.tsv$CPNChange[a] <- 1
                        }
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

                    best.CN <- possible.subclonal.fractions[which.min(abs(mut.info$mutCopyNum / possible.subclonal.fractions - 1))]
                    if (best.CN == 1 | best.CN >= (1 - min.subclonal)) {
                        pyClone.tsv$whichFrac[a] <- c("A,B,C,D")
                        next
                    }
                    var.count <- as.numeric(mut.info$alt_counts)
                    depth.count <- as.numeric(mut.info$alt_counts) + as.numeric(mut.info$ref_counts)
                    expected.prop <- pyClone.tsv$expected.VAF[a] * best.CN

                    # check whether subclonal CN results in clonal mutation
                    # otherwise subclonal CN doesn't explain subclonal MCN
                    if (best.CN != 1 & prop.test(var.count, depth.count * purity, expected.prop / purity, alternative = "less")$p.value > 0.01) {
                        if ((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[1], cellularity = cellularity, CNn = mut.info$normal_cn) / best.CN) > 1) {
                            next
                        }

                        if ((get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[1], cellularity = cellularity, CNn = mut.info$normal_cn) / best.CN) <= 1) {
                            pyClone.tsv$phyloCCF[a] <- mut.info$mutCopyNum / best.CN
                            pyClone.tsv$phyloCCF.0.05[a] <- get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[1], cellularity = cellularity, CNn = mut.info$normal_cn) / best.CN
                            pyClone.tsv$phyloCCF.0.95[a] <- get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw)) + as.numeric(unlist(mut.info$minor_raw)), Vaf = prop.test(var.count, depth.count, expected.prop)$conf.int[2], cellularity = cellularity, CNn = mut.info$normal_cn) / best.CN
                            pyClone.tsv$no.chrs.bearing.mut[a] <- best.CN
                            # pyClone.tsv$whichFrac[a] = c('A')
                            pyClone.tsv$expected.VAF[a] <- expected.prop
                            pyClone.tsv$CPNChange[a] <- 1
                        }
                    }
                }
            }
        }
        # Next, let's deal with potentially amplified mutations
        # convert MCN to subclonal fraction - tricky for amplified mutations
        # test for mutations in more than 1 copy



        p.vals <- sapply(1:nrow(pyClone.tsv), function(var.count, depth, e, i) {
            prop.test(var.count[i],
                depth[i],
                e[i],
                alternative = "greater"
            )$p.value
        },
        var.count = as.numeric(pyClone.tsv$alt_counts),
        depth = as.numeric(pyClone.tsv$alt_counts) + as.numeric(pyClone.tsv$ref_counts),
        e = as.numeric(pyClone.tsv$expected.VAF)
        )

        amplified.muts <- which(p.vals <= 0.05 & pyClone.tsv$mutCopyNum > 1)


        # copy numbers of subclones can only differ by 1 or 0 (as assumed when calling subclones)
        if (length(amplified.muts) > 0) {
            for (a in amplified.muts)
            {
                # print(a)
                mut.info <- pyClone.tsv[a, , drop = FALSE]
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
                    allCNs <- c(1)
                    for (j in 1:max.CN1) {
                        for (k in (j - 1):min(j, max.CN2)) {
                            potential.CN <- j * frac1.mut + k * frac2.mut
                            err <- abs(mut.info$mutCopyNum / potential.CN - 1)
                            # UpperConf = mut.info$
                            if (err < best.err) {
                                # print(potential.CN)
                                pyClone.tsv$no.chrs.bearing.mut[a] <- potential.CN
                                best.err <- err
                                best.CN <- potential.CN
                                allCNs <- c(allCNs, best.CN)
                            }
                        }
                    }




                    # let's just make sure we haven't created a subclonal mutation
                    if (as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN) < 1) {
                        if (prop.test(
                            as.numeric(pyClone.tsv$alt_counts[a]) / 2,
                            round(as.numeric(pyClone.tsv$alt_counts[a]) / as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN)),
                            0.5
                        )$p.value < 0.05 & best.CN > 1) {
                            best.CN <- max(allCNs[!allCNs %in% best.CN])
                        }
                    }


                    pyClone.tsv$phyloCCF[a] <- pyClone.tsv$mutCopyNum[a] / best.CN
                    pyClone.tsv$expected.VAF[a] <- pyClone.tsv$expected.VAF[a] * best.CN
                    pyClone.tsv$phyloCCF.0.05[a] <- pyClone.tsv$phyloCCF.0.05[a] / best.CN
                    pyClone.tsv$phyloCCF.0.95[a] <- pyClone.tsv$phyloCCF.0.95[a] / best.CN
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
                    allCNs <- c(best.CN)
                    for (j in 1:max.CN1) {
                        for (k in (j - 1):min(j, max.CN2)) {
                            potential.CN <- j * frac1.mut + k * frac2.mut
                            err <- abs(mut.info$mutCopyNu / potential.CN - 1)


                            if (err < best.err) {
                                pyClone.tsv$no.chrs.bearing.mut[a] <- potential.CN
                                best.err <- err
                                best.CN <- potential.CN
                                allCNs <- c(allCNs, best.CN)
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


                    # best.err = mut.info$mutCopyNum - 1
                    # best.CN=1

                    for (j in 1:max.CN1) {
                        for (k in (j - 1):min(j, max.CN2)) {
                            potential.CN <- j * frac1.mut + k * frac2.mut
                            err <- abs(mut.info$mutCopyNu / potential.CN - 1)


                            if (err < best.err) {
                                pyClone.tsv$no.chrs.bearing.mut[a] <- potential.CN
                                best.err <- err
                                best.CN <- potential.CN
                                allCNs <- c(allCNs, best.CN)
                            }
                        }
                    }
                    # let's just make sure we haven't created a subclonal mutation
                    if (as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN) < 1) {
                        if (prop.test(
                            as.numeric(pyClone.tsv$alt_counts[a]) / 2,
                            round(as.numeric(pyClone.tsv$alt_counts[a]) / as.numeric(pyClone.tsv$mutCopyNum[a] / best.CN)),
                            0.5
                        )$p.value < 0.05 & best.CN > 1) {
                            best.CN <- max(allCNs[!allCNs %in% best.CN])
                        }
                    }

                    pyClone.tsv$phyloCCF[a] <- pyClone.tsv$mutCopyNum[a] / best.CN
                    pyClone.tsv$expected.VAF[a] <- pyClone.tsv$expected.VAF[a] * best.CN
                    pyClone.tsv$phyloCCF.0.05[a] <- pyClone.tsv$phyloCCF.0.05[a] / best.CN
                    pyClone.tsv$phyloCCF.0.95[a] <- pyClone.tsv$phyloCCF.0.95[a] / best.CN
                }
            }
        }

        # finally, let's sort out 'missing' ones
        pyClone.tsv[pyClone.tsv$alt_counts == 0, "no.chrs.bearing.mut"] <- 0
        pyClone.tsv[pyClone.tsv$alt_counts == 0, "expected.VAF"] <- 0
        pyClone.tsv[pyClone.tsv$alt_counts == 0, "absolute.ccf"] <- 0
        return(pyClone.tsv)
    }
    region.phyloCCF2 <- suppressWarnings(calculate_phylo_ccf_withBH(region,
        complete.mutation.table = as.data.frame(pyclone.table),
        purity = sample.purity
    ))
    out2 <- data.table::as.data.table(region.phyloCCF2)
    data.table::setnames(
        out2, c("expected.VAF", "obs.VAF", "absolute.ccf", "absolute.ccf.0.05", "absolute.ccf.0.95", "phyloCCF", "phyloCCF.0.05", "phyloCCF.0.95"), c("expVAF", "obsVAF", "absCCF", "absCCF_lower", "absCCF_higher", "phyloCCF", "phyloCCF_lower", "phyloCCF_higher")
    )
    common_names <- intersect(names(out), names(out2))
    testthat::expect_equal(
        out[, ..common_names],
        out2[, ..common_names],
        tolerance = sqrt(.Machine$double.eps)
    )
})
