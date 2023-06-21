segRData <- load("../testdata/clonal_dissection/ASCAT/D_LMS025_T1.ascat.seg.RData")
seg.mat.copy <- get(segRData)
seg.mat.copy$CNt <- seg.mat.copy$nMajor + seg.mat.copy$nMinor
seg.mat.copy$CNV <- log2(seg.mat.copy$CNt / seg.mat.copy$Ploidy)

test_that("run_arm_cnv (abs mode) works", {
    testthat::expect_snapshot_value(
        run_arm_cnv(seg.mat.copy,
            cnv_field = "CNt",
            cnv_mode = "abs", ploidy_field = "Ploidy"
        ),
        style = "serialize"
    )
})

test_that("run_arm_cnv (rel mode) works", {
    testthat::expect_snapshot_value(
        run_arm_cnv(seg.mat.copy,
            cnv_field = "CNV",
            cnv_mode = "rel"
        ),
        style = "serialize"
    )
})
