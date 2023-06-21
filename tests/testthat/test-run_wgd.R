segRData <- load("../testdata/clonal_dissection/ASCAT/D_LMS025_T1.ascat.seg.RData")
seg.mat.copy <- get(segRData)

test_that("run_wgd (Craig M. mode) works", {
    testthat::expect_snapshot_value(
        run_wgd(
            seg.mat.copy,
            sample_field = "sample",
            major_cn_field = "nMajor", ploidy_field = "Ploidy"
        ),
        style = "serialize"
    )
})

test_that("run_wgd (Sally M. mode) works", {
    testthat::expect_snapshot_value(
        run_wgd(
            seg.mat.copy,
            sample_field = "sample",
            major_cn_field = "nMajor", minor_cn_field = "nMinor",
            ploidy_field = "Ploidy"
        ),
        style = "serialize"
    )
})
