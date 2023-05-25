testthat::test_that("multiplication works", {
    maf <- readRDS(system.file("extdata", "absolute",
        "run_absolute_example_maf.rds",
        package = "biomisc"
    ))
    testthat::expect_snapshot(
        run_motif_fisher(
            maf[, c("Tumor_Sample_Barcode", "Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele1")],
            ref_genome = "hg38"
        )
    )
})
