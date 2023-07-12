test_that("prepare_snv_signature works well", {
    maf <- readRDS(system.file("extdata", "absolute",
        "run_absolute_example_maf.rds",
        package = "biomisc"
    ))
    ref_genome <- get_genome("hg19")
    testthat::expect_snapshot_value(
        suppressWarnings(prepare_snv_signature(
            maf, "Tumor_Sample_Barcode",
            ref_genome = ref_genome,
            "Chromosome", mut_pos = "Start_position",
            ref_field = "Reference_Allele", alt_field = "Tumor_Seq_Allele2"
        )),
        style = "serialize"
    )
    testthat::expect_snapshot_value(
        suppressWarnings(prepare_snv_signature(
            maf, "Tumor_Sample_Barcode",
            ref_genome = ref_genome,
            "Chromosome", mut_pos = "Start_position",
            ref_field = "Reference_Allele", alt_field = "Tumor_Seq_Allele2",
            method = "genome"
        )),
        style = "serialize"
    )
})
