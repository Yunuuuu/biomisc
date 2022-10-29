.onAttach <- function(libname, pkgname) {
    version <- utils::packageDescription(pkgname, fields = "Version")

    packageStartupMessage(
        "Package: ", pkgname, " (version: ", version, ")\n",
        "Miscellaneous functions for Bioinformatics!"
    )
    invisible()
}

if (getRversion() >= "2.15.1") {
    utils::globalVariables(
        c( # run_absolute
            "Chromosome", "Sample", "Tumor_Sample_Barcode",
            "group_id",
            # get_arm_cytoband
            "chr", "seq_chr", "seq_int", "arm", "chr_arm_order",
            # prepare_pyclone
            "i.ref_counts", "i.var_counts", "x.minor_cn", 
            "x.major_cn", "i.chromosome", "i.pos",
            "x.start_pos", "x.end_pos",
            "mutation_id", "position", "chromosome", "tumour_content",
            "major_cn"
        )
    )
}
