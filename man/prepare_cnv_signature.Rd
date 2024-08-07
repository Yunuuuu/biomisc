% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnv_matrix.R
\name{prepare_cnv_signature}
\alias{prepare_cnv_signature}
\title{Generate CN summary matrix from CN profiles}
\usage{
prepare_cnv_signature(
  seg_data,
  sample_field = NULL,
  chr_field = "chr",
  start_field = "startpos",
  end_field = "endpos",
  major_cn_field = "major_cn",
  minor_cn_field = "minor_cn"
)
}
\arguments{
\item{seg_data}{A \link{data.frame} obeject with segmented chromosome
data.}

\item{sample_field}{A string indicates the sample Id column in seg_data.}

\item{chr_field, start_field, end_field}{A string specifying the column of the
chromosome name, start positions and end positions of the genomic ranges in
seg_data.}

\item{minor_cn_field, major_cn_field}{A string specifying the minor_cn (minor
allele) and major_cn (major allele) column in seg_data.}
}
\description{
Generate CN summary matrix from CN profiles
}
\references{
Christopher D. Steele, Maxime Tarabichi, et al, Undifferentiated
Sarcomas Develop through Distinct Evolutionary Pathways,Cancer Cell,
\url{https://doi.org/10.1016/j.ccell.2019.02.002}
}
