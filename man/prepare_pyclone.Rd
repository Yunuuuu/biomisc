% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_pyclone.R
\name{prepare_pyclone}
\alias{prepare_pyclone}
\title{Prepare data for the input of PyClone or PyClone-vi}
\usage{
prepare_pyclone(
  mut_data,
  cnv_data,
  on_sample = NULL,
  on_chr = "chr",
  mut_pos = "pos",
  ref_field = "ref",
  start_field = "start",
  end_field = "end",
  major_cn_field = "major_cn",
  minor_cn_field = "minor_cn",
  ref_counts_field = "ref_counts",
  var_counts_field = "var_counts",
  purity_field = NULL,
  normal_cn = 2L,
  pyclone_vi = FALSE,
  error_rate = NULL,
  nomatch = NULL
)
}
\arguments{
\item{mut_data}{A data.frame of mutation data.}

\item{cnv_data}{A data.frame of allele-specific CNV data with the top 5
columns containing "chromosome", "start_pos", "end_pos", "major_cn", and
"minor_cn". "sample_id" is optional, details see mut_data. other columns will
be omited. Column names don't matter.}

\item{on_sample}{A string (which can be named), specifying the column in the
sample used for matching mut_data and cnv_data. If NULL, each patient is
considered to contain only one sample.}

\item{on_chr}{A string (can be named), specifying the chromosome column used
to match mut_data and cnv_data.}

\item{mut_pos}{A string indicating the column names in \code{mut_data} that
contains the variants positions in the chromosome.}

\item{ref_field}{A string indicating the column names in \code{mut_data} that
contains the normal reference allele.}

\item{start_field, end_field}{A string indicating the column names in
\code{cnv_data} that contains the start positions and end position of the genomic
ranges.}

\item{major_cn_field, minor_cn_field}{A string indicating the column names in
\code{cnv_data} that contains the major_cn and minor_cn.}

\item{ref_counts_field, var_counts_field}{A string indicating the column names
in \code{mut_data} that contains the ref_counts and var_counts.}

\item{purity_field}{Only used for
\href{https://github.com/Roth-Lab/pyclone-vi}{PyClone-vi}, a string specifying the
tumour content (cellularity) column (can in \code{mut_data} or \code{cnv_data}) of the
sample. Default value is 1.0 if NULL.}

\item{normal_cn}{The copy number of the locus in non-malignant cells. This
should generally be 2 except for sex chromosomes in males.}

\item{pyclone_vi}{a logical value indicates whether prepare data for the
input of PyClone-vi. Details see \url{https://github.com/Roth-Lab/pyclone-vi}}

\item{error_rate}{Only used for
\href{https://github.com/Roth-Lab/pyclone-vi}{PyClone-vi}, sequencing error rate.}

\item{nomatch}{When a row in \code{mut_data} has no match to \code{cnv_data},
nomatch=NA means NA is returned. NULL (or 0 for backward compatibility) means
no rows will be returned for that row of \code{mut_data}.}
}
\description{
Prepare data for the input of PyClone or PyClone-vi
}
\references{
\itemize{
\item \url{https://github.com/Roth-Lab/pyclone}
\item \url{https://github.com/Roth-Lab/pyclone-vi}
\item \url{https://bitbucket.org/sequenzatools/sequenza/src/v2.1.1/R/next.R}
}
}
