% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/facets-wrapper.R
\name{readSnpMatrix}
\alias{readSnpMatrix}
\title{Read in the snp-pileup generated SNP read count matrix file}
\usage{
readSnpMatrix(
  filename,
  skip = 0L,
  err.thresh = Inf,
  del.thresh = Inf,
  perl.pileup = FALSE,
  MandUnormal = FALSE,
  spanT = 0.2,
  spanA = 0.2,
  spanX = 0.2,
  gbuild = "hg19",
  ReferencePileupFile = NULL,
  ReferenceLoessFile = NULL,
  MinOverlap = 0.9,
  useMatchedX = FALSE,
  refX = FALSE,
  unmatched = FALSE,
  donorCounts = FALSE
)
}
\arguments{
\item{filename}{counts file from snp-pileup}

\item{skip}{(character) Skip n number of lines in the input file.}

\item{err.thresh}{(numeric) Error threshold to be used to filter snp-pileup data frame.}

\item{del.thresh}{(numeric) Deletion threshold to be used to filter snp-pileup data frame.}

\item{perl.pileup}{(logical) Is the pileup data generated using perl pileup tool?}

\item{MandUnormal}{(logical) Is CNLR analysis to be peformed using unmatched reference normals?}

\item{spanT}{(numeric) Default span value to be used for loess normalization in tumor sample.}

\item{spanA}{(numeric) Default span value to be used for loess normalization across autosomal chromosomes in the normal sample.}

\item{spanX}{(numeric) Default span value to be used for loess normalization in Chr X in the normal sample.}

\item{gbuild}{(character) Genome build (Default: hg19).}

\item{ReferencePileupFile}{(character) Filepath to an optional snp-pileup generated pileup data of one or more reference normals.}

\item{ReferenceLoessFile}{(character) Filepath to an optional loess data, generated using the facets2n package, of one or more reference normals. The number of normals in this data should match that in the ReferencePileupFile, and should be in the same order.}

\item{MinOverlap}{(numeric) Mininum overlap fraction of loci between a tumor pileup and reference pileup data.}

\item{useMatchedX}{(logical) Is the matched normal to be used for ChrX normalization?}

\item{refX}{(logical) Use matched or reference normal for chrX normalization. excludes unmatched normals, such as pooled references, present in tumor counts matrix.}

\item{unmatched}{(logical) is the tumor being analyzed unmatched}

\item{donorCounts}{(logical) is the counts matrix baseline donor sample(s)}
}
\value{
A dataframe of pileup depth values for Tumor and Matched Normal if MandUnormal is FALSE. Else, a list of data frame with pileup depth values of Tumor, matched Normal, and a best unmatched normal, and the associated span values.
}
\description{
Read in the snp-pileup generated SNP read count matrix file
}
