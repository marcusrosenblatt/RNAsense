% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{MZsox}
\alias{MZsox}
\title{Time resolved RNA seq data for early zygotic development of zebra fish.}
\format{A \link{SummarizedExperiment} object with 15775 rows and 40 columns. Arguments rowData and colData give covariate information on the count data as follows:
\itemize{
  \item{rowData: }{name, name of the gene}
  \item{rowData: }{genename, identifier of gene}
  \item{colData: }{condition, name of condition}
  \item{colData: }{time, measurement time in hours post fertilization}
  \item{colData: }{replicate, identifier of replicate}
}}
\usage{
MZsox
}
\description{
A dataset containing RNA seq count data for two experimental conditions at 8 different time points including replicates.
}
\keyword{datasets}
