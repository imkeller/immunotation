% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AFND_queries.R
\name{build_allele_group}
\alias{build_allele_group}
\title{Building a list of alleles to cover}
\usage{
build_allele_group(allele_selection)
}
\arguments{
\item{allele_selection}{HLA allele for whicht
the allele group should be built.}
}
\value{
list of alleles
}
\description{
\code{build_allele_group} e.g. A*01:01 -> A*01:01:01, 
A*01:01:02, A*01:01:03
}
\examples{
build_allele_group("A*01:01")
}
