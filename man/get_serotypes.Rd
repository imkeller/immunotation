% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRO_interface.R
\name{get_serotypes}
\alias{get_serotypes}
\title{Serotypes}
\usage{
get_serotypes(allele_list, organism = "human", mhc_type)
}
\arguments{
\item{allele_list}{List of allele}

\item{organism}{Organism to be used for MRO lookup. 
If the organism does not match the given allele, a empty object is returned.}

\item{mhc_type}{["MHC-I" or "MHC-II"] MHC class to use for MRO lookup.}
}
\value{
Named list of serotypes, which only contains complexes contained 
in the MRO. 
If no serotype is annoted for a given complex, the list element is NA.
}
\description{
Get the serotypes of the MHC complexes encoded by a list of 
MHC alleles.
}
\examples{
allele_list <- c("A*01:01:01","B*27:01")
get_serotypes(allele_list, mhc_type = "MHC-I")

}
