% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MRO_interface.R
\name{retrieve_chain_lookup_table}
\alias{retrieve_chain_lookup_table}
\title{Retrieve MHC chain lookup table}
\usage{
retrieve_chain_lookup_table(organism)
}
\arguments{
\item{organism}{name of organism (e.g. "human")}
}
\value{
Table containing MHC chain information for the organism. 
It contains chain names, MHC restriction and protein sequence.
}
\description{
Retrieve MHC chain lookup table
}
\examples{
retrieve_chain_lookup_table("mouse")

}
