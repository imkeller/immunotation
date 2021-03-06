% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AFND_queries.R
\name{query_allele_frequencies}
\alias{query_allele_frequencies}
\title{Query allele frequencies}
\usage{
query_allele_frequencies(
  hla_locus = NA,
  hla_selection = NA,
  hla_population = NA,
  hla_country = NA,
  hla_region = NA,
  hla_ethnic = NA,
  hla_sample_size_pattern = NA,
  hla_sample_size = NA,
  standard = "a"
)
}
\arguments{
\item{hla_locus}{HLA locus that will be used for filtering data. A, B, C, 
DPA1, DPB1, DQA1, DQB1, DRB1}

\item{hla_selection}{Allele that will be used for filtering data. 
e.g. A*01:01}

\item{hla_population}{Numeric identifier of the population that will be used 
for filtering. This identifier is defined by the Allele Frequency 
Net Database.}

\item{hla_country}{Country of interest (e.g. Germany, France, ...).}

\item{hla_region}{Geographic region of interest 
(e.g. Europe, North Africa, ...)}

\item{hla_ethnic}{Ethnic origin of interest (e.g. Caucasoid, Siberian, ...)}

\item{hla_sample_size_pattern}{Keyword used to define the filtering for a 
specific population size. e.g. "bigger_than", "equal", "less_than", 
"less_equal_than", "bigger_equal_than"}

\item{hla_sample_size}{Integer number used to define the filtering for a 
specific population size, together with the hla_sample_size_pattern argument.}

\item{standard}{Population standards, as defined in the package vignette. 
"g" - gold, "s" - silver, "a" - all}
}
\value{
data.frame object containing the result of the allele frequency query
}
\description{
Query allele frequencies
}
\examples{

# select frequencies of the A*02:01 allele,
# for gold standard population with more than 10,000 individuals
sel <- query_allele_frequencies(hla_selection = "A*02:01",
hla_sample_size_pattern = "bigger_than", hla_sample_size = 10000,
standard="g")

}
