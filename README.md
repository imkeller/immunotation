# Accessing HLA allele frequencies using `AFNDquery`

## Introduction

### Hyperpolymorphic HLA genes in human populations

MHC (major histocompatibility complex) molecules are cell surface complexes that present peptide antigens to T cells. This processed, also known as MHC-dependent antigen presentation, is essential for the coordination of adaptive immune responses. In humans, MHC molecules are encoded by highly variable HLA (human leukocyte antigen) genes encoded in the hyperpolymorphic HLA locus on chromosome 6. The library of peptide antigens presented on MHC molecules highly depends on the HLA alleles genes present in an individual. The frequencies of individual HLA alleles varies substantially between worldwide human populations.

### Databases for HLA allele frequencies in human populations

The Allele Frequency Net Database (AFND, http://www.allelefrequencies.net/) is a repository for immune gene frequencies in different populations worldwide (Gonzalez-Galarza et al. 2020). In addition to a large collection of HLA alle frequency datasets, the database also contains datasets for allele frequencies of KIR (Killer Cell Immunoglobulin-like Receptor) genes, MIC (Major histocompatibility complex class I chain related) and cytokine genes. The current version of AFNDquery allows the automated R access to the HLA related datasets in AFND.

The HLA frequency datasets in AFND are classified according to the following standards:

| Criteria        | Gold standard           | Silver standard  | Bronze standard |
| ----------- |:-----------------:|:-----------------:|:----------:|
| Allele frequency      | Sum to 1 (± 0.015; or ± 0.05) | Sum to 1 (± 0.015; or ± 0.05) | other |
| Sample size     | >= 50 individuals      |   any | any |
| Resolution | four or more digits      | two or more digits | other |

## Scope of the `AFNDquery` package

The `AFNDquery` package currently supports the following type of queries related to HLA allele frequencies:

1. Query the **HLA allele frequencies** for a given HLA locus, population, sample size (and other available filters). This module returns information contained in the *HLA* section of AFND (allele frequency, population name and id, sample size).

2. Query the **population metainformation** for a given population. This module returns information contained in the *Population* section of the Allele Frequency Net Database (ethnic origin, geographic location, sampling size and date, ...).

## Installation

Install from github using the devtools library:

```{r}
library(devtools)

devtools::install_github("imkeller/AFNDquery")
```

## Querying allele frequencies

**Example 1:** Querying the frequency of allele A*02:01 in all populations with more than 10,000 individuals and fulfilling the gold standard.

```{r}
sel1 <- query_allele_frequencies(hla_selection = "A*02:01", 
                                hla_sample_size_pattern = "bigger_than", 
                                hla_sample_size = 10000, 
                                standard="g")

DT::datatable(sel1)
```

**Example 2:** Querying the frequencies of alleles within the HLA-B locus for population "Peru Lamas City Lama" with population_id 1986. If you are not sure which population_id you should use, we recommend to search your population of interest in the *Population* section of the Allele Frequency Net Database and extract the population id from the url (last 4 digits).

```{r}
sel2 <- query_allele_frequencies(hla_locus = "B", hla_population = 1986)

DT::datatable(sel2)
```

## Querying population metainformation

**Example 3:** Query the metainformation concerning population "Peru Lamas City Lama" (population_id 1986). The webpage concerning the queried information for population Peru Lamas City Lama (1986) can be found here: http://www.allelefrequencies.net/pop6001c.asp?pop_id=1986

```{r}
sel3 <- query_population_detail(1986)

DT::datatable(sel3, options = list(scrollX = TRUE))
```

**Example 4:** Query the metainformation concerning the populations that were listed in the table returned by Example 1

```{r}
sel4 <- query_population_detail(as.numeric(sel1$population_id))

# only select the first 5 columns to discplay in table
DT::datatable(sel4[1:5], options = list(scrollX = TRUE))
```

## References 

1. Gonzalez-Galarza FF, McCabe A, Santos EJ, Jones J, Takeshita LY, Ortega-Rivera ND, Del Cid-Pavon GM, Ramsbottom K, Ghattaoraya GS, Alfirevic A, Middleton D and Jones AR. Allele frequency net database (AFND) 2020 update: gold-standard data classification, open access genotype data and new query tools. Nucleic Acid Research 2020, 48, D783-8.
