# User guide: `immunotation` package

### MHC molecules

MHC (major histocompatibility complex) molecules are a family of diverse cell surface complexes that present antigens to T cells. MHC molecules are divided into three classes (MHC class I, MHC class II, and non-classical MHC), which differ in their protein subunit composition and the types of receptors they can interact with. MHC class I molecules for example consist of one polymorphic $\alpha$-chain and one invariant $\beta$-chain ($\beta$2-microglobulin) and present peptide antigens to T cells that express the MHC-I specific co-receptor CD8. MHC class II molecules are typically composed of one $\alpha$- and one $\beta$-chain, which are both polymorphic. MHC class II molecules present peptide antigens to T cells that express the MHC-II specific co-receptor CD4. 

The repertoire of peptide antigens presented on MHC molecules depends on the sequence of the genes encoded in the MHC locus of an individual. Since the adaptive immune response to an invading pathogen relies on MHC-dependent antigen presentation, a high diversity of MHC genes on population level is beneficial from an evolutionary point of view. Indeed, MHC molecules are polygenic, which means that the MHC locus contains several different genes encoding MHC class I and MHC class II molecules. Moreover, MHC genes are polymorphic, which means that on a population level, there exist multiple variants (alleles) of each gene.

Several experimental techniques exist to identify the different MHC genes, alleles and protein complexes. Protein complexes for example can be identified by binding of subtype-specific anti-MHC antibodies. The resulting information on the protein complex is called MHC serotype. Moreover, MHC genes and alleles can be identified by hybridization with sequence specific probes or by sequencing and mapping to reference databases. However, these techniques often  cover only specific regions of the MHC genes and thus do not allow a complete and unambiguous allele identification.

The [IPD-IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) [1] and [IPD-MHC](https://www.ebi.ac.uk/ipd/mhc/) [2] databases provide a reference of all known MHC genes and alleles in different species. A systematic classification of MHC genes and proteins is provided in the [MHC restriction ontology (MRO)](https://github.com/IEDB/MRO) [3]. The annotation functions in the `immunotation` package use the classification scheme provided by the MRO.

### Hyperpolymorphic HLA genes in human populations

In humans, MHC molecules are encoded by highly variable HLA (human leukocyte antigen) genes of the hyperpolymorphic HLA locus on chromosome 6. Up to date more than 28.000 different alleles have been registered in the [IPD-IMGT/HLA](http://hla.alleles.org/nomenclature/index.html) database [1].

**Nomenclature of HLA genes**

HLA genes and alleles are named according to rules defined by the [WHO Nomenclature Committee for Factors of the HLA System](http://hla.alleles.org/nomenclature/naming.html). The following scheme depicts the components of a complete HLA allele name. In a deprecated naming scheme used before 2010, the components of the naming scheme were not separated by ":".

HLA-(gene)\*(group):(protein):(coding region):(non-coding region)(suffix)

| Term        | Description          | Example  |
| ----------- |:-----------------:|:-----------------:|
| gene        | Sum to 1 (± 0.015; or ± 0.05) | A, B, C, DPA1, DPB1, ... |
| group       | group of HLA alleles with similar protein sequence (protein sequence homology) |  01 |
| protein     | all HLA alleles with same protein sequence  | 01 |
| coding region     | all HLA alleles with same DNA sequence in coding region | 01 |
| non-coding region     | all HLA alleles with same DNA sequence in non-coding region | 01 |
| suffix    | indicates changes expression level (e.g. N - not expressed, L - low surface expression)  | N |

For example:

* *HLA-A\*01:01* and *HLA-A\*01:02* - have a similar but not identical protein sequence
* *HLA-A\*01:01:01* and *HLA-A\*01:01:02* - have the same protein sequence but slightly different DNA sequences in coding region 

**Protein and gene groups**

[*G groups*](http://hla.alleles.org/alleles/g_groups.html) are groups of HLA alleles that have identical *nucleotide* sequences across the exons encoding the peptide binding domains.

[*P groups*](http://hla.alleles.org/alleles/g_groups.html) are groups of HLA alleles that have identical *protein* sequences in the peptide binding domains.

**MAC (Multiple allele codes)**

The National Marrow Donor Program (NMDP) uses [multiple allele codes (MAC)](https://bioinformatics.bethematchclinical.org/hla-resources/allele-codes/allele-code-lists/) to facilitate the reporting and comparison of HLA alleles [4]. MAC represent groups of HLA alleles and are useful when the HLA typing is ambiguous and does not allow to narrow down one single allele from a list of alleles. The `immunotation` packages provides automated access to the [MAC conversion tools](https://hml.nmdp.org/MacUI/) provided by NMDP.

### Variation of HLA alleles across human populations

The frequencies of individual HLA alleles varies substantially between worldwide human populations.
The [Allele Frequency Net Database (AFND)](http://www.allelefrequencies.net/) is a repository for immune gene frequencies in different populations worldwide [5]. In addition to a large collection of HLA allele frequency datasets, the database also contains datasets for allele frequencies of KIR (Killer Cell Immunoglobulin-like Receptor) genes, MIC (Major histocompatibility complex class I chain related) and cytokine genes. The current version of the `immunotation` package allows the automated R access to the HLA related datasets in AFND.

The HLA frequency datasets in AFND are classified according to the following standards:

| Criteria        | Gold standard           | Silver standard  | Bronze standard |
| ----------- |:-----------------:|:-----------------:|:----------:|
| Allele frequency      | Sum to 1 (± 0.015; or ± 0.05) | Sum to 1 (± 0.015; or ± 0.05) | other |
| Sample size     | >= 50 individuals      |   any | any |
| Resolution of allele frequency | four or more digits      | two or more digits | other |

## Scope of the package

The `immunotation` package provides tools for consistent annotation of HLA genes, alleles and protein complexes. The package currently has two main functional modules: 

*1. Conversion and nomenclature functions:*

* access to annotations of MHC loci, genes, protein complexes and serotypes in different species
* conversion between different levels of HLA annotation, namely specific functions to convert between input and output of different immunoinformatics tools
* mapping of alleles and G groups and P groups.
* mapping of HLA allele groups and MAC

*2. Access to HLA allele frequencies:*

* automated access to datasets stored in the [AFND](http://www.allelefrequencies.net/)
* querying and visualization of HLA allele frequencies in human populations worldwide

### Installation

Install `immunotation` from github using the devtools library:

```{r, eval=FALSE}
library(devtools)
devtools::install_github("imkeller/immunotation")
```

**For further details please refer to the user guide located at vignette/immunotation_vignette.Rmd**

## References 

[1] Robinson J, Barker DJ, Georgiou X et al. IPD-IMGT/HLA Database. Nucleic Acids Research (2020)

[2] Maccari G, Robinson J, Ballingall K et al. IPD-MHC 2.0: an improved inter-species database for the study of the major histocompatibility complex. Nucleic Acids Research (2017)

[3] Vita R, Overton JA, Seymour E et al. An ontology for major histocompatibility restriction. J Biomed Semant (2016).

[4] Milius RP, Mack SJ, Hollenbach JA et al. Genotype List String: a grammar for describing HLA and KIR genotyping results in a text string. Tissue Antigens (2013). 

[5] Gonzalez-Galarza FF, McCabe A, Santos EJ at al. Allele frequency net database (AFND) 2020 update: gold-standard data classification, open access genotype data and new query tools. Nucleic Acid Research (2020).
