convert_string <- function(value) {
    # we need this function to format the lists in the sting
    # and replace the NA values in the URL with ""
    if (length(value) > 1) {value <- paste(value, collapse=",")}
    else if(is.na(value)) {""}
    else value
}

# FUNCTIONS FOR ALLELE QUERY

assemble_url_freq <- function(hla_locus,
                              hla_selection,
                              hla_population,
                              hla_sample_size_pattern,
                              hla_sample_size,
                              standard,
                              hla_locus_type = "Classical") {

    freq_url_root <- "http://www.allelefrequencies.net/hla6006a.asp?hla_locus_type=%s&hla_locus=%s&hla_selection=%s&hla_population=%s&hla_sample_size_pattern=%s&hla_sample_size=%s&standard=%s"

    sprintf(freq_url_root, hla_locus_type,
            convert_string(hla_locus),
            convert_string(hla_selection),
            convert_string(hla_population),
            convert_string(hla_sample_size_pattern),
            convert_string(hla_sample_size),
            standard)
}

#' Query allele frequencies
#'
#' @param hla_locus HLA locus that will be used for filtering data. A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1
#' @param hla_selection Allele that will be used for filtering data. e.g. A*01:01
#' @param hla_population Numeric identifier of the population that will be used for filtering. Thie identifier is defined by the Allele Frequency Net Database.
#' @param hla_sample_size_pattern Keyword used to define the filtering for a specific population size. e.g. "bigger_than", "equal", "less_than", "less_equal_than", "bigger_equal_than"
#' @param hla_sample_size Integer number used to define the filtering for a specific population size, together with the hla_sample_size_pattern argument.
#' @param standard Population standards, as defined in the package vignette. "g" - gold, "s" - silver, "a" - all
#'
#' @return data.frame object containing the result of the allele frequency query
#' @export
#'
#' @examples library(AFNDquery)
#'
#' # select frequencies of the A*02:01 allele,
#' # for gold standard population with more than 10,000 individuals
#' sel <- query_allele_frequencies(hla_selection = "A*02:01",
#' hla_sample_size_pattern = "bigger_than", hla_sample_size = 10000,
#' standard="g")
#'
query_allele_frequencies <- function(
    hla_locus = NA,
    hla_selection = NA,
    hla_population = NA,
    hla_sample_size_pattern = NA,
    hla_sample_size = NA,
    standard="a") {
    # check whether input parameters are valid
    verify_parameters(hla_locus,
                      hla_selection,
                      hla_population,
                      hla_sample_size_pattern,
                      hla_sample_size,
                      standard = "a")

    queryurl <- assemble_url_freq(hla_locus,
                             hla_selection,
                             hla_population,
                             hla_sample_size_pattern,
                             hla_sample_size,
                             standard)

    read_complete_allele_freq_table(queryurl)
}

# FUNCTIONS FOR POPULATION QUERY

assemble_url_pop <- function(population_id) {

    pop_url_root <- "http://www.allelefrequencies.net/pop6001c.asp?pop_id=%s"

    sprintf(pop_url_root, population_id)
}

query_single_population_detail <- function(population_id) {
    check_population(population_id)
    pop_url <- assemble_url_pop(population_id)
    read_population_detail(pop_url, population_id)
}

#' Query population metainformation
#'
#' @param population_ids List of numeric identifiers of the population that will be used for filtering. Thie identifier is defined by the Allele Frequency Net Database.
#'
#' @return data.frame object containing the result of the population detail query
#' @export
#'
#' @examples library(AFNDquery)
#' population_detail <- query_population_detail(0001986)
query_population_detail <- function(population_ids) {
    pop_df <- data.frame()
    for(pop_id in population_ids) {
        pop_df <- rbind(pop_df, query_single_population_detail(pop_id))
    }
    pop_df
}
