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

# get the list of valid alleles
retrieve_AFND_valid_allele_list <- function(locus) {
   #check_hla_locus(locus)
    query_alleles_url <- sprintf("http://www.allelefrequencies.net/hla6006c.asp?hla_locus=%s&hla_locus_type=Classical", locus)
    
    html_input <- read_url(query_alleles_url)
    
    rvest_tables <- html_table(html_input, fill = TRUE)
    valid_alleles <- rvest_tables[[3]]$X1
    valid_alleles
}

# store the list of valid alleles somewhere, becasue it takes too long to build it every time

# Update the stored valid allele list 
update_valid_alleles <- function () {
    valid_alleles <- retrieve_AFND_valid_allele_list(locus = "A")
    valid_alleles_B <- retrieve_AFND_valid_allele_list(locus = "B")
    valid_alleles_C <- retrieve_AFND_valid_allele_list(locus = "C")
    valid_alleles_DPA1 <- retrieve_AFND_valid_allele_list(locus = "DPA1")
    valid_alleles_DPB1 <- retrieve_AFND_valid_allele_list(locus = "DPB1")
    valid_alleles_DQA1 <- retrieve_AFND_valid_allele_list(locus = "DQA1")
    valid_alleles_DQB1 <- retrieve_AFND_valid_allele_list(locus = "DQB1")
    valid_alleles_DRB1 <- retrieve_AFND_valid_allele_list(locus = "DRB1")
    
    complete_list <- c(valid_alleles, 
                       valid_alleles_B,
                       valid_alleles_C,
                       valid_alleles_DPA1,
                       valid_alleles_DPB1,
                       valid_alleles_DQA1,
                       valid_alleles_DQB1,
                       valid_alleles_DRB1)
    saveRDS(complete_list, "inst/extdata/AFND_valid_alleles.RData")
    return(complete_list)
}

get_valid_allele_list <- function() {
    file_name <- "/home/katharina/Documents/synced/science_projects/method_development/comp_immuno/immunotation/inst/extdata/AFND_valid_alleles.RData"
    if (file.exists(file_name)) {
        complete_list <- readRDS(file_name)
    } else {
        complete_list <- update_valid_alleles()
    }
    complete_list
}

# convert a allele into all possible alleles contained in this allele
# e.g. A*01:01 -> A*01:01:01, A*01:01:02, A*01:01:03
#' build_allele_group
#'
#' @param allele_selection 
#'
#' @return
#' @export
#'
#' @examples
build_allele_group <- function(allele_selection) {
    # find all alleles that correspond to allele
    # get valid alleles
    valid_alleles <- get_valid_allele_list()
    
    # expand to alleles in the same p group
    p_group_name <- get_P_group(allele_selection)
    alleles_p_group <- get_p_group_members(p_group_name)

     # find alleles corresponding to selection
    valid_p_alleles <- intersect(union(str_replace(p_group_name, "P", ""), alleles_p_group), 
                                       valid_alleles)
    
    p_expand <- sapply(valid_p_alleles, function (X) {valid_alleles[
        # replace * by \\* to make the regexp work properly
        grepl(str_replace(X, "\\*", "\\\\*"), valid_alleles)]
        })
    
    union(valid_p_alleles, union(unlist(p_expand), allele_selection))
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
    
    # maximum 20 alleles can be in the url (last length(hla_selection) because the last chunk missing)
    breaks <- c(seq(from = 0, to = length(hla_selection), by =19),length(hla_selection))
    
    # create empty data frame
    allele_df <- data.frame()
    for (i in seq(from = 1, to = length(breaks)-1)) {
    
        queryurl <- assemble_url_freq(hla_locus,
                                      hla_selection[(breaks[i]+1):(breaks[i+1])],
                                      hla_population,
                                      hla_sample_size_pattern,
                                      hla_sample_size,
                                      standard)
        allele_df <- rbind(read_complete_allele_freq_table(queryurl), allele_df)
    }
    allele_df
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
