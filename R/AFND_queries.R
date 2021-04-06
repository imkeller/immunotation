
convert_string <- function(value) {
    # we need this function to format the lists in the sting
    # and replace the NA values in the URL with ""
    if (length(value) > 1) {value <- paste(value, collapse=",")}
    else if(is.na(value)) {""}
    else value
}

# get the list of valid alleles
retrieve_AFND_valid_allele_list <- function(locus) {
    check_hla_locus(locus)
    query_alleles_url <- sprintf(paste0("http://www.allelefrequencies.net/",
    "hla6006c.asp?hla_locus=%s&hla_locus_type=Classical"), locus)
    
    html_input <- getURL(query_alleles_url, read_method = "html")
    
    rvest_tables <- rvest::html_table(html_input, fill = TRUE)
    valid_alleles <- rvest_tables[[3]]$X1
    valid_alleles
}

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
    
    complete_list <- c(valid_alleles, valid_alleles_B,
        valid_alleles_C, valid_alleles_DPA1, valid_alleles_DPB1,
        valid_alleles_DQA1, valid_alleles_DQB1, valid_alleles_DRB1)
    saveRDS(complete_list, system.file("extdata", "AFND_valid_alleles.RData",
        package = "immunotation"))
    return(complete_list)
}

get_valid_allele_list <- function() {
    file_name <- system.file("extdata", "AFND_valid_alleles.RData",
        package = "immunotation")
    if (file.exists(file_name)) {
        complete_list <- readRDS(file_name)
    } else {
        complete_list <- update_valid_alleles()
    }
    complete_list
}

# convert a allele into all possible alleles contained in this allele
# e.g. A*01:01 -> A*01:01:01, A*01:01:02, A*01:01:03
#' @title  Building a list of alleles to cover 
#' @description \code{build_allele_group} e.g. A*01:01 -> A*01:01:01, 
#' A*01:01:02, A*01:01:03
#'
#' @param allele_selection HLA allele for whicht
#' the allele group should be built.
#'
#' @return list of alleles
#' @export
#'
#' @examples
#' build_allele_group("A*01:01")
build_allele_group <- function(allele_selection) {
    # find all alleles that correspond to allele
    # get valid alleles
    valid_alleles <- get_valid_allele_list()
    
    # expand to alleles in the same p group
    p_group_name <- get_P_group(allele_selection)
    alleles_p_group <- get_p_group_members(p_group_name)
    
    # find alleles corresponding to selection
    valid_p_alleles <- intersect(union(
        stringr::str_replace(p_group_name, "P", ""), 
        alleles_p_group), valid_alleles)
    
    p_expand <- unlist(lapply(valid_p_alleles, function (X) {valid_alleles[
            # replace * by \\* to make the regexp work properly
            grepl(stringr::str_replace(X, "\\*", "\\\\*"), valid_alleles)]
        }))
    
    union(valid_p_alleles, union(unlist(p_expand), allele_selection))
}

# FUNCTIONS FOR ALLELE QUERY

assemble_url_allele_freq <- function(hla_locus,
    hla_selection, hla_population, hla_country, hla_region,
    hla_ethnic, hla_sample_size_pattern, hla_sample_size,
    standard, hla_locus_type = "Classical") {

    freq_url_root <- paste0("http://www.allelefrequencies.net/hla6006a.asp?",
    "hla_locus_type=%s&hla_locus=%s&hla_selection=%s&hla_population=%s&",
    "hla_country=%s&hla_dataset=&hla_region=%s&hla_ethnic=%s&hla_study=&",
    "hla_order=order_1&hla_sample_size_pattern=%s&",
    "hla_sample_size=%s&standard=%s")

    sprintf(freq_url_root, hla_locus_type,
        convert_string(hla_locus), convert_string(hla_selection),
        convert_string(hla_population), convert_string(hla_country),
        convert_string(hla_region), convert_string(hla_ethnic),
        convert_string(hla_sample_size_pattern),
        convert_string(hla_sample_size), standard)
}


#' Query allele frequencies
#'
#' @param hla_locus HLA locus that will be used for filtering data. A, B, C, 
#' DPA1, DPB1, DQA1, DQB1, DRB1
#' @param hla_selection Allele that will be used for filtering data. 
#' e.g. A*01:01
#' @param hla_population Numeric identifier of the population that will be used 
#' for filtering. This identifier is defined by the Allele Frequency 
#' Net Database.
#' @param hla_country Country of interest (e.g. Germany, France, ...).
#' @param hla_region Geographic region of interest 
#' (e.g. Europe, North Africa, ...)
#' @param hla_ethnic Ethnic origin of interest (e.g. Caucasoid, Siberian, ...)
#' @param hla_sample_size_pattern Keyword used to define the filtering for a 
#' specific population size. e.g. "bigger_than", "equal", "less_than", 
#' "less_equal_than", "bigger_equal_than"
#' @param hla_sample_size Integer number used to define the filtering for a 
#' specific population size, together with the hla_sample_size_pattern argument.
#' @param standard Population standards, as defined in the package vignette. 
#' "g" - gold, "s" - silver, "a" - all
#'
#' @return data.frame object containing the result of the allele frequency query
#' @export
#'
#' @examples
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
    hla_country = NA,
    hla_region = NA,
    hla_ethnic = NA,
    hla_sample_size_pattern = NA,
    hla_sample_size = NA,
    standard = "a") {
    
    # check whether input parameters are valid
    verify_parameters(hla_locus = hla_locus, hla_selection = hla_selection,
        hla_population = hla_population, hla_country = hla_country,
        hla_region = hla_region, hla_ethnic = hla_ethnic,
        hla_sample_size_pattern = hla_sample_size_pattern,
        hla_sample_size = hla_sample_size, standard = standard,
        query_type ="allele")
    
    # maximum 20 alleles can be in the url (last length(hla_selection) 
    # because the last chunk missing)
    breaks <- c(seq(from = 0, to = length(hla_selection), by =19),
        length(hla_selection))
    
    # create empty data frame
    allele_df <- data.frame()
    for (i in seq(from = 1, to = length(breaks)-1)) {
    
        queryurl <- assemble_url_allele_freq(hla_locus,
            hla_selection[(breaks[i]+1):(breaks[i+1])], hla_population,
            hla_country, hla_region, hla_ethnic, hla_sample_size_pattern,
            hla_sample_size, standard)
        allele_df <- rbind(read_complete_freq_table(queryurl, type = "allele"),
            allele_df)
    }
    allele_df
}

# FUNCTIONS FOR HAPLOTYPES

assemble_haplotype_string_url_from_allele_list <- function(allele_list) {
    
    # this is the order in which loci need to be passed
    loci_regexp <- c("A\\*","B\\*","C\\*","DRB1\\*",
                        "DPA1\\*","DPB1\\*","DQA1\\*","DQB1\\*")
    # be default set to not include locus
    loci_names <- c("A_not", "B_not", "C_not", "DRB1_not", 
                        "DPA1_not", "DPB1_not", "DQA1_not", "DQB1_not")
    for (i in seq(length(loci_regexp))) {
        locus_name <- allele_list[grepl(loci_regexp[i], allele_list)]
        if (length(locus_name) > 1) {
            stop("In the list of haplotypes for the haplotype assembly, only",
            "one allele per locus may be passed. Following ",
            "entry causes conflict: ", locus_name)
        } else if (length(locus_name) != 0) {
            loci_names[i] <- locus_name
        }
    }
    hla_str <- sprintf(paste0("hla_locus1=%s&hla_locus2=%s&hla_locus3=%s&",
    "hla_locus4=%s&hla_locus5=%s&hla_locus6=%s&hla_locus7=%s&hla_locus8=%s"),
        loci_names[1], loci_names[2], loci_names[3],
        loci_names[4], loci_names[5], loci_names[6],
        loci_names[7], loci_names[8])
    hla_str
}
    

assemble_url_haplotype_freq <- function(hla_selection,
    hla_population, hla_country, hla_region, hla_ethnic,
    hla_sample_size_pattern, hla_sample_size) {
    
    hla_str <- assemble_haplotype_string_url_from_allele_list(hla_selection)
    freq_url_root <- stringr::str_c("http://www.allelefrequencies.net/",
    "hla6003a.asp?", hla_str, "&hla_population=%s&hla_country=%s&hla_dataset=&",
    "hla_region=%s&hla_ethnic=%s&hla_study=&hla_order=order_1",
    "&hla_sample_size_pattern=%s&hla_sample_size=%s&",
    "hla_sample_year_pattern=equal&hla_sample_year=&hla_loci=")
    
    sprintf(freq_url_root, 
        convert_string(hla_population), convert_string(hla_country),
        convert_string(hla_region), convert_string(hla_ethnic),
        convert_string(hla_sample_size_pattern), 
        convert_string(hla_sample_size))
}


#' Query haplotype frequencies
#'
#' @param hla_selection Alleles that will be used to build the haplotype query. 
#' One entry per locus. If no entry for a given locus, the function will search 
#' for haplotypes that do not include specifications for this locus. If any 
#' allele for a given locus should be considered, the list entry should be 
#' "A*" or other locus in same format.
#' @param hla_population Numeric identifier of the population that will be used 
#' for filtering. Thie identifier is defined by the Allele Frequency Net 
#' Database.
#' @param hla_country Country of interest (e.g. Germany, France, ...).
#' @param hla_region Geographic region of interest (e.g. Europe, 
#' North Africa, ...)
#' @param hla_ethnic Ethnic origin of interest (e.g. Caucasoid, Siberian, ...)
#' @param hla_sample_size_pattern Keyword used to define the filtering for a 
#' specific population size. e.g. "bigger_than", "equal", "less_than", 
#' "less_equal_than", "bigger_equal_than"
#' @param hla_sample_size Integer number used to define the filtering for a 
#' specific population size, together with the hla_sample_size_pattern argument.
#'
#' @return data.frame object containing the result of the allele frequency query
#' @export
#'
#' @examples
#' # works only for one haplotype at a time
#' query_haplotype_frequencies(hla_selection = c("A*02:01", "B*", "C*"), 
#' hla_region = "Europe")
#'
# only for one haplotype at a time
query_haplotype_frequencies <- function(
    # this selection is a selection of A,B,C.... alleles that will be assembled 
    # to haplotype
    hla_selection = NA,
    hla_population = NA,
    hla_country = NA,
    hla_region = NA,
    hla_ethnic = NA,
    hla_sample_size_pattern = NA,
    hla_sample_size = NA) {
    
    # check whether input parameters are valid
    verify_parameters(hla_locus = NA,
        hla_selection = hla_selection, hla_population = hla_population,
        hla_country = hla_country, hla_region = hla_region,
        hla_ethnic = hla_ethnic,
        hla_sample_size_pattern = hla_sample_size_pattern,
        hla_sample_size = hla_sample_size, query_type ="haplotype")

    queryurl <- assemble_url_haplotype_freq(hla_selection,
        hla_population, hla_country, hla_region, hla_ethnic,
        hla_sample_size_pattern, hla_sample_size)
    
    allele_df <- read_complete_freq_table(queryurl, type = "haplotype")
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
#' @param population_ids List of numeric identifiers of the population that 
#' will be used for filtering. The identifier is defined by the Allele 
#' Frequency Net Database.
#'
#' @return data.frame object containing the result of the population detail 
#' query
#' @export
#'
#' @examples
#' population_detail <- query_population_detail(0001986)
query_population_detail <- function(population_ids) {
    pop_df <- data.frame()
    for(pop_id in population_ids) {
        pop_df <- rbind(pop_df, query_single_population_detail(pop_id))
    }
    pop_df
}
