#
# FUNCTIONS TO READ ALLELE AND HAPLOTYPE FREQUENCY DATA
#

#' get_nb_pages
#' @param page_tbl html input page
#' @return integer number of pages in the html input
#' @keywords internal
get_nb_pages <- function(page_tbl) {
    stringr::str_extract(page_tbl$X6, "\\d*$")
}

#' parse_allele_freq_html 
#' @description format the allele frequency table and select columns of interest
#' @param allele_freq_table allele frequency table parsed from AFND website
#' @return allele_freq_table reformatted
#' @keywords internal
parse_allele_freq_html <- function(allele_freq_table) {
    colnames(allele_freq_table) <- c("line", "allele", "", "population",
        "perc_individuals_with_allele",
        "allele_frequency", "", "sample_size")
    allele_freq_table[c("allele","population","allele_frequency","sample_size")]
} 

#' parse_haplotype_freq_html 
#' @description format the haplotype frequency table and select 
#' columns of interest
#' @param haplotype_freq_table haplotype frequency table parsed
#' from AFND website
#' @return haplotype_freq_table reformatted
#' @keywords internal
parse_haplotype_freq_html <- function(haplotype_freq_table) {
    colnames(haplotype_freq_table) <- c("line", "haplotype", "", "population", 
        "perc_individuals_with_haplotype",
        "", "sample_size", "")
    haplotype_freq_table[c("haplotype","population",
        "perc_individuals_with_haplotype","sample_size")]
}

#' extract_population_id 
#' @description extract the population ids from the html result
#' @param data html from AFND website
#' @return population ids
#' @keywords internal
extract_population_id <- function(data) {
    pop_hrefs <- rvest::html_attr(
        rvest::html_nodes(
            rvest::html_nodes(data, "td:nth-child(4)"), "a"), "href")
    # on page > 1 there is another href on top and below
    filtered_pop_hrefs <- pop_hrefs[grep("pop6001", pop_hrefs)]
    stringr::str_extract(filtered_pop_hrefs , "\\d\\d\\d\\d$")
}

#' read_complete_freq_table 
#' @param url URL of the website containing frequency table
#' @param type ["allele"|"haplotype"] 
#' @return frequency table
#' @keywords internal
read_complete_freq_table <- function(url, type) {
    html_input <- getURL(url, read_method = "html")
    
    rvest_tables <- rvest::html_table(html_input, fill = TRUE)
    
    output_table <- data.frame()
    # get the number of pages that need to be read
    # for some instances where no entries found,the page number is no indicated 
    # return empty table in this case
    if(length(rvest_tables) > 3) {
        page_nb <- as.numeric(get_nb_pages(rvest_tables[[4]]))
    } else {return(output_table)}
    
    if (length(page_nb) != 0) {
        # for 1-n pages query the allele frequencies
        for (page_id in seq(page_nb)) {
            url_tmp <- paste0(url, "&page=", 
                as.character(page_id), collapse = "")
            html_input_tmp <- getURL(url_tmp, read_method = "html")
            rvest_tables_tmp <- rvest::html_table(html_input_tmp, fill = TRUE)
            # add the population ID
            pop_ids <- extract_population_id(html_input_tmp)
            if (type == "allele") {
                freq_table <- parse_allele_freq_html(rvest_tables_tmp[[5]])
            } else if (type == "haplotype") {
                freq_table <- parse_haplotype_freq_html(rvest_tables_tmp[[5]])
            }
            
            freq_table$population_id <- pop_ids
            
            output_table <- rbind(output_table, freq_table)
        }
    }
    # return the final output table
    output_table
}

#
# FUNCTIONS TO READ POPULATION DETAILS
#

#' extract_population_name
#' @param data html input page
#' @return population name
#' @keywords internal
extract_population_name <- function(data) {
    rvest::html_text(
        rvest::html_nodes(data, xpath = "//td//h1"),
        trim = TRUE)
}

#' extract_population_info
#' @param data html input page
#' @return population information
#' @keywords internal
extract_population_info <- function(data) {
    raw_pop_data <- rvest::html_table(
        rvest::html_nodes(data, ".table04:nth-child(4)"),
        fill = TRUE)[[1]]
    # reformat in a key value manner
    key_names <- stringr::str_replace(
        raw_pop_data$X2[3:nrow(raw_pop_data)], ":",
        "")
    values <- raw_pop_data$X3[3:nrow(raw_pop_data)]
    
    names(values) <- key_names
    values
}

#' extract_sample_info
#' @param data html input page
#' @return sample information
#' @keywords internal
extract_sample_info <- function(data) {
    raw_sample_data <- rvest::html_table(
        rvest::html_nodes(data, ".table04:nth-child(6)"),
        fill = TRUE)[[1]][seq(8),]
    # reformat in a key value manner
    key_names <- stringr::str_replace(
        raw_sample_data$X2[3:nrow(raw_sample_data)], ":",
        "")
    values <- raw_sample_data$X3[3:nrow(raw_sample_data)]
    
    names(values) <- key_names
    values
}

#' read_population_detail
#' @param url url of the page containing population information
#' @param population_id numeric population identifier
#' @return population information
#' @keywords internal
read_population_detail <- function(url, population_id) {
    data <- getURL(url, read_method = "html")
    
    # Population
    population_name <- extract_population_name(data)
    
    # Population data
    pop_data <- extract_population_info(data)
    
    # sample data
    sample_data <- extract_sample_info(data)
    
    # Assemble into a data frame
    # one row with all info in columns
    data.frame(population_id = population_id,
        t(pop_data),
        t(sample_data), stringsAsFactors = FALSE)
}

#
# INPUT VALIDATION FUNCTIONS
#

#' check_hla_locus, stops if input not adequate
#' @param hla_locus locus of hla frequencies to query
#' @return TRUE
#' @keywords internal
check_hla_locus <- function(hla_locus) {
    if (is.na(hla_locus)) {TRUE} else {
        # classical HLA loci
        valid_hla_loci <- c("A","B","C",
                            "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
        if (hla_locus %in% valid_hla_loci) {
            TRUE
        } else {
            stop("The hla_locus ", hla_locus, 
                " does not belong to the list of valid HLA loci: ",
                paste(valid_hla_loci, collapse=",")) 
        }}
}

#' check_hla_selection, stops if input not adequate
#' @param hla_selection HLA alleles used for selection
#' @param query_type ["allele"|"haplotype"] type of AFND query
#' @return TRUE
#' @keywords internal
check_hla_selection <- function(hla_selection, query_type) {
    if (any(is.na(hla_selection))) {TRUE} else {
        # we do not check whether the indicated allele is in the database
        # we just check general formatting
        # hla_selection can be a list of alleles
        if (query_type == "allele") {
            # for alleles at least A*01 level is required
            pattern_match <- grepl("^(\\w*|\\w*\\d)\\*\\d\\d", hla_selection)
        } else if (query_type == "haplotype") {
            # for haplotype A* level is ok to query any HLA-A allele
            pattern_match <- grepl("^(\\w*|\\w*\\d)\\*", hla_selection)  
        } else {stop("query_type must be allele or haplotype")}
        
        if (all(pattern_match)) {
            TRUE
        } else {
            wrong_allele <- hla_selection[!pattern_match]
            stop("The following alleles do not seem to be formated correctly ",
            "(expected format e.g. A*01:01): ",
                paste(wrong_allele, collapse=","))
        }}
}

#' check_population, stops if input not adequate
#' @param hla_population population id
#' @return TRUE
#' @keywords internal
check_population <- function(hla_population) {
    if (is.na(hla_population)) {TRUE} else {
        # We just check formatting
        # we do not check whether this population actually exists
        # format seems to be 4 numbers, but it has to be formatted as numeric
        pattern_match <- vapply(hla_population, is.numeric, 
            FUN.VALUE = logical(1))
        if (all(pattern_match)) {
            TRUE
        } else {
            wrong_pop <- hla_population[!pattern_match]
            stop("The following population IDs do not seem to be formated ",
            "correctly (expected format e.g. 1920): ",
                paste(as.character(wrong_pop), collapse=","))
        }}
}


#' get_valid_geographics
#' @description get a list of valid countries, regions and ethnic origins
#' @return list of valid countries, regions and ethnic origin
#' @keywords internal
get_valid_geographics <- function() {
    url <- "http://www.allelefrequencies.net/hla6006a.asp?"
    html_input <- getURL(url, read_method = "html")
    
    rvest_tables <- rvest::html_table(html_input, fill = TRUE)
    
    # country
    selection_str_1 <- rvest_tables[[3]]$X1[[5]]
    split_selection_str_1 <- stringr::str_split(selection_str_1,
        "(\r\n\t\t\t\t\t)|(\r\n)")[[1]]
    
    valid_countries <- stringr::str_split(split_selection_str_1[[4]],
        pattern = stringr::regex("(?<=[a-z]|\\))(?=[A-Z])"))[[1]]
    
    # region, ethnic
    selection_str_2 <- rvest_tables[[3]]$X1[[6]]
    split_selection_str_2 <- stringr::str_split(selection_str_2,
        "(\r\n\t\t\t\t\t)|(\r\n)")[[1]]
    
    valid_regions <- stringr::str_split(split_selection_str_2[[2]],
        pattern = stringr::regex("(?<=[a-z])(?=[A-Z])"))[[1]]
    
    valid_ethnic <- stringr::str_split(split_selection_str_2[[4]], 
        pattern = stringr::regex("(?<=[a-z])(?=[A-Z])"))[[1]]
    
    list(valid_countries = valid_countries,
        valid_regions = valid_regions,
        valid_ethnic = valid_ethnic)
}

valid_geographics <- get_valid_geographics()

#' check_geographics, stops if input not adequate
#' @param country country used for allele frequency selection
#' @param region geographical region used for allele frequency selection
#' @param ethnic ethical origin used for allele frequency selection
#' @return TRUE
#' @keywords internal
check_geographics <- function(country, region, ethnic) {
    # do not test is its na
    if (!is.na(country)) {
        if(!(country %in% valid_geographics$valid_countries)) {
            stop("Country not part of valid countries list: ", country,
                ". Valid list of countries is: ", 
                paste(valid_geographics$valid_countries, collapse = ", "))
        }
    }
    if (!is.na(region)) {
        if(!(region %in% valid_geographics$valid_regions)) {
            stop("Region not part of valid regions list: ", region, 
                ". Valid list of regions is: ", 
                paste(valid_geographics$valid_regions, collapse = ", "))
        }
    }
    if (!is.na(ethnic)) {
        if(!(ethnic %in% valid_geographics$valid_ethnic)) {
            stop("Ethnic origin not part of valid regions list: ", ethnic, 
                ". Valid list of ethnics origin is: ", 
                paste(valid_geographics$valid_ethnic, collapse = ", "))
        }
    }
}

#' check_sample_size, stops if input not adequate
#' @param hla_sample_size_pattern one of "bigger_than", "equal", "less_than", 
#' "less_equal_than", "bigger_equal_than","different"
#' @param hla_sample_size integer number used for population size
#' @return TRUE
#' @keywords internal
check_sample_size <- function(hla_sample_size_pattern, hla_sample_size) {
    if (is.na(hla_sample_size_pattern) & is.na(hla_sample_size)) {TRUE} else {
        valid_pattern <- c("bigger_than",
            "equal", "less_than",
            "less_equal_than", "bigger_equal_than", "different")
        
        if (hla_sample_size_pattern %in% valid_pattern) {
            if (is.numeric(hla_sample_size)) {TRUE}
            else { stop("hla sample size must be numeric") }
        } else {
            stop("The hla_sample_size_pattern must be one of these options: ",
                paste(valid_pattern, collapse=","))
        }}
}

#' check_standard, stops if input not adequate
#' @param standard one of "a" - all,"s" - silver,"g" - gold
#' @return TRUE
#' @keywords internal
check_standard <- function(standard) {
    # all, silver or gold standard of population data quality
    valid_standard <- c("a","s","g")
    if(standard %in% valid_standard) {
        TRUE
    } else { stop("standard must be one of the following: ",
        paste(valid_standard, collapse=",") )}
}


#' verify_parameters
#'
#' @param hla_locus locus of hla frequencies to query
#' @param hla_selection HLA alleles used for selection
#' @param hla_population population id
#' @param hla_country country used for allele frequency selection
#' @param hla_region geographical region used for allele frequency selection
#' @param hla_ethnic ethical origin used for allele frequency selection
#' @param hla_sample_size_pattern one of "bigger_than", "equal", "less_than", 
#' "less_equal_than", "bigger_equal_than","different"
#' @param hla_sample_size integer number used for population size
#' @param standard  one of "a" - all,"s" - silver,"g" - gold
#' @param query_type "allele" or "haplotype"
#'
#' @return boolean to indicate, whether tests were passed
#' @keywords internal
verify_parameters <- function(hla_locus,
                            # for now we don't check the 
                            hla_selection,
                            hla_population,
                            hla_country,
                            hla_region,
                            hla_ethnic,
                            hla_sample_size_pattern,
                            hla_sample_size,
                            standard = "a",
                            query_type) {
    
    # do the checks only if variables are not na
    if (all(is.na(c(hla_locus, hla_selection, hla_population,
        hla_sample_size_pattern, hla_sample_size)))) {
        stop("You need to at least specify one of the paramters hla_locus, ",
        "hla_selection, hla_population, hla_sample_size_pattern, ",
        "hla_sample_size.")
    } else {
        # do the checks only if variables are not na
        all(check_standard(standard),
            check_sample_size(hla_sample_size_pattern, hla_sample_size),
            check_population(hla_population),
            check_geographics(hla_country, hla_region, hla_ethnic),
            check_hla_selection(hla_selection, query_type),
            check_hla_locus(hla_locus))
    }
}
