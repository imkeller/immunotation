#' get_valid_organisms
#' query the possible species
#' @return a list of organisms that are part of the annotation functions
#' @export
#'
#' @examples
get_valid_organisms <- function() {
    mro.obo$name[organism_children]
}
valid_organisms <- get_valid_organisms()
# get a list of terms
get_name <- function(term) {
    get_term_property(mro.obo, property_name = "name", term, as_names = FALSE)
}

get_name_list <- function(term_list) {
    sapply(term_list, get_name)
}

get_attr <- function(term) {
    get_term_property(mro.obo, property_name = "property_value", term, as_names = FALSE)
}

get_complex_status <- function(id) {
    prop_value <- get_term_property(mro.obo, property_name = "property_value", id, as_names = FALSE)
    # location of MRO:0001984
    loc_comp <- grep("MRO:0001984", prop_value)
    molecule_state <- strsplit(prop_value[loc_comp], split = "\\s\"|\"\\s")[[1]][[2]]
    molecule_state
}

get_complex_chains <- function(id) {
    prop_value <- get_term_property(mro.obo, property_name = "intersection_of", id, as_names = FALSE)
    # location of MRO:0001984
    loc_comp <- grep("BFO:0000051", prop_value)
    chains <- str_extract(prop_value[loc_comp], pattern = "(?<=\\s)(.*)")
    chains
}

get_complex_serotype <- function(id) {
    prop_value <- get_term_property(mro.obo, property_name = "intersection_of", id, as_names = FALSE)
    loc_serotype <- grep("MRO:0000001", prop_value)
    serotype <- str_extract(prop_value[loc_serotype], pattern = "(?<=\\s)(.*)")
    if(length(serotype) == 0) {serotype <- NA}
    serotype
}
