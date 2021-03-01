
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

# extract annotation foe a protein chain

prop_value <- get_term_property(mro.obo, property_name = "property_value", "MRO:0011061", as_names = FALSE)
entries <- strsplit(prop_value[1], split = "\\s\"|\"\\s")

proplist1 <- get_term_property(mro.obo, property_name = "intersection_of", "MRO:0001150", as_names = FALSE)
# location of MRO:0001984
n_comp <- grep("MRO:0001984", proplist1)
entries <- strsplit(proplist1[n_comp], split = "\\s\"|\"\\s")[[1]][[2]]
