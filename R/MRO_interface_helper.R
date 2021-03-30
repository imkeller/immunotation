#
#    Helper functions related to SINGLE CHAINS
#

# build strings for lookup of specific relationships
build_intersection <- function(MRO_ids_queried, type = c("has part", "gene product of")) {
    # find all entries that are "type"
    type_id <- mro.obo$id[mro.obo$name == type]
    # return as string which will be needed for lookup
    str_c(type_id, " ", MRO_ids_queried)
}

# build an intersection of two entry lists
find_intersection_ids <- function(intersection_list, int_list) {
    hits <- sapply(intersection_list, grep, x = int_list)
    if (length(hits[[1]]) == 0) {hits <- NA}
    names(int_list)[hits]
}

# extract the sequence information for a given chain
extract_protein_sequence <- function(chain_list) {
    # build a recursion to extract for all chain ids
    proteins_seqs <- sapply(chain_list,
           function(X) {
               annot <- mro.obo$property_value[mro.obo$id == X][[1]]
               # for some chains no seq is given, this is why we have to test if character(0)
               if(length(annot) > 0) {
                   prop_value <- annot[[3]]
                   return(str_extract(prop_value, "(?<=\\\")[A-Z]*(?=\\\")"))
               } else {return("")}

           }
    )
    proteins_seqs
}


#
#    Helper functions related to PROTEIN COMPLEX
#

# filter the molecules to only get complete molecules
filter_molecules <- function(subcomplex) {
    subcomplex[grep("complete molecule", mro.obo$property_value[subcomplex])]
}

# get all the descendant complexes to a given complex parent
# only return complete molecules
find_descendant_complexes <- function(complex_list) {
    subcomplexes <- sapply(complex_list, function(X) ontologyIndex::get_descendants(roots = X, ontology = mro.obo, exclude_roots = TRUE))
    sapply(subcomplexes, filter_molecules)
}

# find if a molecule complex has complete or partial annotation
get_complex_status <- function(id) {
    prop_value <- ontologyIndex::get_term_property(mro.obo, property_name = "property_value", id, as_names = FALSE)
    # location of MRO:0001984
    loc_comp <- grep("MRO:0001984", prop_value)
    molecule_state <- strsplit(prop_value[loc_comp], split = "\\s\"|\"\\s")[[1]][[2]]
    molecule_state
}

# return chains used by a given complex
get_complex_chains <- function(id) {
    prop_value <- ontologyIndex::get_term_property(mro.obo, property_name = "intersection_of", id, as_names = FALSE)
    # location of MRO:0001984
    loc_comp <- grep("BFO:0000051", prop_value)
    chains <- str_extract(prop_value[loc_comp], pattern = "(?<=\\s)(.*)")
    chains
}

# return serotype used by a given complex
get_complex_serotype <- function(id) {
    prop_value <- ontologyIndex::get_term_property(mro.obo, property_name = "intersection_of", id, as_names = FALSE)
    loc_serotype <- grep("MRO:0000001", prop_value)
    serotype <- str_extract(prop_value[loc_serotype], pattern = "(?<=\\s)(.*)")
    if(length(serotype) == 0) {serotype <- NA}
    serotype
}
