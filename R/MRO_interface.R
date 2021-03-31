# These are the functions used to query information of the MRO
# Helper function that are needed here are defined in MRO_interface_helper.R

#
#    Functions related to ORGANISM
#

#' @title get_valid_organisms
#' @description get the list of organisms that are part of the MRO annotation
#' @return list of organisms
#' @export
#'
#' @examples
#' get_valid_organisms()
#' 
get_valid_organisms <- function() {
    # Organism OBI:0100026
    organism_children <- mro.obo$children$`OBI:0100026`
    mro.obo$name[organism_children]
}

# define a (package internal) variable 
# get a list of all valid organisms 
valid_organisms <- get_valid_organisms()

# Sanity check for the user defined organism input
organism_input <- function(organism) {
    # error when organism wrong
    stopifnot("Organism name invalid" = organism %in% valid_organisms)
    # return the organism id
    names(valid_organisms)[organism == valid_organisms]
}

#
#    Functions related to SINGLE CHAINS
#

# for a given locus (defined by species and mhc-i, mhc-ii....)
# return all children that are proteins
# --> return all protein chains encoded in this locus
extract_mhc_type_table <- function(species_locus, protein_masked_intersections) {
    # get all descendants
    all_species_loci <- ontologyIndex::get_descendants(mro.obo, species_locus)
    # build a list of strings to find the gene products of the above descendants
    intersection_terms <- build_intersection(all_species_loci, type = "gene product of")
    # intersect to find gene product descendants that are proteins
    masked_protein_chains <- find_intersection_ids(intersection_terms, int_list = protein_masked_intersections)
    # get all children of the above defined categories
    all_chains <- mro.obo$children[mro.obo$id %in% masked_protein_chains]
    # assumption that ids are always 11 characters long
    chains <- unlist(all_chains, use.names = FALSE)
    # if there is a result result...
    if (!is.null(chains)) {
        # we assume that the ids are always 11 long
        names(chains) <- substr(names(unlist(all_chains)),1,11)
    }
    return(chains)
}

assemble_chain_lookup_table <- function(organism_id) {
    # get all MHC loci present in this organism
    # MHC locus: MRO:0000004
    # MHC I, MHC II and non-classical MHC
    loci_cat <- mro.obo$children[mro.obo$id == "MRO:0000004"]
    
    # for these 3 loci, get the id of the ones belonging to species of interest
    all_loci <- mro.obo$children[mro.obo$id %in% loci_cat[[1]]]
    
    # all entries belonging to taxon of interest
    # using "in taxon" relationship `RO:0002162`
    entries_for_species <- mro.obo$id[mro.obo$`RO:0002162` == organism_id]
    
    # for every locus select the entry belonging to the species of interest
    species_loci <- sapply(all_loci, intersect, y = entries_for_species, simplify = FALSE)
    ### here eliminate if empty loci exist
    species_loci <- species_loci[!sapply(species_loci, identical, y = character(0))]
    
    # protein masked ontology (to make it possible to intersect later)
    # only entries classified as protein PR:0000001
    protein_masked_intersections <- mro.obo$intersection_of[grep("PR:000000001", mro.obo$intersection_of)]
    protein_mask <- names(protein_masked_intersections)
    
    # for all species_loci, get the table of chains
    tables <- sapply(species_loci, extract_mhc_type_table, 
                     protein_masked_intersections = protein_masked_intersections, simplify = FALSE)
    
    # assemble the results into a dataframe
    unlisted_tables <- unlist(tables)
    names_chains <- names(unlisted_tables)
    chain_table <- data.frame(mhc_type_id = stringr::str_extract(names_chains, pattern = ".+(?=\\.)"),
                              locus_id = stringr::str_extract(names_chains, pattern = "(?<=\\.).+"),
                              chain_id = unlisted_tables)
    
    # filter out above categories such as HLA-DRB1 chain that is child of HLA-DRB chain
    chain_table_filtered <- chain_table[!chain_table$chain_id %in% protein_mask,]
    
    # get names of mhc type
    mhc_type_names <-  mro.obo$name[mro.obo$id %in% unique(chain_table_filtered$mhc_type_id)]
    chain_table_filtered$mhc_type_names <- mhc_type_names[chain_table_filtered$mhc_type_id]
    
    # get names of loci
    locus_names <-  mro.obo$name[mro.obo$id %in% unique(chain_table_filtered$locus_id)]
    chain_table_filtered$locus_names <- locus_names[chain_table_filtered$locus_id]
    
    # get the names of chains
    chain_names <- mro.obo$name[mro.obo$id %in% chain_table_filtered$chain_id]
    chain_names_short <-  stringr::str_extract(chain_names, pattern = ".*(?= chain)")
    names(chain_names_short) <- names(chain_names)
    chain_table_filtered$chain_names <- chain_names_short[chain_table_filtered$chain_id]
    
    # get the chain sequence
    chain_sequences <- extract_protein_sequence(chain_table_filtered$chain_id)
    chain_table_filtered$chain_sequences <- chain_sequences
    
    return(chain_table_filtered)
}

#' @title Retrieve MHC chain lookup table
#'
#' @param organism name of organism (e.g. "human")
#'
#' @return Table containing MHC chain information for the organism. 
#' It contains chain names, MHC restriction and protein sequence.
#' @export
#'
#' @examples
#' retrieve_chain_lookup_table("mouse")
#' 
retrieve_chain_lookup_table <- function(organism) {
    
    # assemble the chain table
    assemble_chain_lookup_table(organism_id = "mouse")
}

#
#    Functions related to PROTEIN COMPLEX
#

# build the table containing mhc_type, complex name, alpha, beta name etc...
# this is done individually for every mhc_type
# complex_sublist_list is a list containing the list of complexes for evey mhc type (list of lists)
extract_mhc_complex_table <- function(complex_sublist_list, mhc_type) {
    # Need to access the name of the list entry
    complex_type <- names(complex_sublist_list)
    complex_sublist <- complex_sublist_list[[complex_type]]
    
    # check if list contains something
    if (length(complex_sublist) != 0) {
        # complex serotype
        complex_serotypes <- sapply(complex_sublist, get_complex_serotype)
        serotypes_names_dict <- mro.obo$name[mro.obo$id %in% unique(complex_serotypes)]
        serotypes_names <- serotypes_names_dict[complex_serotypes]
        # name of protein complex
        complex_names <- mro.obo$name[mro.obo$id %in% complex_sublist]
        # status complete or partial
        complex_status <- sapply(complex_sublist, get_complex_status)
        # chains that compose the locus
        complex_chains <- sapply(complex_sublist, get_complex_chains, simplify = FALSE)
        # alpha chain
        alpha_chains <- sapply(complex_chains, function(x) x[1])
        # dict is used because the chains might be present in the list several times
        alpha_names_dict <- mro.obo$name[mro.obo$id %in% alpha_chains]
        alpha_names <- alpha_names_dict[alpha_chains]
        # beta chain
        beta_chains <- sapply(complex_chains, function(x) x[2])
        
        beta_names_dict <- mro.obo$name[mro.obo$id %in% beta_chains]
        beta_names <- beta_names_dict[beta_chains]
        # build the table 
        complex_table <- data.frame(mhc_type = rep(mhc_type, length(complex_names)),
                                    complex_type = rep(complex_type, length(complex_names)),
                                    complex_name = complex_names,
                                    complex_status = complex_status,
                                    alpha_name = alpha_names,
                                    beta_name = beta_names,
                                    complex_serotype = complex_serotypes,
                                    serotype_name = serotypes_names,
                                    row.names = NULL)
        
        
        return(complex_table)
    } else {
        return(NULL)
    }
}



#' @title Assemble protein complex
#' @description Assemble a table or MHC protein complexes for a given organism.
#' 
#' @param organism_id  Organism for which the lookup should be built (e.g. "human", "mouse", ...). The list of valid 
#' organisms can be found using the function \code{get_valid_organisms}
#'
#' @return a data frame with the MHC complexes annotated in MRO (only completely annotated complexes are returned)
#' @export
#'
#' @examples
#' org_id <- organism_input("mouse")
#' assemble_protein_complex(org_id)
#' 
assemble_protein_complex <- function(organism_id) {
    # MHC class I protein complex: MRO:0001355
    # MHC class II protein complex: MRO:0001356
    # nonclassical MHC protein complex MRO:0001464
    # get all protein complexes belonging to these entries and overlapping with taxon queried
    complexes_mhc_level <- mro.obo$children[mro.obo$id %in% c("MRO:0001355", "MRO:0001356", "MRO:0001464")]
    
    # filter for taxon of interest
    entries_for_species <- mro.obo$id[grep(organism_id, mro.obo$intersection_of)]
    species_complexes_mhc_level <- sapply(complexes_mhc_level, intersect, y = entries_for_species, simplify = FALSE)
    # find all children - locus level
    complexes_locus_level <- mro.obo$children[mro.obo$id %in% species_complexes_mhc_level]
    
    # find the children with entries partial/complete molecule
    ## for every complex locus level find all partial and complete descendant complexes
    complex_list <- sapply(complexes_locus_level, find_descendant_complexes)
    
    # create an empty dataframe for assembling info to come
    complete_table <- data.frame(mhc_type = character(0), complex_type = character(0), complex_name = character(0), 
                                 complex_status = character(0), alpha_name = character(0),  beta_name = character(0), 
                                 complex_serotype= character(0), serotype_name= character(0))
    # loop over mhc class
    for (mhc_type in names(complex_list)) {
        complex_sublist <- complex_list[[mhc_type]]
        
        # complex types
        complex_types <- names(complex_sublist)
        sub_tables <- sapply(complex_types, function(x) extract_mhc_complex_table(complex_sublist[x],
                                                                                  mhc_type = mhc_type), simplify = FALSE)
        int_table <- do.call(rbind, sub_tables)
        complete_table <- rbind(complete_table, int_table)
    }
    
    return(complete_table)
}



# The protein complex table for humans can be assembled before, because it takes quite some time
# get the protein lookup table
#' human_protein_complex_table
#' @details   \code{human_protein_complex_table}: human_protein_complex_table.
#' @export
human_protein_complex_table <- assemble_protein_complex(organism_id = organism_input("human"))


#' @title Serotypes
#' @description Get the serotypes of the MHC complexes encoded by a list of MHC alleles.
#'
#' @param allele_list List of allele
#' @param organism  Organism to be used for MRO lookup. 
#' If the organism does not match the given allele, a empty object is returned.
#' @param mhc_type  ["MHC-I" or "MHC-II"] MHC class to use for MRO lookup.
#'
#' @return Named list of serotypes, which only contains complexes contained in the MRO. 
#' If no serotype is annoted for a given complex, the list element is NA.
#' @export
#'
#' @examples
#' allele_list <- c("A*01:01:01","B*27:01")
#' get_serotypes(allele_list, mhc_type = "MHC-I")
#' 
get_serotypes <- function(allele_list, organism = "human", mhc_type) {
    chain_list <- sapply(allele_list, reformat_allele)
    
    # get the protein lookup table
    if (organism == "human") {
        protein_complexes <- human_protein_complex_table
    } else {
        orgid <- organism_input(organism)
        protein_complexes <- assemble_protein_complex(organism_id = orgid)
    }
    
    # filter the matches to alpha name
    full_chain_list <- sapply(chain_list, function(X) paste0("HLA-", X, " chain"))
    if (mhc_type == "MHC-I") {
        filtered_complexes <- protein_complexes[protein_complexes$alpha_name %in% full_chain_list, ]
    } else if (mhc_type == "MHC-II") {
        filtered_complexes <- protein_complexes[protein_complexes$alpha_name %in% full_chain_list & 
                                                    protein_complexes$beta_name %in% full_chain_list, ]
    } else {stop("mhc_type needs to be MHC-I or MHC-II")}
    
    named_serotype_list <- filtered_complexes$serotype_name
    names(named_serotype_list) <- filtered_complexes$complex_name
    named_serotype_list 
}

