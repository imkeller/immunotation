# get a list of all valid organisms
# Organism OBI:0100026
#' organism_input
#'
#' @param organism 
#'
#' @return
#' @export
#'
#' @examples
#' organism_input("duck")
#' 
organism_input <- function(organism) {
    # error when organism wrong
    stopifnot("Organism name invalid" = organism %in% valid_organisms)
    # return the organism id
    names(valid_organisms)[organism == valid_organisms]
}


# Format the allele list
format_allele_input <- function(allele_list, input_format, organism) {
    if(organism == "human" & input_format == "allele") {
        # for human alleles, extract the first letters + 2 numbers
        chain_level <- str_extract(allele_list, "\\w*\\*\\d*\\:\\d*")
        # add HLA prefix and chain suffix, since needed for querying ontology
        str_c("HLA-", chain_level, " chain")
    }
}



find_MRO_ids <- function(MRO_allele_list) {
    mro.obo$id[mro.obo$name %in% MRO_allele_list]
}

build_intersection <- function(MRO_ids_queried, type = c("has part", "gene product of")) {
    type_id <- mro.obo$id[mro.obo$name == type]
    str_c(type_id, " ", MRO_ids_queried)
}

retrieve_MHC_proteincomplex <- function(alleles = allele_list,
                                        input_format = "allele",
                                        organism = "human") {
    # check if the organism name is valid
    org_id <- organism_input(organism)
    
    # check if the allele list is correctly formatted, get the format needed in MRO
    MRO_allele_list <- format_allele_input(alleles, input_format, organism)
    
    # find the ontology entries corresponding to allele list
    MRO_ids_queried <- find_MRO_ids(MRO_allele_list)
    
    # Get the corresponding entries from the lookup table
    lookup <- assemble_lookup(organism_id = org_id,
                              level = "chain level")
    
    lookup_query <- lookup[lookup$chain_id %in% MRO_ids_queried, ]
    
    # for every locus, build all possible complexes that are fully covered
    
    
    # find the protein complexes containing these chains
    
    # build the intersection of term that indicates part of
    intersection_terms <- build_intersection(MRO_ids_queried)
    
    # this does not work properly
    mro.obo$name[mro.obo$intersection_of %in% intersection_terms]  
    
    # find all complexes, where chain is present
    sapply(intersection_terms, grep, x = mro.obo$intersection_of)
    
}


#' retrieve_chain_lookup_table
#'
#' @param organism 
#'
#' @return
#' @export
#'
#' @examples
#' retrieve_chain_lookup_table("mouse")
#' 
retrieve_chain_lookup_table <- function(organism) {
    # check if the organism name is valid
    org_id <- organism_input(organism)
    assemble_lookup(organism_id = org_id,
                    level = "chain level")
}


find_intersection_ids <- function(intersection_list, int_list) {
    hits <- sapply(intersection_list, grep, x = int_list)
    if (length(hits[[1]]) == 0) {hits <- NA}
    names(int_list)[hits]
}

find_children <- function(node_list) {
    mro.obo$children[mro.obo$id %in% node_list]
}

#species_locus <- "MRO:0000111"
extract_mhc_type_table <- function(species_locus, protein_masked_intersections) {
    all_species_loci <- get_descendants(mro.obo, species_locus)
    intersection_terms <- build_intersection(all_species_loci, type = "gene product of")
    masked_protein_chains <- find_intersection_ids(intersection_terms, int_list = protein_masked_intersections)
    all_chains <- find_children(masked_protein_chains)
    # assumption that ids are always 11 characters long
    chains <- unlist(all_chains, use.names = FALSE)
    if (!is.null(chains)) {
        names(chains) <- substr(names(unlist(all_chains)),1,11)
    }
    return(chains)
}

assemble_lookup <- function(organism_id,
                            level = "chain level") {
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
    
    # protein masked ontology 
    # only entries classified as protein PR:0000001
    protein_masked_intersections <- mro.obo$intersection_of[grep("PR:000000001", mro.obo$intersection_of)]
    
    tables <- sapply(species_loci, extract_mhc_type_table, 
                     protein_masked_intersections = protein_masked_intersections, simplify = FALSE)
    
    unlisted_tables <- unlist(tables)
    names_chains <- names(unlisted_tables)
    chain_table <- data.frame(mhc_type_id = str_extract(names_chains, pattern = ".+(?=\\.)"),
                              locus_id = str_extract(names_chains, pattern = "(?<=\\.).+"),
                              chain_id = unlisted_tables)
    
    protein_mask <- names(protein_masked_intersections)
    # filter out above categories such as HLA-DRB1 chain that is child of HLA-DRB chain
    chain_table_filtered <- chain_table[!chain_table$chain_id %in% protein_mask,]
    
    # get names of mhc type
    mhc_type_names <-  mro.obo$name[mro.obo$id %in% unique(chain_table_filtered$mhc_type_id)]
    chain_table_filtered$mhc_type_names <- mhc_type_names[chain_table_filtered$mhc_type_id]
    
    # get names of locus
    locus_names <-  mro.obo$name[mro.obo$id %in% unique(chain_table_filtered$locus_id)]
    chain_table_filtered$locus_names <- locus_names[chain_table_filtered$locus_id]
    
    # get the names of chains
    chain_names <- mro.obo$name[mro.obo$id %in% chain_table_filtered$chain_id]
    chain_names_short <-  str_extract(chain_names, pattern = ".*(?= chain)")
    names(chain_names_short) <- names(chain_names)
    chain_table_filtered$chain_names <- chain_names_short[chain_table_filtered$chain_id]
    
    return(chain_table_filtered)
}


filter_molecules <- function(subcomplex) {
    subcomplex[grep("complete molecule", mro.obo$property_value[subcomplex])]
}

#complex_list <- complexes_locus_level[[1]]
find_descendant_complexes <- function(complex_list) {
    subcomplexes <- sapply(complex_list, function(X) get_descendants(roots = X, ontology = mro.obo, exclude_roots = TRUE))
    sapply(subcomplexes, filter_molecules )
}

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



#' assemble_protein_complex
#'
#' @param organism_id 
#'
#' @return
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
    #
    #eliminate those that are character(0)
    complete_table <- data.frame(mhc_type = 0, complex_type = 0, complex_name = 0, complex_status = 0, alpha_name = 0,  beta_name = 0, complex_serotype=0, serotype_name=0)
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
#Final table should contain: MHC type, locus, protein complex, chains


