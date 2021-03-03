# naming format
# Genetic
# 
# HLA-A0101
# HLA-A01:01 ok
# HLA-A*01:01 ok
# HLA-A*01:01:01 ok
# HLA-A01 serotype 
# HLA00001 HLA ID
# protein group
# 

#Task 1: 
# arcasHLA presumably mentions gene groups
# convert to MHCpan input
#A*01:01:01 -> A*01:01
#A*03:01:01 -> A*03:01

#' get_mhcpan_input
#'
#' @param name_list 
#'
#' @return Protein chain list as formatted for MHCpan input
#' @export
#'
#' @examples
#' allele_list <- c("A*01:01:01","B*27:01")
#' get_mhcpan_input(allele_list, mhc_class = "MHC-I")
get_mhcpan_input <- function(name_list, mhc_class = c("MHC-I", "MHC-II")) {
    # convert to locus*dd:dd format
    protein_chain_names <- reformat_allele(name_list)
    
    # convert to MHCpan naming scheme
    if(mhc_class == "MHC-I") {
        # convert to MHCpanI format: HLA-A01:01
        # since there is only one chain this is easy
        mhc_name_list <- sapply(protein_chain_names, convert_naming, add_hla = TRUE, add_star = FALSE)
        # check that they are allowed in the MHC input
        valid_in_reference <- sapply(mhc_name_list, function(x) x %in% netmhcI_input_template$netmhc_input)
        if (!all(valid_in_reference)) {warning(paste(mhc_name_list[!valid_in_reference], "not in NetMHCpan input list.\n"))}
        mhc_name_list <- mhc_name_list[valid_in_reference]
    } else if (mhc_class == "MHC-II") {
        # check which chains belong together
        # convert to NetMHCIIpan format: HLA-DPA10103-DPB10101
        mhc_name_list <- build_mhcII_complexes(protein_chain_names)
    } else {stop("mhc_class schould be MHC-I or MHC-II.")}
    mhc_name_list
}

build_mhcII_complexes <- function(protein_chain_names) {
    # check what is present
    valid_in_reference <- sapply(protein_chain_names, function(x) x %in% all_netmhcII_template)
    # check if there are DRA alleles (they are not relevant in the netmhciipan list)
    valid_in_reference <- valid_in_reference | grepl("DRA", protein_chain_names)
    
    if(!all(valid_in_reference)) {
        warning(paste(protein_chain_names[!valid_in_reference], "not in NetMHCIIpan input list.\n"))
    }

    protein_chain_names_valid <- protein_chain_names[valid_in_reference]       
    
    # DPA1 DPB1
    dp_complexes <- assemble_dp_dq(protein_chain_names_valid, type = "DP")
    dp_complexes <- str_replace(dp_complexes, " ", "-")
    # DQA1 DQB1
    dq_complexes <- assemble_dp_dq(protein_chain_names_valid, type = "DQ")
    dq_complexes <- str_replace(dq_complexes, " ", "-")
    # DRA DRB
    drb <- protein_chain_names_valid[grep("DRB", protein_chain_names_valid)]
    drb <- sapply(drb, convert_naming, add_hla = FALSE, add_star = TRUE, remove_colon = TRUE, star = "_")
    c(drb, dq_complexes, dp_complexes)
}

assemble_dp_dq <- function(protein_chain_names, type = c("DP", "DQ")) {
    if (type == "DP") {
        alpha <- "DPA1"; beta <- "DPB1"
    } else if (type == "DQ") {
        alpha <- "DQA1"; beta <- "DQB1"
    } else {stop("type needs to be DP or DQ")}
    
    alpha_names <- protein_chain_names[grep(alpha, protein_chain_names)]
    alpha_names <- sapply(alpha_names, convert_naming, add_hla = TRUE, add_star = FALSE, remove_colon = TRUE)
    beta_names <- protein_chain_names[grep(beta, protein_chain_names)]
    beta_names <- sapply(beta_names, convert_naming, add_hla = FALSE, add_star = FALSE, remove_colon = TRUE)
    dpb_complexes <- unlist(as.list(outer(alpha_names, beta_names, paste)))
    
    dpb_complexes
}

convert_naming <- function(name, add_hla = TRUE, add_star = TRUE, remove_colon = FALSE, star = "*") {
    hla_prefix <- ""
    if (add_hla) {hla_prefix <- "HLA-"}
    if (!add_star) {star <- ""}
    locus_allele <- str_split(name, pattern = "\\*")[[1]]
    locus <- locus_allele[[1]]
    number <- locus_allele[[2]]
    if (remove_colon) {number <- str_replace(number, ":", "")}
    paste0(hla_prefix, locus, star, number)
}

simplify_allele <- function(name_list, from_format = NA, to_format = NA) {
    gene_group_list <- sapply(name_list, derive_gene_group)
    # match to P group
    # manually extract the right formatting
    protein_chain_list <- sapply(gene_group_list, reformat_allele)
    protein_chain_list
}

reformat_allele <- function(allele_name) {
    simplified_allele <- str_extract(allele_name, pattern = "^[:alnum:]+\\*\\d+:\\d+")
    simplified_allele
}

derive_gene_group <- function(allele_name) {
    locus_allele <- str_split(allele_name, pattern = "\\*")[[1]]
    locus <- locus_allele[[1]]
    allele <- locus_allele[[2]]
    
    locus_g_group <- g_group[g_group$locus == locus,]
    # check if can be found in database
    reg_expression_group <- paste0("/",allele,"/|^",allele,"/|/",allele,"$")
    match_position <- grep(reg_expression_group,locus_g_group$g_group)
    
    # if not found, search in gene group
    if (length(match_position) == 0) {
        reg_expression_name <- paste0("^",allele,"G")
        match_position <- grep(reg_expression_name,locus_g_group$g_group_name)
    }
    g_group_name <- unique(locus_g_group$g_group_name[match_position])
    paste0(locus,"*",g_group_name)
}

derive_protein_group <- function(allele_name) {
    locus_allele <- str_split(allele_name, pattern = "\\*")[[1]]
    locus <- locus_allele[[1]]
    allele <- locus_allele[[2]]
    
    locus_p_group <- p_group[p_group$locus == locus,]
    # check if can be found in database
    reg_expression_group <- paste0("/",allele,"/|^",allele,"/|/",allele,"$")
    match_position <- grep(reg_expression_group,locus_p_group$p_group)
    
    # if not found, search in gene group
    if (length(match_position) == 0) {
        reg_expression_name <- paste0("^",allele,"P")
        match_position <- grep(reg_expression_name,locus_p_group$p_group_name)
    }
    g_group_name <- unique(locus_g_group$g_group_name[match_position])
    g_group_name
}

# The protein complex table for humans can be assembled before, because it takes quite some time
# get the protein lookup table
#' mro.obo
#' @details   \code{human_protein_complex_table}: human_protein_complex_table.
#' @import ontologyIndex
#' @export
human_protein_complex_table <- assemble_protein_complex(organism_id = organism_input("human"))

# Serotypes
#' get_serotypes
#'
#' @param allele_list 
#'
#' @return serotype
#' @export
#'
#' @examples
#' allele_list <- c("A*01:01:01","B*27:01")
#' get_serotypes(allele_list, mhc_type = "MHC-I")
#' 
get_serotypes <- function(allele_list, organism = "human", mhc_type = c("MHC-I", "MHC-II")) {
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
    names(named_serotype_list) <- filtered_complexes$alpha_name
    serotypes <- named_serotype_list[full_chain_list]
    serotypes
}

# MHC II
#protein_complexes[grep("HLA-DRB", protein_complexes$beta_name),]
