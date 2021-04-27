#' @title Get format for NetMHCpan tools
#' @description NetMHCpan tools for MHC-peptide binding prediction require 
#' HLA complex names in a specific format. \code{get_mhcpan_input} formats a 
#' list of HLA alleles into a list of NetMHC-formated complexes.
#'
#' @param allele_list list of HLA alles (e.g. c("A*01:01:01","B*27:01"))
#' @param mhc_class ["MHC-I"|"MHC-II"] indicated which NetMHC you want to use.
#'
#' @return protein chain list as formatted for MHCpan input
#' @export
#'
#' @examples
#' allele_list <- c("A*01:01:01","B*27:01")
#' get_mhcpan_input(allele_list, mhc_class = "MHC-I")
get_mhcpan_input <- function(allele_list, mhc_class) {
    # convert to locus*dd:dd format
    protein_chain_names <- reformat_allele(allele_list)
    
    # convert to MHCpan naming scheme
    if(mhc_class == "MHC-I") {
        # convert to MHCpanI format: HLA-A01:01
        # since there is only one chain this is easy
        mhc_name_list <- vapply(protein_chain_names, convert_naming,
            add_hla = TRUE, add_star = FALSE,
            FUN.VALUE = character(1))
        # check that they are allowed in the MHC input
        valid_in_reference <- vapply(mhc_name_list, function(x) x %in%
                netmhcI_input_template$netmhc_input, FUN.VALUE = logical(1))
        if (!all(valid_in_reference)) {stop(mhc_name_list[!valid_in_reference], 
            " not in NetMHCpan input list.\n")}
        mhc_name_list <- mhc_name_list[valid_in_reference]
    } else if (mhc_class == "MHC-II") {
        # check which chains belong together
        # convert to NetMHCIIpan format: HLA-DPA10103-DPB10101
        mhc_name_list <- build_mhcII_complexes(protein_chain_names)
    } else {stop("mhc_class schould be MHC-I or MHC-II.")}
    # names should not appear, because the complexes for mhc-ii are 
    # composed of different elements of the original list
    names(mhc_name_list) <- NULL
    mhc_name_list
}

build_mhcII_complexes <- function(protein_chain_names) {
    # check what is present
    valid_in_reference <- vapply(protein_chain_names, function(x) x %in% 
            all_netmhcII_template, FUN.VALUE = logical(1))
    # check if there are DRA alleles 
    # (they are not relevant in the netmhciipan list)
    valid_in_reference <- valid_in_reference | grepl("DRA", protein_chain_names)
    
    if(!all(valid_in_reference)) {
        stop(protein_chain_names[!valid_in_reference],
            " not in NetMHCIIpan input list.\n")
    }

    protein_chain_names_valid <- protein_chain_names[valid_in_reference]       
    
    # DPA1 DPB1
    dp_complexes <- assemble_dp_dq(protein_chain_names_valid, type = "DP")
    dp_complexes <- stringr::str_replace(dp_complexes, " ", "-")
    # DQA1 DQB1
    dq_complexes <- assemble_dp_dq(protein_chain_names_valid, type = "DQ")
    dq_complexes <- stringr::str_replace(dq_complexes, " ", "-")
    # DRA DRB
    drb <- protein_chain_names_valid[grep("DRB", protein_chain_names_valid)]
    drb <- vapply(drb, convert_naming, add_hla = FALSE,
        add_star = TRUE, remove_colon = TRUE, star = "_",
        FUN.VALUE = character(1))
    c(drb, dq_complexes, dp_complexes)
}

assemble_dp_dq <- function(protein_chain_names, type = c("DP", "DQ")) {
    if (type == "DP") {
        alpha <- "DPA1"; beta <- "DPB1"
    } else if (type == "DQ") {
        alpha <- "DQA1"; beta <- "DQB1"
    } else {stop("type needs to be DP or DQ")}
    
    alpha_names <- protein_chain_names[grep(alpha, protein_chain_names)]
    alpha_names <- vapply(alpha_names, convert_naming, add_hla = TRUE,
        add_star = FALSE, remove_colon = TRUE,
        FUN.VALUE = character(1))
    beta_names <- protein_chain_names[grep(beta, protein_chain_names)]
    beta_names <- vapply(beta_names, convert_naming, add_hla = FALSE,
        add_star = FALSE, remove_colon = TRUE,
        FUN.VALUE = character(1))
    dpb_complexes <- unlist(as.list(outer(alpha_names, beta_names, paste)))
    
    dpb_complexes
}

convert_naming <- function(name, add_hla = TRUE, add_star = TRUE,
            remove_colon = FALSE, star = "*") {
    hla_prefix <- ""
    if (add_hla) {hla_prefix <- "HLA-"}
    if (!add_star) {star <- ""}
    locus_allele <- stringr::str_split(name, pattern = "\\*")[[1]]
    locus <- locus_allele[[1]]
    number <- locus_allele[[2]]
    if (remove_colon) {number <- stringr::str_replace(number, ":", "")}
    paste0(hla_prefix, locus, star, number)
}

simplify_allele <- function(name_list, from_format = NA, to_format = NA) {
    gene_group_list <- vapply(name_list, derive_gene_group,
        FUN.VALUE = character(1))
    # match to P group
    # manually extract the right formatting
    protein_chain_list <- vapply(gene_group_list, reformat_allele,
        FUN.VALUE = character(1))
    protein_chain_list
}

reformat_allele <- function(allele_name) {
    simplified_allele <- stringr::str_extract(allele_name, 
                                pattern = "^[:alnum:]+\\*\\d+:\\d+")
    if (any(is.na(simplified_allele))) {
        # check if maybe HLA formatting is used
        simplified_allele <- stringr::str_extract(allele_name, 
                                pattern = "(?<=HLA-)[:alnum:]+\\*\\d+:\\d+")
        if (any(is.na(simplified_allele))) { 
            stop("Input allele does not have the right formating: ", 
                    allele_name[is.na(simplified_allele)]) }
    }
    simplified_allele
}

allele_name_sanity_check <- function(allele_name) {
    allele_name_old <- allele_name
    hla_string <- grepl("HLA-", allele_name)
    if (hla_string) {
        # remove the HLA prefix
        allele_name <- stringr::str_replace(allele_name, "HLA-", "")
    }
    # check if * : pattern contained
    valid <- grepl("^.+\\*\\d+:\\d+", allele_name)
    if (!valid) {stop("Input allele does not have the right formating: ",
        allele_name_old)}
    allele_name 
}

derive_gene_group <- function(allele_name) {
    allele_name <- allele_name_sanity_check(allele_name)
    locus_allele <- stringr::str_split(allele_name, pattern = "\\*")[[1]]
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
    # if nothing found, keep original allele name
    if (length(match_position) == 0) {
        g_group_return <- allele_name
    } else {
        g_group_name <- unique(locus_g_group$g_group_name[match_position])
        g_group_return  <- paste0(locus,"*",g_group_name)      
    }
    g_group_return 
}

#' @title G groups
#' 
#' @description Get the G groups for a list of HLA alleles. 
#' [G groups](http://hla.alleles.org/alleles/g_groups.html) are groups of 
#' HLA alleles that have identical nucleotide sequences across the exons 
#' encoding the peptide binding domains.
#'
#' @param allele_list List of alleles.
#'
#' @return Named list of G-groups the input alleles belong to.
#' @export
#'
#' @examples
#' allele_list <- c("DQB1*02:02:01", "DQB1*06:09:01")
#' get_G_group(allele_list)
#' 
get_G_group <- function(allele_list) {
    vapply(allele_list, derive_gene_group,
        FUN.VALUE = character(1))
}

derive_protein_group <- function(allele_name) {
    allele_name <- allele_name_sanity_check(allele_name)
    locus_allele <- stringr::str_split(allele_name, pattern = "\\*")[[1]]
    locus <- locus_allele[[1]]
    allele <- locus_allele[[2]]
    
    locus_p_group <- p_group[p_group$locus == locus,]
    # check if can be found in database
    reg_expression_group <- paste0("/",allele,"/|^",allele,"/|/",allele,"$")
    # if the allele is given at a xx:xx:xx level, also allow for further 
    # matching to /01:03:01:xx/|^01:03:01:xx/|/01:03:01:xx
    # this is necessary because a lot of alleles in p group file are given at 
    # the xx:xx:xx:xx level
    if(stringr::str_count(allele_name, ":") > 1 ) {
        add_regexp_group <- paste0("|/",allele,":|^",allele,":")
        reg_expression_group <- paste0(reg_expression_group, add_regexp_group)
    }
    match_position <- grep(reg_expression_group,locus_p_group$p_group)
    
    # if not found, search in gene group
    if (length(match_position) == 0) {
        reg_expression_name <- paste0("^",allele,"P")
        match_position <- grep(reg_expression_name,locus_p_group$p_group_name)
    }
    # if nothing found, keep original allele name
    if (length(match_position) == 0) {
        p_group_return <- allele_name
    } else {
        p_group_name <- unique(locus_p_group$p_group_name[match_position])
        p_group_return  <- paste0(locus,"*",p_group_name)      
    }
    p_group_return
}

# get p group elements
get_p_group_members <- function(p_group_name) {
    p_group_name <- allele_name_sanity_check(p_group_name)
    locus_allele <- stringr::str_split(p_group_name, pattern = "\\*")[[1]]
    locus <- locus_allele[[1]]
    allele <- locus_allele[[2]]
    p_group_selection <- p_group[p_group$locus == locus &
            !is.na(p_group$p_group_name),]
    return_list <- p_group_selection$p_group[
        p_group_selection$p_group_name == allele]
    if(length(return_list) == 0) {return(p_group_name)} else {
        return_list <- stringr::str_split(return_list, pattern = "/")[[1]]
        return_list <- stringr::str_c(locus,"*",return_list)
        return(return_list)}
}

#' @title P groups
#' 
#' @description Get the P groups for a list of HLA alleles. 
#' [P groups](http://hla.alleles.org/alleles/p_groups.html) are groups 
#' of HLA alleles that have identical protein sequences in the peptide binding 
#' domains.
#'
#' @param allele_list list of HLA alleles
#'
#' @return Named list of P-groups the input alleles belong to.
#' @export
#'
#' @examples
#' allele_list <- c("DQB1*02:02:01", "DQB1*06:09:01")
#' get_P_group(allele_list)
#' 
get_P_group <- function(allele_list) {
    vapply(allele_list, derive_protein_group,
        FUN.VALUE = character(1))
}
