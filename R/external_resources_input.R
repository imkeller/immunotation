#
# MRO
#

#' load_mro 
#' @return MRO in ontology_index format
#' @keywords internal
load_mro <- function() {
    obo_file_name <- system.file("extdata", "mro.obo.gz",
        package = "immunotation")
    
    mro.obo <- ontologyIndex::get_OBO(obo_file_name,
        extract_tags="everything")
    
    mro.obo
}

#' @title MRO
#' @description   \code{mro.obo}: MHC Restriction Ontology obo file.
#' @keywords internal
mro.obo <- load_mro()


#
#   Reading URLs
#

#' getURL
#'
#' @param URL  Indicated the url that will be read
#' @param N.TRIES   Integer, how often should the function try to read the URL?
#' @param read_method Method to be used for reading of URL content 
#' ("delim" -> \code{readr::read_delim}, 
#' "lines" -> \code{readr::read_lines}, "html" -> \code{xml2::read_html})
#' @param skip integer indicating how many lines to skip when reading URL 
#' @param delim pattern used for delim 
#' (passed to \code{delim} of read functions)
#' @param col_names list of colnames to use
#'
#' @return returns a the content of the URL. The format of the return object
#' depends on the read_method that was used.
#' @keywords internal
getURL <- function(URL, N.TRIES=2L, 
    read_method = c("delim", "lines", "html"),
    skip = 0, delim = "\t", col_names = TRUE) {
    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))
    while (N.TRIES > 0L) {
        if (read_method == "delim") {
            result <- tryCatch(suppressMessages(readr::read_delim(URL,
                delim = delim, skip = skip, col_names = col_names)),
                error=identity)
        } else if (read_method == "lines") {
            result <- tryCatch(readr::read_lines(URL), error=identity)
        } else if (read_method == "html") {
            result <- tryCatch(xml2::read_html(URL), error=identity)
        } else {stop("Read method unknown")}
        
        if (!inherits(result, "error"))
            break
        N.TRIES <- N.TRIES - 1L
    }
    if (N.TRIES == 0L) {
        stop("'getURL()' failed:",
            "\n  URL: ", URL, "\n  error: ", conditionMessage(result))
    }
    result
}

#
#   NetMHC input
#

# MHC I
# netmhcI_input_template is an internal variable containing list of valid 
# NetMHCpan input alleles
netmhcI_input_template <- getURL(
    URL="https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/allele.list",
    read_method = "delim", delim = "\t",
    col_names = c("netmhc_input", "hla_chain_name", "HLA_gene"))

# MHC II
lines <- getURL(
    URL = paste0("https://services.healthtech.dtu.dk/services/",
    "NetMHCIIpan-4.0/alleles_name.list"),
    read_method = "lines")
lines_rep <- stringr::str_replace_all(lines, "\t+|\\s\\s+", "\t")
netmhcII_input_template <- suppressWarnings(
    suppressMessages(read.delim(textConnection(lines_rep), sep = "\t")))


all_netmhcII_template <- c(netmhcII_input_template$DR, 
    netmhcII_input_template$DP.alpha, netmhcII_input_template$DP.beta,
    netmhcII_input_template$DQ.alpha, netmhcII_input_template$DQ.beta)
# all_netmhcII_template is an internal variable containing list of valid 
# NetMHCIIpan input alleles
all_netmhcII_template <- unique(all_netmhcII_template[!is.na(all_netmhcII_template)])


#
# G and P groups
#

#' get_external_file
#'
#' @param file  Indicated the file that will be read
#' @param skip integer indicating how many lines to skip when reading URL 
#' @param delim pattern used for delim 
#' (passed to \code{delim} of read functions)
#' @param col_names list of colnames to use
#'
#' @return returns a the content of the file. The format of the return object
#' depends on the read_method that was used.
#' @keywords internal
get_external_file <- function(file, skip = 0, delim = "\t", col_names = TRUE) {
    result <- tryCatch(suppressMessages(readr::read_delim(file,
            delim = delim, skip = skip, col_names = col_names)),
            error=identity)
    if (!inherits(result, "error")) {
        return(result)
    }
    else {
        stop("'get_external_file()' failed:",
            "\n  file: ", file, "\n  error: ", conditionMessage(result))
    }
}

# G group
# wget https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt
g_group <- get_external_file(file = system.file("extdata", "hla_nom_g.txt",
        package = "immunotation"),
    skip = 6, delim = ";", 
    col_names = c("locus", "g_group", "g_group_name"))
g_group$locus <- stringr::str_extract(g_group$locus,
    pattern = "[:alnum:]+(?=\\*?)")


# P group
# wget https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt
p_group <- get_external_file(file = system.file("extdata", "hla_nom_p.txt",
        package = "immunotation"),
    skip = 6, delim = ";", 
    col_names = c("locus", "p_group", "p_group_name"))
p_group$locus <- stringr::str_extract(p_group$locus,
    pattern = "[:alnum:]+(?=\\*?)")

