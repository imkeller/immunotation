#### Data ####

load_mro <- function() {
    obo_file_name <- system.file("extdata", "mro_mod.obo", package = "immunotation")
    
    mro.obo <- ontologyIndex::get_OBO(obo_file_name, extract_tags="everything")
    
    mro.obo
}

#' mro.obo
#' @details   \code{mro.obo}: MHC Restriction Ontology obo file.
#' @export
mro.obo <- load_mro()


# NetMHC list of valid chains

# read file from the web
# special function that will be formated according to BC recommendations
# to handle web specific issues like reachability, response time...

#' getURL
#'
#' @param URL   url that will be read
#' @param N.TRIES   number, how often should the function try to read the html
#'
#' @return returns output of read_html(), a html character string.
#' @import rvest
#' @import readr
getURL <- function(URL, N.TRIES=1L, 
                          read_method = c("delim", "lines", "html"),
                          skip = 0, delim = "\t", col_names = TRUE) {
    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))
    while (N.TRIES > 0L) {
        if (read_method == "delim") {
            result <- tryCatch(suppressMessages(readr::read_delim(URL, 
                                                 delim = delim, skip = skip,
                                                 col_names = col_names)), error=identity)
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
             "\n  URL: ", URL,
             "\n  error: ", conditionMessage(result))
    }
    result
}

# MHC I
netmhcI_input_template <- getURL(URL = "https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/allele.list", 
                                read_method = "delim",
                                delim = "\t",
                                col_names = c("netmhc_input", "hla_chain_name", "HLA_gene"))

# MHC II
lines <- getURL(URL = "https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/alleles_name.list",
                read_method = "lines")
lines <- str_replace_all(lines, "\t+|\\s\\s+", "\t")
netmhcII_input_template <- suppressWarnings(
    suppressMessages(readr::read_delim(lines, delim = "\t")))

all_netmhcII_template <- c(netmhcII_input_template$DR, 
                           netmhcII_input_template$`DP alpha`,
                           netmhcII_input_template$`DP beta`,
                           netmhcII_input_template$`DQ alpha`,
                           netmhcII_input_template$`DQ beta`)
all_netmhcII_template <- all_netmhcII_template[!is.na(all_netmhcII_template)]



# G and P groups

# G group
g_group <- getURL(URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt", 
                   read_method = "delim",
                   skip = 6, delim = ";", col_names = c("locus", "g_group", "g_group_name"))
g_group$locus <- str_extract(g_group$locus, pattern = "[:alnum:]+(?=\\*?)")


# P group
p_group <- getURL(URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt", 
                  read_method = "delim",
                  skip = 6, delim = ";", col_names = c("locus", "p_group", "p_group_name"))
p_group$locus <- str_extract(p_group$locus, pattern = "[:alnum:]+(?=\\*?)")

