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
lines <- stringr::str_replace_all(lines, "\t+|\\s\\s+", "\t")
netmhcII_input_template <- suppressWarnings(
    suppressMessages(readr::read_delim(lines, delim = "\t")))

all_netmhcII_template <- c(netmhcII_input_template$DR, 
    netmhcII_input_template$`DP alpha`, netmhcII_input_template$`DP beta`,
    netmhcII_input_template$`DQ alpha`, netmhcII_input_template$`DQ beta`)
# all_netmhcII_template is an internal variable containing list of valid 
# NetMHCIIpan input alleles
all_netmhcII_template <- all_netmhcII_template[!is.na(all_netmhcII_template)]


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


#
# Parse JSON-file arcasHLA output
#

#' parse_json_arcasHLA_output
#' 
#' @param path_to_arcasHLA Indicates the path to file for reading
#' @param path_to_peptide Indicates the path to peptide of interest
#' @param output_bash Indicates output file name for bash script
#' @param mhc_class_to_analyze MHC class to analyze. Possible is MHC-I, MHC-II or both
#' 
#' @return returns list with MHCI/MHCII alleles ready for NetMHCpan
#' @export
#' @importFrom jsonlite fromJSON
#' @import purrr 
#' @import stringr
#' @examples
#' parse_json_arcasHLA_output(path_to_arcasHLA = system.file("extdata", ".", 
#' package = "immunotation"), path_to_peptide = system.file("extdata",   
#' "virus_example.fasta", package = "immunotation"), 
#' output_bash = "example.genotype.sh", 
#' mhc_class_to_analyze = "both")
parse_json_arcasHLA_output <- function(path_to_arcasHLA, path_to_peptide,
                                       output_bash, mhc_class_to_analyze){
    if(!grepl('.json', path_to_arcasHLA)) {
        stop("not JSON-format\n")
    }
    p1 <- jsonlite::fromJSON(txt=path_to_arcasHLA)
    if(mhc_class_to_analyze == "MHC-I"){
        MHCI <- unlist(purrr::map(p1[c('A', 'B', 'C')], 
                                  .f = ~get_mhcpan_input(.x, 'MHC-I')), use.names = FALSE)
        if((!any(grepl('DQA', names(p1))))||(!any(grepl('DQB', names(p1))))) {
            warning("DQA or DQB is missing\n")
        }
        if((!any(grepl('DPB', names(p1))))||(!any(grepl('DPA', names(p1))))) {
            warning("DPA or DPB is missing\n")
        }
        output_xlsI <- paste(stringr::str_sub(tail(unlist(stringr::str_split(path_to_arcasHLA, '/')), 
                                                   n=1), 1, - 6), '_netMHCI.xls',sep='')
        f <- generate_bash_output(output_bash=output_bash, 
                                  path_to_peptide=path_to_peptide, MHCI=MHCI, 
                                  output_xlsI=output_xlsI)
    } else if (mhc_class_to_analyze == "MHC-II") {
        MHCII <- c(get_mhcpan_input(unlist(p1[grepl('DQ', names(p1))]),'MHC-II'),
                   get_mhcpan_input(unlist(p1[grepl('DR', names(p1))]),'MHC-II'),
                   get_mhcpan_input(unlist(p1[grepl('DP', names(p1))]),'MHC-II'))
        output_xlsII <- paste(stringr::str_sub(tail(unlist(stringr::str_split(path_to_arcasHLA, '/')), 
                                                    n=1), 1, - 6), '_netMHCII.xls',sep='')
        f <- generate_bash_output(output_bash=output_bash, 
                                  path_to_peptide=path_to_peptide, 
                                  MHCII=MHCII, output_xlsII=output_xlsII)
    } else if (mhc_class_to_analyze == "both") {
        MHCI <- unlist(purrr::map(p1[c('A', 'B', 'C')], 
                                  .f = ~get_mhcpan_input(.x, 'MHC-I')), use.names = FALSE)
        if((!any(grepl('DQA', names(p1))))||(!any(grepl('DQB', names(p1))))) {
            warning("DQA or DQB is missing\n")
        }
        if((!any(grepl('DPB', names(p1))))||(!any(grepl('DPA', names(p1))))) {
            warning("DPA or DPB is missing\n")
        }
        MHCII <- c(get_mhcpan_input(unlist(p1[grepl('DQ', names(p1))]),'MHC-II'),
                   get_mhcpan_input(unlist(p1[grepl('DR', names(p1))]),'MHC-II'),
                   get_mhcpan_input(unlist(p1[grepl('DP', names(p1))]),'MHC-II'))
        output_xlsI <- paste(stringr::str_sub(tail(unlist(stringr::str_split(path_to_arcasHLA, '/')), 
                                                   n=1), 1, - 6), '_netMHCI.xls',sep='')
        output_xlsII <- paste(stringr::str_sub(tail(unlist(stringr::str_split(path_to_arcasHLA, '/')), 
                                                    n=1), 1, - 6), '_netMHCII.xls',sep='')
        f <- generate_bash_output(output_bash=output_bash, 
                                  path_to_peptide=path_to_peptide, 
                                  MHCI=MHCI, MHCII=MHCII,
                                  output_xlsI=output_xlsI, 
                                  output_xlsII=output_xlsII)
    } else{
        stop(mhc_class_to_analyze, 
             'is not the right parameter for mhc_class_to_analyze')
    }
    return(list(list('MHCI'=MHCI,
                'MHCII'=MHCII)))
}

#
# Jenerate bash-script for netMHC running
#

#' generate_bash_output
#' 
#' @param output_bash Indicate name for output file
#' @param path_to_peptide Indicate path to analyzed peptide
#' @param MHCI List with MHCI alleles
#' @param MHCII List with MHCII alleles
#' @param output_xlsI Name for xls output for MHCI
#' @param output_xlsII Name for xls output for MHCII
#' 
#' @return Generates output bash-script ready to execute from cmd
#' @keywords internal
#' @import stringr
generate_bash_output<- function(output_bash, path_to_peptide, MHCI, MHCII,
                                output_xlsI, output_xlsII){
    if(missing(MHCII)){
        if(!file.exists(output_bash))  {
            file_bash_script<-file(output_bash, 'w')
            writeLines(paste("netMHCpan -f", path_to_peptide, "-BA -xls -a",
                             paste(MHCI, collapse=","),"-xlsfile", output_xlsI,sep=' '), 
                       file_bash_script)
            close(file_bash_script)
        }
        else {
            file_bash_script<-file(output_bash, 'a')
            writeLines(paste("netMHCpan -f", path_to_peptide, "-BA -xls -a",
                             paste(MHCI, collapse=","),"-xlsfile", output_xlsI,sep=' '), 
                       file_bash_script)
            close(file_bash_script)
    }
    } else if(missing(MHCI)){
        if(!file.exists(output_bash))  {
            file_bash_script<-file(output_bash, 'w')
            writeLines(paste("netMHCIIpan -f", path_to_peptide, "-BA -xls -a",
                             paste(MHCII, collapse=","),"-xlsfile", output_xlsII,sep=' '), 
                       file_bash_script)
            close(file_bash_script)
        }
        else {
            file_bash_script<-file(output_bash, 'a')
            writeLines(paste("netMHCIIpan -f", path_to_peptide, "-BA -xls -a",
                             paste(MHCII, collapse=","),"-xlsfile", output_xlsII,sep=' '), 
                       file_bash_script)
            close(file_bash_script)
        }
    } else{
        if(!file.exists(output_bash))  {
            file_bash_script<-file(output_bash, 'w')
            writeLines(paste("netMHCpan -f", path_to_peptide, "-BA -xls -a",
                             paste(MHCI, collapse=","),"-xlsfile", output_xlsI,sep=' '), 
                       file_bash_script)
            writeLines(paste("netMHCIIpan -f", path_to_peptide, "-BA -xls -a",
                             paste(MHCII, collapse=","),"-xlsfile", output_xlsII,sep=' '), 
                       file_bash_script)
            close(file_bash_script)
        }
        else {
            file_bash_script<-file(output_bash, 'a')
            writeLines(paste("netMHCpan -f", path_to_peptide, "-BA -xls -a",
                             paste(MHCI, collapse=","),"-xlsfile", output_xlsI,sep=' '), 
                       file_bash_script)
            writeLines(paste("netMHCIIpan -f", path_to_peptide, "-BA -xls -a",
                             paste(MHCII, collapse=","),"-xlsfile", output_xlsII,sep=' '), 
                       file_bash_script)
            close(file_bash_script)  
        }
    }
}

#
# Process netMHCpan/netMHCIIpan output
#

#' read_netMHCpan_output
#' 
#' @param path_to_file Indicate path to result file
#' @param mhc_class_to_analyze MHC class to analyze. Possible is MHC-I, MHC-II or both
#' 
#' @return returns processed table
#' @export
#' @import readr
#' @import purrr
#' @examples
#' path_to_file <- system.file("extdata", "example1.genotype_netMHCII.xls", package = "immunotation")
#' read_netMHCpan_output(path_to_file = path_to_file, mhc_class_to_analyze = "MHC-I")
read_netMHCpan_output<- function(path_to_file, mhc_class_to_analyze){
    first_line <- strsplit(readr::read_lines(path_to_file, n_max = 1), '\\t+')[[1]]
    if(mhc_class_to_analyze == "MHC-I") {
        updating_column_names <- c(rep(first_line[1], 3), 
                                   unlist(purrr::map(first_line[2:length(first_line)], 
                                                     .f = ~rep(.x,6))))
    } else if (mhc_class_to_analyze == "MHC-II") {
        updating_column_names <- c(rep(first_line[1], 4), 
                                   unlist(purrr::map(first_line[2:length(first_line)], 
                                                     .f = ~rep(.x,5))))
    } else {stop("mhc_class_to_analyze schould be MHC-I or MHC-II.")}
    
    processed_table <- readr::read_delim(path_to_file, delim = "\t", skip=1,col_names = FALSE)
    colnames(processed_table) <- stringr::str_c(unlist(processed_table[1,]), 
                                        sep='.', updating_column_names)
    processed_table <- processed_table[2:nrow(processed_table), ]
    processed_table$NB. <- as.integer(processed_table$NB.)
    processed_table
}

