# MACUI from https://hml.nmdp.org/mac

# ENCODE --------

assemble_encode_curlopts <- function(allele_names) {
    # check that from same locus and right format
    hla_string <- paste(allele_names, collapse="/")
    hla_string
}

#' create_encode_handle
#'
#' @param allele_names  list of HLA alleles
#'
#' @return curl handle
#'
create_encode_handle <- function(allele_names) {
    h <- curl::new_handle()
    opts_string <- assemble_encode_curlopts(allele_names)
    curl::handle_setopt(h, copypostfields = opts_string)
    curl::handle_setheaders(h, "Content-Type" = "text/plain")
}

#' fetch_encoded_MAC
#'
#' @param handle curl handle
#'
#' @return curl handle fetch
#'
fetch_encoded_MAC <- function(handle) {
    res <- curl::curl_fetch_memory("https://hml.nmdp.org/mac/api/encode", 
                             handle = handle)
    res
}

#' @title Encode MAC 
#' @description Encode a list of HLA alleles into multiple allele code (MAC). 
#' The National Marrow Donor Program (NMDP) uses 
#' [MAC](https://bioinformatics.bethematchclinical.org/hla-resources/allele-codes/allele-code-lists/) 
#' to facilitate the reporting and comparison of HLA alleles. MAC represent groups of HLA alleles and are useful when the 
#' HLA typing is ambiguous and does not allow to narrow down one single allele from a list of alleles.
#'
#' @param allele_list list of HLA alleles (e.g. c("A*01:01:01", "A*02:01:01", "A*03:01"))
#'
#' @return encoded MAC 
#' @export
#'
#' @examples
#' allele_list <- c("A*01:01:01", "A*02:01:01", "A*03:01")
#' encode_MAC(allele_list)
#' 
encode_MAC <- function(allele_list) {
    handle <- create_encode_handle(allele_list)
    curl_res <- fetch_encoded_MAC(handle)
    rawToChar(curl_res$content)
}

# DECODE -----

assemble_decode_url <- function(MAC) {
    # check that right format
    sprintf("https://hml.nmdp.org/mac/api/decode?typing=%s&expand=true",
            MAC)
}

fetch_decoded_MAC <- function(MAC) {
    url <- assemble_decode_url(MAC)
    h <- curl::new_handle()
    res <- curl::curl_fetch_memory(url, 
                             handle = h)
    res
}

#' @title Decode MAC
#' @description Decode a multiple allele code (MAC) into a list of HLA alleles.
#' #' The National Marrow Donor Program (NMDP) uses 
#' [MAC](https://bioinformatics.bethematchclinical.org/hla-resources/allele-codes/allele-code-lists/) 
#' to facilitate the reporting and comparison of HLA alleles. MAC represent groups of HLA alleles and are useful when the 
#' HLA typing is ambiguous and does not allow to narrow down one single allele from a list of alleles.
#'
#' @param MAC multiple allele code (e.g. "A*01:ATJNV")
#'
#' @return list of HLA alleles
#' @export
#'
#' @examples
#' MAC <- "A*01:ATJNV"
#' decode_MAC(MAC)
#' 
decode_MAC <- function(MAC) {
    curl_res <- fetch_decoded_MAC(MAC)
    rawToChar(curl_res$content)
}
