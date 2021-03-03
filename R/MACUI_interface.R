# MACUI from https://hml.nmdp.org/mac

# ENCODE --------

assemble_encode_curlopts <- function(allele_names) {
    # check that from same locus and right format
    hla_string <- paste(allele_names, collapse="/")
    hla_string
}

#' create_encode_handle
#'
#' @param allele_names 
#'
#' @return curl handle
#' @import curl
#'
create_encode_handle <- function(allele_names) {
    h <- new_handle()
    opts_string <- assemble_encode_curlopts(allele_names)
    handle_setopt(h, copypostfields = opts_string)
    handle_setheaders(h, "Content-Type" = "text/plain")
}

#' fetch_encoded_MAC
#'
#' @param handle 
#'
#' @return curl handle fetch
#' @import curl
#'
fetch_encoded_MAC <- function(handle) {
    res <- curl_fetch_memory("https://hml.nmdp.org/mac/api/encode", 
                             handle = handle)
    res
}

#' encode_MAC 
#'
#' @param allele_list 
#'
#' @return encoded MAC 
#' @export
#' @import curl
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
    h <- new_handle()
    res <- curl_fetch_memory(url, 
                             handle = h)
    res
}

#' decode_MAC
#'
#' @param MAC 
#'
#' @return decoded MAC
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
