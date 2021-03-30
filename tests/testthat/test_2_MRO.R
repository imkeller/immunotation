context("ORGANISM")

test_that("get_organism() returns a list of organisms", {
    organisms <- get_valid_organisms()
    expect_is(organisms, "character")
    expect_equal(length(organisms), 23)
    expect_true("human" %in% organisms)
    expect_false("NCBITaxon:9606" %in% organisms)
    expect_false("human" %in% names(organisms))
    expect_true("NCBITaxon:9606" %in% names(organisms))
})

test_that("organism_input() distinguishes valid from non-valid organisms", {
    expect_equal(organism_input("human"), "NCBITaxon:9606")
    expect_error(organism_input("marsian"))
    expect_error(organism_input(""))
    expect_error(organism_input())
})


context("CHAIN LOOKUP TABLE")

test_that("retrieve_chain_lookup_table() returns error when arguments are not correct", {
    expect_error(retrieve_chain_lookup_table())
    expect_error(retrieve_chain_lookup_table("marsian"))
})

test_that("retrieve_chain_lookup_table() returns expected result for mouse", {
    m_df <- retrieve_chain_lookup_table("mouse")
    expect_is(m_df, "data.frame")
    
    # check that all loci are present
    expected_loci <- c("MHC class I locus","MHC class II locus","non-classical MHC locus")
    observed_loci <- unique(m_df$mhc_type_names)
    expect_equal(length(intersect(expected_loci, observed_loci)),3) 
})

test_that("retrieve_chain_lookup_table() returns expected result for human", {
    h_df <- retrieve_chain_lookup_table("human")
    expect_is(h_df, "data.frame")
    
    # check that all loci are present
    expected_loci <- c("MHC class I locus","MHC class II locus","non-classical MHC locus")
    observed_loci <- unique(h_df$mhc_type_names)
    expect_equal(length(intersect(expected_loci, observed_loci)),3) 
})

test_that("retrieve_chain_lookup_table() contains sequence information for human", {
    h_df <- retrieve_chain_lookup_table("human")
    expect_true(length(h_df$chain_sequences) > 0)
    expect_true(nchar(h_df$chain_sequences[[1]]) > 10)
})


context("SEROTYPES")

test_that("get_serotypes() returns error when arguments are not correct", {
    expect_error(get_serotypes())
    expect_error(get_serotypes(allele_list = "A*01:01"))
    expect_error(get_serotypes(allele_list = "A*01:01", mhc_type = "foo"))
    expect_error(get_serotypes(allele_list = "A*01:01", mhc_type = "MHC-I", organism = "marsian"))
    expect_error(get_serotypes(allele_list = "foo", mhc_type = "MHC-I"))
    expect_error(get_serotypes(allele_list = "foo", mhc_type = "MHC-II"))
    expect_error(get_serotypes(allele_list = c("A*01:01","foo"), mhc_type = "MHC-II"))
    
    # empty return when allele and organism not matching
    expect_equal(length(get_serotypes(allele_list = "A*01:01", mhc_type = "MHC-I", organism = "mouse")), 0)
    
    # correct query
    expect_equal(get_serotypes(allele_list = "A*01:01", mhc_type = "MHC-I")[[1]], "HLA-A1 serotype")
})
