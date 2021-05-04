context("Reading input file utilities")

test_that("getURL() returns an error on an invalid URL", {
    expect_error(getURL(URL = "https://iamnotanurl", read_method = "delim"))
    expect_error(getURL(URL = "https://iamnotanurl", read_method = "lines"))
    expect_error(getURL(URL = "https://iamnotanurl", read_method = "html"))
})

test_that("getURL() returns an error when the arguments are not matching", {
    expect_error(getURL(URL = "https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/allele.list", 
                        read_method = "iamnotamethod"))
    expect_error(getURL(URL = "https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/allele.list", 
                        read_method = "delim", delim = TRUE))
    expect_warning(getURL(URL = "https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/allele.list", 
                          read_method = "delim", 
        col_names = c("a","b","c","d")))
})


context("External input of MRO ontology")

test_that("load_mro() returns a valid ontology", {
    expect_true(ontologyIndex::check(load_mro()))
})

context("NetMHC valid allele lists")

test_that("netmhcI_input_template has the right column types", {
    # This can be updated
    expect_equal(ncol(netmhcI_input_template), 3)
    expect_is(netmhcI_input_template$netmhc_input, "character")
    expect_is(netmhcI_input_template$hla_chain_name, "character")
})

test_that("netmhcI_input_template has adequate entries", {
    expect_true("HLA-A01:01" %in% netmhcI_input_template$netmhc_input)
    expect_true("HLA-A*01:01" %in% netmhcI_input_template$hla_chain_name)
})

test_that("all_netmhcII_template has adequate entries", {
    # check that there are entries for every expected locus
    patterns_loci <- c("DRB1*", "DRB3*", "DRB4*", "DRB5*", 
                       "DPA1*", "DPB1*", "DQA1*", "DQB1*")
    observed_loci <- unique(stringr::str_extract(all_netmhcII_template, 
                                                 "^\\w+\\*"))
    common_loci <- intersect(patterns_loci, observed_loci)
    expect_equal(length(common_loci), 8)
    
    # check that a exemplary allele is present
    expect_true("DRB1*01:03" %in% all_netmhcII_template)
})

context("G and P groups external sources")

test_that("g_groups has the right column types", {
    expect_equal(ncol(g_group), 3)
    expect_is(g_group$locus, "character")
    expect_is(g_group$g_group, "character")
    expect_is(g_group$g_group_name, "character")
})

test_that("g_groups has all expected loci", {
    patterns_loci <- c("A","B","C","DPA1","DPB1","DQA1",
                       "DQB1","DRB1","DRB3","DRB4","DRB5")
    observed_loci <- unique(g_group$locus)
    common_loci <- intersect(patterns_loci, observed_loci)
    expect_equal(length(common_loci), 11)
})

test_that("g_groups has expected entries", {
    expect_true("01:01:01G" %in% g_group$g_group_name)
    expect_true("01:03:01:01/01:03:01:02/01:03:02/01:287N/01:315" %in% 
                    g_group$g_group)
})

test_that("p_groups has the right column types", {
    expect_equal(ncol(p_group), 3)
    expect_is(p_group$locus, "character")
    expect_is(p_group$p_group, "character")
    expect_is(p_group$p_group_name, "character")
})

test_that("p_groups has all expected loci", {
    patterns_loci <- c("A","B","C","DPA1","DPB1","DQA1",
                       "DQB1","DRB1","DRB3","DRB4","DRB5")
    observed_loci <- unique(p_group$locus)
    common_loci <- intersect(patterns_loci, observed_loci)
    expect_equal(length(common_loci), 11)
})

test_that("p_groups has expected entries", {
    expect_true("01:01P" %in% p_group$p_group_name)
    expect_true("01:03:01:01/01:03:01:02/01:03:02/01:315" %in% p_group$p_group)
})


