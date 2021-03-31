context("Querying allele frequencies")

test_that("query_allele_frequencies() return error when arguments not adequate", {
    expect_error(query_allele_frequencies())
    expect_error(query_allele_frequencies(hla_locus = NA, hla_selection = NA, hla_population = NA,
                                          hla_country = NA, hla_region = NA, hla_ethnic = NA,
                                          hla_sample_size_pattern = NA, hla_sample_size = NA,  standard = "a"))
    expect_error(query_allele_frequencies(hla_locus = "foo", hla_selection = NA, hla_population = NA,
                                          hla_country = NA, hla_region = NA, hla_ethnic = NA,
                                          hla_sample_size_pattern = NA, hla_sample_size = NA,  standard = "a"))
    expect_error(query_allele_frequencies(hla_locus = NA, hla_selection = "foo", hla_population = NA,
                                          hla_country = NA, hla_region = NA, hla_ethnic = NA,
                                          hla_sample_size_pattern = NA, hla_sample_size = NA,  standard = "a"))
    expect_error(query_allele_frequencies(hla_locus = NA, hla_selection = NA, hla_population = "foo",
                                          hla_country = NA, hla_region = NA, hla_ethnic = NA,
                                          hla_sample_size_pattern = NA, hla_sample_size = NA,  standard = "a"))
    expect_error(query_allele_frequencies(hla_locus = NA, hla_selection = NA, hla_population = NA,
                                          hla_country = "foo", hla_region = NA, hla_ethnic = NA,
                                          hla_sample_size_pattern = NA, hla_sample_size = NA,  standard = "a"))
    expect_error(query_allele_frequencies(hla_locus = NA, hla_selection = NA, hla_population = NA,
                                          hla_country = NA, hla_region = "foo", hla_ethnic = NA,
                                          hla_sample_size_pattern = NA, hla_sample_size = NA,  standard = "a"))
    expect_error(query_allele_frequencies(hla_locus = NA, hla_selection = NA, hla_population = NA,
                                          hla_country = NA, hla_region = NA, hla_ethnic = "foo",
                                          hla_sample_size_pattern = NA, hla_sample_size = NA,  standard = "a"))
    expect_error(query_allele_frequencies(hla_locus = NA, hla_selection = NA, hla_population = NA,
                                          hla_country = NA, hla_region = NA, hla_ethnic = NA,
                                          hla_sample_size_pattern = "foo", hla_sample_size = NA,  standard = "a"))
    expect_error(query_allele_frequencies(hla_locus = NA, hla_selection = NA, hla_population = NA,
                                          hla_country = NA, hla_region = NA, hla_ethnic = NA,
                                          hla_sample_size_pattern = NA, hla_sample_size = "foo",  standard = "a"))
    expect_error(query_allele_frequencies(hla_locus = NA, hla_selection = NA, hla_population = NA,
                                          hla_country = NA, hla_region = NA, hla_ethnic = NA,
                                          hla_sample_size_pattern = NA, hla_sample_size = NA,  standard = "foo"))
})

test_that("query_allele_frequencies() returns expected results", {
    sel1 <- query_allele_frequencies(hla_selection = "A*02:01", 
                                     hla_sample_size_pattern = "bigger_than", 
                                     hla_sample_size = 10000, 
                                     standard="g")
    expect_equal(ncol(sel1), 5)
})

context("Building overall allele groups")

test_that("build_allele_group() return error when arguments not adequate", {
    expect_error(build_allele_group())
    expect_error(build_allele_group("HLA01:01"))
    expect_error(build_allele_group("foo"))
    expect_error(build_allele_group("*"))
    expect_error(build_allele_group("DRB1*01"))
})

test_that("build_allele_group() returns expected results", {
    # this is not yet working
    expect_error(build_allele_group("A*01"))
    
    l1 <- build_allele_group("A*01:01")
    expect_is(l1, "character")
    expect_true(length(l1) > 3)
    
    l2 <- build_allele_group("A*01:01:01:01")
    expect_is(l2, "character")
    expect_true(length(l2) > 3)
    
    l3 <- build_allele_group("HLA-A*01:01")
    expect_is(l3, "character")
    expect_true(length(l3) > 3)
    
    l4 <- build_allele_group("DRB1*01:01:01")
    expect_is(l4, "character")
    expect_true(length(l4) > 3)
    
    l5 <- build_allele_group("HLA-DRB1*01:01:01")
    expect_is(l5, "character")
    expect_true(length(l5) > 3)
})
    
context("Querying haplotype frequency")

test_that("build_allele_group() return error when arguments not adequate", {
    expect_error(query_haplotype_frequencies())
    expect_error(query_haplotype_frequencies(hla_selection = NA, hla_population = NA, hla_country = NA,
                                    hla_region = NA, hla_ethnic = NA, hla_sample_size_pattern = NA,
                                    hla_sample_size = NA))
    
    # not behaving as expected
    expect_error(query_haplotype_frequencies(hla_selection = "foo", hla_population = NA, hla_country = NA,
                                             hla_region = NA, hla_ethnic = NA, hla_sample_size_pattern = NA,
                                             hla_sample_size = NA))  
    expect_error(query_haplotype_frequencies(hla_selection = "", hla_population = NA, hla_country = NA,
                                             hla_region = NA, hla_ethnic = NA, hla_sample_size_pattern = NA,
                                             hla_sample_size = NA)) 
    expect_error(query_haplotype_frequencies(hla_selection = "HLA-A", hla_population = NA, hla_country = NA,
                                             hla_region = NA, hla_ethnic = NA, hla_sample_size_pattern = NA,
                                             hla_sample_size = NA)) 
    
    expect_error(query_haplotype_frequencies(hla_selection = NA, hla_population = "foo", hla_country = NA,
                                             hla_region = NA, hla_ethnic = NA, hla_sample_size_pattern = NA,
                                             hla_sample_size = NA))
    expect_error(query_haplotype_frequencies(hla_selection = NA, hla_population = NA, hla_country = "foo",
                                             hla_region = NA, hla_ethnic = NA, hla_sample_size_pattern = NA,
                                             hla_sample_size = NA))
    expect_error(query_haplotype_frequencies(hla_selection = NA, hla_population = NA, hla_country = NA,
                                             hla_region = "foo", hla_ethnic = NA, hla_sample_size_pattern = NA,
                                             hla_sample_size = NA))
    expect_error(query_haplotype_frequencies(hla_selection = NA, hla_population = NA, hla_country = NA,
                                             hla_region = NA, hla_ethnic = "foo", hla_sample_size_pattern = NA,
                                             hla_sample_size = NA))
    expect_error(query_haplotype_frequencies(hla_selection = NA, hla_population = NA, hla_country = NA,
                                             hla_region = NA, hla_ethnic = NA, hla_sample_size_pattern = "foo",
                                             hla_sample_size = NA))
    expect_error(query_haplotype_frequencies(hla_selection = NA, hla_population = NA, hla_country = NA,
                                             hla_region = NA, hla_ethnic = NA, hla_sample_size_pattern = NA,
                                             hla_sample_size = "foo"))
})

test_that("build_allele_group() returns expected results", {
    df <- query_haplotype_frequencies(hla_selection = c("A*02:01", "B*", "C*"),
                                hla_region = "Europe")
    expect_equal(ncol(df), 5)
    
    expect_is(query_haplotype_frequencies(hla_selection = "A*01:01", hla_population = NA, hla_country = NA,
                                             hla_region = NA, hla_ethnic = NA, hla_sample_size_pattern = NA,
                                             hla_sample_size = NA), "data.frame")  
})


context("Plotting allele frequencies")

test_that("plot_allele_frequency() return error when arguments not adequate", {
    expect_error(plot_allele_frequency())
    
    # not behaving as expected
    expect_error(plot_allele_frequency(data.frame()))
})

test_that("plot_allele_frequency() returns expected results", {
    sel <- query_allele_frequencies(hla_selection = "A*02:01", 
                                    hla_sample_size_pattern = "bigger_than", 
                                    hla_sample_size = 10000, 
                                    standard="g")
    expect_is(plot_allele_frequency(sel), "ggplot")
})

context("Querying population details")

test_that("query_population_detail() return error when arguments not adequate", {
    expect_error(query_population_detail())
    expect_error(query_population_detail("foo"))
    expect_error(query_population_detail("3725"))
    
    # not behaving as expected
    expect_error(query_population_detail(1))
    expect_error(query_population_detail(c(1111)))
})

test_that("query_population_detail() returns expected results", {
    df <- query_population_detail(3725)
    expect_equal(ncol(df), 13)
    expect_equal(nrow(df), 1)
    expect_equal(df$population_id, 3725)
})