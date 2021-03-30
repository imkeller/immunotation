context("P and G groups")

test_that("get_G_group() and get_P_group() return error when arguments not adequate", {
    get_G_group("foo")
    get_P_group("foo")
})

test_that("get_G_group() and get_P_group() return expected results", {
    # this should return something else
    get_G_group("HLA*A*01:01:01:01")
    # this should return something else
    get_P_group("foo*")
    get_P_group("*")
    
    
    
    
    get_G_group("A*01:01:01:01")
    get_P_group("A*01:01:01:01")
})


context("Get MHC input from allele list")

test_that("get_mhcpan_input() returns error when arguments not adequate", {
    expect_error(get_mhcpan_input())
    expect_error(get_mhcpan_input(allele_list = "", mhc_class = "MHC-I"))
    expect_error(get_mhcpan_input(mhc_class = "MHC-I"))
    expect_error(get_mhcpan_input(name_list = "A*01:01", mhc_class = "foo"))
})

test_that("get_mhcpan_input() returns error when alleles not in the list", {
    # this should be an error message
    get_mhcpan_input("A*01", mhc_class = "MHC-I")
    get_mhcpan_input("A*01", mhc_class = "MHC-II")
    
    # invert mhc types
    expect_error(get_mhcpan_input("A*01:01", mhc_class = "MHC-II"))
    expect_error(get_mhcpan_input("DRB1*01:01", mhc_class = "MHC-I"))
})

test_that("get_mhcpan_input() returns expected results", {
    expect_equal(get_mhcpan_input("A*01:01", mhc_class = "MHC-I"), "HLA-A01:01")
})


context("MAC conversion")

test_that("encode_MAC() returns error when arguments not adequate", {
    expect_error(encode_MAC())
    expect_equal(encode_MAC(""), "glstring is empty")
    expect_error(encode_MAC("A*01:01", "A*01:02"))
    
    # this should return an error warning
    expect_equal(encode_MAC("foo"), "A service failures was logged.  Contact administrator for more information")
})

test_that("encode_MAC() returns expected results", {
    expect_equal(encode_MAC("A*01:01"), "A*01:01")
    expect_equal(encode_MAC(c("A*01:01", "A*01:02")), "A*01:AB")
})

test_that("decode_MAC() returns error when arguments not adequate", {
    expect_error(decode_MAC())
    expect_equal(decode_MAC(""), "Empty allele designation")
    expect_error(decode_MAC("A*01:01", "A*01:02"))
    
    # this should return an error warning
    expect_equal(decode_MAC("foo"), "No asterisk in allele designation: 'foo'")
})

test_that("decode_MAC() returns expected results", {
    expect_is(decode_MAC("A*01:01"), "character")
    expect_is(decode_MAC("A*01:AB"), "character")
})