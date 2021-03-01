# MHC Restriction ontology
obo_file_name <- system.file("extdata", "mro_mod.obo", package="immunontology")
obo_file_name <- "/home/katharina/Documents/synced/science_projects/method_development/comp_immuno/MHC_ontology/immunontology/inst/extdata/mro_mod.obo"

mro.obo <- ontologyIndex::get_OBO(obo_file_name, extract_tags="everything")

ontologyIndex::check(mro.obo)



# NetMHC list of valid chains
# MHC I
lines_mhcI <- readr::read_delim("https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/allele.list", delim = "\t",
                                col_names = c("netmhc_input", "hla_chain_name", "HLA_gene"))

# MHC II
# STILL A PROBLEM
# The mouse entry is separated by 
lines <- readr::read_lines("https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/alleles_name.list")
colnames <- str_split(lines[1], pattern = "\\s\\s++")[[1]][1:5]
lines <- lines[-1]
lines_split <- str_split(lines, pattern = "\t|\\s\\s++")

# remove the mouse entries in line 1-3
for (i in seq.int(1,2)) {
    lines_split[[i]] <- lines_split[[i]][1:5]
}




# AFND
# Not sure if necessary
pop_url_root <- "http://www.allelefrequencies.net/pop6001c.asp?pop_id=%s"

freq_url_root <- "http://www.allelefrequencies.net/hla6006a.asp?hla_locus_type=%s&hla_locus=%s&hla_selection=%s&hla_population=%s&hla_sample_size_pattern=%s&hla_sample_size=%s&standard=%s"