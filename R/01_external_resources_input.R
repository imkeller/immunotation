#### Data ####

load_mro <- function() {
    obo_file_name <- "/home/katharina/Documents/synced/science_projects/method_development/comp_immuno/immunotation/inst/extdata/mro_mod.obo"
    
    mro.obo <- ontologyIndex::get_OBO(obo_file_name, extract_tags="everything")
    
    #ontologyIndex::check(mro.obo)
    
    mro.obo
}

#' mro.obo
#' @details   \code{mro.obo}: MHC Restriction Ontology obo file.
#' @export
mro.obo <- load_mro()








# NetMHC list of valid chains
# MHC I
netmhcI_input_template <- readr::read_delim("https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/allele.list", delim = "\t",
                                col_names = c("netmhc_input", "hla_chain_name", "HLA_gene"))

# MHC II
# STILL A PROBLEM
# The mouse entry is separated by 
lines <- readr::read_lines("https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/alleles_name.list")
lines <- str_replace_all(lines, "\t+|\\s\\s+", "\t")
netmhcII_input_template <- read_delim(lines, delim = "\t")

all_netmhcII_template <- c(netmhcII_input_template$DR, 
                           netmhcII_input_template$`DP alpha`,
                           netmhcII_input_template$`DP beta`,
                           netmhcII_input_template$`DQ alpha`,
                           netmhcII_input_template$`DQ beta`)
all_netmhcII_template <- all_netmhcII_template[!is.na(all_netmhcII_template)]



# AFND
# Not sure if necessary
pop_url_root <- "http://www.allelefrequencies.net/pop6001c.asp?pop_id=%s"

freq_url_root <- "http://www.allelefrequencies.net/hla6006a.asp?hla_locus_type=%s&hla_locus=%s&hla_selection=%s&hla_population=%s&hla_sample_size_pattern=%s&hla_sample_size=%s&standard=%s"


# G and P groups

# G group
g_group <- readr::read_delim("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt", skip = 6,
                      delim = ";", col_names = c("locus", "g_group", "g_group_name"))
g_group$locus <- str_extract(g_group$locus, pattern = "[:alnum:]+(?=\\*?)")

unique(g_group$locus)


# P group
p_group <- readr::read_delim("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt", skip = 6,
                             delim = ";", col_names = c("locus", "p_group", "p_group_name"))
p_group$locus <- str_extract(p_group$locus, pattern = "[:alnum:]+(?=\\*?)")

unique(p_group$locus)

