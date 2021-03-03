
assemble_frequency_lat_long_df <- function(global_alleles) {

    global_alleles <- global_alleles %>%
        dplyr::mutate(sample_size = str_remove_all(sample_size, pattern = ",") %>% as.integer(),
                      allele_frequency = as.double(allele_frequency)) %>%
        dplyr::filter(allele_frequency != 0.000) %>%
        #dplyr::filter(sample_size > 200) %>%
        tidyr::drop_na() %>%
        tidyr::as_tibble()
    
    population_details <- query_population_detail(unique(as.numeric(global_alleles$population_id)))
    
    pop_long_lat <- population_details %>%
        separate(Latitude, into = c("lata", "mina", "dira")) %>%
        separate(Longitude, into = c("longo", "mino", "diro")) %>%
        dplyr::mutate(lat = ifelse(dira == "N", as.numeric(lata) + as.numeric(mina)/60,
                                   -as.numeric(lata) - as.numeric(mina)/60)) %>%
        dplyr::mutate(long = ifelse(diro == "E", as.numeric(longo) + as.numeric(mino)/60,
                                    -as.numeric(longo) - as.numeric(mino)/60))
    
    df_lat_long <- merge(global_alleles, pop_long_lat, by = "population_id")
}

#' plot_allele_frequency
#'
#' @param allele_frequency 
#'
#' @return ggplot
#' @export
#'
#' @examples
#' 
#' # select frequency of given allele
#' sel_allele_freq <- query_allele_frequencies(hla_selection = "A*02:01", 
#' hla_sample_size_pattern = "bigger_than", 
#' hla_sample_size = 10000, standard="g")
#' 
#' plot_allele_frequency(sel_allele_freq)
#' 
plot_allele_frequency <- function(allele_frequency) {
    
    df_lat_long <- assemble_frequency_lat_long_df(allele_frequency)
    
    WorldData <- map_data("world") %>% filter(region != "Antarctica") %>% fortify
    
    ggplot() +
        geom_map(data = WorldData, map = WorldData,
                 # removing x and y here makes the plot scale to the dots
                 aes( group = group, map_id=region), fill = "darkgrey") +
        geom_point(data = df_lat_long,
                   aes(x = as.numeric(long), y = as.numeric(lat), color = allele_frequency),
                   size = 5) + theme(legend.position=c(0.14,0.23)) +
        # scale x and y so that its uniform for all plots
        xlim(min(WorldData$long)+30, max(WorldData$long)-10) +
        ylim(min(WorldData$lat)+10, max(WorldData$lat)-10) +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank())
}
