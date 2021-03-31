
assemble_frequency_lat_long_df <- function(global_alleles) {
    
    global_alleles$sample_size <- as.integer(stringr::str_remove_all(global_alleles$sample_size, pattern = ","))
    global_alleles$allele_frequency <- as.double(global_alleles$allele_frequency)
    
    global_alleles <- global_alleles[global_alleles$allele_frequency != 0.000,]
    global_alleles <- tidyr::drop_na(global_alleles)
    
    population_details <- query_population_detail(unique(as.numeric(global_alleles$population_id)))
    
    population_details <- tidyr::separate(population_details, Latitude, into = c("lata", "mina", "dira"))
    population_details$lat <- ifelse(population_details$dira == "N", 
                                     as.numeric(population_details$lata) + as.numeric(population_details$mina)/60,
                                    -as.numeric(population_details$lata) - as.numeric(population_details$mina)/60)
    
    population_details <- tidyr::separate(population_details, Longitude, into = c("longo", "mino", "diro"))
    population_details$long <- ifelse(population_details$diro == "E", 
                                     as.numeric(population_details$longo) + as.numeric(population_details$mino)/60,
                                     -as.numeric(population_details$longo) - as.numeric(population_details$mino)/60)
    
    df_lat_long <- merge(global_alleles, population_details, by = "population_id")
    df_lat_long
}

#' @title Plotting allele frequencies
#' @description \code{plot_allele_frequency} Generate a World map displaying the frequency of a given table of HLA alleles. Use the function
#' \link{query_allele_frequencies} to generate a table with allele frequencies.
#'
#' @param allele_frequency returned by \link{query_allele_frequencies}
#'
#' @return ggplot2 object displaying the allele frequencies on a world map.
#' @export
#' @import maps
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
    
    WorldData <- ggplot2::map_data("world")[ggplot2::map_data("world")$region != "Antarctica",]
    
    ggplot2::ggplot() +
        ggplot2::geom_map(data = WorldData, map = WorldData,
                 # removing x and y here makes the plot scale to the dots
                 ggplot2::aes(group = group, map_id=region), fill = "darkgrey") +
        ggplot2::geom_point(data = df_lat_long,
                 ggplot2::aes(x = as.numeric(long), y = as.numeric(lat), color = allele_frequency),
                 size = 5) +
        # scale x and y so that its uniform for all plots
        ggplot2::xlim(min(WorldData$long)+30, max(WorldData$long)-10) +
        ggplot2::ylim(min(WorldData$lat)+10, max(WorldData$lat)-10) +
        ggplot2::theme(legend.position=c(0.14,0.23),
              axis.title.x=ggplot2::element_blank(), 
              axis.title.y=ggplot2::element_blank())
}
