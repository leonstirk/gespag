#' see the grid overlay for a gdf created with a given resolution
#' @param gdf data.frame
#' @param map ggmap
#' @param fill_value character
#' @return ggplot
#' @export
seeGrid <- function(gdf, map, fill_value = 'ln_sale_price') {

  ## Visualise grid
  grid_data <- as.data.frame(t(matrix(unlist(by(gdf, gdf$cell_id, function(x) {
    c(x$cell_lon[1], x$cell_lat[1], mean(as.numeric(x[[fill_value]])), count(x))
  })), nrow = 4)))

  names(grid_data) <- c('lon', 'lat', 'fill', 'count')

  return(
    ggmap::ggmap(map) +
      ggplot2::geom_tile(data = grid_data, aes(x = lon, y = lat, fill = fill, alpha = count)) +
      viridis::scale_fill_viridis()
  )
}

#' see the shp overlay for a gdf created with a given resolution
#' @param gdf data.frame
#' @param shp sp
#' @param map ggmap
#' @param fill_value character
#' @return ggplot
#' @export
seeShp <- function(gdf, shp, map, fill_value = 'ln_sale_price') {

  shp_data <- gdf %>% group_by(cell_id) %>%
    summarise(fill = mean(as.numeric(get(fill_value))), count = length(as.numeric(get(fill_value))))
  shp <- st_as_sf(shp)
  shp$id <- seq(1:nrow(shp)) %>% as.character %>% as.factor
  plt <- dplyr::left_join(shp , shp_data, by=c("id"="cell_id"))

  quartz()
  return(
    ggmap::ggmap(map) +
      ggplot2::geom_sf(data = stats::na.omit(plt), aes(fill = fill, alpha = count), lwd = 0.2, inherit.aes = F) +
      viridis::scale_fill_viridis() +
      ggplot2::scale_alpha_continuous()
  )
}


#' Plot the organism evolution
#' @param organism list
#' @return plot
#' @export
seeEvolution <- function(organism) {

  evolution <- organism[['evolution']]
  theoretical_max <- organism[['other']][['theoretical_max']]
  baseline <- organism[['other']][['baseline']]

  plot(evolution$train ~ evolution$generation,
       type = 'l',
       col = 'red',
       ylim = c(baseline, theoretical_max),
       ylab = 'r2',
       xlab = 'generation')
  with(evolution, lines(test ~ generation, col = 'blue', lty = 2))
  with(evolution, lines(validate ~ generation, col = 'green', lty = 2))
  abline(h = baseline)
  abline(h = theoretical_max)
}
