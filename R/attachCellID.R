#' Attach grid variables to dataset
#' @param df data.frame
#' @param make_grid list
#' @param lon_colname character
#' @param lat_colname character
#' @return data.frame
#' @export
## df <- df %>% attachGridVars('lon', 'lat', 10)
attachGridVars <- function(df, make_grid, lon_colname, lat_colname) {
  dat <- df %>% dplyr::select(lon_colname, lat_colname)
  origin <- make_grid$origin
  cellsize <- make_grid$cellsize
  cell_id <- dat %>% cellXY(origin, cellsize) %>% cantorPair() %>% as.factor()
  cell_coords <- dat %>% cellXY(origin, cellsize) %>% cellCoords(origin, cellsize)
  out <- list('cell_id' = cell_id, 'cell_lon' = cell_coords[,1], 'cell_lat' = cell_coords[,2]) %>% cbind(df) %>% dplyr::as_tibble()
  return(out)
}

#' Make grid object
#' @param df data.frame
#' @param lon_colname character
#' @param lat_colname character
#' @param resolution integer
#' @param square boolean
#' @return list
#' @export
makeGrid <- function(df, lon_colname, lat_colname, resolution, square = FALSE) {
  dat <- df %>% dplyr::select(lon_colname, lat_colname)
  origin <- dat %>% sapply(min)
  cellsize <- dat %>% sapply(., function(x) { diff(range(x))/resolution })
  if(square) { cellsize[which(cellsize == max(cellsize))] <- min(cellsize) }
  return(list('origin' = origin, 'cellsize' = cellsize))
}

#' Create hex binnings as an sf object
#' @param df data.frame
#' @param lon_colname character
#' @param lat_colname character
#' @param resolution integer
#' @param crs CRS
#' @return SpatialPolygons
#' @export
makeHex <- function(df, lon_colname, lat_colname, resolution, crs) {
  df_sf <- sf::st_as_sf(x = as.data.frame(df), coords = c(lon_colname, lat_colname), crs = crs)
  cellsize <- diff(range(df[,lon_colname]))/resolution
  a <- sf::st_bbox(df_sf) %>% matrix(2,2)
  a[,2] <- a[,2] + cellsize
  a[,1] <- a[,1] - cellsize
  rownames(a) <- c('lon', 'lat')
  colnames(a) <- c('min', 'max')
  hex_poly <- sf::st_make_grid(df_sf, cellsize = cellsize, square = FALSE)
  return(hex_poly)
}

 #' Attach custom shapefile cell ids to data
#' @param df data.frame
#' @param shp_poly SpatialPolygons
#' @param lon_colname character
#' @param lat_colname character
#' @return data.frame
#' @export
attachShpVars <- function(df, shp_poly, lon_colname, lat_colname) {
  crs <- sf::st_crs(shp_poly)
  df_sf <- sf::st_as_sf(x = as.data.frame(df), coords = c(lon_colname, lat_colname), crs = crs)
  names(df_sf)[which(names(df_sf) == lon_colname)] <- 'lon'
  names(df_sf)[which(names(df_sf) == lat_colname)] <- 'lat'
  df$cell_id <- unlist(lapply(sf::st_intersects(df_sf, shp_poly), function(x) { x[1] })) %>% as.character %>% as.factor
  return(df %>% dplyr::as_tibble())
}

