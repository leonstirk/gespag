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
#' @return list
#' @export
makeGrid <- function(df, lon_colname, lat_colname, resolution) {
  dat <- df %>% dplyr::select(lon_colname, lat_colname)
  origin <- dat %>% sapply(min)
  cellsize <- dat %>% sapply(., function(x) { diff(range(x))/resolution })
  return(list('origin' = origin, 'cellsize' = cellsize))
}

#' Create hex binnings as a SpatialPolygonsDataFrame
#' @param df data.frame
#' @param lon_colname character
#' @param lat_colname character
#' @param resolution integer
#' @param crs CRS
#' @return SpatialPolygons
#' @export
makeHex <- function(df, lon_colname, lat_colname, resolution, crs) {
  sp::coordinates(df) <- ~ lon + lat
  sp::proj4string(df) <- crs
  cellsize <- diff(range(df@coords[,lon_colname]))/resolution
  df@bbox[,'max'] <- df@bbox[,'max'] + cellsize
  df@bbox[,'min'] <- df@bbox[,'min'] - cellsize
  hex_pts <- sp::spsample(df, type = 'hexagonal', cellsize = cellsize, bb = df@bbox)
  hex_poly <- sp::HexPoints2SpatialPolygons(hex_pts)
  return(hex_poly)
}

#' Attach hexbins
#' @param df data.frame
#' @param make_hex SpatialPolygons
#' @param lon_colname character
#' @param lat_colname character
#' @return data.frame
#' @export
attachHexVars <- function(df, make_hex, lon_colname, lat_colname) {
  df_sp <- df
  crs <- make_hex@proj4string
  names(df_sp)[which(names(df_sp) == lon_colname)] <- 'lon'
  names(df_sp)[which(names(df_sp) == lat_colname)] <- 'lat'
  sp::coordinates(df_sp) <- ~ lon + lat
  sp::proj4string(df_sp) <- crs
  df$cell_id <- sp::over(df_sp, make_hex) %>% as.character() %>% as.factor()
  return(df %>% dplyr::as_tibble())
}
 #' Attach custom shapefile cell ids to data
#' @param df data.frame
#' @param shp_poly SpatialPolygons
#' @param poly_id_colname character
#' @param lon_colname character
#' @param lat_colname character
#' @return data.frame
#' @export
attachShpVars <- function(df, shp_poly, poly_id_colname, lon_colname, lat_colname) {
  df_sp <- df
  crs <- shp_poly@proj4string
  names(df_sp)[which(names(df_sp) == lon_colname)] <- 'lon'
  names(df_sp)[which(names(df_sp) == lat_colname)] <- 'lat'
  sp::coordinates(df_sp) <- ~ lon + lat
  sp::proj4string(df_sp) <- crs
  ptinpoly <- sp::over(df_sp, shp_poly)
  df$cell_id <- ptinpoly[[poly_id_colname]] %>% as.factor()
  return(df %>% dplyr::as_tibble())
}

